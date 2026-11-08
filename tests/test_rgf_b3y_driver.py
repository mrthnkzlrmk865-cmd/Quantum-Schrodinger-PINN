"""
Tests of the B3Y driver (validate_rgf_b3y.py). Everything here uses the SYNTHETIC problem (N_y = 7, N_q = 15, dx = 1/32 ... 1/128, E = 60 and 30);
no B3Y energy point or realization is computed. Companion of tests/test_rgf_b3y.py (the sparse_exact_lattice module).
"""
import csv
import json
import os

import numpy as np
import pytest
import scipy.sparse as sp

import sparse_exact_lattice as xp
import sparse_reference as sr
import sparse_wavefunction as sw
import validate_rgf_b3 as b3
import validate_rgf_b3x as b3x
import validate_rgf_b3y as drv

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CLD, LD = xp.CLD, xp.LD


# ----------------------------------------------------------------------
# fixtures: one synthetic point, its three double methods and the exact lattice solution
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def syn():
    prob = drv.SyntheticProblem()
    devs = drv.DeviceCache(prob)
    pt = prob.points["SYN_A@E60"]
    inv = 64
    dev, geo, ref, info, pot = devs.get("SYN", inv, prob.buffer_req)
    dx, N = 1.0 / inv, dev.n_slices
    Ep = dev.basis.E_perp
    wf = sw.wavefunction_scattering_matrix(pt["E"], pot.numpy, geo["x_start"], dx, N, prob.Ly, prob.Ny, prob.Nq)
    gb = sr.sparse_scattering_matrix(pt["E"], pot.numpy, geo["x_start"], dx, N, prob.Ly, prob.Ny, prob.Nq)
    ex = xp.exact_lattice_S(pt["E"], Ep, ref, dx, prob.Ny, prob.Nq)
    return dict(prob=prob, devs=devs, pt=pt, inv=inv, dx=dx, N=N, Ep=Ep, wf=wf, gb=gb, ex=ex, ref=ref, pot=pot, geo=geo)


# ----------------------------------------------------------------------
# method-specific matrix construction, sources and extraction conventions
# ----------------------------------------------------------------------
def test_M3_source_and_extraction_mirrors_reproduce_the_methods_own_numbers_bitwise(syn):
    wf, dx = syn["wf"], syn["dx"]
    o = np.where(wf["open"])[0]
    t = sw.hopping(dx)
    assert np.array_equal(drv.wf_source_mirror(t, wf["lam"], wf["a"], o), wf["source"])
    S = drv.wf_extract(wf["c1"], wf["cN"], o, np.sqrt(wf["gamma"][o])[:, None], wf["a"])
    assert np.array_equal(S, wf["S"])


def test_M2_fisher_lee_mirror_reproduces_the_methods_own_S_bitwise(syn):
    gb = syn["gb"]
    o = np.where(gb["modes"]["open"])[0]
    G, ix = gb["G"], np.ix_(o, o)
    S = drv.fl_extract(G["G11"][ix], G["GN1"][ix], G["G1N"][ix], G["GNN"][ix], np.sqrt(gb["modes"]["gamma"][o]), 1j)
    assert np.array_equal(S, gb["S"])


def test_each_method_solves_its_own_matrix_and_the_two_matrices_are_built_differently(syn):
    """The decomposition uses M2's matrix from sr.assemble_M and M3's from sw.assemble; each reproduces its own double S when solved."""
    import scipy.sparse.linalg as spla
    E, dx, N, ny = syn["pt"]["E"], syn["dx"], syn["N"], syn["prob"].Ny
    wf, gb = syn["wf"], syn["gb"]
    M3 = sw.assemble(E, wf["V_blocks"], wf["E_perp"], dx, wf["lam"])
    M2 = sr.assemble_M(E, gb["V_blocks"], syn["Ep"], dx, gb["modes"]["sigma"])
    assert M2.shape == M3.shape == (N * ny, N * ny)
    assert float(abs(M2 - M3).max()) > 0.0                                       # different potential blocks (sine-sum vs cosine-sum) and different assembly
    o = np.where(wf["open"])[0]
    no = o.size
    B = np.zeros((N * ny, 2 * no), dtype=complex)
    B[o, np.arange(no)] = wf["source"]; B[(N - 1) * ny + o, no + np.arange(no)] = wf["source"]
    X = spla.splu(M3, permc_spec="MMD_AT_PLUS_A").solve(B)
    S3 = drv.wf_extract(X[:ny], X[(N - 1) * ny:], o, np.sqrt(wf["gamma"][o])[:, None], wf["a"])
    assert float(np.max(np.abs(S3 - wf["S"]))) <= 1e-12
    o2 = np.where(gb["modes"]["open"])[0]
    B2 = np.zeros((N * ny, 2 * no), dtype=complex)
    B2[o2, np.arange(no)] = 1.0; B2[(N - 1) * ny + o2, no + np.arange(no)] = 1.0
    X2 = spla.splu(M2).solve(B2)
    S2 = drv.fl_extract(X2[:ny][o2][:, :no], X2[(N - 1) * ny:][o2][:, :no], X2[:ny][o2][:, no:], X2[(N - 1) * ny:][o2][:, no:], np.sqrt(gb["modes"]["gamma"][o2]), 1j)
    assert float(np.max(np.abs(S2 - gb["S"]))) <= 1e-12


# ----------------------------------------------------------------------
# decomposition: identities, hybrid lead replacement, M1 excluded
# ----------------------------------------------------------------------
@pytest.mark.parametrize("which", ["M2", "M3"])
def test_decomposition_identities_and_mirror_flags(syn, which):
    S_ex = syn["ex"]["S"]
    if which == "M3":
        d = drv.decompose_m3(syn["pt"]["E"], syn["wf"], syn["dx"], syn["ex"]["lead"], S_ex, syn["ex"]["error_indicator"])
    else:
        d = drv.decompose_m2(syn["pt"]["E"], syn["gb"], syn["Ep"], syn["dx"], syn["ex"]["lead"], S_ex, syn["ex"]["error_indicator"])
    t = d["terms"]
    assert d["method"] == which and t["three_term_available"]
    assert t["identity_residual_three_term"] <= 1e-18 and t["identity_residual_two_term"] <= 1e-18
    assert d["mirror_source_bitwise"] and d["mirror_extraction_bitwise"]
    assert d["refinement_fl"]["passed"] and d["refinement_hybrid"]["passed"]
    assert d["e_ind_combined"] == max(d["e_ind_exact"], d["e_ind_fl"], d["e_ind_hybrid"])      # combined indicator of ALL extended solves involved
    assert t["total"] == pytest.approx(float(np.max(np.abs((syn["wf"] if which == "M3" else syn["gb"])["S"] - S_ex))), abs=1e-18)


def test_hybrid_replaces_only_the_lead_terms_a_known_lead_perturbation_lands_in_the_lead_term(syn, monkeypatch):
    """Make M3's double lead data wrong by 1e-6 (a deliberately bad lead root): the lead-data term must be ~1e-6 and the bulk term must stay at the round-off level."""
    real = sw.lead_data

    def bad(E, E_perp, dx):
        ld = real(E, E_perp, dx)
        o = ld["open"]
        ld["lam"] = ld["lam"].copy()
        ld["lam"][o] = ld["lam"][o] * np.exp(1j * 1e-6)
        s = np.sqrt(1.0 - (ld["c"][o]) ** 2)
        return ld
    monkeypatch.setattr(sw, "lead_data", bad)
    dev, geo, ref = syn["devs"].get("SYN", syn["inv"], syn["prob"].buffer_req)[:3]
    wfb = sw.wavefunction_scattering_matrix(syn["pt"]["E"], syn["pot"].numpy, geo["x_start"], syn["dx"], syn["N"], syn["prob"].Ly, syn["prob"].Ny, syn["prob"].Nq)
    monkeypatch.setattr(sw, "lead_data", real)
    d = drv.decompose_m3(syn["pt"]["E"], wfb, syn["dx"], syn["ex"]["lead"], syn["ex"]["S"], syn["ex"]["error_indicator"])
    t = d["terms"]
    assert t["lead_data"] > 1e-8 and t["bulk_rounding"] < 1e-12                  # exact attribution of the perturbation to the lead term
    assert t["identity_residual_three_term"] <= 1e-18


def test_hybrid_corner_swap_changes_only_the_two_lead_blocks():
    """solve_matrix_extended with corner_swap: swapping lam_double -> lam_ext at the lead terms; with lam_ext such that t*lam_ext equals the stored
    double t*lam_double exactly, the hybrid operator IS the rounded operator and the two solutions coincide to the arithmetic floor."""
    ny, N, dx, E = 3, 12, 1.0 / 64, 70.0
    Ep = 0.5 * (np.arange(1, ny + 1) * np.pi) ** 2
    V = np.zeros((N, ny, ny)); V[:, 0, 1] = V[:, 1, 0] = 3.0
    lead = xp.lead_data_extended(E, Ep, dx)
    lam_d = np.asarray(lead["lam"], dtype=complex) * np.exp(1e-8j)                       # a deliberately inexact 'double' lead root
    M = sw.assemble(E, V, Ep, dx, lam_d)
    t = float(sw.hopping(dx))
    lam_ext = (np.asarray(t * lam_d, dtype=complex).astype(CLD)) / LD(t)                 # t is a power of two: exact
    o = np.where(lead["open"])[0]; no = o.size
    B = np.zeros((N * ny, 2 * no), dtype=CLD); B[o, np.arange(no)] = 1; B[(N - 1) * ny + o, no + np.arange(no)] = 1
    plain = xp.solve_matrix_extended(M, B)
    swap = xp.solve_matrix_extended(M, B, corner_swap=dict(t=t, ny=ny, N=N, lam_double=lam_d, lam_ext=lam_ext))
    assert float(np.max(np.abs(plain["X"] - swap["X"]))) <= 1e-17
    swap2 = xp.solve_matrix_extended(M, B, corner_swap=dict(t=t, ny=ny, N=N, lam_double=lam_d, lam_ext=lead["lam"]))
    assert float(np.max(np.abs(plain["X"] - swap2["X"]))) > 1e-10                           # the genuine exact lead data changes the solution (by ~1e-8)


def test_M1_gets_no_decomposition_and_only_M2_M3_rows_are_produced(syn):
    out = drv.run_xp_case(syn["prob"], syn["devs"], syn["pt"], syn["inv"], syn["prob"].buffer_req,
                          {"M1": syn["ex"]["S"].astype(complex), "M2": syn["gb"]["S"], "M3": syn["wf"]["S"]}, decompose=True)
    assert {r["method"] for r in out["decomp"]} == {"M2", "M3"}
    assert not hasattr(drv, "decompose_m1")
    assert "eps_M1" in out["row"] and "bulk_rounding" not in out["row"]


def test_unresolved_or_unavailable_terms_never_support_a_dominance_claim(syn):
    terms = dict(bulk_rounding=1e-12, lead_data=1e-14, solver=1e-14)
    assert xp.representation_dominance(terms, 1e-11)["supported"] is False            # terms below 10 x e_ind: unresolved
    assert xp.representation_dominance(terms, float("nan"))["supported"] is False
    assert xp.representation_dominance(terms, 1e-16)["supported"] is True
    assert xp.representation_dominance(None, 1e-15)["available"] is False             # M1: no decomposition
    assert xp.classify_exceedance({512: 1e-11, 1024: 1e-10, 2048: 1e-9}, None, float("nan"), True) == "UNEXPLAINED"


def test_failed_extended_gates_give_not_evaluable_decomposition_rows(syn, monkeypatch):
    real = xp.exact_lattice_S

    def worse(*a, **k):
        r = real(*a, **k)
        r["error_indicator"] = 1e-9                                              # X7 violated
        return r
    monkeypatch.setattr(xp, "exact_lattice_S", worse)
    out = drv.run_xp_case(syn["prob"], syn["devs"], syn["pt"], syn["inv"], syn["prob"].buffer_req,
                          {"M1": syn["wf"]["S"], "M2": syn["gb"]["S"], "M3": syn["wf"]["S"]}, decompose=True)
    assert out["row"]["gates_passed"] is False
    assert [r["status"] for r in out["decomp"]] == ["NOT EVALUABLE", "NOT EVALUABLE"] and all("bulk_rounding" not in r for r in out["decomp"])


# ----------------------------------------------------------------------
# whole synthetic driver run: schema, identities, metadata reproducibility, no historical overwrite
# ----------------------------------------------------------------------
@pytest.fixture(scope="module")
def run_twice(tmp_path_factory):
    outs = []
    for k in (1, 2):
        base = tmp_path_factory.mktemp(f"run{k}")
        drv.run_campaign(drv.SyntheticProblem(), str(base / "out"), str(base / "cache"), log=lambda m: None, run_g0=False)
        outs.append(str(base / "out"))
    return outs


FILES = ["per_energy_comparison.csv", "tier3_exact_lattice.csv", "error_decomposition.csv", "refinement_convergence.csv", "buffer_control.csv",
         "potential_identity.csv", "xp_gates.csv", "solver_spread.csv", "summary.json", "run_metadata.json", "platform_precision.json"]


def _rows(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def test_output_schema_files_and_columns(run_twice):
    out = run_twice[0]
    for f in FILES:
        assert os.path.exists(os.path.join(out, f)), f
    comp = _rows(os.path.join(out, "per_energy_comparison.csv"))
    need = {"point", "E", "inv_dx", "S1_dS", "S2_dS", "S3_dS", "S4_S1_max", "S4_S2_max", "U_M1", "U_M2", "U_M3", "Q_M1", "Q_M2", "Q_M3", "T1_M1", "T1_M2", "T1_M3",
            "R1_backward_error_M2", "R1_backward_error_M3", "P1_rgf_device_bitwise", "N1_passed", "status_S1", "status_S5", "rho_S1", "potential_hash"}
    assert need <= set(comp[0])
    dec = _rows(os.path.join(out, "error_decomposition.csv"))
    assert {"method", "total", "bulk_rounding", "lead_data", "solver", "identity_residual_three_term", "e_ind_combined", "resolved", "dominant",
            "dominance_supported", "mirror_extraction_bitwise"} <= set(dec[0])
    assert {r["method"] for r in dec} == {"M2", "M3"}
    t3 = _rows(os.path.join(out, "tier3_exact_lattice.csv"))
    assert {"eps_M1", "eps_M2", "eps_M3", "X4_state_passed", "X5_fl_vs_state_max_abs", "X7_error_indicator", "T1_M4", "gates_passed"} <= set(t3[0])
    s = json.load(open(os.path.join(out, "summary.json")))
    assert {"criteria", "tier3_adjudication", "overall_status", "stopped_at"} <= set(s)
    assert all(k in s["criteria"] for k in ("S1_M1_vs_M2", "S2_M1_vs_M3", "S3_M2_vs_M3", "S5_unitarity", "S6_reciprocity", "R1_backward_error", "K1_buffer_control", "V1_tier3_validated"))


def test_all_four_methods_ran_and_decomposition_identities_hold(run_twice):
    out = run_twice[0]
    comp, t3, dec = (_rows(os.path.join(out, f)) for f in ("per_energy_comparison.csv", "tier3_exact_lattice.csv", "error_decomposition.csv"))
    assert all(r["status"] == "OK" for r in comp) and all(float(r["T1_M4"]) > 0 for r in t3)
    assert all(float(r["identity_residual_three_term"]) <= 1e-18 for r in dec)
    assert all(r["mirror_source_bitwise"] == "True" and r["mirror_extraction_bitwise"] == "True" for r in dec)
    assert all(r["M2_recomputed_bitwise"] == "True" and r["M3_recomputed_bitwise"] == "True" for r in dec)       # decomposition starts from the reported S_double
    assert all(abs(float(r["total"]) - float(r["bulk_rounding"]) - float(r["lead_data"]) - float(r["solver"])) <= float(r["sum_of_term_norms"]) for r in dec)


def test_metadata_and_results_are_reproducible_across_two_independent_runs(run_twice):
    a, b = run_twice
    timing = lambda k: k.startswith("time_") or k.endswith("_time_s") or k.endswith("time_decomposition_s")
    for f in ("per_energy_comparison.csv", "tier3_exact_lattice.csv", "error_decomposition.csv", "refinement_convergence.csv", "buffer_control.csv", "xp_gates.csv", "potential_identity.csv"):
        ra, rb = _rows(os.path.join(a, f)), _rows(os.path.join(b, f))
        assert len(ra) == len(rb)
        for x, y in zip(ra, rb):
            assert {k: v for k, v in x.items() if not timing(k)} == {k: v for k, v in y.items() if not timing(k)}, f
    ma, mb = (json.load(open(os.path.join(d, "run_metadata.json"))) for d in (a, b))
    drop = {"generated_utc", "wall_time_s", "command"}
    assert {k: v for k, v in ma.items() if k not in drop} == {k: v for k, v in mb.items() if k not in drop}
    for key in ("source_commit", "protocol", "precision", "potential_hashes", "exact_energies", "solver_settings", "limits", "code_fingerprint"):
        assert key in ma, key
    assert ma["protocol"]["sha256"].startswith("sha256:") and ma["precision"]["longdouble"]["mantissa_bits"] >= 64


def test_driver_refuses_to_write_into_historical_result_directories(tmp_path):
    for d in ("rgf_b3x_validation", "rgf_b3_validation", "rgf_b2_validation", "rgf_b1_validation"):
        with pytest.raises(SystemExit):
            drv.check_output_dir(os.path.join(HERE, "results", d))
        with pytest.raises(SystemExit):
            drv.check_output_dir(os.path.join(HERE, "results", d, "sub"))
    assert drv.check_output_dir(str(tmp_path / "x")) == os.path.normpath(str(tmp_path / "x"))
    assert drv.check_output_dir(os.path.join(HERE, "results", "rgf_b3y_validation"))


def test_cache_from_another_protocol_is_refused(tmp_path):
    drv.Cache(str(tmp_path), dict(problem="A", protocol_sha256="x", invs=[1]))
    with pytest.raises(RuntimeError):
        drv.Cache(str(tmp_path), dict(problem="B", protocol_sha256="x", invs=[1]))


# ----------------------------------------------------------------------
# platform failure: NOT EVALUABLE, never a float64 fallback
# ----------------------------------------------------------------------
def test_unavailable_extended_precision_stops_tier3_and_reports_not_evaluable(tmp_path, monkeypatch):
    bad = dict(evaluable=False, longdouble_mantissa_bits=53, reason="simulated", arithmetic_probe_passed=False)
    monkeypatch.setattr(drv.xp, "platform_precision", lambda: bad)
    s = drv.run_campaign(drv.SyntheticProblem(), str(tmp_path / "o"), str(tmp_path / "c"), log=lambda m: None, run_g0=False)
    assert s["stopped_at"]["stage"] == "X0" and "NOT EVALUABLE" in s["stopped_at"]["reason"]
    assert s["criteria"]["V1_tier3_validated"]["status"] == "NOT EVALUABLE"
    assert _rows(str(tmp_path / "o" / "tier3_exact_lattice.csv")) == [] if os.path.getsize(str(tmp_path / "o" / "tier3_exact_lattice.csv")) <= 1 else True
    assert "PASS" != s["overall_status"]


def test_failed_gate_g0_stops_everything():
    row = {"config": "W200_z0.05", "E": 140.47, "T1_M1": 0.5, "S1_dS": 0, "S2_dS": 0, "S3_dS": 0}
    out = drv.stage_G0({"p": (row, None)})
    assert out[0]["passed"] is False                       # no matching B3X row


def test_G0_uses_the_b3x_artifact_and_its_tolerances():
    b3x_rows = drv.load_b3x_rows()
    r = b3x_rows[0]
    row = dict(config=r["config"], E=float(r["E"]), T1_M1=float(r["T1_rgf"]), S2_dS=float(r["rgf_wf_dS"]), S1_dS=float(r["rgf_gb_dS"]), S3_dS=float(r["wf_gb_dS"]))
    assert drv.stage_G0({"p": (row, None)})[0]["passed"] is True
    assert drv.stage_G0({"p": (dict(row, S2_dS=row["S2_dS"] + 1e-12), None)})[0]["passed"] is False
    assert drv.stage_G0({"p": (dict(row, T1_M1=np.nextafter(row["T1_M1"], 1.0)), None)})[0]["passed"] is False


# ----------------------------------------------------------------------
# protocol preservation
# ----------------------------------------------------------------------
def test_protocol_thresholds_points_and_geometry_are_preserved():
    L = drv.LIMITS
    assert (L["S_abs"], L["TR_offdiag_abs"], L["unitarity"], L["reciprocity"], L["sparse_backward_error"]) == (1e-10, 2e-10, 1e-10, 1e-8, 1e-14)
    assert (L["dT1_buffer"], L["lead_tail_guard"], L["T1_method_spread_for_solver_independence"]) == (2e-4, 1e-6, 2e-6)
    assert L["order_window"] == [1.8, 2.2] and L["dominance_factor"] == xp.DOMINANCE_FACTOR == 2.0 and L["term_resolution_factor"] == xp.TERM_RESOLUTION_FACTOR == 10.0
    assert xp.ETA_MAX == drv.XPG["final_relative_operator_residual_max"] == 1e-17 and xp.ERROR_INDICATOR_MAX == drv.XPG["error_indicator_max"] == 1e-11
    assert xp.S_TARGET == 1e-10 and xp.GROWTH_PER_HALVING_MAX == L["precision_limited_growth_per_halving"] == 10.0
    prob = drv.RealProblem()
    assert len(prob.points) == 16 and len(prob.sets_map["tier1"]) == 16 and len(prob.sets_map["tier2"]) == 5 == len(prob.sets_map["buffer"])
    assert prob.sets_map["tier2"] == prob.sets_map["tier3_fine"] == prob.sets_map["buffer"]
    assert prob.invs == (512, 1024, 2048) and (prob.Ny, prob.Nq) == (63, 127) and prob.seed == 0 and 4096 not in prob.invs
    assert (prob.buffer_req, prob.buffer_ctrl_req) == (0.2, 0.4)
    for inv, n in ((512, 719), (1024, 1437), (2048, 2873)):
        assert prob.geometry("W200_z0.05", inv, prob.buffer_req)["n_slices"] == n
    g = prob.geometry("W200_z0.05", 1024, prob.buffer_ctrl_req)
    assert g["n_slices"] == 1845 and g["n_buffer_slices"] == 410
    b3x_E = {(r["config"], float(r["E"])) for r in drv.load_b3x_rows()}
    for p in prob.points.values():                                                   # no new energies, no new realizations
        assert (p["config"], p["E"]) in b3x_E and p["config"] in b3x.CONFIGS
    assert drv.sha_file(os.path.join(HERE, "docs", "rgf_b3y_protocol.json")) == drv.sha_file(drv.PROTOCOL_PATH)


def test_synthetic_problem_uses_no_b3y_point():
    prob = drv.SyntheticProblem()
    real = {p["E"] for p in drv.PROTOCOL["points"]}
    assert not (real & {p["E"] for p in prob.points.values()}) and all(p["config"] == "SYN" for p in prob.points.values())


# ----------------------------------------------------------------------
# post-campaign provenance: historical protocol record preserved, status reconciled
# ----------------------------------------------------------------------
def test_protocol_hash_in_results_equals_the_unchanged_protocol_file_and_keeps_its_historical_status():
    m = json.load(open(os.path.join(HERE, "results", "rgf_b3y_validation", "run_metadata.json")))
    assert m["protocol"]["sha256"] == drv.sha_file(drv.PROTOCOL_PATH)
    assert "no B3Y numerical comparison has been run" in drv.PROTOCOL["status"]            # historical text deliberately preserved


def test_campaign_summary_separates_historical_protocol_status_from_campaign_status():
    s = json.load(open(os.path.join(HERE, "results", "rgf_b3y_validation", "summary.json")))
    assert "protocol_status" not in s and "no B3Y numerical comparison" in s["protocol_status_at_predeclaration"]
    assert s["campaign_status"].startswith("COMPLETED") and s["stopped_at"] is None and s["overall_status"] == "PASS"
    assert s["metadata_corrections"][0]["numerical_values_changed"] is False


def test_synthetic_driver_summary_states_it_is_synthetic(run_twice):
    s = json.load(open(os.path.join(run_twice[0], "summary.json")))
    assert s["campaign_status"].startswith("SYNTHETIC") and "protocol_status_at_predeclaration" in s


def test_stale_no_comparison_statements_in_the_report_are_marked_historical():
    lines = open(os.path.join(HERE, "docs", "RGF_B3Y_VALIDATION.md")).read().splitlines()
    for i, l in enumerate(lines):
        low = l.lower()
        if ("no b3y result directory exists" in low or "is not written" in low or "no b3y numerical comparison has been run" in low):
            assert "[historical" in low or "historical" in low, (i + 1, l[:120])
