"""B3X: extended independent validation of the RGF (validate_rgf_b3x.py, sparse_wavefunction.py,
docs/RGF_B3X_VALIDATION.md, docs/rgf_b3x_protocol.json).

Part 1 tests the independent wave-function reference on problems with known answers; part 2 the
driver and the predeclared protocol; part 3 the committed artifacts (skipped if absent).
Thresholds are checked for preservation only, never changed.
"""
import ast
import copy
import csv
import hashlib
import json
import math
import os

import numpy as np
import pytest

import smatrix_spec as sm
import sparse_reference as sr
import sparse_wavefunction as sw
import validate_rgf_b3 as b3
import validate_rgf_b3x as v
from stationary import Device, potential_blocks, scattering_matrix, transverse_basis

RES = "results/rgf_b3x_validation"
SUMMARY = os.path.join(RES, "summary.json")
HAVE = os.path.exists(SUMMARY) and os.path.exists(os.path.join(RES, "per_energy_comparison.csv"))
SMALL = dict(dx=1 / 128, Ny=15, Nq=31, buffer=0.2)


def _toy(n_slices=10, ny=6, dx=0.05, seed=5, vmax=40.0):
    rng = np.random.default_rng(seed)
    A = rng.uniform(-vmax, vmax, (n_slices, ny, ny))
    V = 0.5 * (A + np.transpose(A, (0, 2, 1)))
    basis = transverse_basis(1.0, ny)
    return V, basis, dx, Device(xs=dx * np.arange(n_slices), dx=dx, basis=basis, V=V, V_lead_max=0.0)


# ----------------------------------------------------------------------
# Part 1: the wave-function reference
# ----------------------------------------------------------------------
def test_wavefunction_reference_is_independent_of_rgf_and_fisher_lee():
    src = open("sparse_wavefunction.py").read()
    tree = ast.parse(src)
    imported = {n.names[0].name for n in ast.walk(tree) if isinstance(n, ast.Import)} | \
               {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    assert not imported & {"stationary", "smatrix_spec", "sparse_reference", "validate_rgf_b3", "validate_rgf_b3x"}
    code = "\n".join(l for l in src.splitlines() if not l.strip().startswith("#"))
    body = code.split('"""', 2)[2]                                   # drop the module docstring
    for forbidden in ("fisher_lee", "assemble_S", "lead_modes", "rgf_endpoint_blocks", "G11", "GN1"):
        assert forbidden not in body, forbidden


def test_retarded_root_closed_form_matches_the_complex_energy_root():
    Ep = transverse_basis(1.0, 63).E_perp
    for E in (13.3, 88.5134, 140.474, 204.93, 371.356, 649.7):
        ld = sw.lead_data(E, Ep, 1 / 512)
        eta = sw.retarded_root_eta(E, Ep, 1 / 512)
        assert np.max(np.abs(ld["lam"] - eta)) <= v.LIMITS["retarded_root_agreement"]
        assert np.all(np.abs(ld["lam"][ld["open"]].imag) > 0) and np.allclose(np.abs(ld["lam"][ld["open"]]), 1.0, atol=1e-15)
        cl = ld["lam"][~ld["open"]]
        assert np.all(cl.imag == 0) and np.all(np.abs(cl) < 1)
        assert np.all(ld["gamma"][~ld["open"]] == 0.0)
        k_dx = np.arccos(ld["c"][ld["open"]])
        assert np.allclose(ld["gamma"][ld["open"]], 2 * sw.hopping(1 / 512) * np.sin(k_dx), rtol=1e-13)
        # the retarded open root is e^{+ i k dx} (outgoing), not its conjugate
        assert np.allclose(ld["lam"][ld["open"]].imag, np.sin(k_dx), rtol=1e-13)
    with pytest.raises(ValueError):
        sw.lead_data(Ep[2], Ep, 1 / 512)                             # exactly at a threshold


def test_cosine_sum_potential_blocks_equal_the_matrix_product_quadrature():
    rng = np.random.default_rng(1)
    ny, nq = 12, 29
    V = rng.uniform(-50, 50, (4, nq))
    basis = transverse_basis(1.0, ny, nq)
    explicit = np.empty((4, ny, ny))
    for j in range(4):
        for n in range(1, ny + 1):
            for m in range(1, ny + 1):
                explicit[j, n - 1, m - 1] = np.sum(V[j] * basis.X[:, n - 1] * basis.X[:, m - 1]) * basis.weight
    cos_blocks = sw.potential_blocks_cosine(V, ny, nq)
    assert np.max(np.abs(cos_blocks - explicit)) <= 1e-13 * 50
    assert np.max(np.abs(cos_blocks - potential_blocks(V, basis))) <= 1e-13 * 50
    assert np.max(np.abs(cos_blocks - np.transpose(cos_blocks, (0, 2, 1)))) <= 1e-13 * 50
    # y-uniform potential: V_j = V(x_j) I (channels decouple)
    u = np.tile(np.array([3.0, -7.0, 11.0])[:, None], (1, nq))
    ub = sw.potential_blocks_cosine(u, ny, nq)
    assert np.max(np.abs(ub - u[:, :1, None] * np.eye(ny))) <= 1e-13 * 11


def test_clean_guide_amplitudes_and_phase():
    N, ny, dx, E = 12, 4, 0.1, 30.0
    Ep = transverse_basis(1.0, ny).E_perp
    out = sw.wavefunction_S_from_blocks(E, np.zeros((N, ny, ny)), Ep, dx)
    lam = out["lam"][out["open"]]
    assert np.allclose(np.diag(out["t_RL"]), lam ** (N - 1), atol=1e-11) and np.max(np.abs(out["r_L"])) < 1e-11
    assert np.allclose(np.diag(out["t_LR"]), lam ** (N - 1), atol=1e-11)


def test_matches_a_green_function_free_transfer_matrix_1d():
    dx, N, E = 0.05, 40, 30.0
    rng = np.random.default_rng(11)
    Vj = np.where(np.arange(N) % 9 < 4, rng.uniform(5, 25, N), 0.0)
    Ep = np.array([0.5 * np.pi ** 2])
    out = sw.wavefunction_S_from_blocks(E, Vj.reshape(N, 1, 1), Ep, dx)
    t = sw.hopping(dx)
    lam = sw.lead_data(E, Ep, dx)["lam"][0]
    psi = {N + 1: lam, N: 1.0 + 0j}
    for j in range(N, 0, -1):
        psi[j - 1] = ((2 * t + Ep[0] + Vj[j - 1] - E) * psi[j] - t * psi[j + 1]) / t
    a, b = np.linalg.solve(np.array([[1.0, 1.0], [1 / lam, lam]]), np.array([psi[1], psi[0]]))
    assert abs(out["t_RL"][0, 0] - 1 / a) < 1e-9 and abs(out["r_L"][0, 0] - b / a) < 1e-9
    assert 0.05 < abs(1 / a) ** 2 < 0.95


def test_random_multichannel_agrees_with_rgf_and_greenblock_reference_and_is_unitary_reciprocal():
    V, basis, dx, dev = _toy()
    E = 60.0                                                         # 3 open + 3 closed channels
    w = sw.wavefunction_S_from_blocks(E, V, basis.E_perp, dx)
    r = scattering_matrix(E, dev)
    g = sr.sparse_S_from_blocks(E, V, basis.E_perp, dx)
    assert w["n_open"] == r["n_open"] == g["n_open"] == 3
    assert np.max(b3.s_difference(w["S"], r["S"])) < 1e-12 and np.max(b3.s_difference(w["S"], g["S"])) < 1e-12
    assert b3.unitarity(w["S"])["U"] < 1e-13 and b3.reciprocity(w["S"], 3)["Q"] < 1e-13
    assert w["residual_after_refinement"] < 1e-14
    assert w["closed_to_open_weight_slice1"] > 0.1                   # the evanescent content is real
    cc = v.closed_channel_checks(E, r, w, basis.E_perp, dx=dx)      # closed-channel amplitudes from the RGF Green blocks
    assert cc["closed_amp_diff_rgf_vs_wf"] <= 1e-10 * cc["closed_amp_scale"] and cc["closed_amp_scale"] > 0
    assert cc["gamma_closed_max_rgf"] == cc["gamma_closed_max_wf"] == 0.0 and cc["open_masks_equal"]


def test_wrong_source_sign_is_caught_by_unitarity():
    V, basis, dx, _ = _toy()
    good = sw.wavefunction_S_from_blocks(60.0, V, basis.E_perp, dx)
    bad = sw.wavefunction_S_from_blocks(60.0, V, basis.E_perp, dx, _source_sign=-1.0)
    assert b3.unitarity(good["S"])["U"] < 1e-13 and b3.unitarity(bad["S"])["U"] > 1e-3
    assert np.max(b3.s_difference(bad["S"], scattering_matrix(60.0, _toy()[3])["S"])) > 1e-3


# ----------------------------------------------------------------------
# Part 2: protocol and driver
# ----------------------------------------------------------------------
def test_predeclared_protocol_preserved():
    with open("docs/rgf_b3x_protocol.json") as f:
        pr = json.load(f)
    assert pr["limits"] == dict(S_agreement=1e-10, T_R_agreement=2e-10, offdiag_transmission_agreement=2e-10, unitarity=1e-10,
                                reciprocity=1e-8, dT1_refinement=2e-4, potential_block_agreement_relative=1e-12, retarded_root_agreement=1e-9)
    assert v.LIMITS == pr["limits"]
    assert len(pr["criteria"]) == 14 and set(pr["limit_justification"]) == set(pr["limits"])
    assert sorted(v.CONFIGS) == sorted(["W100_z0.02", "W100_z0.05", "W100_z0.2", "W200_z0.02", "W200_z0.05", "W200_z0.2"])
    assert {(c["W"], c["zeta"]) for c in v.CONFIGS.values()} == {(W, z) for W in (100.0, 200.0) for z in (0.02, 0.05, 0.2)}
    assert v.SEED == 0 and v.BASE == dict(dx=1 / 512, Ny=63, Nq=127, buffer=0.2)
    assert v.VARIANTS["Ny127_Nq254"]["Nq"] == v.VARIANTS["Nq"]["Nq"] == 254 and v.VARIANTS["dx"]["dx"] == 1 / 1024
    assert v.VARIANTS["buffer"]["buffer"] == 0.4 and v.LADDER == (1 / 256, 1 / 512, 1 / 1024, 1 / 2048)
    # B1 / B2 / B3 thresholds untouched
    import validate_rgf_b1 as b1
    import validate_rgf_b2 as b2
    assert b1.LIMITS == dict(T_diag=1e-12, mixing=1e-12, R=1e-12, U=1e-10, phase=1e-10, Q=1e-8)
    assert b2.LIMITS == dict(dT=2e-4, mixing=1e-12, U=1e-10, Q=1e-8)
    assert b3.LIMITS == dict(S_agreement=1e-10, U=1e-10, Q=1e-8, dT1=2e-4)


@pytest.mark.parametrize("cid", sorted(v.CONFIGS))
def test_configuration_seed_and_energies(cid):
    cfg, pot = v.config_potential(cid)
    p = v.CONFIGS[cid]
    assert cfg.disorder_strength == p["W"] and cfg.disorder_corr_length == p["zeta"] and tuple(cfg.disorder_region) == (1.9, 2.9)
    assert pot.disorder.seed == 0 and pot.disorder.n_modes == 512
    E = b3.case_energies(cfg)
    assert len(E) == 5 and [e["n_open_exact"] for e in E] == [4, 5, 6, 7, 8]
    assert all(e["rel_dist_to_threshold"] >= 1e-3 for e in E)
    assert np.allclose([e["E"] for e in E], [0.5 * ((20 + j / (math.sqrt(2) * 0.2)) ** 2 + math.pi ** 2) for j in (-2, -1, 0, 1, 2)], rtol=1e-14)


def test_same_realization_and_buffer_geometry():
    for cid in ("W200_z0.02", "W100_z0.2"):
        cfg, pot = v.config_potential(cid)
        dev, geo, ref, info = v.build_rgf(cid, cfg, pot, SMALL)
        assert info["bitwise_equal"] and info["max_abs_diff"] == 0.0 and info["shape"] == [dev.n_slices, SMALL["Nq"]]
        out = sw.wavefunction_scattering_matrix(b3.case_energies(cfg)[1]["E"], pot.numpy, geo["x_start"], SMALL["dx"], dev.n_slices,
                                                1.0, SMALL["Ny"], SMALL["Nq"])
        assert np.array_equal(out["V_samples"], ref) and np.array_equal(out["xs"], dev.xs)
        assert np.max(np.abs(out["V_blocks"] - dev.V)) <= 1e-12 * np.max(np.abs(ref))      # C2 on a live device
    cfgA = v.config_cfg("W100_z0.05")
    g1, g2, g3 = (b3.geometry(cfgA, 1 / 512, 0.2), b3.geometry(cfgA, 1 / 1024, 0.2), b3.geometry(cfgA, 1 / 512, 0.4))
    assert g1["buffer_left"] == g2["buffer_left"] == 103 / 512 and g3["buffer_left"] == 205 / 512 and cfgA.disorder_region == (1.9, 2.9)
    # W = 200 is exactly twice W = 100 for the same seed and zeta (same draw); different zeta are different potentials
    _, pa = v.config_potential("W100_z0.05"); _, pb = v.config_potential("W200_z0.05")
    xs = np.linspace(1.5, 3.3, 7)[:, None] * np.ones((1, 5)); ys = np.ones((7, 1)) * np.linspace(0.1, 0.9, 5)[None, :]
    assert np.array_equal(pb.numpy(xs, ys), 2 * pa.numpy(xs, ys))


def test_pipeline_on_a_small_device_all_three_methods_agree():
    cfg, pot = v.config_potential("W200_z0.02")
    dev, geo, ref, info = v.build_rgf(("W200_z0.02"), cfg, pot, SMALL)
    E = b3.case_energies(cfg)[2]["E"]
    rg = scattering_matrix(E, dev)
    wf = sw.wavefunction_scattering_matrix(E, pot.numpy, geo["x_start"], SMALL["dx"], dev.n_slices, 1.0, SMALL["Ny"], SMALL["Nq"])
    gb = sr.sparse_scattering_matrix(E, pot.numpy, geo["x_start"], SMALL["dx"], dev.n_slices, 1.0, SMALL["Ny"], SMALL["Nq"])
    c = v.compare_S(rg["S"], wf["S"], rg["n_open"])
    assert c["dS"] <= v.LIMITS["S_agreement"] and c["dP_all"] <= v.LIMITS["T_R_agreement"] and c["dP_offdiag_trans"] <= v.LIMITS["offdiag_transmission_agreement"]
    assert np.max(b3.s_difference(rg["S"], gb["S"])) <= v.LIMITS["S_agreement"]
    mm = v.method_metrics(wf["S"], rg["n_open"])
    assert mm["U"] <= v.LIMITS["unitarity"] and mm["Q"] <= v.LIMITS["reciprocity"]
    cc = v.closed_channel_checks(E, rg, wf, dev.basis.E_perp, dx=SMALL["dx"])
    assert cc["n_closed"] == SMALL["Ny"] - rg["n_open"] >= 1 and cc["gamma_closed_max_rgf"] == 0.0


def test_compare_S_sees_phase_and_reflection_and_offdiagonal_errors():
    n = 3
    S = np.eye(2 * n, dtype=complex)
    S2 = S.copy(); S2[0, 0] = np.exp(1j * 3e-10)                      # phase-only error in r_L
    c = v.compare_S(S, S2, n)
    assert c["dS"] > 1e-10 and c["dP_all"] < 1e-15 and c["dP_refl"] < 1e-15            # probabilities blind, S not
    S3 = S.copy(); S3[n + 1, 0] = 1e-4                                  # off-diagonal transmission error (t_RL, m != n)
    c3 = v.compare_S(S, S3, n)
    assert c3["dP_offdiag_trans"] == pytest.approx(1e-8, rel=1e-6) and c3["dP_trans"] == pytest.approx(1e-8, rel=1e-6)
    S4 = S.copy(); S4[1, 1] = 1.0 - 1e-5                               # reflection-block diagonal error
    assert v.compare_S(S, S4, n)["dP_refl"] > 1e-6


# ----------------------------------------------------------------------
# Part 3: committed artifacts
# ----------------------------------------------------------------------
def _csv(name):
    with open(os.path.join(RES, name + ".csv")) as f:
        return [{c: b3._conv(x) for c, x in row.items()} for row in csv.DictReader(f)]


def _load_R():
    return dict(energy=_csv("per_energy_comparison"), identity=_csv("potential_identity"), refine=_csv("refinement_convergence"),
                ladder=_csv("dx_ladder_diagnostic"), geometry={})


@pytest.mark.skipif(not HAVE, reason="B3X artifacts not generated")
def test_artifacts_complete():
    for n in ("per_energy_comparison per_configuration_summary refinement_convergence dx_ladder_diagnostic potential_identity").split():
        assert os.path.exists(os.path.join(RES, n + ".csv")), n
    assert os.path.exists(os.path.join(RES, "run_metadata.json"))
    en = _csv("per_energy_comparison")
    assert len(en) == 30 and {r["config"] for r in en} == set(v.CONFIGS)
    for c in v.CONFIGS:
        assert len([r for r in en if r["config"] == c]) == 5
    assert len(_csv("per_configuration_summary")) == 6 and len(_csv("dx_ladder_diagnostic")) == 30
    ref = [r for r in _csv("refinement_convergence") if r["role"] == "criterion"]
    assert len(ref) == 120 and {r["name"] for r in ref} == {"dx", "Nq", "buffer", "Ny_fixed_Nq254"}
    for r in en:
        assert r["S_dim"] == 2 * r["n_open_rgf"] and r["matrix_order"] == r["n_slices"] * 63 and r["N_q"] == 127 and r["dx"] == 1 / 512


@pytest.mark.skipif(not HAVE, reason="B3X artifacts not generated")
def test_every_criterion_is_evaluated_and_the_summary_is_the_summarize_of_the_csv():
    with open(SUMMARY) as f:
        s = json.load(f)
    R = _load_R()
    s2 = v.summarize(copy.deepcopy(R))
    assert len(s["criteria_passed"]) == 14 and s["all_criteria_evaluated"] is True
    assert s2["criteria_passed"] == s["criteria_passed"] and s2["overall_B3X_status"] == s["overall_B3X_status"]
    for k in s["criteria_passed"]:
        if "value" in s[k]:
            assert s[k]["value"] == s2[k]["value"]
            if "per_case_max" in s[k]:
                assert s[k]["per_case_max"] == s2[k]["per_case_max"] and set(s[k]["per_case_max"]) == set(v.CONFIGS)
            assert s[k]["passed"] == (s[k]["value"] <= s[k]["limit"])
    assert s["failed_criteria"] == [k for k, x in s["criteria_passed"].items() if not x] and s["b3x_full_pass"] == (not s["failed_criteria"])
    assert s["overall_B3X_status"] == ("PASS" if s["b3x_full_pass"] else "FAIL: " + ", ".join(s["failed_criteria"]))
    assert s["limits"] == v.LIMITS
    # recompute headline numbers from the raw per-energy rows
    en = R["energy"]
    assert s["C5_S_agreement_RGF_vs_wavefunction"]["value"] == max(r["rgf_wf_dS"] for r in en)
    assert s["C9_unitarity_all_methods"]["value"] == max(max(r["U_rgf"], r["U_wf"], r["U_gb"]) for r in en)
    for r in en:
        assert r["U_rgf"] == max(r["U_StS_rgf"], r["U_SSt_rgf"])


@pytest.mark.skipif(not HAVE, reason="B3X artifacts not generated")
def test_classification_responds_to_mutated_values_exactly_at_thresholds():
    R0 = _load_R()
    base = v.summarize(copy.deepcopy(R0))
    # (a) agreement exactly at the limit passes, a hair above fails (C5), and the overall status then contains C5
    R = copy.deepcopy(R0)
    for r in R["energy"]:
        r["rgf_wf_dS"] = 0.0
    R["energy"][0]["rgf_wf_dS"] = 1e-10
    assert v.summarize(copy.deepcopy(R))["C5_S_agreement_RGF_vs_wavefunction"]["passed"]
    R["energy"][0]["rgf_wf_dS"] = 1e-10 * (1 + 1e-9)
    s = v.summarize(copy.deepcopy(R))
    assert not s["C5_S_agreement_RGF_vs_wavefunction"]["passed"] and "C5_S_agreement_RGF_vs_wavefunction" in s["failed_criteria"] and not s["b3x_full_pass"]
    # (b) open-channel count mismatch -> C3; a nonzero closed-channel Gamma -> C4; bad unitarity -> C9; bad reciprocity -> C10
    for mut, crit in ((lambda r: r.update(n_open_wf=r["n_open_wf"] + 1), "C3_open_channel_count"),
                      (lambda r: r.update(gamma_closed_max_rgf=1e-30), "C4_closed_channels"),
                      (lambda r: r.update(U_wf=2e-10), "C9_unitarity_all_methods"),
                      (lambda r: r.update(Q_gb=2e-8), "C10_reciprocity_all_methods"),
                      (lambda r: r.update(rgf_gb_dP_offdiag_trans=3e-10), "C8_offdiag_transmission_agreement"),
                      (lambda r: r.update(rgf_wf_dR1=3e-10), "C7_T_R_agreement")):
        R = copy.deepcopy(R0)
        mut(R["energy"][3])
        s = v.summarize(copy.deepcopy(R))
        assert not s[crit]["passed"] and crit in s["failed_criteria"], crit
    # (c) refinement rows: exactly 2e-4 passes, above fails; an artificial all-pass R gives PASS (logic only)
    R = copy.deepcopy(R0)
    for r in R["refine"]:
        if r["role"] == "criterion":
            r["abs_delta_T1"], r["passed"] = 1e-9, True
    assert v.summarize(copy.deepcopy(R))["C11_dx_refinement"]["passed"]
    row = [r for r in R["refine"] if r["name"] == "dx"][0]
    row["abs_delta_T1"], row["passed"] = 2e-4, True
    assert v.summarize(copy.deepcopy(R))["C11_dx_refinement"]["passed"]
    row["abs_delta_T1"], row["passed"] = 2e-4 * (1 + 1e-9), False
    assert not v.summarize(copy.deepcopy(R))["C11_dx_refinement"]["passed"]
    # the real data are not altered by any of this
    assert v.summarize(copy.deepcopy(R0))["criteria_passed"] == base["criteria_passed"]


@pytest.mark.skipif(not HAVE, reason="B3X artifacts not generated")
def test_metadata_protocol_hash_and_provenance():
    with open(os.path.join(RES, "run_metadata.json")) as f:
        m = json.load(f)
    for k in ("generated_utc", "command", "source_commit", "working_tree_dirty", "python", "numpy", "scipy", "platform", "machine", "random_seed",
              "configurations", "energies", "base_discretization", "variants", "methods", "limits", "code_fingerprint", "code_files", "protocol", "geometry"):
        assert k in m, k
    assert m["random_seed"] == 0 and sorted(m["configurations"]) == sorted(v.CONFIGS) and len(m["energies"]) == 30
    assert m["protocol"]["sha256"] == v.sha_file(v.PROTOCOL_PATH) and m["protocol"]["commit"]
    from provenance import code_fingerprint
    fp = code_fingerprint("rgf_b3x")
    assert set(m["code_files"]) == set(fp["code_files"])
    for f in ("stationary.py", "smatrix_spec.py", "sparse_reference.py", "sparse_wavefunction.py", "potential.py", "disorder.py", "config.py"):
        assert m["code_files"][f] == fp["code_files"][f], f"{f} changed after the artifacts were generated"
    if "post_run_tests" in m:
        assert all({"exit_status", "command", "summary"} <= set(r) for r in m["post_run_tests"]["runs"])


@pytest.mark.skipif(not HAVE, reason="B3X artifacts not generated")
def test_per_configuration_summary_and_cross_check_consistent():
    pcs = _csv("per_configuration_summary")
    en = _csv("per_energy_comparison")
    for row in pcs:
        rr = [r for r in en if r["config"] == row["config"]]
        assert row["max_dS_rgf_wf"] == max(r["rgf_wf_dS"] for r in rr) and row["n_energies"] == 5
        assert row["pass_S_rgf_wf"] == (row["max_dS_rgf_wf"] <= 1e-10) and row["pass_dx"] == (row["max_abs_dT1_dx"] <= 2e-4)
        assert row["pass_Ny"] == (row["max_abs_dT1_Ny_fixed_Nq254"] <= 2e-4)
    if os.path.exists(os.path.join(RES, "cross_check_b3.csv")) and os.path.exists("results/rgf_b3_validation/rgf_sparse_comparison.csv"):
        xc = _csv("cross_check_b3")
        assert len(xc) == 20 and max(r["abs_diff_rgf"] for r in xc) <= 1e-12          # same RGF, same energies as the existing B3


@pytest.mark.skipif(not HAVE, reason="B3X artifacts not generated")
def test_one_row_reproduces_live():
    r0 = [r for r in _csv("per_energy_comparison") if r["config"] == "W200_z0.02"][2]
    cfg, pot = v.config_potential("W200_z0.02")
    dev, geo, _, _ = v.build_rgf("W200_z0.02", cfg, pot, v.BASE)
    E = r0["E"]
    rg = scattering_matrix(E, dev)
    wf = sw.wavefunction_scattering_matrix(E, pot.numpy, geo["x_start"], v.BASE["dx"], dev.n_slices, 1.0, 63, 127)
    assert abs(float(np.max(b3.s_difference(rg["S"], wf["S"]))) - r0["rgf_wf_dS"]) < 1e-13
    assert abs(b3.T1(rg["S"], rg["n_open"]) - r0["T1_rgf"]) < 1e-13
