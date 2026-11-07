"""B3Y Tier 3: the exact-lattice reference (sparse_exact_lattice.py), SYNTHETIC problems only.

Nothing here computes a B3Y point (no disorder realization, no B3X failed point, no B3Y comparison, no benchmark).
Synthetic energies are deliberately different from every B3Y energy; the clean guide, the y-uniform smooth barrier,
the single-channel chain and the small dense-mixing toys have exact (mpmath) or independently known answers.

Gates exercised: X0 platform, X1 lead data vs mpmath, X2 clean guide, X3 smooth barrier vs a single-channel transfer
matrix, X4 refinement contraction and operator residual, X5 scattering-state vs Fisher-Lee, X6 preconditioner-ordering
stability, X7 convergence/error indicator. The B3Y-energy instances of the gates and X8 belong to the later,
separately authorized numerical stage.
"""
import ast
import hashlib
import json
import math
import os
import types

import mpmath as mp
import numpy as np
import pytest
import scipy.sparse as sp

import sparse_exact_lattice as xp

PROTOCOL = "docs/rgf_b3y_protocol.json"
with open(PROTOCOL) as _f:
    PR = json.load(_f)

E_SYN = {"A": 95.0, "B": 150.25, "C": 263.5}          # 4, 5, 7 open channels; none is a B3Y energy


def eperp(ny):
    return 0.5 * (np.arange(1, ny + 1) * np.pi / 1.0) ** 2          # identical expression to the transverse basis of the solvers


def geometry(inv):
    nb = 103 * inv // 512
    return 1.0 / inv, 2 * nb + inv + 1, 1.9 - nb / inv


def barrier(xs):                                                    # the B2 smooth barrier (V_b = 150, width 0.1, edge 0.01, centre 2.4)
    return 150.0 * 0.5 * (np.tanh((xs - 2.35) / 0.01) - np.tanh((xs - 2.45) / 0.01))


def to_mp(z):
    with mp.workdps(50):
        return xp._to_mp(z)


# ----------------------------------------------------------------------
# independence, constants, protocol preservation
# ----------------------------------------------------------------------
def test_module_imports_nothing_from_the_existing_solvers():
    tree = ast.parse(open("sparse_exact_lattice.py").read())
    imported = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            imported |= {a.name.split(".")[0] for a in n.names}
        elif isinstance(n, ast.ImportFrom):
            imported.add((n.module or "").split(".")[0])
    assert imported <= {"math", "platform", "sys", "time", "mpmath", "numpy", "scipy", "__future__"}, imported
    names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)} | {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
    for forbidden in ("stationary", "smatrix_spec", "sparse_reference", "sparse_wavefunction", "rgf_endpoint_blocks", "fisher_lee",
                      "lead_modes", "scattering_matrix", "wavefunction_scattering_matrix", "sparse_scattering_matrix"):
        assert forbidden not in names, forbidden


def test_constants_equal_the_predeclared_protocol_values():
    g, lim = PR["xp_gates"], PR["limits"]
    assert xp.MIN_MANTISSA_BITS == g["platform_min_mantissa_bits"] == 64
    assert xp.LEAD_REL_TOL == g["lead_data_vs_mpmath_rel"] and xp.SCALAR_ABS_TOL == g["scalar_problem_vs_mpmath_abs"]
    assert xp.CONTRACTION_MAX == g["refinement_contraction_max"] and xp.MAX_ITERATIONS == g["refinement_max_iterations"]
    assert xp.ETA_MAX == g["final_relative_operator_residual_max"] and xp.FL_WF_TOL == g["fisher_lee_vs_wavefunction_forms_abs"]
    assert xp.ORDERING_TOL == g["preconditioner_ordering_stability_abs"] and xp.ERROR_INDICATOR_MAX == g["error_indicator_max"] == 1e-11
    assert xp.VALIDATION_ABS == g["validation_vs_B3X_refs_at_dx_1_512_abs"] and xp.MPMATH_DPS == g["mpmath_dps"]
    assert xp.GROWTH_PER_HALVING_MAX == lim["precision_limited_growth_per_halving"] == 10.0
    assert xp.DOMINANCE_FACTOR == lim["dominance_factor"] == 2.0 and xp.TERM_RESOLUTION_FACTOR == lim["term_resolution_factor"] == 10.0
    assert xp.S_TARGET == lim["S_abs"] == 1e-10                      # the absolute S criterion is preserved


def test_protocol_revision_1_clarifications_present_and_numbers_preserved():
    assert "uncertainty_max" not in PR["xp_gates"] and "error_indicator_note" in PR["xp_gates"]
    assert "rigorous uncertainty bound" in PR["xp_gates"]["error_indicator_note"]
    assert "64 mantissa bits" in PR["xp_gates"]["platform_requirement"] and "NOT EVALUABLE" in PR["xp_gates"]["platform_requirement"]
    assert "never substituted" in PR["xp_gates"]["platform_requirement"] and "Windows" in PR["xp_gates"]["platform_requirement"]
    assert [r["revision"] for r in PR["revision_history"]] == [0, 1]
    assert PR["limits"]["S_abs"] == 1e-10 and "DIAGNOSTIC ONLY" in PR["rho_diagnostic"]["status"]
    for k in ("PRECISION-CONSISTENT", "UNEXPLAINED"):
        assert k in PR["status_vocabulary"]
    assert "operational" in PR["limits"]["precision_limited_growth_per_halving_note"] and "NOT a mathematical proof" in PR["limits"]["precision_limited_growth_per_halving_note"]
    assert "NO mechanism claim" in PR["tier3_adjudication"]["step3"] and "UNEXPLAINED" in PR["tier3_adjudication"]["step5"]
    doc = open("docs/RGF_B3Y_VALIDATION.md").read()
    assert "e_ind" in doc and "NOT EVALUABLE" in doc and "64 mantissa bits" in doc and "not a mathematical proof" in doc
    assert doc.count("u_XP") == 1                                      # only the sentence announcing the rename
    # earlier numbers untouched (B3X limits)
    b = json.load(open("docs/rgf_b3x_protocol.json"))["limits"]
    assert (b["S_agreement"], b["T_R_agreement"], b["unitarity"], b["reciprocity"], b["dT1_refinement"]) == (1e-10, 2e-10, 1e-10, 1e-8, 2e-4)


def _sha(p):
    return "sha256:" + hashlib.sha256(open(p, "rb").read()).hexdigest()


def _tree(d):
    h = hashlib.sha256()
    n = 0
    for root, _, files in sorted(os.walk(d)):
        for f in sorted(files):
            p = os.path.join(root, f)
            h.update((os.path.relpath(p, d) + ":" + _sha(p) + "\n").encode())
            n += 1
    return "sha256:" + h.hexdigest(), n


def test_protected_solver_files_and_prior_artifacts_are_unchanged():
    for f, h in PR["protected_files_sha256"].items():
        assert _sha(f) == h, f"{f} changed"
    for name, rec in PR["protected_results_trees"].items():
        digest, n = _tree(os.path.join("results", name))
        assert (digest, n) == (rec["sha256"], rec["n_files"]), f"results/{name} changed"
    for f, h in PR["b3x_artifacts_sha256"].items():
        assert _sha(os.path.join("results/rgf_b3x_validation", f)) == h, f


def test_synthetic_energies_are_not_b3y_points():
    pts = [p["E"] for p in PR["points"]]
    for E in E_SYN.values():
        assert min(abs(E - p) for p in pts) > 1e-3


# ----------------------------------------------------------------------
# X0: platform precision
# ----------------------------------------------------------------------
def test_X0_platform_precision_here_and_simulated_double_longdouble(monkeypatch):
    info = xp.platform_precision()
    assert info["longdouble_mantissa_bits"] >= 64 and info["evaluable"] and info["machine"] and info["system"]
    xp.require_extended_precision()
    real = np.finfo

    def fake53(dt):                                                  # Windows/MSVC or Apple arm64: longdouble == float64
        return types.SimpleNamespace(nmant=52, eps=2.220446049250313e-16) if dt is np.longdouble else real(dt)

    monkeypatch.setattr(xp.np, "finfo", fake53)
    bad = xp.platform_precision()
    assert not bad["evaluable"] and bad["longdouble_mantissa_bits"] == 53 and "53" in bad["reason"]
    with pytest.raises(xp.ExtendedPrecisionUnavailable, match="NOT EVALUABLE"):
        xp.require_extended_precision()
    Ep = eperp(3)
    V = np.zeros((10, 3, 3))
    for call in (lambda: xp.lead_data_extended(60.0, Ep, 0.05), lambda: xp.potential_blocks_extended(np.zeros((4, 7)), 3, 7),
                 lambda: xp.UnmergedOperator(60.0, Ep, V, 0.05, np.ones(3)),
                 lambda: xp.exact_lattice_S_from_blocks(60.0, Ep, V, 0.05)):
        with pytest.raises(xp.ExtendedPrecisionUnavailable):         # never silently falls back to double
            call()

    def fake113(dt):
        return types.SimpleNamespace(nmant=112, eps=1.9e-34) if dt is np.longdouble else real(dt)

    monkeypatch.setattr(xp.np, "finfo", fake113)
    assert xp.platform_precision()["evaluable"]                      # IEEE quad (aarch64 Linux) is acceptable


# ----------------------------------------------------------------------
# X1: lead data vs mpmath
# ----------------------------------------------------------------------
@pytest.mark.parametrize("inv", [512, 1024, 2048])
@pytest.mark.parametrize("key", ["A", "C"])
def test_X1_lead_roots_and_gamma_against_mpmath(inv, key):
    Ep = eperp(63)
    r = xp.lead_data_vs_mpmath(E_SYN[key], Ep, 1.0 / inv)
    assert r["passed"] and max(r["rel_lambda"], r["rel_gamma"]) <= 1e-17
    ld = xp.lead_data_extended(E_SYN[key], Ep, 1.0 / inv)
    assert int(ld["open"].sum()) == int(np.sum(Ep < E_SYN[key]))
    assert np.all(ld["gamma"][~ld["open"]] == 0) and np.all(ld["lam"][~ld["open"]].imag == 0) and np.all(np.abs(ld["lam"][~ld["open"]]) < 1)
    assert float(np.max(np.abs(np.abs(ld["lam"][ld["open"]]) - 1))) <= 1e-18


def test_X1_threshold_energy_is_rejected_and_double_formula_loses_precision():
    Ep = eperp(63)
    with pytest.raises(ValueError):
        xp.lead_data_extended(float(Ep[3]), Ep, 1 / 512)
    # why the module exists: Gamma from c = 1 - x in double has a relative error ~ eps/x that grows as 1/dx^2
    E, dx = E_SYN["A"], 1 / 2048
    ref = xp.lead_data_mpmath(E, Ep, dx)
    n = 3                                                           # channel 4, nearest threshold
    c = 1.0 - (E - Ep[n]) * dx * dx
    gam_double = 2 * (1 / (2 * dx * dx)) * math.sqrt(1.0 - c * c)
    with mp.workdps(50):
        err_double = float(abs(mp.mpf(gam_double) / ref["gamma"][n] - 1))
    ext = xp.lead_data_vs_mpmath(E, Ep, dx)
    assert err_double > 1e-12 and ext["rel_gamma"] <= 1e-17


# ----------------------------------------------------------------------
# potential blocks and the operator
# ----------------------------------------------------------------------
def test_extended_potential_blocks_against_direct_mpmath_sine_sums_and_y_uniform():
    ny, nq = 5, 11
    rng = np.random.default_rng(3)
    Vs = rng.uniform(-40, 40, (3, nq))
    B = xp.potential_blocks_extended(Vs, ny, nq)
    assert B.dtype == np.longdouble and np.max(np.abs(B - np.transpose(B, (0, 2, 1)))) <= 1e-18 * 40
    with mp.workdps(50):
        for (j, n, m) in ((0, 1, 1), (1, 2, 4), (2, 5, 3)):
            ref = sum(mp.mpf(float(Vs[j, l - 1])) * 2 / (nq + 1) * mp.sin(mp.pi * n * l / (nq + 1)) * mp.sin(mp.pi * m * l / (nq + 1)) for l in range(1, nq + 1))
            got = mp.mpf(np.format_float_positional(B[j, n - 1, m - 1], precision=24, unique=False))
            assert abs(got - ref) <= 1e-17 * 40
    u = np.tile(np.array([3.0, -7.0])[:, None], (1, nq))
    ub = xp.potential_blocks_extended(u, ny, nq)
    assert np.max(np.abs(ub - np.array([3.0, -7.0], dtype=np.longdouble)[:, None, None] * np.eye(ny, dtype=np.longdouble))) <= 1e-17 * 7
    with pytest.raises(ValueError):
        xp.potential_blocks_extended(Vs, 6, nq)                    # n_quad < 2 n_modes


def test_unmerged_operator_differs_from_the_merged_double_matrix_at_the_ulp_of_2t():
    dx, N, ny, E = 1 / 2048, 40, 3, 88.5
    Ep = eperp(ny)
    rng = np.random.default_rng(7)
    A = rng.uniform(-30, 30, (N, ny, ny)); V = 0.5 * (A + A.transpose(0, 2, 1))
    lead = xp.lead_data_extended(E, Ep, dx)
    op = xp.UnmergedOperator(E, Ep, V, dx, lead["lam"])
    M = xp.assemble_preconditioner(E, Ep, V, dx, lead["lam"])
    X = (rng.standard_normal((N, ny, 2)) + 1j * rng.standard_normal((N, ny, 2))).astype(xp.CLD)
    y_unmerged = op.apply(X)
    y_rounded = xp.apply_matrix_extended(M, X.reshape(N * ny, 2)).reshape(N, ny, 2)
    d = float(np.max(np.abs(y_unmerged - y_rounded)))
    scale = float(np.max(np.abs(X)))
    ulp = np.spacing(2 * float(xp.hopping(dx)))
    assert 1e-13 * scale < d < 4 * ulp * scale * 4               # visible, and bounded by a few ulps of the merged diagonal (9.3e-10 here)
    # in extended precision the unmerged operator applied to a smooth field is accurate: A x for x = const is P x + t*(ghost terms) only
    ones = np.ones((N, ny, 1), dtype=xp.CLD)
    Y = op.apply(ones)
    interior = Y[1:-1, :, 0]
    expected = (op.P[1:-1].astype(xp.CLD) @ np.ones(ny, dtype=xp.CLD))
    assert np.max(np.abs(interior - expected)) <= 1e-18 * 100        # second difference of a constant vanishes exactly


def test_apply_matrix_extended_is_the_exact_product_of_a_rounded_matrix():
    rng = np.random.default_rng(5)
    n = 30
    M = sp.random(n, n, density=0.2, random_state=5, format="csr", dtype=np.float64) + 1j * sp.random(n, n, density=0.2, random_state=6, format="csr")
    M = (M + sp.identity(n)).tocsr()
    X = (rng.standard_normal((n, 2)) + 1j * rng.standard_normal((n, 2))).astype(xp.CLD)
    Y = xp.apply_matrix_extended(M, X, chunk_rows=7)
    with mp.workdps(50):
        for i in (0, 11, 29):
            for c in range(2):
                ref = sum(mp.mpc(float(M[i, j].real), float(M[i, j].imag)) * to_mp(X[j, c]) for j in M[i].indices)
                assert abs(to_mp(Y[i, c]) - ref) <= 1e-17 * (1 + abs(ref))


# ----------------------------------------------------------------------
# helpers for exact-lattice runs
# ----------------------------------------------------------------------
def run_xp(E, V_samples, inv, ny=63, nq=127, **kw):
    dx, N, _ = geometry(inv)
    return xp.exact_lattice_S(E, eperp(ny), V_samples, dx, ny, nq, **kw)


def check_gates(res, tag):
    g = res["refinement_gate"]
    assert g["passed"], (tag, g)                                                      # X4
    assert res["green_refinement_gate"]["passed"], (tag, res["green_refinement_gate"])
    assert res["fl_vs_state_max_abs"] <= xp.FL_WF_TOL, (tag, res["fl_vs_state_max_abs"])   # X5
    assert res["error_indicator"] <= xp.ERROR_INDICATOR_MAX, (tag, res["error_indicator"])  # X7
    assert res["eta_final"] <= xp.ETA_MAX


# ----------------------------------------------------------------------
# X2: clean guide
# ----------------------------------------------------------------------
@pytest.mark.parametrize("inv,key", [(512, "A"), (512, "C"), (1024, "A"), (2048, "A")])
def test_X2_clean_guide_matches_mpmath_phase_and_has_no_reflection_or_mixing(inv, key):
    dx, N, _ = geometry(inv)
    Ep = eperp(63)
    res = run_xp(E_SYN[key], np.zeros((N, 127)), inv, condition=(inv != 2048) or key == "A")
    check_gates(res, f"clean {inv} {key}")
    ref = xp.lead_data_mpmath(E_SYN[key], Ep, dx)
    no = res["n_open"]
    assert no == int(np.sum(Ep < E_SYN[key]))
    t, r = res["t_RL"], res["r_L"]
    with mp.workdps(50):
        for n in range(no):
            exact = ref["lam"][n] ** (N - 1)
            assert abs(to_mp(t[n, n]) - exact) <= xp.SCALAR_ABS_TOL
            assert abs(to_mp(res["t_LR"][n, n]) - exact) <= xp.SCALAR_ABS_TOL
    off = t - np.diag(np.diag(t))
    assert float(np.max(np.abs(off))) <= xp.SCALAR_ABS_TOL and float(np.max(np.abs(r))) <= xp.SCALAR_ABS_TOL and float(np.max(np.abs(res["r_R"]))) <= xp.SCALAR_ABS_TOL
    U = float(np.max(np.abs(res["S"].conj().T @ res["S"] - np.eye(2 * no))))
    assert U <= 1e-15


# ----------------------------------------------------------------------
# X3: smooth y-uniform barrier vs single-channel transfer matrix; X6: ordering stability
# ----------------------------------------------------------------------
@pytest.mark.parametrize("inv,key", [(512, "B"), (1024, "B"), (2048, "A")])
def test_X3_smooth_barrier_matches_single_channel_transfer_matrix(inv, key):
    dx, N, x0 = geometry(inv)
    xs = x0 + dx * np.arange(N)
    Vx = barrier(xs)
    Ep = eperp(63)
    res = run_xp(E_SYN[key], np.tile(Vx[:, None], (1, 127)), inv)
    check_gates(res, f"barrier {inv} {key}")
    Vb = res["V_blocks"]
    assert float(np.max(np.abs(Vb - Vb[:, np.arange(63), np.arange(63)][:, :, None] * np.eye(63)[None]))) <= 1e-17 * 150   # channels decouple
    no = res["n_open"]
    assert no == int(np.sum(Ep < E_SYN[key]))
    with mp.workdps(50):
        for n in range(no):
            r_ex, t_ex = xp.mp_single_channel_amplitudes(E_SYN[key], Ep[n], Vx, dx)
            assert abs(to_mp(res["t_RL"][n, n]) - t_ex) <= xp.SCALAR_ABS_TOL, (inv, n)
            assert abs(to_mp(res["r_L"][n, n]) - r_ex) <= xp.SCALAR_ABS_TOL, (inv, n)
    t = res["t_RL"]
    assert float(np.max(np.abs(t - np.diag(np.diag(t))))) <= xp.SCALAR_ABS_TOL                   # no channel mixing
    S = res["S"]
    assert float(np.max(np.abs(S.conj().T @ S - np.eye(2 * no)))) <= 1e-15
    assert float(np.max(np.abs(S - S.T))) <= 1e-15
    assert res["kappa1"] > 1e3 and np.isfinite(res["error_indicator"])


def test_X6_preconditioner_ordering_does_not_change_the_extended_result():
    dx, N, x0 = geometry(512)
    Vx = barrier(x0 + dx * np.arange(N))
    Ep = eperp(63)
    Vb = xp.potential_blocks_extended(np.tile(Vx[:, None], (1, 127)), 63, 127)
    r = xp.ordering_stability(E_SYN["B"], Ep, Vb, dx)
    assert r["passed"] and r["max_abs_difference"] <= xp.ORDERING_TOL and r["orderings"] == ["COLAMD", "MMD_AT_PLUS_A"]


def test_independent_green_solve_and_reuse_agree():
    dx, N, x0 = geometry(512)
    Vx = barrier(x0 + dx * np.arange(N))
    Vb = xp.potential_blocks_extended(np.tile(Vx[:, None], (1, 127)), 63, 127)
    a = xp.exact_lattice_S_from_blocks(E_SYN["B"], eperp(63), Vb, dx, independent_green_solve=True, condition=False)
    b = xp.exact_lattice_S_from_blocks(E_SYN["B"], eperp(63), Vb, dx, independent_green_solve=False, condition=False)
    assert float(np.max(np.abs(a["S_fisher_lee"] - b["S_fisher_lee"]))) <= 1e-15 and "green_refinement_gate" in a and "green_refinement_gate" not in b


# ----------------------------------------------------------------------
# dense mixing (multichannel): exact mpmath dense solve, and the other solvers on a toy
# ----------------------------------------------------------------------
@pytest.mark.parametrize("dx,N", [(1 / 512, 24), (1 / 2048, 30)])
def test_dense_channel_mixing_against_an_exact_mpmath_dense_solve(dx, N):
    ny, E = 3, 25.0                                                  # 2 open, 1 closed channel
    Ep = eperp(ny)
    rng = np.random.default_rng(11)
    A = rng.uniform(-20, 20, (N, ny, ny)); V = 0.5 * (A + A.transpose(0, 2, 1))
    res = xp.exact_lattice_S_from_blocks(E, Ep, V, dx)
    check_gates(res, f"dense {dx}")
    S_mp = xp.mp_dense_S(E, Ep, V, dx)
    assert S_mp.shape == res["S"].shape == (4, 4)
    with mp.workdps(50):
        d = max(abs(to_mp(res["S"][i, k]) - S_mp[i, k]) for i in range(4) for k in range(4))
    assert d <= xp.SCALAR_ABS_TOL
    assert float(np.max(np.abs(res["t_RL"] - np.diag(np.diag(res["t_RL"]))))) > 1e-6            # genuine mixing is present


def test_dense_toy_agrees_with_the_rgf_and_both_sparse_references_without_importing_them_in_the_module():
    import smatrix_spec  # noqa: F401  (the test, not the module, imports the existing solvers)
    import sparse_reference as sr
    import sparse_wavefunction as sw
    from stationary import Device, scattering_matrix, transverse_basis
    N, ny, dx, E = 10, 6, 0.05, 60.0
    rng = np.random.default_rng(5)
    A = rng.uniform(-40, 40, (N, ny, ny)); V = 0.5 * (A + A.transpose(0, 2, 1))
    basis = transverse_basis(1.0, ny)
    assert np.array_equal(basis.E_perp, eperp(ny))                    # bitwise identical E_perp doubles
    res = xp.exact_lattice_S_from_blocks(E, basis.E_perp, V, dx)
    S = res["S"].astype(np.complex128)
    rg = scattering_matrix(E, Device(xs=dx * np.arange(N), dx=dx, basis=basis, V=V, V_lead_max=0.0))["S"]
    wf = sw.wavefunction_S_from_blocks(E, V, basis.E_perp, dx)["S"]
    gb = sr.sparse_S_from_blocks(E, V, basis.E_perp, dx)["S"]
    for other in (rg, wf, gb):
        assert float(np.max(np.abs(S - other))) <= 1e-12
    assert float(np.max(np.abs(S.conj().T @ S - np.eye(6)))) <= 1e-14 and float(np.max(np.abs(S - S.T))) <= 1e-14
    assert res["n_open"] == 3 and res["S_fisher_lee"].shape == res["S"].shape


# ----------------------------------------------------------------------
# error decomposition with a stand-in double-precision method (single channel, fine dx)
# ----------------------------------------------------------------------
def _naive_lead(E, Ep, dx):
    t = 1.0 / (2 * dx * dx)
    c = 1.0 - (E - Ep) * dx * dx
    s = math.sqrt(1.0 - c * c)
    return complex(c, s), 2 * t * s


def _extract_two_columns(x, gamma, a):
    """x: (N, 2) solution columns [left source, right source] -> S (2x2) with the wave-function formulas; works in any precision."""
    sg = np.sqrt(gamma)
    r_L = sg * (x[0, 0] - a); t_LR = sg * x[0, 1]; t_RL = sg * x[-1, 0]; r_R = sg * (x[-1, 1] - a)
    return np.array([[r_L, t_LR], [t_RL, r_R]])


@pytest.mark.parametrize("inv", [512, 2048])
def test_error_decomposition_separates_representation_from_solver_error(inv):
    dx, N, x0 = geometry(inv)
    xs = x0 + dx * np.arange(N)
    Vx = 5.0 * 0.5 * (np.tanh((xs - 2.35) / 0.01) - np.tanh((xs - 2.45) / 0.01)) + 3.0 * 0.5 * (np.tanh((xs - 2.6) / 0.02) - np.tanh((xs - 2.75) / 0.02))
    Ep = np.array([0.5 * (4 * np.pi) ** 2]); E = 88.5134459632352                   # channel-4-like threshold; a synthetic single channel
    V = Vx.reshape(N, 1, 1)
    ex = xp.exact_lattice_S_from_blocks(E, Ep, V, dx)
    check_gates(ex, f"chain {inv}")
    S_exact = ex["S"]
    with mp.workdps(50):                                                              # the exact answer itself, from mpmath
        r_ex, t_ex = xp.mp_single_channel_amplitudes(E, Ep[0], Vx, dx)
        assert abs(to_mp(S_exact[1, 0]) - t_ex) <= xp.SCALAR_ABS_TOL and abs(to_mp(S_exact[0, 0]) - r_ex) <= xp.SCALAR_ABS_TOL
    # ---- the stand-in double method: merged matrix, double lead data, double LU, double extraction ----
    lam_d, gam_d = _naive_lead(E, Ep[0], dx)
    t = float(xp.hopping(dx))
    M = xp.assemble_preconditioner(E, Ep, V, dx, np.array([lam_d]))
    import scipy.sparse.linalg as spla
    lu = spla.splu(M.tocsc())
    a_d = gam_d ** -0.5
    src_d = -t * (1 / lam_d - lam_d) * a_d
    Bd = np.zeros((N, 2), complex); Bd[0, 0] = src_d; Bd[-1, 1] = src_d
    x_d = lu.solve(Bd)
    S_double = _extract_two_columns(x_d, gam_d, a_d)
    # ---- solutions of the ROUNDED matrix to extended precision ----
    B_fl = Bd.astype(xp.CLD)
    fl = xp.solve_matrix_extended(M, B_fl)
    assert fl["refinement_gate"]["passed"]
    S_fl = _extract_two_columns(fl["X"], xp.LD(gam_d), xp.LD(a_d))                     # the method's DOUBLE lead data, evaluated exactly
    lead = ex["lead"]
    hyb = xp.solve_matrix_extended(M, xp.build_sources(lead, N, 1).reshape(N, 2),
                                   corner_swap=dict(t=t, ny=1, N=N, lam_double=np.array([lam_d]), lam_ext=lead["lam"]))
    assert hyb["refinement_gate"]["passed"]
    a_x = (lead["gamma"][0] ** xp.LD(-0.5))
    S_hyb = _extract_two_columns(hyb["X"], lead["gamma"][0], a_x)                      # rounded bulk, EXACT lead data
    dec = xp.error_decomposition(S_double, S_exact, S_fl, S_hyb)
    assert dec["identity_residual"] <= 1e-18                                           # the telescoping identity
    assert dec["total"] <= 1e-9 and np.isfinite(dec["bulk_rounding"])
    assert abs(dec["total"] - float(np.max(np.abs(S_double - S_exact)))) <= 1e-18
    terms = dict(bulk_rounding=dec["bulk_rounding"], lead_data=dec["lead_data"], solver=dec["solver"])
    dom = xp.representation_dominance(terms, ex["error_indicator"])
    assert dom["available"] and dom["supported"] == (dom["resolved"] and dom["dominant"])
    if inv == 2048:
        assert dec["total"] > 1e-12                                                    # the representation phenomenon is visible at fine dx


# ----------------------------------------------------------------------
# X4 gate logic, e_ind, classification rules (pure functions)
# ----------------------------------------------------------------------
def test_X4_gate_logic_detects_slow_nonconverging_and_unconverged_refinement():
    def res(etas, dns, converged=True, eta_final=1e-20):
        h = [dict(k=i, eta=e, correction_inf=d, x_inf=1.0, contraction=(dns[i] / dns[i - 1]) if i else float("nan")) for i, (e, d) in enumerate(zip(etas, dns))]
        return dict(history=h, converged=converged, n_iterations=len(h), eta_final=eta_final)

    good = xp.refinement_gate(res([1e-16, 1e-20], [1e-14, 1e-19]))
    assert good["passed"] and good["worst_contraction"] == pytest.approx(1e-5)
    slow = xp.refinement_gate(res([1e-12, 1e-13, 1e-20], [1e-9, 5e-10, 1e-19]))        # contraction 0.5 while above the target
    assert not slow["passed"] and slow["worst_contraction"] == pytest.approx(0.5)
    assert not xp.refinement_gate(res([1e-16, 1e-20], [1e-14, 1e-19], converged=False))["passed"]
    assert not xp.refinement_gate(res([1e-16, 1e-20], [1e-14, 1e-19], eta_final=2e-17))["passed"]          # operator residual above 1e-17
    assert not xp.refinement_gate(res([1e-16] * 11, [1e-14] * 11))["passed"]                                # more than 10 iterations
    boundary = xp.refinement_gate(res([1e-16, 1e-20], [1e-14, 1e-15], eta_final=1e-17))                    # contraction exactly 0.1, eta exactly 1e-17
    assert boundary["passed"]


def test_error_indicator_is_the_maximum_of_the_iterate_difference_and_kappa_times_eta():
    dx, N, x0 = geometry(512)
    V = np.zeros((40, 3, 3))
    res = xp.exact_lattice_S_from_blocks(60.0, eperp(3), V, 0.05)
    d = res["history"][-1]["dS_vs_previous"]
    assert res["error_indicator"] == pytest.approx(max(d, res["kappa1"] * res["eta_final"]), rel=1e-12)
    assert "not a rigorous" in res["error_indicator_note"]


def test_classification_rules():
    ok = {512: 2.0 ** -40, 1024: 2.0 ** -38, 2048: 10 * 2.0 ** -38}   # growth exactly 4, then exactly 10 (the limit; binary-exact values): consistent
    assert xp.precision_consistency(ok)["passed"] and xp.precision_consistency(ok)["growth_per_halving"] == pytest.approx([4.0, 10.0])
    assert not xp.precision_consistency({512: 5e-12, 1024: 6.0001e-11})["passed"]               # growth 12.0002 > 10
    assert xp.precision_consistency({512: 5e-12, 1024: 4.9999e-11})["passed"]                   # growth 9.9998 <= 10
    assert xp.precision_consistency({512: 1.0000001e-10})["passed"] is False                    # eps(1/512) just above the target
    assert xp.precision_consistency({512: 1e-10})["passed"]                                      # eps(512) exactly at the target passes
    assert not xp.precision_consistency({1024: 1e-11, 2048: 2e-11})["passed"]                    # no 1/512 anchor
    terms = dict(bulk_rounding=2.0e-11, lead_data=1.0e-11, solver=1.0e-11)
    assert xp.representation_dominance(terms, 1e-13)["supported"]                                # exactly factor 2 passes
    assert not xp.representation_dominance(dict(bulk_rounding=1.99e-11, lead_data=1e-11, solver=1e-11), 1e-13)["supported"]
    assert not xp.representation_dominance(terms, 1e-12)["supported"]                            # 1e-11 is exactly 10 x e_ind: not strictly greater -> unresolved
    assert xp.representation_dominance(terms, 9.9e-13)["supported"]
    none = xp.representation_dominance(None, 1e-13)
    assert not none["available"] and not none["supported"]                                       # M1: no decomposition, no claim
    # full classification: never PASS; never claims dominance without the decomposition
    assert xp.classify_exceedance(ok, terms, 1e-13, claim_representation_dominance=False) == "PRECISION-CONSISTENT"
    assert xp.classify_exceedance(ok, terms, 1e-13, True) == "PRECISION-LIMITED (measured; representation-dominated)"
    assert xp.classify_exceedance(ok, None, 1e-13, True) == "UNEXPLAINED"                        # M1 attempting a mechanism claim
    assert xp.classify_exceedance(ok, dict(bulk_rounding=1e-11, lead_data=9e-12, solver=9e-12), 1e-13, True) == "UNEXPLAINED"   # decomposition does not support it
    assert xp.classify_exceedance({512: 5e-12, 1024: 1e-10}, terms, 1e-13, False) == "UNEXPLAINED"   # growth 20 per halving
    assert xp.classify_exceedance({512: 5e-12, 1024: 1e-10}, terms, 1e-13, True) == "UNEXPLAINED"
    for s in (xp.classify_exceedance(ok, terms, 1e-13, c) for c in (False, True)):
        assert s != "PASS"


def test_error_decomposition_identity_on_arbitrary_matrices():
    z = np.zeros((4, 4), dtype=complex)
    d = xp.error_decomposition(3e-11 * (z + 1), z, 2e-11 * (z + 1), 0.5e-11 * (z + 1))      # S_double, S_exact, S_fl, S_hybrid
    assert d["total"] == pytest.approx(3e-11, rel=1e-12) and d["identity_residual"] <= 1e-25
    assert d["bulk_rounding"] == pytest.approx(0.5e-11, rel=1e-12)           # S_hybrid - S_exact
    assert d["lead_data"] == pytest.approx(1.5e-11, rel=1e-12)               # S_fl - S_hybrid
    assert d["solver"] == pytest.approx(1.0e-11, rel=1e-9)                   # S_double - S_fl
    assert d["representation_total"] == pytest.approx(2e-11, rel=1e-12)
    assert d["bulk_rounding"] + d["lead_data"] + d["solver"] == pytest.approx(d["total"], rel=1e-9)
