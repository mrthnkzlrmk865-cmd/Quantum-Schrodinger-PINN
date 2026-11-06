"""Benchmark B3: RGF vs independent direct sparse Green-function reference on
fixed disorder realizations (validate_rgf_b3.py, sparse_reference.py,
docs/RGF_B3_VALIDATION.md).

Part 1 tests the machinery on small live problems (including the sparse
solver on synthetic cases with known answers). Part 2 checks the committed
artifacts in results/rgf_b3_validation (skipped if absent). Predeclared
thresholds are checked for preservation only, never changed.
"""
import ast
import csv
import json
import math
import os

import numpy as np
import pytest

import smatrix_spec as sm
import sparse_reference as sr
import validate_rgf_b3 as b3
from disorder import disorder_from_config, make_disorder
from potential import Potential
from stationary import Device, scattering_matrix, transverse_basis

RES = "results/rgf_b3_validation"
SUMMARY = os.path.join(RES, "summary.json")
HAVE = os.path.exists(SUMMARY) and os.path.exists(os.path.join(RES, "rgf_sparse_comparison.csv"))
SMALL = dict(dx=1 / 128, Ny=15, Nq=31, buffer=0.2)        # fast device for machinery tests


# ----------------------------------------------------------------------
# protocol, configuration, energies, thresholds
# ----------------------------------------------------------------------
def test_predeclared_protocol_preserved():
    assert b3.LIMITS == dict(S_agreement=1e-10, U=1e-10, Q=1e-8, dT1=2e-4)
    assert b3.CASES == {"A": dict(seed=0, W=100.0, Ld=1.0, zeta=0.02), "B": dict(seed=0, W=100.0, Ld=1.0, zeta=0.05),
                        "C": dict(seed=0, W=100.0, Ld=1.0, zeta=0.20), "D": dict(seed=0, W=200.0, Ld=1.0, zeta=0.05)}
    assert b3.BASE == dict(dx=1 / 512, Ny=63, Nq=127, buffer=0.2)
    assert b3.VARIANTS["dx"] == dict(b3.BASE, dx=1 / 1024)
    assert b3.VARIANTS["Nq"]["Nq"] == 2 * b3.BASE["Nq"] and b3.VARIANTS["buffer"]["buffer"] == 0.4
    assert b3.VARIANTS["Ny127"]["Ny"] == 127 and b3.VARIANTS["Ny255"]["Ny"] == 255
    assert b3.KX_OFFSETS == (-2, -1, 0, 1, 2) and b3.THRESHOLD_REL == 1e-3
    # B1 / B2 thresholds untouched
    import validate_rgf_b1 as b1
    import validate_rgf_b2 as b2
    assert b1.LIMITS == dict(T_diag=1e-12, mixing=1e-12, R=1e-12, U=1e-10, phase=1e-10, Q=1e-8)
    assert b2.LIMITS == dict(dT=2e-4, mixing=1e-12, U=1e-10, Q=1e-8)


def test_pass_fail_at_exact_thresholds():
    for lim in (1e-10, 1e-8, 2e-4):
        assert b3.passes(lim, lim)                               # equality passes
        assert b3.passes(0.0, lim) and b3.passes(lim * (1 - 1e-12), lim)
        assert not b3.passes(lim * (1 + 1e-12), lim)
        assert not b3.passes(float("nan"), lim) and not b3.passes(float("inf"), lim)
    rows = [dict(case="A", kx=1.0, E=2.0, v=1e-10), dict(case="B", kx=1.0, E=2.0, v=1e-10 * (1 + 1e-9))]
    c = b3._crit(rows, "v", 1e-10)
    assert c["per_case_pass"] == {"A": True, "B": False} and not c["passed"] and c["argmax"]["case"] == "B"
    assert b3._crit(rows[:1], "v", 1e-10)["passed"]
    assert not b3._crit([dict(case="A", kx=1.0, E=2.0, v=float("nan"))], "v", 1e-10)["passed"]


@pytest.mark.parametrize("case", list("ABCD"))
def test_configuration_and_seed(case):
    p = b3.CASES[case]
    cfg, pot = b3.case_potential(case)
    assert cfg.disorder_strength == p["W"] and cfg.disorder_corr_length == p["zeta"]
    assert tuple(cfg.disorder_region) == (1.9, 2.9) and cfg.Ly == 1.0
    d = pot.disorder
    assert d.seed == 0 and d.strength == p["W"] and d.corr_length == p["zeta"] and d.n_modes == 512
    ref = make_disorder(p["W"], p["zeta"], 0, n_modes=512, region=(1.9, 2.9), edge=cfg.disorder_edge)
    assert np.array_equal(d.qx, ref.qx) and np.array_equal(d.phi, ref.phi)
    # same seed = same random draw rescaled by zeta (B2/B3 pairing property)
    other = disorder_from_config(b3.case_config("B"), seed=0)
    assert np.allclose(d.qx * p["zeta"], other.qx * 0.05, rtol=1e-14)


@pytest.mark.parametrize("case", list("ABCD"))
def test_energies_follow_packet_convention(case):
    cfg = b3.case_config(case)
    E = b3.case_energies(cfg)
    sk = 1 / (math.sqrt(2) * cfg.sigma_x)
    assert len(E) == 5
    for e, j in zip(E, (-2, -1, 0, 1, 2)):
        assert e["kx"] == cfg.kx0 + j * sk
        assert e["E"] == 0.5 * (e["kx"] ** 2 + math.pi ** 2)
        assert e["rel_dist_to_threshold"] >= 1e-3                # threshold rule would not apply
        assert e["n_open_exact"] == 4 + (0, 1, 2, 3, 4)[e["j"] + 2]
    assert [e["n_open_exact"] for e in E] == [4, 5, 6, 7, 8]


def test_threshold_rule_rejects_energy_on_threshold(monkeypatch):
    cfg = b3.case_config("A")
    monkeypatch.setattr(b3, "KX_OFFSETS", (0,))
    monkeypatch.setattr(b3, "THRESHOLD_REL", 0.5)                # absurdly wide: must raise, never silently drop
    with pytest.raises(RuntimeError):
        b3.case_energies(cfg)


def test_geometry_buffers_and_nested_grids():
    cfg = b3.case_config("B")
    g = {(dx, b): b3.geometry(cfg, dx, b) for dx in (1 / 512, 1 / 1024) for b in (0.2, 0.4)}
    base, fine, wide = g[(1 / 512, 0.2)], g[(1 / 1024, 0.2)], g[(1 / 512, 0.4)]
    assert base["n_slices"] == 719 and fine["n_slices"] == 1437 and wide["n_slices"] == 923
    for (dx, b), v in g.items():
        assert v["buffer_left"] >= b and v["buffer_left"] == pytest.approx(math.ceil(b * 512) / 512, abs=1e-15)
        assert v["x_start"] + (v["n_slices"] - 1) * dx == pytest.approx(v["x_end"], abs=1e-12)
        assert v["x_start"] < 1.9 < 2.9 < v["x_end"]
    xs = lambda v, dx: v["x_start"] + dx * np.arange(v["n_slices"])
    # nested lattices: every base slice is also a slice of the refined / wider device
    assert np.min(np.abs(xs(fine, 1 / 1024)[None, :] - xs(base, 1 / 512)[:, None]), axis=1).max() < 1e-12
    assert np.min(np.abs(xs(wide, 1 / 512)[None, :] - xs(base, 1 / 512)[:, None]), axis=1).max() < 1e-12
    assert fine["buffer_left"] == base["buffer_left"]             # a dx refinement does not change the buffer length
    assert wide["buffer_left"] == pytest.approx(205 / 512) and base["buffer_left"] == pytest.approx(103 / 512)
    # the physical disorder region does not depend on the buffer
    assert cfg.disorder_region == (1.9, 2.9)


# ----------------------------------------------------------------------
# same physical realization
# ----------------------------------------------------------------------
@pytest.mark.parametrize("case", list("ABCD"))
def test_same_realization_bitwise(case):
    cfg, pot = b3.case_potential(case)
    dev, geo, info, ref = b3.make_rgf_device(case, cfg, pot, SMALL)
    assert info["bitwise_equal"] and info["coords_bitwise_equal"] and info["max_abs_diff"] == 0.0
    assert info["shape"] == [dev.n_slices, SMALL["Nq"]] and info["potential_hash"] == b3.sha(ref)
    # a completely fresh Potential/disorder object gives the same numbers
    pot2 = Potential(cfg, disorder_from_config(cfg, seed=0))
    X, Y = np.meshgrid(dev.xs, dev.basis.y, indexing="ij")
    assert np.array_equal(pot2.numpy(X, Y), ref)
    # the sparse reference samples the same realization through its own coordinates
    out = sr.sparse_scattering_matrix(b3.case_energies(cfg)[0]["E"], pot.numpy, geo["x_start"], SMALL["dx"],
                                      dev.n_slices, 1.0, SMALL["Ny"], SMALL["Nq"])
    assert np.array_equal(out["V_samples"], ref) and np.array_equal(out["xs"], dev.xs)
    # and it is a real disorder realization (not zero), different for different (W, zeta)
    assert ref.std() > 1.0


def test_realizations_differ_between_cases_and_comparison_detects_changes():
    refs = {}
    for c in "ABD":
        cfg, pot = b3.case_potential(c)
        refs[c] = b3.make_rgf_device(c, cfg, pot, SMALL)[3]
    assert b3.sha(refs["A"]) != b3.sha(refs["B"]) != b3.sha(refs["D"])
    a = refs["B"].copy()
    assert b3.compare_arrays(a, refs["B"])["bitwise_equal"]
    a[3, 4] *= 1 + 1e-15
    c = b3.compare_arrays(a, refs["B"])
    assert not c["bitwise_equal"] and 0 < c["max_abs_diff"] < 1e-10


# ----------------------------------------------------------------------
# sparse reference: construction, independence, synthetic correctness
# ----------------------------------------------------------------------
def _toy(n_slices=7, ny=4, dx=0.1, seed=3, vmax=20.0):
    rng = np.random.default_rng(seed)
    A = rng.uniform(-vmax, vmax, (n_slices, ny, ny))
    V = 0.5 * (A + np.transpose(A, (0, 2, 1)))
    E_perp = 0.5 * (np.arange(1, ny + 1) * np.pi) ** 2
    return V, E_perp, dx


def test_sparse_matrix_construction():
    V, Ep, dx = _toy()
    N, ny = V.shape[:2]
    E = 20.0
    modes = sm.lead_modes(E, Ep, dx)
    M = sr.assemble_M(E, V, Ep, dx, modes["sigma"])
    t = sm.hopping(dx)
    D = np.zeros((N * ny, N * ny), dtype=complex)               # explicit dense block-tridiagonal M = E - H - Sigma
    for j in range(N):
        s = slice(j * ny, (j + 1) * ny)
        H = 2 * t * np.eye(ny) + np.diag(Ep) + V[j]
        D[s, s] = E * np.eye(ny) - H
        if j in (0, N - 1):
            D[s, s] -= np.diag(modes["sigma"])
        if j + 1 < N:
            D[s, (j + 1) * ny:(j + 2) * ny] = t * np.eye(ny)
            D[(j + 1) * ny:(j + 2) * ny, s] = t * np.eye(ny)
    # not bitwise: the diagonal is summed as (E - 2t) - E_perp - V here and E - (2t + E_perp + V) above
    assert np.max(np.abs(M.toarray() - D)) <= 8 * np.finfo(float).eps * np.max(np.abs(D))
    assert np.array_equal(M.toarray() != 0, D != 0)                  # identical sparsity pattern
    assert M.nnz <= N * ny * ny + 2 * (N - 1) * ny              # block-tridiagonal sparsity, no fill
    assert (M - M.T).nnz == 0                                    # complex symmetric (spec §2.7.8)


def test_sparse_solver_independent_of_rgf(monkeypatch):
    tree = ast.parse(open("sparse_reference.py").read())
    imported = {n.names[0].name for n in ast.walk(tree) if isinstance(n, ast.Import)} | \
               {n.module for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    assert "stationary" not in imported and "validate_rgf_b3" not in imported
    assert "rgf_endpoint_blocks" not in open("sparse_reference.py").read()
    import stationary
    boom = lambda *a, **k: (_ for _ in ()).throw(AssertionError("RGF recursion was called by the sparse reference"))
    monkeypatch.setattr(stationary, "rgf_endpoint_blocks", boom)
    V, Ep, dx = _toy()
    out = sr.sparse_S_from_blocks(20.0, V, Ep, dx)               # must not touch the RGF
    assert out["S"].shape == (2 * out["n_open"], 2 * out["n_open"])
    with pytest.raises(AssertionError):                           # and the patch does hit the RGF solver itself
        Dv = Device(xs=np.arange(7) * dx, dx=dx, basis=transverse_basis(1.0, 4), V=V, V_lead_max=0.0)
        scattering_matrix(20.0, Dv)


def test_sparse_vs_transfer_matrix_1d():
    """Ny = 1 chain: the sparse S equals a Green-function-free backward
    transfer-matrix wave function (complex amplitudes, convention P)."""
    dx, N, E = 0.05, 40, 30.0
    rng = np.random.default_rng(11)
    Vj = np.where(np.arange(N) % 9 < 4, rng.uniform(5, 25, N), 0.0)
    Ep = np.array([0.5 * np.pi ** 2])
    out = sr.sparse_S_from_blocks(E, Vj.reshape(N, 1, 1), Ep, dx)
    t = sm.hopping(dx)
    lam = sm.lead_modes(E, Ep, dx)["lam"][0]
    psi = {N + 1: lam, N: 1.0 + 0j}                              # right lead: d lam^{j-N}, d := 1
    for j in range(N, 0, -1):                                    # -t psi_{j-1} + (2t + Ep + V_j - E) psi_j - t psi_{j+1} = 0
        psi[j - 1] = ((2 * t + Ep[0] + Vj[j - 1] - E) * psi[j] - t * psi[j + 1]) / t
    # left lead: psi_j = a lam^{j-1} + b lam^{-(j-1)}  (convention P: a, b referred to slice 1)
    A = np.array([[1.0, 1.0], [1 / lam, lam]])
    a, b = np.linalg.solve(A, np.array([psi[1], psi[0]]))
    t_ref, r_ref = 1 / a, b / a
    assert abs(out["t_RL"][0, 0] - t_ref) < 1e-9 and abs(out["r_L"][0, 0] - r_ref) < 1e-9
    assert abs(abs(t_ref) ** 2 + abs(r_ref) ** 2 - 1) < 1e-9
    assert 0.05 < abs(t_ref) ** 2 < 0.95                          # non-trivial scatterer
    assert out["relative_residual"] < 1e-12


def test_sparse_clean_guide_and_residual():
    N, ny, dx, E = 12, 4, 0.1, 30.0
    Ep = 0.5 * (np.arange(1, ny + 1) * np.pi) ** 2
    out = sr.sparse_S_from_blocks(E, np.zeros((N, ny, ny)), Ep, dx)
    lam = out["modes"]["lam"][out["modes"]["open"]]
    assert np.allclose(np.diag(out["t_RL"]), lam ** (N - 1), atol=1e-11)
    assert np.max(np.abs(out["r_L"])) < 1e-11 and out["relative_residual"] < 1e-13


def test_sparse_equals_rgf_on_random_multichannel_device():
    V, Ep, dx = _toy(n_slices=10, ny=6, dx=0.05, seed=5, vmax=40.0)
    basis = transverse_basis(1.0, 6)
    assert np.array_equal(basis.E_perp, Ep)
    dev = Device(xs=0.0 + dx * np.arange(10), dx=dx, basis=basis, V=V, V_lead_max=0.0)
    E = 60.0                                                     # 3 open, 3 closed channels
    r = scattering_matrix(E, dev)
    s = sr.sparse_S_from_blocks(E, V, Ep, dx)
    assert r["n_open"] == s["n_open"] == 3 and r["S"].shape == (6, 6)
    assert np.max(b3.s_difference(r["S"], s["S"])) < 1e-10
    assert np.max(np.abs(np.abs(r["t_RL"][:, 0]) ** 2)) > 1e-4   # genuine channel mixing present
    assert np.sum(np.abs(sm.split_S(s["S"], 3)["t_RL"][1:, 0]) ** 2) > 1e-6


def test_full_complex_S_comparison_sees_imaginary_and_phase_errors():
    S = np.eye(4, dtype=complex)
    d1 = b3.s_difference(S, S + 2e-10j)                           # purely imaginary error
    assert d1.max() == pytest.approx(2e-10) and not b3.passes(d1.max(), 1e-10)
    S2 = S.copy(); S2[0, 0] = np.exp(1j * 3e-10)                  # phase-only error, |S| unchanged
    assert not b3.passes(b3.s_difference(S, S2).max(), 1e-10) and abs(S2[0, 0]) == pytest.approx(1.0, abs=1e-15)
    Z = np.zeros((4, 4), dtype=complex)                           # zero background: the difference is exactly representable
    assert b3.passes(b3.s_difference(Z, Z + 1e-10).max(), 1e-10)
    assert not b3.passes(b3.s_difference(Z, Z + 1e-10 * (1 + 1e-9)).max(), 1e-10)
    with pytest.raises(ValueError):
        b3.s_difference(np.eye(3), np.eye(4))


# ----------------------------------------------------------------------
# unitarity and reciprocity diagnostics (synthetic S)
# ----------------------------------------------------------------------
def _sym_unitary(n, seed=1):
    rng = np.random.default_rng(seed)
    Z = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    Q, _ = np.linalg.qr(Z)
    return Q @ Q.T                                               # unitary and symmetric


def test_unitarity_both_products():
    S = _sym_unitary(6)
    u = b3.unitarity(S)
    assert u["U"] < 1e-14 and set(u) == {"U_StS", "U_SSt", "U"}
    u2 = b3.unitarity(S * (1 + 1e-9))
    assert u2["U_StS"] == pytest.approx(2e-9, rel=0.2) and u2["U_SSt"] == pytest.approx(2e-9, rel=0.2)
    assert not b3.passes(u2["U"], 1e-10)
    S3 = S.copy(); S3[1, 4] += 3e-10                              # a single-entry violation is seen by both products
    assert b3.unitarity(S3)["U_StS"] > 1e-10 and b3.unitarity(S3)["U_SSt"] > 1e-10


def test_reciprocity_full_blocks_probability_and_report_only():
    n_open = 3
    S = _sym_unitary(2 * n_open, seed=2)
    q = b3.reciprocity(S, n_open)
    assert q["Q"] < 1e-14 and q["Q_prob"] < 1e-14
    assert q["same_dir_RL_report_only"] > 1e-3                    # same-direction symmetry is NOT implied
    assert q["Q"] == max(q["Q_full"], q["Q_rL"], q["Q_rR"], q["Q_t"], q["Q_prob"])
    for blk, key, sl in (("rL", "Q_rL", (slice(0, 3), slice(0, 3))), ("rR", "Q_rR", (slice(3, 6), slice(3, 6))),
                         ("t", "Q_t", (slice(3, 6), slice(0, 3)))):
        S2 = S.copy()
        S2[sl[0], sl[1]] = S2[sl[0], sl[1]] + 1e-6 * np.triu(np.ones((3, 3)), 1)   # break one block's symmetry
        q2 = b3.reciprocity(S2, n_open)
        assert q2[key] > 1e-7 and q2["Q_full"] > 1e-7
        assert not b3.passes(q2["Q"], 1e-8)
    # probability form catches a change of |t_LR| that matches no symmetric partner
    S3 = S.copy(); S3[0, 4] *= 1.001
    assert b3.reciprocity(S3, n_open)["Q_prob"] > 1e-4
    # index mapping: P_RL[m, n] (in left n, out right m) = P_LR[n, m] (in right m, out left n)
    p = sm.split_S(S, n_open)
    assert np.allclose(np.abs(p["t_RL"]) ** 2, (np.abs(p["t_LR"]) ** 2).T, atol=1e-14)


def test_T1_is_total_transmission_of_first_left_channel():
    S = _sym_unitary(6, seed=4)
    p = sm.split_S(S, 3)
    assert b3.T1(S, 3) == pytest.approx(float(np.sum(np.abs(S[3:, 0]) ** 2)), rel=1e-14)
    assert b3.T1(S, 3) == pytest.approx(1 - float(np.sum(np.abs(p["r_L"][:, 0]) ** 2)), abs=1e-14)


# ----------------------------------------------------------------------
# live small-case run of the B3 pipeline on one realization
# ----------------------------------------------------------------------
def test_pipeline_small_device_agreement_unitarity_reciprocity():
    cfg, pot = b3.case_potential("D")
    dev, geo, info, ref = b3.make_rgf_device("D", cfg, pot, SMALL)
    for e in b3.case_energies(cfg)[::2]:
        r = scattering_matrix(e["E"], dev)
        s = sr.sparse_scattering_matrix(e["E"], pot.numpy, geo["x_start"], SMALL["dx"], dev.n_slices, 1.0,
                                        SMALL["Ny"], SMALL["Nq"])
        assert r["n_open"] == s["n_open"] and r["S"].shape == s["S"].shape == (2 * r["n_open"],) * 2
        assert np.max(b3.s_difference(r["S"], s["S"])) <= 1e-10
        for S in (r["S"], s["S"]):
            assert b3.unitarity(S)["U"] <= 1e-10 and b3.reciprocity(S, r["n_open"])["Q"] <= 1e-8
        assert b3.T1(r["S"], r["n_open"]) == pytest.approx(b3.T1(s["S"], s["n_open"]), abs=1e-10)
        assert 0 < b3.T1(r["S"], r["n_open"]) < 1


def test_transverse_refinement_helpers_change_one_parameter_only():
    base = b3.BASE
    for name, v in b3.VARIANTS.items():
        changed = {k for k in base if v[k] != base[k]}
        expected = {"dx": {"dx"}, "Nq": {"Nq"}, "buffer": {"buffer"}, "Ny127": {"Ny", "Nq"},
                    "Ny127_Nq254": {"Ny", "Nq"}, "Ny255": {"Ny", "Nq"}}[name]
        assert changed == expected, name
    assert b3.VARIANTS["Ny127"]["Nq"] == 2 * 127 + 1          # default N_q = 2 N_y + 1 at each N_y
    assert b3.VARIANTS["Ny127_Nq254"]["Nq"] == b3.VARIANTS["Nq"]["Nq"]   # exploratory: same N_q as the Nq run


# ----------------------------------------------------------------------
# committed artifacts
# ----------------------------------------------------------------------
def _csv(name):
    with open(os.path.join(RES, name + ".csv")) as f:
        return list(csv.DictReader(f))


def _num(rows, k):
    return np.array([float(r[k]) for r in rows])


@pytest.mark.skipif(not HAVE, reason="B3 artifacts not generated")
def test_artifacts_complete_and_row_counts():
    for n in ("rgf_sparse_comparison unitarity reciprocity closed_channel_convergence dx_convergence "
              "quadrature_convergence buffer_convergence potential_identity").split():
        assert os.path.exists(os.path.join(RES, n + ".csv")), n
    assert os.path.exists(os.path.join(RES, "summary.json")) and os.path.exists(os.path.join(RES, "run_metadata.json"))
    comp = _csv("rgf_sparse_comparison")
    assert len(comp) == 20 and {r["case"] for r in comp} == set("ABCD")
    for c in "ABCD":
        assert len([r for r in comp if r["case"] == c]) == 5
    assert len(_csv("unitarity")) == 40 and len(_csv("reciprocity")) == 40
    assert len(_csv("dx_convergence")) == len(_csv("quadrature_convergence")) == len(_csv("buffer_convergence")) == 20
    closed = _csv("closed_channel_convergence")
    crit = [r for r in closed if r["role"] == "criterion"]
    assert len(crit) == 20 and len([r for r in closed if r["role"].startswith("exploratory")]) == 20
    cond = [r for r in closed if r["role"].startswith("conditional")]
    # N_y = 255 was run iff the zeta = 0.02 case failed at N_y = 127, and only for that case
    A_failed = any(r["passed"] == "False" for r in crit if r["case"] == "A")
    assert bool(cond) == A_failed and {r["case"] for r in cond} <= {"A"}
    assert len(cond) == (10 if A_failed else 0)


@pytest.mark.skipif(not HAVE, reason="B3 artifacts not generated")
def test_energies_dimensions_and_pass_flags_recomputed():
    with open(os.path.join(RES, "run_metadata.json")) as f:
        m = json.load(f)
    comp = _csv("rgf_sparse_comparison")
    for c in "ABCD":
        es = [e for e in m["energies"] if e["case"] == c]
        assert len(es) == 5
        cfg = b3.case_config(c)
        assert [float(r["E"]) for r in comp if r["case"] == c] == [e["E"] for e in es] == \
               [e["E"] for e in b3.case_energies(cfg)]
    for r in comp:
        assert int(r["S_dim"]) == 2 * int(r["n_open"]) and int(r["matrix_order"]) == int(r["n_slices"]) * 63
        assert (float(r["max_abs_dS"]) <= 1e-10) == (r["passed"] == "True")
        assert float(r["max_abs_dS"]) == max(float(r[k]) for k in
                                             ("max_abs_dS_rL", "max_abs_dS_tLR", "max_abs_dS_tRL", "max_abs_dS_rR"))
    for fn, key in (("unitarity", "U"), ("reciprocity", "Q")):
        lim = 1e-10 if key == "U" else 1e-8
        for r in _csv(fn):
            assert (float(r[key]) <= lim) == (r["passed"] == "True")
    for r in _csv("unitarity"):
        assert float(r["U"]) == max(float(r["U_StS"]), float(r["U_SSt"]))
    for fn in ("closed_channel_convergence", "dx_convergence", "quadrature_convergence", "buffer_convergence"):
        for r in _csv(fn):
            assert float(r["abs_delta_T1"]) == abs(float(r["T1_variant"]) - float(r["T1_reference"]))
            assert (float(r["abs_delta_T1"]) <= 2e-4) == (r["passed"] == "True")
            assert int(r["n_open_reference"]) >= 4
    for r in _csv("closed_channel_convergence"):
        assert int(r["closed_reference"]) == int(r["Ny_reference"]) - int(r["n_open_reference"]) >= 1
        assert int(r["closed_variant"]) == int(r["Ny_variant"]) - int(r["n_open_variant"]) >= 1


@pytest.mark.skipif(not HAVE, reason="B3 artifacts not generated")
def test_summary_agrees_with_csv():
    with open(SUMMARY) as f:
        s = json.load(f)
    comp, uni, rec = _csv("rgf_sparse_comparison"), _csv("unitarity"), _csv("reciprocity")
    assert s["limits"] == b3.LIMITS
    assert s["S_agreement_RGF_vs_sparse"]["value"] == _num(comp, "max_abs_dS").max()
    for meth in ("RGF", "sparse"):
        u = [r for r in uni if r["method"] == meth]
        assert s[f"unitarity_{meth}"]["value"] == _num(u, "U").max()
        q = [r for r in rec if r["method"] == meth]
        assert s[f"reciprocity_{meth}"]["full"]["value"] == _num(q, "Q_full").max()
        assert s[f"reciprocity_{meth}"]["probability"]["value"] == _num(q, "Q_prob").max()
        assert s[f"reciprocity_{meth}"]["report_only_same_direction"]["max_RL"] == _num(q, "same_dir_RL_report_only").max()
    for key, fn, role in (("closed_channel_Ny63_vs_Ny127", "closed_channel_convergence", "criterion"),
                          ("dx_1_512_vs_1_1024", "dx_convergence", "criterion"),
                          ("Nq_127_vs_254", "quadrature_convergence", "criterion"),
                          ("buffer_0.2_vs_0.4", "buffer_convergence", "criterion")):
        rows = [r for r in _csv(fn) if r["role"] == role]
        assert s[key]["value"] == _num(rows, "abs_delta_T1").max()
        assert s[key]["passed"] == bool(_num(rows, "abs_delta_T1").max() <= 2e-4)
        failed = [r for r in rows if r["passed"] == "False"]
        assert len(s[key]["rows_failed"]) == len(failed)
        for c in "ABCD":
            assert s[key]["per_case_max"][c] == _num([r for r in rows if r["case"] == c], "abs_delta_T1").max()
    ident = _csv("potential_identity")
    assert s["potential_identity"]["passed"] == all(r["bitwise_equal"] == "True" and r["coords_bitwise_equal"] == "True"
                                                    for r in ident)
    cp = s["criteria_passed"]
    assert s["all_B3_criteria_passed"] == all(cp.values())
    assert s["overall_B3_status"] == ("PASS" if all(cp.values()) else "FAIL (see criteria_passed)")
    for k in cp:                                                   # per-criterion flag equals value <= limit
        v = s[k]
        if "value" in v:
            assert cp[k] == (v["value"] <= v["limit"])


@pytest.mark.skipif(not HAVE, reason="B3 artifacts not generated")
def test_potential_identity_rows_and_hashes_reproduce():
    ident = _csv("potential_identity")
    for c in "ABCD":
        rows = [r for r in ident if r["case"] == c]
        assert {r["variant"] for r in rows} >= {"base", "dx", "Nq", "buffer", "Ny127", "Ny127_Nq254"}
        assert any("sparse" in r["evaluator"] for r in rows)
        base = [r for r in rows if r["variant"] == "base"]
        assert len({r["potential_hash"] for r in base}) == 1       # RGF device and sparse reference sampled identical arrays
        # recompute the base hash from a direct, independent call
        cfg = b3.case_config(c)
        geo = b3.geometry(cfg, 1 / 512, 0.2)
        _, _, ref = b3.direct_reference_samples(c, geo["x_start"], 1 / 512, geo["n_slices"], 1.0, 127)
        assert b3.sha(ref) == base[0]["potential_hash"] and json.loads(base[0]["shape"]) == list(ref.shape)
        assert float(base[0]["max_abs_diff"]) == 0.0
    assert len({r["potential_hash"] for r in ident if r["variant"] == "base"}) == 4


@pytest.mark.skipif(not HAVE, reason="B3 artifacts not generated")
def test_metadata_complete_and_source_hashes_current():
    with open(os.path.join(RES, "run_metadata.json")) as f:
        m = json.load(f)
    for k in ("generated_utc", "command", "source_commit", "working_tree_dirty", "python", "numpy", "scipy",
              "machine", "cases", "energies", "base_discretization", "variants", "sparse_reference",
              "rgf_implementation", "shared_conventions", "code_fingerprint", "code_files", "limits", "geometry"):
        assert k in m, k
    assert sorted(m["cases"]) == list("ABCD") and len(m["energies"]) == 20
    assert "splu" in m["sparse_reference"] and "no import of stationary" in m["sparse_reference"]
    from provenance import code_fingerprint
    fp = code_fingerprint("rgf_b3")
    assert set(m["code_files"]) == set(fp["code_files"])
    # solver and spec files must be the ones the run used (protected files unchanged since the run)
    for f in ("stationary.py", "smatrix_spec.py", "sparse_reference.py", "validate_rgf_b3.py", "potential.py", "disorder.py"):
        assert m["code_files"][f] == fp["code_files"][f], f"{f} changed after the artifacts were generated"
    if "post_run_tests" in m:
        assert all("exit_status" in r and "command" in r and "summary" in r for r in m["post_run_tests"]["runs"])


@pytest.mark.skipif(not HAVE, reason="B3 artifacts not generated")
def test_one_comparison_row_reproduces_live():
    r0 = [r for r in _csv("rgf_sparse_comparison") if r["case"] == "C"][2]
    cfg, pot = b3.case_potential("C")
    dev, geo, _, _ = b3.make_rgf_device("C", cfg, pot, b3.BASE)
    E = float(r0["E"])
    r = scattering_matrix(E, dev)
    s = sr.sparse_scattering_matrix(E, pot.numpy, geo["x_start"], b3.BASE["dx"], dev.n_slices, 1.0, 63, 127)
    assert abs(float(np.max(b3.s_difference(r["S"], s["S"]))) - float(r0["max_abs_dS"])) < 1e-13
    assert abs(b3.T1(r["S"], r["n_open"]) - float(r0["T1_RGF"])) < 1e-13
    assert abs(b3.T1(s["S"], s["n_open"]) - float(r0["T1_sparse"])) < 1e-13
