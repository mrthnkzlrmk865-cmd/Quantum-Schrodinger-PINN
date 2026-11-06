"""Benchmark B2: RGF vs independent 1D ODE for the y-uniform smooth barrier
(validate_rgf_b2.py, docs/RGF_B2_VALIDATION.md).

Machinery tests run live on small cases; artifact tests check the committed
results/rgf_b2_validation files (skipped if absent). Thresholds are the
predeclared ones of docs/PHASE2_METHOD_SELECTION.md §3 B2 and are only
checked for preservation here, never changed.
"""
import csv
import json
import os

import numpy as np
import pytest

import validate_rgf_b2 as b2
from config import waveguide_config
from ode1d import potential_barrier_callable
from potential import Potential

RES = "results/rgf_b2_validation"
CSV = os.path.join(RES, "rgf_ode_comparison.csv")
SUMMARY = os.path.join(RES, "summary.json")
META = os.path.join(RES, "run_metadata.json")
HAVE = os.path.exists(CSV) and os.path.exists(SUMMARY)


# ----------------------------------------------------------------------
# thresholds and parameters are the predeclared ones
# ----------------------------------------------------------------------
def test_predeclared_thresholds_and_parameters_preserved():
    assert b2.LIMITS == dict(dT=2e-4, mixing=1e-12, U=1e-10, Q=1e-8)
    assert b2.PROVISIONAL == dict(richardson_dT=2e-5, order_window=(1.8, 2.2))
    assert b2.ODE_LIMITS == dict(sharp_T_vs_formula=1e-8, flux=1e-9, smooth_rtol_change=1e-9,
                                 smooth_range_change=1e-10)
    assert b2.TD_LIMITS == dict(interval=1e-3, P_in=1e-4, quad_halving=2e-4)
    assert b2.VB_LIST == (150.0, 250.0) and b2.DXS == (1 / 256, 1 / 512, 1 / 1024) and b2.BINDING_DX == 1 / 512
    assert b2.NY == 63 and b2.THRESHOLD_REL == 1e-3
    assert np.array_equal(b2.KX, np.linspace(4.09, 35.91, 60))
    from ode1d import DEFAULT_ATOL, DEFAULT_RTOL
    assert DEFAULT_RTOL == 1e-12 and DEFAULT_ATOL == 1e-12
    cfg = waveguide_config(barrier_height=250.0)
    assert (cfg.barrier_width, cfg.barrier_edge, cfg.barrier_center, cfg.disorder_strength) == (0.1, 0.01, 2.4, 0.0)


# ----------------------------------------------------------------------
# energy convention and threshold rule
# ----------------------------------------------------------------------
def test_total_energy_conversion():
    Ep = b2.e_perp()
    assert Ep[0] == 0.5 * np.pi ** 2
    grid = b2.energy_grid(Ep)
    assert len(grid) == 60
    for g, kx in zip(grid, b2.KX):
        assert g["kx"] == kx and g["eps_original"] == 0.5 * kx ** 2
        assert g["E_total"] == g["eps"] + Ep[0]                    # E = eps + E_perp_1 (not E = eps)
        if not g["adjusted"]:
            assert g["eps"] == g["eps_original"]


def test_threshold_rule_deterministic_and_grid_preserved():
    Ep = b2.e_perp()
    grid = b2.energy_grid(Ep)
    # the 60 B2 energies are all >= 1e-3 E from every threshold: none is moved
    assert not any(g["adjusted"] for g in grid)
    assert min(g["rel_dist"] for g in grid) >= 1e-3
    assert grid == b2.energy_grid(Ep)                              # deterministic
    # a point too close to a threshold is moved to E_n (1 +/- 2e-3) on its own side
    for n in (2, 5, 9):
        En = Ep[n - 1]
        for E, side in ((En * (1 + 3e-4), 1), (En * (1 - 3e-4), -1), (En, 1)):
            a = b2.adjust_energy(E, Ep)
            assert a["adjusted"] and a["nearest_threshold_n"] == n
            assert a["E"] == En * (1 + side * 2e-3)
            assert np.all(np.abs(a["E"] - Ep) >= 1e-3 * a["E"])
            assert a["n_open"] == int(np.sum(Ep < E)) + (1 if E == En else 0)
            assert b2.adjust_energy(E, Ep) == a
    far = b2.adjust_energy(100.0, Ep)
    assert not far["adjusted"] and far["E"] == 100.0


# ----------------------------------------------------------------------
# same potential, bitwise
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Vb", b2.VB_LIST)
def test_rgf_samples_are_the_potential_barrier_callable_bitwise(Vb):
    cfg = waveguide_config(barrier_height=Vb)
    pot = Potential(cfg)
    dev, chk = b2.rgf_device(pot, 1 / 256)
    assert chk["samples_bitwise_equal_barrier_numpy"]
    assert dev.n_slices == 513 and dev.xs[0] == 1.5 and dev.xs[-1] == 3.5
    # the ODE callable is the same function: identical values at the RGF slices
    V = potential_barrier_callable(cfg)
    assert np.array_equal(np.array([V(x) for x in dev.xs]), pot.barrier_numpy(dev.xs))
    # clean leads: the barrier tail on the first lead slices is negligible
    assert chk["V_lead_max"] < 1e-30 and max(chk["V_at_device_ends"]) < 1e-30


def test_recording_wrapper_detects_a_different_potential():
    cfg = waveguide_config(barrier_height=150.0)

    class Shifted(Potential):
        def barrier_numpy(self, X):
            return super().barrier_numpy(X) * (1 + 1e-15)

    pot = Potential(cfg)
    rec = b2.RecordingPotential(Shifted(cfg))
    from stationary import build_device
    dev = build_device(rec, 1.5, 3.5, 1 / 256, 1.0, 63)
    Vs = rec.calls[0][1]
    assert not np.array_equal(Vs[:, 0], pot.barrier_numpy(dev.xs))


@pytest.mark.parametrize("Vb", b2.VB_LIST)
def test_y_uniform_potential_blocks_are_diagonal(Vb):
    pot = Potential(waveguide_config(barrier_height=Vb))
    dev, chk = b2.rgf_device(pot, 1 / 256)
    Vmax = chk["max_V"]
    assert Vmax > 0.99 * Vb
    assert chk["max_offdiag_V_block"] <= 1e-13 * Vmax          # V_j = V(x_j) I up to rounding
    assert chk["max_diag_V_block_minus_Vx"] <= 1e-13 * Vmax
    assert np.array_equal(dev.V, np.transpose(dev.V, (0, 2, 1)))


# ----------------------------------------------------------------------
# RGF: decoupling, unitarity, reciprocity, agreement (live, coarse grid)
# ----------------------------------------------------------------------
@pytest.mark.parametrize("Vb,i", [(150.0, 0), (250.0, 18), (250.0, 59)])
def test_rgf_point_mixing_unitarity_reciprocity_and_ode_agreement(Vb, i):
    cfg = waveguide_config(barrier_height=Vb)
    pot = Potential(cfg)
    dev, _ = b2.rgf_device(pot, 1 / 256)
    g = b2.energy_grid()[i]
    r = b2.rgf_point(g["E_total"], dev)
    assert r["channel_mixing"] <= b2.LIMITS["mixing"]
    assert r["channel_mixing_all_channels"] <= b2.LIMITS["mixing"]
    assert r["unitarity_residual"] <= b2.LIMITS["U"]
    assert r["reciprocity_residual"] <= b2.LIMITS["Q"]
    assert abs(r["T_RGF"] - r["T_RGF_11"]) <= 1e-12               # all transmission stays in channel 1
    assert abs(r["T_RGF"] - r["T_RGF_R_to_L_ch1"]) <= 1e-10       # both incidence directions
    o = b2.ode_point(g["eps"], cfg, b2.ode_setup(cfg))
    # at dx = 1/256 (not the binding grid) agreement to the lattice-dispersion level only
    assert abs(r["T_RGF"] - o["T_ODE"]) < 1e-2
    assert o["ODE_flux_error"] <= b2.ODE_LIMITS["flux"]


# ----------------------------------------------------------------------
# convergence order and Richardson
# ----------------------------------------------------------------------
def test_observed_order():
    for p in (1.0, 2.0, 2.7):
        h = np.array([1 / 256, 1 / 512, 1 / 1024])
        e = 3.0 * h ** p
        assert b2.observed_order(e[0], e[1]) == pytest.approx(p, rel=1e-12)
        assert b2.observed_order(-e[1], -e[2]) == pytest.approx(p, rel=1e-12)
    assert np.isnan(b2.observed_order(1e-3, 1e-10))               # below the reference floor
    assert np.isnan(b2.observed_order(1e-12, 1e-3))
    assert b2.observed_order(1e-4, 1e-4) == 0.0
    assert b2.in_window(2.0) and b2.in_window(1.8) and b2.in_window(2.2)
    assert not b2.in_window(1.79) and not b2.in_window(float("nan"))
    T_exact, c = 0.4, 5.0
    Tc, Tf = T_exact + c * (1 / 512) ** 2, T_exact + c * (1 / 1024) ** 2
    assert b2.richardson(Tc, Tf) == pytest.approx(T_exact, abs=1e-15)


def test_energy_average_weights_and_interval():
    cfg = waveguide_config()
    k = b2.kx_window(cfg.kx0, cfg.sigma_x)
    assert k[0] == pytest.approx(4.0901, abs=1e-4) and k[-1] == pytest.approx(35.9099, abs=1e-4)
    assert k.size == 1274 and abs(np.diff(k).mean() - 0.025) < 1e-4
    kf = b2.kx_window(cfg.kx0, cfg.sigma_x, b2.DKX / 2)
    assert np.allclose(kf[::2], k, atol=1e-12, rtol=0)
    av, normP = b2.energy_average(k, np.ones_like(k), cfg.kx0, cfg.sigma_x)
    # the +/- 4.5 sigma_k window (spec B4) holds 1 - erfc(4.5/sqrt 2) of the weight
    from scipy.special import erfc
    assert normP == pytest.approx(1 - erfc(4.5 / np.sqrt(2)), abs=1e-9)
    assert av == normP                                            # clean guide: T_pkt = int P (= 1 - 6.8e-6)
    av2, _ = b2.energy_average(k, k, cfg.kx0, cfg.sigma_x)
    assert av2 / normP == pytest.approx(cfg.kx0, rel=1e-9)        # mean of P is kx0
    assert b2.interval_distance(0.5, 0.4, 0.6) == 0.0
    assert b2.interval_distance(0.3, 0.4, 0.6) == pytest.approx(0.1)
    assert b2.interval_distance(0.7, 0.4, 0.6) == pytest.approx(0.1)


# ----------------------------------------------------------------------
# committed artifacts
# ----------------------------------------------------------------------
def _rows():
    with open(CSV) as f:
        return list(csv.DictReader(f))


@pytest.mark.skipif(not HAVE, reason="B2 artifacts not generated")
def test_csv_complete_and_columns():
    rows = _rows()
    assert len(rows) == 2 * 60 * 3
    need = ("Vb width edge kx epsilon E_total dx N_y N_q number_of_slices number_of_open_channels T_RGF R_RGF "
            "T_ODE R_ODE absolute_T_error relative_T_error observed_convergence_order unitarity_residual "
            "reciprocity_residual channel_mixing min_rcond ODE_flux_error ODE_nfev").split()
    for k in need:
        assert k in rows[0]
    for Vb in b2.VB_LIST:
        for dx in b2.DXS:
            rr = [r for r in rows if float(r["Vb"]) == Vb and float(r["dx"]) == dx]
            assert [float(r["kx"]) for r in rr] == list(b2.KX)
            for r in rr:
                assert float(r["E_total"]) == float(r["epsilon"]) + 0.5 * np.pi ** 2
                assert float(r["absolute_T_error"]) == abs(float(r["T_RGF"]) - float(r["T_ODE"]))
                assert int(r["N_y"]) == 63 and int(r["N_q"]) >= 126
                assert int(r["number_of_slices"]) == round(2 / dx) + 1


@pytest.mark.skipif(not HAVE, reason="B2 artifacts not generated")
def test_summary_agrees_with_csv():
    rows = _rows()
    with open(SUMMARY) as f:
        s = json.load(f)
    for Vb in b2.VB_LIST:
        pb = s["per_barrier"][str(Vb)]
        for dx in b2.DXS:
            rr = [r for r in rows if float(r["Vb"]) == Vb and float(r["dx"]) == dx]
            err = np.array([float(r["absolute_T_error"]) for r in rr])
            v = pb["per_dx"][f"1/{round(1 / dx)}"]
            assert v["max_abs_T_error"] == err.max()
            assert np.isclose(v["rms_abs_T_error"], np.sqrt(np.mean(err ** 2)), rtol=1e-12)
            assert v["argmax"]["kx"] == float(rr[int(np.argmax(err))]["kx"])
            assert v["max_unitarity"] == max(float(r["unitarity_residual"]) for r in rr)
            assert v["max_channel_mixing"] == max(float(r["channel_mixing"]) for r in rr)
        assert pb["binding_dx_max_abs_T_error"] == pb["per_dx"]["1/512"]["max_abs_T_error"]
        # per-energy orders in the CSV are recomputable from the errors
        r512 = [r for r in rows if float(r["Vb"]) == Vb and float(r["dx"]) == 1 / 512]
        r256 = [r for r in rows if float(r["Vb"]) == Vb and float(r["dx"]) == 1 / 256]
        for a, b in zip(r256, r512):
            p = b2.observed_order(float(a["signed_T_error"]), float(b["signed_T_error"]))
            q = float(b["observed_convergence_order"])
            assert (np.isnan(p) and np.isnan(q)) or p == pytest.approx(q, rel=1e-12)
    c = s["criteria"]
    assert c["A_binding_dT_dx512"]["limit"] == 2e-4 and c["C_channel_mixing"]["limit"] == 1e-12
    assert c["D_unitarity"]["limit"] == 1e-10 and c["E_reciprocity"]["limit"] == 1e-8
    assert c["A_binding_dT_dx512"]["pass"] == (max(c["A_binding_dT_dx512"]["value"].values()) <= 2e-4)
    for k in ("C_channel_mixing", "D_unitarity", "E_reciprocity"):
        col = {"C_channel_mixing": "channel_mixing", "D_unitarity": "unitarity_residual",
               "E_reciprocity": "reciprocity_residual"}[k]
        assert c[k]["value"] == max(float(r[col]) for r in rows)
        assert c[k]["pass"] == (c[k]["value"] <= c[k]["limit"])
    assert s["potential_checks"]["all_samples_bitwise_equal"]
    if s["overall_B2_status"] == "PASS":
        assert s["rgf_vs_ode_all_binding_pass"] and s["td_energy_averaged_comparison"]["all_pass"]


@pytest.mark.skipif(not HAVE, reason="B2 artifacts not generated")
def test_metadata_complete():
    with open(META) as f:
        m = json.load(f)
    for k in ("generated_utc", "command", "source_commit", "working_tree_dirty", "python", "numpy", "scipy",
              "machine", "parameters", "limits", "code_fingerprint"):
        assert k in m
    assert m["limits"]["binding"] == b2.LIMITS


@pytest.mark.skipif(not HAVE, reason="B2 artifacts not generated")
def test_binding_row_reproduces_live():
    """The worst dx = 1/512 row of V_b = 250 recomputed from scratch."""
    rows = [r for r in _rows() if float(r["Vb"]) == 250.0 and float(r["dx"]) == 1 / 512]
    w = max(rows, key=lambda r: float(r["absolute_T_error"]))
    cfg = waveguide_config(barrier_height=250.0)
    dev, _ = b2.rgf_device(Potential(cfg), 1 / 512)
    r = b2.rgf_point(float(w["E_total"]), dev)
    o = b2.ode_point(float(w["epsilon"]), cfg, b2.ode_setup(cfg))
    assert abs(r["T_RGF"] - float(w["T_RGF"])) <= 1e-13
    assert o["T_ODE"] == float(w["T_ODE"])
