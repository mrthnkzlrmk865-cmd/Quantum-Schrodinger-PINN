"""Validation of the independent 1D ODE scattering reference (ode1d.py).

Criteria from docs/PHASE2_METHOD_SELECTION.md, B2 (unchanged):
  (a) sharp rectangular barrier: |T_ODE - T_formula| <= 1e-8
  (b) smooth barrier: rtol 1e-10 -> 1e-12 changes T by <= 1e-9;
      doubling the clean integration range changes T by <= 1e-10
Additional criteria, fixed before the first measurement:
  flux conservation |R + T - 1| <= 1e-9; complex amplitudes vs exact
  plane-wave matching <= 1e-8.
The ODE parameters are the specification defaults (DOP853,
rtol = atol = 1e-12); they are not tuned per test.
"""
import numpy as np
import pytest

from analytic import rectangular_barrier_matching, rectangular_barrier_transmission
from config import waveguide_config
from ode1d import (clean_interval, potential_barrier_callable, rectangular_barrier, shift_reference,
                   solve_scattering)

# 1D longitudinal energies spanned by the B4 packet window: kx in [4.09, 35.91]
EPS_PACKET = np.linspace(0.5 * 4.09**2, 0.5 * 35.91**2, 20)


def _near(V0):
    return [V0 * (1 - 1e-3), V0 * (1 - 1e-6), V0, V0 * (1 + 1e-6), V0 * (1 + 1e-3)]


# ----------------------------------------------------------------------
# (a) sharp rectangular barrier
# ----------------------------------------------------------------------
@pytest.mark.parametrize("V0,w", [(150.0, 0.1), (250.0, 0.1), (250.0, 0.3)])
def test_sharp_barrier_vs_formula_and_matching(V0, w):
    V, bp = rectangular_barrier(V0, 0.0, w)
    energies = np.concatenate([EPS_PACKET, _near(V0)])
    for eps in energies:
        res = solve_scattering(eps, V, 0.0, w, bp)
        T_formula = rectangular_barrier_transmission(np.array([eps]), V0, w)[0]
        t_m, r_m = rectangular_barrier_matching(eps, V0, w)
        assert abs(res["T"] - T_formula) <= 1e-8, (eps, res["T"], T_formula)
        assert abs(res["t"] - t_m) <= 1e-8 and abs(res["r"] - r_m) <= 1e-8, eps
        assert res["flux_error"] <= 1e-9, eps


def test_exact_threshold_uses_analytic_limit():
    """At eps = V0 the ODE is regular (phi'' = 0 inside the barrier). It is
    compared with the analytic limit T = 1/(1 + V0 w^2/2), not with the
    singular textbook expression."""
    for V0, w in [(150.0, 0.1), (250.0, 0.1), (250.0, 0.3)]:
        V, bp = rectangular_barrier(V0, 0.0, w)
        res = solve_scattering(V0, V, 0.0, w, bp)
        assert abs(res["T"] - 1.0 / (1.0 + V0 * w**2 / 2.0)) <= 1e-8


def test_tunnelling_regime_below_barrier():
    """Deep tunnelling (thick barrier, low energy): small T must still be
    accurate in absolute and relative terms."""
    V0, w = 250.0, 0.3
    V, bp = rectangular_barrier(V0, 0.0, w)
    for eps in (5.0, 20.0, 60.0):
        res = solve_scattering(eps, V, 0.0, w, bp)
        Tf = rectangular_barrier_transmission(np.array([eps]), V0, w)[0]
        assert Tf < 1e-3
        assert abs(res["T"] - Tf) <= 1e-8 and abs(res["T"] / Tf - 1) <= 1e-6


def test_sharp_barrier_self_convergence_and_lead_range():
    V0, w = 250.0, 0.1
    V, bp = rectangular_barrier(V0, 0.0, w)
    for eps in (60.0, 249.0, 251.0, 400.0):
        base = solve_scattering(eps, V, 0.0, w, bp)
        tight = solve_scattering(eps, V, 0.0, w, bp, rtol=1e-13, atol=1e-13)
        loose = solve_scattering(eps, V, 0.0, w, bp, rtol=1e-10, atol=1e-10)
        assert abs(tight["T"] - base["T"]) <= 1e-9 and abs(loose["T"] - base["T"]) <= 1e-9
        wide = solve_scattering(eps, V, -0.5, w + 0.5, bp)
        assert abs(wide["T"] - base["T"]) <= 1e-10
        back = shift_reference(wide, 0.0, w)           # same reference planes as base
        assert abs(back["t"] - base["t"]) <= 1e-9 and abs(back["r"] - base["r"]) <= 1e-9


def test_clean_lead_is_transparent():
    V = lambda x: 0.0
    for eps in (1.0, 200.0, 600.0):
        res = solve_scattering(eps, V, 0.0, 0.7)
        assert abs(res["T"] - 1) <= 1e-12 and res["R"] <= 1e-20
        assert abs(res["t"] - np.exp(1j * res["k"] * 0.7)) <= 1e-10


def test_input_validation():
    V, bp = rectangular_barrier(250.0, 0.0, 0.1)
    with pytest.raises(ValueError):
        solve_scattering(0.0, V, 0.0, 0.1, bp)
    with pytest.raises(ValueError):
        solve_scattering(10.0, V, 0.1, 0.0, bp)


# ----------------------------------------------------------------------
# (b) smooth barrier of the existing Potential class
# ----------------------------------------------------------------------
def _smooth_setup(Vb):
    cfg = waveguide_config(barrier_height=Vb)
    V = potential_barrier_callable(cfg)
    x1 = cfg.barrier_center - cfg.barrier_width / 2
    x2 = cfg.barrier_center + cfg.barrier_width / 2
    a, b = clean_interval(V, x1, x2, tol=1e-12)
    return cfg, V, x1, x2, a, b


def test_smooth_barrier_uses_the_potential_callable():
    from potential import Potential
    cfg, V, x1, x2, a, b = _smooth_setup(250.0)
    xs = np.linspace(a, b, 101)
    ref = Potential(cfg.replace(disorder_strength=0.0)).barrier_numpy(xs)
    assert np.array_equal(np.array([V(x) for x in xs]), ref)    # bitwise identical
    assert cfg.barrier_edge > 0 and abs(V(a)) <= 1e-12 and abs(V(b)) <= 1e-12


@pytest.mark.parametrize("Vb", [150.0, 250.0])
def test_smooth_barrier_criteria(Vb):
    cfg, V, x1, x2, a, b = _smooth_setup(Vb)
    for eps in np.concatenate([EPS_PACKET[::2], _near(Vb)]):
        base = solve_scattering(eps, V, a, b)
        assert base["flux_error"] <= 1e-9
        loose = solve_scattering(eps, V, a, b, rtol=1e-10, atol=1e-10)
        assert abs(loose["T"] - base["T"]) <= 1e-9, eps                  # B2(b) tolerance criterion
        wide = solve_scattering(eps, V, a - (x1 - a), b + (b - x2))       # clean distance doubled
        assert abs(wide["T"] - base["T"]) <= 1e-10, eps                  # B2(b) range criterion


def test_smooth_barrier_approaches_sharp_barrier_as_edge_shrinks():
    """Consistency (informative): as the tanh edge width e -> 0 the smooth
    barrier tends to the sharp one, and so does T. The difference must
    decrease monotonically with e."""
    Vb, w, eps = 250.0, 0.1, 150.0
    Tsharp = rectangular_barrier_transmission(np.array([eps]), Vb, w)[0]
    diffs = []
    for e in (0.01, 0.003, 0.001):
        cfg = waveguide_config(barrier_height=Vb, barrier_edge=e)
        V = potential_barrier_callable(cfg)
        x1 = cfg.barrier_center - w / 2
        x2 = cfg.barrier_center + w / 2
        a, b = clean_interval(V, x1, x2, tol=1e-12, step=0.002)
        diffs.append(abs(solve_scattering(eps, V, a, b)["T"] - Tsharp))
    assert diffs[0] > diffs[1] > diffs[2]
