"""Analytic rectangular-barrier transmission, including the threshold eps = V0."""
import numpy as np
import pytest

from analytic import rectangular_barrier_transmission

V0, W = 250.0, 0.1
T_LIMIT = 1.0 / (1.0 + V0 * W**2 / 2.0)


def T(e):
    return rectangular_barrier_transmission(np.atleast_1d(np.asarray(e, float)), V0, W)


def _textbook(e):
    """The textbook piecewise expressions (singular at eps = V0)."""
    if e > V0:
        q = np.sqrt(2 * (e - V0))
        return 1 / (1 + V0**2 * np.sin(q * W) ** 2 / (4 * e * (e - V0)))
    k = np.sqrt(2 * (V0 - e))
    return 1 / (1 + V0**2 * np.sinh(k * W) ** 2 / (4 * e * (V0 - e)))


def test_exactly_at_threshold_uses_continuous_limit():
    assert T(V0)[0] == pytest.approx(T_LIMIT, abs=1e-15)


@pytest.mark.parametrize("delta", [1e-12, 1e-9, 1e-6, 1e-3])
def test_immediately_below_and_above_threshold(delta):
    lo, hi = T([V0 - delta, V0 + delta])
    # continuous, approaching the limit linearly from both sides
    assert lo < T_LIMIT < hi
    slope = (hi - lo) / (2 * delta)
    assert abs(lo - T_LIMIT) <= 2 * abs(slope) * delta + 1e-14
    assert abs(hi - T_LIMIT) <= 2 * abs(slope) * delta + 1e-14


def test_derivative_finite_and_consistent_across_threshold():
    s = [(T(V0 + d)[0] - T(V0 - d)[0]) / (2 * d) for d in (1e-2, 1e-3, 1e-4)]
    assert np.allclose(s, s[-1], rtol=1e-3)


def test_identical_to_textbook_form_away_from_threshold():
    es = np.concatenate([np.linspace(5.0, 600.0, 400), [V0 * (1 - 1e-3), V0 * (1 + 1e-3)]])
    new = T(es)
    old = np.array([_textbook(e) for e in es])
    assert np.max(np.abs(new - old)) < 1e-14


def test_limits_and_bounds():
    es = np.linspace(1.0, 2000.0, 1000)
    t = T(es)
    assert np.all((t > 0) & (t <= 1))
    assert T(1e-6)[0] < 1e-6                     # deep tunnelling, low energy
    assert T(1e6)[0] > 0.999                     # far above the barrier


def test_rejects_non_propagating_energy():
    with pytest.raises(ValueError):
        T(0.0)
