r"""
Independent 1D stationary scattering reference (validation benchmark B2,
docs/PHASE2_METHOD_SELECTION.md).

Solves  -1/2 phi''(x) + V(x) phi(x) = eps phi(x)   (hbar = m = 1, eps > 0)
for a real potential V that vanishes in the clean leads x <= x_left and
x >= x_right. It integrates the ODE directly with scipy.integrate.solve_ivp:
no lattice, no Green functions, no code shared with the time-dependent or
the planned stationary (RGF) solver.

Scattering boundary conditions (incidence from the left, k = sqrt(2 eps)):
    x <= x_left :  phi = e^{ik(x - x_left)} + r e^{-ik(x - x_left)}
    x >= x_right:  phi = t e^{ik(x - x_right)}
The amplitudes are referred to the planes x_left (incident and reflected)
and x_right (transmitted). Both leads have V = 0, so the flux-normalized
probabilities are T = |t|^2, R = |r|^2, and R + T = 1 for real V.
Moving a reference plane multiplies t or r by a unimodular phase
(shift_reference); T and R do not change.

Method (backward shooting, which is exact for a linear ODE):
  1. Start at x_right with the pure outgoing solution phi = 1,
     phi' = i k (i.e. a provisional t' = 1).
  2. Integrate phi'' = 2 (V(x) - eps) phi from x_right to x_left as the
     first-order system (phi, phi'), segment by segment between the
     breakpoints (discontinuities of V), carrying phi and phi' across each
     breakpoint (both are continuous for a finite step in V).
  3. At x_left decompose into plane waves:
         A = (phi + phi'/(ik))/2   (incident),  B = (phi - phi'/(ik))/2   (reflected).
  4. Normalize to unit incidence: t = 1/A, r = B/A.

Within a segment, V is evaluated at x clipped into the segment's interior,
so the adaptive integrator never samples V on the wrong side of a
discontinuity (this matters only for sharp potentials).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

DEFAULT_RTOL = 1e-12   # docs/PHASE2_METHOD_SELECTION.md, B2 (rtol = atol = 1e-12)
DEFAULT_ATOL = 1e-12


def solve_scattering(eps, V, x_left, x_right, breakpoints=(), rtol=DEFAULT_RTOL, atol=DEFAULT_ATOL,
                     method="DOP853"):
    """Transmission/reflection of a left-incident plane wave.

    V           callable V(x) -> float (scalar x), real; V = 0 outside [x_left, x_right]
    breakpoints positions of discontinuities of V inside (x_left, x_right)
    Returns dict(t, r, T, R, k, flux_error=|T + R - 1|, nfev, ...)."""
    eps = float(eps)
    if eps <= 0:
        raise ValueError("eps must be > 0 (propagating incident wave)")
    if not x_right > x_left:
        raise ValueError("need x_right > x_left")
    k = np.sqrt(2.0 * eps)
    pts = sorted({float(x_left), float(x_right)} | {float(b) for b in breakpoints if x_left < b < x_right})
    y = np.array([1.0 + 0j, 1j * k])          # outgoing wave at x_right
    nfev = 0
    for xb, xa in zip(pts[::-1][:-1], pts[::-1][1:]):    # segments right -> left: [xa, xb]
        span = xb - xa
        lo, hi = xa + 1e-12 * span, xb - 1e-12 * span

        def rhs(x, u, lo=lo, hi=hi):
            v = float(V(min(max(x, lo), hi)))
            return [u[1], 2.0 * (v - eps) * u[0]]

        sol = solve_ivp(rhs, (xb, xa), y, method=method, rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(f"ODE integration failed on [{xa}, {xb}]: {sol.message}")
        y = sol.y[:, -1]
        nfev += sol.nfev
    phi, dphi = y
    A = 0.5 * (phi + dphi / (1j * k))
    B = 0.5 * (phi - dphi / (1j * k))
    t, r = 1.0 / A, B / A
    T, R = abs(t) ** 2, abs(r) ** 2
    return dict(t=t, r=r, T=T, R=R, k=k, flux_error=abs(T + R - 1.0), nfev=nfev,
                x_left=float(x_left), x_right=float(x_right), rtol=rtol, atol=atol, method=method)


def shift_reference(res, new_left, new_right):
    """Amplitudes referred to other reference planes (V must vanish between
    the old and new planes). Incident and reflected waves are referred to
    new_left, the transmitted wave to new_right:
        t' = t exp(ik[(new_right - x_right) + (x_left - new_left)]),
        r' = r exp(2ik (x_left - new_left))."""
    k = res["k"]
    t = res["t"] * np.exp(1j * k * ((new_right - res["x_right"]) + (res["x_left"] - new_left)))
    r = res["r"] * np.exp(2j * k * (res["x_left"] - new_left))
    return dict(res, t=t, r=r, x_left=new_left, x_right=new_right)


# ----------------------------------------------------------------------
# Potentials
# ----------------------------------------------------------------------
def rectangular_barrier(V0, x1, x2):
    """Sharp barrier V = V0 on (x1, x2), 0 elsewhere. Returns (V, breakpoints)."""
    def V(x):
        return V0 if x1 < x < x2 else 0.0
    return V, (x1, x2)


def clean_interval(V, x_lo, x_hi, tol, step=0.01, max_iter=100000):
    """Expand [x_lo, x_hi] outward in steps of `step` until |V| <= tol at
    both ends and at the next step beyond them. For a smooth barrier this
    gives an interval outside which the neglected tails are below tol.

    Raises ValueError for invalid inputs, and RuntimeError if no clean end
    point is found within max_iter steps on either side (e.g. a potential
    that does not decay to |V| <= tol)."""
    x_lo, x_hi, tol, step = float(x_lo), float(x_hi), float(tol), float(step)
    if not (np.isfinite(x_lo) and np.isfinite(x_hi)) or x_hi < x_lo:
        raise ValueError(f"need finite x_lo <= x_hi, got [{x_lo}, {x_hi}]")
    if not (np.isfinite(tol) and tol >= 0):
        raise ValueError(f"tol must be finite and >= 0, got {tol}")
    if not (np.isfinite(step) and step > 0):
        raise ValueError(f"step must be finite and > 0, got {step}")
    if int(max_iter) != max_iter or max_iter < 1:
        raise ValueError(f"max_iter must be a positive integer, got {max_iter}")
    a, b = x_lo, x_hi
    for _ in range(int(max_iter)):
        if abs(V(a)) <= tol and abs(V(a - step)) <= tol:
            break
        a -= step
    else:
        raise RuntimeError(f"clean_interval: |V| > tol = {tol:g} still at x = {a:g} after {max_iter} steps "
                           f"of {step:g} to the left; the potential does not decay to tol")
    for _ in range(int(max_iter)):
        if abs(V(b)) <= tol and abs(V(b + step)) <= tol:
            break
        b += step
    else:
        raise RuntimeError(f"clean_interval: |V| > tol = {tol:g} still at x = {b:g} after {max_iter} steps "
                           f"of {step:g} to the right; the potential does not decay to tol")
    return a, b


def potential_barrier_callable(cfg):
    """The barrier of the existing Potential class (the same callable the
    planned RGF solver must sample; B2). Disorder is not included."""
    from potential import Potential
    pot = Potential(cfg.replace(disorder_strength=0.0))
    return lambda x: float(pot.barrier_numpy(np.asarray(x, dtype=float)))
