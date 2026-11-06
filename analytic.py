r"""
Exact solutions used to validate the numerical solvers.

1. Free Gaussian packet in unbounded 2D space (hbar = m = 1). For
   psi0 = (pi sx sy)^{-1/2} exp(-xi^2/(2 sx^2) - eta^2/(2 sy^2)) e^{i(kx x + ky y)},
   xi = x - x0, eta = y - y0, the solution factorizes into 1D packets

       phi(xi, t) = (pi s^2)^{-1/4} alpha^{-1/2}
                    exp( [-xi^2/(2 s^2) + i k xi - i k^2 t/2] / alpha ) e^{i k x0},
       alpha = 1 + i t / s^2.

   In a box with Dirichlet walls this is accurate as long as the packet
   has negligible amplitude at the walls (checked by the caller).

2. Dirichlet box eigenmode  psi = (2/sqrt(Lx Ly)) sin(m pi x/Lx) sin(n pi y/Ly) e^{-iEt},
   E = ((m pi/Lx)^2 + (n pi/Ly)^2)/2, valid for V = 0 (or V = const, with
   E -> E + V).
"""

import numpy as np


def free_gaussian_1d(x, t, x0, s, k):
    alpha = 1.0 + 1j * t / s**2
    xi = x - x0
    return ((np.pi * s**2) ** -0.25 * alpha**-0.5
            * np.exp((-xi**2 / (2 * s**2) + 1j * k * xi - 0.5j * k**2 * t) / alpha)
            * np.exp(1j * k * x0))


def free_gaussian_2d(cfg, X, Y, t):
    assert cfg.packet == "gaussian"
    return (free_gaussian_1d(X, t, cfg.x0, cfg.sigma_x, cfg.kx0)
            * free_gaussian_1d(Y, t, cfg.y0, cfg.sigma_y, cfg.ky0))


def box_eigenmode(X, Y, t, Lx, Ly, m=1, n=1, V0=0.0):
    E = 0.5 * ((m * np.pi / Lx) ** 2 + (n * np.pi / Ly) ** 2) + V0
    return (2.0 / np.sqrt(Lx * Ly) * np.sin(m * np.pi * X / Lx)
            * np.sin(n * np.pi * Y / Ly) * np.exp(-1j * E * t))


def _sinc_like(z, hyperbolic):
    """sin(z)/z (hyperbolic=False) or sinh(z)/z (hyperbolic=True), z >= 0,
    evaluated without cancellation near z = 0 (value 1 at z = 0)."""
    z = np.asarray(z, dtype=float)
    out = np.ones_like(z)
    small = z < 1e-4
    zs = z[small]
    sgn = 1.0 if hyperbolic else -1.0
    out[small] = 1.0 + sgn * zs**2 / 6.0 + zs**4 / 120.0
    zl = z[~small]
    out[~small] = (np.sinh(zl) if hyperbolic else np.sin(zl)) / zl
    return out


def rectangular_barrier_transmission(eps, V0, w):
    """Exact transmission probability of a sharp rectangular barrier of
    height V0 > 0 and width w for a particle of energy eps > 0 (hbar = m = 1):

        eps > V0:  T = [1 + V0^2 sin^2(q w)  / (4 eps (eps - V0))]^-1,  q = sqrt(2 (eps - V0))
        eps < V0:  T = [1 + V0^2 sinh^2(k w) / (4 eps (V0 - eps))]^-1,  k = sqrt(2 (V0 - eps))
        eps = V0:  T = [1 + V0 w^2 / 2]^-1                               (continuous limit)

    Evaluated in the algebraically identical form
        T = [1 + V0^2 w^2 s^2 / (2 eps)]^-1,
        s = sin(q w)/(q w)  (eps > V0),  sinh(k w)/(k w)  (eps < V0),  1  (eps = V0),
    which follows from sin^2(qw) = q^2 w^2 (sin(qw)/(qw))^2 and q^2 = 2(eps - V0).
    This form has no singular intermediate expression, so it is continuous
    and accurate through eps = V0.

    Standard textbook result. In this project it is checked against an
    independent exact plane-wave matching calculation
    (tests/test_smatrix_spec.py::test_rectangular_barrier_formula_vs_matching,
    tests/test_analytic.py)."""
    eps = np.asarray(eps, dtype=float)
    if np.any(eps <= 0):
        raise ValueError("eps must be > 0 (propagating incident wave)")
    d = eps - V0
    z = np.sqrt(2.0 * np.abs(d)) * w
    # explicit masks: each branch is evaluated only where it applies (an
    # unused sinh branch would overflow for large above-barrier z > ~710)
    s = np.ones_like(z)                      # eps == V0: s = 1 (continuous limit)
    above, below = d > 0, d < 0
    s[above] = _sinc_like(z[above], hyperbolic=False)
    s[below] = _sinc_like(z[below], hyperbolic=True)
    return 1.0 / (1.0 + V0**2 * w**2 * s**2 / (2.0 * eps))


def rectangular_barrier_matching(eps, V0, w):
    """Exact complex amplitudes (t, r) for a sharp rectangular barrier on
    [0, w], by matching psi and psi' at x = 0 and x = w (continuum plane
    waves, hbar = m = 1). Convention (same as ode1d.solve_scattering with
    x_left = 0, x_right = w):
        x <= 0: e^{ikx} + r e^{-ikx};   x >= w: t e^{ik(x - w)}.
    Inside: A e^{iqx} + B e^{-iqx} with q = sqrt(2(eps - V0)) (imaginary
    below the barrier); exactly at eps = V0, C + D x (q = 0), so there is no
    singular case. Independent of the closed-form T expression."""
    eps = float(eps)
    if eps <= 0:
        raise ValueError("eps must be > 0")
    k = np.sqrt(2.0 * eps)
    if eps == V0:
        # unknowns r, C, D, t:  1 + r = C;  ik(1 - r) = D;  C + D w = t;  D = ik t
        M = np.array([[-1, 1, 0, 0], [1j * k, 0, 1, 0], [0, 1, w, -1], [0, 0, 1, -1j * k]], dtype=complex)
        rhs = np.array([1, 1j * k, 0, 0], dtype=complex)
        r, C, D, t = np.linalg.solve(M, rhs)
        return complex(t), complex(r)
    q = np.sqrt(2.0 * (eps - V0) + 0j)
    e_p, e_m = np.exp(1j * q * w), np.exp(-1j * q * w)
    M = np.array([
        [-1, 1, 1, 0],
        [1j * k, 1j * q, -1j * q, 0],
        [0, e_p, e_m, -1],
        [0, 1j * q * e_p, -1j * q * e_m, -1j * k],
    ], dtype=complex)
    rhs = np.array([1, 1j * k, 0, 0], dtype=complex)
    r, A, B, t = np.linalg.solve(M, rhs)
    return complex(t), complex(r)
