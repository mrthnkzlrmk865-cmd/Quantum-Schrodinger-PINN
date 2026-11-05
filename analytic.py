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


def rectangular_barrier_transmission(eps, V0, w):
    """Exact transmission probability of a sharp rectangular barrier of
    height V0 > 0 and width w for a particle of energy eps > 0 (hbar = m = 1):

        eps > V0:  T = [1 + V0^2 sin^2(q w)  / (4 eps (eps - V0))]^-1,  q = sqrt(2 (eps - V0))
        eps < V0:  T = [1 + V0^2 sinh^2(k w) / (4 eps (V0 - eps))]^-1,  k = sqrt(2 (V0 - eps))

    Standard textbook result. In this project it is checked against an
    independent exact plane-wave matching calculation
    (tests/test_smatrix_spec.py::test_rectangular_barrier_formula_vs_matching).
    Undefined at eps == V0 (use the limit or avoid that energy)."""
    eps = np.asarray(eps, dtype=float)
    out = np.empty_like(eps)
    above = eps > V0
    q = np.sqrt(2.0 * (eps[above] - V0))
    out[above] = 1.0 / (1.0 + V0**2 * np.sin(q * w) ** 2 / (4.0 * eps[above] * (eps[above] - V0)))
    k = np.sqrt(2.0 * (V0 - eps[~above]))
    out[~above] = 1.0 / (1.0 + V0**2 * np.sinh(k * w) ** 2 / (4.0 * eps[~above] * (V0 - eps[~above])))
    return out
