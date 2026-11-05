"""
Physical observables computed from a (batch of) wavefunction(s) on the
solver grid. All functions accept psi with shape (..., nx, ny) and return
arrays with the leading shape (...). The grid spacings are dx, dy.

Definitions (rho = |psi|^2, integrals approximated by grid sums):

  norm           N    = int rho
  region prob.   P_A  = int_A rho
  centre of mass <r>  = int r rho / N
  MSD            <|r - r0|^2> = int |r - r0|^2 rho / N   (r0 = initial centre)
  participation  PR   = (int rho)^2 / int rho^2     [units of area]
  ratio          IPR  = 1 / PR

PR is the effective area occupied by the probability distribution: for a
uniform density over an area A it equals A, and for a Gaussian
rho ~ exp(-r^2/s^2) it equals 2 pi s^2. It is meaningful for the closed-box
spreading experiment. In the open waveguide it is not used, because
probability leaves the domain.
"""

from __future__ import annotations

import numpy as np


def density(psi):
    return np.abs(psi) ** 2


def norm(psi, dx, dy):
    return density(psi).sum(axis=(-2, -1)) * dx * dy


def region_probability(psi, mask, dx, dy):
    return np.where(mask, density(psi), 0.0).sum(axis=(-2, -1)) * dx * dy


def centre_of_mass(psi, X, Y, dx, dy):
    rho = density(psi)
    n = rho.sum(axis=(-2, -1))
    return (rho * X).sum(axis=(-2, -1)) / n, (rho * Y).sum(axis=(-2, -1)) / n


def mean_square_displacement(psi, X, Y, r0):
    rho = density(psi)
    d2 = (X - r0[0]) ** 2 + (Y - r0[1]) ** 2
    return (rho * d2).sum(axis=(-2, -1)) / rho.sum(axis=(-2, -1))


def participation_ratio(psi, dx, dy):
    rho = density(psi)
    num = (rho.sum(axis=(-2, -1)) * dx * dy) ** 2
    den = (rho**2).sum(axis=(-2, -1)) * dx * dy
    return num / den


def inverse_participation_ratio(psi, dx, dy):
    return 1.0 / participation_ratio(psi, dx, dy)


def uniform_box_msd(Lx, Ly, r0):
    """MSD of a uniform distribution over [0,Lx]x[0,Ly] about r0
    (the 'fully spread / ergodic' reference value)."""
    x0, y0 = r0
    mx = (Lx**2) / 3 - Lx * x0 + x0**2
    my = (Ly**2) / 3 - Ly * y0 + y0**2
    return mx + my


def ergodic_box_pr(Lx, Ly):
    """PR of the long-time-average density of an ergodic (random-wave)
    state in a box with Dirichlet walls. For random superpositions the
    density fluctuates (Porter-Thomas / random-wave statistics give
    <rho^2> = 3 <rho>^2 for real waves and 2 <rho>^2 for complex waves),
    so even a fully delocalized state has PR ~ A/2 to A/3, not A.
    We return A/2 (complex random wave). This sets the scale only."""
    return Lx * Ly / 2.0
