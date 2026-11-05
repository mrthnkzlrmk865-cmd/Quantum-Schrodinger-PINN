r"""
Gaussian-correlated random potential, defined in the continuum and
evaluated identically in NumPy (reference solver) and PyTorch (PINN).

Definition
----------
For a seed s, draw M independent wave vectors and phases

    q_j ~ N(0, I_2 / zeta^2),      phi_j ~ Uniform[0, 2 pi),   j = 1..M

and define the field

    V_dis(x, y) = W * sqrt(2/M) * sum_j cos(q_j . r + phi_j) * S(x)

where S(x) is an optional smooth window that restricts the disorder to a
slab x_a < x < x_b (S = 1 everywhere if no slab is set):

    S(x) = 1/2 [ tanh((x - x_a)/delta) - tanh((x - x_b)/delta) ].

This is the spectral-representation ("random phase") method for stationary
random fields (Shinozuka & Deodatis, Appl. Mech. Rev. 44, 191 (1991)).
Averaging over the ensemble of (q_j, phi_j) gives, for any M,

    <V_dis(r)> = 0
    <V_dis(r) V_dis(r')> = W^2 E_q[cos(q.(r - r'))] = W^2 exp(-|r - r'|^2 / (2 zeta^2))

inside the slab, i.e. the field has RMS amplitude W and a Gaussian
correlation function with correlation length zeta. Its power spectrum is
S(q) = 2 pi W^2 zeta^2 exp(-q^2 zeta^2 / 2). For large M the one-point
distribution tends to a Gaussian (central limit theorem). The sample
variance of one realization fluctuates around W^2; we do NOT re-normalize
each realization, so the ensemble statistics stay exactly as stated.

Why this representation
-----------------------
* It is a smooth analytic function of (x, y), so it can be evaluated at
  arbitrary PINN collocation points and is differentiable in PyTorch.
* One seed defines one physical field independent of any grid, so grid
  refinement studies compare the SAME potential, and the PINN and the
  reference solver see the SAME realization.
* The correlation length is an explicit, physically meaningful parameter.

The module also contains the lowest-order (Born / Fermi golden rule)
scattering rates for this disorder in an unbounded 2D system. They are
approximate benchmarks, valid only for weak disorder (k * l >> 1), and are
used only for comparison with the numerical results.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import torch
from scipy.special import ive


@dataclass(frozen=True)
class DisorderRealization:
    strength: float                # W (RMS amplitude, ensemble value)
    corr_length: float             # zeta
    qx: np.ndarray                 # (M,)
    qy: np.ndarray                 # (M,)
    phi: np.ndarray                # (M,)
    region: Optional[Tuple[float, float]] = None
    edge: float = 0.02
    seed: int = 0

    @property
    def n_modes(self) -> int:
        return self.qx.shape[0]

    # ------------------------------------------------------------------
    # Window
    # ------------------------------------------------------------------
    def window_numpy(self, x: np.ndarray) -> np.ndarray:
        if self.region is None:
            return np.ones_like(x)
        a, b = self.region
        d = self.edge
        return 0.5 * (np.tanh((x - a) / d) - np.tanh((x - b) / d))

    def window_torch(self, x: torch.Tensor) -> torch.Tensor:
        if self.region is None:
            return torch.ones_like(x)
        a, b = self.region
        d = self.edge
        return 0.5 * (torch.tanh((x - a) / d) - torch.tanh((x - b) / d))

    # ------------------------------------------------------------------
    # Field evaluation
    # ------------------------------------------------------------------
    def numpy(self, X: np.ndarray, Y: np.ndarray, chunk: int = 64) -> np.ndarray:
        """Evaluate on arrays of any (equal) shape, float64."""
        X = np.asarray(X, dtype=np.float64)
        Y = np.asarray(Y, dtype=np.float64)
        if self.strength == 0.0:
            return np.zeros_like(X)
        out = np.zeros_like(X)
        for s in range(0, self.n_modes, chunk):
            qx = self.qx[s:s + chunk]
            qy = self.qy[s:s + chunk]
            ph = self.phi[s:s + chunk]
            out += np.cos(X[..., None] * qx + Y[..., None] * qy + ph).sum(axis=-1)
        out *= self.strength * math.sqrt(2.0 / self.n_modes)
        return out * self.window_numpy(X)

    def torch(self, x: torch.Tensor, y: torch.Tensor, chunk: int = 128) -> torch.Tensor:
        """Differentiable evaluation; x, y of shape (N, 1) (any float dtype).
        Computed in the dtype of x; parameters are cast from float64."""
        if self.strength == 0.0:
            return torch.zeros_like(x)
        dt, dev = x.dtype, x.device
        qx = torch.as_tensor(self.qx, dtype=dt, device=dev)
        qy = torch.as_tensor(self.qy, dtype=dt, device=dev)
        ph = torch.as_tensor(self.phi, dtype=dt, device=dev)
        out = torch.zeros_like(x)
        for s in range(0, self.n_modes, chunk):
            arg = x * qx[s:s + chunk] + y * qy[s:s + chunk] + ph[s:s + chunk]
            out = out + torch.cos(arg).sum(dim=-1, keepdim=True)
        out = out * (self.strength * math.sqrt(2.0 / self.n_modes))
        return out * self.window_torch(x)


def make_disorder(strength: float, corr_length: float, seed: int,
                  n_modes: int = 512, region=None, edge: float = 0.02) -> DisorderRealization:
    """Draw one realization. The random numbers depend only on
    (seed, n_modes); the field scales with corr_length and strength, so
    the same seed at different zeta or W gives a 'rescaled' copy of the same
    random draw (useful for paired comparisons)."""
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), 20241105]))
    g = rng.standard_normal((n_modes, 2))
    phi = rng.uniform(0.0, 2.0 * np.pi, size=n_modes)
    q = g / corr_length
    return DisorderRealization(
        strength=float(strength), corr_length=float(corr_length),
        qx=q[:, 0].copy(), qy=q[:, 1].copy(), phi=phi,
        region=None if region is None else (float(region[0]), float(region[1])),
        edge=float(edge), seed=int(seed),
    )


def disorder_from_config(cfg, seed: Optional[int] = None) -> DisorderRealization:
    return make_disorder(
        cfg.disorder_strength, cfg.disorder_corr_length,
        cfg.disorder_seed if seed is None else seed,
        n_modes=cfg.disorder_n_modes, region=cfg.disorder_region,
        edge=cfg.disorder_edge,
    )


# ----------------------------------------------------------------------
# Theory: correlation function and Born-approximation scattering rates
# ----------------------------------------------------------------------
def correlation_function(r, strength: float, corr_length: float):
    """Target two-point function C(r) = W^2 exp(-r^2 / (2 zeta^2))."""
    return strength**2 * np.exp(-np.asarray(r) ** 2 / (2.0 * corr_length**2))


def born_rates(strength, corr_length, k):
    r"""Fermi-golden-rule rates for a particle with wave number k (energy
    k^2/2, hbar = m = 1) in an unbounded 2D Gaussian-correlated potential:

        1/tau_s  = (1/2pi) \int_0^{2pi} S(2k sin(theta/2)) d theta
                 = 2 pi W^2 zeta^2 e^{-a} I_0(a)
        1/tau_tr = (1/2pi) \int S(.) (1 - cos theta) d theta
                 = 2 pi W^2 zeta^2 e^{-a} [I_0(a) - I_1(a)],     a = k^2 zeta^2

    Returns dict with rates, scattering mean free path l_s = k tau_s and
    transport mean free path l_tr = k tau_tr (velocity = k).
    Valid only for k * l >> 1 (weak disorder). The rates are not monotonic
    in zeta at fixed W: they grow like zeta^2 for k zeta << 1 (white-noise
    limit) and fall like 1/zeta for k zeta >> 1 (smooth disorder scatters
    mainly forward).
    """
    a = (np.asarray(k) * np.asarray(corr_length)) ** 2
    pref = 2.0 * np.pi * np.asarray(strength) ** 2 * np.asarray(corr_length) ** 2
    gs = pref * ive(0, a)
    gtr = pref * (ive(0, a) - ive(1, a))
    with np.errstate(divide="ignore"):
        ls = np.asarray(k) / gs
        ltr = np.asarray(k) / gtr
    return dict(rate_s=gs, rate_tr=gtr, l_s=ls, l_tr=ltr, a=a)
