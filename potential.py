"""
Potential energy V(x, y) for the 2D infinite well / tunneling scenario.

Implemented twice with an identical mathematical definition:
  - `potential_torch`: differentiable, used inside the PINN's autograd loss.
  - `potential_numpy`: used by the classical split-step reference solver.

Keeping two implementations (instead of forcing one to serve both) avoids
autograd-vs-numpy shape headaches while guaranteeing they encode the same V.
A unit test (tests/test_potential.py) checks that they agree pointwise.
"""

import numpy as np
import torch

from config import cfg


def _disorder_field_numpy(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """
    Smooth random potential roughness, used only in SCENARIO == "disordered".

    Built by band-limiting white noise in Fourier space so the disorder has
    a controllable spatial correlation length (`cfg.disorder_correlation_length`)
    instead of being pixel-to-pixel uncorrelated. This mimics the kind of
    static energetic disorder used in Environment-Assisted Quantum Transport
    (ENAQT) models of chromophore networks, where site energies fluctuate
    smoothly rather than randomly at every point.
    """
    rng = np.random.default_rng(cfg.disorder_seed)
    n = x.shape[0]
    noise = rng.normal(size=(n, n))
    # Gaussian smoothing via FFT convolution
    kx = np.fft.fftfreq(n, d=1.0 / n)
    ky = np.fft.fftfreq(n, d=1.0 / n)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    sigma_k = 1.0 / (2 * np.pi * cfg.disorder_correlation_length * n)
    kernel = np.exp(-(KX**2 + KY**2) * (sigma_k**2) / 2.0)
    field = np.fft.ifft2(np.fft.fft2(noise) * kernel).real
    field -= field.mean()
    std = field.std()
    if std > 1e-12:
        field /= std
    return field * cfg.disorder_strength


def potential_torch(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    V = torch.zeros_like(x)
    if cfg.SCENARIO in ("tunneling", "disordered"):
        mask = (x > cfg.barrier_center - cfg.barrier_width / 2) & (
            x < cfg.barrier_center + cfg.barrier_width / 2
        )
        V = torch.where(mask, torch.full_like(x, cfg.barrier_height), V)
    if cfg.SCENARIO == "disordered":
        # NOTE: the disorder field is not differentiable w.r.t. (x, y) in a
        # physically meaningful way when built from a fixed random grid, so
        # the disordered scenario is intended for the *reference solver*
        # comparison (see experiment_disorder.py), not for PINN training
        # directly. Raise loudly rather than silently returning a wrong V.
        raise NotImplementedError(
            "Disordered potential is only implemented for the numpy "
            "reference solver (potential_numpy). See experiment_disorder.py."
        )
    return V


def potential_numpy(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    V = np.zeros_like(X)
    if cfg.SCENARIO in ("tunneling", "disordered"):
        mask = (X > cfg.barrier_center - cfg.barrier_width / 2) & (
            X < cfg.barrier_center + cfg.barrier_width / 2
        )
        V = np.where(mask, cfg.barrier_height, V)
    if cfg.SCENARIO == "disordered":
        V = V + _disorder_field_numpy(X, Y)
    return V
