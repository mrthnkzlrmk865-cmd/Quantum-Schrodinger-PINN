"""
Initial wave packet psi(x, y, 0), one definition evaluated by both the
PINN (torch, hard IC constraint) and the reference solver (numpy).

"gaussian" (closed boxes, PINN problem):
    psi0 = A exp(-(x-x0)^2/(2 sx^2) - (y-y0)^2/(2 sy^2)) exp(i(kx0 x + ky0 y)),
    A = 1/sqrt(pi sx sy)  so that  int |psi0|^2 dx dy = 1 on the whole plane.

"waveguide_mode" (open waveguide):
    psi0 = (pi sx^2)^(-1/4) exp(-(x-x0)^2/(2 sx^2)) exp(i kx0 x)
           * sqrt(2/Ly) sin(n pi y / Ly),
    i.e. a Gaussian packet in the n-th transverse Dirichlet mode; it
    satisfies psi = 0 on y = 0, Ly exactly and is normalized on the strip.

On a finite box the Gaussian tails at the walls are not exactly zero. For
all presets the analytic norm on the box differs from 1 by < 1e-8; tests
check this. Momentum-space width: |psi(k)|^2 ~ exp(-k^2 s^2), i.e. the
standard deviation of k_x is 1/(sqrt(2) sx).

(Historical note: an earlier version used an unnormalized envelope with
peak 1, which contradicted the norm = 1 loss target. See README.)
"""

from __future__ import annotations

import math

import numpy as np
import torch


def _amplitude_and_phase(cfg, x, y, lib):
    if cfg.packet == "gaussian":
        A = 1.0 / math.sqrt(math.pi * cfg.sigma_x * cfg.sigma_y)
        env = A * lib.exp(-((x - cfg.x0) ** 2) / (2 * cfg.sigma_x**2)
                          - ((y - cfg.y0) ** 2) / (2 * cfg.sigma_y**2))
        phase = cfg.kx0 * x + cfg.ky0 * y
    elif cfg.packet == "waveguide_mode":
        A = (math.pi * cfg.sigma_x**2) ** (-0.25) * math.sqrt(2.0 / cfg.Ly)
        env = A * lib.exp(-((x - cfg.x0) ** 2) / (2 * cfg.sigma_x**2)) \
            * lib.sin(cfg.mode_n * math.pi * y / cfg.Ly)
        phase = cfg.kx0 * x
    else:
        raise ValueError(f"unknown packet type {cfg.packet!r}")
    return env, phase


def initial_condition_numpy(cfg, X, Y):
    """Complex float64 array."""
    env, phase = _amplitude_and_phase(cfg, np.asarray(X, float), np.asarray(Y, float), np)
    return env * np.exp(1j * phase)


def initial_condition_torch(cfg, x, y):
    """Returns (u0, v0) = (Re psi0, Im psi0)."""
    env, phase = _amplitude_and_phase(cfg, x, y, torch)
    return env * torch.cos(phase), env * torch.sin(phase)
