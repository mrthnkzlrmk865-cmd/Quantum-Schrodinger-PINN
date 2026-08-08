"""Initial wave packet, shared by model.py (hard IC constraint) and
pinn_core.py (used to build the reference IC values, and previously the
soft IC loss term)."""

import math

import torch

from config import cfg

# Normalization constant so that integral(u0^2 + v0^2) dx dy = 1 over the
# whole plane, i.e. A such that A^2 * pi * sigma^2 = 1.
#
# BUG THIS FIXES: the original (and an earlier draft of this rewrite) used
# an unnormalized envelope with peak amplitude 1.0, whose total probability
# is only pi*sigma^2 (~0.011 for sigma=0.06) -- not 1. That silently made
# the "IC" and the "normalization = 1" training targets mutually
# contradictory: the network was being pushed to both reproduce the
# (small-probability) unnormalized packet exactly AND have total
# probability 1 at every time including t=0. This was the root cause of a
# training collapse that survived several other fixes (importance
# sampling, Fourier features, hard IC/BC constraints) -- see README
# "Debugging note" for the full story of how this was tracked down.
_NORM_CONST = 1.0 / (cfg.sigma * math.sqrt(math.pi))


def initial_condition_torch(x: torch.Tensor, y: torch.Tensor):
    envelope = _NORM_CONST * torch.exp(
        -((x - cfg.x0) ** 2 + (y - cfg.y0) ** 2) / (2 * cfg.sigma**2)
    )
    phase = cfg.kx0 * x + cfg.ky0 * y
    u0 = envelope * torch.cos(phase)
    v0 = envelope * torch.sin(phase)
    return u0, v0
