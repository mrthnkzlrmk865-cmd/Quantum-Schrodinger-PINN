"""
Real potential V(x, y) = V_barrier(x) + V_dis(x, y), defined ONCE and
evaluated by both solvers:

  * ``Potential.torch(x, y)``  differentiable, used in the PINN residual.
  * ``Potential.numpy(X, Y)``  float64, used by the split-step solver.

Both evaluate the same formulas with the same parameters (including the
same disorder realization), and tests/test_potential.py checks pointwise
agreement. The complex absorbing potential is NOT part of V; it belongs
to the open-boundary reference solver only (reference_solver.py).

Barrier
-------
    V_b(x) = V0 * 1/2 [tanh((x - x1)/e) - tanh((x - x2)/e)],
    x1,2 = barrier_center -/+ barrier_width/2, e = barrier_edge.
With barrier_edge = 0 the barrier is the sharp strip x1 < x < x2 (the
original definition). The smooth version is preferred for convergence
studies because the potential then does not depend on where grid points
fall relative to the edge.
"""

from __future__ import annotations

import numpy as np
import torch

from disorder import DisorderRealization, disorder_from_config


class Potential:
    def __init__(self, cfg, disorder: DisorderRealization | None = None):
        self.V0 = float(cfg.barrier_height)
        self.xc = float(cfg.barrier_center)
        self.w = float(cfg.barrier_width)
        self.edge = float(cfg.barrier_edge)
        if disorder is None and cfg.disorder_strength != 0.0:
            disorder = disorder_from_config(cfg)
        self.disorder = disorder

    # ------------------------------------------------------------------
    def barrier_numpy(self, X):
        X = np.asarray(X, dtype=np.float64)
        if self.V0 == 0.0:
            return np.zeros_like(X)
        x1, x2 = self.xc - self.w / 2, self.xc + self.w / 2
        if self.edge > 0:
            return self.V0 * 0.5 * (np.tanh((X - x1) / self.edge) - np.tanh((X - x2) / self.edge))
        return np.where((X > x1) & (X < x2), self.V0, 0.0)

    def barrier_torch(self, x):
        if self.V0 == 0.0:
            return torch.zeros_like(x)
        x1, x2 = self.xc - self.w / 2, self.xc + self.w / 2
        if self.edge > 0:
            return self.V0 * 0.5 * (torch.tanh((x - x1) / self.edge) - torch.tanh((x - x2) / self.edge))
        return torch.where((x > x1) & (x < x2), torch.full_like(x, self.V0), torch.zeros_like(x))

    # ------------------------------------------------------------------
    def numpy(self, X, Y):
        V = self.barrier_numpy(X)
        if self.disorder is not None:
            V = V + self.disorder.numpy(X, Y)
        return V

    def torch(self, x, y):
        V = self.barrier_torch(x)
        if self.disorder is not None:
            V = V + self.disorder.torch(x, y)
        return V

    __call__ = torch


# Backwards-compatible functional API (default config) -----------------
def potential_numpy(X, Y, cfg=None):
    from config import cfg as default_cfg
    return Potential(cfg or default_cfg).numpy(X, Y)


def potential_torch(x, y, cfg=None):
    from config import cfg as default_cfg
    return Potential(cfg or default_cfg).torch(x, y)
