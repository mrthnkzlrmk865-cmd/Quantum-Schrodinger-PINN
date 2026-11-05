"""
PINN mapping (x, y, t) -> (u, v), the real and imaginary parts of
psi = u + i v, with the initial and boundary conditions built into the
output (hard constraints):

    u(x,y,t) = u0(x,y) + A * B(x,y) * tau * N_u(x/Lx, y/Ly, tau)
    v(x,y,t) = v0(x,y) + A * B(x,y) * tau * N_v(x/Lx, y/Ly, tau)

  tau = t/T in [0, 1]                    normalized time
  B   = 16 (x/Lx)(1-x/Lx)(y/Ly)(1-y/Ly)  zero on all walls, max 1
  A   = max |psi0|                       output scale, so N is O(1)
  u0, v0                                 the exact initial condition
                                         (initial_condition.py)

At t = 0, u = u0 and v = v0 exactly. On the walls u = u0(wall),
v = v0(wall), which is < 4e-5 for the PINN preset (the packet is 5 sigma
from the nearest wall). So the Dirichlet condition holds up to that
tolerance; tests/test_pinn.py checks it.

Inputs are normalized to [0, 1] before the (optional) Fourier embedding,
so the same network hyperparameters can be used for any domain size or
time window.

Why hard constraints: a soft-penalty version (separate IC/BC losses)
reproducibly collapsed to the trivial psi ~ 0 solution in this project's
history (README "Debugging history").
"""

import math

import torch
import torch.nn as nn

from embeddings import FourierFeatureEmbedding
from initial_condition import initial_condition_torch


class SchrodingerPINN(nn.Module):
    def __init__(self, cfg):
        super().__init__()
        self.cfg = cfg
        if cfg.use_fourier:
            self.embedding = FourierFeatureEmbedding(3, cfg.fourier_num_freq, cfg.fourier_max_freq)
            in_features = self.embedding.out_dim
        else:
            self.embedding = None
            in_features = 3
        layers = [nn.Linear(in_features, cfg.hidden_units), nn.Tanh()]
        for _ in range(cfg.hidden_layers - 1):
            layers += [nn.Linear(cfg.hidden_units, cfg.hidden_units), nn.Tanh()]
        layers += [nn.Linear(cfg.hidden_units, 2)]
        self.net = nn.Sequential(*layers)
        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)
        if cfg.packet == "gaussian":
            self.amp = 1.0 / math.sqrt(math.pi * cfg.sigma_x * cfg.sigma_y)
        else:
            self.amp = (math.pi * cfg.sigma_x**2) ** (-0.25) * math.sqrt(2.0 / cfg.Ly)

    def raw(self, xn, yn, tau):
        inp = torch.cat([xn, yn, tau], dim=1)
        if self.embedding is not None:
            inp = self.embedding(inp)
        out = self.net(inp)
        return out[:, 0:1], out[:, 1:2]

    def forward(self, x, y, t):
        c = self.cfg
        xn, yn, tau = x / c.Lx, y / c.Ly, t / c.T
        nu, nv = self.raw(xn, yn, tau)
        u0, v0 = initial_condition_torch(c, x, y)
        bump = 16.0 * xn * (1 - xn) * yn * (1 - yn)
        g = self.amp * bump * tau
        return u0 + g * nu, v0 + g * nv
