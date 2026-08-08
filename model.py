"""
Fully-connected network mapping (x, y, t) -> (u, v), the real and imaginary
parts of the wavefunction Psi(x, y, t) = u + i*v.
"""

import torch
import torch.nn as nn

from config import cfg
from embeddings import FourierFeatureEmbedding
from initial_condition import initial_condition_torch


class SchrodingerPINN(nn.Module):
    """
    Outputs (u, v) with the initial condition and Dirichlet boundary
    condition satisfied BY CONSTRUCTION, instead of via soft penalty losses.

    Parametrization
    ---------------
        u(x,y,t) = u0(x,y) + B(x,y) * t * NN_u(x,y,t)
        v(x,y,t) = v0(x,y) + B(x,y) * t * NN_v(x,y,t)

    where u0, v0 is the exact analytic initial condition, and
    B(x,y) = 16*x*(L-x)*y*(L-y)/L^4 is a smooth bump that is exactly zero on
    all four walls of the box (normalized so its peak value is 1).

    At t=0, the t-multiplied term vanishes, so u=u0, v=v0 EXACTLY -- no IC
    loss term needed. At the walls, B=0 for all t, so u=u0(wall), v=v0(wall)
    -- not exactly 0, but the Gaussian packet (sigma=0.06) is centered well
    inside the box, so u0/v0 at any wall are numerically negligible
    (~1e-6 or smaller for this configuration).

    Why: a soft-penalty version of this network (plain MLP + separate IC/BC
    MSE loss terms) reproducibly collapsed to the trivial near-zero solution
    during training -- the network could reach a deceptively low combined
    loss without ever learning the packet's actual shape (see README
    "Debugging note" for the full diagnosis, which also covers a spectral
    bias issue fixed separately via FourierFeatureEmbedding below). Hard-
    constraining the IC removes the optimization pathway that let this
    happen: the initial amplitude can no longer be optimized away to zero.
    """

    def __init__(
        self,
        in_dim=3,
        out_dim=2,
        hidden=cfg.hidden_units,
        layers=cfg.hidden_layers,
        num_frequencies=6,
    ):
        super().__init__()
        self.embedding = FourierFeatureEmbedding(in_dim, num_frequencies=num_frequencies)
        net = [nn.Linear(self.embedding.out_dim, hidden), nn.Tanh()]
        for _ in range(layers - 1):
            net += [nn.Linear(hidden, hidden), nn.Tanh()]
        net += [nn.Linear(hidden, out_dim)]
        self.net = nn.Sequential(*net)

        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def _raw_forward(self, x: torch.Tensor, y: torch.Tensor, t: torch.Tensor):
        inp = torch.cat([x, y, t], dim=1)
        emb = self.embedding(inp)
        out = self.net(emb)
        return out[:, 0:1], out[:, 1:2]

    def forward(self, x: torch.Tensor, y: torch.Tensor, t: torch.Tensor):
        nn_u, nn_v = self._raw_forward(x, y, t)
        u0, v0 = initial_condition_torch(x, y)
        bump = 16.0 * x * (cfg.L - x) * y * (cfg.L - y) / (cfg.L**4)
        u = u0 + bump * t * nn_u
        v = v0 + bump * t * nn_v
        return u, v
