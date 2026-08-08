"""
Fourier feature embedding, used to counteract the "spectral bias" of plain
MLPs (Rahaman et al., 2019: neural networks preferentially fit low-frequency
functions first and are slow/unable to fit high-frequency ones at typical
training budgets).

Root-cause context for this project
------------------------------------
The initial wave packet carries momentum kx0=40, giving it a spatial phase
cos(kx0*x), sin(kx0*x) that oscillates roughly every 2*pi/40 ~ 0.157 units --
about 2-3 oscillations within one packet width (sigma=0.06). A raw
(x, y, t) -> Tanh-MLP struggles to represent this (verified empirically: a
pure supervised fit of the network to the initial condition alone plateaus
around MSE ~0.01 and does not improve with more steps -- see README
"Debugging note"). Mapping inputs through sinusoids at multiple frequencies
before the MLP (Tancik et al., 2020, "Fourier Features Let Networks Learn
High Frequency Functions in Low Dimensional Domains") gives the network
direct access to the relevant frequency content and removes this bottleneck.
"""

import torch
import torch.nn as nn


class FourierFeatureEmbedding(nn.Module):
    """
    Maps each input coordinate to [coordinate, sin(2*pi*f_i*coordinate),
    cos(2*pi*f_i*coordinate) for f_i in frequencies], concatenated across
    coordinates.

    Frequency range: linearly spaced from 1 to max_freq, NOT octave-spaced
    (1,2,4,8,...). Octave spacing was tried first and reached very high
    frequencies (up to 128) that are unnecessary for this problem (the
    packet's phase oscillates at kx0/(2*pi) ~= 6.4 cycles, and its envelope,
    sigma=0.06, has meaningful spectral content up to roughly 1/sigma ~= 17
    cycles) -- but those high frequencies caused d^2/dx^2 of the embedding
    to blow up by a factor of (2*pi*f)^2 (~6.4e6 at f=128), which exploded
    the PDE-residual loss once it passed through the network's second
    derivatives via autograd (see README "Debugging note"). Linear spacing
    up to a physically-motivated max_freq avoids that.
    """

    def __init__(self, in_dim: int, num_frequencies: int = 6, max_freq: float = 10.0):
        super().__init__()
        freqs = torch.linspace(1.0, max_freq, num_frequencies)
        self.register_buffer("freqs", freqs)
        self.in_dim = in_dim
        self.out_dim = in_dim * (1 + 2 * num_frequencies)

    def forward(self, coords: torch.Tensor) -> torch.Tensor:
        # coords: (N, in_dim)
        proj = coords.unsqueeze(-1) * self.freqs.view(1, 1, -1) * 2 * torch.pi  # (N, in_dim, F)
        sin_feat = torch.sin(proj)
        cos_feat = torch.cos(proj)
        feats = torch.cat(
            [coords.unsqueeze(-1), sin_feat, cos_feat], dim=-1
        )  # (N, in_dim, 1+2F)
        return feats.reshape(coords.shape[0], -1)
