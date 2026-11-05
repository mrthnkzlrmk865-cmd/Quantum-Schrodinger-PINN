"""
Fourier feature embedding (Tancik et al., NeurIPS 2020), used to
counteract the spectral bias of plain MLPs (Rahaman et al., ICML 2019):
networks fit low-frequency content first and are slow to fit high-
frequency content.

Each (normalized, in [0, 1]) input coordinate c is mapped to
    [c, sin(2 pi f_i c), cos(2 pi f_i c)]  for f_i = linspace(1, max_freq, F).

Frequencies are linearly (not octave) spaced and capped at a physically
motivated max_freq. Very high frequencies make the second derivatives of
the embedding scale like (2 pi f)^2, which blows up the PDE residual at
initialization (this happened in an earlier version; see README
"Debugging history"). For the PINN problem the packet's phase has kx0/2pi
~ 2.4 cycles per unit length and its envelope has spectral content up to a
few cycles, so max_freq = 6 covers the relevant range. This is a modelling
choice, and the PINN study compares it against no embedding.
"""

import torch
import torch.nn as nn


class FourierFeatureEmbedding(nn.Module):
    def __init__(self, in_dim: int, num_frequencies: int = 6, max_freq: float = 6.0):
        super().__init__()
        freqs = torch.linspace(1.0, max_freq, num_frequencies)
        self.register_buffer("freqs", freqs)
        self.in_dim = in_dim
        self.out_dim = in_dim * (1 + 2 * num_frequencies)

    def forward(self, coords: torch.Tensor) -> torch.Tensor:
        proj = coords.unsqueeze(-1) * self.freqs.view(1, 1, -1) * (2 * torch.pi)
        feats = torch.cat([coords.unsqueeze(-1), torch.sin(proj), torch.cos(proj)], dim=-1)
        return feats.reshape(coords.shape[0], -1)
