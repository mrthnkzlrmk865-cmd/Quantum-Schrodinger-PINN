import numpy as np
import torch

from config import pinn_config, waveguide_config, closed_box_config
from initial_condition import initial_condition_numpy, initial_condition_torch
from reference_solver import build_grid


def _grid_norm(cfg, n=400):
    x = (np.arange(n) + 0.5) * cfg.Lx / n
    y = (np.arange(n) + 0.5) * cfg.Ly / n
    X, Y = np.meshgrid(x, y, indexing="ij")
    return np.sum(np.abs(initial_condition_numpy(cfg, X, Y)) ** 2) * cfg.Lx * cfg.Ly / n**2


def test_normalized_on_box():
    for cfg in (pinn_config(), waveguide_config(), closed_box_config(1.0)):
        assert abs(_grid_norm(cfg) - 1.0) < 1e-6


def test_torch_numpy_ic_agree():
    cfg = pinn_config()
    X, Y, _, _ = build_grid(cfg, 31, 31)
    p = initial_condition_numpy(cfg, X, Y)
    u, v = initial_condition_torch(cfg, torch.tensor(X), torch.tensor(Y))
    assert np.allclose(u.numpy(), p.real, atol=1e-12) and np.allclose(v.numpy(), p.imag, atol=1e-12)


def test_ic_small_on_walls():
    cfg = pinn_config()
    s = np.linspace(0, 1, 201)
    z = np.zeros_like(s)
    o = np.ones_like(s)
    wall = np.concatenate([initial_condition_numpy(cfg, a, b) for a, b in ((z, s), (o, s), (s, z), (s, o))])
    assert np.abs(wall).max() < 5e-5
    cfg = waveguide_config()
    assert np.abs(initial_condition_numpy(cfg, s * cfg.Lx, z)).max() < 1e-12  # exact transverse zero
