"""Potential: torch/numpy agreement, disorder reproducibility and statistics."""
import numpy as np
import pytest
import torch

from config import pinn_config, waveguide_config
from disorder import born_rates, make_disorder
from potential import Potential


@pytest.mark.parametrize("edge", [0.0, 0.01])
def test_torch_numpy_pointwise_agreement(edge):
    cfg = pinn_config(barrier_height=150.0, barrier_center=0.6, barrier_edge=edge,
                      disorder_strength=56.0, disorder_corr_length=0.08, disorder_seed=3)
    pot = Potential(cfg)
    rng = np.random.default_rng(0)
    x = rng.uniform(0, 1, 5000)
    y = rng.uniform(0, 1, 5000)
    vn = pot.numpy(x, y)
    vt = pot.torch(torch.tensor(x[:, None]), torch.tensor(y[:, None])).numpy()[:, 0]
    assert np.max(np.abs(vn - vt)) < 1e-9 * max(1.0, np.abs(vn).max())
    # float32 path used in PINN training
    vt32 = pot.torch(torch.tensor(x[:, None], dtype=torch.float32),
                     torch.tensor(y[:, None], dtype=torch.float32)).numpy()[:, 0]
    if edge > 0:  # sharp step can flip at the edge in float32
        assert np.max(np.abs(vn - vt32)) < 1e-3 * np.abs(vn).max()


def test_same_seed_same_field_any_resolution():
    d1 = make_disorder(50.0, 0.05, seed=7)
    d2 = make_disorder(50.0, 0.05, seed=7)
    d3 = make_disorder(50.0, 0.05, seed=8)
    x = np.linspace(0, 1, 33)
    X, Y = np.meshgrid(x, x, indexing="ij")
    assert np.array_equal(d1.numpy(X, Y), d2.numpy(X, Y))
    assert not np.allclose(d1.numpy(X, Y), d3.numpy(X, Y))
    # a coarse grid samples exactly the same continuum field as a fine grid
    xf = np.linspace(0, 1, 65)
    Xf, Yf = np.meshgrid(xf, xf, indexing="ij")
    assert np.allclose(d1.numpy(Xf, Yf)[::2, ::2], d1.numpy(X, Y), atol=1e-12)


def test_disorder_is_differentiable():
    d = make_disorder(10.0, 0.1, seed=1)
    x = torch.tensor([[0.3]], dtype=torch.float64, requires_grad=True)
    y = torch.tensor([[0.4]], dtype=torch.float64)
    V = d.torch(x, y)
    (g,) = torch.autograd.grad(V.sum(), x)
    h = 1e-6
    fd = (d.numpy(np.array(0.3 + h), np.array(0.4)) - d.numpy(np.array(0.3 - h), np.array(0.4))) / (2 * h)
    assert abs(g.item() - float(fd)) < 1e-5 * max(1, abs(float(fd)))


def test_ensemble_mean_variance_and_correlation():
    """Ensemble over seeds reproduces <V>=0, <V^2>=W^2 and
    C(r)=W^2 exp(-r^2/(2 zeta^2)) within sampling error."""
    W, zeta = 1.0, 0.1
    r = np.array([0.0, 0.05, 0.1, 0.2])
    pts = np.array([0.5 + r, np.full_like(r, 0.5)])
    vals = np.array([make_disorder(W, zeta, seed=s, n_modes=256).numpy(pts[0], pts[1])
                     for s in range(3000)])
    v0 = vals[:, 0]
    assert abs(v0.mean()) < 4 * W / np.sqrt(3000)
    assert abs(v0.var() - W**2) < 0.1
    C = (vals * v0[:, None]).mean(axis=0)
    target = W**2 * np.exp(-r**2 / (2 * zeta**2))
    assert np.max(np.abs(C - target)) < 0.08


def test_slab_window_confines_disorder():
    cfg = waveguide_config(disorder_strength=100.0)
    pot = Potential(cfg)
    a, b = cfg.disorder_region
    y = np.linspace(0.01, 0.99, 50)
    for x in (a - 0.3, b + 0.3):
        assert np.max(np.abs(pot.numpy(np.full_like(y, x), y))) < 1e-6 * 100


def test_born_rates_limits():
    W, k = 10.0, 20.0
    # white-noise limit: rate_tr ~ rate_s ~ 2 pi W^2 zeta^2
    r = born_rates(W, 1e-4, k)
    assert np.isclose(r["rate_s"], 2 * np.pi * W**2 * 1e-8, rtol=1e-3)
    # non-monotonic in zeta at fixed W
    z = np.linspace(0.005, 0.5, 400)
    rt = born_rates(W, z, k)["rate_tr"]
    i = np.argmax(rt)
    assert 0 < i < len(z) - 1
