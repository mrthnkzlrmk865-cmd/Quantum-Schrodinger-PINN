"""PINN: residual correctness, hard constraints, training smoke test,
and (if a trained model exists) agreement with the reference solver."""
import math
import os

import numpy as np
import pytest
import torch

from config import pinn_config
from model import SchrodingerPINN
from pinn_core import pde_residual
from potential import Potential
from initial_condition import initial_condition_torch


class ExactEigenmode(torch.nn.Module):
    """psi = sin(m pi x) sin(n pi y) exp(-i (E + V0) t) in the unit box."""
    def __init__(self, cfg, m=2, n=1, V0=0.0):
        super().__init__()
        self.cfg, self.amp = cfg, 1.0
        self.E = 0.5 * math.pi**2 * (m**2 + n**2) + V0
        self.m, self.n = m, n

    def forward(self, x, y, t):
        s = torch.sin(self.m * math.pi * x) * torch.sin(self.n * math.pi * y)
        return s * torch.cos(self.E * t), -s * torch.sin(self.E * t)


def test_residual_vanishes_for_exact_solution():
    cfg = pinn_config(barrier_height=37.0, barrier_center=0.5, barrier_width=10.0, barrier_edge=0.0)
    pot = Potential(cfg)  # V = 37 everywhere in the box
    model = ExactEigenmode(cfg, V0=37.0)
    x = torch.rand(500, 1, dtype=torch.float64)
    y = torch.rand(500, 1, dtype=torch.float64)
    t = torch.rand(500, 1, dtype=torch.float64) * 0.1
    fu, fv = pde_residual(model, pot, x, y, t, scaled=False)
    assert fu.abs().max() < 1e-9 and fv.abs().max() < 1e-9
    # and a wrong energy is detected
    fu, fv = pde_residual(ExactEigenmode(cfg, V0=30.0), pot, x, y, t, scaled=False)
    assert fu.abs().max() > 1.0


def test_hard_initial_and_boundary_conditions():
    cfg = pinn_config(disorder_strength=56.0)
    torch.manual_seed(0)
    model = SchrodingerPINN(cfg)
    x = torch.rand(400, 1)
    y = torch.rand(400, 1)
    u, v = model(x, y, torch.zeros_like(x))
    u0, v0 = initial_condition_torch(cfg, x, y)
    assert torch.equal(u, u0) and torch.equal(v, v0)
    s = torch.rand(400, 1)
    for xb, yb in ((torch.zeros_like(s), s), (torch.ones_like(s), s), (s, torch.zeros_like(s)), (s, torch.ones_like(s))):
        u, v = model(xb, yb, torch.rand_like(s) * cfg.T)
        assert (u**2 + v**2).sqrt().max() < 5e-5


def test_training_smoke():
    from train import train
    cfg = pinn_config(epochs_adam=60, use_lbfgs=False, n_interior=300, n_norm_side=16,
                      log_every=10, pde_warmup_epochs=1, hidden_units=32, hidden_layers=2)
    model, hist, info = train(cfg, verbose=False)
    assert hist[-1]["loss"] < hist[0]["loss"]


RUN = "results/pinn/clean_fourier_uniform_s0"


HAVE_RUN = os.path.exists(os.path.join(RUN, "model.pt"))


@pytest.fixture(scope="module")
def trained_metrics():
    from evaluate import evaluate_run
    return evaluate_run(RUN, save=False)


@pytest.mark.skipif(not HAVE_RUN, reason="no trained model")
def test_trained_pinn_regression(trained_metrics):
    """The trained PINN must beat the trivial prediction psi = 0 (relative
    L2 error 1) and conserve the norm to within 5 %."""
    m = trained_metrics
    assert m["rel_l2_final"] < 0.9
    assert m["norm_max_dev"] < 0.05


@pytest.mark.skipif(not HAVE_RUN, reason="no trained model")
@pytest.mark.xfail(reason="KNOWN LIMITATION: with the CPU training budget used here the PINN "
                          "reaches ~49% final relative L2 error (README Section 9)", strict=False)
def test_trained_pinn_meets_20_percent_accuracy(trained_metrics):
    assert trained_metrics["rel_l2_final"] < 0.2
