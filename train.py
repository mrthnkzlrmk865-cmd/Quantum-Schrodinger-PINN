"""
Train the Schrodinger PINN.

Two-phase optimization, following the standard PINN recipe (Raissi et al.,
2019): Adam for fast, robust early progress, then L-BFGS for high-precision
convergence once the loss landscape is closer to a local basin. The original
draft only ran Adam, which typically plateaus at a higher residual than a
PINN with an L-BFGS finishing phase.

Usage:
    python train.py
"""

import time

import torch

from config import cfg
from model import SchrodingerPINN
from pinn_core import sample_points, total_loss

torch.manual_seed(cfg.seed)


def train_adam(model, optimizer, scheduler):
    interior, ic, bc, norm_pts = sample_points()
    for epoch in range(1, cfg.epochs_adam + 1):
        # Linear warm-up of the PDE loss weight (see config.py docstring).
        warmup_frac = min(1.0, epoch / max(1, cfg.pde_warmup_epochs))
        weights = dict(cfg.loss_weights)
        weights["pde"] = cfg.loss_weights["pde"] * warmup_frac

        optimizer.zero_grad()
        loss, parts = total_loss(model, interior, ic, bc, norm_pts, weights=weights)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), cfg.grad_clip_norm)
        optimizer.step()
        scheduler.step()

        if epoch % cfg.resample_every == 0:
            interior, ic, bc, norm_pts = sample_points()

        if epoch % 200 == 0 or epoch == 1:
            print(
                f"[Adam {epoch:6d}/{cfg.epochs_adam}] total={loss.item():.6f} "
                f"pde={parts['pde']:.6f} ic={parts['ic']:.6f} "
                f"bc={parts['bc']:.6f} norm={parts['norm']:.6f} "
                f"(pde_w={weights['pde']:.3f})"
            )
    return model


def train_lbfgs(model):
    interior, ic, bc, norm_pts = sample_points()

    optimizer = torch.optim.LBFGS(
        model.parameters(),
        lr=1.0,
        max_iter=cfg.epochs_lbfgs,
        history_size=50,
        line_search_fn="strong_wolfe",
    )

    iteration = {"n": 0}

    def closure():
        optimizer.zero_grad()
        loss, parts = total_loss(model, interior, ic, bc, norm_pts)
        loss.backward()
        iteration["n"] += 1
        if iteration["n"] % 50 == 0:
            print(
                f"[L-BFGS {iteration['n']:5d}] total={loss.item():.6f} "
                f"pde={parts['pde']:.6f} ic={parts['ic']:.6f} "
                f"bc={parts['bc']:.6f} norm={parts['norm']:.6f}"
            )
        return loss

    optimizer.step(closure)
    return model


def train():
    model = SchrodingerPINN().to(cfg.device)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.lr_adam)
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=1000, gamma=0.5)

    t0 = time.time()
    model = train_adam(model, optimizer, scheduler)
    print(f"Adam phase done in {time.time() - t0:.1f}s")

    if cfg.use_lbfgs:
        t1 = time.time()
        model = train_lbfgs(model)
        print(f"L-BFGS phase done in {time.time() - t1:.1f}s")

    torch.save(model.state_dict(), "schrodinger_pinn.pt")
    print("Model saved -> schrodinger_pinn.pt")
    return model


if __name__ == "__main__":
    print(f"Device: {cfg.device} | Scenario: {cfg.SCENARIO}")
    train()
