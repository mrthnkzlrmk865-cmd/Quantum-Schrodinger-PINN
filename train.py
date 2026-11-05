"""
Train the Schrodinger PINN: Adam (with a linear warm-up of the PDE weight
and step learning-rate decay), then optionally L-BFGS (Liu & Nocedal 1989)
on a fixed set of collocation points, following the usual PINN recipe
(Raissi et al. 2019).

Usage:
    python train.py                                  # clean potential
    python train.py --strength 56 --zeta 0.08 --seed 0 --out results/pinn/run
Outputs (in --out directory): model.pt, config.json, history.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time

import numpy as np
import torch

from config import Config, pinn_config
from model import SchrodingerPINN
from pinn_core import norm_quadrature_points, sample_interior, total_loss
from potential import Potential


def set_seed(seed: int):
    torch.manual_seed(seed)
    np.random.seed(seed)


def train(cfg: Config, potential: Potential | None = None, verbose=True):
    set_seed(cfg.seed)
    gen = torch.Generator(device=cfg.device).manual_seed(cfg.seed)
    potential = potential or Potential(cfg)
    model = SchrodingerPINN(cfg).to(cfg.device)
    xq, yq, w = norm_quadrature_points(cfg, cfg.device)

    def norm_pts():
        times = torch.rand(cfg.n_norm_times, generator=gen, device=cfg.device) * cfg.T
        return xq, yq, w, times

    history = []
    t_start = time.time()

    opt = torch.optim.Adam(model.parameters(), lr=cfg.lr_adam)
    sched = torch.optim.lr_scheduler.StepLR(opt, step_size=cfg.lr_decay_every, gamma=cfg.lr_decay_gamma)
    interior = sample_interior(cfg, model, potential, gen)
    npts = norm_pts()
    for epoch in range(1, cfg.epochs_adam + 1):
        weights = dict(cfg.loss_weights)
        weights["pde"] *= min(1.0, epoch / max(1, cfg.pde_warmup_epochs))
        opt.zero_grad()
        loss, parts = total_loss(model, potential, interior, npts, weights, cfg.causal_eps, cfg.causal_bins)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), cfg.grad_clip_norm)
        opt.step()
        sched.step()
        if epoch % cfg.log_every == 0 or epoch == 1:
            history.append(dict(phase="adam", epoch=epoch, loss=loss.item(), **parts,
                                wall_time=time.time() - t_start))
            if verbose:
                print(f"[Adam {epoch:5d}] loss={loss.item():.3e} pde={parts['pde']:.3e} "
                      f"norm={parts['norm']:.3e} wmin={parts['causal_wmin']:.2f} "
                      f"({time.time()-t_start:.0f}s)", flush=True)
        every = cfg.rad_resample_every if cfg.sampling == "rad" else cfg.resample_every
        if epoch % every == 0 and epoch < cfg.epochs_adam:
            interior = sample_interior(cfg, model, potential, gen)
            npts = norm_pts()
    t_adam = time.time() - t_start

    if cfg.use_lbfgs and cfg.epochs_lbfgs > 0:
        interior = sample_interior(cfg, model, potential, gen, n=cfg.lbfgs_points_factor * cfg.n_interior)
        npts = norm_pts()
        opt = torch.optim.LBFGS(model.parameters(), lr=1.0, max_iter=cfg.epochs_lbfgs,
                                history_size=50, line_search_fn="strong_wolfe",
                                tolerance_grad=1e-9, tolerance_change=1e-12)
        it = {"n": 0}

        def closure():
            opt.zero_grad()
            loss, parts = total_loss(model, potential, interior, npts, cfg.loss_weights,
                                     cfg.causal_eps, cfg.causal_bins)
            loss.backward()
            it["n"] += 1
            if it["n"] % 50 == 0:
                history.append(dict(phase="lbfgs", epoch=cfg.epochs_adam + it["n"], loss=loss.item(),
                                    **parts, wall_time=time.time() - t_start))
                if verbose:
                    print(f"[L-BFGS {it['n']:4d}] loss={loss.item():.3e} pde={parts['pde']:.3e} "
                          f"norm={parts['norm']:.3e}", flush=True)
            return loss

        opt.step(closure)
    total_time = time.time() - t_start
    # generalization check: residual on fresh uniform points never used in training
    hold = torch.Generator(device=cfg.device).manual_seed(10_000 + cfg.seed)
    from pinn_core import _uniform, pde_residual
    xh, yh, th = _uniform(cfg, 20000, hold, cfg.device)
    fu, fv = pde_residual(model, potential, xh, yh, th)
    info = dict(train_time_s=total_time, adam_time_s=t_adam,
                final_loss=history[-1]["loss"] if history else None,
                final_train_pde=history[-1]["pde"] if history else None,
                holdout_pde=float(((fu**2).mean() + (fv**2).mean()).detach()))
    return model, history, info


def save_run(out_dir, cfg, model, history, info):
    os.makedirs(out_dir, exist_ok=True)
    torch.save(model.state_dict(), os.path.join(out_dir, "model.pt"))
    cfg.save_json(os.path.join(out_dir, "config.json"))
    with open(os.path.join(out_dir, "train_info.json"), "w") as f:
        json.dump(info, f, indent=2)
    if history:
        with open(os.path.join(out_dir, "history.csv"), "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(history[0].keys()))
            wr.writeheader()
            wr.writerows(history)


def load_run(out_dir, device="cpu"):
    cfg = Config.load_json(os.path.join(out_dir, "config.json")).replace(device=device)
    model = SchrodingerPINN(cfg).to(device)
    model.load_state_dict(torch.load(os.path.join(out_dir, "model.pt"), map_location=device))
    model.eval()
    return cfg, model


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--strength", type=float, default=0.0)
    ap.add_argument("--zeta", type=float, default=0.08)
    ap.add_argument("--disorder-seed", type=int, default=0)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--epochs-adam", type=int, default=None)
    ap.add_argument("--epochs-lbfgs", type=int, default=None)
    ap.add_argument("--no-fourier", action="store_true")
    ap.add_argument("--sampling", default="uniform")
    ap.add_argument("--T", type=float, default=None)
    ap.add_argument("--causal-eps", type=float, default=0.0)
    ap.add_argument("--no-lbfgs", action="store_true")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--out", default="results/pinn/default")
    a = ap.parse_args()
    over = dict(disorder_strength=a.strength, disorder_corr_length=a.zeta,
                disorder_seed=a.disorder_seed, seed=a.seed,
                use_fourier=not a.no_fourier, sampling=a.sampling)
    if a.epochs_adam is not None:
        over["epochs_adam"] = a.epochs_adam
    if a.epochs_lbfgs is not None:
        over["epochs_lbfgs"] = a.epochs_lbfgs
    if a.T is not None:
        over["T"] = a.T
    over["causal_eps"] = a.causal_eps
    if a.no_lbfgs:
        over["use_lbfgs"] = False
    cfg = pinn_config(**over)
    torch.set_num_threads(a.threads or max(1, os.cpu_count() or 1))
    model, hist, info = train(cfg)
    save_run(a.out, cfg, model, hist, info)
    print(f"saved -> {a.out}  ({info['train_time_s']:.0f}s)")
