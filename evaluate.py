"""
Evaluate a trained PINN against the DST split-step reference solution of
the SAME configuration (same potential, same disorder realization, same
initial condition), on the reference grid.

Metrics at each of cfg.n_frames times:
  rel_l2       ||psi_PINN - psi_ref||_2 / ||psi_ref||_2   (complex wavefunction;
               sensitive to the global phase)
  rel_l2_rho   same for the probability density rho = |psi|^2
               (insensitive to phase errors)
  norm_pinn    int |psi_PINN|^2      (ideal: 1)
  com_x_err    |<x>_PINN - <x>_ref|  (centre-of-mass error)
  p_right_err  |P(x > 0.5)_PINN - P(x > 0.5)_ref|

Reference: the PINN error is meaningful only if the reference error is much
smaller. validate_solver.py shows the reference for the PINN problem
(127^2, dt = 1e-5) has relative L2 error ~1e-8 vs the analytic free packet
and ~1e-6 dt-self-convergence error with disorder.

Usage:
    python evaluate.py results/pinn/<run_dir>
"""

from __future__ import annotations

import argparse
import csv
import json
import os

import numpy as np
import torch

from observables import centre_of_mass, region_probability
from potential import Potential
from reference_solver import SplitStepDST


def pinn_on_grid(model, X, Y, times, device="cpu"):
    x = torch.tensor(X.reshape(-1, 1), dtype=torch.float32, device=device)
    y = torch.tensor(Y.reshape(-1, 1), dtype=torch.float32, device=device)
    out = []
    with torch.no_grad():
        for t in times:
            u, v = model(x, y, torch.full_like(x, float(t)))
            out.append((u + 1j * v).cpu().numpy().reshape(X.shape).astype(np.complex128))
    return np.array(out)


def compare(cfg, psi_pinn, psi_ref, X, Y, times):
    dx, dy = X[1, 0] - X[0, 0], Y[0, 1] - Y[0, 0]
    rows = []
    mask = X > 0.5 * cfg.Lx
    for i, t in enumerate(times):
        p, r = psi_pinn[i], psi_ref[i]
        rho_p, rho_r = np.abs(p) ** 2, np.abs(r) ** 2
        cp = centre_of_mass(p, X, Y, dx, dy)[0]
        cr = centre_of_mass(r, X, Y, dx, dy)[0]
        rows.append(dict(
            t=float(t),
            rel_l2=float(np.linalg.norm(p - r) / np.linalg.norm(r)),
            rel_l2_rho=float(np.linalg.norm(rho_p - rho_r) / np.linalg.norm(rho_r)),
            norm_pinn=float(rho_p.sum() * dx * dy),
            norm_ref=float(rho_r.sum() * dx * dy),
            com_x_err=float(abs(cp - cr)),
            p_right_err=float(abs(region_probability(p, mask, dx, dy) - region_probability(r, mask, dx, dy))),
        ))
    return rows


def evaluate_run(run_dir, nx=None, save=True):
    from train import load_run

    cfg, model = load_run(run_dir)
    if nx is not None:
        cfg = cfg.replace(nx=nx, ny=nx)
    solver = SplitStepDST(cfg, Potential(cfg))
    ref = solver.run()
    times = ref["times"]
    psi_ref = ref["frames"][:, 0]
    psi_pinn = pinn_on_grid(model, solver.X, solver.Y, times)
    rows = compare(cfg, psi_pinn, psi_ref, solver.X, solver.Y, times)
    rel = np.array([r["rel_l2"] for r in rows])
    rel_rho = np.array([r["rel_l2_rho"] for r in rows])
    norms = np.array([r["norm_pinn"] for r in rows])
    summary = dict(
        rel_l2_mean=float(rel.mean()), rel_l2_final=float(rel[-1]), rel_l2_max=float(rel.max()),
        rel_l2_rho_mean=float(rel_rho.mean()), rel_l2_rho_final=float(rel_rho[-1]),
        norm_max_dev=float(np.max(np.abs(norms - 1.0))),
        com_x_err_max=float(max(r["com_x_err"] for r in rows)),
        p_right_err_max=float(max(r["p_right_err"] for r in rows)),
    )
    if save:
        with open(os.path.join(run_dir, "eval_timeseries.csv"), "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)
        with open(os.path.join(run_dir, "eval_summary.json"), "w") as f:
            json.dump(summary, f, indent=2)
        np.savez_compressed(os.path.join(run_dir, "eval_fields.npz"), X=solver.X, Y=solver.Y,
                            times=times, psi_pinn=psi_pinn[::5], psi_ref=psi_ref[::5],
                            V=solver.V[0])
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    a = ap.parse_args()
    s = evaluate_run(a.run_dir)
    for k, v in s.items():
        print(f"{k:>18s}: {v:.4e}")
