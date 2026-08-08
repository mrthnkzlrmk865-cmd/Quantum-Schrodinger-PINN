"""
Evaluate a trained PINN against:
  (1) its own probability conservation over time (does integral(|Psi|^2)
      stay near 1, or did training collapse toward the trivial zero
      solution?), and
  (2) the classical DST split-step reference solution (reference_solver.py),
      via an L2 error curve over time.

Usage:
    python train.py                 # produces schrodinger_pinn.pt
    python reference_solver.py      # produces reference_solution.npz
    python evaluate.py              # produces figures/ + printed report
"""

import os

import numpy as np
import torch

from config import cfg
from model import SchrodingerPINN

FIG_DIR = "figures"


def load_model(path="schrodinger_pinn.pt"):
    model = SchrodingerPINN().to(cfg.device)
    model.load_state_dict(torch.load(path, map_location=cfg.device))
    model.eval()
    return model


def evaluate_on_reference(model, ref_path="reference_solution.npz"):
    data = np.load(ref_path)
    X, Y, frame_times = data["X"], data["Y"], data["frame_times"]
    psi_ref = data["psi_real"] + 1j * data["psi_imag"]  # (n_frames, n_grid, n_grid)

    n_grid = X.shape[0]
    dx = X[1, 0] - X[0, 0]

    x_flat = torch.tensor(X.reshape(-1, 1), dtype=torch.float32, device=cfg.device)
    y_flat = torch.tensor(Y.reshape(-1, 1), dtype=torch.float32, device=cfg.device)

    pinn_prob_total = []
    l2_errors = []
    rel_l2_errors = []

    with torch.no_grad():
        for i, t_val in enumerate(frame_times):
            t_flat = torch.full_like(x_flat, float(t_val))
            u, v = model(x_flat, y_flat, t_flat)
            u = u.cpu().numpy().reshape(n_grid, n_grid)
            v = v.cpu().numpy().reshape(n_grid, n_grid)

            prob_density = u**2 + v**2
            total_prob = float(np.sum(prob_density) * dx * dx)
            pinn_prob_total.append(total_prob)

            psi_pinn = u + 1j * v
            diff = psi_pinn - psi_ref[i]
            l2 = float(np.sqrt(np.sum(np.abs(diff) ** 2) * dx * dx))
            ref_norm = float(np.sqrt(np.sum(np.abs(psi_ref[i]) ** 2) * dx * dx))
            l2_errors.append(l2)
            rel_l2_errors.append(l2 / max(ref_norm, 1e-12))

    return {
        "frame_times": frame_times,
        "pinn_prob_total": np.array(pinn_prob_total),
        "l2_errors": np.array(l2_errors),
        "rel_l2_errors": np.array(rel_l2_errors),
    }


def make_plots(results):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(FIG_DIR, exist_ok=True)
    t = results["frame_times"]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t, results["pinn_prob_total"], label="PINN: integral |Psi|^2")
    ax.axhline(1.0, color="gray", linestyle="--", label="ideal (probability conserved)")
    ax.set_xlabel("t")
    ax.set_ylabel("Total probability")
    ax.set_title("Conservation check: total probability vs time")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "conservation_check.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(t, results["rel_l2_errors"])
    ax.set_xlabel("t")
    ax.set_ylabel("Relative L2 error vs. reference")
    ax.set_title("PINN vs. classical (DST split-step) reference solution")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG_DIR, "l2_error_vs_reference.png"), dpi=150)
    plt.close(fig)


def report(results):
    prob = results["pinn_prob_total"]
    rel_err = results["rel_l2_errors"]

    print("=== Conservation check ===")
    print(f"  total probability: min={prob.min():.4f} max={prob.max():.4f} "
          f"(ideal = 1.0 at every t)")
    print(f"  max deviation from 1: {np.max(np.abs(prob - 1.0)):.4f}")
    if np.max(np.abs(prob - 1.0)) > 0.2:
        print("  WARNING: probability is not well conserved. The model may be "
              "under-trained, or collapsing toward / diverging from the trivial "
              "zero solution. Consider increasing loss_weights['norm'] or "
              "epochs, and re-checking the loss curve.")

    print("\n=== Accuracy vs. classical reference (DST split-step) ===")
    print(f"  relative L2 error: mean={rel_err.mean():.4f} "
          f"final(t=T)={rel_err[-1]:.4f} max={rel_err.max():.4f}")

    print(f"\nFigures written to {FIG_DIR}/")


if __name__ == "__main__":
    model = load_model()
    results = evaluate_on_reference(model)
    report(results)
    make_plots(results)
