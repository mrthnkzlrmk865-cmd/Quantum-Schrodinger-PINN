"""
Original-contribution experiment: quantum transport through a barrier with
smooth, correlated spatial disorder, compared to the clean barrier.

Motivation
-----------
This is the piece that moves the project beyond "PINN applied to a textbook
tunneling problem" (see README "Scope and original contribution"). It asks
a question inspired by Environment-Assisted Quantum Transport (ENAQT) in
photosynthetic light-harvesting complexes: static energetic disorder along
a transport pathway can, counter-intuitively, *increase* transmission
through constructive interference effects, rather than only hindering it
via Anderson localization. Most ENAQT literature studies this with open
quantum system master equations (Lindblad/HEOM) and a *dephasing*
environment; this experiment asks a narrower, complementary question with
a much simpler, fully unitary (Schrodinger-only) model: does adding smooth,
correlated *static* disorder around a tunneling barrier change transmission,
and does it do so systematically (helps/hurts) or unpredictably
(seed-dependent)?

This uses the classical DST split-step solver (reference_solver.py), not
the PINN -- the PINN's job in this project is to reproduce known unitary
dynamics accurately; the disorder *sweep* itself is intentionally run with
the solver already validated for probability conservation to machine
precision (see tests/test_basic.py), so any transmission trend we see here
is physical, not a training artifact.

Usage:
    python experiment_disorder.py
"""

import numpy as np

from config import cfg
import reference_solver as rs


def transmitted_probability(psi: np.ndarray, X: np.ndarray, dx: float) -> float:
    """Probability mass to the right of the barrier center (x > barrier_center)."""
    mask = X > cfg.barrier_center
    return float(np.sum(np.abs(psi[mask]) ** 2) * dx * dx)


def run_scenario(scenario: str, seed: int = 0):
    cfg.SCENARIO = scenario
    cfg.disorder_seed = seed
    X, Y, frame_times, frames = rs.propagate()
    dx = X[1, 0] - X[0, 0]
    trans = [transmitted_probability(p, X, dx) for p in frames]
    total_prob = [rs.total_probability(p, dx) for p in frames]
    return frame_times, np.array(trans), np.array(total_prob)


def run_disorder_sweep(n_seeds: int = 8):
    """Compare final transmission for the clean barrier vs. several
    independent disorder realizations, to see if disorder systematically
    helps/hurts or is just seed-noise."""
    cfg.eval_n_grid = 96
    cfg.eval_n_steps = 3000

    t_clean, trans_clean, prob_clean = run_scenario("tunneling")
    final_clean = trans_clean[-1]

    finals_disordered = []
    for seed in range(n_seeds):
        _, trans_d, prob_d = run_scenario("disordered", seed=seed)
        finals_disordered.append(trans_d[-1])
        drift = abs(prob_d[-1] - 1.0)
        assert drift < 1e-6, f"seed {seed}: probability not conserved (drift={drift})"

    finals_disordered = np.array(finals_disordered)
    return final_clean, finals_disordered, t_clean, trans_clean


def make_plot(final_clean, finals_disordered, t_clean, trans_clean):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    ax = axes[0]
    ax.axhline(final_clean, color="black", linestyle="--", label="clean barrier")
    ax.scatter(
        range(len(finals_disordered)), finals_disordered,
        color="tab:orange", label="disordered (per seed)"
    )
    ax.axhline(
        finals_disordered.mean(), color="tab:orange", linestyle=":",
        label=f"disordered mean ({finals_disordered.mean():.3f})"
    )
    ax.set_xlabel("disorder realization (seed)")
    ax.set_ylabel("Transmitted probability at t=T")
    ax.set_title("Final transmission: clean vs. disordered")
    ax.legend(fontsize=8)

    ax = axes[1]
    cfg.SCENARIO = "tunneling"
    t_d, trans_d0, _ = run_scenario("disordered", seed=0)
    ax.plot(t_clean, trans_clean, label="clean barrier")
    ax.plot(t_d, trans_d0, label="disordered (seed 0)")
    ax.set_xlabel("t")
    ax.set_ylabel("Transmitted probability")
    ax.set_title("Transmission over time (one realization)")
    ax.legend(fontsize=8)

    fig.tight_layout()
    fig.savefig("figures/disorder_transmission.png", dpi=150)
    print("Saved -> figures/disorder_transmission.png")


if __name__ == "__main__":
    import os

    os.makedirs("figures", exist_ok=True)

    print("Running clean-vs-disordered transmission sweep "
          "(this uses the classical reference solver, not the PINN)...")
    final_clean, finals_disordered, t_clean, trans_clean = run_disorder_sweep(n_seeds=8)

    print(f"\nClean barrier, final transmitted probability: {final_clean:.4f}")
    print(f"Disordered barrier, {len(finals_disordered)} seeds: "
          f"mean={finals_disordered.mean():.4f} std={finals_disordered.std():.4f}")
    print(f"Per-seed values: {np.round(finals_disordered, 4).tolist()}")

    delta = finals_disordered.mean() - final_clean
    direction = "increases" if delta > 0 else "decreases"
    print(f"\nOn average, this disorder configuration {direction} transmission "
          f"by {abs(delta):.4f} ({abs(delta)/final_clean*100:.1f}% relative to clean), "
          f"with realization-to-realization std {finals_disordered.std():.4f}.")
    if finals_disordered.std() > abs(delta):
        print("NOTE: the seed-to-seed spread is larger than the mean shift -- "
              "read this as 'disorder matters, direction is not yet a clean, "
              "statistically confident trend' rather than a firm effect. "
              "More seeds and/or a scan over disorder_strength would be the "
              "natural next step before claiming a directional ENAQT-like effect.")

    make_plot(final_clean, finals_disordered, t_clean, trans_clean)
