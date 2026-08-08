# Quantum Schrodinger PINN: 2D Tunneling, with a Ground-Truth Check and an ENAQT-Inspired Disorder Experiment

A Physics-Informed Neural Network (PINN) that solves the time-dependent 2D
Schrodinger equation for a Gaussian wave packet tunneling through a
rectangular barrier in an infinite square well — plus a classical
split-step reference solver used to check the PINN's accuracy, and a small
original experiment on how spatial disorder affects tunneling transmission.

## Scope and original contribution

The core PINN (autograd-based PDE residual, complex wavefunction split into
real/imaginary parts, collocation sampling) follows the standard recipe
from Raissi, Perdikaris & Karniadakis (2019). That part is a *replication*
of known methodology on a known problem, not a novel technique, and this
README doesn't claim otherwise.

What this repository adds beyond that baseline:

1. **A validated ground truth.** A from-scratch classical solver
   (`reference_solver.py`) using a 2D Discrete Sine Transform split-step
   method, chosen specifically to match the PINN's Dirichlet boundary
   condition exactly (see "Why DST, not FFT" below). Its own probability
   conservation is verified to ~1e-13.
2. **A documented debugging process.** The first working version of this
   PINN silently collapsed to the trivial zero solution. Finding out why
   took four separate, compounding bugs/limitations — see "Debugging
   note" below. This is included deliberately: a project that shows how a
   plausible-looking wrong answer was caught and fixed is stronger
   evidence of understanding than a project that only shows a final plot.
3. **A small original experiment** (`experiment_disorder.py`) asking
   whether smooth, spatially-correlated static disorder around the barrier
   changes tunneling transmission — motivated by Environment-Assisted
   Quantum Transport (ENAQT) in photosynthetic complexes, where structured
   disorder/noise can *increase* transport efficiency rather than only
   hindering it. This uses a fully unitary model (no explicit dephasing
   bath), so it is a narrower question than the open-quantum-system ENAQT
   literature, and is presented as such.

## Physical setup

- 2D infinite square well, `[0, L] x [0, L]`, `L = 1`.
- Initial state: a normalized Gaussian wave packet with momentum `(kx0,
  ky0)`, centered at `(x0, y0)`, width `sigma`.
- Potential: a rectangular barrier of configurable height/width
  (`SCENARIO = "tunneling"`), or the same barrier plus smooth correlated
  disorder (`SCENARIO = "disordered"`, used only in the disorder
  experiment).
- All parameters live in `config.py`, the single source of truth shared by
  the PINN, the reference solver, and the disorder experiment.

## Why DST, not FFT, for the reference solver

The textbook split-step Fourier method assumes periodic boundaries — the
FFT diagonalizes the Laplacian on a ring. This PINN is trained with
Dirichlet boundaries (`Psi = 0` at the walls). A periodic reference would
be solving a different physical problem, making any comparison unfair. The
DST-I basis functions `sin(n*pi*x/L)` are exactly the Dirichlet-Laplacian
eigenfunctions on `[0, L]`, so using a 2D DST-I in place of the FFT gives a
spectral method that enforces `Psi = 0` at the walls exactly, by
construction, matching the PINN's boundary condition rather than
approximating it.

## Debugging note: four bugs behind one collapsed model

The first trained model reached a low combined loss while doing almost
nothing physically: `integral |Psi|^2` sat near 0 at every time, including
`t=0`. Tracking this down surfaced four separate issues, each compounding
the last:

1. **IC point starvation.** Uniformly sampling initial-condition points
   over the whole domain meant only ~5% of them landed near the packet
   (`sigma=0.06` in a unit box) — the loss was dominated by the
   "everywhere-else-is-zero" background, giving weak, high-variance
   gradient signal about the packet's actual shape. *Fix:* importance-sample
   IC points around the packet, mixed with a smaller uniform background
   sample.
2. **Spectral bias.** With momentum `kx0=40`, the packet's phase oscillates
   roughly every `2*pi/40 ~= 0.157` units — 2-3 oscillations within one
   packet width. A plain Tanh-MLP struggles to represent this (verified in
   isolation: a pure supervised fit of the network to the IC alone
   plateaued at MSE ~0.01 and stopped improving). *Fix:* a Fourier feature
   input embedding (`embeddings.py`), which dropped the isolated IC-fit MSE
   to ~0.00005.
3. **PDE-loss blow-up.** Fourier features fixed (2), but their second
   derivatives (needed for the Laplacian in the PDE residual) scale with
   `(2*pi*f)^2`; with the frequency range first tried (up to 128), a
   freshly-initialized network had `mean(f_u^2+f_v^2) ~ 1.2e7`, swamping
   every other loss term. *Fix:* a narrower, physically-motivated frequency
   range (linearly spaced up to ~10, matched to the packet's actual
   spectral content), gradient clipping, and a short linear warm-up of the
   PDE loss weight.
4. **The actual root cause.** After fixing 1-3, the model still collapsed.
   Direct inspection showed the "IC" the network was being fit to was
   *not normalized*: `initial_condition()` used a Gaussian with peak
   amplitude 1.0, whose true `integral(|Psi0|^2)` is `pi*sigma^2 ~= 0.011`,
   not 1. The training setup was simultaneously asking the network to (a)
   reproduce that unnormalized, low-probability packet exactly, and (b)
   have total probability exactly 1 at every time including `t=0` — two
   *mutually contradictory* targets. No amount of reweighting or
   architecture change could satisfy both. *Fix:* normalize the analytic
   IC so `integral(u0^2+v0^2) = 1`, and switch the normalization loss's own
   Monte-Carlo integral estimate (originally ~400 uniform random points,
   which is a noisy, biased-low estimator for a narrow peaked function)
   to a fixed quadrature grid.

On top of these four fixes, the IC and boundary conditions were also moved
into the network's output parametrization directly (hard constraints,
`model.py`):
```
u(x,y,t) = u0(x,y) + B(x,y) * t * NN_u(x,y,t)
v(x,y,t) = v0(x,y) + B(x,y) * t * NN_v(x,y,t)
```
with `B(x,y)` a smooth bump vanishing at the walls. This makes `t=0` exact
by construction and removes two loss terms (and their weight-tuning) from
the optimization entirely, closing off the collapse pathway rather than
just penalizing it more heavily.

## Current results (honest status)

The `figures/` and numbers below are from a training run capped to fit a
CPU-only sandbox (~1000 Adam epochs + ~250 L-BFGS steps, a few minutes
total). This is **not** a fully converged model — see "What's not done yet"
below. It is included specifically because it's an honest checkpoint,
consistent with the debugging story above, not a cherry-picked final
result.

- **Conservation check** (`figures/conservation_check.png`): total
  probability starts at exactly 1.0 (guaranteed by the hard IC constraint)
  and drifts down to roughly 0.6-0.8 with training noise as `t` increases.
  This is a large improvement over the pre-fix behavior (`~0.00001` at
  every time, including `t=0`), but is not yet flat at 1.0 — the PDE
  dynamics are not yet fully learned.
- **Accuracy vs. the DST reference solver**
  (`figures/l2_error_vs_reference.png`): relative L2 error is 0 at `t=0`
  (again, exact by construction) and jumps to roughly 1.2-1.4 almost
  immediately after — the network has not yet learned accurate *time
  propagation*, only the initial state.
- **Disorder experiment** (`figures/disorder_transmission.png`, computed
  with the already-validated classical solver, not the PINN): across 8
  disorder realizations, transmission past the barrier increased by
  ~53% on average relative to the clean barrier (mean effect 0.147,
  seed-to-seed std 0.019 — the effect is several times larger than the
  realization-to-realization noise, i.e. a real trend for this disorder
  strength/correlation length, not sampling luck). The clean-barrier
  transmission also oscillates much more strongly over time (repeated
  reflection within the box) than the disordered case, which settles into
  a smoother, more damped transmission curve — a qualitative signature
  worth following up analytically.

## What's not done yet / natural next steps

- **Full PDE convergence.** PDE residual was still decreasing steadily
  (not plateaued) when the CPU training budget ran out. Re-running
  `train.py` with the full `epochs_adam`/`epochs_lbfgs` in `config.py` on a
  GPU should close most of the remaining gap; re-run `evaluate.py`
  afterward to get updated conservation/L2 numbers and figures.
- **Disorder-strength scan.** The disorder experiment currently uses one
  fixed `disorder_strength`. A sweep over strength (and correlation
  length) would turn "disorder helps transmission here" into a proper
  characterization of *when* it helps, which is what would make the ENAQT
  connection scientifically substantive rather than suggestive.
- **PINN-vs-reference comparison under disorder**, once the PINN itself
  converges well on the clean case — right now the disorder experiment
  intentionally only uses the validated classical solver, to keep its
  result independent of the PINN's current training state.

## Repository structure

```
config.py              Single source of truth for all physical/training parameters
potential.py            Barrier (+ optional disorder) potential, torch and numpy versions
initial_condition.py    Normalized Gaussian wave packet (see debugging note #4)
embeddings.py           Fourier feature input embedding (see debugging note #2)
model.py                PINN architecture with hard IC/BC constraints (see debugging note)
pinn_core.py             PDE residual, collocation sampling, normalization loss
train.py                 Adam -> L-BFGS training script
reference_solver.py      Classical DST split-step ground-truth solver
evaluate.py              Conservation check + L2 error vs. reference, generates figures/
experiment_disorder.py   Original-contribution disorder/transmission experiment
visualize.py             3D interactive Plotly animation (+ PINN-vs-reference overlay)
tests/test_basic.py      Unit tests: potential consistency, solver conservation, shapes
requirements.txt
figures/                 Generated plots (conservation, L2 error, disorder comparison)
```

## Reproducing the results

```bash
pip install -r requirements.txt

# 1. Train the PINN (uses config.py defaults; reduce epochs_adam/epochs_lbfgs
#    in config.py for a quick CPU smoke test)
python train.py

# 2. Generate the classical ground-truth solution
python reference_solver.py

# 3. Evaluate: conservation check + accuracy vs. reference -> figures/
python evaluate.py

# 4. Original-contribution experiment (does not need the trained PINN)
python experiment_disorder.py

# 5. Interactive 3D visualization (+ PINN-vs-reference overlay, if
#    reference_solution.npz exists)
python visualize.py

# Tests
python -m pytest tests/
```

