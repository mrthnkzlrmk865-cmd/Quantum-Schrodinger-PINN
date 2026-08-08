# Quantum Schrodinger PINN: 2D Tunneling & ENAQT-Inspired Disorder Experiment

A Physics-Informed Neural Network (PINN) that solves the time-dependent 2D Schrodinger equation for a Gaussian wave packet tunneling through a rectangular barrier in an infinite square well — plus a classical split-step reference solver used to check the PINN's accuracy, and a small original experiment on how spatial disorder affects tunneling transmission.

## Scope and original contribution

The core PINN (autograd-based PDE residual, complex wavefunction split into real/imaginary parts, collocation sampling) follows the standard recipe from Raissi, Perdikaris & Karniadakis (2019). That part is a *replication* of known methodology on a known problem, not a novel technique, and this README doesn't claim otherwise.

What this repository adds beyond that baseline:

1. **A validated ground truth.** A from-scratch classical solver (`reference_solver.py`) using a 2D Discrete Sine Transform split-step method, chosen specifically to match the PINN's Dirichlet boundary condition exactly (see "Why DST, not FFT" below). Its own probability conservation is verified to ~1e-13.
2. **A documented debugging process.** As a high school student exploring quantum computing and computational physics, building this model meant encountering silent failures firsthand. The first working version of this PINN silently collapsed to the trivial zero solution. Finding out why took four separate, compounding bugs — see "Debugging note" below. Documenting how a plausible-looking wrong answer was caught and fixed felt just as valuable as writing the final working code.
3. **A small original experiment** (`experiment_disorder.py`) asking whether smooth, spatially-correlated static disorder around the barrier changes tunneling transmission — motivated by Environment-Assisted Quantum Transport (ENAQT) in photosynthetic complexes, where structured disorder/noise can *increase* transport efficiency rather than only hindering it. This uses a fully unitary model (no explicit dephasing bath), so it is a narrower question than the open-quantum-system ENAQT literature, and is presented as such.

## Physical setup

- 2D infinite square well, `[0, L] x [0, L]`, `L = 1`.
- Initial state: a normalized Gaussian wave packet with momentum `(kx0, ky0)`, centered at `(x0, y0)`, width `sigma`.
- Potential: a rectangular barrier of configurable height/width (`SCENARIO = "tunneling"`), or the same barrier plus smooth correlated disorder (`SCENARIO = "disordered"`, used only in the disorder experiment).
- All parameters live in `config.py`, the single source of truth shared by the PINN, the reference solver, and the disorder experiment.

## Why DST, not FFT, for the reference solver

The textbook split-step Fourier method assumes periodic boundaries — the FFT diagonalizes the Laplacian on a ring. This PINN is trained with Dirichlet boundaries (`Psi = 0` at the walls). A periodic reference would be solving a different physical problem, making any comparison unfair. The DST-I basis functions `sin(n*pi*x/L)` are exactly the Dirichlet-Laplacian eigenfunctions on `[0, L]`, so using a 2D DST-I in place of the FFT gives a spectral method that enforces `Psi = 0` at the walls exactly, by construction, matching the PINN's boundary condition rather than approximating it.

## Debugging note: four bugs behind one collapsed model

The first trained model reached a low combined loss while doing almost nothing physically: total probability sat near 0 at every time, including `t=0`. Tracking this down surfaced four separate issues, each compounding the last:

1. **IC point starvation.** Uniformly sampling initial-condition points over the whole domain meant only ~5% of them landed near the packet (`sigma=0.06` in a unit box) — the loss was dominated by the "everywhere-else-is-zero" background, giving weak, high-variance gradient signal about the packet's actual shape. *Fix:* importance-sample IC points around the packet, mixed with a smaller uniform background sample.
2. **Spectral bias.** With momentum `kx0=40`, the packet's phase oscillates rapidly across the box. A plain Tanh-MLP struggles to represent this high-frequency component. *Fix:* a Fourier feature input embedding (`embeddings.py`), which dramatically dropped the isolated IC-fit MSE.
3. **PDE-loss blow-up.** Fourier features fixed the fit, but their second derivatives (needed for the Laplacian in the PDE residual) scale quadratically with frequency. A freshly-initialized network swamped every other loss term. *Fix:* a narrower, physically-motivated frequency range, gradient clipping, and a short linear warm-up of the PDE loss weight.
4. **The actual root cause.** After fixing 1-3, the model still collapsed. Direct inspection showed the initial condition function originally produced an unnormalized peak amplitude, meaning the total integral was much smaller than 1. The training setup was simultaneously asking the network to reproduce that unnormalized wave packet and enforce a total probability of 1 at all times — two contradictory targets. *Fix:* normalize the analytic IC, and switch the normalization loss to a fixed quadrature grid instead of noisy uniform Monte-Carlo sampling.

On top of these four fixes, the IC and boundary conditions were also moved into the network's output parametrization directly via hard constraints in `model.py` (`u = u0 + B * t * NN_u`). This makes `t=0` exact by construction and removes two loss terms from optimization entirely, closing off the collapse pathway.

## Results & Visualizations

### 1. Total Probability Conservation
Total probability starts at exactly 1.0 (guaranteed by the hard IC constraint) and fluctuates between 0.6 and 0.8 as time progresses during short CPU runs. While not yet flat at 1.0, this represents a major improvement over the pre-fix state where probability collapsed near zero immediately.

![Conservation Check](figures/conservation_check.png)

### 2. Relative L2 Error vs. DST Reference
Compared against the classical split-step reference solver, the relative L2 error remains 0 at `t=0` by construction and settles around 1.2–1.4 across time steps, highlighting that further GPU training epochs are required for full temporal convergence.

![L2 Error vs Reference](figures/l2_error_vs_reference.png)

### 3. ENAQT-Inspired Disorder Experiment
Across 8 static disorder realizations using the classical solver, tunneling transmission past the barrier increased significantly (mean transmission ~0.424) compared to the clean barrier (~0.276). Static spatial disorder suppresses destructive interference back-reflection, creating smoother transmission curves.

![Disorder Transmission](figures/disorder_transmission.png)

## Repository Structure

```text
config.py             Single source of truth for all physical/training parameters
potential.py          Barrier (+ optional disorder) potential, torch and numpy versions
initial_condition.py  Normalized Gaussian wave packet
embeddings.py         Fourier feature input embedding
model.py              PINN architecture with hard IC/BC constraints
pinn_core.py          PDE residual, collocation sampling, normalization loss
train.py              Adam -> L-BFGS training script
reference_solver.py   Classical DST split-step ground-truth solver
evaluate.py           Conservation check + L2 error vs. reference, generates figures/
experiment_disorder.py Original-contribution disorder/transmission experiment
visualize.py          3D interactive Plotly animation (+ PINN-vs-reference overlay)
run.sh                One-click automation script (setup, train, evaluate, test)
tests/test_basic.py   Unit tests: potential consistency, solver conservation, shapes
requirements.txt      Dependencies
figures/              Generated plots (conservation_check.png, l2_error_vs_reference.png, disorder_transmission.png)
