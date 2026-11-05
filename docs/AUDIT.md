# Audit of the original implementation (commit `42cea06`)

This is the audit done **before** any changes. Every number below came from
running the original code at commit `42cea06` (check out that commit to
reproduce). Environment: CPU only, PyTorch 2.x, NumPy/SciPy, float64
reference solver.

Severity: **[C]** critical (invalidates a claimed result), **[M]** major,
**[m]** minor.

---

## 1. Mathematical inconsistencies

| # | Finding | Severity |
|---|---------|----------|
| 1.1 | README states a **1D** TDSE, a 1D initial condition `(1/πσ²)^{1/4}`, and **periodic** boundary conditions plus a periodic BC loss. The code solves a **2D** problem with **Dirichlet** walls and enforces IC/BC by a hard constraint, not by losses. The README loss does not include the normalization term the code actually uses. | [M] |
| 1.2 | The PDE residuals in `pinn_core.py` (`f_u = u_t + ½∇²v − Vv`, `f_v = v_t − ½∇²u + Vu`) are **correct** for `i∂ψ/∂t = −½∇²ψ + Vψ`. | OK |
| 1.3 | IC normalization `A = 1/(σ√π)` is correct for a 2D Gaussian (`∫A²e^{-r²/σ²} = A²πσ² = 1`). | OK |
| 1.4 | Hard constraint `u = u₀ + B(x,y)·t·NN` gives `u = u₀(wall) ≠ 0` on the walls (≈3.5×10⁻⁵ at x = 0). Negligible here, but the BC is not satisfied exactly as the docstring implies. | [m] |

## 2. Incorrect physical assumptions

| # | Finding | Severity |
|---|---------|----------|
| 2.1 | **"Tunneling" scenario is mostly over-barrier scattering.** E₀ = k₀²/2 = 800 > V_b = 600. With σ = 0.06 the momentum spread is σ_k = 1/(√2σ) ≈ 11.8, so only the components with k < √1200 ≈ 34.6 (≈ 32 % of the packet) are below the barrier. | [M] |
| 2.2 | **"Transmission" is measured in a closed box after ~40 wall traversals.** v = k₀ = 40, L = 1, T = 1: the packet hits the barrier at t ≈ 0.0055 and then bounces ~40 times between the walls. `P(x > 0.55)` at t = T is a snapshot of a strongly fluctuating quantity (the original right-hand panel oscillates between 0.15 and 0.6). It is not a transmission coefficient. | [C] |
| 2.3 | The headline "disorder increases transmission (0.424 vs 0.277)" compares **one time snapshot of the clean case** against the disorder average. The disordered value ≈ 0.42 is close to the area fraction to the right of the barrier (0.45), i.e. consistent with disorder simply **randomizing (ergodizing) the probability over the closed box**. There is no evidence for disorder-enhanced transport. | [C] |
| 2.4 | ENAQT / photosynthetic motivation does not match the model: there is no environment, dephasing, Lindblad or HEOM dynamics. The model is unitary with static disorder. | [M] |

## 3. Numerical instability / accuracy

| # | Finding | Severity |
|---|---------|----------|
| 3.1 | **The reference solver is not converged at the settings used.** At t = T, N = 128: 4000 vs 16000 time steps differ by **51 %** (relative L2). N = 127 vs N = 255 (16000 steps, compared on the shared grid points) differ by **177 %**. Even 16000 vs 64000 steps differ by 3.1 %. Long-time dynamics in a chaotic-ish closed box with a sharp barrier amplify small discretization errors. The "ground truth" was therefore not ground truth. | [C] |
| 3.2 | The sharp (discontinuous) barrier sits at a different place relative to the grid at every resolution, so spatial convergence studies cannot isolate discretization error. | [M] |
| 3.3 | The PINN has to represent ≈127 temporal phase cycles (E₀T/2π) and 40 wall traversals, while the Fourier features only cover 1–10 cycles per unit input. This is a representation mismatch, not a tuning problem. Loss scales were not normalized (residual ≈ E₀·|ψ| ~ 10⁴, so the residual MSE is ~10⁷ at initialization, as the code comments note). | [M] |

## 4. PINN vs reference solver

| # | Finding | Severity |
|---|---------|----------|
| 4.1 | **The saved model `schrodinger_pinn.pt` has failed.** Re-evaluating it against the reference solver gives a relative L2 error of **1.25–1.45 for every t > 0**. This is worse than predicting ψ = 0, which gives exactly 1.0. Total probability drops to 0.70–0.89. The README presents these plots as "verifying unitarity" and as a "numerical benchmark validation". | [C] |
| 4.2 | The 20 MB `wavefunction_animation.html` (served by `index.html` on GitHub Pages) shows this failed model's output as if it were the wavefunction. | [M] |

## 5. Normalization

| # | Finding | Severity |
|---|---------|----------|
| 5.1 | The PINN IC is normalized analytically and the reference IC is normalized on the grid. They agree to ~10⁻⁶. OK, but the two IC implementations are duplicated (`initial_condition.py`, `reference_solver.initial_wavepacket`) and could drift apart. | [m] |
| 5.2 | **The disorder field is re-standardized for each realization** (zero sample mean, unit sample std). This changes the ensemble statistics and makes the definition depend on the grid. | [m] |

## 6. Boundary conditions

| # | Finding | Severity |
|---|---------|----------|
| 6.1 | The DST-I solver enforces ψ = 0 exactly on Dirichlet walls. The PINN enforces it approximately (1.4). Consistent in practice. | OK |
| 6.2 | The disorder field is built with an FFT, so it is **periodic** on the grid even though the physical domain has hard walls. This is harmless but undocumented. | [m] |

## 7. Disorder implementation

| # | Finding | Severity |
|---|---------|----------|
| 7.1 | **The correlation length is wrong and depends on the resolution.** The Fourier kernel uses `σ_k = 1/(2πℓn)` in integer-wavenumber units. The resulting real-space correlation length is ≈ 1/(4π²ℓn) ≈ 1 grid spacing, not ℓ = 0.03. Measured lag where the correlation falls below e^{-1/2}: 0.0206 (N = 96), 0.0155 (N = 128), 0.0078 (N = 256), i.e. ∝ dx. Every disorder result is effectively grid-scale white noise whose properties change when the grid is refined. | [C] |
| 7.2 | The random field is drawn on the solver grid with `default_rng(seed)`, so the **same seed gives a different physical potential at every resolution**. Convergence testing is impossible. | [C] |
| 7.3 | `potential_torch` raises `NotImplementedError` for disorder, so the PINN cannot be compared with the reference solver on any disordered landscape. | [M] |

## 8. Reproducibility

| # | Finding | Severity |
|---|---------|----------|
| 8.1 | `experiment_disorder.py` mutates the global `cfg` (`SCENARIO`, `eval_n_grid`, `eval_n_steps`, `disorder_seed`). Results depend on the call order. | [M] |
| 8.2 | `config.py` seeds torch at import time (hidden global side effect). Configs and seeds are not saved with outputs. Raw data are not saved (only PNGs). | [M] |
| 8.3 | The README lists `run.sh`, `tests/test_basic.py`, `figures/` and `REPORT.md`; **none of these exist**. The tracked PNGs are in the repository root. `wavefunction_animation.html` is tracked although `.gitignore` lists it. | [M] |
| 8.4 | `visualize.py`'s docstring names a nonexistent `visualize_wavefunction.py`. | [m] |

## 9. Statistics

| # | Finding | Severity |
|---|---------|----------|
| 9.1 | 8 realizations, one disorder strength, one correlation length; mean and std only; no standard error, confidence intervals, convergence-in-N_real check or distribution. The clean reference is a single non-converged snapshot (2.2, 3.1). | [M] |

## 10. Unsupported claims (README / comments)

* "Original contribution … ENAQT-inspired": no ENAQT mechanism is modeled (2.4).
* "Verifies that the neural network enforces quantum unitarity": contradicted by the figure (≈30 % loss of norm).
* "High-precision reference solver": not converged (3.1).
* "Visualizes spatial density scattering, … interference patterns" for `disorder_transmission.png`: the figure shows two scalar time series.
* "probability conserved to machine precision … so any transmission trend we see here is physical": conservation of the norm says nothing about discretization accuracy or about the meaning of the observable.

## 11. Missing tests

The test suite referenced in the README does not exist. Nothing checks the
potential consistency between torch and numpy, the PDE residual, the IC, the
BC, solver accuracy or convergence, disorder statistics, or PINN/reference
agreement.

## What is sound and will be preserved

* The DST-I Strang split-step scheme is mathematically correct for Dirichlet
  walls (verified below against exact eigenmode evolution and the analytic
  free Gaussian packet). It is kept and generalized.
* The real/imaginary residual formulation, the hard IC constraint, the Fourier
  feature embedding, Adam→L-BFGS training and the normalization loss are kept
  as options.
* The debugging history in the README (four PINN failure modes) is kept as an
  honest engineering record.

---

## Prioritized implementation plan

1. **Shared, resolution-independent, differentiable disorder** (`disorder.py`).
   Gaussian-correlated random field as a random Fourier series, evaluated
   identically in NumPy and PyTorch. One seed gives one physical field at any
   resolution, for both solvers.
2. **One potential and one initial-condition definition** used by both solvers.
   Pointwise-agreement tests.
3. **Generalize the DST solver**: rectangular domain (waveguide), optional complex
   absorbing potential (CAP) with exact absorbed-flux bookkeeping, so
   transmission and reflection become well-defined asymptotic scattering
   probabilities instead of closed-box snapshots.
4. **Validate the solver**: exact eigenmode phase, analytic free Gaussian,
   Δt² convergence, spatial convergence, CAP spurious reflection, norm
   conservation, convergence of the observables under disorder.
5. **Physics experiments with statistics**: sweeps over the disorder strength,
   the correlation length, the momentum, the slab length and the barrier
   geometry. Many realizations, bootstrap CIs, a convergence-in-N_real check,
   raw CSV/NPZ with configs and seeds. Compare with the Born-approximation mean
   free path and quasi-1D localization scales. A localization feasibility study
   (participation ratio and mean-square displacement in closed boxes of two
   sizes, transmission vs length).
6. **PINN on a feasible problem**: the PINN uses the shared potential, with
   input normalization and residual scaling. Compare clean vs disordered,
   strength, correlation length, Fourier features on/off, and uniform vs
   residual-adaptive sampling, against the validated reference solver.
7. **Documentation**: rewrite the README, verify every reference, remove ENAQT,
   state the limitations.
