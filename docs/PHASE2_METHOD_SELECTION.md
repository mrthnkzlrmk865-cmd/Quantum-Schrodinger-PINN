# Phase 2 preparation: reproducibility closure, method selection, validation and statistics plan

Date: 2026-11-05. Starting point: `origin/main` at `3c42f7a` (merge of PR #2,
Phase 1). Main scientific question of Phase 2:

> How does the spatial correlation length ζ of a Gaussian-correlated random
> potential affect the transmission of a quantum wave packet through a
> finite disordered section of a 2D hard-wall waveguide?

Nothing in this document changes the Hamiltonian, the physical model, the
committed baseline results or the conclusions in the README. No new
parameter sweeps were run. No independent solver was implemented (see
§2.6). Anderson localization is not claimed. No new physical discovery is
claimed.

Contents
1. Task A: PINN thread count and environment compatibility (implemented)
2. Task B: independent fixed-energy validation method (selection)
3. Task C: minimal validation benchmark with pre-defined criteria
4. Task D: paired statistical analysis of the correlation-length sweep (performed)
5. Implementation plan
6. Required tests and validation criteria (checklist)
7. Remaining scientific limitations
8. References

---

## 1. Task A: PINN cache thread count and environment compatibility

### 1.1 Problem
Phase 1 showed that PINN training depends on the PyTorch intra-op thread
count (2 threads: final rel. L2 0.48850; 4 threads: 0.48987, max |Δw| = 0.011),
but `torch_threads` was only recorded as metadata. A model trained with 2
threads could have been silently reused for a 4-thread request.

### 1.2 Changes (implemented and tested)
* `experiment_pinn.pinn_spec(cfg, torch_threads)`: the spec now contains
  `execution.torch_threads`, so the thread count is part of the cache
  fingerprint.
* `run_one(..., threads=None)` applies `--threads` **before** the cache
  spec is built and validated. The spec uses the *effective* count
  `torch.get_num_threads()`, so a run without `--threads` is also pinned to
  the count actually used. A guard checks that the count did not change
  during training.
* `provenance.py` (schema 2): an **environment fingerprint** for each
  experiment kind, built from the library versions that determine the
  floating-point results, the Python version and the CPU architecture:

  | kind | fields |
  |---|---|
  | transport, localization, validation | python, numpy, scipy, machine |
  | pinn | python, numpy, scipy, torch, machine |

  `check_cache(..., strict_env=True)` refuses a cache from a different
  environment. By default it warns, like a source-code change.
* New flags on all experiment drivers: `--strict-env`, and `--strict`
  (= `--strict-code` + `--strict-env`). **Runs whose results enter a report
  must use `--strict`.**
* Schema-1 (Phase 1) records have no environment fingerprint and no
  `machine` field. The check derives both from the versions and the
  platform string that were recorded, so Phase 1 results stay verifiable.
* Existing results were not modified. The two Phase 1 PINN reproduction
  runs have no thread count in their spec, so they can no longer match any
  request. This is intended (a test asserts it), and they stay on disk as
  records.

### 1.3 Environment review: what is and is not in the fingerprint
| Candidate | Decision | Evidence / reason |
|---|---|---|
| NumPy, SciPy versions | **included** | they implement all reference-solver arithmetic (FFT, exp, sums) |
| PyTorch version (PINN only) | **included** | kernels and reduction order differ between versions; the reference solver does not use it (a test checks that a torch change does not affect transport caches) |
| Python version | **included** | floating-point formatting and `math` functions; cheap to include |
| CPU architecture | **included** | x86_64 vs ARM changes SIMD paths and libm |
| FFT worker / CPU count (reference solver) | **excluded** | measured: wavefunctions are **bit-identical** for 1, 2 and 4 FFT workers (new regression test) |
| PyTorch thread count | **in the PINN spec** (not the environment) | measured to change training results (Phase 1) |
| BLAS/LAPACK vendor, CPU micro-architecture (AVX level), OS kernel | **recorded only** (`platform` string) | cannot be read portably and reliably; listed as a residual risk (§7) |
| GPU / device | **excluded, recorded** | all results are CPU. A future GPU run should get its own spec field |

Verified outcomes: `python -m pytest -q tests` → **78 passed, 1 xfailed**
(the xfail is the pre-existing, documented PINN 20 % accuracy target).
There are 17 new tests, including
`test_pinn_cache_trained_with_2_threads_not_reused_with_4`,
`test_pinn_thread_count_applied_before_cache_validation`,
`test_pinn_default_threads_use_effective_count`,
`test_environment_change_warns_or_fails_in_strict_mode`,
`test_reference_solver_bit_identical_across_fft_workers`, and the 6 tests in
`tests/test_paired.py`.

Observation: the code fingerprint hashes whole driver files. Adding the
`--strict-env` option to `experiment_disorder.py` therefore changed the
transport code fingerprint, so the Phase 1 results now get a "source code
changed" warning, and `--strict-code` refuses them. This is conservative by
design, but it means that **non-numerical edits to a driver also invalidate
strict reuse**. See §7.

---

## 2. Task B: independent fixed-energy validation method

### 2.1 The target problem (stationary form)
Same Hamiltonian as the time-dependent (TD) solver, ħ = m = 1:

  H = −½(∂²ₓ + ∂²ᵧ) + V(x, y),  0 < y < L_y,  ψ(x, 0) = ψ(x, L_y) = 0,

with V = V_dis·S(x) (+ optional barrier). V is negligible outside
a − 0.2 < x < b + 0.2: the window tail is ≈ 2e^{−2d/δ} ≤ 4×10⁻⁹ of W at
d = 0.2 with δ = 0.02. Outside that region the guide is clean, which
defines the asymptotic **leads**. Transverse modes are
χₙ(y) = √(2/L_y) sin(nπy/L_y) with energies Eₙ^⊥ = (nπ/L_y)²/2. At total
energy E the open channels are those with E > Eₙ^⊥, with
kₙ = √(2(E − Eₙ^⊥)). For the baseline (k_x0 = 20) there are N = 6 open
channels. The quantity to compute is the transmission from incoming mode 1:
T₁(E) = Σₘ T_{m←1}(E), with R₁(E) defined likewise and T₁ + R₁ = 1 (unitarity).

In the TD solver the CAPs replace the semi-infinite leads. A stationary
method treats the leads exactly (self-energies or mode matching), so it
tests the CAP approximation as well as the propagation.

### 2.2 Candidates compared

| Criterion | (G) Direct Green's function / sparse solve of the whole region | (RGF) Recursive Green's function, slice by slice | (MM) Mode matching / coupled-channel stationary scattering | (TMM) Transfer matrix |
|---|---|---|---|---|
| Same Hamiltonian and Dirichlet walls | yes (any discretization) | yes | yes (transverse sine basis = exact Dirichlet eigenfunctions) | yes |
| Leads / open channels | exact lead self-energies Σ_L, Σ_R; transmission via Fisher–Lee | same (Σ only at the two ends) | exact: plane waves e^{±ikₙx} per channel | exact at the ends |
| Evanescent (closed) channels | included automatically, stable | included, **stable** (only bounded quantities are propagated) | must be kept for convergence; plain shooting is unstable; needs a log-derivative/R-matrix or S-matrix propagator | **unstable** without stabilization: closed channels grow like e^{κₙL}. For n = 20, κ ≈ 59 gives e^{59} ≈ 10²⁵ over L_d = 1, which destroys all precision. Needs QR or S-matrix composition, which is effectively RGF |
| Cost per energy (L_d = 1, dx = 1/512 ≈ 614 slices, N_y = 63 modes) | sparse LU of ≈ 4×10⁴ unknowns with dense 63×63 blocks: about the same as RGF; more memory | ≈ 614 × O(N_y³) ≈ 10⁹ flop ≈ 0.1–1 s; memory ≈ 40 MB | similar to RGF when a stable propagator is used | ≈ RGF when stabilized |
| Numerical stability | excellent | excellent (standard for disordered wires) | good only with a stable propagator | poor unless stabilized |
| Uses exactly the same disorder realization | yes: evaluate `Potential(cfg, disorder_from_config(cfg, seed))` at the grid | yes, same | yes, via V_nm(x) = ∫χₙVχₘ dy | yes |
| Independence from the split-step TD solver | high: energy domain, linear solve, exact leads, no CAP, different x-discretization | high, same | high | high |
| Implementation risk | low (scipy.sparse) but memory scales with the whole system | low–moderate; textbook algorithm | moderate (propagator) | high (stabilization) |
| Literature status for disordered wires | standard (e.g. Kwant) | standard since Lee & Fisher 1981, MacKinnon 1985 | standard for waveguides (Szafer & Stone 1989; Usuki et al. 1995) | standard, with stabilization (Kramer & MacKinnon 1993) |

### 2.3 Recommendation: mode-space recursive Green's function (primary)

**Primary method: RGF in a hybrid basis.** The transverse direction uses
the sine basis χₙ(y) (n = 1…N_y, the same exact Dirichlet eigenfunctions
as the DST solver). The propagation direction x uses second-order finite
differences, and both ends get exact discrete-lead self-energies.
Transmission follows from the Fisher–Lee relation T_{m←1} = |t_{m1}|²
with t = i√v_R G_{N1} √v_L (equivalently Tr[Γ_R G Γ_L G†]).

Reasons:
1. **Exact Hamiltonian and boundary conditions.** The sine basis satisfies
   the Dirichlet walls exactly and gives exact transverse energies, so the
   channel thresholds Eₙ^⊥ are not shifted. A pure finite-difference
   transverse grid would shift them by ≈ (nπ dy)²/12, about 0.7 % for
   n = 6 at dy = 1/64. That matters for near-threshold channels.
2. **Exact treatment of the open and closed channels** through the lead
   self-energies Σₙ = −t e^{i kₙ dx}, t = 1/(2dx²), with
   cos(kₙdx) = 1 − (E − Eₙ^⊥)dx². For closed channels the root with
   |e^{ikₙdx}| < 1 is used. Leads and scattering region share the same
   discretization, so a clean guide transmits exactly (benchmark B1).
3. **Stable** for the long evanescent decay that kills a plain transfer
   matrix.
4. **Cheap**: one realization and one energy cost ≲ 1 s, so the ≈ 1300
   energies needed for energy averaging (§3, B4) take ≈ 20 min.
5. **Same realization**: V(x_j, y) is evaluated from the identical
   `Potential` object. V_nm(x_j) is obtained by sine-transforming
   V(x_j, y_l) on a y-grid finer than the TD grid, so the stationary result
   samples the continuum potential rather than copying the TD solver's
   collocation.
6. **Independent of the split-step solver** in every step that could hide
   an error: energy domain instead of time domain, linear solve instead of
   exponential splitting, finite differences in x instead of spectral, exact
   leads instead of CAPs. What is shared on purpose is only the physical
   definition (H, V, transverse basis).

**Internal cross-check (not an independent method):** a direct sparse solve
(G) of the same discretized system must agree with the RGF to 10⁻¹⁰. This
catches coding errors in the recursion, but it does not test the physics or
the discretization.

Rejected as primary:
* **Transfer matrix**: unstable with the closed channels that are needed
  for convergence at small ζ (§2.4). The stabilized forms add complexity
  and gain nothing over RGF.
* **Continuum mode matching / log-derivative propagator**: a good
  method, but it needs a carefully implemented stable propagator (Johnson
  1973), which is more new code to validate than RGF.
* **Whole-system sparse Green's function**: kept as the internal coding
  cross-check.

Optional, later: if the package `kwant` can be installed, a pure
tight-binding Kwant calculation of the same system is a third,
externally written implementation. Its transverse discretization differs
(finite differences), so it must agree only within that discretization
error. This is not a dependency.

### 2.4 Discretization requirements (estimated; to be verified by B3)
* **Longitudinal dispersion error.** The FD lead dispersion is
  E = Eₙ^⊥ + (1 − cos k dx)/dx² instead of k²/2, a relative error in k of
  about (k dx)²/24. That is ≈ 6×10⁻⁵ at dx = 1/512 for k ≈ 20 (35 at the
  top of the packet's spectrum). Plan: dx = 1/512, check against 1/1024,
  and Richardson-extrapolate.
* **Number of transverse modes N_y.** The disorder couples modes n and m
  through wave numbers |q_y| ≈ |n − m|π/L_y up to a few × 1/ζ. For
  ζ = 0.02 that means q ≈ 150, i.e. coupling to n ≈ 50. Plan: N_y = 63 (as
  in the TD grid), checked against 127; ζ = 0.02 is the stringent case.
* **Scattering-region buffer:** 0.2 on each side of the slab (window tail
  ≤ 4×10⁻⁹ W), checked against 0.4.

### 2.5 How transmission is obtained independently
For each energy E: build the lead self-energies, run the RGF recursion
through the slices, and compute t_{m1} from G_{N1}. T₁(E), R₁(E) and the
full matrices come out directly, with no time propagation, no CAP, and no
wave packet. The wave-packet comparison (§3, B4) is a separate step that
weights T₁(E) with the TD packet's known momentum distribution (§3, B4).

### 2.6 Feasibility verdict on implementing the full solver now
The method and the benchmark are well defined, and the cost is small.
I did **not** implement it in this task, because doing it safely requires:
(i) settling the analytic reference formulas against the original sources.
The smooth-step formula in Landau & Lifshitz could not be checked against
the book text in this session (§3.2).
(ii) Writing the benchmark criteria down and committing them **before** the
code exists. This document does that.
Implementing the solver is the first step of Phase 2 (§5).

---

## 3. Task C: minimal validation benchmark (criteria fixed in advance)

All thresholds below are fixed now, before any stationary-solver result
exists. If a criterion fails, the cause is investigated and reported. The
threshold is not changed and parameters are not tuned to pass.

### B1. Clean waveguide (analytic: T_{m←1} = δ_{m1}, R = 0)
* Energies: 50 values spanning the 1st to 6th channel window, avoiding
  |E − Eₙ^⊥| < 10⁻³ E (thresholds are singular in any discretization).
* Pass: |T_{1←1} − 1| ≤ 10⁻¹², Σ_{m≠1}T_{m←1} ≤ 10⁻¹², R ≤ 10⁻¹²,
  |Σ T + Σ R − 1| ≤ 10⁻¹⁰.
* TD side (already known): the clean-guide TD packet has T = 0.99989 at
  k_x0 = 20 (V6/V5). The deviation from 1 is the CAP spurious return
  (≈ 10⁻⁴), a known TD error.

### B2. Smooth rectangular barrier, uniform in y
The barrier is the existing `Potential` barrier (tanh edges, edge width e > 0).
Because V does not depend on y, the channels decouple exactly:
T_{m←1} = 0 for m ≠ 1, and T_{1←1}(E) = T_1D(E − E₁^⊥) for the 1D problem
−½φ'' + V_b(x)φ = εφ.
* Independent 1D reference: direct ODE integration (scipy `solve_ivp`,
  rtol = atol = 10⁻¹²) from a pure outgoing wave on the right back to the
  left lead. This shares no code with either solver. The integrator itself
  is validated first:
  * (a) a sharp rectangular barrier, approached by e → 0, against the
    closed form T = [1 + V₀² sin²(k′w)/(4ε(ε − V₀))]⁻¹ (sinh form for
    ε < V₀);
  * (b) the single smooth (logistic) step V₀/(1 + e^{−x/α}), whose exact
    reflection coefficient is given in Landau & Lifshitz Vol. 3, §25
    (problems). **The formula, its parameter convention and the problem
    number must be checked against the book before use.** If it cannot be
    checked, (b) is dropped and (a) plus self-convergence of the ODE is
    used.
* Heights: V_b ∈ {150, 250} (below and above E_x = 200), as in the
  existing barrier sweep, width 0.1, edge 0.01. 60 energies across the
  packet spectrum.
* Pass (RGF vs ODE): |ΔT| ≤ 2×10⁻⁴ at dx = 1/512; observed dx-convergence
  order between 1.8 and 2.2; Richardson-extrapolated |ΔT| ≤ 2×10⁻⁵;
  off-diagonal T_{m≠1←1} ≤ 10⁻¹²; unitarity ≤ 10⁻¹⁰.
* Pass (TD vs energy-averaged ODE, as in B4): |ΔT| ≤ 1×10⁻³ for both barrier
  heights.

### B3. One fixed disorder realization (internal convergence of the RGF)
Realizations: seed 0 at W = 100, L_d = 1, ζ ∈ {0.02, 0.05, 0.20}, plus seed 0
at W = 200, ζ = 0.05 (strong scattering, kℓ_tr ≈ 2.6). 5 test energies
each.
* Pass: |ΔT₁(E)| ≤ 2×10⁻⁴ for each of these refinements: dx 1/512 → 1/1024,
  N_y 63 → 127, y-quadrature ×2, buffer 0.2 → 0.4. Unitarity
  |T₁ + R₁ − 1| ≤ 10⁻¹⁰. Reciprocity |T_{m←n} − T_{n←m}| ≤ 10⁻⁸ (the
  Hamiltonian is real and symmetric). RGF vs direct sparse solve ≤ 10⁻¹⁰.

### B4. Fixed realization: stationary vs time-dependent, with the energy width accounted for
Same realizations as B3, TD solver at the production settings.

**The energy-width correction is not optional.** The TD packet is
φ(x) ∝ e^{−(x−x₀)²/2σₓ²}e^{ik₀x} in mode 1, with momentum distribution

  P(kₓ) = (σₓ/√π) exp[−(kₓ − k₀)²σₓ²],  E(kₓ) = (kₓ² + π²/L_y²)/2.

For σₓ = 0.2 the standard deviation of kₓ is 1/(√2σₓ) ≈ 3.5, so the
packet covers E ≈ 205 ± 70 (1σ). For linear scattering, with the packet
initially entirely in the left lead and moving right, the transmitted
probability at t → ∞ is exactly the energy average

  T_pkt = ∫ P(kₓ) T₁(E(kₓ)) dkₓ

(standard wave-packet scattering theory). The cross terms between
different energies vanish asymptotically. Conditions, all checked for the
baseline: initial overlap with the slab ≈ ½erfc(3.5) ≈ 4×10⁻⁷; weight of
kₓ < 0 ≈ ½erfc(4) ≈ 8×10⁻⁹; weight of kₓ < 5, which is slow and poorly
absorbed by the CAP, ≈ 10⁻⁵.

The stationary result is therefore compared as ⟨T⟩_E = ∫P T₁ dkₓ, **never
as T₁(E₀) at the central energy.**

* Energy grid: T₁(E) of a single realization fluctuates on the
  correlation-energy scale E_c ~ ħ/τ, where τ is the longest relevant dwell
  time. Measured from the committed strength sweep (∫P_in dt and the
  late-time exponential decay of ⟨P_in(t)⟩): mean dwell 0.078 (W = 100) and
  0.075 (W = 200); late-time decay time 0.17 (W = 100) and 0.19 (W = 200).
  So E_c ≈ 5–6. Plan: uniform in kₓ over k₀ ± 4.5/(√2σₓ) with
  δkₓ = 0.025 (≈ 1270 energies; δE = kₓδkₓ ≤ 0.9 ≈ E_c/6 over the whole
  window). Convergence check: halving δkₓ changes ⟨T⟩_E by ≤ 2×10⁻⁴.
  Estimated cost ≈ 20 min per realization. (An earlier draft estimate of
  E_c ~ 10 from the mean dwell time alone was too optimistic.)
* TD run: same realization, t_end extended until P_in(t_end) ≤ 10⁻⁴
  (the remaining probability is reported as an uncertainty bound,
  T ∈ [T_TD, T_TD + P_in]).
* Error budget (estimated from the Phase 1 validation): TD Δt 2×10⁻⁵, dx
  1.4×10⁻⁵, CAP ≈ 10⁻⁴, P_in ≤ 10⁻⁴. RGF dx ≈ 10⁻⁴, energy quadrature
  ≤ 2×10⁻⁴. Total ≲ 5×10⁻⁴.
* **Pass: |T_TD − ⟨T⟩_E| ≤ 1×10⁻³ and |R_TD − ⟨R⟩_E| ≤ 1×10⁻³** for every
  B4 realization. This is ≈ 2× the error budget and ≈ 10× below the
  ensemble standard error of ⟨T⟩ (≈ 10⁻²), so passing means the two
  methods agree far better than the statistics need.
* B4b (narrow packet, informative, not pass/fail): a TD run with σₓ = 0.8
  (σ_k ≈ 0.88) in a longer domain. It is compared (i) with the
  energy-averaged stationary result using the *same* criterion, and (ii)
  with the monochromatic T₁(E₀), only to **quantify** how far a naive
  monochromatic comparison would be off.

Phase 2 transport results are considered independently validated only if
B1–B4 all pass.

---

## 4. Task D: paired statistical analysis of the correlation-length sweep (performed)

### 4.1 Design
Realization i of every ζ point uses disorder seed i. For a given seed,
`make_disorder` uses the same Gaussian draws and phases at every ζ
(q = g/ζ), so the seeds are **paired**. Verified on the data: identical
seed arrays (0–47) at all 7 ζ, and specs that differ only in
`disorder_corr_length` (`paired_analysis.load_sweep` refuses anything
else). Analysis unit: the per-seed difference
dᵢ = Σⱼ wⱼ Tᵢ(ζⱼ). Bootstrap: 20 000 resamples of the seed index, shared by
all contrasts (seed 20261105). Percentile CIs: 95 % and Bonferroni 98.33 %
for the 3 primary contrasts. Paired t and Wilcoxon signed-rank p-values
with Holm adjustment over the primary family. Effect sizes:
d_z = mean(d)/sd(d), with a bootstrap CI, and d_av = mean(d)/mean(sd).

Primary contrasts:
* **C1** T(0.05) − T(0.02);
* **C2** T(0.05) − T(0.20);
* **C3** curvature K = [T(0.03) + T(0.07)]/2 − T(0.05).

Secondary: the two sides of C3; the bootstrap distribution of the grid
argmin; the vertex of a parabola in ln ζ through 0.03, 0.05, 0.07.

Transparency note: these contrasts were chosen after the *unpaired* results
had been seen (README §9.3). The paired analysis is a confirmatory
re-analysis of the same 48 realizations, not an independent replication.

Input: `results/reproduction/transport/corr_length/` (fingerprinted,
bit-identical to the committed baseline). Command: `python paired_analysis.py`.
Outputs: `results/analysis/paired_corr_length.{json,csv}`,
`figures/paired_corr_length.png`.

### 4.2 Results (n = 48 seeds, W = 100, L_d = 1, k ≈ 20.25)

| Contrast | mean ΔT | paired 95 % CI | Bonferroni 98.3 % CI | unpaired Welch 95 % CI | r (pairs) | d_z [95 % CI] | Holm p (t / Wilcoxon) |
|---|---|---|---|---|---|---|---|
| C1 T(.05)−T(.02) | −0.176 | [−0.198, −0.154] | [−0.203, −0.149] | [−0.200, −0.152] | +0.14 | −2.25 [−2.89, −1.86] | 1×10⁻¹⁹ / 2×10⁻¹⁴ |
| C2 T(.05)−T(.20) | −0.240 | [−0.284, −0.193] | [−0.294, −0.182] | [−0.290, −0.190] | +0.18 | −1.49 [−2.19, −1.07] | 3×10⁻¹³ / 6×10⁻¹¹ |
| C3 curvature K | +0.046 | [0.023, 0.068] | [0.018, 0.073] | – | – | +0.58 [0.27, 1.00] | 2×10⁻⁴ / 7×10⁻⁵ |
| C3a T(.03)−T(.05) | +0.051 | [0.029, 0.073] | [0.024, 0.078] | [0.024, 0.078] | +0.31 | +0.64 | secondary |
| C3b T(.07)−T(.05) | +0.041 | [0.006, 0.074] | **[−0.002, 0.081]** | [0.006, 0.075] | −0.01 | +0.34 | secondary |

* Grid argmin of ⟨T⟩: ζ = 0.05 in 98.9 % of bootstrap replicates, ζ = 0.07 in 1.1 %.
* Parabola vertex (model-dependent, 3 points): kζ* = 0.95, 95 % CI
  [0.88, 1.12]. For comparison, the minimum of the Born transport mean
  free path is at kζ = 0.889, and the minimum of the energy-averaged
  Born + diffusion estimate of T is at kζ = 0.883. These lie at the lower
  edge of the CI.

### 4.3 Interpretation (limited to what the analysis supports)
* With seed-wise pairing and correction for three primary comparisons,
  ⟨T⟩ at ζ = 0.05 is lower than at ζ = 0.02 and at ζ = 0.20, with large
  effect sizes, and ζ = 0.05 lies below the mean of its two grid
  neighbours (C3).
* **The minimum is resolved only coarsely.** The difference to ζ = 0.07 is
  small (d_z = 0.34), and its multiplicity-corrected interval includes
  zero. The data locate the minimum between ζ ≈ 0.03 and 0.07, i.e.
  kζ ≈ 0.6–1.4. A finer ζ grid is needed for a sharper statement.
* Pairing helps little here: the per-seed correlation of T across ζ is only
  r ≈ −0.01 to +0.31, because the same draw rescaled in ζ gives a quite
  different scattering landscape. The paired and unpaired intervals are
  therefore similar. The paired analysis is still the correct one.
* These results describe this geometry (N = 6 channels, L_d = 1,
  W/E₀ ≈ 0.49, energy-averaged packet) and the solver. Their numerical
  correctness is what the B1–B4 benchmark is designed to confirm
  independently.

---

## 5. Implementation plan (Phase 2, in order)

1. **Reference formulas.** Check the smooth-step formula against Landau &
   Lifshitz, and write `analytic_scattering.py` (rectangular barrier, smooth
   step) plus `ode1d.py` (independent 1D integrator). Tests: integrator vs
   formulas to 10⁻⁸.
2. **`stationary.py`: mode-space RGF.** Lead modes and self-energies, the
   recursion, Fisher–Lee transmission and reflection matrices. A direct
   sparse-solve variant serves as the internal cross-check. It uses the
   `Potential`/`disorder` objects unchanged.
3. **`validate_stationary.py`**: runs B1–B4 with `provenance` metadata and
   `--strict`, writes `results/stationary_validation/` and a pass/fail table
   against the criteria in §3. Estimated cost: B1–B3 a few minutes. B4 is
   4 realizations × (≈ 20 min RGF energy grid, ≈ 40 min with the halving
   check, + ≈ 10 min TD with extended t_end) ≈ 3–3.5 h. B4b ≈ 40 min.
4. **Report** in `docs/PHASE2_VALIDATION.md`: every criterion with its
   measured value, including any failure.
5. Only after B1–B4 pass: a **finer ζ grid** around the minimum
   (e.g. ζ ∈ {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.07}, 48–96 paired
   seeds, `--strict`). Its pre-registered contrasts are the location of
   the minimum (vertex CI) and C3-type curvature. Optionally, a
   stationary-solver ζ sweep at fixed energy to separate the
   energy-averaging effect.
6. Phase-2 science items from README §9.11 (longer slabs, channel-number
   dependence) only after step 5.

## 6. Required tests and validation criteria (checklist)

Implemented now (passing):
- [x] PINN cache rejects a different torch thread count; the thread count is applied before validation; the effective count is used by default
- [x] Environment fingerprint; strict refusal; warning by default; schema-1 derivation; torch irrelevant for transport
- [x] Reference solver bit-identical across FFT worker counts
- [x] Paired analysis: CI coverage, pairing benefit on correlated data, curvature/vertex, Holm, refusal of unpaired or inconsistent inputs, loads the reproduced data

Required for the stationary solver (to be written with it):
- [ ] `test_lead_modes`: open/closed classification; |e^{ikdx}| < 1 for closed modes; FD dispersion relation
- [ ] `test_clean_guide_unit_transmission` (B1 thresholds)
- [ ] `test_unitarity_and_reciprocity` (random weak potential, 10⁻¹⁰ / 10⁻⁸)
- [ ] `test_rgf_equals_direct_sparse` (10⁻¹⁰)
- [ ] `test_y_uniform_barrier_decouples_channels` (off-diagonal ≤ 10⁻¹²)
- [ ] `test_ode1d_vs_rectangular_formula` and (if the source is verified) `test_ode1d_vs_smooth_step` (10⁻⁸)
- [ ] `test_rgf_barrier_vs_ode1d` (2×10⁻⁴, order ≈ 2)
- [ ] `test_same_realization` (V on the RGF slices equals `Potential.numpy` at the same points, bitwise)
- [ ] `test_energy_average_weights` (P(kₓ) normalized; reproduces T_pkt = 1 for the clean guide)
- [ ] B3/B4 as slow, marked integration tests (run by `validate_stationary.py`, not by the default suite)

## 7. Remaining scientific limitations

1. **No independent validation yet.** All transport results still come from
   one solver family (split-step + CAP), validated by convergence,
   analytic free/eigenmode tests and conservation. B1–B4 are designed to
   close this gap.
2. **Energy averaging.** Every TD result is an average over a broad band
   (E ≈ 205 ± 70). Statements about "transmission at kζ ≈ 1" refer to the
   central k. The stationary solver will make fixed-energy statements
   possible.
3. **Coarse ζ grid at the minimum.** The minimum is located only to
   kζ ≈ 0.6–1.4. The ζ = 0.07 side is not separated after multiplicity
   correction.
4. **Re-analysis, not replication.** The paired contrasts were chosen
   after the unpaired results were seen. A finer, pre-registered sweep with
   new seeds (e.g. 1000–1095) would be an independent confirmation.
5. **Single geometry.** N = 6 channels, L_d = 1, one disorder strength for
   the ζ sweep. Born-approximation comparisons are order-of-magnitude
   (unbounded 2D, no channel structure).
6. **Finite observation time and CAP.** P_in(t_end) up to 10⁻³; CAP return
   ≈ 10⁻⁴ at k ≈ 20 and 1.3 % at k_x0 = 12.
7. **No localization claim.** Nothing here addresses localization. The
   ζ-dependence concerns ⟨T⟩ in a diffusive/quasi-ballistic regime
   (kℓ_tr ≈ 10–50 for W = 100).
8. **Reproducibility residuals.** BLAS vendor and CPU micro-architecture are
   recorded but not fingerprinted. Bit-identity has been shown on one
   machine only. The code fingerprint is per file, so non-numerical driver
   edits trigger it (conservative, sometimes inconvenient). The PINN study
   (17 of 18 runs) has not been re-run with a recorded, fixed thread
   count.

## 8. References (new in this document; others in docs/LITERATURE.md)

* P. A. Lee, D. S. Fisher, "Anderson localization in two dimensions," *Phys. Rev. Lett.* **47**, 882 (1981). doi:10.1103/PhysRevLett.47.882
* D. S. Fisher, P. A. Lee, "Relation between conductivity and transmission matrix," *Phys. Rev. B* **23**, 6851 (1981). doi:10.1103/PhysRevB.23.6851
* A. MacKinnon, "The calculation of transport properties and density of states of disordered solids," *Z. Phys. B* **59**, 385 (1985). doi:10.1007/BF01328846
* S. Datta, *Electronic Transport in Mesoscopic Systems*, Cambridge University Press (1995).
* A. Szafer, A. D. Stone, "Theory of quantum conduction through a constriction," *Phys. Rev. Lett.* **62**, 300 (1989). doi:10.1103/PhysRevLett.62.300
* Usuki et al., "Numerical analysis of ballistic-electron transport in magnetic fields by using a quantum point contact and a quantum wire," *Phys. Rev. B* **52**, 8244 (1995). doi:10.1103/PhysRevB.52.8244 (author order differs between indexes; check on the APS page before citing)
* B. R. Johnson, "The multichannel log-derivative method for scattering calculations," *J. Comput. Phys.* **13**, 445 (1973). doi:10.1016/0021-9991(73)90049-1
* C. W. Groth, M. Wimmer, A. R. Akhmerov, X. Waintal, "Kwant: a software package for quantum transport," *New J. Phys.* **16**, 063065 (2014). doi:10.1088/1367-2630/16/6/063065
* L. D. Landau, E. M. Lifshitz, *Quantum Mechanics: Non-Relativistic Theory* (Course of Theoretical Physics Vol. 3), §25 and its problems (smooth-step reflection coefficient; to be checked against the text before use).
* J. R. Taylor, *Scattering Theory: The Quantum Theory of Nonrelativistic Collisions*, Wiley (1972): wave-packet description of scattering.
* B. Kramer, A. MacKinnon, *Rep. Prog. Phys.* **56**, 1469 (1993) (transfer-matrix methods; already in docs/LITERATURE.md).
