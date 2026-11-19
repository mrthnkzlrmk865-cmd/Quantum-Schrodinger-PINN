# B4 Implementation Readiness Audit

**Scope.** A read-only audit of whether the planned benchmark B4 can be implemented and run. B4 compares stationary
(RGF) transmission and reflection with time-dependent (TDSE) wave-packet transmission and reflection for the same
disorder realization.

**Protocol source.** [PHASE2_METHOD_SELECTION.md](PHASE2_METHOD_SELECTION.md) §3 (B4), together with §2.3–§2.5, §5,
§6 and §9.

**What this audit did not do:**

- It did not modify any source code, protocol, threshold, result or report.
- It did not run any B4 computation, RGF energy sweep or TDSE propagation.

**What was run:**

- Inspection commands.
- An existing unit-test subset: 50 passed and 1 xfailed (the documented B1-F xfail), in 50 s.
- Configuration and hash comparisons. These generate disorder *parameters* only; no field is propagated.
- A read-only analysis of committed TD files.

No new benchmark was run. Every runtime figure below is either a measured value quoted from a committed artifact or
document, with its source given, or an arithmetic extrapolation from such a value, labelled **ESTIMATE**.

**Status vocabulary in this report:**

| Category | Meaning |
|---|---|
| **[V]** | existing validated result |
| **[L]** | existing but limited numerical evidence |
| **[A]** | unvalidated assumption |
| **[M]** | missing implementation |
| **[P]** | proposed future work |

---

## 0. Summary

- **B4 is not implemented.**
  - No B4 driver, configuration, protocol file, test, provenance entry or result directory exists.
  - The components it needs exist separately:
    - the single-energy RGF S-matrix;
    - the split-step TDSE solver;
    - the shared continuum disorder field;
    - the packet-weight, window, trapezoid, halving and interval helpers (inside the B2 driver).
  - B4 itself has never been executed.
- **The same-realization requirement can be met with existing code.** The RGF (via the B3 path) and the committed TD
  runs use bitwise-identical disorder parameters for all four B4 configurations (§7).
- **Three issues must be settled before any B4 scientific run:**
  1. **RGF resolution policy.** The B4 error budget assumes an RGF dx error of about 10⁻⁴. B3/B3X/B3Y measured
     extrapolated dx = 1/512 errors up to 8.0×10⁻⁴ for B4 configuration 4 (W = 200, ζ = 0.05) and 2.4×10⁻⁴ for
     configuration 2. They did so only at 5 energies per configuration, covering E ≈ 88–371 of the B4 window
     E ≈ 13–650.
  2. **TD stopping rule.** No committed seed-0 TD run of the four B4 realizations reaches P_in ≤ 10⁻⁴: at t = 1 the
     values are 3.2×10⁻⁴ to 2.2×10⁻³. Late-time P_in(t) is **non-monotone** in three of the four (§6.5). "Extend
     until P_in ≤ 10⁻⁴" is therefore not yet an operational rule.
  3. **CAP spurious return.** It was measured only for mode-1 packets. Disorder scatters into higher transverse
     channels with small longitudinal wave number, for which the measured mode-1 return reaches 10⁻²–3×10⁻² (§6.7).
- **Overall classification: READY FOR LIMITED PILOT.**
  - The pilot is a set of labelled TD-only and RGF-only *diagnostics* (§9), **not** B4.
  - The full B4 campaign is **NOT READY**.

---

## 1. Repository-wide B4 implementation discovery

**Searches performed:**

- `grep -rIl "\bB4\b"` over the whole repository (excluding `.git`), followed by inspection of every code hit;
- function and class listings of all stationary, TD, disorder, configuration and driver modules;
- searches for `energy_average`, `packet_weights`, `kx_window`, `interval_distance` and `P_in`;
- listing of `results/` for any B4 directory;
- `experiments/`, which contains no "B4" string.

**Findings:**

- **No file implements B4.** There is no `validate_stationary.py` (§5 step 3 of the protocol), no
  `validate_rgf_b4.py`, no `docs/rgf_b4_protocol.json`, no `tests/test_rgf_b4.py`, and no
  `results/stationary_validation/` or `results/rgf_b4_validation/`.
- **"B4" in code occurs only in comments** that cite the B4 packet convention: `validate_rgf_b2.py`,
  `validate_rgf_b3.py`, `validate_ode1d.py`, `tests/test_ode1d.py` and `tests/test_rgf_b2.py`.
- **The protocol's own test checklist (§6) leaves the B4 tests unchecked:** `test_energy_average_weights`,
  `test_b4_interval_criterion`, and "B3/B4 as slow, marked integration tests".
- **Classification: NOT IMPLEMENTED; never executed.** Implementation status was determined from the code and
  `results/`, not from documentation.

**Components relevant to B4:**

| File | Function / class | What it does | Reusable for B4? | Missing or incompatible |
|---|---|---|---|---|
| `stationary.py` | `transverse_basis`, `potential_blocks`, `build_device`, `rgf_endpoint_blocks`, `scattering_matrix` | Flux-normalized S over the open channels at **one** energy, for a device built from a potential callable; dx, N_y, N_q and the device extent are arguments; reports min rcond | **Yes** (core RGF) | no energy loop; no packet weighting; dx has no default (must be chosen by the caller); `min_closed = 1` check only |
| `smatrix_spec.py` | `lead_modes`, `fisher_lee`, `assemble_S`, `split_S`, `caroli_transmission`, `unitarity_residual`, `reciprocity_residual` | Specified conventions (lead roots, Σ, Γ, Fisher–Lee, S ordering, checks) | **Yes** | none |
| `validate_rgf_b3.py` | `case_config`, `case_potential`, `geometry`, `buffer_length`, `make_rgf_device`, `T1`, `unitarity`, `reciprocity`, `Recorder`, `direct_reference_samples` | The four B4 realizations (cases A–D are exactly the B4 configurations); nested-grid geometry with a dx-independent buffer (103/512); bitwise same-realization check; T₁ = Σₘ\|(t_RL)ₘ₁\|² | **Yes** (import or copy) | no R₁ helper (R₁ = Σₘ\|(r_L)ₘ₁\|² from `split_S` is trivial); driver-local functions, not a library |
| `validate_rgf_b2.py` | `packet_weights`, `kx_window`, `energy_average`, `interval_distance`, `adjust_energy`, `run_td` | The B4 weighting P(kₓ) = (σₓ/√π)exp[−(kₓ−k₀)²σₓ²]; uniform kₓ window k₀ ± 4.5σ_k with δkₓ = 0.025 (1274 points); trapezoid; halving by a δkₓ/2 grid whose even points equal the coarse grid; interval distance; threshold rule; and a complete **TD vs energy-averaged 1D-ODE** comparison for the y-uniform barrier | **Partly** | written for the single-channel ODE, not the multichannel RGF; no handling of threshold-moved points inside the quadrature; module-private to a driver |
| `sparse_reference.py`, `sparse_wavefunction.py`, `sparse_exact_lattice.py` | `sparse_scattering_matrix`, `wavefunction_scattering_matrix`, `exact_lattice_S` | Independent references for the **same lattice problem** (B3/B3X/B3Y) | Optional spot checks | 3.8–15 s per energy (B3Y II.1); not suitable for a 1274-point sweep |
| `reference_solver.py` | `SplitStepDST` (`run` with an `observables` callback), `build_grid`, `cap_profile` | Strang split-step in a 2D DST-I basis with a quadratic CAP; exact absorbed-probability bookkeeping per side | **Yes** (TD core) | fixed `t_end`; no adaptive stop |
| `transport.py` | `run_transport` | T(t) = absorbed_right + P(x > b); R(t) = absorbed_left + P(x < a); **P_in = 1 − T − R** | **Yes**, or call `SplitStepDST.run` directly | P_in is a residual, not a direct sum; no norm or gap-region output; fixed t_end |
| `initial_condition.py` | `initial_condition_numpy` (`packet = "waveguide_mode"`) | Gaussian in x × mode-1 sine, analytically normalized | **Yes** | none |
| `disorder.py` / `potential.py` | `make_disorder`, `disorder_from_config`, `DisorderRealization.numpy`, `Potential` | Continuum Gaussian-correlated field from (seed, M = 512); evaluated at any points; tanh-windowed slab | **Yes** | none |
| `config.py` | `waveguide_config` | Production TD settings: L_x = 3.9, nx = 249 (dx = 3.9/250 = 0.0156), ny = 63, dt = 2×10⁻⁴, T = 1, CAP width 0.6, strength 300, power 2, x₀ = 1.2, σₓ = 0.2, k_x0 = 20 | **Yes** | T = 1 is too short for P_in ≤ 10⁻⁴ (§6.5) |
| `validate_solver.py` | V3, V4, V5, V6 | TD dt, dx, CAP and conservation validation | Evidence only | V4 covers W = 200 at ζ = 0.02 and 0.05 at t = 1; V5 covers mode-1 packets only |
| `analyze_transport.py` | `est_energy_averaged` | Born-rate estimate averaged over 401 points | **No** | a theory estimate, not an RGF average |
| `provenance.py` | per-benchmark entries (e.g. `rgf_b3`) | Code fingerprints and metadata | Pattern reusable | no B4 entry |

---

## 2. Reconstruction of the B4 protocol (from PHASE2_METHOD_SELECTION.md §3)

**What the protocol states:**

- **Realizations:** "Same realizations as B3, TD solver at the production settings."
  - B3: seed 0, L_d = 1, W = 100 with ζ ∈ {0.02, 0.05, 0.20}, and W = 200 with ζ = 0.05.
- **Packet:** φ ∝ e^{−(x−x₀)²/2σₓ²}e^{ik₀x} in mode 1, with momentum distribution
  P(kₓ) = (σₓ/√π)exp[−(kₓ−k₀)²σₓ²] and E(kₓ) = (kₓ² + π²/L_y²)/2.
- **Stationary quantity:** ⟨T⟩_E = ∫P(kₓ)T₁(E(kₓ))dkₓ, with T₁(E) = Σₘ T^{L→R}_{m←1}(E) summed over the open
  channels; ⟨R⟩_E likewise. The protocol states explicitly that it is "never … T₁(E₀) at the central energy alone".
- **Energy grid:**
  - uniform in **kₓ** over k₀ ± 4.5/(√2σₓ) = [4.090, 35.910], δkₓ = 0.025 ("≈ 1270 energies");
  - E ∈ [13.30, 649.70], with 1 to 11 open channels and 10 thresholds inside the window;
  - grid points within 10⁻³E of a threshold are moved off it.
- **Convergence check:** halving δkₓ must change ⟨T⟩_E by ≤ 2×10⁻⁴ (*kept*).
- **TD side:**
  - the same realization, extended until P_in(t_end) ≤ 10⁻⁴; this is **a pass condition**, not only a stopping rule;
  - the TD result is an interval: 𝒯 = [T_TD, T_TD + P_in], ℛ = [R_TD, R_TD + P_in].
- **Pass, per realization:**
  - P_in(t_end) ≤ 10⁻⁴; **and**
  - dist(⟨T⟩_E, 𝒯) ≤ 10⁻³; **and**
  - dist(⟨R⟩_E, ℛ) ≤ 10⁻³.
- **Error budget (§3):**

  | Source | Budgeted error |
  |---|---|
  | TD Δt | 2×10⁻⁵ |
  | TD dx | 1.4×10⁻⁵ |
  | CAP | ≈ 10⁻⁴ |
  | P_in | ≤ 10⁻⁴ |
  | RGF dx | ≈ 10⁻⁴ |
  | energy quadrature | ≤ 2×10⁻⁴ |
  | **total** | ≲ 5×10⁻⁴ |

- **B4b (informative, not pass/fail):** σₓ = 0.8 in a longer domain, compared with the average and with T₁(E₀).
- **Validation rule:** Phase 2 transport is "independently validated only if B1–B4 all pass".

**Discrepancies between the task description and the repository** (reported; the protocol is unchanged):

| # | Task description | Repository protocol | Note |
|---|---|---|---|
| D1 | "≈ 1270 energy points" | "≈ 1270" | the existing helper gives **1274** points with δkₓ = 0.024996 (`kx_window`; also used in B2) |
| D2 | "repeat with half the original **energy** step" | halve the **kₓ** step (the grid is uniform in kₓ, so δE = kₓδkₓ is not uniform) | criterion: change of ⟨T⟩_E ≤ 2×10⁻⁴ |
| D3 | "per-realization T and R differences ≤ 10⁻³" | **interval** distance to [T_TD, T_TD + P_in] and [R_TD, R_TD + P_in] ≤ 10⁻³, **plus** P_in ≤ 10⁻⁴ as a pass condition | a plain \|ΔT\| would be wrong by up to P_in |
| D4 | — | the RGF resolution for B4 is **not stated** in §3 B4 | §2.4 plans "dx = 1/512, check against 1/1024, and Richardson-extrapolate"; the budget assumes an RGF dx error ≈ 10⁻⁴. B3/B3X recommend establishing dx convergence per B4 realization and using dx ≤ 1/1024 or a verified Richardson extrapolation (recommendation only; protocol not changed) |
| D5 | — | "TD solver at the production settings" | production T = 1 does not satisfy the P_in condition (§6.5); "extended t_end" is part of the protocol |
| D6 | — | B4b exists | informative only |

**How the packet's energy distribution is represented.**

- The protocol and the B2 helper use the exact analytic |φ̃(kₓ)|² of the initial Gaussian. The TD kinetic operator is
  spectral, so E(kₓ) is exact on the TD side.
- The RGF is evaluated at the same total E, not at a lattice wave number.
- A single RGF calculation at E₀ is **not** used. B2 measured why: the monochromatic value differs from the average
  by 0.024–0.026 for the barrier, which is 30–70× the tolerance.

---

## 3. Implications of B1–B3Y (statuses unchanged)

Three different questions are kept apart:

- **(a) solver agreement** on the same lattice problem;
- **(b) accuracy of the discretized lattice problem**, i.e. convergence in dx, N_y, N_q and buffer;
- **(c) convergence to the continuum Schrödinger problem.**

| Benchmark | Status | What it established | What it did not establish |
|---|---|---|---|
| B1 | 5 of 6 configurations pass; **B1-F fails** the provisional phase criterion (N = 3483) | lattice algebra, leads and Fisher–Lee in a clean guide; E ∈ [4.93, 241.8] only | potential matrix elements, mode mixing, energies above 241.8 |
| B2 | **PASS** | RGF vs ODE for a y-uniform barrier at 60 energies over the whole B4 window (ε ∈ [8.4, 644.8]); **the only executed TD vs energy-averaged comparison** (barrier): distances 3.2×10⁻⁴ and 7.0×10⁻⁴ | channel coupling, closed channels as intermediates, y-dependent V; TD vs RGF |
| B3 | **does not pass** (dx criterion) | (a) RGF = sparse to 4.9×10⁻¹²; N_y, N_q, buffer converged; **dx 1/512 → 1/1024 change 5.97×10⁻⁴** (case D) | energies outside the 5 per case; dx beyond 1/1024 (case D exploratory only) |
| B3X | **does not pass** (C11 dx criterion) | (a) with a second independent reference; dx ladder 1/256 … 1/2048 (RGF only) at 30 points; observed order 1.96–2.14 | the references at dx ≠ 1/512; other energies; the continuum |
| B3Y | **PASS** (its own predeclared criteria) | (a)+(b): four methods, including an extended-precision exact-lattice solution, reproduce the same complex S (≤ 9.1×10⁻¹¹) and the same T₁ dx dependence (spread ≤ 3.3×10⁻¹¹; order 1.996–2.006) at 16 points (dx 1/512, 1/1024) and at 5 points at 1/2048 | **not (c)**; not the B4 energy range; one realization; a physical transport result neither claimed nor implied |

**Coverage of the B4 configurations by existing convergence evidence.**

- All existing disorder evidence is at the **5 energies kₓ = k₀ + jσ_k, j = −2…2**, i.e.
  E ∈ {88.51, 140.47, 204.93, 281.90, 371.36}.
- B4 integrates over kₓ ∈ [4.09, 35.91], i.e. E ∈ [13.3, 649.7].
- Under P(kₓ), the region kₓ < 12.93 or kₓ > 27.07 (outside ±2σ_k) carries 2(1 − Φ(2)) = **4.6 % of the packet
  weight**, and it has **no** disorder convergence evidence at all.
  - This is not small. For example, a T₁ error of 10⁻² confined to those tails would contribute ≈ 5×10⁻⁴ to ⟨T⟩_E.
  - The energies there (1–3 open channels below E ≈ 88; 9–11 above E ≈ 371, up to 650) are also outside the
    tested range of N_y convergence and outside B1's energy range.

**Evidence per B4 configuration** (each uses 5 energies at E ≈ 88–371 only):

| B4 configuration (B3 case) | dx 1/512→1/1024 change (B3/B3X C11) | extrapolated dx = 1/512 error (B3X diagnostic) | extrapolated dx = 1/1024 error | dx = 1/2048 data | multi-method agreement at refined dx |
|---|---|---|---|---|---|
| 1: W = 100, ζ = 0.02 (A) | 7.0×10⁻⁵ (pass) | ≤ 9.4×10⁻⁵ | ≤ 2.3×10⁻⁵ | RGF only (B3X ladder) | none (M2/M3 at 1/512 only) |
| 2: W = 100, ζ = 0.05 (B) | 1.81×10⁻⁴ (pass, near the limit) | ≤ **2.4×10⁻⁴** | ≤ 6.1×10⁻⁵ | RGF only | B3Y control point E = 371.36, at 1/512 and 1/1024 |
| 3: W = 100, ζ = 0.20 (C) | 3.4×10⁻⁵ (pass) | ≤ 4.6×10⁻⁵ | ≤ 1.1×10⁻⁵ | RGF only | none |
| 4: W = 200, ζ = 0.05 (D) | **5.97×10⁻⁴ (FAIL)** | ≤ **8.0×10⁻⁴** | ≤ **2.0×10⁻⁴** | RGF (all 5 energies, B3X) and **M1–M4 at 3 energies** (B3Y Tier 2/3) | B3Y: 1/512, 1/1024, 1/2048 at E = 140.47, 204.93, 371.36 |

**Reading the table:**

- **[L]** The "extrapolated" errors are Richardson ESTIMATES from three or four lattice grids. B3Y states that they
  are not continuum statements.
- **[L]** The sign of ΔT₁ varies between energies (case D: +2.8×10⁻⁴, −5.97×10⁻⁴, …). The size of the dx error of the
  **energy average** ⟨T⟩_E is therefore **unknown**: it may partly cancel or not. No existing computation measures it.
- **[A]** Assuming that the 5-point evidence represents the whole window (in particular E < 88 and E > 371, with up to
  11 open channels and 10 threshold cusps) is an **unvalidated assumption**.
- **[V]** B3Y establishes only that the dx dependence is a property of the finite-difference lattice problem,
  reproduced by independent solvers. It does **not** validate any B4 energy point beyond its 16 or 5 points, and it
  is not a physical transport result.

---

## 4. RGF numerical-resolution readiness

**Resolution parameters:**

- **Default dx.** `stationary.build_device` has **no default**; dx is a required argument. The B3 base is 1/512
  (`validate_rgf_b3.BASE`), and B3X and B3Y used 1/512, 1/1024 and 1/2048.
- **dx is configurable** and validated: x_end − x_start must be an integer multiple of dx.
- **Slice count.** The geometry `validate_rgf_b3.geometry` uses a buffer fixed on the 1/512 lattice (103/512 = 0.2012)
  at every dx, so N = (b − a)/dx + 2·(0.2012/dx) + 1:

  | dx | 1/512 | 1/1024 | 1/2048 |
  |---|---|---|---|
  | N (L_d = 1) | **719** | **1437** | **2873** |

  These counts are confirmed in B3Y II.1. The grids are nested: the 1/512 slices are a subset of the finer ones.
- **Same realization at every resolution.**
  - The potential is a continuum callable evaluated at the slice positions x_j and the quadrature points
    y_l = l/(N_q + 1).
  - Nothing is interpolated or resampled, so changing dx changes **only** the sampling, not the realization.
  - Verified: B3 bitwise identity, 28/28; B3Y P1–P4 and cross-resolution identity on coincident slices.
- **Multi-resolution energy sweep.** It can be run with existing functions (a loop over E per device), but no such
  sweep code exists **[M]**.
- **Transverse resolution.** N_y = 63, N_q = 127 or 254.
  - It is converged to ≤ 10⁻⁷ (N_y) and ≤ 1.9×10⁻⁶ (N_q) **at E ≤ 371** **[L]**.
  - At E up to 650 there are 11 open and 52 closed channels at N_y = 63; convergence there has not been tested **[A]**.

**Clean-lead guard.** The lead-tail guard of 10⁻⁶ was met in B3Y at 1/2048. However, B3X notes that W200_z0.02 is
within 4 % of the guard at 1/2048. That is not a B4 configuration, but it shows that a finer dx may need a larger
buffer.

**Feasibility of comparing 1/512, 1/1024 and 1/2048.**

- Measured per-energy RGF cost: 0.4 / 0.8 / 1.6 s at 1/512 / 1/1024 / 1/2048 (B3Y II.1, M1, one BLAS thread). B1 §3
  gives 0.5 s at N = 718.
- **ESTIMATE** per realization, using 2547 energies (the δkₓ/2 grid, which contains the 1274-point grid):

  | dx | Time per realization |
  |---|---|
  | 1/512 | ≈ 17 min |
  | 1/1024 | ≈ 34 min |
  | 1/2048 | ≈ 68 min |

  - All three resolutions for all four realizations: ≈ 8 h serial.
  - The energies are independent, so the sweep parallelizes across cores. `OPENBLAS_NUM_THREADS=1` is required for
    this speed (`stationary.py` docstring: 16–22 s per energy with multithreaded BLAS).
  - The cost does not depend on E (it is O(N·N_y³)).
- **Memory, ESTIMATE from array shapes:** at 1/2048 the real potential blocks take 2873 × 63² × 8 B ≈ 91 MB, and the
  H blocks take another ≈ 91 MB. That is well within the 16 GB of the validation machine.
- **No benchmark was run for this audit**, because the existing measurements suffice.

**Resolution policy.** The existing evidence does **not** support declaring any single dx sufficient for all four
configurations across the B4 window:

- **Configurations 1 and 3:** dx = 1/512 is probably within the budget (≤ 9.4×10⁻⁵ estimated), but only at 5 energies.
- **Configuration 2:** 1/512 exceeds the 10⁻⁴ budget at the measured points (≈ 2.4×10⁻⁴); 1/1024 is within it
  (≈ 6×10⁻⁵).
- **Configuration 4:** 1/1024 is at about 2×10⁻⁴ (2× the budget); 1/2048 is estimated at ≈ 5×10⁻⁵ (B3X/B3Y second-order
  scaling, not measured as an error).

**Additional validation required before choosing a policy:**

- **(i)** A predeclared rule for the B4 RGF resolution: for example, compute ⟨T⟩_E and ⟨R⟩_E at two or three dx for
  every realization and report the change of the **energy average**; Richardson extrapolation only inside the
  predeclared order window.
- **(ii)** dx-ladder evidence at energies **outside** E ≈ 88–371.
- **(iii)** N_y spot checks at the top of the window.

---

## 5. Energy-integration audit

**What the existing helpers in `validate_rgf_b2.py` do:**

- **Weights.** `packet_weights` returns the analytic P(kₓ), whose integral over all kₓ is 1. It is **not
  renormalized** on the window.
  - `energy_average` returns ∫P on the grid: 1 − 6.8×10⁻⁶ for ±4.5σ_k (test
    `test_energy_average_weights_and_interval`, which passed in this audit).
  - The truncation therefore biases ⟨T⟩ and ⟨R⟩ by at most 6.8×10⁻⁶ (B2 reports this as `window_truncation_bound`).
- **Grid.** `kx_window` gives a uniform kₓ grid with 1274 points; the halving grid has 2547 points and its even
  points equal the coarse grid (asserted to 10⁻¹²).
- **Integrand coverage.** T₁ must be evaluated at **every** grid point; nothing is interpolated in B2 **[V]**.
- **Limits.** k₀ ± 4.5σ_k covers the packet spectrum. The weight of kₓ < 0 is ≈ 8×10⁻⁹, and of kₓ < 5 ≈ 10⁻⁵
  (protocol).
- **Halving check.** Implemented for the ODE average (B2: change 2.3×10⁻¹⁰) **[V for the barrier]**. Not implemented
  for an RGF average **[M]**.
- **Reflection.** The same function applies to R₁. B2 does so for the ODE **[V]**. For the RGF,
  R₁ = Σₘ|(r_L)ₘ₁|² is not yet coded **[M]**.
- **Units.** ħ = m = 1; E(kₓ) = (kₓ² + π²/L_y²)/2 with L_y = 1 is the same expression in the protocol, B2 and B3.
  - The RGF evaluates S at the total E. The TD packet's transverse energy is exactly π²/2 (an exact sine mode on the
    DST grid), and both solvers place the thresholds Eₙ^⊥ = (nπ)²/2 exactly, because both use exact sine transverse
    modes.
  - **Consistent.**

**Issues found:**

1. **Threshold-moved points [A].** `adjust_energy` moves a point within 10⁻³E of a threshold to Eₙ(1 ± 2×10⁻³), but
   the trapezoid weight stays at the original kₓ.
   - On the δkₓ = 0.025 grid, δE = kₓδkₓ ≤ 0.9, while 10⁻³E ≈ 0.02–0.65, so some points will be moved. Their number
     has not been computed in this audit.
   - The protocol assigns the resulting error, and the √-cusp error of T₁ at the ten thresholds, to the halving check.
     Whether halving detects a cusp error reliably is an **unvalidated assumption**: the trapezoid error near a √ cusp
     converges more slowly than O(δkₓ²).
2. **No interpolation is used**, so interpolation adds no error. The cost is a full RGF solve per point.
3. **The definitions differ, so a small integration error is not physical equivalence:**
   - The RGF gives probabilities for ideal semi-infinite leads at fixed E.
   - The TD gives absorbed probability plus the probability in the gap regions (x > 2.9 or x < 1.9) at t_end, in a
     finite box with CAPs and hard end walls.
   - The two agree only asymptotically, and only if the CAP returns nothing (§6.7) and the gap regions hold no
     probability that will later re-enter (§6.5).
   - The interval criterion covers P_in only, **not** the gap or CAP contents.

---

## 6. TDSE readiness

**6.1 Same realization.** `run_transport` builds `Potential(cfg, disorder_from_config(cfg, seed))` and evaluates it on
the TD grid. This is the same callable the RGF samples (§7) **[V]**.

**6.2 Discretization consistency.** The TD and RGF solvers discretize the same continuum potential in different ways.
This is intended, since B4 compares two discretizations of one continuum problem:

| | TD | RGF |
|---|---|---|
| x | spectral (DST-I), dx = 0.0156 (nx = 249) | finite differences, 1/512 or finer |
| y | collocation, ny = 63 (dy = 1/64) | DST-I quadrature, N_q = 127 or 254 (dy = 1/128 or 1/255) |
| leads | CAP | exact |

**TD spatial convergence:**

- **[L]** V4: dx 1/64 → 1/128 changes T by 3.7×10⁻⁶ (W = 200, ζ = 0.02) and 1.3×10⁻⁵ (W = 200, ζ = 0.05), for seeds
  0 and 1 at t = 1, with dt = 10⁻⁴.
- Not measured for W = 100 at ζ = 0.20 or 0.02, nor at the extended t_end.
- **[L]** B2 measured a TD dx error of 2.6×10⁻⁴–5.7×10⁻⁴ for the 0.01-wide barrier edges, which the B4 budget
  (1.4×10⁻⁵) underestimates for sharp features. The disorder slab edges have width 0.02.

**6.3 Initial packet.** `packet = "waveguide_mode"`: a Gaussian with σₓ = 0.2 at x₀ = 1.2, times √2 sin(πy), times
e^{i20x}, analytically normalized. Its |φ̃(kₓ)|² is exactly the protocol's P(kₓ), so the **TD and RGF energy
distributions match by construction [V]**.

Two small TD-side biases are **not** in the error budget:

- **(a)** The initial |ψ|² inside the left CAP (x < 0.6) is ½erfc(3) ≈ 1.1×10⁻⁵ (analytic). It is absorbed and booked
  as R.
- **(b)** The overlap with the slab is ≈ 4×10⁻⁷ (protocol).

Both are ≪ 10⁻³ **[L, analytic]**.

**6.4 Measurement.** T = absorbed_right + P(x > 2.9); R = absorbed_left + P(x < 1.9); **P_in = 1 − T − R**, a residual.

- P_in is defined on the nominal slab [1.9, 2.9]; the tanh window tails extend slightly beyond it.
- Because P_in is a residual, any bookkeeping or normalization error is silently folded into P_in. The solver's
  bookkeeping is exact to 1.6×10⁻¹³ (V6), but `run_transport` does not record the norm.

**6.5 P_in as a stopping or pass condition (committed data, seed 0, production settings, t_end = 1):**

| Configuration | T_TD | R_TD | P_in(1) | ESTIMATE t for P_in = 10⁻⁴ (exponential fit, t ≥ 0.7) | P_in(t) monotone for t ≥ 0.7? |
|---|---|---|---|---|---|
| 1: W100, ζ0.02 | 0.728707 | 0.270294 | 9.99×10⁻⁴ | ≈ 1.7 | yes |
| 2: W100, ζ0.05 | 0.511014 | 0.488668 | 3.17×10⁻⁴ | ≈ 1.3 (fit over the last 0.1 is nearly flat; unreliable) | **no** (0.900: 3.52e-4 → 0.975: 4.12e-4 → 1.000: 3.17e-4) |
| 3: W100, ζ0.20 | 0.886358 | 0.111464 | 2.18×10⁻³ | ≈ 2.1–2.9 | **no** (re-rise at 0.850) |
| 4: W200, ζ0.05 | 0.328163 | 0.671025 | 8.12×10⁻⁴ | ≈ 1.7–2.5 | **no** (re-rises at 0.800–0.825 and 0.925) |

Sources: `results/transport/corr_length/zeta{0.02,0.05,0.2}.npz` and `results/transport/strength/W200.npz`, all
48 seeds, seed 0 column. Frames are spaced Δt = 0.025.

**Findings:**

- **[V]** None of the four satisfies P_in ≤ 10⁻⁴ at the production t_end = 1. **Every B4 TD run must be extended.**
- **[L]** P_in(t) is not monotone for one realization.
  - The cause cannot be established from the repository. Candidates are slow (near-threshold or high-transverse-mode)
    components crossing the slab, re-entry of CAP-returned waves, and trapping in resonances.
  - Consequences:
    1. A "first time P_in ≤ 10⁻⁴" rule can stop at a point after which probability re-enters.
    2. The exponential extrapolations above are crude ESTIMATES, not predictions.
    3. P_in ≤ 10⁻⁴ does not bound the probability still in the gap regions or inside the CAP layers.
- **[M]** No adaptive stop exists. A rule must be predeclared, for example a fixed t_max per configuration, plus a
  requirement that P_in ≤ 10⁻⁴ holds over a trailing window, and reporting of the gap and CAP contents.

**6.6 Conservation.** `SplitStepDST` books the absorbed probability exactly (norm + absorbed = 1 to 1.6×10⁻¹³; V6)
**[V]**. A B4 driver should record the norm and compute P_in **directly** (Σ|ψ|² over the slab) as well as by residual.
This is possible today through the `observables` callback, without modifying any file.

**6.7 Boundary bias (CAP).**

- **[V/L]** V5 (clean guide, mode-1 packets) measured the spurious CAP return as a function of the packet's
  longitudinal wave number:

  | k_x0 | 10 | 12 | 16 | 20 | 25 | 30 |
  |---|---|---|---|---|---|---|
  | spurious return | 3.1×10⁻² | 1.3×10⁻² | 1.3×10⁻³ | 1.0×10⁻⁴ | 1.6×10⁻⁴ | 5.2×10⁻⁴ |

- W_cap depends on x only, so channels decouple inside the clean CAP region, and the return depends on each channel's
  **longitudinal** wave number kₓ,ₙ = √(2E − (nπ)²).
- Disorder converts mode-1 flux into higher channels. At E = 205, channel n = 5 has kₓ ≈ 12.8 and n = 6 has
  kₓ ≈ 7.4. In the V5 data such kₓ values have a 10⁻²-level return.
- **[A]** The B4 budget's CAP ≈ 10⁻⁴ is therefore an **unvalidated assumption** for disordered scattering. Depending
  on the (unmeasured) fraction of flux in low-kₓ channels, the CAP return could reach the 10⁻³ tolerance.
- A B4 failure could then reflect a TD boundary artifact rather than a stationary–TD disagreement.
- The returned probability is not covered by the P_in interval.

**6.8 Truncation.** Hard walls at x = 0 and x = L_x lie behind the CAPs. Anything not absorbed is reflected back. A
longer t_end gives more time for absorption, but also for CAP-returned waves to re-enter the slab (§6.5).

**6.9 Cost.**

- **Measured:** a single-seed production run took 5.8–6.2 s for t_end = 1.5 (`results/rgf_b2_validation/td_energy_average.json`,
  live reruns), and a dx = 1/128 run took 45–46 s for the same t_end.
- **ESTIMATE:** one B4 TD run at t_end ≈ 2–3 takes ≈ 10–15 s at production resolution and ≈ 1–1.5 min at dx = 1/128.
- The protocol's "≈ 10 min TD" is conservative. TD cost is not a constraint.

---

## 7. Same-realization requirement

**Checks, all from existing code and data:**

| Item | Result |
|---|---|
| Seed handling | `make_disorder` uses `SeedSequence([seed, 20241105])`; the parameters (q, φ) depend only on (seed, M); ζ and W scale them deterministically |
| Disorder generation | identical function for both solvers (`disorder_from_config`) |
| Spatial correlation | continuum spectral representation (M = 512 cosines); no grid in the definition |
| Interpolation / resampling | none; both solvers call `DisorderRealization.numpy` at their own points |
| Normalization | none per realization (by design); identical for both |
| Geometry | slab (1.9, 2.9), edge 0.02, L_y = 1 in both. The RGF device is (1.9 − 0.2012, 2.9 + 0.2012) with clean ideal leads; the TD box is [0, 3.9] with CAPs on [0, 0.6] and [3.3, 3.9] |

**Configuration and hash check (this audit; no propagation).**

- For each B4 configuration, `waveguide_config(W, ζ)` (the B3 path) and the configuration stored in the committed TD
  file give the same disorder parameters, compared by sha256 over (qx, qy, φ):

  | Configuration | Hash prefix |
  |---|---|
  | A | `f022790469f9d24a` |
  | B | `f47346b43bf5a752` |
  | C | `2122ece8feba83ac` |
  | D | `f47346b43bf5a752` |

- No solver-relevant configuration field differs.
- The prefixes match the hashes recorded in RGF_B3_VALIDATION.md §5.

**Conclusion.** **The same physical realization is available to both methods [V]. Not a blocking issue.**

- A B4 driver must still **assert** it at run time: equality of the configurations and disorder hashes, bitwise RGF
  slice samples (as in B3), and a recorded hash of the TD-grid potential **[M]**.
- B and D are the same draw at two amplitudes, and A, B and C the same draw rescaled. The four B4 realizations are not
  independent samples (B3 §5). This is acceptable for a per-realization method comparison.

---

## 8. B1-F and other known limitations

**B1-F (status unchanged: FAIL of the provisional 10⁻¹⁰ phase criterion).**

- **What failed.** A phase error of t of 1.54×10⁻¹⁰ (genuine RGF rounding, confirmed at high precision) on a
  **clean**, 3483-slice device (L_d = 3, dx = 1/1024), worst near thresholds. |t| was exact to ≤ 3.2×10⁻¹³ there.
- **Mechanism (B1 §6).** Rounding of ε·2t_x per slice, with a bound ∝ N·ε/sin(kdx).
- **Geometry and length for B4.** L_d = 1, with N = 719 / 1437 / 2873 slices at 1/512 / 1/1024 / 1/2048.
  - At 1/512 and 1/1024 the devices are shorter than B1-F, and equal to or longer than B1 configurations that passed.
  - **At dx = 1/2048, N = 2873 is close to B1-F's 3483**, and the bound scales as N/(k·dx). An ESTIMATE from the bound
    gives ≈ 1.6× B1-F's bound, so phase errors of a few ×10⁻¹⁰ are possible near thresholds.
- **Direct effect on B4.**
  - B4 uses probabilities |t|² and |r|², not phases.
  - In a disordered device, rounding-level phase errors inside the device perturb interference only at the same
    relative level, so the expected effect on T₁ is ~10⁻¹⁰, which is 10⁷ below the B4 tolerance.
  - **[L]** This is supported by B3Y: four independent methods, including the extended-precision lattice solution,
    agree on the full complex S to ≤ 9.1×10⁻¹¹ at 1/2048 for 5 disordered points.
  - **[A]** Near-threshold energies at 1/2048 (which B4 contains and B3Y did not test) remain untested.
- **Verdict.** Not a blocking issue for B4's probability criteria. Not dismissed: a B4 driver should record min rcond
  (already returned) and unitarity U per energy, and flag near-threshold points.

**Other limitations carried into B4:**

- B1 tested only E ≤ 241.8, whereas B4 reaches 650. B2 covers the whole window, but only for decoupled channels.
- N_y has been tested only at E ≤ 371.
- Bit-identity across BLAS thread counts has been shown only for the clean guide (B1 §8).
- All measurements were made on one machine.

---

## 9. Final readiness assessment

| Component | Status | Evidence | Blocking issue? | Required action |
|---|---|---|---|---|
| B4 driver, configuration, protocol file, tests, outputs | NOT IMPLEMENTED | none found anywhere (§1) | **Yes** | write them (P0.2–P0.4) |
| Single-energy RGF S-matrix (`stationary.scattering_matrix`) | READY | B1 (5/6), B2 PASS, B3/B3X solver agreement ≤ 6×10⁻¹², B3Y ≤ 9.1×10⁻¹¹; unit tests pass | No | reuse unchanged |
| RGF energy sweep over 1274 / 2547 points | NOT IMPLEMENTED | — | **Yes** | part of the driver; cost ESTIMATE §4 |
| RGF dx accuracy over the B4 window (lattice → continuum) | NOT VALIDATED | B3/B3X/B3Y: 5 energies per configuration (E 88–371); configuration 4 dx = 1/512 error ≈ 8×10⁻⁴, configuration 2 ≈ 2.4×10⁻⁴ (estimates) | **Yes** (configurations 2 and 4) | predeclared resolution policy (P0.1); dx evidence outside E 88–371 (P1.1) |
| Transverse truncation N_y = 63 / N_q | PARTIALLY READY | converged ≤ 1.9×10⁻⁶ at E ≤ 371 | No (P1) | spot checks at E 371–650 (P1.2) |
| Packet weights, window, trapezoid, halving grid | PARTIALLY READY | `validate_rgf_b2` helpers; unit test passes; used in B2 | No | reuse or copy into a B4 module; multichannel use untested |
| Threshold rule and quadrature near the 10 thresholds | PARTIALLY READY | `adjust_energy` tested; cusp/moved-point quadrature effect untested | No | test and report the moved points; predeclare their weighting (P0.1, P0.3) |
| Interval criterion | PARTIALLY READY | `interval_distance` tested (B2) | No | B4 test with P_in logic (P0.3) |
| R₁ from the RGF | NOT IMPLEMENTED | trivial from `split_S` | No | in the driver, with a test |
| TD solver at production settings | READY | V1–V7; B2 bitwise rerun | No | extend t_end |
| TD stopping rule (P_in ≤ 10⁻⁴) | NOT IMPLEMENTED | committed seed-0 P_in(1) = 3.2×10⁻⁴–2.2×10⁻³; non-monotone tails (§6.5) | **Yes** | predeclared rule (P0.1); direct P_in, gap and CAP contents (P0.2) |
| TD resolution for the four realizations | PARTIALLY READY | V3 (W200, t = 1), V4 (W200, ζ 0.02/0.05, t = 1) | No (P1) | informative dx = 1/128 and dt = 10⁻⁴ reruns at the B4 t_end (P1.3) |
| CAP spurious return under mode conversion | NOT VALIDATED | V5 mode-1 only: 10⁻²-level return at kₓ ≈ 10–12 | **Yes** (for interpreting a FAIL) | TD-only CAP diagnostic (P1.4; may be part of the pilot) |
| T/R definition consistency (RGF vs TD) | PARTIALLY READY | the same E(kₓ), P(kₓ) and thresholds; TD includes gap regions and CAP bookkeeping | No | record gap and CAP contents (P0.2) |
| Probability conservation monitoring | PARTIALLY READY | solver bookkeeping exact (V6); `run_transport` does not record the norm; P_in is a residual | No | observables callback (P0.2) |
| Same realization (RGF ↔ TD) | READY | hashes identical (§7); B3 bitwise sampling test | No | runtime assertion in the driver (P0.2) |
| B1-F impact on B4 | NOT VALIDATED | probabilities unaffected at measured points (B3Y); near-threshold points at 1/2048 untested | No | record rcond and U; flag near-threshold points |
| Provenance and metadata for B4 | NOT IMPLEMENTED | `provenance.py` has per-benchmark entries | No (needed for a reportable run) | additive entry (P0.2) |
| B4b (σₓ = 0.8, informative) | NOT IMPLEMENTED | — | No | P2.1 |
| Compute resources | READY | per-energy and TD costs measured; ESTIMATE ≈ 8 h serial for three dx × four realizations | No | parallelize across energies |

**Overall B4 readiness: READY FOR LIMITED PILOT.**

- **The limited pilot is not B4 and must not be reported as B4.** It consists of predeclared, separately labelled
  diagnostics that use existing components, written to a separate directory (e.g. `results/b4_pilot/`):
  - **(i)** TD-only long-time runs of the four seed-0 realizations, recording P_in(t) (direct and residual), the gap
    and CAP contents and the norm. Purpose: design the stopping rule.
  - **(ii)** a TD-only CAP sensitivity check (§6.7);
  - **(iii)** RGF-only dx ladders (1/512, 1/1024, 1/2048) at a **small, predeclared** set of energies outside
    E ≈ 88–371, to decide the resolution policy.
- The pilot must not compute ⟨T⟩_E vs TD comparisons for the four realizations. Such numbers could be mistaken for B4
  results and would contaminate the predeclaration.
- The full four-configuration campaign is NOT READY until every P0 item is resolved.

---

## 10. Proposed next steps (none implemented)

### P0 — before any B4 scientific run

**P0.1 Protocol amendment, fixed before any B4 number** (new `docs/rgf_b4_protocol.json` and
`docs/RGF_B4_VALIDATION.md` Part I).

- **What it must fix in advance:**
  - **(a) RGF resolution policy** per configuration, for example compute ⟨T⟩_E and ⟨R⟩_E at 1/1024 and 1/2048,
    report the change in the average, and use a Richardson estimate only if the per-energy order lies in [1.8, 2.2];
  - **(b) TD stopping rule**: t_max per configuration, the frame spacing, and the requirement that P_in ≤ 10⁻⁴ hold
    over a trailing window; which frame defines T_TD, R_TD and P_in;
  - **(c) treatment of threshold-moved points** in the quadrature;
  - **(d) reporting** of the gap and CAP contents;
  - **(e) separation** of pilot and official outputs.
- **Why.** The §3 B4 text leaves the RGF dx and the stopping procedure unspecified (D4, D5). The existing evidence
  shows that both matter at the 10⁻³ level.
- **Verified by:** a test that the driver's constants equal the JSON values, as B2 and B3 do
  (`test_predeclared_*_preserved`).
- **Changes the B4 protocol?** It specifies items the protocol leaves open. No threshold changes: 10⁻³, 10⁻⁴ and
  2×10⁻⁴ are all kept.
- **Amendment required?** **Yes**, as a documented, dated amendment committed before results, following the B3X/B3Y
  precedent.

**P0.2 B4 driver** (new `validate_rgf_b4.py`; additive entry in `provenance.py`).

- **Functions:**
  - realization setup via `validate_rgf_b3.case_config` and `case_potential`;
  - nested geometry via `geometry`;
  - the 1274 / 2547 kₓ grid via `kx_window`, with `adjust_energy` per point;
  - `scattering_matrix` per energy, giving T₁, R₁, U, Q, min rcond and n_open;
  - per-point caching;
  - TD via `SplitStepDST.run` with an observables callback for the direct slab, gap and CAP probabilities and the
    norm;
  - the interval criterion;
  - run-time assertions of configuration and disorder-hash equality and bitwise RGF slice samples;
  - the TD-grid potential hash;
  - `OPENBLAS_NUM_THREADS=1`;
  - provenance metadata;
  - a fail-closed output directory.
- **Why.** Nothing exists.
- **Verified by:** P0.3.
- **Protocol change?** No (it implements §3 together with P0.1).
- **Note.** No existing source file needs to change: the TD observables hook already exists. The `provenance.py`
  entry is additive, the same pattern as `rgf_b3`.

**P0.3 Tests** (new `tests/test_rgf_b4.py`):

- **(a)** weights: ∫P = 1 − erfc(4.5/√2); mean k₀; halving grid nested;
- **(b)** clean guide: the RGF T₁ ≡ 1 gives ⟨T⟩_E = ∫P and ⟨R⟩_E ≈ 0 — the unchecked `test_energy_average_weights`
  item of §6;
- **(c)** `test_b4_interval_criterion` with P_in, including both endpoints;
- **(d)** the threshold-moved points are listed, and each stays on its own side;
- **(e)** the same realization on the RGF slices and the TD grid, and a mutation (a different seed or ζ) is detected;
- **(f)** the stopping-rule logic on synthetic P_in(t) curves, including the non-monotone case;
- **(g)** direct P_in equals the residual P_in to rounding on a short run;
- **(h)** the predeclared constants are preserved;
- **(i)** the slow integration runs are marked and excluded from the default suite.

**Protocol change?** No.

### P1 — before the full four-configuration campaign

- **P1.1 RGF dx evidence outside E ≈ 88–371** for the four realizations (driver option or pilot; RGF only).
  - **Why:** the B4 window carries 4.6 % of its weight there, with no evidence.
  - **Verified by:** observed orders and consecutive changes reported per energy.
  - **Protocol change:** no; a diagnostic under P0.1.
- **P1.2 N_y 63 → 127 (fixed N_q) spot checks at E ≈ 450–650.**
  - **Why:** untested above E = 371, where there are up to 11 open channels.
  - **Verified by:** |ΔT₁| reported against the existing 2×10⁻⁴ value.
  - **Protocol change:** no.
- **P1.3 Informative TD refinements at the B4 t_end** (dx = 1/128 and dt = 10⁻⁴, the four realizations).
  - **Why:** V4 covers only W = 200 at t = 1, and B2 showed larger TD dx errors for sharp features.
  - **Protocol change:** no (informative, like the B2 refinements).
- **P1.4 TD CAP sensitivity:** the four realizations rerun with a wider or stronger CAP (and a correspondingly longer
  box), plus mode-n or low-kₓ clean-guide return measurements.
  - **Why:** the CAP budget for mode-converted flux is unvalidated (§6.7).
  - **Verified by:** the change of T_TD and R_TD reported; it must be ≪ 10⁻³ for B4 to be interpretable.
  - **Protocol change?** Not if it is informative. If production CAP settings must change for B4, then **yes**, an
    amendment, because the protocol says "production settings".
- **P1.5 Resource plan:** parallel energy evaluation, per-point cache, a single-thread BLAS check of bit-identity on a
  disordered device.
  - **Why:** bit-identity across thread counts has been shown only for the clean guide.
  - **Protocol change:** no.

### P2 — useful, not blocking

- **P2.1 B4b** (σₓ = 0.8, longer domain; informative per the protocol).
- **P2.2 Move the energy-average helpers into a shared module** (new file; the B2 driver stays untouched).
- **P2.3 Sparse-reference (M2/M3) spot checks** at a few B4 energies outside the B3 set, at the chosen dx.
- **P2.4 Optional Kwant cross-check** (protocol §2.3).
- **P2.5 Update the §6 checklist and README status** after B4 exists. This is a documentation change; it is not
  done here.

None of the P2 items changes the protocol.

---

## 11. Scientific-integrity classification

**1. Existing validated results [V]:**

- B2 PASS, including its TD vs energy-averaged ODE comparison for the y-uniform barrier.
- B3Y PASS: solver agreement, and the lattice dx dependence reproduced by four methods.
- The RGF–sparse agreement in B3 and B3X.
- The TD bookkeeping (V6).
- The same disorder parameters for the RGF and TD paths (§7).
- B1 5/6 and B1-F FAIL, B3 FAIL (dx) and B3X FAIL (C11): all statuses unchanged.

**2. Existing but limited evidence [L]:**

- dx and N_y convergence at 5 energies per configuration (E 88–371).
- The Richardson error estimates.
- TD dx and dt convergence for W = 200 at t = 1 only.
- CAP return for mode-1 packets only.
- The committed seed-0 P_in(t) up to t = 1, non-monotone.
- B3Y full-S agreement at 1/2048 for 5 points.

**3. Unvalidated assumptions [A]:**

- RGF accuracy at E < 88 and E > 371.
- The B4 budget's RGF dx error of 10⁻⁴ for configurations 2 and 4.
- CAP return ≈ 10⁻⁴ under mode conversion.
- That P_in ≤ 10⁻⁴ is reached and stays reached.
- That halving detects the threshold-cusp quadrature error.
- N_y = 63 adequacy at E up to 650.
- Negligible B1-F-type rounding at near-threshold points at 1/2048.

**4. Missing implementation [M]:**

- The B4 driver, protocol file and tests.
- The RGF energy sweep and R₁.
- The TD adaptive or predeclared stopping.
- Direct P_in, gap and CAP recording.
- The B4 provenance entry.
- B4b.

**5. Proposed future work [P]:** §10, P0–P2.

**What this audit does not establish.** It establishes no B4 result, no TD–RGF agreement for any disordered
realization, no continuum convergence, and no physical transport statement. It does not imply anything about Anderson
localization. Phase 2 transport results remain **not independently validated**, because B4 has not been performed.
