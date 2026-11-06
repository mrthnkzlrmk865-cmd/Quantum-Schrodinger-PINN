# Stationary RGF solver: implementation and B1 (clean waveguide) validation

Phase 2, implementation stage 2. Date: 2026-11-06. Base commit: `bce6587`.
Files: `stationary.py` (solver), `tests/test_stationary.py`,
`validate_rgf_b1.py`, `results/rgf_b1_validation/`.

**Scope.** Only benchmark B1 (clean, matched waveguide) of
[PHASE2_METHOD_SELECTION.md](PHASE2_METHOD_SELECTION.md) §3 is performed.
B1 is an internal-consistency test: the clean guide has an exact analytic
answer *on the lattice*. **It is not an independent physical validation.**
It cannot detect errors that only appear with a non-zero potential. B2
(barrier vs the 1D ODE), B3 (disorder realization) and B4 (vs the
time-dependent solver) are **not implemented**. No disorder simulation was
run, and no baseline result was changed.

**Outcome in one line:** all B1 criteria pass on the 50-energy spec grid in
5 of 6 tested configurations. The **provisional** phase criterion fails on
the longest device (dx = 1/1024, L_d = 3, N = 3483): 1.77×10⁻¹⁰ against
10⁻¹⁰ (§6). A reproducible high-precision audit of all six configurations
(§6a) confirms that this is a genuine RGF phase error (1.54×10⁻¹⁰ against
the exact lattice phase), not a reference artifact. The threshold was not
relaxed. **B1 therefore does not fully pass.**

Reproduce:
```bash
python -m pytest -q tests/test_stationary.py       # ~40 s
python validate_rgf_b1.py                          # ~11 min on 4 CPU cores; --quick: ~1 min
python audit_rgf_phase.py                          # high-precision phase audit, ~10 min (§6a)
python -m pytest -q tests/test_phase_audit.py
```
Both set `OPENBLAS_NUM_THREADS=1` if it is unset (§7).

---

## 1. Mathematical formulation (implemented contract: spec §2.7)
Units ħ = m = 1. Total energy E. Transverse Dirichlet modes
χₙ(y) = √(2/L_y) sin(nπy/L_y), with Eₙ^⊥ = (nπ/L_y)²/2, n = 1…N_y. All N_y
modes, open and closed, are kept.

* **Slices:** x_j = x_start + (j − 1)dx, j = 1…N. The device spans the
  transport region (slab ± ≥ 0.2 buffer, rounded up to whole dx). The leads
  are clean lattices continuing beyond both ends. `build_device` raises an
  error if |V| > 10⁻⁶ on the first lead slice on either side.
* **Hamiltonian blocks:** H_jj = 2t_x I + diag(E^⊥) + V_j and
  H_{j,j±1} = −t_x I, with t_x = 1/(2dx²).
* **Potential blocks:** V_j = Xᵀ diag(V(x_j, y_l)) X · L_y/(N_q + 1), with
  X_{ln} = χₙ(y_l), y_l = l·L_y/(N_q + 1), l = 1…N_q, and N_q = 2N_y + 1
  (≥ 2N_y as required). V(x, y) is the existing `Potential(cfg).numpy`
  callable; the test checks that the samples are bitwise identical. The
  product is symmetrized, (B + Bᵀ)/2, which changes it only at rounding
  level, so the blocks are exactly symmetric.
* **Retarded Green function:** G = M⁻¹ with M = E − H_D − Σ_L ⊕ Σ_R (Σ_L on
  block (1,1), Σ_R on block (N,N)). Its off-diagonal blocks are
  M_{j,j±1} = +t_x I.

## 2. Lead self-energy: derivation and branch choice
For channel n of a clean lead, the recurrence
−t_x c_{j−1} + (2t_x + Eₙ^⊥ − E)c_j − t_x c_{j+1} = 0 with c_j = λ^j gives

  λ + 1/λ = 2cₙ,  cₙ = 1 − (E − Eₙ^⊥)dx².

The surface Green function of the semi-infinite lead satisfies
g = [E − (2t_x + Eₙ^⊥) − t_x²g]⁻¹. Its two roots are −λ/t_x and −1/(t_xλ).
The retarded root is the one analytic for E + i0⁺, i.e. |λ| < 1:
gₙ = −λₙ/t_x. With coupling τ = −t_x I, Σₙ = τ†gτ = **−t_x λₙ**, the same
for both leads.

* **Open channel** (|cₙ| < 1): λₙ = e^{+ikₙdx} with kₙdx = arccos cₙ ∈ (0, π).
  This is the outgoing wave with positive group velocity
  vₙ = sin(kₙdx)/dx, since E + i0 gives |λ| → 1⁻. Then
  Im Σₙ = −t_x sin(kₙdx) < 0 and Γₙ = i(Σ − Σ†) = 2t_x sin(kₙdx) = vₙ/dx.
* **Closed channel** (|cₙ| > 1, here cₙ > 1): λₙ = cₙ − √(cₙ² − 1) ∈ (0, 1).
  This is the solution that decays away from the device. Σₙ = −t_xλₙ is
  **real and negative**, and Γₙ = 0.
* **Closed channels are kept in Σ.** They carry no flux, but their real
  self-energy shifts the evanescent modes, and those couple to the open
  channels through V. Test `test_closed_channel_self_energy_is_included_and_matters`
  shows that dropping the closed-channel Σ changes G_N1 by more than 10⁻³
  relative on a random device.
* All of this comes from `smatrix_spec.lead_modes`; the solver contains no
  second copy. The branch selection is checked in `test_smatrix_spec.py`
  (roots of the quadratic at E + 10⁻⁹i; Im Σ < 0).

## 3. Recursive algorithm
Block Gaussian elimination from the left. With A_j = E − H_jj − Σ_L δ_{j1} − Σ_R δ_{jN}:

  D_1 = A_1,  D_j = A_j − t_x² g_{j−1},  g_j = D_j⁻¹,
  G_{j1} = −t_x g_j G_{j−1,1},  G_{1j} = −t_x G_{1,j−1} g_j,
  G_{11} ← G_{11} + t_x² G_{1,j−1} g_j G_{j−1,1},  G_{NN} = g_N.

These are the Schur-complement identities of the inverse of a block
tridiagonal matrix (left-connected Green functions).

* Each D_j is LU-factorized once (`scipy.linalg.lu_factor`). g_j and
  g_jG_{j−1,1} come from one `lu_solve` with stacked right-hand sides.
* The full device matrix is never inverted. The dense inverse is used only
  in tests, on small random toys, as a coding check.
* G_{1N} is computed by its **own** recursion, not set to G_{N1}ᵀ, so S = Sᵀ
  is a genuine numerical test.
* The reciprocal 1-norm condition number of every D_j is estimated
  (LAPACK `gecon`) and returned (`min_rcond`). A singular D_j or a
  non-finite block raises an error.
* Nothing is clipped or renormalized.
* Cost: O(N·N_y³) per energy. Measured (1 BLAS thread) per energy: 0.5 s
  (N_y = 63, N = 718), 2.3 s (N_y = 127, N = 718), 4.7 s (N_y = 127, N = 1435),
  2.4 s (N_y = 63, N = 3483).

Coding check against a dense inverse (random real-symmetric toy devices:
15 slices, 6 modes, 3 seeds × 3 energies): all four blocks agree to
≤ 1×10⁻¹² relative (measured ≤ 3×10⁻¹⁵), and G_{1N} = G_{N1}ᵀ to 10⁻¹².

## 4. Fisher–Lee implementation and phase convention
`smatrix_spec.fisher_lee` and `assemble_S` are called unchanged
(convention P, reference planes at slices 1 and N, [outgoing, incoming]):

  t_RL = iΓ^{1/2}G_{N1}Γ^{1/2},  t_LR = iΓ^{1/2}G_{1N}Γ^{1/2},
  r_L = iΓ^{1/2}G_{11}Γ^{1/2} − I,  r_R = iΓ^{1/2}G_{NN}Γ^{1/2} − I,

  S = [[r_L, t_LR], [t_RL, r_R]],

restricted to the open channels (identical leads). Incidence from both
directions is computed. Reciprocity is checked as S = Sᵀ, i.e.
t_LR = t_RLᵀ. Same-direction symmetry is **not** assumed: a random toy
device violates it by > 10⁻³, as expected. `scattering_matrix` requires at
least one closed channel in the basis (otherwise open channels would be
truncated) and rejects threshold energies (via `lead_modes`).

## 5. B1: analytic expectations (clean guide, V = 0)
A clean device with exact leads is the infinite clean lattice guide.
Everything is diagonal in modes, for every channel n including closed ones:

  G_{11,nn} = 1/(t_x(λₙ − 1/λₙ)),  G_{N1,nn} = λₙ^{N−1} G_{11,nn}.

For an open channel G_{11,nn} = −i/Γₙ. Hence, under convention P,
**t_nn = λₙ^{N−1}** (free propagation from slice 1 to slice N) and r = 0.
So S = [[0, D], [D, 0]] with D = diag(λ^{N−1}). **S is not the identity.**
The diagonal entries t_nn are unimodular phases e^{i kₙ(N−1)dx}.

## 6. B1 test configurations and measured results
Device: the transport geometry, i.e. slab 1.9 … 1.9 + L_d with ≥ 0.2
buffers on each side, L_y = 1, from `Potential(waveguide_config(...))` with
no disorder and no barrier (V = 0 exactly). Energies: the spec grid, 50
values spanning the 1st–6th channel windows (E ∈ [4.93, 241.8]) with
|E − Eₙ^⊥| ≥ 10⁻³E. This includes energies with 1 open channel (E < 19.74)
and with 2–6 open channels. Every energy also has ≥ 57 closed channels
(N_y = 63) or ≥ 121 (N_y = 127). Both incidence directions are checked.

Thresholds are from the spec, **unchanged**. The reciprocity threshold
10⁻⁸ is the existing B2/B3 value; B1 itself has none.

| Configuration | N | \|T_nn−1\| ≤1e-12 | mixing ≤1e-12 | R_n ≤1e-12 | U ≤1e-10 | phase ≤1e-10 (prov.) | Q ≤1e-8 | min rcond | B1 |
|---|---|---|---|---|---|---|---|---|---|
| N_y 63, dx 1/512, L_d 1 | 718 | 1.5e-13 | 0 | 9.6e-24 | 1.5e-13 | 1.9e-11 | 4.4e-15 | 3.8e-3 | **pass** |
| N_y 127, dx 1/512, L_d 1 | 718 | 1.5e-13 | 0 | 9.6e-24 | 1.5e-13 | 1.9e-11 | 4.4e-15 | 1.8e-3 | **pass** |
| N_y 63, dx 1/1024, L_d 1 | 1435 | 4.2e-13 | 0 | 1.4e-22 | 4.2e-13 | 7.3e-11 | 5.9e-15 | 3.8e-3 | **pass** |
| N_y 127, dx 1/1024, L_d 1 | 1435 | 4.2e-13 | 0 | 1.4e-22 | 4.2e-13 | 7.3e-11 | 5.9e-15 | 1.9e-3 | **pass** |
| N_y 63, dx 1/512, L_d 3 | 1742 | 2.8e-13 | 0 | 7.3e-23 | 2.8e-13 | 4.6e-11 | 6.4e-15 | 3.8e-3 | **pass** |
| N_y 63, dx 1/1024, L_d 3 | 3483 | 6.4e-13 | 0 | 5.8e-22 | 6.4e-13 | **1.77e-10** | 8.5e-15 | 3.8e-3 | **fail (phase)** |

(Maxima over the 50 energies; exact values per configuration are in
`results/rgf_b1_validation/summary.json`, per-energy values in `b1_grid.csv`.)

**Additional B1 checks (tests, N_y = 63, N = 718):**
* At E = 10, 100 and 230 (1, 4 and 6 open channels), all criteria hold.
* Off-diagonal entries of G_11 and G_N1 are ≤ 10⁻¹² × max.
* Closed-channel diagonals match the infinite-chain values to ≤ 10⁻⁹
  relative; open-channel diagonals to ≤ 10⁻¹⁰ × max. In an exploratory run,
  G_N1 and G_11 matched the analytic blocks to 4.1×10⁻¹¹ and 1.4×10⁻¹¹
  relative to the largest entry.
* The open-channel S is identical (≤ 10⁻¹²) for N_y = 15, 63 and 127.

**The failed criterion (investigated, threshold not relaxed).** The
predeclared B1 phase quantity on the 3483-slice device is 1.7685×10⁻¹⁰
(> 10⁻¹⁰). The reproducible high-precision audit in §6a shows that this is a
**genuine RGF phase error**, not an artifact of the double-precision
reference. Against the converged high-precision lattice phase the RGF itself
deviates by up to 1.5425×10⁻¹⁰ (8 of 424 rows exceed 10⁻¹⁰). |t| is still
exact to ≤ 3.2×10⁻¹³ there.

> *Correction to the previous version of this section.* It stated the
> high-precision result for configuration F without a committed script or
> artifact that reproduced it (the committed `run_phase_exact()` evaluated
> only configuration A). It also said the double reference deviates "by
> only 2–3×10⁻¹¹" and that the dx scaling of the error matched a factor ×4.
> The committed audit (§6a) reproduces the 1.54×10⁻¹⁰ value. However:
> * the double reference deviates by up to 6.1×10⁻¹¹ over all F rows (the
>   2–3×10⁻¹¹ values held only at the worst-error energies);
> * the dx scaling is not a clean ×4 (details below).

**Mechanism (supported as a bound, not as an exact scaling law).** The
diagonal of D_j contains 2t_x ≈ 10⁶ (dx = 1/1024), while E − Eₙ^⊥ is small
near a threshold. Rounding of order ε_mach·2t_x gives a per-slice
wavenumber error of order ε/sin(kdx), which accumulates over N slices. The
audit data (per row, RGF vs high-precision reference) show:
* in every configuration, |Δφ| ≤ 0.89·N·ε_mach/sin(kdx) (median ratio
  0.21–0.23);
* at fixed dx the error grows exactly linearly with the number of slices:
  per row, E/A has median ×2.43 (IQR 2.41–2.43) for N 718 → 1742;
* halving dx does **not** multiply the per-row error by a fixed factor.
  Per row, C/A has median ×2.6 (IQR 1.2–5.6) and F/A median ×7.0 (IQR
  3.0–13.6), against the bound's ×4 and ×9.7. A different dx changes the
  rounding of every Hamiltonian entry, so only the bound scales
  predictably. The earlier "×3.9 vs ×4" agreement was based on maxima of
  the B1 quantity and was partly coincidental.

Rescaling H by 1/t_x does not remove the cancellation; a different
discretization would be needed, and the specification fixes this one. The
strict xfail `test_long_device_provisional_phase_criterion` records the
failure (the suite flags it if it changes). Note: the B3 plan (spec §3)
uses L_d = 1 with dx = 1/512 and 1/1024, and those configurations pass.

## 6a. Reproducible high-precision phase audit (all six configurations)
Script `audit_rgf_phase.py`. Committed artifacts in
`results/rgf_b1_validation/`:
* `phase_audit_hp.csv`: 2544 rows, i.e. 6 configurations × 50 energies ×
  open channels × 2 directions;
* `phase_audit_hp_summary.json`;
* `phase_audit_hp_metadata.json`: source commit `d5a97a3`, clean working
  tree, Python, NumPy, SciPy and mpmath versions, precision, the exact
  command, timestamp, code fingerprint.

The earlier `phase_exact.csv` (configuration A, L→R only) is left
unchanged.

**What is compared.** For the clean guide, convention P gives
t_nn = λₙ^{N−1} exactly (§5). Three quantities are kept apart:

| Quantity | Definition | Role |
|---|---|---|
| (i) B1 phase quantity | \|t_nn − λ_double^{N−1}\|: double-precision RGF vs double-precision reference (λ from `lead_modes`) | the **predeclared B1 criterion** (≤ 10⁻¹⁰), unchanged |
| (ii) RGF vs high-precision reference | Δφ = arg(t_nn·e^{−iφ_ref}) | the RGF's actual numerical phase error |
| (iii) double reference vs high-precision reference | arg(λ_double^{N−1}·e^{−iφ_ref}) | the error of the double-precision reference itself |

Magnitude errors ||t| − 1| are reported separately.

**High-precision reference.** φ_ref = (N − 1)·θ with θ = arccos(1 − δ),
δ = (E − Eₙ^⊥)·dx², evaluated in mpmath at 100 significant digits.
* Taken exactly as given (double values are exact binary rationals): E (the
  B1 grid value), Eₙ^⊥ (the same double array the RGF uses), dx and N. The
  reference is therefore the exact phase *of the lattice problem the RGF was
  given*. It is not the continuum phase. Differences from the continuum are
  discretization or parameter rounding, not RGF error.
* No double-precision intermediate of the RGF or of `lead_modes` (λ, k,
  angles) enters it.

**Phase-difference convention.** Δφ = arg(t·e^{−iφ_ref}) ∈ (−π, π], the
principal value of φ_RGF − φ_ref. It is independent of 2π branches. This
matters: the unwrapped φ_ref reaches 73 rad (≈ 23π) in configuration F, so
naive subtraction of principal values would fail across the ±π cut. A test
covers that case. It is valid while the true difference is below π; here it
is ~10⁻¹⁰. The mpmath reduction into (−π, π] is done at an explicit
precision, not the ambient one; a test checks this.

**Precision and conditioning of the reference.** Near a threshold forming
1 − δ loses about 11 of the 100 digits, and arccos amplifies relative error
by ~1/√(2δ) ≲ 10⁶; N − 1 ≤ 3482 adds about 3.5 more digits. Measured over
all rows:

| Check | Max over all rows |
|---|---|
| dps 30 vs dps 100 | ≤ 1.1×10⁻²⁹ |
| dps 50 vs dps 100 | ≤ 1.1×10⁻⁴⁹ |
| arccos form vs cancellation-free 2·arcsin(√(δ/2)) (dps 100) | ≤ 9.1×10⁻¹⁰⁰ |

So the reference resolves the 10⁻¹⁰ effects with ≥ 19 orders of magnitude
to spare. It distinguishes RGF error, double-reference error and branch
effects unambiguously.

**Results (all six configurations; each 50 energies, 424 rows):**

| Config | N | (ii) RGF vs HP: max \|Δφ\| | RMS | at E, channel | (iii) double ref vs HP: max | max \|\|t\|−1\| | (i) B1 quantity max | B1 phase criterion |
|---|---|---|---|---|---|---|---|---|
| A (63, 1/512, 1) | 718 | 2.360e-11 | 3.25e-12 | 79.247, n=4 | 2.35e-11 | 7.4e-14 | 1.878e-11 | pass |
| B (127, 1/512, 1) | 718 | 2.360e-11 | 3.25e-12 | 79.247, n=4 | 2.35e-11 | 7.4e-14 | 1.878e-11 | pass |
| C (63, 1/1024, 1) | 1435 | 6.358e-11 | 1.35e-11 | 125.692, n=5 | 2.52e-11 | 2.1e-13 | 7.289e-11 | pass |
| D (127, 1/1024, 1) | 1435 | 6.358e-11 | 1.35e-11 | 125.692, n=5 | 2.52e-11 | 2.1e-13 | 7.289e-11 | pass |
| E (63, 1/512, 3) | 1742 | 5.765e-11 | 7.92e-12 | 79.247, n=4 | 5.70e-11 | 1.4e-13 | 4.560e-11 | pass |
| **F (63, 1/1024, 3)** | 3483 | **1.5425e-10** | 3.29e-11 | 125.692, n=5 | 6.12e-11 | 3.2e-13 | **1.7685e-10** | **fail** |

(Exact values are in `phase_audit_hp_summary.json`.)
* N_y does not matter for the clean guide: A ≡ B and C ≡ D bitwise.
* The recomputed B1 quantity agrees with `b1_grid.csv` to rounding (test).
* Applying the same 10⁻¹⁰ value to (ii) as a *diagnostic* gives the same
  verdicts: F fails (1.54×10⁻¹⁰), the others pass.

**Interpretation.**
1. **Configuration F fails B1 on its genuine RGF error, not on reference
   error.** At F's worst row the RGF error is +1.542×10⁻¹⁰ and the double
   reference error is −2.26×10⁻¹¹. They add up to the B1 quantity
   1.769×10⁻¹⁰, but the RGF error alone already exceeds 10⁻¹⁰.
2. **The predeclared B1 phase quantity (i) is a comparison of two
   double-precision computations.** It can over- or under-state the RGF
   error. It can even *mask* it when both share the same rounding error. At
   E's worst row the RGF and the double reference deviate from the exact
   phase by +5.77×10⁻¹¹ and +5.70×10⁻¹¹, so (i) reports only 7×10⁻¹³. In
   these data this never changes a verdict. It is a methodological
   limitation of criterion (i) as declared, which is kept unchanged.
3. Phase wrapping plays no role in any verdict: all |Δφ| are ~10⁻¹⁰ ≪ π.

**Reproducibility.** Rerunning configurations A and F reproduced their 848
rows identically (text-exact). Regenerating all six from commit `d5a97a3`
gave a CSV byte-identical to an earlier full run. Rows are deterministic on
this machine and library set (`OPENBLAS_NUM_THREADS=1`).

Reproduce: `python audit_rgf_phase.py` (≈ 10 min, all six) or
`python audit_rgf_phase.py --configs F --out DIR`; then
`python -m pytest -q tests/test_phase_audit.py`.

## 7. Numerical stability study (measured)
* **Conditioning.** min rcond(D_j) over the spec grid is 1.8–3.8×10⁻³ in all
  configurations (smaller for N_y = 127). It does not degrade with device
  length.
* **Number of transverse modes.** N_y = 15, 63, 127 give the same
  open-channel S (≤ 10⁻¹²). Errors are identical for N_y = 63 and 127.
* **Device length.** Magnitude errors grow mildly with N (|T−1|: 1.5e-13 at
  N = 718 to 6.4e-13 at N = 3483). Unitarity and reciprocity stay at
  ≤ 6.4e-13 and ≤ 8.5e-15. The phase error grows as described in §6.
* **Near (not at) thresholds** (outside the B1 grid; measured, not
  pass/fail; N_y = 63, N = 718, E = Eₙ^⊥(1 ± δ), n = 1…7):

  | δ | max \|T_nn−1\| | max U | max R | max phase err. | Q |
  |---|---|---|---|---|---|
  | 10⁻³ | 1.3e-13 | 1.3e-13 | 1.4e-20 | 1.2e-10 | 2.7e-15 |
  | 10⁻⁵ | 1.6e-13 | 1.6e-13 | 1.7e-17 | 4.1e-9 | 4.4e-15 |
  | 10⁻⁷ | 1.1e-13 | 1.1e-13 | 4.2e-16 | 2.1e-8 | 8.3e-15 |
  | 10⁻⁹ | rejected by `lead_modes` (\|cₙ\| within 10⁻¹² of 1), as specified | | | | |

  Magnitudes, unitarity and reciprocity stay at the 10⁻¹³ level up to the
  threshold. Only the phase degrades, as ε/sin(kdx). min rcond near the
  thresholds is 4.9×10⁻⁶: worse, but far from singular. 13 of the 52
  near-threshold points (all δ = 10⁻⁹) are rejected.
* **BLAS threads (performance, measured).** With the default multithreaded
  OpenBLAS (4 threads), one energy at N = 718 took 16–22 s, because the
  63×63 BLAS calls suffered quantized ~4–16 ms thread wake-up latency. With
  `OPENBLAS_NUM_THREADS=1` (or 2) it took 0.4–0.6 s. **For the clean guide
  the S matrices for 1, 2 and 4 threads are bit-identical** (E = 30, 100,
  200; N_y = 63, N = 718). The test session and the validation script
  therefore set `OPENBLAS_NUM_THREADS=1` (setdefault).
* **A refuted hypothesis (recorded for transparency).** I first attributed
  the slowness to subnormal floating-point intermediates. These are indeed
  ~40× slower in an isolated product, but an explicit flush of tiny
  entries removed only 22 entries, did not speed anything up, and left S
  bit-identical. Per-slice timings were also erratic and independent of
  magnitude. That flush step was therefore removed; the solver modifies no
  values.

## 8. Numerical limitations
* **B1 is internal consistency only.** The clean guide checks the lattice
  algebra, lead self-energies, recursion, Fisher–Lee normalization, phase
  convention and reciprocity bookkeeping. It does **not** test the
  potential-matrix elements in transport, mode mixing by a potential, the
  agreement with the continuum (B2), or with the time-dependent solver (B4).
* **The predeclared B1 phase quantity compares two double-precision
  computations** (§6a). It can mask or exaggerate the RGF error by up to
  ~6×10⁻¹¹ here. The high-precision audit is the measurement of the
  RGF's own phase error.
* **Provisional phase criterion fails for N = 3483** (§6). The phase error
  scales as L/(k·dx²) near thresholds. Phase-sensitive applications on long
  fine grids need this taken into account. The error in |t| is ≤ 6.4×10⁻¹³
  throughout.
* The 10⁻¹² magnitude thresholds have their smallest margin on the longest
  device: |T_nn − 1| = 6.4×10⁻¹³, about 1.6×.
* Threshold energies within ~10⁻¹²/dx² of Eₙ^⊥ are rejected by design.
* Bit-identity across BLAS thread counts was verified only for the clean
  guide and one device size. It is not established for disordered devices.
* All measurements are on one machine and library set (environment recorded
  in `results/rgf_b1_validation/run_metadata.json`).
