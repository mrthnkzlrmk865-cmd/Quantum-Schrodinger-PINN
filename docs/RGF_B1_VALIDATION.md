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
10⁻¹⁰ (§6). The threshold was not relaxed.

Reproduce:
```bash
python -m pytest -q tests/test_stationary.py       # ~40 s
python validate_rgf_b1.py                          # ~11 min on 4 CPU cores; --quick: ~1 min
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

**The failed criterion (investigated, threshold not relaxed).** The phase
error on the 3483-slice device is a genuine RGF accuracy limit, not an
artifact of the double-precision reference λ^{N−1}. Against the exact
lattice phase evaluated in 50-digit arithmetic (mpmath), the RGF deviates by
up to 1.54×10⁻¹⁰, while the double reference deviates by only 2–3×10⁻¹¹.
|t| is still exact (||t| − 1| ≤ 2×10⁻¹³).

Mechanism: the diagonal of D_j contains 2t_x ≈ 10⁶ (dx = 1/1024), while the
physical combination E − Eₙ^⊥ is small near a threshold. Rounding of order
ε_mach·2t_x gives a per-slice wavenumber error δ(kdx) ≈ ε/sin(kdx), which
accumulates over N slices. The predicted scaling N/sin(kdx) ∝ L/(k dx²)
matches the measurements: ×3.9 for dx/2 (predicted ×4), ×2.4 for L_d = 3
(predicted ×2.4), ×9.4 for both (predicted ×9.7).

The error is concentrated in channels with sin(kdx) < 0.01. In the base
configuration, comparing against the exact phase: ≤ 2.4×10⁻¹¹ for
sin(kdx) < 0.01, ≤ 8×10⁻¹² for 0.01–0.03, ≤ 4×10⁻¹² for 0.03–0.1
(`phase_exact.csv`). Rescaling H by 1/t_x does not remove the cancellation;
a different discretization would be needed. The specification fixes this
one. The strict xfail `test_long_device_provisional_phase_criterion`
records the failure (the suite will flag it if it changes). Note: the B3
plan (spec §3) uses L_d = 1 with dx = 1/512 and 1/1024, and those
configurations pass.

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
