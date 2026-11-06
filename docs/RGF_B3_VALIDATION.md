# Benchmark B3: RGF vs an independent direct sparse Green-function reference on fixed disorder realizations

Phase 2, benchmark B3 of [PHASE2_METHOD_SELECTION.md](PHASE2_METHOD_SELECTION.md)
§3. Date: 2026-11-06.

* Production code commit: `d3f9291` (`sparse_reference.py`, `validate_rgf_b3.py`,
  `tests/test_rgf_b3.py`, additive `rgf_b3` entry in `provenance.py`). The B3 artifacts were
  generated from that commit with a clean working tree.
* Exploratory diagnosis script: `explore_b3_dx.py` (commit `a0616cf`).
* Artifacts: `results/rgf_b3_validation/`.
* Unchanged: `stationary.py`, `smatrix_spec.py`, `ode1d.py`, every B1/B2 threshold and
  artifact, all baseline results. No predeclared threshold was changed.

Status words used below:
**PASS** / **FAIL**: a predeclared criterion, evaluated as specified.
**EXPLORATORY**: computed to understand a result, not part of the predeclared protocol,
never replaces a criterion.
**NOT TESTED**: nothing was computed.

## 0. Outcome

**B3 does not pass: 9 of 10 criterion groups pass, and one fails.**

| Criterion | Result | Limit | Status |
|---|---|---|---|
| RGF vs sparse, full complex S | max 4.9×10⁻¹² | 10⁻¹⁰ | **PASS** |
| Unitarity U (RGF / sparse) | 4.9×10⁻¹⁴ / 2.6×10⁻¹³ | 10⁻¹⁰ | **PASS** |
| Reciprocity Q: full, r_L, r_R, t_LR vs t_RLᵀ, probabilities (RGF / sparse) | 2.6×10⁻¹⁴ / 2.1×10⁻¹³ | 10⁻⁸ | **PASS** |
| N_y 63 → 127: max \|ΔT₁\| | 1.8×10⁻⁶ | 2×10⁻⁴ | **PASS** |
| N_q 127 → 254: max \|ΔT₁\| | 1.9×10⁻⁶ | 2×10⁻⁴ | **PASS** |
| Buffer 0.2012 → 0.4004: max \|ΔT₁\| | 7.6×10⁻¹¹ | 2×10⁻⁴ | **PASS** |
| Same realization (bitwise potential identity) | 28 of 28 checks bitwise equal | bitwise | **PASS** |
| **dx 1/512 → 1/1024: max \|ΔT₁\|** | **5.97×10⁻⁴** (case D, E = 140.47) | **2×10⁻⁴** | **FAIL** |

The failure is real and well characterized (§7.2, §8.2). It is a discretization effect of the
specified finite-difference scheme in the strongest-disorder case (W = 200). It is **not** a solver
error: the sparse reference agrees with the RGF to 10⁻¹² at dx = 1/1024 as well.
Consequence: at the base resolution dx = 1/512, T₁ is not converged to 2×10⁻⁴ for the
strong-disorder realization (§8.2, §9).

B3 is a numerical-consistency validation on four fixed realizations. It makes no statement
about disorder-averaged transport, about the ζ dependence, or about localization.

---

## 1. Objective and scope
Test, for one fixed disorder realization per configuration, that
1. the RGF solver and a differently-implemented Green-function solution of the same discretized
   problem give the same complex S-matrix;
2. S is unitary and reciprocal, as the physics requires;
3. T₁ is stable under the numerical refinements the specification lists.

Not within scope: the physical model, the Fisher–Lee convention, or the comparison with the
time-dependent solver (B4). A pass of the first item validates the Green-function algebra and its
numerical implementation, **not** the physical model (§4).

## 2. Mathematical formulation (taken from the repository, spec §2.7)
Units ħ = m = 1, hard-wall strip of width L_y = 1.
* **Basis and Hamiltonian.** Slices x_j = x_start + (j−1)dx, j = 1…N. Transverse modes
  χₙ(y) = √(2/L_y) sin(nπy/L_y), n = 1…N_y, E⊥ₙ = (nπ/L_y)²/2, t = 1/(2dx²):
  H_jj = 2t·I + diag(E⊥) + V_j, H_{j,j±1} = −t·I.
* **Potential blocks.** (V_j)ₙₘ = Σ_l χₙ(y_l)V(x_j, y_l)χₘ(y_l)·L_y/(N_q+1), y_l = lL_y/(N_q+1),
  l = 1…N_q, N_q ≥ 2N_y (DST-I quadrature).
* **Leads.** Channel n: λ + 1/λ = 2cₙ, cₙ = 1 − (E − E⊥ₙ)dx². Open (|cₙ| < 1): λ = e^{ikₙdx},
  retarded root, Γₙ = 2t sin(kₙdx) = vₙ/dx. Closed: λ real, |λ| < 1, Γₙ = 0. Self-energy
  Σ = −t·diag(λ), identical on both sides.
* **Device Green function.** G = M⁻¹, M = E·I − H_D − Σ_L ⊕ Σ_R (Σ on the first and last slice
  blocks). Needed blocks: G₁₁, G_N1, G₁N, G_NN.
* **Fisher–Lee, convention P** (all indices [outgoing, incoming], open channels only):
  t_RL = iΓ_R^{1/2}G_N1Γ_L^{1/2}, t_LR = iΓ_L^{1/2}G₁NΓ_R^{1/2},
  r_L = iΓ_L^{1/2}G₁₁Γ_L^{1/2} − I, r_R = iΓ_R^{1/2}G_NNΓ_R^{1/2} − I.
  S = [[r_L, t_LR], [t_RL, r_R]] (rows: outgoing L, R; columns: incoming L, R).
* **Quantities.**
  * T₁(E) = Σₘ|(t_RL)ₘ₁|², the total transmission for incidence from the left in channel 1.
  * U = max over both products |SᴴS − I| and |SSᴴ − I|.
  * Q = max|S − Sᵀ|, which implies the block relations r_L = r_Lᵀ, r_R = r_Rᵀ, t_LR = t_RLᵀ.
  * **Probability reciprocity (index mapping).** T^{L→R}_{m←n} = |t_RL[m,n]|² (in: left n, out:
    right m) equals T^{R→L}_{n←m} = |t_LR[n,m]|² (in: right m, out: left n). So the array
    comparison is `|t_RL|² − (|t_LR|²)ᵀ`; the transpose is the index mapping. The same-direction
    asymmetry |t_RL[m,n]|² vs |t_RL[n,m]|² is **report-only** (spec §2.7.8): it does not follow
    from time reversal, and measures 0.49 here.

## 3. Shared conventions between the two methods
Both methods represent the same discretized physical problem. They share by construction:
* the lattice Hamiltonian above (t, E⊥, the DST-I quadrature definition of V_j);
* the clean-lead roots, Σ and Γ: `smatrix_spec.lead_modes`;
* the Fisher–Lee amplitudes and the S assembly: `smatrix_spec.fisher_lee`, `assemble_S`;
* the channel ordering (open channels ascending) and the left/right block orientation;
* the potential callable `Potential(cfg, disorder_from_config(cfg, seed)).numpy`.

Agreement therefore says nothing about whether these conventions are the physically right ones.
They were verified separately against a Green-function-free transfer-matrix calculation
(`tests/test_smatrix_spec.py`).

## 4. Computational independence
| Step | RGF (`stationary.py`) | Sparse reference (`sparse_reference.py`) |
|---|---|---|
| Transverse basis, y points | `transverse_basis` | own `transverse_data` (own sine matrix) |
| Sampling of V(x_j, y_l) | `build_device` | own coordinates and sampling |
| V_j blocks | einsum + symmetrization | own per-slice matmul, no symmetrization |
| Device matrix | blocks used one at a time | **one** assembled `scipy.sparse` matrix M of order N·N_y |
| Solution method | left-to-right block recursion (Schur complements, LAPACK `getrf`/`getrs` on 63×63 blocks) | **one** SuperLU factorization (`splu`, COLAMD) and direct triangular solves for 2N_y boundary right-hand sides |
| Boundary blocks | recursion for G₁₁, G_N1, G₁N, G_NN | columns of X = M⁻¹B; G₁N is a separate solve, not G_N1ᵀ |
| Conditioning check | reciprocal condition numbers of the D_j | relative residual of M X = B |
| Dense inverse | none | none |
| Lead data, Fisher–Lee, S | shared (§3) | shared (§3) |

Enforced by tests:
* `sparse_reference.py` does not import `stationary` and never refers to the RGF recursion;
* with `stationary.rgf_endpoint_blocks` monkey-patched to raise, the sparse reference still runs,
  and the same patch breaks the RGF itself;
* the assembled sparse M equals an explicit dense block-tridiagonal M, to 8 machine epsilons
  (not bitwise: the diagonal is summed in a different order);
* on synthetic problems the sparse solver reproduces
  (i) a Green-function-free backward transfer-matrix wave function (1D chain, complex amplitudes
  t and r to 10⁻⁹, non-trivial scatterer with T = 0.28);
  (ii) the exact clean-guide S;
  (iii) the RGF on a random six-channel device with genuine channel mixing.

This is an independent numerical implementation of the Green-function calculation, **not** an
independent physical model.

## 5. Configurations, energies, discretization
Predeclared configurations (all seed 0, L_d = 1; disorder region (1.9, 2.9), edge width 0.02,
512 random modes):

| Case | W | ζ | disorder parameter hash (sha256, first 20 hex) | potential hash at dx = 1/512 |
|---|---|---|---|---|
| A | 100 | 0.02 | `f022790469f9d…` | `6f66a2ad277cfe4…` |
| B | 100 | 0.05 | `f47346b43bf5a…` | `07717723b62164f…` |
| C | 100 | 0.20 | `2122ece8feba8…` | `fe7636039b79b13…` |
| D | 200 | 0.05 | `f47346b43bf5a…` | `be2d1b6953e6a20…` |

**B and D are the same random draw at two amplitudes**: the same seed and ζ give identical random
parameters, and V_D = 2·V_B bitwise. Likewise A, B and C are the same draw rescaled by ζ. The four
cases are therefore not independent samples.

**Energies** (five per case; the specification names the packet-energy convention but no
list, so this rule was fixed before any B3 result existed): kx = kx₀ + j σ_k, j = −2…2,
σ_k = 1/(√2 σ_x) = 3.5355, E = (kx² + π²/L_y²)/2.

| j | kx | E | open channels | nearest threshold |
|---|---|---|---|---|
| −2 | 12.929 | 88.51 | 4 | n = 4, 10.8 % away |
| −1 | 16.464 | 140.47 | 5 | n = 5, 12.2 % |
| 0 | 20.000 | 204.93 | 6 | n = 6, 13.3 % |
| +1 | 23.536 | 281.90 | 7 | n = 8, 12.0 % |
| +2 | 27.071 | 371.36 | 8 | n = 9, 7.6 % |

The threshold rule of spec B1/B4 (|E − Eₙ| ≥ 10⁻³E) moves none of them. No energy, channel, seed or
configuration was dropped.

**Base discretization:** dx = 1/512, N_y = 63, N_q = 2N_y + 1 = 127, so the device has
N = 719 slices and the sparse matrix has order 45 297. The S-matrices are 8×8 … 16×16.
**Geometry.** Slices lie on the lattice 1.9 + k·dx. The device extends n_b·dx beyond the disorder
region on each side, with the buffer length fixed physically to ⌈buffer·512⌉/512, i.e. 0.2012
(103/512) and 0.4004 (205/512). The same length is used at every dx. The disorder region (1.9, 2.9)
is the same for every buffer, and the devices are nested on one lattice. The potential on the first
lead slices is ≤ 6.3×10⁻⁷ (guard in `build_device`: 10⁻⁶, not loosened).

*A flaw found and fixed before the production run.* My first geometry rounded the buffer up at each
dx separately, so the dx = 1/1024 device had a 0.2002 buffer instead of 0.2012, and a unit test of
the nesting property failed. A dx refinement must not also change the buffer, so the buffer is now
defined once on the coarse lattice. No result had been produced with the earlier geometry.

**Refinements** (one parameter at a time, against the base):
dx = 1/1024; N_q = 254; buffer 0.4004; N_y = 127 with the default N_q = 255.
*Note:* the N_y criterion changes N_q along with N_y, because N_q ≥ 2N_y is required. The pure
transverse truncation (N_y = 63 vs 127 at fixed N_q = 254) was computed separately and is
labelled EXPLORATORY (§8.1).

## 6. Predeclared criteria
All from spec §3 B3 (values unchanged): max|S_RGF − S_sparse| ≤ 10⁻¹⁰ (complex, full S);
U ≤ 10⁻¹⁰ (both products, both methods); Q ≤ 10⁻⁸ (full, three blocks, probabilities);
|ΔT₁| ≤ 2×10⁻⁴ for N_y 63→127 (N_y = 255 only if the ζ = 0.02 case fails), dx, N_q and buffer
refinements; potential samples bitwise identical to a direct call. Same-direction transmission
symmetry is report-only.

## 7. Results

### 7.1 Same realization (PASS)
For every device built (28 checks: 4 cases × {base, dx, N_q, buffer, N_y = 127 (two variants)} for
the RGF, plus the sparse reference's own sampling for each base device), the sampled array V(x_j, y_l)
equals, **bitwise**, the array from a direct call of
`Potential(cfg, disorder_from_config(cfg, seed=0)).numpy` on independently built coordinates,
with a freshly drawn disorder object. Maximum difference 0.0, and the coordinate arrays are also
bitwise equal. Shapes are (719, 127) at the base, hashes are in `potential_identity.csv`.
Potential ranges at the base: A [−383, 380], B [−432, 346], C [−203, 194], D [−863, 692].

### 7.2 RGF vs sparse reference, full complex S (PASS)
20 comparisons (4 cases × 5 energies), S dimensions 8–16 (2·N_open), 719 slices.

| Case | max \|S_RGF − S_sparse\| over the 5 energies |
|---|---|
| A | 2.8×10⁻¹² |
| B | **4.9×10⁻¹²** (E = 88.51, in r_R) |
| C | 2.7×10⁻¹² |
| D | 2.5×10⁻¹² |

Margin ≥ 20× against 10⁻¹⁰. Every block (r_L, t_LR, t_RL, r_R) agrees to ≤ 4.9×10⁻¹². The sparse
solves have relative residuals ≤ 8.6×10⁻¹⁶ (LU fill 5.8–5.9×10⁶ non-zeros for 2.9×10⁶ in M). The
smallest RGF reciprocal condition number of the D_j is 5.7×10⁻³ (case D, E = 88.5). T₁ from both
methods agrees to ≤ 6.1×10⁻¹³.

### 7.3 Unitarity (PASS)
Maximum over all 20 points; both products are evaluated and reported.

| Method | max \|SᴴS − I\| | max \|SSᴴ − I\| | limit |
|---|---|---|---|
| RGF | 4.85×10⁻¹⁴ | 4.05×10⁻¹⁴ | 10⁻¹⁰ |
| sparse | 2.59×10⁻¹³ | 2.00×10⁻¹³ | 10⁻¹⁰ |

### 7.4 Reciprocity (PASS)
| Method | full max\|S − Sᵀ\| | r_L | r_R | t_LR − t_RLᵀ | probabilities | limit |
|---|---|---|---|---|---|---|
| RGF | 1.88×10⁻¹⁴ | 1.08×10⁻¹⁴ | 1.88×10⁻¹⁴ | 1.65×10⁻¹⁴ | 2.56×10⁻¹⁴ | 10⁻⁸ |
| sparse | 2.12×10⁻¹³ | 1.04×10⁻¹³ | 2.12×10⁻¹³ | 1.85×10⁻¹³ | 1.83×10⁻¹³ | 10⁻⁸ |

The RGF computes G₁N independently of G_N1, so S = Sᵀ is a genuine test of the recursion. Same-direction
transmission asymmetry (report-only): max 0.49 for both methods. This is expected, since it is not
a consequence of time reversal and **not** a criterion.

### 7.5 Numerical refinements of T₁ (|ΔT₁| ≤ 2×10⁻⁴)
Maximum over the 5 energies of each case:

| Refinement | A (ζ 0.02) | B (ζ 0.05) | C (ζ 0.20) | D (W 200) | Status |
|---|---|---|---|---|---|
| N_y 63 → 127 (N_q 127 → 255) | 1.8×10⁻⁶ | 4.0×10⁻⁷ | 4.7×10⁻⁸ | 1.6×10⁻⁶ | **PASS** |
| N_q 127 → 254 | 1.9×10⁻⁶ | 4.0×10⁻⁷ | 4.7×10⁻⁸ | 1.6×10⁻⁶ | **PASS** |
| buffer 0.2012 → 0.4004 | 2.7×10⁻¹¹ | 7.6×10⁻¹¹ | 4.7×10⁻¹¹ | 7.0×10⁻¹¹ | **PASS** |
| **dx 1/512 → 1/1024** | 7.0×10⁻⁵ | **1.81×10⁻⁴** | 3.4×10⁻⁵ | **5.97×10⁻⁴** | **FAIL** (D) |

* **Closed channels** retained: 55–59 at N_y = 63 and 119–123 at N_y = 127 (open channels 4–8).
* **N_y = 255:** the rule requires it only if the ζ = 0.02 case fails at N_y = 127. Case A passes
  (1.8×10⁻⁶), so it was **not run** (NOT TESTED).
* The N_y change is dominated by the accompanying change of N_q (§8.1).

dx refinement, all 20 points (T₁ at dx = 1/512, 1/1024, change):

| Case | E | T₁(1/512) | ΔT₁ = T₁(1/1024) − T₁(1/512) | status |
|---|---|---|---|---|
| A | 88.51 / 140.47 / 204.93 / 281.90 / 371.36 | 0.42634 / 0.83118 / 0.73164 / 0.70983 / 0.81891 | −1.0×10⁻⁶ / 1.6×10⁻⁵ / 7.0×10⁻⁵ / 1.1×10⁻⁵ / 2.1×10⁻⁵ | all pass |
| B | same | 0.42933 / 0.34864 / 0.49970 / 0.63517 / 0.74025 | −9.8×10⁻⁵ / 2.7×10⁻⁵ / −5.6×10⁻⁵ / −2.0×10⁻⁵ / **1.81×10⁻⁴** | all pass |
| C | same | 0.54684 / 0.85415 / 0.90949 / 0.97371 / 0.99586 | −3.4×10⁻⁵ / −9.2×10⁻⁶ / −7.5×10⁻⁶ / −9.6×10⁻⁶ / −3.3×10⁻⁶ | all pass |
| D | same | 0.41483 / 0.30208 / 0.10372 / 0.23938 / 0.39468 | 2.1×10⁻⁵ / **−5.97×10⁻⁴** / **2.80×10⁻⁴** / −1.27×10⁻⁴ / **−2.80×10⁻⁴** | **3 of 5 fail** |

(`dx_convergence.csv` has the full-precision values.)

## 8. Numerical convergence analysis

### 8.1 Transverse truncation and quadrature
* Raising N_q alone (127 → 254) changes T₁ by up to 1.9×10⁻⁶ (case A).
* The criterion run (N_y 63→127 with N_q 127→255) differs from that by at most 9.3×10⁻⁸ per
  energy, so the N_y entry of the criterion is almost entirely the N_q change.
* **EXPLORATORY:** the pure truncation effect (N_y 63 vs 127, both at N_q = 254) is ≤ 9.5×10⁻⁸
  (case A, E = 371.36), and ≤ 1.3×10⁻⁸ for B, C and D.
* This shows that the spec's expectation that ζ = 0.02 needs the most transverse modes holds only
  weakly here: with 55–59 closed channels retained, the ζ = 0.02 realization is converged in N_y
  to < 10⁻⁷ at these five energies.

### 8.2 The dx failure: EXPLORATORY diagnosis (case D only)
Script `explore_b3_dx.py`, files `exploratory_dx_case_D.csv` and `…_summary.json`. This is **not**
part of the predeclared protocol and does not replace the criterion. Case D was chosen because it
is the failing case, and only it was extended (cases A–C were not).

RGF T₁ at dx = 1/256, 1/512, 1/1024, 1/2048 (same N_y = 63, N_q = 127, buffer):

| E | T₁(1/512) | T₁(1/1024) | T₁(1/2048) | Δ(512→1024) | Δ(1024→2048) | orders |
|---|---|---|---|---|---|---|
| 88.51 | 0.414827 | 0.414848 | 0.414853 | 2.1×10⁻⁵ | 5.2×10⁻⁶ | 1.98, 1.99 |
| 140.47 | 0.302077 | 0.301481 | 0.301332 | −5.97×10⁻⁴ | −1.49×10⁻⁴ | 2.01, 2.00 |
| 204.93 | 0.103725 | 0.104004 | 0.104074 | 2.80×10⁻⁴ | 7.0×10⁻⁵ | 1.99, 2.00 |
| 281.90 | 0.239380 | 0.239254 | 0.239222 | −1.27×10⁻⁴ | −3.2×10⁻⁵ | 2.03, 2.01 |
| 371.36 | 0.394676 | 0.394396 | 0.394327 | −2.80×10⁻⁴ | −7.0×10⁻⁵ | 2.02, 2.01 |

* The observed order is **1.98–2.03** at all five energies, i.e. in the window where Richardson
  extrapolation is permitted (applied here for diagnosis only). The error is the lattice-dispersion
  error of the second-order scheme and converges cleanly.
* With the extrapolated limit, the dx = 1/512 error of T₁ is up to **8.0×10⁻⁴** (4× the 2×10⁻⁴ limit)
  and the dx = 1/1024 error is up to 2.0×10⁻⁴, i.e. just at the limit.
* **The sparse reference at dx = 1/1024 agrees with the RGF** to 8.8×10⁻¹³ in T₁ and 4.9×10⁻¹² in
  S. Both methods see the same dx dependence; it belongs to the discretization, not to a numerical
  error of the RGF.
* *Plausible cause (a consistent explanation, not proved):* the lattice dispersion error
  ≈ (k dx)²/24 grows with the local wave number, and for W = 200 the potential reaches −863, so the
  local k = √(2(E − V)) reaches ≈ 45. The phase error accumulated across the slab is then
  ~10⁻² rad, which changes T₁ at the 10⁻⁴ level.
* For case B the largest Δ(512→1024) is 1.81×10⁻⁴. **If** the second-order behaviour measured in D
  also holds there (not tested), the dx = 1/512 error would be ≈ 2.4×10⁻⁴, also above 2×10⁻⁴.
  This is an extrapolation and is **not** established.

## 9. Limitations
* **Fixed realizations, five energies each.** B3 says nothing about ensemble statistics, the ζ
  dependence, or localization. B and D are the same draw at two amplitudes (§5).
* **Convergence tested only against the specified one-step refinements.**
  * N_y = 255: NOT TESTED (rule not triggered).
  * dx finer than 1/1024: only case D, EXPLORATORY.
  * dx behaviour of A, B, C beyond the criterion: NOT TESTED.
* **The energy set is my choice** within the specified convention (§5), fixed before any result.
  A different set might give a different largest dx change. The failing criterion should not be
  read as the worst case.
* **Shared conventions** (§3): the comparison cannot detect an error in the lead formulas or the
  Fisher–Lee convention; those rest on `tests/test_smatrix_spec.py` and B1/B2.
* **The sparse reference was run at the base configuration only** (plus case D at dx = 1/1024,
  exploratory). It was not run for the refinements. The refinements use the RGF alone, but the sparse
  agreement at the base and at dx = 1/1024 gives no indication of an RGF error.
* **Direct sparse LU is practical only for moderate sizes** (order 4.5×10⁴, ≈ 5 s per energy); it
  was chosen for transparency, not speed.
* **B1 is still not fully passed** (its provisional phase criterion fails on a 3483-slice device).
  The B3 devices have 719–1437 slices, and the RGF–sparse S agreement (5×10⁻¹²) shows no sign of
  accumulated round-off at this size.
* Bit-identity was shown on one machine (x86_64, 4 cores, OpenBLAS with one thread).

**Implication for B4 (recommendation, protocol not changed).** The B4 error budget in the spec
assumes an RGF dx error ≈ 10⁻⁴. For W = 200 at dx = 1/512 this is not met (≈ 8×10⁻⁴ at one energy).
Before B4 the dx convergence should be established for each B4 realization, and either dx = 1/1024
(≈ 2×10⁻⁴) or a verified Richardson extrapolation should be used. No such change has been made here.

## 10. Reproducibility
```bash
python validate_rgf_b3.py                  # 5.8 min on 4 cores; writes results/rgf_b3_validation/
python explore_b3_dx.py                    # EXPLORATORY dx diagnosis (case D), 1.5 min
python validate_rgf_b3.py --record-tests   # runs the B3 and the full test suites, records the outcome in run_metadata.json
python -m pytest -q tests/test_rgf_b3.py
```
* A second full run into a separate directory reproduced **every CSV field identically** (the sparse
  wall-clock column excluded). The same machine class was used.
* `run_metadata.json` records: commit `d3f9291`, clean tree, Python 3.13.18, NumPy 2.5.3,
  SciPy 1.18.1, x86_64, the four configurations and all 20 energies, the discretization of
  every variant, per-device geometry and potential hashes, the sparse method
  (`scipy.sparse.linalg.splu`, SuperLU with COLAMD), the RGF implementation
  (`stationary.rgf_endpoint_blocks`), sha256 hashes of the source files (including
  `stationary.py` and `smatrix_spec.py`), and the test results below.

## 11. Tests (actual commands and exit statuses)
| Command | Exit status | Result |
|---|---|---|
| `python -m pytest -q tests/test_rgf_b3.py` | 0 | 35 passed |
| `python -m pytest -q tests` (full suite) | 0 | 224 passed, 2 xfailed (the two known expected failures) |

Both were run by `validate_rgf_b3.py --record-tests` after the artifacts existed. Before that, with
the artifacts absent, the live tests gave 28 passed and 6 skipped. The tests cover:
* configuration and seed selection; energies; geometry and nesting;
* bitwise same-realization sampling, including detection of a 10⁻¹⁵ perturbation;
* sparse matrix construction and independence from the RGF recursion;
* sparse solver correctness on synthetic problems (1D transfer matrix, clean guide, random multichannel);
* the full complex S comparison (imaginary-only and phase-only errors are detected);
* both unitarity products and the full, block-wise and probability reciprocity (each block checked individually);
* pass/fail handling exactly at the thresholds (equality passes, +10⁻¹² relative fails, NaN fails);
* all four refinements, one parameter at a time;
* artifact integrity: row counts, summary vs CSV, recomputed pass flags, hashes, metadata, current source hashes,
  one comparison row recomputed live, and the exploratory results' consistency with the main run.

The B3 test file contains 35 tests. Two earlier test failures during development were bugs in my own
tests (an over-strict bitwise comparison and a too-weak scatterer) and were fixed in the tests, not by
changing any threshold or solver.
