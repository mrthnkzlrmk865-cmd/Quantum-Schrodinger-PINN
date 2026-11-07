# B3X: extended independent validation of the RGF under fixed disorder realizations

Phase 2, extension of benchmark B3 ([RGF_B3_VALIDATION.md](RGF_B3_VALIDATION.md), unchanged).
Date: 2026-11-07. Machine-readable protocol: [rgf_b3x_protocol.json](rgf_b3x_protocol.json).

**Part I (this part) was written and committed before any B3X number was computed.** Part II (results)
is appended after the run.

## Part I: Protocol

### I.1 Scope and relation to B3
B3 (four configurations, 20 energies) already exists, with one failing criterion (dx). B3X is a
separate campaign with the full 2 × 3 matrix requested for this extension:

| W | ζ = 0.02 | ζ = 0.05 | ζ = 0.20 |
|---|---|---|---|
| 100 | W100_z0.02 | W100_z0.05 | W100_z0.2 |
| 200 | W200_z0.02 (new) | W200_z0.05 | W200_z0.2 (new) |

Seed 0, L_d = 1, 5 energies each (30 points), dx = 1/512, N_y = 63, N_q = 127. It writes only to
`results/rgf_b3x_validation/`; no B1, B2 or B3 artifact is touched.

**Not independent samples.** One seed is shared by all six cases: different ζ rescale one random draw,
and W = 200 is exactly twice W = 100 (bitwise). B3X validates numerical consistency of fixed
realizations. It says nothing about ensemble statistics, ζ dependence or localization.

### I.2 Problem (taken from the repository, spec §2.7)
ħ = m = 1, strip of width L_y = 1 with hard walls; transverse modes χₙ = √(2/L_y) sin(nπy/L_y),
E⊥ₙ = (nπ/L_y)²/2; t = 1/(2dx²). On slices x_j = x₀ + (j−1)dx, j = 1…N, in the mode basis
n = 1…N_y: (H c)_j = −t c_{j−1} + (2t + E⊥ + V_j) c_j − t c_{j+1}, with
(V_j)ₙₘ = ∫χₙ V(x_j, y) χₘ dy evaluated by the DST-I quadrature y_l = l L_y/(N_q+1), l = 1…N_q, N_q ≥ 2N_y.
Leads (j ≤ 0, j ≥ N+1) are clean.

### I.3 Independent derivation of the reference formulation ("wave-function reference")
The reference computes the S-matrix **from the scattering wave function**, not from Fisher–Lee
formulas on Green-function blocks.

*Lead modes.* In a clean lead channel n: −t c_{j−1} + (2t + E⊥ₙ − E)c_j − t c_{j+1} = 0, c_j = λʲ,
so λ + 1/λ = 2cₙ, cₙ = 1 − (E − E⊥ₙ)dx². Open (|cₙ| < 1): λ = cₙ ± i√(1−cₙ²) = e^{±ikₙdx}. For
E → E + i0⁺ the root with |λ| < 1 is the retarded one; linearizing in the imaginary part of E gives
dλ = −η dx² e^{iθ}/sin θ (radially inward for θ ∈ (0, π)), so the retarded open root is e^{+ikₙdx}, a wave
moving away from the device. Closed (|cₙ| > 1): both roots are real, the retarded one is the one
with |λ| < 1. The wave carries the lattice flux J = 2t Im(c*_j c_{j+1}) = |a|² 2t sin(kₙdx) =: |a|² Γₙ, which
defines the flux normalization of an amplitude.

*Scattering state for incidence in left channel n.* In the left lead,
c_j = a λₙ^{j−1}eₙ + Σₘ bₘ λₘ^{−(j−1)} eₘ (incident + outgoing, all channels m including evanescent ones),
and in the right lead c_j = Σₘ dₘ λₘ^{j−N} eₘ. At the device ends the lead solution continues onto
slices 1 and N: c₁ = a eₙ + b, c_N = d. The equation at slice 1 reads
(E − E⊥ − V₁ − 2t)c₁ + t c₀ + t c₂ = 0. With c₀ = Λ(c₁ − a eₙ) + λₙ⁻¹ a eₙ this becomes
(E − H₁₁ + tΛ)c₁ + t c₂ = −t(λₙ⁻¹ − λₙ) a eₙ. The same elimination at slice N for the right lead
(c_{N+1} = Λ c_N) gives (E − H_NN + tΛ)c_N + t c_{N−1}. Therefore the scattering state solves the sparse system

  M c = s,  M = E I − H_D + t Λ at the two end blocks (Σ = −tΛ),  s₁ = −t(λₙ⁻¹ − λₙ) a eₙ, s_j = 0 otherwise.

For an open channel s₁ = i Γₙ a eₙ. The source is computed from this formula, not from a Γ shortcut.
Unit incident flux fixes a = Γₙ^{−1/2}. The flux-normalized amplitudes (open m only) are
r_{mn} = Γₘ^{1/2}(c₁ − a eₙ)ₘ, t_{mn} = Γₘ^{1/2}(c_N)ₘ. Incidence from the right is the mirror problem with the
source at slice N and r_R = Γ^{1/2}(c_N − a eₙ), t_LR = Γ^{1/2} c₁. The S-matrix is assembled in the order
[left open channels ascending, right open channels ascending], rows = outgoing, columns = incoming.

*Why this is independent.* It never forms a Green-function block, never uses Γ^{1/2} G Γ^{1/2} − I
(including the "−I" term) and never calls `smatrix_spec`. Agreement with the RGF therefore tests the
Fisher–Lee convention and the recursion together, against a direct physical definition of the amplitudes.

### I.4 What is shared and what is independent
| Item | RGF (`stationary.py` + `smatrix_spec.py`) | Old B3 reference (`sparse_reference.py`) | B3X wave-function reference (`sparse_wavefunction.py`) |
|---|---|---|---|
| Lattice Hamiltonian (formulas of §I.2) | yes | yes | yes (re-implemented) |
| Potential callable | `Potential(...).numpy` | same | same |
| Lead roots, Σ, Γ | `smatrix_spec.lead_modes` | same function | **own closed form**, retarded sign checked against the E+iη root |
| Transverse basis | `transverse_basis` | own sine matrix | own; V_j by cosine sums g_k = Σ v_l cos(πkl/(N_q+1)), V_nm = (g_{|n−m|} − g_{n+m})/(N_q+1) |
| Linear algebra | block recursion (Schur complements) | one SuperLU factorization, COLAMD | SuperLU with MMD_AT_PLUS_A ordering + one step of iterative refinement |
| S from | Fisher–Lee on G blocks | `fisher_lee`, `assemble_S` on G blocks | **own wave-function decomposition** (§I.3) |

Shared by all three: the physical problem definition (§I.2) and the potential realization. A mistake in
that definition would not be detected by B3X (it is checked in B1, B2 and the spec tests).

### I.5 Pitfalls identified in the audit, and how each is handled
1. *Name collision:* B3 exists → B3X, separate directory, no B3 file modified.
2. *Shared conventions in the old reference:* addressed by the new reference (§I.4).
3. *Retarded root and source sign:* derived (§I.3) and tested against the E + iη root; unitarity is a second check.
4. *Closed channels:* kept in the device basis; Σ real, Γ = 0 exactly; reported open/closed counts must equal the analytic count from the thresholds.
5. *N_y vs N_q:* N_q ≥ 2N_y couples them; the N_y criterion is evaluated at fixed N_q = 254 (single variable), as in B3 Revision 1.
6. *Grid nesting:* buffer length fixed physically (ceil(buffer·512)/512) so a dx or buffer refinement changes one thing only.
7. *Clean-lead guard (10⁻⁶):* W200_z0.02 has 8.3×10⁻⁷ at buffer 0.2012 (pre-flight geometry check, no transmission computed). The guard is not loosened.
8. *Short ζ vs sampling:* ζ = 0.02 is ≈ 2.6 quadrature points per ζ in y at N_q = 127; the N_q criterion is the safeguard.
9. *dx convergence:* B3 already failed this criterion at W = 200. It is predeclared at the spec value and will be reported as found. The diagnostic dx ladder (1/256 … 1/2048) is labelled diagnostic and cannot rescue a failed criterion.
10. *Potential identity:* "bitwise" only where the evaluation path is identical (RGF sampling vs a direct call). The two independent algebraic routes to V_j differ by rounding and are compared with a relative tolerance.

### I.6 Acceptance criteria (fixed now; values and justifications in `rgf_b3x_protocol.json`)
| ID | Quantity | Limit | Why this tolerance |
|---|---|---|---|
| C1 | same realization: sampled V bitwise equal to a direct call | bitwise | identical evaluation path |
| C2 | potential blocks, matrix-product route vs cosine-sum route | 10⁻¹² · max\|V\| | rounding-level difference of two exact formulas |
| C3 | open-channel count equal in all methods and equal to the analytic count | exact | integer |
| C4 | closed channels: Γ = 0 exactly, λ real with \|λ\| < 1; ≥ 1 closed channel retained | exact | definitions |
| C5 | max\|S_RGF − S_wavefunction\| (complex, full S) | 10⁻¹⁰ | spec value; ~20× above observed 5×10⁻¹² |
| C6 | max\|S_RGF − S_greenblock\| | 10⁻¹⁰ | same |
| C7 | transmission and reflection agreement: max\|\|S_ij\|² − \|S'_ij\|²\|, and \|T1 − T1'\|, \|R1 − R1'\| | 2×10⁻¹⁰ | ≤ 2 max\|dS\| |
| C8 | off-diagonal channel transmission agreement (all m ≠ n) | 2×10⁻¹⁰ | same |
| C9 | unitarity max\|SᴴS − I\| and max\|SSᴴ − I\|, all methods | 10⁻¹⁰ | spec value |
| C10 | reciprocity: full \|S − Sᵀ\|, blocks r_L, r_R, t_LR vs t_RLᵀ, probabilities \|t_RL\|² vs (\|t_LR\|²)ᵀ, all methods | 10⁻⁸ | spec value |
| C11 | \|ΔT1\|, dx 1/512 → 1/1024 | 2×10⁻⁴ | spec value (B4 budget) |
| C12 | \|ΔT1\|, N_q 127 → 254 | 2×10⁻⁴ | same |
| C13 | \|ΔT1\|, buffer 0.2012 → 0.4004 | 2×10⁻⁴ | same |
| C14 | \|ΔT1\|, N_y 63 → 127 at fixed N_q = 254 (N_y = 255 only if a ζ = 0.02 case fails) | 2×10⁻⁴ | same |

Report-only: same-direction transmission asymmetry; dx ladder orders and Richardson estimates
(applied only where the measured orders lie in [1.8, 2.2]); consistency with the existing B3 numbers
for the four overlapping configurations. **B3X passes only if every criterion C1–C14 is evaluated and met.**

### I.7 Outputs
`results/rgf_b3x_validation/`: per-energy comparison CSV, per-configuration summary CSV,
convergence CSVs, `summary.json`, `run_metadata.json` (commit, Python and library versions, parameters,
seed, protocol hash and commit, source hashes, test commands and results). Reproduction commands are
in Part II.

---

## Part II: Results

Protocol commit `6cee803` (Part I and `rgf_b3x_protocol.json`, committed before any B3X number existed).
Code commit used for the production run: `5d907cd` (clean working tree). Artifacts:
`results/rgf_b3x_validation/`. Nothing in B1, B2 or B3 was modified.

Status words: **PASS** / **FAIL** = a predeclared criterion evaluated as specified; **DIAGNOSTIC** =
computed to understand a result, not a criterion; **NOT TESTED** = nothing was computed.

### II.0 Outcome

**B3X does not pass: 13 of 14 predeclared criteria pass, and C11 (longitudinal refinement) fails.**
All 14 criteria were evaluated on all 30 energy points (6 configurations × 5 energies).

| ID | Criterion | Max observed (where) | Limit | Status |
|---|---|---|---|---|
| C1 | same realization, bitwise (54 checks) | 0.0 difference | bitwise | PASS |
| C2 | potential blocks, two independent routes | 2.3×10⁻¹⁵ relative | 10⁻¹² | PASS |
| C3 | open-channel count (3 methods + analytic) | 0 mismatches (4–8 open) | exact | PASS |
| C4 | closed channels (Γ = 0, real λ, \|λ\| < 1, ≥ 1 retained) | Γ_closed = 0.0, max \|λ\| = 0.985, 55–59 retained | exact | PASS |
| C5 | max\|S_RGF − S_wave-function\| | 5.3×10⁻¹² (W200_z0.2, E = 88.5) | 10⁻¹⁰ | PASS |
| C6 | max\|S_RGF − S_Green-block\| | 5.0×10⁻¹² (W200_z0.2, E = 88.5) | 10⁻¹⁰ | PASS |
| C7 | transmission and reflection agreement | 2.6×10⁻¹² (W200_z0.2, E = 88.5) | 2×10⁻¹⁰ | PASS |
| C8 | off-diagonal channel transmission agreement | 1.1×10⁻¹² (W200_z0.2, E = 88.5) | 2×10⁻¹⁰ | PASS |
| C9 | unitarity, both products, 3 methods | 1.2×10⁻¹² (wave-function ref., W100_z0.2) | 10⁻¹⁰ | PASS |
| C10 | reciprocity: full, blocks, probabilities, 3 methods | 8.8×10⁻¹³ (wave-function ref.) | 10⁻⁸ | PASS |
| **C11** | **\|ΔT1\|, dx 1/512 → 1/1024** | **9.17×10⁻⁴** (W200_z0.2, E = 88.5) | **2×10⁻⁴** | **FAIL** |
| C12 | \|ΔT1\|, N_q 127 → 254 | 6.7×10⁻⁶ (W200_z0.02, E = 140.5) | 2×10⁻⁴ | PASS |
| C13 | \|ΔT1\|, buffer 0.2012 → 0.4004 | 1.3×10⁻¹⁰ (W200_z0.2, E = 140.5) | 2×10⁻⁴ | PASS |
| C14 | \|ΔT1\|, N_y 63 → 127 at fixed N_q = 254 | 6.9×10⁻⁷ (W200_z0.02, E = 140.5) | 2×10⁻⁴ | PASS |

`summary.json`: `overall_B3X_status = "FAIL: C11_dx_refinement"`, `b3x_full_pass = false`.
The N_y = 255 rule (run only if a ζ = 0.02 configuration fails C14) was not triggered: NOT TESTED.

What this does and does not show: the RGF reproduces the complex S-matrix of two differently implemented
sparse references to 5×10⁻¹², is unitary and reciprocal to 10⁻¹³ or better, and is insensitive to N_y,
N_q and the buffer at the 10⁻⁶ level. At the base resolution dx = 1/512 the transmission of the strong-disorder
configurations is **not** converged to 2×10⁻⁴. That is a property of the specified finite-difference
scheme, not of the RGF (§II.4). B3X validates numerical consistency of six fixed realizations. It
says nothing about disorder-averaged transport, ζ dependence or localization.

### II.1 Configurations, energies, discretization
Six configurations (W ∈ {100, 200} × ζ ∈ {0.02, 0.05, 0.20}), seed 0, L_d = 1; energies
E = 88.51, 140.47, 204.93, 281.90, 371.36 (kx = kx₀ + jσ_k, j = −2…2; 4, 5, 6, 7, 8 open channels; the nearest
threshold is ≥ 7.6 % away, so the threshold rule moved nothing). Base: dx = 1/512, N_y = 63, N_q = 127, 719
slices, sparse matrix order 45 297, S-matrices 8×8 … 16×16. Device buffers 0.2012 (103/512) and 0.4004
(205/512). The first-lead-slice potential is ≤ 8.3×10⁻⁷ at the base; for W200_z0.02 it is 9.1×10⁻⁷ in the
dx = 1/1024 criterion run and 9.6×10⁻⁷ in the dx = 1/2048 diagnostic run (the first lead slice moves closer to the
slab as dx shrinks): **within 4 % of the 10⁻⁶ guard**, which was not loosened.
The six cases are one random draw (W = 200 is exactly 2 × W = 100; ζ rescales the draw), so they are not independent samples.

### II.2 RGF vs the two references (C5–C8)
| Configuration | max\|ΔS\| RGF–wave-function | RGF–Green-block | wave-function–Green-block | max T/R discrepancy | max off-diagonal T discrepancy | max off-diagonal T itself |
|---|---|---|---|---|---|---|
| W100_z0.02 | 2.7×10⁻¹² | 2.8×10⁻¹² | 6.0×10⁻¹³ | 9.6×10⁻¹³ | 8.1×10⁻¹³ | 0.51 |
| W100_z0.05 | 4.9×10⁻¹² | 4.9×10⁻¹² | 5.1×10⁻¹³ | 1.4×10⁻¹² | 6.4×10⁻¹³ | 0.37 |
| W100_z0.2 | 2.8×10⁻¹² | 2.7×10⁻¹² | 8.0×10⁻¹³ | 1.4×10⁻¹² | 6.9×10⁻¹³ | 0.59 |
| W200_z0.02 | 3.8×10⁻¹² | 4.5×10⁻¹² | 1.0×10⁻¹² | 1.6×10⁻¹² | 4.2×10⁻¹³ | 0.26 |
| W200_z0.05 | 2.5×10⁻¹² | 2.5×10⁻¹² | 4.3×10⁻¹³ | 1.1×10⁻¹² | 5.1×10⁻¹³ | 0.23 |
| W200_z0.2 | 5.3×10⁻¹² | 5.0×10⁻¹² | 1.1×10⁻¹² | 2.6×10⁻¹² | 1.1×10⁻¹² | 0.37 |

All blocks (r_L, t_LR, t_RL, r_R) agree at this level. Channel mixing is strong (off-diagonal transmission up to
0.59), so the agreement of the off-diagonal elements is a real test. The linear-solve backward errors are
≤ 9×10⁻¹⁶ (sparse relative residual, after one refinement step for the wave-function reference); the smallest RGF reciprocal condition number is 5.7×10⁻³.

### II.3 Unitarity, reciprocity, channels (C3, C4, C9, C10)
| Method | max\|SᴴS − I\| | max\|SSᴴ − I\| | Q full | r_L | r_R | t_LR − t_RLᵀ | probabilities |
|---|---|---|---|---|---|---|---|
| RGF | 4.9×10⁻¹⁴ | 4.1×10⁻¹⁴ | 1.9×10⁻¹⁴ | 1.8×10⁻¹⁴ | 1.9×10⁻¹⁴ | 1.7×10⁻¹⁴ | 2.6×10⁻¹⁴ |
| wave-function reference | 1.2×10⁻¹² | 9.8×10⁻¹³ | 8.8×10⁻¹³ | 7.8×10⁻¹³ | 4.5×10⁻¹³ | 8.8×10⁻¹³ | 6.6×10⁻¹³ |
| Green-block reference | 3.6×10⁻¹³ | 5.0×10⁻¹³ | 2.1×10⁻¹³ | 2.1×10⁻¹³ | 2.1×10⁻¹³ | 1.8×10⁻¹³ | 1.8×10⁻¹³ |

Same-direction transmission asymmetry (report-only, not a criterion): up to 0.49.
* **Channels:** 4–8 open and 55–59 closed channels are retained. The open-channel count agrees in all
  three methods and with the analytic threshold count. Γ of closed channels is exactly 0; closed λ are real
  with \|λ\| ≤ 0.985; the closed-form retarded root agrees with the E + iη root to 4.5×10⁻¹³ and with
  `smatrix_spec.lead_modes` to 2.2×10⁻¹⁵.
* **Closed-channel content (report-only; added after the protocol was fixed):** evanescent components carry up to 16 %
  of the open-channel amplitude norm at the first device slice, so closed channels are not negligible near the
  device. The RGF's closed-channel amplitudes (from its Green blocks) equal the reference wave function's to
  9×10⁻¹⁵ (amplitude scale 2.4×10⁻³).
* The unitarity check is sensitive to a sign error in the reference: flipping the source sign gives U ≈ 6 on the random six-channel test device.

### II.4 Convergence (C11–C14)
Maximum |ΔT1| per configuration:

| Configuration | dx 1/512→1/1024 (C11) | N_q 127→254 (C12) | buffer 0.2012→0.4004 (C13) | N_y 63→127, N_q = 254 (C14) |
|---|---|---|---|---|
| W100_z0.02 | 7.0×10⁻⁵ | 1.9×10⁻⁶ | 2.7×10⁻¹¹ | 9.5×10⁻⁸ |
| W100_z0.05 | 1.81×10⁻⁴ | 4.0×10⁻⁷ | 7.6×10⁻¹¹ | 1.8×10⁻⁹ |
| W100_z0.2 | 3.4×10⁻⁵ | 4.7×10⁻⁸ | 4.7×10⁻¹¹ | 5.8×10⁻¹¹ |
| W200_z0.02 | **2.58×10⁻⁴ (FAIL)** | 6.7×10⁻⁶ | 4.5×10⁻¹¹ | 6.9×10⁻⁷ |
| W200_z0.05 | **5.97×10⁻⁴ (FAIL)** | 1.6×10⁻⁶ | 7.0×10⁻¹¹ | 1.3×10⁻⁸ |
| W200_z0.2 | **9.17×10⁻⁴ (FAIL)** | 5.4×10⁻⁶ | 1.3×10⁻¹⁰ | 1.2×10⁻⁸ |

**C11 fails for 5 of 30 energy points** (all W = 200): W200_z0.02 at E = 281.9 (2.58×10⁻⁴); W200_z0.05 at
E = 140.5, 204.9 and 371.4 (5.97×10⁻⁴, 2.80×10⁻⁴, 2.80×10⁻⁴); W200_z0.2 at E = 88.5 (9.17×10⁻⁴).
The change grows with disorder strength and with the potential's steepness. N_y, N_q and the buffer are
converged by 2–5 orders of magnitude.

**DIAGNOSTIC: dx ladder (1/256, 1/512, 1/1024, 1/2048; not a criterion).**
* The observed order is **1.96–2.14** at all 30 points, so the error is second order, and a Richardson
  estimate is permitted by the predeclared window [1.8, 2.2] at all 30 points.
* The extrapolated error of the T1 at dx = 1/512 is up to 9.4×10⁻⁵ (W100_z0.02), 2.4×10⁻⁴ (W100_z0.05),
  4.6×10⁻⁵ (W100_z0.2), 3.4×10⁻⁴ (W200_z0.02), 8.0×10⁻⁴ (W200_z0.05) and **1.2×10⁻³ (W200_z0.2)**.
  At dx = 1/1024 it is up to 2.3×10⁻⁵, 6.1×10⁻⁵, 1.1×10⁻⁵, 8.6×10⁻⁵, 2.0×10⁻⁴ and **3.0×10⁻⁴**.
* **C11 compares consecutive grids; it does not bound the error.** For a second-order scheme the dx = 1/512
  error is about 4/3 of the consecutive change. W100_z0.05 passes C11 (1.81×10⁻⁴) yet its extrapolated
  dx = 1/512 error is 2.4×10⁻⁴. For W200_z0.2, even dx = 1/1024 is estimated to be 3.0×10⁻⁴ off, so a 1/1024 run
  would meet the consecutive-change criterion only marginally and an error target of 2×10⁻⁴ not at all.
  By second-order scaling dx = 1/2048 would be about 7.6×10⁻⁵ (extrapolated, not measured).
* **Likely cause (consistent with the data, not proved):** the lattice-dispersion error ≈ (k dx)²/24 grows with
  the local wave number, and the stronger disorder reaches deeper potentials.

**Implication for B4 (recommendation; the protocol was not changed).** The B4 budget assumes an RGF dx error of
≈ 10⁻⁴. That is not met at dx = 1/512 for W = 200 (and is marginal for W100_z0.05). B4 should establish dx
convergence per realization and use dx ≤ 1/1024 or a verified Richardson extrapolation, with the error, not the
consecutive change, as the target.

### II.5 Cross-check with the existing B3 (report-only)
For the four overlapping configurations (20 points) the RGF T1 equals the existing B3 values exactly (difference
0.0) and differs from the existing sparse reference's T1 by ≤ 6.1×10⁻¹³. The B3 artifacts were not modified.

### II.6 Limitations
* Fixed realizations, seed 0, five energies each. Not statistics. The six cases are one draw rescaled.
* **Shared by all three methods:** the problem definition (lattice Hamiltonian, DST-I quadrature, sine basis) and the
  potential callable. An error in that definition would not be detected here (B1, B2 and the spec tests cover it).
* The two references were run at the base configuration only; all refinements and the dx ladder use the RGF alone.
  The agreement at the base resolution gives no sign of an RGF error, but the references were not run at the finer dx.
* C11 is a consecutive-change criterion, not an error bound (§II.4).
* N_y = 255 NOT TESTED (not triggered). dx finer than 1/2048 NOT TESTED.
* The energy set is the project's convention, chosen before any result; other energies might give larger dx changes.
* The clean-lead tail of W200_z0.02 is within 4 % of its guard at dx = 1/2048 (9.6×10⁻⁷ vs 10⁻⁶); a finer dx would need a larger buffer.
* B1 is still not a complete pass (its longest configuration has a genuine phase discrepancy, confirmed
  independently at 100-digit precision); the B3X devices have 719–2873 slices.
* Items added after the protocol was fixed, all report-only: the closed-channel amplitude comparison and weights (§II.3).

### II.7 Reproducibility and tests
```bash
python validate_rgf_b3x.py                 # ~10 min on 4 cores; writes results/rgf_b3x_validation/
python validate_rgf_b3x.py --record-tests  # runs the B3X and the full suites, records them in run_metadata.json
python -m pytest -q tests/test_rgf_b3x.py
```
A second full run into a separate directory reproduced every field of the six CSV files identically (the
LU fill-count column was excluded from the comparison) and gave the same criterion outcomes. Same machine class.

`run_metadata.json` records: code commit, clean-tree flag, protocol file hash and commit, Python 3.13.18,
NumPy 2.5.3, SciPy 1.18.1, platform, all six configurations (seed 0, disorder parameter hashes), all 30 energies,
every variant, per-device geometry and potential hashes, the three methods, source-file hashes, and the test results.

| Command | Exit status | Result |
|---|---|---|
| `python -m pytest -q tests/test_rgf_b3x.py` | 0 | 23 passed |
| `python -m pytest -q tests` (full suite) | 0 | 252 passed, 2 xfailed (the two known expected failures) |

Both were recorded by `--record-tests` after the artifacts existed (no CI system was involved). Before the artifacts
existed, the live tests gave 17 passed and 6 skipped. The 23 tests cover: independence of the wave-function reference
(import and source inspection), retarded-root selection against the complex-energy root, cosine-sum vs matrix-product
potential blocks, clean-guide phase, a Green-function-free 1D transfer-matrix check, agreement with both other methods
on a random multichannel device, detection of a wrong source sign, preservation of every predeclared threshold (B1, B2,
B3 and B3X), configuration/seed/energy selection, same-realization bitwise sampling, buffer geometry, sensitivity of the
S comparison to phase/reflection/off-diagonal errors, artifact completeness, the summary equalling `summarize()` of the
CSVs, classification at the exact thresholds on mutated data, metadata/protocol-hash integrity, and a live recomputation
of one row. One test bug found during development (a missing key in a bookkeeping check) was fixed in the test.
