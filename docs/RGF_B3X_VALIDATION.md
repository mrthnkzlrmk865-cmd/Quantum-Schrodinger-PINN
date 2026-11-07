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
