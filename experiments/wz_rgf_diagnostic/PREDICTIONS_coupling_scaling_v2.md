# Design and predictions: W/Z coupling-scaling experiment v2 (independent; preregistration document)

Revision 3.4 (document-only; on top of Revision 3.3 `4d471b0`; earlier: Revision 3.2 `f16e7c2`, Revision 3.1 `ae05c0d`). Revision 3 (`d40960b`) implements the corrections of the
independent R0 design review of `f101224` (verdict: CONDITIONAL PASS). Earlier revisions:
- Revision 2 (`ef4bf12`, `e52cd51`, `f101224`): the certificate audit;
- Revision 1 (`b62ed42`): the scientific audit of `01bf018`.

This is a design only:
- nothing in this document has been implemented or run;
- no code exists for it, no matrix of it has been generated, and no previous file is modified.

Relation to the previous experiment (`PREDICTIONS_coupled_blocks.md`, design `7d69260`, implementation `63e046d`, results
`e3c67a6`):
- Its decision, **DESIGN INSUFFICIENT**, is final and is not re-scored, re-fitted or reinterpreted here.
- Its post-hoc 3-point slopes (0.985–0.989) are **not** evidence for anything in this design and are not used to choose any
  threshold.
- Its matrices are not part of this experiment.
- This experiment **cannot** attribute the previous 0.78–0.79 slopes to any specific cause. Those configurations are
  excluded, and the previous fitted quantity (a block maximum with s-dependent normalization) differs from v2's. v2 tests
  only whether the candidate mechanisms occur in a new population.

Standing status:
- B1 configuration F still FAILS the unchanged 1e-10 criterion.
- W/Z is experimental and NOT validated.
- No production change is authorized.

## N. Notation (Revision 3; used throughout)

| symbol | meaning |
|---|---|
| s | coupling multiplier, V_j(s) = s·V_j(1) (unchanged definition) |
| L(s), R_W(s), R_orig(s) | represented lattices: exact; W/Z-represented; production-represented |
| class R, class H, class D | potential classes. Always written with the word "class" |
| M(s), M^(j)(s) | full device matrix and the j-slice prefix matrix (left lead only), physical units |
| A0, A1 | for an affine matrix A(s) = A0 + s·A1 |
| M1 … M16, M4L, M4R | mutant labels. Never matrices |
| Δ_j | W/Z variable (E·I − diag(E_perp) − V_j)·dx² |
| Λ, λ_n, γ_n | lead roots; γ_n = 1 − δ_n is the lead parameter |
| Σ_L, Σ_R | lead self-energies. A bare Σ is a summation |
| Ĝ = t·G | dimensionless Green function; Ĝ^(j) is the prefix inverse (M^(j)/t)⁻¹ |
| χ_j, ψ_j | left and right prefix chain blocks Ĝ^(j)_j1 and Ĝ^(j)_1j |
| θ_j | ψ_{j−1} P_j χ_{j−1} |
| P_j, W_j, Z | W/Z propagator, W/Z variable, closing matrix |
| SC_j | Schur complement of slice j in M^(j)/t, so that W_j = I + SC_j |
| A_mid, A_e | exact dyadic midpoint matrix and enclosure error, A = A_mid + A_e |
| ε_A | validated upper bound of ‖A_e‖_F |
| Y_κ | approximate inverse at certification point κ |
| Ξ_κ | exact residual I − A_mid(s_κ)·Y_κ |
| ω_κ | validated bound of ‖I − A(s_κ)Y_κ‖₂ |
| a_κ, h_κ, b, r_κ | validated bounds of ‖A(s_κ)⁻¹‖₂, of ‖A(s_κ)‖₂ and of ‖A1‖₂; certified radius |
| ū(·), ℓ | rigorous rational upper square root, with ℓ = 256 bits |
| ξ_n | radicand of the lead root |
| r_0 | common Taylor radius (§2.4 item 3) |
| l̄(·) | rigorous rational lower square root (§2.4.1 step 4) |
| η_A, β̃_i, E_m | per-matrix reference-error bound; per-block norm bound; polynomial error bound (§5.6) |
| U_m(s), G_lo(s) | upper bound on D*_true and lower bound on G_true, rigorous by construction; validation pending T1–T6 and R0 (§5.6) |
| D*_true(s), G_true(s) | exact-lattice R_W oracle discrepancy and ‖Ĝ_full‖₂ (not computed; enclosed) (§5.6) |
| SNR_ideal(s) | idealised exact-level ratio D*_true/(10⁻¹³·G_true) (§5.6) |
| S_test | the 43 preregistered design-grid points with s ≤ 2⁻¹⁹ (§5.6) |
| THRESHOLD | the exact rational `Fraction(1, 10**10)` (§5.6, T6.7) |
| n_tot | matrix dimension N·N_y |
| D_m(s), D*_m(s) | double and exact (oracle) discrepancy of the designated block X_m, on R_W |
| η_X, ν_X, ν_X,m | observed correct-code noise; SNR noise level; mutant-specific apparatus tolerance level |
| Mag_m(s) | magnitude scale of the oracle monomials of mutant m (§4.3) |
| τ_id | identity threshold (§4.7) |
| C_k, q, μ, μ_est | Taylor coefficients; asymptotic order; first-order remainder constant and its estimate |
| T_m, deg_m, α, β, ζ | majorant constants (§5.2) |
| φ(s) | local log-log slope |
| dev(s) | crossover deviation |
| Ω_asym | weak-coupling window |
| ϑ_i | OLS weights |
| p̂, p̂* | double and exact OLS slopes |
| Fact 1–7, Hyp1–4, Asm1–3, Gap1–6 | facts, hypotheses, assumptions, previous gaps |
| set XP, set HO | exploratory and held-out configuration sets |
| WZ solver, orig solver, Ref(·) | solvers under test and the direct-inverse reference |

## R. Revision log

Revision 1 (vs `01bf018`):

| # | change | audit item |
|---|---|---|
| R1.1 | Removed the restriction q ∈ {1, 2, 3}. q is any finite integer ≥ 1 (identity theorem), identified from Taylor coefficients computed by truncated Neumann-series algebra with a precision-scaling zero test, up to K_max = 12 | 1 |
| R1.2 | A non-integer finite-scale slope estimate is not INVALID. INVALID requires a contradiction between exact computations or with a theorem. Unresolved order is DESIGN INSUFFICIENT | 1 |
| R1.3 | Added a continuous-interval certificate (Neumann-radius covering of [0, 1]). Its inverse-norm bound was validated only from R2.1 on. Analytic claims are restricted to the certified interval; R_W statements are pointwise | 2 |
| R1.4 | Error quantities are classified as rigorous/validated, empirical, heuristic or fixed criterion (§4.6). Exceeding a heuristic bound alone is an anomaly (DESIGN INSUFFICIENT), never FAIL | 3 |
| R1.5 | FAIL is split into FAIL-W (W/Z accuracy, the 1e-10 criterion) and FAIL-D (a diagnostic hypothesis falsified; not evidence of a W/Z defect) | 3 |
| R1.6 | Specification-derived oracle; clean-room transcription; review gate R0 | 4 |
| R1.7 | Every mutant's exact code change is documented (§4.2) | 4 |
| R1.8 | Claims are separated by class (the physical scope is class R plus the class-D controls; class H is a separate stress-test outcome) | 5 |
| R1.9 | Sampling is conditional on the acceptance filter | 6 |
| R1.10 | Legacy diagnostic set removed | scope |
| R1.11 | Fixed asymptote points replaced by Taylor coefficients | 1 |
| R1.12 | Coverage counted per class | 5 |

Revision 2 (vs `b62ed42`):

| # | change |
|---|---|
| R2.1 | Validated approximate-inverse certificate with an exactly computed residual. It replaces an mp value inflated by (1 + 1e-30), which was not a proof |
| R2.2 | Every bound is stated in the spectral norm, computed through ‖·‖₂ ≤ ‖·‖_F with exact Frobenius norms |
| R2.3 | Special case b = 0 |
| R2.4 | Certificate failure and cost policy; no silent empirical fallback |
| R2.5 | Mandatory certificate self-tests; R0 extended to the certificate code |
| R2.6 | Taylor disk from the validated certificate. The oracles (except M10) are polynomials in inverse blocks, so no series inversion is needed |
| R2.7 | Validated reference residuals |
| R2.8 | The zero test is labelled non-validated |
| R2.9 | Ratio-test tail bound, with proof |
| R2.10 | λ-enclosure radius ≤ 2⁻²⁵⁵ |
| R2.11 | ū defined piecewise, ū(0) = 0 |
| R2.12 | "ε_A = 0 iff A_e = 0" replaced by the implication ε_A = 0 ⟹ A_e = 0 |

Revision 3 (vs `f101224`; R0 review items E1–E17). **No sample count, grid, fitting window, mutant definition, threshold
value or the 1e-10 criterion is changed.** Every change to a decision rule is listed explicitly here and in §6.

| # | review item | change | rule change? |
|---|---|---|---|
| R3.1 | E1 | Coverage uses certified points only. An uncertified point is discarded, never given radius 0 (§2.4 item 3) | no (clarification) |
| R3.2 | E2 | Common Taylor radius r_0 := min over all matrices of r_0,A = 1/(2β), β = max_A a_A,0·b_A (§2.4, Fact 3, §5.2) | no |
| R3.3 | E3 | Lead perturbation A_e specified separately for the full M (N ≥ 2), for N = 1, and for prefix matrices (§2.4.1) | no |
| R3.4 | E5 | The validated Y is the **exact** dyadic value of the mpmath result, with no further rounding (§2.4.1, §4.7). The original wording, which also asserted a mantissa width equal to the working precision, was removed in Rev3.1-4 | no |
| R3.5 | E6 | Self-test T5 (exact-rational assembly, floating-point intermediates prohibited) | adds a precondition |
| R3.6 | E9 | ε_A explanation corrected. The one-way implication is kept. The example "A_e = 0 with ε_A > 0" is shown to be unreachable in this design and is not claimed to be realizable | no |
| R3.7 | E10 | M9 product written as an explicit ordered product (§4.2) | no |
| R3.8 | E11 | Rationale of §4.7 item 3 corrected: diagonal rounding is O(ε), independent of s, with effect O(ε/gap) | no |
| R3.9 | E4 | The two conflicting tail-comparison criteria of §5.2 (ii) are replaced by one criterion | **yes**: the redundant "relative 1e-40" is removed; the precision allowance's reference scale is defined (Mag_m) |
| R3.10 | E7 | Mutant-specific apparatus tolerance ν_X,m. **Double-vs-oracle and double-vs-transcription excess becomes a mutant-apparatus anomaly (DESIGN INSUFFICIENT), no longer INVALID.** Oracle-vs-transcription (mp vs mp) disagreement remains INVALID. The first-deviation check (also mp vs mp), previously listed without an explicit outcome, is now explicitly INVALID | **yes** |
| R3.11 | E8 | Agreement checks compare matrices directly, not discrepancy norms | **yes** (stricter) |
| R3.12 | E12 | A correct-code error reliably measured above 1e-10 is FAIL-W regardless of the judgeability gate. An error above 1e-10 that cannot be measured reliably → DESIGN INSUFFICIENT. The 1e-10 value is unchanged | **yes** (stricter) |
| R3.13 | E13 | Defined: identity threshold τ_id, meaningful-candidate pair, judgeable case, reliably measured | no |
| R3.14 | E14 | R_W is the lattice for D*, φ_exact, p̂*, SNR* and dev(s); ‖C_q‖ comes from L. The noise-trend window, previously "over Ω_asym" (undefined for a per-configuration control), is defined as all grid s ≤ 2⁻¹², per configuration and block, i.e. the mutant-independent part of Ω_asym | clarification; the window choice is flagged for R0 |
| R3.15 | E15 | Energy rounding is half-to-even in exact rationals. A maximum redraw-attempt count, with its consequence | **adds a rule**; the cap value is flagged as not derived (§9) |
| R3.16 | E16 | Notation disambiguated (§N). Labels renamed: R1.x/R2.x/R3.x, Fact/Hyp/Asm/Gap, sets XP/HO, WZ/orig solver, Ref(·) | no |
| R3.17 | E17 | Predicted risks added: q ≥ 2 pairs and zero-test scaling (§9) | no |
| R3.18 | review §F | Self-tests T1–T4 fully specified; T5 added (§2.4.2) | adds preconditions |
| R3.19 | (consequence of R3.12) | W_j references are computed in the validated Schur form W_j = I + SC_j, which is algebraically identical to I + (Ĝ^(j)_jj)⁻¹ (§4.2, §4.7) | no |
| R3.20 | E2 | s = 0 must itself be certified for every matrix, since the Taylor disk is built from the s = 0 certificate (previously implicit) | filter condition made explicit |

Revision 3.1 (vs `d40960b`). **No sample count, grid, fitting window, mutant, SNR threshold, zero-test threshold, numerical
tolerance or the 1e-10 criterion is changed.**

| # | change | rule change? |
|---|---|---|
| Rev3.1-1 | The series-vs-direct check (§5.2 (ii)) no longer triggers INVALID through the non-validated B_prec. On failure there is an independent dps-200 rerun; persistent failure → DESIGN INSUFFICIENT ("unresolved series inconsistency"). B_tail is unchanged (validated) | **yes**: INVALID → escalation then DESIGN INSUFFICIENT (§5.2, §4.6, §6) |
| Rev3.1-2 | Explicit binary64 import boundary (exact dyadic import of existing doubles; rational assembly only afterwards; direct float inputs to assembly rejected) (§2.4.1, T5) | no (clarification) |
| Rev3.1-3 | T5 wording: binary64 addition rounds the represented s-term; exact assembly preserves it. Replaces "would be lost" | no |
| Rev3.1-4 | Y precision wording: Y is the exact dyadic value of the mpmath result without subsequent rounding. The mantissa-width claim is removed (§2.4.1 step 1, R3.4) | no |
| Rev3.1-5 | Cleanup: §2.4 items 4 and 5 and the certificate failure policy moved into their structural position before §2.4.1; the duplicate "Validity on the interval" paragraph and a duplicate lead-enclosure bullet in §4.7 removed; R3.19/R3.20 ordered | no |
| Rev3.1-6 | §5.6: feasibility analysis of regression validity for q ≥ 2, with no rule change; the predicted risk is retained (§9 item 12) | no |

Revision 3.2 (vs `ae05c0d`). **No grid, SNR threshold, fitting window, point count, sample count, threshold, mutant or
decision rule is changed.**

| # | change | rule change? |
|---|---|---|
| Rev3.2-1 | §5.6 (C-q1) now applies to the SNR-qualified points actually included in Ω_asym (at least 8 points, span ≥ 100, SNR*_m ≥ 1e3 at each included point), not to every grid point of [2⁻¹⁹, 2⁻¹²]. The necessary consequence (a qualified point at s ≤ 2⁻¹⁹) is retained. The point-count bullet is corrected accordingly | no (the analysis is aligned with the existing rules) |
| Rev3.2-2 | §5.6 classification "Possible but unlikely" renamed "Likelihood not assessed"; no likelihood estimate is made | no |
| Rev3.2-3 | §5.2 (ii) escalation wording: B_tail validated, B_prec non-validated, the dps-200 rerun is an independent numerical consistency check, persistent disagreement is DESIGN INSUFFICIENT and not proof of a contradiction, no third precision level | no (wording; the rule is as in Rev3.1-1) |
| Rev3.2-4 | **Reviewer decisions confirmed**, as listed below | no (confirmation of existing values) |

**Reviewer decisions confirmed (Rev3.2-4):**
- the redraw cap of 10,000 attempts (§2.4 item 5);
- Mag_m as the reference scale of B_prec and of the mutant-specific heuristic tolerance τ_abs,m and ν_X,m, with their
  non-validated status preserved (§4.3, §5.2 (ii));
- the noise-trend window: all grid s ≤ 2⁻¹², without SNR filtering (§5.5; R3.14);
- the zero-test thresholds unchanged (§5.2);
- the q ≥ 2 feasibility risk retained (§5.6, §9 item 12);
- dps 200 as the final escalation level of the series check (§5.2 (ii)).

Revision 3.3 (vs `f16e7c2`). This is a mathematical consistency correction of the §5.6 diagnostic test only. **No
threshold, grid, sample count, mutant, fitting window, SNR condition, q classification, majorant formula or decision rule
is changed.**

| # | change | rule change? |
|---|---|---|
| Rev3.3-1 | §5.6 majorant test: it may conclude infeasibility only if r_0 ≥ 2⁻¹⁹ (so the span-derived necessary point s ≤ 2⁻¹⁹ lies inside the certified Taylor disk) and the evaluated grid subset is verified nonempty. If r_0 < 2⁻¹⁹ → INCONCLUSIVE, not infeasible. This replaces the unjustified "every grid s ≤ min(2⁻¹⁹, r_0)" conclusion | no (the diagnostic test has no decision role) |
| Rev3.3-2 | §5.6 heading "Mathematically impossible" renamed "Majorant-based infeasibility under the stated conditions"; the caveats are unchanged | no |

Revision 3.4 (vs `4d471b0`; final-audit corrections C1–C5). **No threshold, grid point, mutant definition, sampling rule or
decision rule is changed.**

| # | change | rule change? |
|---|---|---|
| Rev3.4-1 | §5.6: the L-majorant infeasibility claim is replaced by the pointwise certificate on R_W (U_m, G_lo, exact comparison; INCONCLUSIVE otherwise). It no longer depends on the 1e-30 restriction, §4.7 item 3, r_0 or the zero test. The L-majorant is retained as "L-lattice estimate only", with no R_W claim. The conclusion is restricted to the 43 preregistered design-grid points (C1) | no (diagnostic, no decision role) |
| Rev3.4-2 | √N_y replaced by ū(N_y); l̄ defined in §2.4.1 step 4 (the §4.7 S-reference bullet now refers to it); the threshold is the exact rational `Fraction(1, 10**10)` | no |
| Rev3.4-3 | Exact-level (SNR_ideal) and rule-level claims separated (C5); eligibility versus certificate inequality distinguished with respect to the zero test (C2); the diagnostic's non-decision role stated | no (clarification) |
| Rev3.4-4 | Self-test T6 added (§2.4.2), including C3 (exact threshold type) and C4 (supplementary irrational case reported as a numerical cross-check only). Precondition references updated from T1–T5 to T1–T6 (§2.4.1, §4.3, §8.1, §8.2, §9 item 7) | adds a precondition |
| Rev3.4-5 | Notation table and §4.6 row added. The certificate's status is "Designed; validation pending T1–T6 and R0" | no |

## 0. Inspection of the previous experiment: verified definitions and reproducibility gaps

Verified from `coupled_blocks.py` and `wz_mutants.py` at `63e046d`:

| item | exact definition found |
|---|---|
| coupling parameter s | `V_j(s) = v0 · s · M_j` with v0 = 4. **s multiplies the entire potential**, including the diagonal and identity terms. Values used: s ∈ {1, 2⁻⁴, 2⁻⁸, 2⁻¹²} |
| W/Z input | Δ_j = (E·I − diag(E_perp) − V_j(s))·dx² in double; leads from `lead_from_delta(fl(E − E_perp)·dx²)`; dx = 1/8, t = 32, L_y = 1 |
| mutant discrepancy D_mut | max over the four endpoint blocks of ‖X_mut − X_ref‖_F / ‖Ĝ_full(s)‖₂, against the R_W reference |
| slope | `np.polyfit(ln s, ln D_mut, 1)` over the 4 scales |
| "D ≫ τ_end" | D ≥ 100·τ_end |
| mutants | `wz_mutants.MUTANTS` (M1–M16, M4L, M4R). Pure order mutants M1, M2, M3, M5, M6, M7, M8, M15. `MIN_N` per mutant |

Previous gaps, documented and not repaired retroactively:
- Gap1. The fit procedure and "≫" existed only in code.
- Gap2. The maximum over blocks can switch block between s values.
- Gap3. The normalization depended on s.
- Gap4. One fit spanned both coupling regimes, with no crossover rule.
- Gap5. The commuting and non-commuting parts of V both scale with s.
- Gap6. No noise-free signal was computed.

v2 keeps **the same definition of s** (V_j(s) = s·V_j(1)) and the same W/Z definitions.

## 1. Research questions, mathematical facts, hypotheses

| RQ | question | answered by |
|---|---|---|
| RQ1 | Is the deviation near s = 1 a crossover rather than numerical noise (in the new population)? | exact vs double local slopes; a preregistered crossover definition (§5.4) |
| RQ2 | Does the asymptotic law emerge at weak coupling? | the exact order q from Taylor coefficients (§5.2); the double slope in Ω_asym (§5.3) |
| RQ3 | Is the scaling stable across independent matrices? | 120 held-out configurations; every meaningful class-R pair must satisfy the rules |
| RQ4 | Can ordering errors be distinguished from noise? | the oracle signal vs measured noise; SNR rules; controls (class D, M12, M9) |
| RQ5 | Can the criterion be fixed without tuning to previous slopes? | all thresholds derived from precision arguments or algebra; non-derived values are flagged (§9) |

Lattice structure:
- L(s) is the exact lattice with V(s) = sV(1) and exact inputs.
- Every prefix and full matrix of L is **affine** in s: M_L(s) = A0 + s·A1, with A1 = −blockdiag(V_j(1)) (physical units).
- The leads do not depend on s.
- R_W(s) (double Δ_fl(s) taken as exact) is **not** affine in s, because of diagonal rounding.

**Mathematical facts** (validity conditions, not hypotheses):
- **Fact 1 (rationality and analyticity).**
  - On L(s), every reference block, every intermediate of the correct algorithm and every mutant output is rational in s,
    and analytic wherever all inverted matrices are invertible.
  - The correct algorithm inverts W_1 − I = SC_1, I − W_j = −SC_j and Z = SC_N (Schur complements of the prefix matrices
    M^(j)/t and of M/t). Since det M^(j) = det M^(j−1)·det(t·SC_j), these are invertible wherever **every prefix matrix
    M^(j) (j = 1..N−1) and the full M are invertible**.
  - All mutants except M10 invert only these matrices, their negatives or their transposes.
  - M10 inverts Z − (Λ − I), which is checked only pointwise.
  - **The grid conditioning check does not prove invertibility between grid points.** Interval claims rely only on the
    certificate of §2.4.
- **Fact 2.** At s = 0 every Δ_j and Λ is diagonal. Pure order mutants coincide exactly with the correct algorithm.
- **Fact 3 (order of vanishing).** Let F_m(s) be the discrepancy of the designated block X_m (mutant minus exact reference)
  on L(s). On the common certified disk |s| ≤ r_0 (§2.4 item 3), F_m(s) = Σ_{k≥0} s^k C_k.
  - For pure order mutants, C_0 = 0 (Fact 2).
  - By the identity theorem, either F_m ≡ 0, or a finite smallest k = q with C_q ≠ 0 exists. Then
    ‖F_m(s)‖_F = s^q ‖C_q‖_F (1 + O(s)) as s → 0⁺.
  - **No upper bound on q is proven.** q is identified computationally and can be unresolved.
- **Fact 4.** In class D, F_m ≡ 0 for every pure order mutant.
- **Fact 5.** M12 is algebraically equivalent to the correct code: F_m ≡ 0.
- **Fact 6.** For M9, C_0 = 0: at s = 0 the clean device has W_j = w*, so the stale and current P coincide (on L and on R_W,
  because Δ_fl(0) = δ_lead,fl).
- **Fact 7.** For M4, M4L, M4R, M10, M13, M14 and M16, generically C_0 ≠ 0 (q = 0).

**Empirical hypotheses** (tested on held-out data, per class; all exact quantities are on R_W, §4.4):
- **Hyp1 (attribution):** at every grid point resolvable for both neighbours (SNR* ≥ 1e3), |φ_double − φ_exact| ≤ 0.01.
- **Hyp2 (asymptotic emergence):** in Ω_asym, |p̂ − q| ≤ 0.02 and |p̂ − p̂*| ≤ 0.005.
- **Hyp3 (W/Z accuracy):** every correct-code accuracy quantity is ≤ 1e-10 at every grid point (the fixed criterion).
- **Hyp4 (distinguishability):** the double results of class-D order mutants and of M12 stay within the noise tolerance.
  Structural mutants are detected wherever resolvable.

Descriptive only: the distribution of q, the crossover scale s_× and the regime model.

## 2. Configurations, generation and the conditioning filter

### 2.1 Fixed lattice

- dx = 1/8, t = 32, L_y = 1.
- E_perp,n from `stationary.transverse_basis(1.0, N_y)`.
- Clean leads.
- W/Z and production definitions as in `PREDICTIONS_coupled_blocks.md` §1 and §4.

### 2.2 Dimensions and classes

- N_y ∈ {2, 3, 4} and N ∈ {2, 3, 4, 6}. N = 1 has no order mutant on its code path. All matrices are ≤ 24 × 24.
- Classes:
  - **class R:** real symmetric. This is the documented physical scope of `wz_rgf.py`.
  - **class H:** complex Hermitian. This is an algebraic stress test only.
  - **class D:** diagonal (control).

### 2.3 Random matrix generation (exact dyadic values)

Root entropies (fixed): set XP **202611110001**; set HO **202611110002**; bootstrap **202611110003**.

Per configuration index idx (manifest order of §2.6), the generator is
`numpy.random.Generator(numpy.random.PCG64(numpy.random.SeedSequence(entropy=ROOT, spawn_key=(idx, attempt))))`.

Fixed draw order:
1. `n_open = gen.integers(1, N_y + 1)`.
2. `u = gen.random()`.
3. For j = 1..N:
   - `Jr_j = gen.integers(-8, 9, size=(N_y, N_y))`;
   - `Ji_j = gen.integers(-8, 9, size=(N_y, N_y))` (always drawn).

Matrix construction:
- B_j = (Jr_j + i·Ji_j·[class H])/8.
- V_j(1) = v0·(B_j + B_jᴴ)/2 with **v0 = 8**. Entries are multiples of 1/2 in [−8, 8].
- Class D: zero the off-diagonal entries.
- V_j(s) = s·V_j(1), exact for dyadic s.

v0 = 8 places s = 1 in a strong-coupling regime (a design choice for RQ1, not tuned).

### 2.4 Acceptance filter (deterministic; on failure, increment `attempt` and redraw everything)

1. **Energy rule.**
   - The window is E ∈ (E_perp,n_open, E_upper), with E_upper = E_perp,n_open+1, or E_perp,1 + 64·(2 − δ_min) if
     n_open = N_y.
   - E = rhe(64·(E_lo + u·(E_upper − E_lo)))/64, evaluated in **exact rational arithmetic** from the double values of
     E_lo, E_upper and u. rhe is round-half-to-even to the nearest integer.
   - Every channel must satisfy |δ_n| ≥ δ_min = 0.02 and δ_n ≤ 2 − δ_min, where δ_n = (E − E_perp,n)dx².
2. **Grid conditioning** (direct inverses only; empirical filter). At every grid point s, on both L(s) and R_W(s), each of the
   following must be ≤ 1e3: cond₂(M), cond₂(M^(j)), cond₂(SC_j) (from prefix inverses), cond₂(Z), and cond₂(Z − (Λ − I))
   (the M10 matrix).
3. **Continuous-interval invertibility certificate on [0, 1]** (validated, on L(s); method in §2.4.1).
   - Matrices: every A(s) ∈ {M_L(s), M_L^(j)(s), j = 1..N−1}, each affine, A(s) = A(s_κ) + (s − s_κ)A1.
   - At a certification point s_κ, §2.4.1 either returns validated rational bounds a_κ ≥ ‖A(s_κ)⁻¹‖₂, h_κ ≥ ‖A(s_κ)‖₂ and
     b ≥ ‖A1‖₂, or reports s_κ as **uncertified**.
   - **Neumann step** (operator 2-norm; real or complex s): if |s − s_κ|·a_κ·b < 1, then A(s) is invertible and
     ‖A(s)⁻¹‖₂ ≤ a_κ/(1 − |s − s_κ| a_κ b).
   - **Half radius:** r_κ = 1/(2a_κ b), an exact rational (a_κ > 0 because ω_κ < 1 implies Y_κ ≠ 0). On |s − s_κ| ≤ r_κ:
     ‖A(s)⁻¹‖₂ ≤ 2a_κ, ‖A(s)‖₂ ≤ h_κ + r_κ b, and cond₂(A(s)) ≤ 2a_κ(h_κ + r_κ b).
   - **b = 0** (equivalently A1 = 0, §2.4.1 step 7):
     - A(s) ≡ A(s_κ) for all s, and r_κ := +∞ is used only in the coverage comparisons, never multiplied by b.
     - One certified point certifies all of [0, 1], with bounds a_κ and h_κ.
     - The expression 1/(2a_κ b) is never evaluated.
   - **Certification points and coverage (certified points only):**
     - Initial points: s = 0, the 81 grid points (§3) and s = 1, all exact dyadics.
     - Every point is evaluated by §2.4.1. **Uncertified points are discarded. They are never assigned a radius, in
       particular not r = 0, which would wrongly assert invertibility at that point.**
     - Let s_(1) < … < s_(K) be the certified points. Coverage of [0, 1] holds iff all of:
       - s_(1) − r_(1) ≤ 0;
       - s_(K) + r_(K) ≥ 1;
       - s_(i+1) − s_(i) ≤ r_(i) + r_(i+1) for every consecutive certified pair.
       These are exact rational comparisons, with +∞ handled symbolically.
     - Every uncovered gap (between consecutive certified points, or at an end) is bisected at its dyadic midpoint. New
       uncertified midpoints are discarded likewise. The limits are a depth of 30 per gap and 2000 certification points per
       matrix.
   - **Accept only if [0, 1] is covered for every matrix A.** Then every listed matrix is invertible for **every** s ∈ [0, 1],
     with the recorded validated bounds. Record the validated continuous cond bound and the minimum r_κ.
   - **Common Taylor radius (Fact 3, §5.2).** For each matrix A, its s = 0 certificate (which must be certified; otherwise
     0 has no Taylor disk and the configuration fails the filter, even if 0 is covered by a neighbouring point; logged as
     R3.20) gives r_0,A = 1/(2a_A,0 b_A), or +∞ if b_A = 0. Define
     **r_0 := min_A r_0,A = 1/(2β)** with β := max_A a_A,0 b_A (r_0 = +∞ iff β = 0).
     - Every matrix entering any oracle is invertible on the complex disk |s| ≤ r_0, and its Neumann series converges there
       absolutely, with |s|·a_A,0 b_A ≤ 1/2.
   - **Not certified continuously:** R_W(s) (pointwise checks only) and the M10 matrix (pointwise only).
   - **Failure and cost policy:**
     - If [0, 1] is not covered for some matrix within the depth and point budget, the configuration fails the acceptance
       filter and is redrawn. This is recorded.
     - If the certificate is infeasible for the manifest as a whole, no empirical fallback is allowed. A new design revision
       is required.
4. **Prohibited selection.** No selection on commutators, mutant signals, slopes or any solver output.
5. **Redraw cap.**
   - At most **attempt = 0 … 9999** (10,000 attempts) per configuration index. This value is a design choice, not derived;
     see §9.
   - If the cap is reached for any idx, manifest generation **aborts**. No partial manifest is committed, no index is
     skipped and no sample count is reduced.
   - The experiment is then **not run**, and its outcome is **DESIGN INSUFFICIENT (pre-run: the preregistered sample cannot
     be drawn under the filter)**. Any continuation requires a new design revision.

#### 2.4.1 Validated method for the inverse-norm bound

**Exact representation of the inputs; the binary64 import boundary.**
- **Source-value import (the only place binary64 values are accepted).**
  - Existing binary64 values are converted once, by an import function, to their **exact dyadic rational value** (e.g.
    `fractions.Fraction(x)`, which is exact for finite binary64 x). Complex binary64 values are imported as two exact
    rationals. NaN and ±∞ are rejected.
  - Imported values: E, E_perp,n (from `stationary.transverse_basis`), the manifest's E and V_j(1), the R_W arrays Δ_fl, the
    production arrays `E*eye − H[j]` and σ_code, and the mpmath results that define Y_κ (§2.4.1 step 1).
- **Rational assembly (after import).**
  - All assembly, arithmetic and comparison uses exact rational arithmetic only. **No floating-point arithmetic is allowed
    after conversion.**
  - The rational assembly functions accept only `int` or `Fraction` (complex values as pairs of them) and **reject direct
    float inputs**.
- Every imported value is a dyadic rational: E, E_perp,n, V_j(1), s, dx² = 2⁻⁶ and t = 32.
- Every entry of A(s) is an exact dyadic rational, except the lead terms tλ_n (from −Σ = tλ). Here γ_n = 1 − δ_n is exact and:
  - open channel (|γ_n| < 1): λ_n = γ_n + i√ξ_n with ξ_n = 1 − γ_n²;
  - closed channel (γ_n > 1): λ_n = γ_n − √ξ_n with ξ_n = γ_n² − 1.
  - γ_n < −1 and |γ_n| = 1 are excluded by the energy rule, so ξ_n > 0 is an exact rational.
- **Enclosure.**
  - Bounds (exact `math.isqrt`): lo_n = isqrt(⌊ξ_n·4^ℓ⌋)/2^ℓ ≤ √ξ_n ≤ hi_n = (isqrt(⌈ξ_n·4^ℓ⌉) + 1)/2^ℓ.
  - Width: 2^−ℓ ≤ hi_n − lo_n ≤ 2^(1−ℓ), because isqrt(n + 1) ≤ isqrt(n) + 1.
  - The midpoint μ̃_n = (lo_n + hi_n)/2 replaces √ξ_n, giving λ̃_n, with a rigorous radius
    e_n := (hi_n − lo_n)/2 ∈ [2^−(ℓ+1), 2^−ℓ] ⊂ (0, 2⁻²⁵⁵] on the imaginary (open) or real (closed) part.
- **Lead perturbation, per matrix type** (A = A_mid + A_e; A_e is diagonal and s-independent):
  - full M, N ≥ 2: slice 1 and slice N each contribute A_e,nn = t(λ_n − λ̃_n), with |A_e,nn| ≤ t·e_n.
    ε_A := ū(2·Σ_n (t e_n)²).
  - full M, N = 1: one block receives both leads, A_e,nn = 2t(λ_n − λ̃_n). ε_A := ū(Σ_n (2t e_n)²).
  - prefix M^(j), 1 ≤ j ≤ N − 1: only slice 1 (left lead). ε_A := ū(Σ_n (t e_n)²).
  - In every case ‖A_e‖_F ≤ ε_A.
  - The lead roots are those of the lattice concerned: δ_L for L; δ_lead,fl for R_W.
  - R_orig has no enclosed lead (σ_code are doubles), so A_e = 0 and ε_A = ū(0) = 0.
- A1 is exactly dyadic and s-independent.

**Approximate inverse with an exactly computed residual.**
1. **Y_κ** is the approximate inverse of A_mid(s_κ) computed by mpmath (dps 60; dps 120 on retry).
   - Y_κ is defined as **the exact dyadic rational represented by each returned mpmath entry** (real and imaginary parts),
     imported **without any subsequent rounding**.
   - So the validated object is exactly the matrix used in every subsequent comparison (§4.7). How Y_κ was obtained, and
     its precision, are irrelevant to rigor.
2. **Ξ_κ** = I − A_mid(s_κ)·Y_κ, computed **exactly** in dyadic integer arithmetic (Python integers or `fractions.Fraction`;
   no floating point).
3. Exact rationals: f_Ξ = ‖Ξ_κ‖_F² (complex entries: re² + im²), f_Y = ‖Y_κ‖_F², f_mid = ‖A_mid(s_κ)‖_F² and
   f_1 = ‖A1‖_F².
4. **Upper square root** for rationals f ≥ 0, piecewise: ū(0) := 0; ū(f) := (isqrt(⌈f·4^ℓ⌉) + 1)/2^ℓ for f > 0.
   - For f > 0: with N_f = ⌈f·4^ℓ⌉ ≥ f·4^ℓ, isqrt(N_f) + 1 > √N_f ≥ 2^ℓ√f, so ū(f) > √f > 0.
   - So ū(f) ≥ √f for all f ≥ 0, and ū(f) = 0 iff f = 0. The branch is decided by an exact comparison with 0.
   - **Lower square root:** l̄(f) := isqrt(⌊f·4^ℓ⌋)/2^ℓ for rationals f ≥ 0 (l̄(0) = 0) satisfies 0 ≤ l̄(f) ≤ √f.
5. **ω_κ** := ū(f_Ξ) + ε_A·ū(f_Y). Since I − A(s_κ)Y_κ = Ξ_κ − A_e·Y_κ,
   ‖I − A(s_κ)Y_κ‖₂ ≤ ‖Ξ_κ‖_F + ‖A_e‖_F‖Y_κ‖_F ≤ ω_κ.
6. If ω_κ < 1 (an exact comparison):
   - A(s_κ)Y_κ = I − R′ with ‖R′‖₂ ≤ ω_κ < 1, so A(s_κ) is invertible, and A(s_κ)⁻¹ = Y_κ(I − R′)⁻¹.
   - Then ‖A(s_κ)⁻¹‖₂ ≤ ‖Y_κ‖₂/(1 − ω_κ) ≤ **a_κ := ū(f_Y)/(1 − ω_κ)**.
   - This uses the operator 2-norm (‖I‖₂ = 1). Frobenius norms are used only as upper bounds of 2-norms.
7. **h_κ** := ū(f_mid) + ε_A ≥ ‖A(s_κ)‖₂, and **b** := ū(f_1) ≥ ‖A1‖₂.
   - Since A1 is exact, f_1 = 0 iff A1 = 0, so **b = 0 iff A1 = 0**.
   - b is computed per matrix: a prefix whose slices 1..j all have V(1) = 0 takes the b = 0 branch even if the full matrix
     does not.
   - Zero-branch consequences:
     - f_Ξ = 0 gives ω_κ = ε_A·ū(f_Y), which is valid.
     - f_Y = 0 means Y_κ = 0, so Ξ_κ = I, ω_κ ≥ √n_tot ≥ 1, and the point is uncertified.
     - f_mid = 0 would give Ξ_κ = I and an uncertified point, so this branch never produces a bound (it is vacuous).
   - **On ε_A:**
     - ε_A = 0 ⟹ A_e = 0, because ‖A_e‖_F ≤ ε_A. **The converse is not claimed.** ε_A is built from the enclosure radii, not
       from the actual error.
     - In this design the converse situation cannot occur. Every enclosed lead has e_n ≥ 2^−(ℓ+1) > 0, and A_e ≠ 0 as well:
       - λ̃_n = λ_n would require √ξ_n to equal the midpoint (2k + 1)/2^(ℓ+1) for an integer k, i.e.
         ξ_n = (2k + 1)²/4^(ℓ+1), with a denominator of 2^514;
       - but ξ_n is built from E (multiples of 1/64), the doubles E_perp,n and dx² = 2⁻⁶, so its denominator is at most
         about 2^130;
       - and when √ξ_n is itself dyadic (e.g. γ_n = 0, ξ_n = 1), the midpoint misses by 2^−(ℓ+1).
     - So ε_A > 0 and A_e ≠ 0 for L and R_W, while ε_A = 0 and A_e = 0 for R_orig.
     - Every use of ε_A needs only ‖A_e‖_F ≤ ε_A.
   - No division by ū(·) occurs. The only divisions are by 1 − ω_κ > 0 and by 2a_κ b with a_κ, b > 0.
8. If ω_κ ≥ 1: retry once with Y_κ at dps 120. If it is still ≥ 1, s_κ is **uncertified** and is discarded (§2.4 item 3).

**Rounding coverage.** Steps 2–7 are exact rational arithmetic. The only approximations are the two square-root enclosures,
which are rigorous by construction. No floating-point value enters an inequality except through Y_κ, which is used as an
exact dyadic matrix. Rigor rests on:
- (i) the inequalities above;
- (ii) Python's exact integer arithmetic and `math.isqrt`;
- (iii) implementation correctness, addressed by T1–T6 and R0.

**Alternative considered.**
- Directed-rounding interval arithmetic (`mpmath.iv`, real `ivmpf` and complex `ivmpc`) would also be valid.
- It is not selected, because:
  - it depends on library directed rounding for complex operations;
  - wrapping inflates enclosures with n_tot;
  - it still needs a validated interval-matrix norm.
- It may serve as an optional cross-check of a_κ, with no decision role.


#### 2.4.2 Certificate, assembly and pointwise-certificate self-tests (T1–T6)

All six tests are **deterministic**: fixed exact inputs, exact arithmetic, no random draws. They must all pass **before the
manifest is generated**. Any failure means the certificate code is not used and the manifest is not generated. Each test is
part of the R0 review.

**T1 — validated inverse bound.**
- **Property.** For an exact matrix A with known inverse:
  - ω < 1 implies a ≥ ‖A⁻¹‖₂;
  - h ≥ ‖A‖₂ and b ≥ ‖A1‖₂;
  - ω ≥ 1 yields "uncertified" and **no** bound.
- **Inputs (exact):**
  - (a) A = diag(2, −3, 1/2);
  - (b) A = 2·Perm, with Perm the 4 × 4 cyclic permutation;
  - (c) A = [[1+i, 1/2],[0, 1−i]];
  - (d) Y deliberately poor: Y = A⁻¹/2 for case (a), so Ξ = I/2 and ‖Ξ‖_F = √3/2 < 1, giving a valid but loose bound; and
    Y = A⁻¹/2 for case (b), n_tot = 4, ‖Ξ‖_F = 1, which must be uncertified;
  - (e) Y = 0, which must be uncertified;
  - (f) the full-M and prefix lead-enclosure paths, using an L-lattice matrix at s = 0 with N_y = 2 and N ∈ {1, 2, 6} (N = 1
    tests the 2t·e_n rule), plus a 24 × 24 case (N_y = 4, N = 6).
- **Independent reference:**
  - for (a)–(e): exact closed forms of ‖A⁻¹‖₂ and ‖A‖₂. For the 2 × 2 case, the squared singular values satisfy
    σ²_max = (F + √(F² − 4|det|²))/2 with F = ‖·‖_F². The check a ≥ σ_max is then done exactly via
    (2a² − F)² ≥ F² − 4|det|² with 2a² − F ≥ 0.
  - for (f): an mpmath SVD at dps 300, accepted only as a numerical reference with margin 10^−290·σ_ref. This is the
    identity-threshold convention 10^−(dps−10).
- **Success:** every returned bound ≥ the reference (exactly, or beyond the stated margin for (f)); every case predicted
  uncertified is reported uncertified; ω, a and h are computed as exact rationals (checked by type).
- **Failure:** any bound below the reference; any certification with ω ≥ 1; any floating-point value inside an inequality.
- **Edge cases:** exact inverse (Ξ = 0, so ω = ε_A·ū(f_Y)); ε_A = 0 (no leads); N = 1.

**T2 — no certification across singularities.**
- **Property.** No certified interval or disk contains a point where the matrix is singular. For a family singular on
  [0, 1], coverage is reported failed and the configuration is rejected.
- **Inputs (exact affine families):**
  - (a) A(s) = diag(1 − 2s, 1): singular at 1/2, a dyadic point that is also a bisection midpoint;
  - (b) A(s) = diag(1 − 3s, 1): singular at 1/3, non-dyadic;
  - (c) A(s) = diag(s, 1): singular at 0;
  - (d) A(s) = diag(1 − s, 1): singular at 1;
  - (e) 24 × 24: A(s) = I − 3s·e_1e_1ᵀ, singular at 1/3;
  - (f) A(s) = [[1 − 2s, 2⁻¹⁰],[−2⁻¹⁰, 1 − 2s]]: invertible on the real line, singular at s = 1/2 ± i·2⁻¹¹;
  - (g) a constant singular matrix (b = 0, A0 = diag(0, 1)).
- **Independent reference:** the exact singular points. For each certified point check r_κ < |s_κ − s_sing| exactly. For the
  complex singular points (f), compare squares: r_κ² < (s_κ − 1/2)² + 2⁻²².
- **Success:**
  - (a)–(e) and (g): coverage fails and no certified interval contains the singular point;
  - (f): every certified disk excludes both complex singular points. Coverage of the real interval may succeed, and the
    Taylor radius r_0 at s = 0 must satisfy r_0 < |s_sing|.
- **Failure:** any certified interval or disk that contains a singular point, or coverage reported for (a)–(e) or (g).
- **Edge cases:** a singular point coinciding with an initial grid point (a), with an endpoint (c, d), or with no dyadic
  point at all (b, e).

**T3 — lead-root enclosure.**
- **Property:**
  - lo_n² ≤ ξ_n ≤ hi_n² (exact);
  - 2^−ℓ ≤ hi_n − lo_n ≤ 2^(1−ℓ);
  - correct branch: open λ is imaginary-part enclosed with |λ| = 1, checked exactly as γ_n² + lo_n² ≤ 1 ≤ γ_n² + hi_n²;
    closed λ ∈ (0, 1);
  - e_n > 0;
  - A_e ≠ 0 in every case (§2.4.1 step 7).
- **Inputs (exact rational δ):**
  - 1/50, −1/50, 99/50 (= 2 − δ_min), 1 (γ = 0, ξ = 1: exact dyadic square, open branch), 1/2, −1/2;
  - −1/4 (γ = 5/4, ξ = 9/16, √ξ = 3/4: exact dyadic square, closed branch). An open-branch ξ = 1/4 is unreachable for
    rational γ (it would need γ² = 3/4), so it is not used;
  - δ_n values for one fixed N_y = 4 configuration specified in the test file by explicit integers.
- **Independent reference:** the exact rational inequalities above (no floating point). Supplementary: an mpmath dps-300
  evaluation of λ_n lies in the enclosure.
- **Success:** all exact inequalities hold; the branch is correct; e_n > 0; λ̃_n ≠ λ_n.
- **Failure:** any violated inequality, a wrong branch, e_n = 0, or a midpoint equal to the root.
- **Edge cases:** exact dyadic squares, ξ_n close to 0 (δ = ±0.02), the band edge δ = 1.98.

**T4 — the b = 0 path.**
- **Property:**
  - ū(0) = 0 exactly, and ū(f)² > f for f > 0;
  - a matrix with A1 = 0 is certified on [0, 1] from a single point, with r = +∞;
  - neither 1/(2a_κ b) nor r_κ·b is evaluated when b = 0;
  - a prefix matrix with b = 0 is handled while the full matrix has b ≠ 0.
- **Inputs:**
  - ū at f ∈ {0, 1, 2, 1/4, 2⁻⁶⁰⁰};
  - A(s) = diag(2, 3) with A1 = 0;
  - an N = 2, N_y = 2 L-lattice configuration with V_1(1) = 0 and V_2(1) ≠ 0 (fixed integers).
- **Independent reference:** exact arithmetic. **Instrumentation:** the implementation raises an exception if 1/(2a_κ b) or
  r_κ·b is evaluated with b = 0. The test asserts that no exception occurs and that the b = 0 flag is set for the right
  matrices.
- **Success:** all of the above hold.
- **Failure:** ū(0) ≠ 0, an instrumented exception, a wrong b = 0 flag, or more than one certification point used for a
  constant matrix.
- **Edge cases:** f = 0; tiny f > 0; prefix-only b = 0.

**T5 — exact-rational assembly** (A_mid(s), A1, and the R_W and R_orig reference matrices).
- **Property (two separate stages):**
  - **(a) Import:** the import function maps every finite binary64 value to its exact dyadic rational. The round trip
    float(rational) == x holds, and the denominator is a power of two. NaN and ±∞ are rejected.
  - **(b) Rational assembly:** the certificate and validated-reference matrices are assembled from imported rationals only:
    - L: E, E_perp,n, V(1), s, t, dx² and the enclosure midpoints;
    - R_W: the imported Δ_fl;
    - R_orig: the imported production arrays `E*eye − H[j]` and σ_code.
  - **No floating-point arithmetic occurs after import.**
- **Inputs:**
  - three fixed configurations (N_y = 2, 3, 4; N = 1, 2, 6), given by explicit integer arrays Jr_j and Ji_j and explicit
    dyadic E;
  - s ∈ {0, 1, 3/4, 2⁻⁴⁰, 3·2⁻⁴²}.
- **Independent reference:** a second assembly routine, written separately from the definition M = E − H − Σ_L ⊕ Σ_R using
  `fractions.Fraction` only.
- **Import checks:**
  - binary64 inputs including 0.1, −0.0, 2⁻¹⁰⁷⁴ (subnormal), 1.7976931348623157e308, π (as binary64) and one complex128
    value: their exact rationals and round trips;
  - NaN and +∞ must be rejected.
- **Type checks:** the rational assembly functions must raise when given a `float` directly, even one that is exactly
  representable such as 0.5. Import is the only accepted path for binary64 values.
- **Exactness check:** A_mid(2⁻⁴⁰) − A_mid(0) = −2⁻⁴⁰·blockdiag(V(1)) exactly.
  - Rationale: in binary64, forming E − E_perp − s·V_aa rounds each intermediate sum to the nearest double. This perturbs
    the represented s-term by up to half an ulp of the sum, which matters when the term is small, as at the smallest grid s.
  - Exact rational assembly preserves the represented term exactly.
- **Success:** import exact with correct rejections; entrywise exact equality with the independent assembly; no float
  object in any post-import intermediate; direct float input rejected; the exactness check holds.
- **Failure:** any import inexactness, any mismatch, any float intermediate after import, or acceptance of a direct float
  input by an assembly function.
- **Edge cases:** N = 1 (both leads on one block); prefix matrices; s = 0; the smallest grid s.

**T6 — pointwise infeasibility certificate (§5.6).** Not implemented or executed. No pass is claimed.
- **Independent reference implementation.** A separate module sharing no code with the certificate module, providing:
  - (i) exact `Fraction` Gaussian elimination for A⁻¹;
  - (ii) the §4.2 mutant recursions (code-change column) executed in exact rationals;
  - (iii) exact positive-semidefiniteness tests by rational LDLᴴ, deciding ‖X‖₂ ≤ c exactly as c²I − XᴴX ⪰ 0, and the strict
    counterpart.
- **Exact test lattices** (true inverse exact, because the lead roots are rational):
  - closed channels with λ ∈ {1/2, 1/3, 2/5} (e.g. γ = 5/4 gives λ = 1/2);
  - open channels from Pythagorean pairs (e.g. γ = 3/5 gives λ = 3/5 + 4i/5). Here √ξ = 4/5 is not dyadic, so A_e ≠ 0 is
    exercised while the truth is known exactly;
  - N_y ∈ {1, 2}, N ∈ {2, 3, 4} (N = 4 for M5);
  - dyadic V(1) with nonzero off-diagonal entries.
- **T6.1 (polynomial evaluation).** For M1, M2, M3, M5, M6, M7, M8 and M15, at every N where the mutant is on path, F_m
  evaluated on the exact blocks must equal, **exactly**, the exact mutant recursion's designated block minus the exact
  reference block. This also checks the polynomial forms. Any mismatch is a failure.
- **T6.2 (reference-error propagation).** Inputs: Y = A⁻¹ + Δ with Δ = 2⁻⁴⁰·e₁e₁ᵀ (rank one, so ‖Δ‖₂ = ‖Δ‖_F exactly), and
  Y = A⁻¹. Success: η_A ≥ t‖Δ‖₂ exactly, and E_m ≥ ‖F_m(B) − F_m(B̃)‖₂ by the exact PSD test. Failure: any bound violated.
- **T6.3 (full and prefix certificates).** For every prefix j = 1..N−1 and the full matrix: ω < 1, and
  a_A = ū(f_Y)/(1 − ω) with a_A²I − A⁻ᴴA⁻¹ ⪰ 0 exactly. Failure: a violation, or a certificate returned with ω ≥ 1.
- **T6.4 (lead-enclosure propagation).**
  - On the rational-λ and Pythagorean lattices, A_e = A_true − A_mid is known exactly. Success: ε_A² ≥ ‖A_e‖_F², and
    ω²I − RᴴR ⪰ 0 for R = I − A_true·Y. Failure: any violation.
  - **Supplementary irrational-lead case (C4).** This case is checked against an mpmath dps-300 reference and is reported in
    a separate field labelled **"numerical cross-check (not validated)"**.
    - Disagreement is blocking: it must be investigated before the manifest freeze.
    - Agreement is never reported or cited as an exact validation result, and does not contribute to any "validated" claim.
    - The exact T6.4 checks are those on the rational-λ and Pythagorean lattices only.
- **T6.5 (U_m against the oracle).** On the T6.1 lattices at several dyadic s: U_m ≥ 0 and U_m² ≥ D*_true² (exact).
  Failure: any violation.
- **T6.6 (G_lo against G).** On the T6.1 full matrices: G_lo²I − ĜᴴĜ is **not** positive definite (strict exact test), hence
  G_lo ≤ ‖Ĝ‖₂. Failure: G_lo²I − ĜᴴĜ ≻ 0.
- **T6.7 (exact comparison and threshold type; C3).**
  - **Threshold.** The threshold is the module constant THRESHOLD, constructed from integers as `Fraction(1, 10**10)`.
    The test asserts:
    - `THRESHOLD == Fraction(1, 10**10)`;
    - `THRESHOLD != Fraction(1e-10)`, which guards against importing the binary64 approximation of 1e-10;
    - `type(THRESHOLD) is Fraction`.
  - **Operand types.** The comparison function accepts only `int` or `Fraction` operands. Passing a `float`,
    `numpy.floating`, or a `Fraction` built from a float in place of THRESHOLD must raise or be detected.
  - **Comparison cases:**
    - (a) D*_true far below 10⁻¹⁰·G gives "infeasible at point";
    - (b) D*_true far above gives "not infeasible";
    - (c) the synthetic exact equality U = THRESHOLD·G gives "not infeasible" (strict inequality).
  - **Failure:** any wrong classification, any float accepted, or a threshold not equal to `Fraction(1, 10**10)`.
- **T6.8 (INCONCLUSIVE paths).** Each of the following must return INCONCLUSIVE or not applicable, **never** infeasible and
  never feasible:
  - Y = 0;
  - Y = A⁻¹/2 with n_tot ≥ 4 (ω ≥ 1);
  - a missing prefix certificate;
  - a test-only coarse enclosure (ℓ = 4, a harness parameter; the production ℓ = 256 is unchanged) making ω ≥ 1;
  - a near-threshold lattice making U_m ≥ THRESHOLD·G_lo;
  - M10 or an unsupported mutant requested;
  - an evaluated set different from S_test, including an empty one.
- **Edge cases:**
  - N_y = 1;
  - class D (F_m ≡ 0 exactly from exact blocks);
  - the j = 2 terms of M5 vanish exactly;
  - M8 at N = 2 (identical by algebra) is not applicable;
  - S_test on the design grid has exactly 43 points.

### 2.5 Conditional sampling (explicit limitation)

The accepted configurations are draws **conditional on passing §2.4**.

Consequences:
- Every claim refers to that conditional population only.
- Poorly conditioned, near-singular, resonant, bound-state-adjacent and near-threshold (|δ_n| < 0.02) devices are
  systematically excluded. Nothing can be inferred about them.

Recorded and reported, never used to choose thresholds:
- the rejection log;
- the rejection rate per stratum (N_y, class, N);
- the validated minimum r_κ per configuration;
- the validated continuous cond bound per configuration.

### 2.6 Configuration sets

Index order: N_y ∈ (2, 3, 4) → class ∈ (R, H, D) → N ∈ (2, 3, 4, 6) → draw d.

| set | class R | class H | class D | total |
|---|---|---|---|---|
| XP (exploratory) | 3 × 4 × 2 = 24 | 24 | 3 × 4 × 1 = 12 | 60 |
| HO (held-out) | 3 × 4 × 4 = 48 | 48 | 3 × 4 × 2 = 24 | 120 |

## 3. Coupling grid

- **Grid:** s ∈ {2⁻ᵏ : k = 0..40} ∪ {3·2⁻⁽ᵏ⁺²⁾ : k = 0..39}. That is 81 dyadic points (about 6.6 per decade). s = 1 is a
  diagnostic point.
- **Lower edge:**
  - ν_X ≥ τ_abs ≥ 1e-13·‖Ĝ_full‖₂ (§4.5). For an order-1 signal with |C_1| ≤ 1 (Asm2, used only for sizing), SNR* < 10
    before 2⁻⁴⁰. The grid therefore extends past the detectability limit.
  - At 2⁻⁴⁰ the diagonal s-part of Δ_fl is still represented, with relative rounding ~1e-3. R_W accounts for this exactly.
- **Asymptotic order:** from Taylor coefficients (§5.2), not from grid points.
- **Crossover:** §5.4. No grid choice uses previous slopes.

## 4. Quantities, references, mutant oracle and numerical controls

### 4.1 Solvers and exact computations per (configuration, s)

- **WZ solver:** `wz_rgf.wz_endpoint_blocks` and `wz_fisher_lee` (read-only).
  - Plus the instrumented copy `wz_mutants.wz_instrumented` (unchanged from `63e046d`), which must be bit-identical at the
    four endpoint blocks.
  - The double mutants come from `wz_mutants.py`, unmodified.
- **orig solver:** production `stationary.scattering_matrix(E, dev, min_closed=0)` on a directly built `Device`. Accuracy is
  reported only; no decision weight.
- **Ref(·):** the direct mpmath inverse of the full block-tridiagonal matrix, assembled exactly (T5), with Σ_L and Σ_R added
  separately.
  - Built for L, R_W and R_orig at dps 60, with escalation to dps 120.
  - Prefix references: the same routine without Σ_R.
  - No recursion code is used.
- **Oracle D*_oracle:** the closed forms of §4.2 from the exact prefix and full inverses, on R_W (primary) and on L.
- **Transcription:** a clean-room mpmath transcription of every mutant (§4.3), on R_W and L.

### 4.2 Exact algebraic change of every mutant, and the specification-derived oracle

Notation (exact, from direct inverses of R_W or L, in Ĝ units):
- Ĝ^(j)_ab is a block of (M^(j)/t)⁻¹, and Ĝ_ab is an endpoint block of (M/t)⁻¹.
- P_j = −Ĝ^(j)_jj, and P_1 := (I − W_1)⁻¹ = −Ĝ^(1)_11.
- χ_j = Ĝ^(j)_j1, ψ_j = Ĝ^(j)_1j, θ_j = ψ_{j−1}P_jχ_{j−1}.
- Z⁻¹ = Ĝ_NN.
- **W_j = I + SC_j**, with SC_1 = A_11/t (slice 1 including the left lead) and SC_j = A_jj/t − Ĝ^(j−1)_{j−1,j−1} for j ≥ 2.
  This is algebraically identical to I + (Ĝ^(j)_jj)⁻¹, and is the form used for validated references (§4.7).
- Σ_j sums over the loop slices j = 2..N−1.

Correct code (`wz_mutants.py`, `mutant=None`, identical to `wz_rgf.wz_endpoint_blocks`):
- init: `W = 2Δ_1 + Lm`; `G11 = solve(W − I, I)`; `Gj1 = G1j = G11`.
- loop: `W = 2Δ_j + solve(I − W, W)`; `P = solve(I − W, I)`; `G11 = G11 − G1j@P@Gj1`; `Gj1 = P@Gj1`; `G1j = G1j@P`.
- closing: `Z = 2Δ_N + Lm + solve(I − W, W)`; `Zi = solve(Z, I)`; `G11 = G11 + G1j@Zi@Gj1`; `GN1 = −Zi@Gj1`;
  `G1N = −G1j@Zi`; `GNN = Zi`.

| mutant | exact code change (`wz_mutants.py`) | designated block X_m | oracle: exact mutant value of X_m |
|---|---|---|---|
| M1 | `Gj1 = Gj1@P` | Ĝ_N1 | −Z⁻¹ Ĝ^(1)_11 P_2 P_3 ⋯ P_{N−1} (indices increasing left to right) |
| M2 | `G1j = P@G1j` | Ĝ_1N | −(P_{N−1} P_{N−2} ⋯ P_2 Ĝ^(1)_11) Z⁻¹ |
| M3 | M1 and M2 | Ĝ_N1 | as M1 |
| M4 | chains updated first, then `G11 = G11 − G1j@P@Gj1` | Ĝ_11 | Ĝ_11 − Σ_j ψ_{j−1}(P_j³ − P_j)χ_{j−1} |
| M4L | `Gj1 = P@Gj1` before the G11 update; `G1j` after | Ĝ_11 | Ĝ_11 − Σ_j ψ_{j−1}(P_j² − P_j)χ_{j−1} |
| M4R | `G1j = G1j@P` before the G11 update; `Gj1` after | Ĝ_11 | as M4L |
| M5 | `G11 = G11 − Gj1@P@G1j` | Ĝ_11 | Ĝ_11 − Σ_j (χ_{j−1}P_jψ_{j−1} − θ_j). The j = 2 term is 0 |
| M6 | `GN1 = −(Gj1@Zi)` | Ĝ_N1 | −χ_{N−1} Z⁻¹ |
| M7 | `G1N = −(Zi@G1j)` | Ĝ_1N | −Z⁻¹ ψ_{N−1} |
| M8 | `G11 = G11 + Gj1@Zi@G1j` (closing) | Ĝ_11 | Ĝ_11 + χ_{N−1}Z⁻¹ψ_{N−1} − ψ_{N−1}Z⁻¹χ_{N−1}. Zero change at N = 2 |
| M9 | `P = solve(I − W_old, I)` before the W update | Ĝ_N1 | −Z⁻¹·(P_{N−2} P_{N−3} ⋯ P_2 P_1)·Ĝ^(1)_11: the ordered product Π_{k=N−2 down to 1} P_k, indices decreasing left to right. For N = 3 it is the single factor P_1, giving −Z⁻¹P_1Ĝ^(1)_11 |
| M10 | `Z = 2Δ_N + solve(I − W, W)` (Lm dropped) | Ĝ_NN | (Z − (Λ − I))⁻¹ |
| M12 | `solve((I−W).T, W.T).T` for (I − W)⁻¹W (loop and Z) | Ĝ_11 | Ĝ_11 (algebraic identity) |
| M13 | `G11 = G11 + G1j@P@Gj1` (loop) | Ĝ_11 | Ĝ_11 + 2Σ_j θ_j |
| M14 | init `G11 = solve(I − W, I)` | Ĝ_11 | Ĝ_11 − 2Ĝ^(1)_11 (and Ĝ_N1 → −Ĝ_N1, Ĝ_1N → −Ĝ_1N) |
| M15 | M6 + M7 | Ĝ_N1 | as M6 |
| M16 | closing: GN1, G1N first; `G11 = G11 + G1N@Zi@GN1` | Ĝ_11 | Ĝ^(N−1)_11 + ψ_{N−1}Z⁻³χ_{N−1} |
| M11 | N = 1 boundary `2Δ_1 + Lm` | — | off-path for every N ≥ 2 |

The oracle values follow from the code changes by exact algebra, using the identities between recursion intermediates and
prefix-reference quantities (verified at ≤ 1e-50, §4.7 item 5). The oracle uses **no recursion code**.

### 4.3 Independent mutant reference, apparatus accounting, and shared-error checks

Three implementations of each mutant output:
1. the double mutant (`wz_mutants.py`);
2. the clean-room mpmath transcription;
3. the oracle (§4.2).

**Clean-room transcription and review gate R0:**
- The transcription is written from §4.2's "exact code change" column and the correct-code listing only, without reading
  `wz_mutants.py`. Otherwise it is written by an independent implementer. Authorship is recorded.
- **R0 (mandatory before any solver run)** is an independent reviewer's sign-off in `REVIEW_coupling_scaling_v2.md`
  (identity, date, commit) covering:
  - (a) §4.2's change column against `wz_mutants.py`;
  - (b) each oracle formula, by independent derivation;
  - (c) the oracle implementation;
  - (d) the transcription;
  - (e) the certificate, validated-residual and pointwise-certificate (§5.6) code, and T1–T6.
  - Running without a committed R0 sign-off → INVALID.
- The design-level R0 review of `f101224` (CONDITIONAL PASS) verified (a) and (b) against the document. It does not replace
  the implementation-level sign-off of (c)–(e).

**Mutant-specific magnitude scale and apparatus tolerance** (heuristic, same form as the inherited τ_end):
- Mag_m(s) := Σ over the monomials of the oracle formula of X_m of |coefficient|·Π‖factor‖₂. The factors are exact
  reference blocks on R_W at s.
- τ_abs,m(s) := max(τ_abs(s), c·n_tot·ε·max(cond₂(M), K_Schur)·Mag_m(s)), with c = 10 and K_Schur as in τ_end
  (`PREDICTIONS_coupled_blocks.md` §5).
- ν_X,m(s) := max(η_X(s), τ_abs,m(s), ε·‖X_ref(s)‖_F).
- No new numerical constant is introduced. The factor 10 and c = 10 are inherited.

**Agreement checks (matrices compared directly; R3.11):**

| check | quantity | tolerance | on violation |
|---|---|---|---|
| oracle vs transcription (mp vs mp), designated block | ‖X_m^transcription − X_m^oracle‖_F | ≤ 1e-50·‖Ĝ_full‖₂ | **INVALID** (global) |
| double vs oracle, designated block, at every grid point with SNR*_m ≥ 10 | ‖X_m^double − X_m^oracle‖_F | ≤ 10·ν_X,m | **mutant-apparatus anomaly** → the pair is excluded from detection claims, and DESIGN INSUFFICIENT for its scope |
| double vs transcription, all four endpoint blocks, every grid point | max over blocks of ‖X^double − X^transcription‖_F | ≤ 10·ν_X,m | mutant-apparatus anomaly (as above) |

Rationale for the classification (R3.10):
- A double-vs-oracle excess can come from legitimate rounding amplification in a structural mutant (e.g. P³, Z⁻³) that
  ν_X,m underestimates, so it must not invalidate the whole experiment.
- It can also come from an implementation error in `wz_mutants.py`. Either way it blocks PASS (DESIGN INSUFFICIENT), so no
  unexplained disagreement can pass.
- Oracle-vs-transcription disagreement involves no double rounding and indicates an apparatus error → INVALID.

**First-deviation check:** the transcription's first deviating intermediate equals the prefix-level closed form at the slice
where the mutant first acts (e.g. M1: Ĝ^(1)_11 P_2 in place of P_2 Ĝ^(1)_11), within 1e-50·‖Ĝ_full‖₂. Violation → INVALID.

**Residual risk:** an error in §4.2's change column itself, shared by all three implementations, is detectable only by R0
item (a).

### 4.4 Fitted error quantity and lattice

- D_m(s) = ‖X_m^double(s) − X_m^ref(s)‖_F (absolute, Ĝ units).
- D*_m(s) = ‖X_m^oracle(s) − X_m^ref(s)‖_F.
- **All on R_W(s)**, the lattice the double code sees. **R_W is the lattice for D*, φ_exact, p̂*, SNR* and dev(s).**
- L is used only for the Taylor coefficients C_k, q and ‖C_q‖, and for control §4.7 item 3.
- All four blocks are recorded as secondary data.

### 4.5 Noise and SNR (unchanged rule; notation updated)

- η_X(s) = ‖X^correct,double − X^ref‖_F (empirical, R_W).
- ν_X(s) = max(η_X, τ_abs, ε·‖X_ref‖_F), with τ_abs = τ_end·‖Ĝ_full‖₂ (heuristic).
- SNR*_m(s) = D*_m(s)/ν_X(s).

| SNR*_m | meaning |
|---|---|
| ≥ 10 | detection-eligible |
| ≥ 1e3 | regression-usable |
| < 10 | **no claim; never interpreted as harmless** |

- **Detection:** D_m ≥ 10·ν_X **and** the double-vs-oracle matrix agreement of §4.3 holds.

### 4.6 Classification of every error quantity and threshold

| quantity | nature | role |
|---|---|---|
| 1e-10 accuracy criterion (endpoint, prefix, W_j, P_j, ΔS) | **fixed criterion**, unchanged | FAIL-W when reliably measured (§6) |
| judgeability gate τ_end ≤ 1e-11 | heuristic | decides only whether heuristic-bound comparisons are made. A non-judgeable case without a reliably measured > 1e-10 error → DESIGN INSUFFICIENT |
| τ_end, τ_pre, τ_W, τ_P, τ_S, U/Q/symmetry bounds, τ_abs,m | heuristic | exceedance → anomaly → DESIGN INSUFFICIENT; never FAIL |
| reference accuracy ‖A⁻¹ − Y‖₂ ≤ ‖Y‖₂·ω/(1 − ω) (§2.4.1 applied to the reference Y) | **validated** | ω/(1 − ω) > 1e-55 → escalation (dps 120), then INVALID |
| certificate bounds (§2.4) | **validated**, conditional on T1–T5 and R0 | acceptance filter |
| pointwise infeasibility certificate U_m, G_lo (§5.6) | **Designed; validation pending T1–T6 and R0.** | reported only; no decision role |
| Taylor tail majorant B_tail (§5.2 (ii)) | **validated** (from validated a, b) | series check; failure after dps-200 escalation → DESIGN INSUFFICIENT |
| B_prec precision allowance (§5.2 (ii)) | non-validated (mp rounding) | series check only; never triggers INVALID by itself |
| zero test (§5.2) | non-validated numerical classification | order identification |
| dps-100 and `lu_solve` cross-checks | empirical | INVALID on disagreement |
| η_X | empirical | noise |
| SNR levels, detection factor, apparatus factor 10 | design constants (inherited) | detection / anomaly |
| slope tolerances 0.01, 0.02, 0.005; window adequacy 0.02 | derived from Asm1 plus empirical noise | FAIL-D; window adequacy → DESIGN INSUFFICIENT |
| zero-test thresholds | design constants (not derived; §9) | order identification |
| conditioning ≤ 1e3; δ_min = 0.02; redraw cap 10,000 | design constants (the cap is not derived; §9) | filter |

### 4.7 Represented-lattice and reference controls (violation → INVALID unless stated)

1. Off-diagonal Δ_fl(s) entries are bitwise −s·V_off(1)·dx². The diagonal rounding of the s-part is recorded.
2. e_glob(Ref(R_W) − Ref(L)) and e_glob(Ref(R_orig) − Ref(L)) are ≤ 10·cond₂(M)·ε.
3. |D*_m,R_W(s) − D*_m,L(s)|/D*_m,L(s) ≤ 1e-6 wherever D*_m,L ≥ 1e-30.
   - **Rationale (corrected, R3.8).** R_W differs from L by the rounding of fl(E − E_perp) and of fl(·) − s·V_aa on the
     diagonal of Δ_fl. This perturbation is O(ε·|Δ|) in absolute size and **does not shrink with s**.
   - It is diagonal, so it commutes with the s = 0 structure. It changes the commutator-driven discrepancy only through the
     diagonal spectrum, by a relative O(ε/gap), where gap is the smallest separation of the unperturbed diagonal values.
   - The off-diagonal parts are exact (item 1).
   - The tolerance 1e-6 is far above the expected O(1e-14). The estimate is heuristic and the check is empirical.
4. **Validated references.**
   - (a) For the full and prefix systems of L and R_W at **every** grid point, the reference inverse is the exact dyadic Y of
     §2.4.1 step 1 (no further rounding; R3.4), and ω/(1 − ω) ≤ 1e-55.
     - R_orig is validated at every grid point if the runtime allows; otherwise at k ≡ 0 (mod 4) only, with the remaining
       points labelled empirical (no decision weight).
     - **Validated W_j and P_j references:**
       - P_j = −(Ĝ^(j)_jj) is a block of a validated inverse.
       - W_j = I + SC_j (§4.2) uses only the exact A_jj/t (with the enclosed lead term for j = 1) and a block of the
         validated (j−1)-prefix inverse.
       - The validated bound of each therefore follows from §2.4.1 plus the enclosure radius.
     - **Validated S reference:** S_ref = 2i√(s_n s_n′)·Ĝ_nn′ − δ_nn′ on open channels, with s_n = Im λ_n = √ξ_n enclosed
       in [lo_n, hi_n] (§2.4.1). Computed with the midpoints and the reference Y, its error bound is
       |ΔS_ref,nn′| ≤ 2·ū(hi_n hi_n′)·err_Ĝ,nn′ + 2·|Ĝ^Y_nn′|^up·(ū(hi_n hi_n′) − l̄(lo_n lo_n′)), evaluated in exact rationals.
       - err_Ĝ is the validated block error of §2.4.1 (t·‖Y‖₂·ω/(1 − ω) bounds every entry).
       - l̄ is the rigorous lower square root of §2.4.1 step 4.
       - |·|^up is an exact rational upper bound of the modulus: ū of re² + im².
     - If any of these bounds is not computed at a point, measurements against that reference at that point are **not
       reliable** (§6).
   - (b) dps-100 agreement within 1e-55 at k ≡ 0 (mod 4), s = 1, 2⁻⁴⁰ and on the duplicate set.
   - (c) `mp.lu_solve` vs `mp.inverse` at s ∈ {1, 2⁻¹², 2⁻⁴⁰}.
   - Escalation to dps 120 on failure; INVALID if it persists.
5. Prefix and closing identities hold at ≤ 1e-50·‖Ĝ_full‖₂.
6. The instrumented copy is bit-identical to `wz_endpoint_blocks` at every execution.
7. Determinism: the first 5 configurations of XP and of HO, over the full grid with `OPENBLAS_NUM_THREADS=1`, are run twice
   and must be bit-identical.
8. Manifest hash and R0 sign-off: §8.

**Identity threshold:** τ_id(s; dps) := 10^−(dps−10)·‖Ĝ_full(s)‖₂ (the convention of the original design). A computed matrix
discrepancy ≤ τ_id counts as **zero within precision**. It is used in §5.1 (ordering-zero) and §5.2 (iii).

## 5. Scaling analysis (all procedures fixed before data collection)

### 5.1 Meaningfulness classification per (configuration, mutant); exact rules only

A **meaningful-candidate pair** is a (configuration, mutant) pair with:
- a pure order mutant (M1, M2, M3, M5, M6, M7, M8, M15);
- a class R or class H configuration;
- the mutant on path (N ≥ `MIN_N[m]`);
- not identical_by_algebra.

The classification is assigned before the order classification.

| category | rule |
|---|---|
| off_path | N < `MIN_N[m]`, or m = M11 |
| identical_by_algebra | M5 at N = 3; M8 at N = 2 |
| equivalent | M12 |
| class_D_order | a pure order mutant in class D |
| ordering-meaningful | a meaningful-candidate pair with resolved finite order q ≥ 1 (§5.2) |
| ordering-zero | a meaningful-candidate pair whose C_0..C_K_max all test zero **and** whose D*_m(s) ≤ τ_id at every grid point; excluded from detection claims and counted against coverage |
| order-unresolved | a meaningful-candidate pair whose zero test is ambiguous, or whose order exceeds K_max |
| perturbation-control | M9 |
| structural | M4, M4L, M4R, M10, M13, M14, M16 |

Identical, equivalent, class_D_order and ordering-zero pairs are **never** counted as detections.

Expected meaningful-candidate pairs:
- per N: N = 2 → 3, N = 3 → 7, N = 4 → 8, N = 6 → 8, i.e. 26 per (N_y, class, draw);
- set HO: class R 312, class H 312;
- set XP: class R 156, class H 156.

### 5.2 Identification of the asymptotic order

**Taylor coefficients (on L).**
- Each prefix and full inverse is expanded by the Neumann series A(s)⁻¹ = Σ_k (−s)^k (A0⁻¹A1)^k A0⁻¹, truncated at
  K_max = 12.
- Every oracle except M10 is a **polynomial** in these blocks (Z⁻³ = Ĝ_NN³; P_1 = −Ĝ^(1)_11; W_j appears in no oracle), so
  only Cauchy products are needed.
- M10's series is not computed (its matrix is not certified; M10 is structural and outside the scaling decision).
- The coefficients are computed in floating mp arithmetic at dps 100 and dps 200, so they are approximate.
- Convergence holds on |s| ≤ r_0 = 1/(2β) (§2.4 item 3).

**Zero test (non-validated classification; thresholds unchanged).** With z_k(dps) = ‖C_k‖_F / ‖Ĝ_full(0)‖₂:
- C_k is **nonzero** if z_k(200) ≥ 1e-60 and |z_k(100) − z_k(200)| ≤ 1e-20·z_k(200);
- C_k is **zero** if z_k(100) ≤ 1e-80 and z_k(200) ≤ 1e-170;
- otherwise C_k is **unresolved**.
- q = the smallest k with C_k nonzero, provided every C_j (j < k) tests zero.
- An unresolved C_j (j ≤ q), or all C_0..C_12 zero while D*_m is not ≤ τ_id on the grid, makes the pair
  **order-unresolved** → DESIGN INSUFFICIENT for the pair. It is never INVALID because the order is large or the coefficients
  are tiny.
- (Predicted risk: §9 item 11.)

**Consistency checks** ((i) and (iii): genuine contradictions → INVALID; (ii): escalation, then DESIGN INSUFFICIENT):
- (i) **Pure order mutants and M9:** C_0 must test zero (Fact 2, Fact 6). Nonzero → INVALID. Unresolved → DESIGN INSUFFICIENT
  for the pair.
- Criterion (ii) below is a consistency **check** whose failure, after escalation, leads to DESIGN INSUFFICIENT (not
  INVALID).
- (ii) **Series vs direct evaluation (single criterion; R3.9).**
  - The check point s* = 2⁻⁶⁰ if 2⁻⁶⁰ ≤ r_0, otherwise s* = r_0·2⁻²⁰.
  - At dps 100, compute the truncated series S_K(s*) = Σ_{k ≤ K_max} C_k s*^k and the direct oracle discrepancy F_dir(s*)
    from mp inverses of L(s*).
  - **Criterion.** The comparison passes if ‖S_K(s*) − F_dir(s*)‖_F ≤ B_tail(s*) + B_prec(s*). B_tail and B_prec are
    defined below.
  - **Status of the two terms.** **B_tail is validated** (a rigorous majorant from validated bounds, evaluated in exact
    rationals). **B_prec is non-validated** (an mp rounding allowance). The comparison itself uses mp quantities that are not
    validated (S_K and F_dir), so a failed comparison is not a proof of anything.
  - **Escalation (Rev3.1-1; wording Rev3.2-3).** Exceeding the criterion never by itself triggers INVALID.
    - If the dps-100 comparison fails, perform an **independent higher-precision numerical consistency check**:
      - F_dir(s*) is recomputed from fresh mpmath inverses of L(s*) at dps 200;
      - S_K(s*) is recomputed from the dps-200 coefficients (already computed for the zero test);
      - B_prec is re-evaluated at dps = 200.
    - If the dps-200 comparison passes, the discrepancy is attributed to precision and recorded.
    - If it also fails, the pair is classified **unresolved series inconsistency → DESIGN INSUFFICIENT** for its scope, and
      it is reported for investigation. A persistent disagreement is **not** treated as proof of a mathematical
      contradiction, and is not INVALID.
    - **dps 200 is the final escalation level. No third precision level is specified** (confirmed, Rev3.2-4).
    - B_tail is unchanged by the escalation.
    - **B_tail (validated majorant).**
      - Each block series of a matrix A in Ĝ units has coefficients with ‖·‖₂ ≤ t·a_A,0·(a_A,0 b_A)^k. With the convention
        0⁰ = 1, b_A = 0 gives zero coefficients for k ≥ 1.
      - Write F_m as a sum of monomials of degree ≤ deg_m in such blocks, with T_m = the sum of the absolute values of the
        scalar coefficients (read off §4.2, including the reference block; e.g. M13 has coefficient 2 per loop term).
      - With α = max(1, max_A t·a_A,0) and x = |s*|β ≤ 1/2, the coefficients are majorized by T_m·α^deg_m·(1 − x)^(−deg_m).
      - The tail beyond K_max is bounded by B_tail = T_m·α^deg_m·C(K_max + deg_m, deg_m − 1)·x^(K_max+1)/(1 − ζx), with
        ζ = (K_max + deg_m + 1)/(K_max + 2).
      - Proof: the coefficients of (1 − x)^(−d) are C(k + d − 1, d − 1), and their consecutive-term ratio (k + d)/(k + 1)·x
        decreases in k, so for k ≥ K_max + 1 it is ≤ ζx.
      - Requires ζx < 1, checked exactly. If it fails, s* is halved until it holds.
      - If β = 0, then x = 0 and B_tail = 0.
      - Evaluated in exact rationals.
    - **B_prec (non-validated precision allowance)** := 10^−(dps−10)·Mag_m(s*) at the dps of the comparison (100, or 200 on
      escalation). Mag_m(s*) (§4.3) is the
      magnitude of the oracle monomials before cancellation, i.e. the scale on which mp rounding acts. This defines the
      reference scale of the previous "10^(−(dps−10)) relative" allowance.
  - The previous additional "relative 1e-40" requirement is removed as conflicting (R3.9).
- (iii) **Class D and M12:** every coefficient must test zero and D*_m ≤ τ_id at every grid point (Fact 4, Fact 5). Otherwise
  → INVALID.

**Finite-scale slopes** are never used to identify q and are never INVALID.

### 5.3 Weak-coupling window and regression (unchanged rules)

**Window:** Ω_asym = {grid s ≤ 2⁻¹²} ∩ {SNR*_m(s) ≥ 1e3} (per pair, on R_W).
- Edge rationale: the local slope deviates from q by ≈ μs/(1 + μs). Under **Asm1: μ ≤ 25** (heuristic), a model error ≤ 0.01
  requires s ≤ 4e-4, so the edge is 2⁻¹².
- μ_est = ‖C_{q+1}‖/‖C_q‖ (from L) is reported. Since |μ_eff| ≤ μ_est (Cauchy–Schwarz), it is an upper estimate. It is not
  used to move the window.

**Regression:** unweighted OLS of ln D_m (double) and ln D*_m (exact, R_W) on ln s over all Ω_asym points.
- **Validity:** ≥ 8 points; ln-span ≥ ln 100; no removal except by the SNR rule.
- **Window adequacy** (exact R_W data; otherwise DESIGN INSUFFICIENT for the pair): |p̂* − q| ≤ 0.02, and every exact centred
  φ in Ω_asym is within 0.02 of q.
- **Hyp2:** |p̂ − q| ≤ 0.02 and |p̂ − p̂*| ≤ 0.005.
- **Uncertainty:** the deterministic interval p̂ ± Σ_i |ϑ_i|·|ln(D_i/D*_i)| (ϑ_i the OLS weights); the OLS SE is descriptive.
  Across configurations: median, min and max, plus a seeded bootstrap (10⁴ resamples) of the median.
- **Outliers:** none are removed.

### 5.4 Local slopes, crossover and regimes (RQ1; exact quantities on R_W)

- φ(s_k) = ln(D(s_{k+1})/D(s_{k−1})) / ln(s_{k+1}/s_{k−1}), exact and double. The one-sided φ(1) uses s = 1 and 0.5.
- **Hyp1:** |φ_double − φ_exact| ≤ 0.01 wherever both neighbours have SNR*_m ≥ 1e3, including s ∈ [2⁻⁴, 1].
- **Crossover:**
  - dev(s) = D*_m,R_W(s)/(‖C_q‖_L s^q) − 1;
  - s_× = the largest grid s with |dev(s′)| ≤ 0.1 for all grid s′ ≤ s;
  - secondary: the largest s with |φ_exact − q| ≤ 0.05 below it.
- **Regime models** (descriptive): ln D* = a + q ln s + ln|1 + Σ_{i=1}^{mreg} b_i s^i| for mreg = 0, 1, 2. Select the smallest
  mreg with max |log residual| ≤ 0.01; otherwise "strongly nonlinear in range". No decision depends on this.
- **RQ1 answer rule:** if Hyp1 holds, any departure of φ from q near s = 1 is algebraic (a crossover). Otherwise it is at least
  partly numerical. This applies to the new population only (see the header on the previous slopes).

### 5.5 Controls

- **class_D_order and M12:** the double D_m ≤ 10·ν_X at every s (Hyp4). (The matrix agreement of §4.3 applies in addition.)
- **M9:** q reported (expected ≥ 1 in all classes including D).
- **Structural mutants:** detected at every s with SNR*_m ≥ 10 (Hyp4).
- **Noise-trend control (window defined, R3.14):**
  - For each configuration and each of the four endpoint blocks X of the **correct** code: OLS slope of
    ln max(η_X(s), ε‖X_ref(s)‖_F) on ln s over **all grid points s ≤ 2⁻¹²** (57 points; no SNR filter, because the SNR
    filter is mutant-specific and the correct-code noise is mutant-independent).
  - It must lie in [−0.25, 0.25]; otherwise → DESIGN INSUFFICIENT.
  - The s ≤ 2⁻¹² edge is the mutant-independent part of Ω_asym.

### 5.6 Feasibility of regression validity for q ≥ 2 pairs (analysis only; no rule, grid, threshold or window change)

Notation: G(s) := ‖Ĝ_full(s)‖₂ (on R_W).

**Requirements that a pair must meet** (all existing rules):
- Ω_asym ⊆ {grid s ≤ 2⁻¹²}, with SNR*_m ≥ 1e3 at every point used;
- at least 8 points;
- ln-span ≥ ln 100, i.e. max(Ω_asym)/min(Ω_asym) ≥ 100;
- window adequacy (exact local slopes within 0.02 of q; |p̂* − q| ≤ 0.02).

**Exact consequences of the definitions:**
- **Span ⟹ a small resolvable point.** max(Ω_asym) ≤ 2⁻¹², so the span requires some point s_min ≤ 2⁻¹²/100 ≈ 2.44e-6. The
  largest grid value ≤ 2.44e-6 is 2⁻¹⁹ ≈ 1.91e-6 (3·2⁻²⁰ ≈ 2.86e-6 is too large). So SNR*_m ≥ 1e3 must hold at some grid
  s ≤ 2⁻¹⁹.
- **Noise floor.** ν_X(s) ≥ τ_abs(s) = τ_end·G(s) ≥ 1e-13·G(s), because τ_end ≥ 1e-13 by definition.
- **Necessary condition (exact).** D*_m(s) ≥ 1e3·ν_X(s) ≥ 1e-10·G(s) at some grid s ≤ 2⁻¹⁹.
- **Point count.** The ≥ 8 requirement concerns the **SNR-qualified points actually included in Ω_asym**, not every grid
  point.
  - The interval [2⁻¹⁹, 2⁻¹²] contains 15 grid points (8 powers of two and 7 points of the form 3·2⁻ᵏ), and the full
    admissible range s ≤ 2⁻¹² contains 57.
  - So 8 qualified points with the required span can exist without every grid point in an interval being qualified.

**Leading-order translation** (approximate; valid only where the window is adequate, i.e. D*_m(s) ≈ ‖C_q‖ s^q).
Evaluating at the most favourable admissible point s = 2⁻¹⁹:

    ‖C_q‖ ≳ 1e-10 · G · 2^(19q) · (ν_X/(1e-13·G))

| q | factor 1e-10·2^(19q) | resulting requirement |
|---|---|---|
| 1 | ≈ 5.2e-5 | ‖C_1‖ ≳ 5.2e-5·G |
| 2 | ≈ 27.5 | ‖C_2‖ ≳ 27.5·G |
| 3 | ≈ 1.4e7 | ‖C_3‖ ≳ 1.4e7·G |

Each further order multiplies the requirement by 2¹⁹ ≈ 5.2e5. The factor ν_X/(1e-13·G) ≥ 1 enlarges the requirement
whenever the actual noise level exceeds its floor.

**Window adequacy imposes a second, opposing condition.** The local-slope deviation is ≈ μ_eff·s/(1 + μ_eff·s), so adequacy
at s = 2⁻¹² requires approximately |μ_eff| ≲ 0.0204·2¹² ≈ 84, where μ_eff = Re⟨C_q, C_{q+1}⟩/‖C_q‖². A large ‖C_q‖ is
therefore admissible only if ‖C_{q+1}‖ does not grow correspondingly (|μ_eff| ≤ μ_est = ‖C_{q+1}‖/‖C_q‖).

**Classification.**
- **Pointwise infeasibility certificate on R_W (pure order mutants M1, M2, M3, M5, M6, M7, M8, M15; reported only, no
  decision role).** The bounds are rigorous by construction. **Designed; validation pending T1–T6 and R0.**
  - **Polynomial forms.** For each supported mutant, F_m is the §4.2 oracle minus the designated reference block. Blocks common
    to both cancel symbolically. Coefficient magnitudes are 1:
    - M1, M3: −Ĝ_NN·Ĝ^(1)_11·P_2⋯P_{N−1} − Ĝ_N1;
    - M2: −P_{N−1}⋯P_2·Ĝ^(1)_11·Ĝ_NN − Ĝ_1N;
    - M5: −Σ_{j=2}^{N−1}(χ_{j−1}P_jψ_{j−1} − ψ_{j−1}P_jχ_{j−1});
    - M6, M15: −χ_{N−1}Ĝ_NN − Ĝ_N1;
    - M7: −Ĝ_NNψ_{N−1} − Ĝ_1N;
    - M8: χ_{N−1}Ĝ_NNψ_{N−1} − ψ_{N−1}Ĝ_NNχ_{N−1}.
  - **Certificate quantities at each grid s** (all exact rationals). For every R_W matrix A (the full matrix and each prefix
    used), with reference Y_A and ω_A < 1 (§4.7 item 4(a); the R_W lead enclosure is included via ε_A):
    - η_A := t·ū(f_Y,A)·ω_A/(1 − ω_A);
    - for each block B̃_i of t·Y_A, β̃_i := ū(‖B̃_i‖_F²), and η_i := η_A for the matrix containing the block;
    - E_m := Σ_μ |c_μ|·[Π_{i∈μ}(β̃_i + η_i) − Π_{i∈μ} β̃_i];
    - **U_m(s) := ū(‖F_m(B̃)‖_F²) + ū(N_y)·E_m**. Then U_m(s) ≥ D*_true(s) := ‖F_m(B)‖_F, with B the exact blocks of the R_W
      lattice;
    - **G_lo(s) := max(t/h_W, t·l̄(max_ab |Y_full,ab|²) − η_full)**, with h_W := ū(f_mid,W) + ε_A,W of the full R_W matrix and
      l̄ from §2.4.1 step 4. Then 0 < G_lo(s) ≤ G_true(s) := ‖Ĝ_full(s)‖₂ of the exact R_W lattice.
  - **Conclusion rule (C1).** Let S_test be the set of **preregistered design-grid points (§3) with s ≤ 2⁻¹⁹**. Under the
    fixed grid this is exactly 43 points: 2⁻ᵏ for k = 19..40 and 3·2⁻ᵏ for k = 21..41.
    - The implementation verifies that the evaluated set equals this preregistered set exactly (dyadic comparison). Any
      deviation, including an empty set, gives INCONCLUSIVE.
    - If U_m(s) < THRESHOLD·G_lo(s) at **every** s ∈ S_test, the pair is **infeasible at the exact level on the design
      grid**. THRESHOLD is the exact rational `Fraction(1, 10**10)` = 10⁻¹⁰, constructed from integers, and the comparison
      is an exact rational comparison. Otherwise the pair is INCONCLUSIVE.
    - INCONCLUSIVE includes: any uncertified or missing reference; ω ≥ 1; a bound that is too large; an unsupported mutant
      (including M10). INCONCLUSIVE is never read as feasibility.
    - No statement is made about non-grid values of s. None is needed, because Ω_asym consists of design-grid points only
      (§5.3).
  - **What this proves (exact level; C5).** Define, for the exact R_W lattice, the idealised exact-level ratio
    SNR_ideal(s) := D*_true(s)/(10⁻¹³·G_true(s)). This is the exact oracle discrepancy divided by the definitional minimum
    of the noise level, evaluated with exact quantities.
    - The certificate proves SNR_ideal(s) < 10³ at every s ∈ S_test, i.e. at every preregistered design-grid point with
      s ≤ 2⁻¹⁹.
    - Equivalently: even with exact evaluation and the noise level at its definitional minimum, the exact profile has no
      design-grid point at s ≤ 2⁻¹⁹ that could anchor the span required by §5.3.
    - The relevant constant is 10⁻¹⁰ = 10³ × 10⁻¹³, exactly.
    - This statement does not involve, and makes no statement about, the recorded ν_X, η_X, G or SNR*_m.
  - **What this does not prove.**
    - It does not validate any recorded mpmath value (D*_m, ν_X, η_X, SNR*_m, G).
    - It does not show that the recorded rule-level outcome is free of numerical evaluation error.
    - It says nothing about the double-precision solvers or about L.
    - In particular, a recorded SNR*_m ≥ 10³ at some s ∈ S_test does not contradict the certificate mathematically. It
      indicates that the recorded value differs from the idealised exact-level quantity, and it is reported as a diagnostic
      note.
  - **Rule-level outcome.** Regression validity in the run is determined solely by §5.1–§5.3 from the recorded quantities,
    and §6 uses only those. This diagnostic has no decision role: giving it one would change the preregistered rules.
    A disagreement between this exact-level result and the recorded outcome is reported as a diagnostic note only.
  - **Inputs and eligibility (C2).**
    - The certificate inequality and its implication use **no** Taylor coefficient, no r_0, no zero-test output and no
      L-to-R_W comparison. They hold for any supported pure order mutant pair, whatever its order.
    - **Eligibility** for reporting under §5.6 is restricted to pairs classified ordering-meaningful with q ≥ 2 by §5.2,
      whose zero test is non-validated. A misclassified q can change which pairs are reported. It cannot affect the validity
      of an issued certificate.
  - **Dependencies.** The certificate is conditional on:
    - the correctness of the §4.2 oracle formulas and of the polynomial forms above (R0 (b); checked exactly by T6.1 on the
      T6 test lattices only);
    - implementation correctness (T1–T6, R0);
    - §4.7 item 4(a) holding at every s ∈ S_test.
- **L-lattice estimate only (supplementary; no claim about R_W; no decision role).**
  - The validated L-lattice majorant of §5.2 (ii), D*_m,L(s) ≤ B_q(s) := T_m α^deg_m Σ_{k≥q} C(k + deg_m − 1, deg_m − 1)(βs)^k,
    holds for s ≤ r_0. It is conditional on the zero-test classification C_0..C_{q−1} = 0, which is not validated.
  - With h the validated bound on ‖M_L(s)‖₂, it may be reported whether B_q(s) < 10⁻¹⁰·t/h at every s ∈ S_test.
    - This is reported only if r_0 ≥ 2⁻¹⁹ and the evaluated set equals S_test; otherwise INCONCLUSIVE.
    - It is a statement about the L lattice only. **No conclusion about R_W, or about regression validity, is drawn from
      it.** The L-to-R_W transfer is not available: §4.7 item 3 is empirical and silent where D*_m,L < 1e-30, and t/h
      bounds G on L, not on R_W.
- **Likelihood not assessed.** The design assumes and predicts no coefficient distribution, so no likelihood estimate is
  made for any q. In particular, no case is classified as "possible but unlikely".
- **Possible under explicit coefficient conditions.** Regression validity requires, at leading order, both:
  - (C-q1) there is a set Ω of grid points with s ≤ 2⁻¹² such that ‖C_q‖·s^q ≳ 1e3·ν_X(s) holds **at every point of Ω**
    (the leading-order form of SNR*_m ≥ 1e3, which is exactly the inclusion condition of Ω_asym), with |Ω| ≥ 8 and
    max Ω / min Ω ≥ 100.
    - Points of the interval that are not in Ω are simply not included. They need not be qualified.
    - Necessarily min Ω ≤ 2⁻¹⁹, so at least one SNR-qualified point lies at s ≤ 2⁻¹⁹ (the necessary condition above).
    - For q = 2 this gives ‖C_2‖ ≳ 27.5·G at the noise floor, evaluated at that point;
  - (C-q2) |μ_eff| ≲ 84, together with exact window adequacy.
  Beyond leading order, these conditions are checked directly by the existing rules on the exact R_W profile.

**Consequence for the predicted risk (§9 item 12, retained).** A q ≥ 2 pair that does not meet (C-q1) and (C-q2) is not
regression-valid. Under the existing rule that every ordering-meaningful pair must be regression-valid, such a pair yields
DESIGN INSUFFICIENT. The risk is unchanged. The analysis makes it quantitative and gives a per-pair infeasibility test, with
no decision role.

## 6. Preregistered decision rules

Decisions use **set HO only**. Set XP is reported and can never yield PASS.

Two outcomes:
- **Physical-scope decision:** class R plus the class-D controls.
- **Algebraic stress-test outcome:** class H.

The class-H outcome never enters the physical-scope decision and must not be cited as validation within the documented scope
of `wz_rgf.py`. INVALID applies globally.

Precedence: **INVALID > FAIL-W > FAIL-D > DESIGN INSUFFICIENT > PASS**.

Definitions:
- **Judgeable case:** a (configuration, s) point with τ_end ≤ 1e-11. Judgeability governs only the heuristic-bound
  comparisons (anomalies).
- **Reliably measured** (for a correct-code quantity at a point): the validated reference bound of §4.7 item 4 (a) is
  available for that quantity at that point. The quantity is reliably measured to exceed 1e-10 when
  **(measured error) − (validated reference error bound) > 1e-10**.

**INVALID** (global): any of
- a §4.7 control failure;
- an oracle-vs-transcription or first-deviation failure (§4.3);
- a missing R0 sign-off;
- a manifest hash mismatch;
- a §5.2 consistency contradiction (C_0 nonzero where a theorem requires zero; class-D or M12 coefficients or grid
  discrepancies not zero within τ_id). The series-vs-direct check (ii) is **not** an INVALID condition (Rev3.1-1);
- a determinism mismatch.

**FAIL-W** (W/Z accuracy; per scope): any correct-code WZ-solver endpoint, prefix, W_j or P_j error, or ΔS, **reliably
measured** to exceed **1e-10**, **regardless of judgeability** (R3.12). For the physical scope this covers class R and
class D. This is the only W/Z failure criterion; the value 1e-10 is unchanged.

**FAIL-D** (a diagnostic hypothesis falsified; **not** evidence of a W/Z defect): any of
- Hyp1 violated for an ordering-meaningful pair;
- Hyp2 violated for a regression-valid, window-adequate ordering-meaningful pair;
- Hyp4 violated (a class-D order mutant or M12 above 10·ν_X; a structural mutant undetected where SNR*_m ≥ 10).

**DESIGN INSUFFICIENT** (per scope, unless marked global): any of
- an anomaly (a heuristic bound exceeded in a judgeable case while ≤ 1e-10);
- a non-judgeable case (without a reliably measured > 1e-10 error);
- **a correct-code error measured above 1e-10 that is not reliably measured** (validated reference bound unavailable) (R3.12);
- **a mutant-apparatus anomaly** (§4.3; R3.10);
- an order-unresolved meaningful-candidate pair;
- an ordering-meaningful pair that is not regression-valid or not window-adequate;
- a failed noise-trend control;
- coverage below requirement: ordering-meaningful pairs fewer than 90 % of the expected count (class R: 312; class H: 312),
  or fewer than 10 per pure order mutant;
- an unresolved C_0 (§5.2 (i));
- an unresolved series inconsistency (§5.2 (ii) after the dps-200 escalation; Rev3.1-1);
- **(global, pre-run)** manifest generation aborted at the redraw cap (§2.4 item 5).

**PASS** (per scope): none of the above.

Physical-scope PASS claim:

> For 48 independently generated held-out class-R devices (N_y ≤ 4, N ≤ 6, equal lead widths, clean leads), conditional on
> the §2.4 acceptance filter, the correct W/Z recursion meets the 1e-10 criterion and its heuristic error estimates at
> 81 coupling values, and every meaningful ordering mutant's discrepancy follows its exactly identified asymptotic order
> down to the resolvable noise floor, with class-D controls and the M12 negative control behaving as required.

Class-H PASS claim: the same statement, as an **algebraic stress test** only.

Interpretation:
- SNR*_m < 10 carries no information about harm.
- Class-D, M12, identical and ordering-zero non-detections are not detections.
- FAIL-D says nothing about W/Z correctness.
- No outcome says anything about the previous experiment's slopes.

## 7. Scope

**Can establish** (conditional on §2.4, in the new population):
- whether a strong-coupling crossover explains slope departures near s = 1;
- whether the exactly identified asymptotic order emerges in double precision;
- whether ordering errors are separable from rounding noise down to a quantified floor;
- whether the correct W/Z code meets the 1e-10 criterion across a dense coupling range.

Class R results refer to the documented physical scope. Class H results are algebraic stress tests. Class D results are
controls.

**Cannot establish:**
- general W/Z validation;
- long-device stability (N ≤ 6);
- near-threshold behaviour;
- poorly conditioned, resonant or near-singular devices;
- unequal lead widths;
- disordered systems;
- anything about B1 configuration F;
- any change to, or explanation of, the previous experiment's DESIGN INSUFFICIENT decision or its 0.78–0.79 slopes;
- validation within the physical scope from class-H data.

## 8. Reproducibility and protocol

### 8.1 Planned files (none exist yet)

All under `experiments/wz_rgf_diagnostic/`:
- `coupling_scaling_v2_manifest.py`;
- `coupling_scaling_v2_certify.py` (§2.4.1, the §5.6 pointwise certificate, and T1–T6; T6's independent reference
  implementation lives in a separate module);
- `coupling_scaling_v2.py` (driver);
- `wz_mutants_oracle.py`;
- `wz_mutants_mp.py` (clean-room transcription);
- `REVIEW_coupling_scaling_v2.md` (R0 sign-off);
- directory `results_coupling_scaling_v2/`:
  - `manifest_XP.json`, `manifest_HO.json`: per configuration, the set, idx, attempt, N_y, class, N, E, n_open, the
    integer arrays Jr_j and Ji_j, V_j(1) as exact rationals, the rejection log, the certificate (points, radii, cond bound,
    r_0), and the versions;
  - `manifest.sha256`;
  - `results_XP.json`, `results_HO.json`;
  - `points_XP.csv`, `points_HO.csv`: per (configuration, s), the metrics of both solvers, conditioning, commutator norms,
    the reference checks, and per mutant D_m, D*_m (oracle and transcription, R_W and L) for all four blocks, η_X, ν_X,
    ν_X,m, SNR*_m, the agreement results and the category;
  - `pairs_XP.csv`, `pairs_HO.csv`: per pair, category, C_0..C_12 norms at both precisions, the zero-test results, q, μ_est,
    p̂, p̂*, the interval, OLS SE, points, span, adequacy, s_× (both), the regime model, φ(1) exact and double, and the
    Hyp1/Hyp2/Hyp4 flags;
  - `coupling_scaling_v2_report.md` (separate class R, class H and class D sections);
  - `previous_files.sha256`.

### 8.2 Phases

1. This design is committed.
2. **Implementation and manifest freeze:**
   - Write the code.
   - **T1–T6 must pass before the manifest is generated.**
   - Generate the manifests (filter and certificate only; redraw cap §2.4 item 5) and commit them with `manifest.sha256`
     before any solver or reference run.
   - Complete R0 (items (a)–(e)) and commit its sign-off.
3. **Exploratory run (set XP):**
   - Publish the report.
   - At most one amendment may follow, limited to bug fixes, feasibility and added diagnostics. No threshold, window, grid,
     SNR level, tolerance, coverage, K_max, cap or decision rule may change; any such change requires a new experiment ID.
   - An amendment touching the oracle, the transcription or the certificate requires a new R0 sign-off.
4. **Confirmatory run (set HO)**, then both decisions by §6.

### 8.3 Preserving the previous experiment

- Before any run, record the sha256 of every existing file under `experiments/wz_rgf_diagnostic/` (excluding
  `__pycache__`) in `previous_files.sha256`. Verify it after every phase; any change → INVALID.
- Code imports `wz_rgf.py` and `wz_mutants.py` read-only and writes only to `results_coupling_scaling_v2/`.

### 8.4 Independent reproduction

1. Fresh clone.
2. Set `OPENBLAS_NUM_THREADS=1`.
3. Verify `manifest.sha256` (the frozen manifest is authoritative).
4. Run the driver.
5. Bitwise identity on the same build. On another build, values may differ at rounding level and the decisions must be
   unchanged.
6. Record the versions.

Cost: uncertain; to be measured in the exploratory phase. Escalations to dps 120 are expected to be common (§9 item 13).
Parallelization is allowed only over configurations, with fixed output ordering.

## 9. Assumptions, non-derived values and predicted risks

Assumptions:
- **Asm1 (heuristic):** μ ≤ 25 for the window edge. Checked by window adequacy.
- **Asm2 (heuristic):** |C_1| ≤ 1, used only for grid sizing.
- **Asm3 (heuristic):** the inherited τ bounds and their mutant-specific form τ_abs,m. Exceedance is an anomaly, never FAIL.

Values that are **not derived** from the protocol (flagged as required; not silently justified):
1. **Zero-test thresholds** (1e-60, 1e-80, 1e-170, 1e-20). Inherited design constants.
2. **K_max = 12.** Inherited search depth.
3. **Redraw cap 10,000** (R3.15). Introduced in Revision 3 because a finite cap was required. No derivation exists.
   **Confirmed by the reviewer (Rev3.2-4).** Its only effect is to abort, never to change the sample.
4. **Bisection depth 30 and 2000 points per matrix.** Inherited feasibility constants.
5. **The use of Mag_m(s*) as the reference scale of B_prec, and of Mag_m in τ_abs,m** (R3.9, R3.10). Definitional choices
   in the form of the inherited heuristic, with no new numerical constant. **Accepted by the reviewer (Rev3.2-4).** Both
   remain non-validated.

Remaining ambiguities and limitations:
6. **The certificate is on L(s) only.** R_W statements are pointwise. L uses the double E and E_perp,n as exact values.
7. **The certificate is validated conditional on implementation correctness** (T1–T6, R0). A bug missed by both would
   invalidate it.
8. **Non-validated steps that remain:** the zero test, B_prec, the dps-100 and `lu_solve` cross-checks, the grid conditioning
   filter, and all heuristic bounds.
9. **Validated statements available:**
   - nonzero coefficients could be proven by a validated enclosure (optional, no decision role);
   - only C_0 = 0 for order mutants and M9 is a theorem;
   - vanishing of C_1 … C_{q−1} for q ≥ 2 rests on the zero test.
10. **M10's matrix** is checked pointwise only.
11. **Predicted risk — zero-test scaling (R3.17).**
    - z_k is normalized by ‖Ĝ_full(0)‖₂ only. The coefficient scale grows like β^k, and within the filter β^12 can reach
      about 1e36.
    - The rounding residue of a truly zero coefficient at dps 100 (about 1e-100·β^k) can then exceed 1e-80, classifying it
      "unresolved". The pair becomes order-unresolved → DESIGN INSUFFICIENT.
    - This affects only pairs with q ≥ 2, where C_1, …, C_{q−1} must test zero; C_0 and C_1 are at safe scales.
    - Not corrected here, because the remedy (normalizing by the validated majorant) changes a threshold. It requires a new
      decision.
12. **Predicted risk — q ≥ 2 pairs (R3.17; quantified in §5.6, Rev3.1-6).**
    - Ω_asym is fixed at s ≤ 2⁻¹² with SNR*_m ≥ 1e3, which needs D* ≳ 1e-10·‖Ĝ‖.
    - For q = 2, D* ≈ |C_2|s², so two decades of span require |C_2| on the order of 10 or more. Such pairs are expected to be
      **not regression-valid**.
    - Since every ordering-meaningful pair must be regression-valid, **a single q ≥ 2 pair is predicted to produce DESIGN
      INSUFFICIENT**. This prediction is retained deliberately; it is unfavourable.
13. **Escalation frequency.** At dps 60, ‖Ξ‖_F can be close to the 1e-55 limit for n_tot = 24 and cond ≤ 1e3. dps-120
    escalation will be frequent (cost, not validity).
14. **Strictness choices** (preregistered, not derived): coverage ≥ 90 %; 100 % of ordering-meaningful pairs regression-valid
    and window-adequate; anomaly → DESIGN INSUFFICIENT; mutant-apparatus anomaly → DESIGN INSUFFICIENT.
15. **v0 = 8** positions the crossover. A different strength moves s_× but not q.
16. **Determinism** is required only within the same build.
17. **The previous slopes.** v2 cannot attribute the previous 0.78–0.79 slopes. The previous configurations are excluded and
    the previous fitted quantity differs. This limitation is deliberate (scope hygiene, R1.10). Re-evaluating the previous
    configurations as a descriptive set would require an explicit decision.
