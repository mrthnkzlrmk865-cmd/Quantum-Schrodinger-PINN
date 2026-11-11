# Design and predictions: W/Z coupling-scaling experiment v2 (independent; preregistration document)

Revision 2: focused mathematical audit of the continuous-interval certificate and of every validated-bound claim, on top of
Revision 1 (`b62ed42`), itself an audit of `01bf018`. Design only.
- Nothing in this document has been implemented or run.
- No code exists for it, no matrix of it has been generated, and no previous file is modified.

Relation to the previous experiment (`PREDICTIONS_coupled_blocks.md`, design `7d69260`, implementation `63e046d`, results
`e3c67a6`):
- Its decision, **DESIGN INSUFFICIENT**, is final and is not re-scored, re-fitted or reinterpreted here.
- Its post-hoc 3-point slopes (0.985–0.989) are **not** evidence for anything in this design and are not used to choose any
  threshold.
- Its matrices are not part of this experiment.

Standing status:
- B1 configuration F still FAILS the unchanged 1e-10 criterion.
- W/Z is experimental and NOT validated.
- No production change is authorized.

## R. Revision log (Revision 1 vs `01bf018`)

| # | change | audit item |
|---|---|---|
| D1 | Removed the restriction q ∈ {1, 2, 3}. q is any finite integer ≥ 1 (identity theorem). It is identified from Taylor coefficients computed by truncated Neumann-series algebra, with a precision-scaling zero test, up to K_max = 12 | 1 |
| D2 | A non-integer finite-scale slope estimate is no longer INVALID. INVALID now requires a contradiction between two exact computations or with a theorem. Unresolved order (ambiguous zero test, or order > K_max) is DESIGN INSUFFICIENT | 1 |
| D3 | The grid conditioning check no longer stands in for analyticity. Added a continuous-interval certificate (a Neumann-radius covering of [0, 1]) for all matrices that are inverted. Its inverse-norm bound was not yet validated; corrected in Revision 2, C1. It is part of the acceptance filter. Analytic claims are restricted to the certified interval; R_W statements are pointwise | 2 |
| D4 | Error quantities are classified as rigorous, empirical, heuristic or fixed criterion (§4.6). Exceeding a heuristic τ bound alone is an **anomaly** (it blocks PASS and yields DESIGN INSUFFICIENT), never FAIL. The fixed 1e-10 criterion is unchanged and remains the only W/Z-accuracy FAIL criterion | 3 |
| D5 | FAIL is split into FAIL-W (W/Z accuracy, the 1e-10 criterion) and FAIL-D (a diagnostic/scaling hypothesis falsified; explicitly **not** evidence of a W/Z defect) | 3 |
| D6 | Added a specification-derived **oracle**: a closed-form mutant output for every mutant, built only from exact prefix and full inverses (§4.2). It is compared three ways with the mp transcription and the double mutants. Added a clean-room transcription and a mandatory independent review gate R0 (§4.3) | 4 |
| D7 | Every mutant's exact algebraic change is documented (§4.2 table) | 4 |
| D8 | Claims are separated by class. The physical-scope decision uses class R plus the class D controls. Class H receives its own "algebraic stress-test" outcome, which never enters the physical-scope decision (§6, §7) | 5 |
| D9 | Sampling is stated as conditional on the preregistered filter (energy margins, grid conditioning, continuous certificate). Rejection statistics are recorded, and generalization limits are stated (§2.5, §7) | 6 |
| D10 | Legacy diagnostic set G **removed**, to avoid any appearance of revisiting the previous experiment | scope hygiene |
| D11 | The exact asymptote points on a fixed grid (2⁻⁹⁰…2⁻¹⁰⁰) are replaced by Taylor coefficients (D1). Finite-s exact profiles are kept only as consistency checks | 1 |
| D12 | Coverage requirements are counted per class (R: 312 expected held-out ordering pairs) | 5 |

Revision 2 changes (vs `b62ed42`). Mutant definitions, sample counts, grid, thresholds, fitting window and decision rules are
**unchanged**.

| # | change | reason |
|---|---|---|
| C1 | The certificate's inverse-norm bound is replaced. The earlier version used an mp approximation inflated by (1 + 1e-30), which is not a proof. The new method (§2.4.1) is a validated **approximate inverse with an exactly computed residual**: all inputs are dyadic, and the only irrational inputs (lead roots λ) are enclosed by exact integer square-root bounds | the earlier bound was not validated |
| C2 | Every bound is stated for the spectral norm ‖·‖₂ and computed through the rigorous inequality ‖·‖₂ ≤ ‖·‖_F. Frobenius norms are evaluated exactly, then given a rigorous rational upper square root | norm consistency; ‖I‖_F ≠ 1 invalidates a Frobenius-only Neumann bound |
| C3 | Special case b = 0 (A1 = 0, so A(s) is constant): radius = +∞; certification reduces to one point | the radius formula is undefined at b = 0 |
| C4 | Failure and cost policy for the certificate: refine X, bisect within a budget, otherwise reject the configuration through the existing filter. **No silent fallback to an empirical certificate**; if the method is infeasible overall, a new design revision is required | audit request |
| C5 | Mandatory certificate self-tests, including negative tests across a known singularity, before the manifest is frozen. Review gate R0 is extended to the certificate and validated-residual code | rigor depends on implementation correctness |
| C6 | Taylor disk (§5.2) now uses the validated certificate at s = 0. Correction: every oracle formula except M10 is a **polynomial** in prefix and full inverse blocks, so no power-series inversion is needed. The earlier text wrongly invoked series inversion of W_j and Z⁻³ (Z⁻³ = (Ĝ_NN)³; W_j never appears in any oracle). M10's series is not computed. The tail bound in §5.2 (ii) is now a validated Cauchy-product majorant | mathematical correction |
| C7 | Reference accuracy (§4.6, §4.7) now uses the same validated exact-residual bound instead of an mp-evaluated residual | the residual was not validated |
| C8 | §4.6 classification updated. The zero test of §5.2 is explicitly a **non-validated numerical classification** | do not call it rigorous |
| C9 | The binomial tail bound in §5.2 (ii) is replaced by a ratio-test bound with a stated proof, α ≥ 1, and T counting coefficient magnitudes. The earlier /(1 − x)^m tail factor was asserted without proof | mathematical correction |
| C12 | §2.4.1 step 7: the claimed equivalence "ε_A = 0 iff ΔA = 0" is replaced by the implication ε_A = 0 ⟹ ΔA = 0. The converse is not guaranteed (enclosure radii stay positive). No other statement relied on the equivalence | logical correctness |
| C11 | The upper square root ū is defined piecewise, with ū(0) = 0 exactly. The earlier formula returned 2⁻ᵐ at f = 0, making the b = 0 branch unreachable and contradicting "b = 0 iff f_1 = 0". The upper-bound guarantee is unchanged for f > 0. The zero-case effects on ρ, a_k, n_k and ε_A are checked (§2.4.1 steps 4 and 7), and the Taylor-disk convention for b = 0 is stated | internal consistency |
| C10 | The λ-enclosure radius is stated as ≤ 2⁻²⁵⁵, a correct worst case. The historical D3 entry is annotated: its rigor claim is superseded by C1 | accuracy of claims |

## 0. Inspection of the previous experiment: verified definitions and reproducibility gaps

Verified from `coupled_blocks.py` and `wz_mutants.py` at `63e046d`:

| item | exact definition found |
|---|---|
| coupling parameter s | `V_j(s) = v0 · s · M_j` with v0 = 4 (`coupled_blocks.V`). **s multiplies the entire potential**, including the diagonal and identity terms. Values used: s ∈ {1, 2⁻⁴, 2⁻⁸, 2⁻¹²} |
| W/Z input | Δ_j = (E·I − diag(E_perp) − V_j(s))·dx² in double; leads from `lead_from_delta(fl(E − E_perp)·dx²)`; dx = 1/8, t = 32, L_y = 1 |
| mutant discrepancy D_mut | max over the four endpoint blocks of ‖X_mut − X_ref‖_F / ‖Ĝ_full(s)‖₂, against C(R_W) |
| slope | `np.polyfit(ln s, ln D_mut, 1)` over the 4 scales |
| "D ≫ τ_end" | D ≥ 100·τ_end |
| mutants | `wz_mutants.MUTANTS` (M1–M16, M4L, M4R). Pure order mutants M1, M2, M3, M5, M6, M7, M8, M15. `MIN_N` per mutant |

Previous gaps, documented and not repaired retroactively:
- G1. The fit procedure and "≫" existed only in code.
- G2. The maximum over blocks can switch block between s values.
- G3. The normalization depended on s.
- G4. One fit spanned the strong- and weak-coupling regimes, with no crossover rule.
- G5. The commuting and non-commuting parts of V both scale with s.
- G6. No noise-free mutant signal was computed.

v2 keeps **the same definition of s** (V_j(s) = s·V_j(1)) and the same W/Z definitions.

## 1. Research questions, mathematical facts, hypotheses

| RQ | question | answered by |
|---|---|---|
| RQ1 | Is the deviation near s = 1 a crossover rather than numerical noise? | exact vs double local slopes, and a preregistered crossover definition (§5.4) |
| RQ2 | Does the asymptotic law emerge at weak coupling? | the exact order q from Taylor coefficients (§5.2), and the double slope in a preregistered window (§5.3) |
| RQ3 | Is the scaling stable across independent matrices? | 120 held-out configurations; every meaningful class-R pair must satisfy the rules |
| RQ4 | Can ordering errors be distinguished from noise? | the oracle signal D* vs measured noise; SNR rules; controls (class D, M12, M9) |
| RQ5 | Can the criterion be fixed without tuning to previous slopes? | all thresholds derived in §4–§5 from precision arguments or algebra |

Notation:
- L(s) is the exact lattice with V(s) = sV(1) and exact inputs. Every block of M_L(s) is **affine** in s:
  M_L(s) = A0 + s·A1, with A1 = −blockdiag(V_j(1)) (physical units).
- The leads do not depend on s.
- R_W(s) is the represented W/Z lattice (double Δ_fl(s) taken as exact). Because of diagonal rounding it is **not** affine in s.

**Mathematical facts** (provable; used as validity conditions, not as hypotheses):
- **F1 (rationality and analyticity).**
  - On L(s), every reference block, every intermediate of the correct algorithm and every mutant output is a rational
    function of s. This holds because each is built from finitely many sums, products and inverses of matrices that are
    affine or rational in s.
  - Each such function is analytic on any real interval or complex disk on which all inverted matrices are invertible.
  - The correct algorithm inverts W_1 − I, I − W_j, Z_N, i.e. the Schur complements of the prefix matrices M^(j) and of M.
    By the determinant identity det M^(j) = det M^(j−1)·det D_j, these are all invertible wherever **every prefix matrix
    M^(j) (j = 1..N−1) and the full M are invertible**.
  - All mutants except M10 invert only these same matrices, or their negatives and transposes.
  - M10 inverts Z_N − (Λ − I), which is not a Schur complement of M; its invertibility is checked only pointwise (§2.4).
  - **The grid conditioning check alone does not prove invertibility between grid points.** Continuous-interval claims rely
    only on the certificate of §2.4.
- **F2.** At s = 0 every Δ_j and Λ is diagonal. Pure order mutants therefore coincide exactly with the correct algorithm.
- **F3 (order of vanishing).** Let F_m(s) be the discrepancy of the designated block X_m (mutant minus exact reference) on
  L(s). On the certified disk around 0, F_m(s) = Σ_{k≥0} s^k C_k.
  - For pure order mutants, C_0 = 0 (F2).
  - By the identity theorem, either F_m ≡ 0, or there is a **finite** smallest k = q with C_q ≠ 0. Then
    ‖F_m(s)‖_F = s^q ‖C_q‖_F (1 + O(s)) for real s → 0⁺, and the exact local slope tends to q.
  - **No upper bound on q is proven.** The commutator structure gives C_1 = 0 whenever the off-diagonal part of V(1) is
    zero, but C_1 may also vanish by symmetry in other cases. q is identified computationally (§5.2) and can be unresolved.
- **F4.** In class D, F_m ≡ 0 for every pure order mutant.
- **F5.** M12 is algebraically equivalent to the correct code, so F_m ≡ 0.
- **F6.** For M9, C_0 = 0: at s = 0 the clean device has W_j = w*, so the stale and current P coincide. M9 is not an
  ordering mutant and scales with s in every class, including D.
- **F7.** For M4, M4L, M4R, M10, M13, M14 and M16, generically C_0 ≠ 0 (q = 0).

**Empirical hypotheses** (tested on held-out data, per class):
- **H1 (attribution):** at every grid point resolvable for both neighbours (SNR ≥ 1e3), the double local slope equals the
  exact local slope to within 0.01.
- **H2 (asymptotic emergence):** in the preregistered window, the double OLS slope is within 0.02 of q and within 0.005 of
  the exact-data slope over the same points.
- **H3 (W/Z accuracy):** at every grid point, every correct-code accuracy quantity is ≤ 1e-10 (the fixed criterion).
- **H4 (distinguishability):** the double results of class-D order mutants and of M12 stay within the noise tolerance.
  Structural mutants are detected wherever resolvable.

Descriptive only: the distribution of q (expected mostly q = 1), the crossover scale s_×, and the regime model selected.

## 2. Configurations, generation and the conditioning filter

### 2.1 Fixed lattice

- dx = 1/8, t = 32, L_y = 1.
- E_perp,n from `stationary.transverse_basis(1.0, N_y)`.
- Clean leads.
- W/Z and production definitions as in `PREDICTIONS_coupled_blocks.md` §1 and §4.

### 2.2 Dimensions and classes

- N_y ∈ {2, 3, 4} and N ∈ {2, 3, 4, 6} (N = 1 has no order mutant on its code path). All matrices are ≤ 24 × 24.
- Classes:
  - **R:** real symmetric. This is the documented physical scope of `wz_rgf.py`.
  - **H:** complex Hermitian. Time-reversal is broken; this is outside the documented scope and is an algebraic stress test
    only.
  - **D:** diagonal (control).

### 2.3 Random matrix generation (exact dyadic values)

Root entropies are fixed now: set X (exploratory) **202611110001**; set C (held-out) **202611110002**; bootstrap
**202611110003**.

Per configuration index k (manifest order of §2.6), the generator is
`numpy.random.Generator(numpy.random.PCG64(numpy.random.SeedSequence(entropy=ROOT, spawn_key=(k, attempt))))`.

Fixed draw order:
1. `n_open = gen.integers(1, N_y + 1)`.
2. `u = gen.random()`.
3. For j = 1..N:
   - `R_j = gen.integers(-8, 9, size=(N_y, N_y))`;
   - `I_j = gen.integers(-8, 9, size=(N_y, N_y))` (always drawn).

Matrix construction:
- B_j = (R_j + i·I_j·[class H])/8.
- V_j(1) = v0·(B_j + B_jᴴ)/2 with **v0 = 8**. Entries are multiples of 1/2 in [−8, 8].
- Class D: zero the off-diagonal entries.
- V_j(s) = s·V_j(1), which is exact for dyadic s.

v0 = 8 places s = 1 in a strong-coupling regime (‖V dx²‖₂ ≲ 0.5, comparable to the δ spacings). This is a design choice for
RQ1, not tuned to any result.

### 2.4 Acceptance filter (deterministic; on failure, increment `attempt` and redraw everything)

1. **Energy rule:**
   - E ∈ (E_perp,n_open, E_upper), with E_upper = E_perp,n_open+1, or E_perp,1 + 64·(2 − δ_min) if n_open = N_y.
   - E = round(64·(E_lo + u·(E_upper − E_lo)))/64.
   - Every channel must satisfy |δ_n| ≥ δ_min = 0.02 and δ_n ≤ 2 − δ_min.
2. **Grid conditioning** (direct inverses only; never a recursion or solver). At every grid point s, on both L(s) and R_W(s),
   each of the following must be ≤ 1e3:
   - cond₂(M), cond₂(M^(j)), cond₂(I − W_j) (from prefix inverses), cond₂(Z_N) (= cond₂ of [(M/t)⁻¹]_NN⁻¹);
   - cond₂(Z_N − (Λ − I)), the M10 matrix.
3. **Continuous-interval invertibility certificate on [0, 1] (validated, on L(s); method in §2.4.1).**
   - For each A(s) ∈ {M_L(s), M_L^(j)(s), j = 1..N−1}, all affine: A(s) = A(s_k) + (s − s_k)A1.
   - At each certification point s_k, §2.4.1 yields validated rational upper bounds:
     - a_k ≥ ‖A(s_k)⁻¹‖₂;
     - n_k ≥ ‖A(s_k)‖₂;
     - b ≥ ‖A1‖₂.
   - **Neumann step** (operator 2-norm, valid for real and complex s): if |s − s_k|·a_k·b < 1, then A(s) is invertible and
     ‖A(s)⁻¹‖₂ ≤ a_k / (1 − |s − s_k| a_k b).
   - Coverage uses the **half radius** r_k = 1/(2 a_k b), computed exactly as a rational. On |s − s_k| ≤ r_k this gives the
     validated bounds:
     - ‖A(s)⁻¹‖₂ ≤ 2a_k;
     - ‖A(s)‖₂ ≤ n_k + r_k b;
     - cond₂(A(s)) ≤ 2a_k(n_k + r_k b).
   - **b = 0** (equivalently A1 = 0, i.e. V_j(1) = 0 on every slice of that matrix; see §2.4.1 step 7): A(s) ≡ A(s_k) for all
     s. Then r_k := +∞ (used only in the coverage comparisons, never multiplied by b), one
     validated point certifies all of [0, 1], and the bounds are a_k and n_k with no s-dependence. The formula 1/(2a_k b) is
     never evaluated.
   - **Certification points:**
     - Start with s = 0, the 81 grid points (§3) and s = 1, all exact dyadics.
     - Sort the points. Coverage of [0, 1] holds iff 0 and 1 are covered and every consecutive pair satisfies
       s_{k+1} − s_k ≤ r_k + r_{k+1}. These are exact rational comparisons.
     - Every uncovered gap is bisected at its dyadic midpoint, to a depth of at most 30 per gap. The total budget is 2000
       certification points per matrix A.
   - **Accept only if [0, 1] is covered for every A.**
     - Then, for **every** s ∈ [0, 1] (not only grid points), every matrix in the list is invertible, with the recorded
       validated bounds. The closed intervals [s_k − r_k, s_k + r_k] cover the compact interval [0, 1].
     - Record the validated continuous condition-number bound and the validated minimum r_k.
   - **Taylor disk:** the s = 0 certificate also certifies the complex disk |s| ≤ r_0 = 1/(2 a_0 b) around 0 (or all of ℂ if
     b = 0). That disk is the analytic neighbourhood used in §5.2.
   - **Not certified continuously:** R_W(s), which is not affine (diagonal rounding) and has pointwise checks only; and the
     M10 matrix Z − (Λ − I) (pointwise only, §2.4 item 2).
   - The grid conditioning check (item 2) remains an **empirical** filter on R_W and L. It is not a certificate.

#### 2.4.1 Validated method for the inverse-norm bound (selected), and the alternative considered

**Exact representation of the inputs.**
- Every double is a dyadic rational m·2^e with m an integer: E, E_perp,n, V_j(1) (multiples of 1/2), s (dyadic grid and
  bisection midpoints), dx² = 2⁻⁶ and t = 32.
- So every entry of A(s) is an exact dyadic rational, **except** the lead terms −tλ_n. The lead roots λ_n are determined by
  the exact rational c_n = 1 − δ_n, δ_n = (E − E_perp,n)dx²:
  - open channel (|c_n| < 1): λ_n = c_n + i√(1 − c_n²);
  - closed channel (c_n > 1): λ_n = c_n − √(c_n² − 1). The design excludes c_n < −1 via δ_n ≤ 2 − δ_min.
- **Enclosure of λ.** For a rational q > 0 and m = 256 bits:
  - lo = isqrt(⌊q·4^m⌋)/2^m ≤ √q ≤ hi = (isqrt(⌈q·4^m⌉) + 1)/2^m, using exact Python integer `math.isqrt` (floor square root).
  - Set λ̃_n to the dyadic midpoint of the resulting enclosure of λ_n, with a rigorous radius e_λ,n = (hi − lo)/2 ≤ 2⁻²⁵⁵
    on the affected (real or imaginary) part.
- The enclosure splits A(s) = A_c(s) + ΔA:
  - A_c(s) is an exact dyadic matrix, with λ_n replaced by λ̃_n;
  - ΔA is nonzero only on the diagonal of slices 1 and N, with |ΔA_nn| ≤ t·e_λ,n per lead. For N = 1 both leads act on one
    block: |ΔA_nn| ≤ 2t·e_λ,n.
  - Validated bound: ‖ΔA‖_F ≤ ε_A := ū(Σ of the squared entry bounds), with ū from step 4 below. If no lead term is enclosed
    (R_orig, whose σ_code entries are exact doubles), the sum is 0 and ε_A = 0 exactly, consistent with ΔA = 0 in §4.7.
  - ΔA does not depend on s, and A1 is exactly dyadic.

**Approximate inverse with an exactly computed residual.**
1. X: compute an approximate inverse of A_c(s_k) in mpmath at dps 60 (dps 120 on retry). Round every entry to a dyadic with
   a 256-bit mantissa. X is then an **exact** dyadic matrix. How X was obtained is irrelevant to rigor.
2. R = I − A_c(s_k)·X: computed **exactly** in dyadic integer arithmetic (Python integers or `fractions.Fraction`). There is
   **no rounding at all**.
3. Exact rationals: f_R = ‖R‖_F² = Σ |R_ab|² (complex entries: re² + im²), f_X = ‖X‖_F², f_c = ‖A_c(s_k)‖_F² and
   f_1 = ‖A1‖_F².
4. Rigorous upper square roots for rationals f ≥ 0 (m = 256), defined **piecewise**:
   - ū(0) := 0 (exact);
   - ū(f) := (isqrt(⌈f·4^m⌉) + 1)/2^m for f > 0.

   Guarantees:
   - **f > 0:** with M = ⌈f·4^m⌉ ≥ f·4^m, isqrt(M) = ⌊√M⌋, so isqrt(M) + 1 > √M ≥ 2^m√f. Hence ū(f) > √f > 0.
   - **f = 0:** ū(0) = 0 = √0.
   - So ū(f) ≥ √f for every f ≥ 0, and **ū(f) = 0 iff f = 0**. The f = 0 branch returns the exact value, so no upper-bound
     guarantee is weakened.
   - Every f to which ū is applied (f_R, f_X, f_c, f_1 and the sum in ε_A) is an exact non-negative rational, so the branch is
     decided by an exact comparison with 0.
5. ρ := ū(f_R) + ε_A·ū(f_X). Since I − A(s_k)X = R − ΔA·X, this gives
   ‖I − A(s_k)X‖₂ ≤ ‖R‖_F + ‖ΔA‖_F‖X‖_F ≤ ρ (using ‖·‖₂ ≤ ‖·‖_F and submultiplicativity).
6. If ρ < 1 (an exact rational comparison):
   - A(s_k)X = I − R′ with ‖R′‖₂ ≤ ρ < 1, so A(s_k)X is invertible, and A(s_k) is invertible because it is square.
   - A(s_k)⁻¹ = X(I − R′)⁻¹, so ‖A(s_k)⁻¹‖₂ ≤ ‖X‖₂/(1 − ρ) ≤ a_k := ū(f_X)/(1 − ρ), exact rational.
   - This uses the **operator 2-norm**, for which ‖I‖₂ = 1 and ‖(I − R′)⁻¹‖₂ ≤ 1/(1 − ρ). Frobenius norms appear only as
     upper bounds of 2-norms; a Frobenius-only Neumann bound would be wrong, since ‖I‖_F = √n.
7. n_k := ū(f_c) + ε_A ≥ ‖A(s_k)‖₂, and b := ū(f_1) ≥ ‖A1‖₂. A1 is an exact dyadic matrix, so f_1 = ‖A1‖_F² = 0 iff A1 = 0.
   By the piecewise definition, **b = 0 iff A1 = 0**: the b = 0 branch of §2.4 item 3 is reachable exactly when A(s) is
   constant. For A1 ≠ 0, b > 0, so r_k = 1/(2a_k b) is a finite positive rational.
   - b is computed per matrix A ∈ {M_L, M_L^(j)}. A prefix matrix whose slices 1..j all have V(1) = 0 takes the b = 0 branch
     even if the full matrix does not.
   - Effect of the zero branch on the other uses:
     - f_R = 0 (X exactly inverse to A_c) gives ρ = ε_A·ū(f_X), which is valid.
     - f_X = 0 means X = 0, so R = I, ρ ≥ 1 and the point is uncertified. a_k is never formed with f_X = 0 and ρ < 1.
     - f_c = 0 gives n_k = ε_A ≥ ‖ΔA‖₂ = ‖A(s_k)‖₂, which is valid.
     - ε_A = 0 ⟹ ΔA = 0, because ‖ΔA‖_F ≤ ε_A. **The converse is not guaranteed.** ε_A is built from the enclosure radii
       e_λ,n, which are positive for every enclosed lead term (hi − lo ≥ 2⁻ᵐ). So ε_A > 0 can hold even when the actual
       error ΔA happens to be 0, for example if λ_n is a dyadic rational equal to its midpoint λ̃_n.
       - In this design ε_A = 0 occurs only when no lead term is enclosed (R_orig).
       - Every use of ε_A needs only the upper-bound property ‖ΔA‖_F ≤ ε_A, so a conservative positive ε_A is always
         valid; it merely enlarges ρ and n_k slightly.
     - No division by ū(·) occurs anywhere; the only divisions are by 1 − ρ > 0 and by 2a_k b with b > 0.
8. If ρ ≥ 1: retry once with X recomputed at dps 120. If that fails too, s_k is uncertified, and its neighbours must cover it
   through bisection.

**Rounding coverage.**
- Steps 2–7 are exact rational arithmetic. The only approximations are the two square-root enclosures (inputs λ and outputs
  ū), and both are rigorous by construction from the exact integer floor square root.
- No floating-point result enters any inequality except through X, which is used as an exact dyadic matrix.
- Rigor therefore rests on:
  - (i) the mathematical inequalities above;
  - (ii) correctness of Python's exact integer arithmetic and `math.isqrt`;
  - (iii) correctness of the implementation, which is addressed by the self-tests and review below.

**Validity on the interval.** Given validated a_k, b and n_k, the Neumann step of item 3 holds exactly for every real or
complex s with |s − s_k| ≤ r_k. The coverage condition is checked in exact rationals. The bound therefore holds for every s in
[0, 1], not merely at evaluated points.

**Alternative considered: directed-rounding interval arithmetic.**
- mpmath provides an interval context (`mpmath.iv`, with real `ivmpf` and complex `ivmpc` intervals using outward rounding).
- An interval Gaussian elimination, or a Krawczyk-type enclosure of A(s_k)⁻¹, would also be valid. Its disadvantages:
  - (a) its rigor depends on the correctness of the library's directed rounding for every operation, including complex
    ones, which is harder to audit;
  - (b) wrapping and dependency effects inflate the enclosure with n (≤ 24) and could reject certifiable points;
  - (c) it still needs a validated norm of an interval matrix.
- The exact-residual method needs only exact integer arithmetic plus a scalar square-root enclosure. It is selected.
- Interval arithmetic may be used as an optional independent cross-check of a_k. It has no decision role.

**Self-tests (mandatory, at implementation time, before the manifest freeze; failure → the certificate code is not used):**
- (T1) Exactly invertible matrices with known inverses (diagonal, permutation, 2×2 complex): a_k ≥ the true ‖A⁻¹‖₂.
- (T2) Negative test: an affine family with a known singular point, e.g. A(s) = diag(1 − 2s, 1) and a 24 × 24 analogue. The
  covering must **not** be achieved across s = 1/2. Any certified interval [s_k − r_k, s_k + r_k] that contains s = 1/2 is a
  fatal error.
- (T3) λ enclosure: the true root lies in the enclosure, checked against an independent dps-300 evaluation.
- (T4) b = 0 path: ū(0) = 0 exactly; a matrix with A1 = 0 is certified on [0, 1] from one point; no expression 1/(2a_k b)
  or r_k·b is evaluated.
- Review gate R0 (§4.3) is extended to the certificate and validated-residual code.

**Failure and cost policy:**
- If [0, 1] is not covered for some matrix within the depth and point budget, the **configuration fails the acceptance
  filter** and is redrawn (§2.4); this is recorded.
- If the validated certificate is infeasible for the manifest as a whole (systematic runtime or memory failure), the
  experiment **must not** fall back to an empirical certificate under this design. A new design revision is required, which
  would downgrade F1/F3 interval claims to grid-point claims.
4. **Prohibited selection.** No selection on commutators, mutant signals, slopes or any solver output.

### 2.5 Conditional sampling (explicit limitation)

The accepted configurations are draws from the generator distribution **conditional on passing §2.4**.

Consequences:
- Every claim of this experiment refers to that conditional population only.
- Devices that are poorly conditioned, near-singular, resonant, near a bound state, or near a channel threshold
  (|δ_n| < 0.02) are **systematically excluded**. Nothing can be inferred about W/Z behaviour, or about mutant detectability,
  for them.
- In such devices the Schur complements can be nearly singular, rounding can be amplified far beyond the bounds used here,
  and the analytic structure (F1, F3) can have nearby poles.

Recorded for transparency (reported, never used to choose thresholds):
- the rejection log, with the rule that rejected and the attempt count;
- the rejection rate per stratum (N_y, class, N);
- for each accepted configuration, the certified minimum distance to singularity, min over A of r_k at the most
  constraining point, and the certified continuous condition-number bound.

A high rejection rate in a stratum narrows the population that the stratum represents. This is reported, but it is not a
decision criterion.

### 2.6 Configuration sets

Index order: N_y ∈ (2, 3, 4) → class ∈ (R, H, D) → N ∈ (2, 3, 4, 6) → draw d.

| set | R | H | D | total |
|---|---|---|---|---|
| X (exploratory) | 3 × 4 × 2 = 24 | 24 | 3 × 4 × 1 = 12 | 60 |
| C (held-out) | 3 × 4 × 4 = 48 | 48 | 3 × 4 × 2 = 24 | 120 |

## 3. Coupling grid

- **Grid:** s ∈ {2⁻ᵏ : k = 0..40} ∪ {3·2⁻⁽ᵏ⁺²⁾ : k = 0..39}.
  - 81 dyadic points from 1 to 2⁻⁴⁰ (about 6.6 per decade), so s·V(1) is exact in double.
  - s = 1 is a diagnostic point.
- **Lower edge, justified by precision:**
  - The SNR noise level satisfies N_X ≥ τ_abs ≥ 1e-13·‖Ĝ_full‖₂ (§4.5).
  - For an order-1 signal with |C_1| ≤ 1, SNR < 10 before 2⁻⁴⁰. (|C_1| ≤ 1 is assumption A2; it is used only for this
    sizing argument, never in a decision.)
  - The grid therefore extends past the detectability limit, so the transition into the noise floor is observed.
  - At 2⁻⁴⁰ the diagonal s-part of Δ_fl is still represented, with relative rounding ~1e-3. This is accounted for exactly
    in R_W.
- **Asymptotic order** is not taken from grid points. It comes from Taylor coefficients (§5.2).
- **Crossover:** §5.4. No grid choice uses previous slopes.

## 4. Quantities, references, mutant oracle and numerical controls

### 4.1 Solvers and exact computations per (configuration, s)

- **A:**
  - `wz_rgf.wz_endpoint_blocks` and `wz_fisher_lee` (read-only).
  - The instrumented copy `wz_mutants.wz_instrumented`, unchanged from `63e046d`, must be bit-identical at the four endpoint
    blocks.
  - The double mutants come from `wz_mutants.py`, unmodified.
- **B:** production `stationary.scattering_matrix(E, dev, min_closed=0)` on a directly built `Device`. Accuracy is reported
  only; it carries no decision weight.
- **C (reference):**
  - The direct mpmath inverse of the full block-tridiagonal matrix, assembled from its definition, with Σ_L and Σ_R added
    separately.
  - Built for L(s), R_W(s) and R_orig(s) at dps 60.
  - Prefix references from the same routine without Σ_R.
  - No recursion code is used.
- **D\*_oracle (primary exact mutant signal):** the closed forms of §4.2, evaluated from the exact prefix and full inverses
  of R_W(s) and of L(s).
- **D\*_transcription:** a clean-room mpmath transcription of every mutant (§4.3), run on R_W(s) and L(s).

### 4.2 Exact algebraic change of every mutant, and the specification-derived oracle

Notation (all exact, from direct inverses of R_W(s) or L(s) in Ĝ = tG units):
- Ĝ^(j)_ab is a block of the prefix inverse (M^(j))⁻¹·t.
- P_j = −Ĝ^(j)_jj, with P_1 := (I − W_1)⁻¹ = −Ĝ^(1)_11.
- W_j = I + (Ĝ^(j)_jj)⁻¹.
- Z⁻¹ = Ĝ_NN of the full inverse.
- Abbreviations: L_j := Ĝ^(j)_j1, R_j := Ĝ^(j)_1j, T_j := Ĝ^(j−1)_{1,j−1} P_j Ĝ^(j−1)_{j−1,1}, Ĝ_ab := full endpoint blocks.
- A sum Σ_j runs over the loop slices j = 2..N−1.

Correct code (from `wz_mutants.py`, `mutant=None`, identical to `wz_rgf.wz_endpoint_blocks`):
- init: `W = 2Δ_1 + Lm`; `G11 = solve(W − I, I)`; `Gj1 = G1j = G11`.
- loop: `W = 2Δ_j + solve(I − W, W)`; `P = solve(I − W, I)`; `G11 = G11 − G1j@P@Gj1`; `Gj1 = P@Gj1`; `G1j = G1j@P`.
- closing: `Z = 2Δ_N + Lm + solve(I − W, W)`; `Zi = solve(Z, I)`; `G11 = G11 + G1j@Zi@Gj1`; `GN1 = −Zi@Gj1`;
  `G1N = −G1j@Zi`; `GNN = Zi`.

| mutant | exact code change (`wz_mutants.py`) | designated block X_m | oracle: exact mutant value of X_m |
|---|---|---|---|
| M1 | `Gj1 = Gj1@P` | Ĝ_N1 | −Z⁻¹ Ĝ^(1)_11 P_2 P_3 ⋯ P_{N−1} |
| M2 | `G1j = P@G1j` | Ĝ_1N | −(P_{N−1} ⋯ P_2 Ĝ^(1)_11) Z⁻¹ |
| M3 | both M1 and M2 | Ĝ_N1 | as M1 |
| M4 | chains updated first, then `G11 = G11 − G1j@P@Gj1` | Ĝ_11 | Ĝ_11 − Σ_j Ĝ^(j−1)_{1,j−1}(P_j³ − P_j)Ĝ^(j−1)_{j−1,1} |
| M4L | `Gj1 = P@Gj1` before the G11 update; `G1j` after | Ĝ_11 | Ĝ_11 − Σ_j Ĝ^(j−1)_{1,j−1}(P_j² − P_j)Ĝ^(j−1)_{j−1,1} |
| M4R | `G1j = G1j@P` before the G11 update; `Gj1` after | Ĝ_11 | as M4L |
| M5 | `G11 = G11 − Gj1@P@G1j` | Ĝ_11 | Ĝ_11 − Σ_j (Ĝ^(j−1)_{j−1,1} P_j Ĝ^(j−1)_{1,j−1} − T_j). The j = 2 term is 0 |
| M6 | `GN1 = −(Gj1@Zi)` | Ĝ_N1 | −L_{N−1} Z⁻¹ |
| M7 | `G1N = −(Zi@G1j)` | Ĝ_1N | −Z⁻¹ R_{N−1} |
| M8 | `G11 = G11 + Gj1@Zi@G1j` (closing) | Ĝ_11 | Ĝ_11 + L_{N−1}Z⁻¹R_{N−1} − R_{N−1}Z⁻¹L_{N−1}. Equals 0 change at N = 2 |
| M9 | `P = solve(I − W_old, I)` before the W update | Ĝ_N1 | −Z⁻¹ P_{N−2} ⋯ P_2 P_1 Ĝ^(1)_11 (stale P_j → P_{j−1}) |
| M10 | `Z = 2Δ_N + solve(I − W, W)` (Lm dropped) | Ĝ_NN | (Z − (Λ − I))⁻¹ |
| M12 | `solve((I−W).T, W.T).T` for (I − W)⁻¹W (loop and Z) | Ĝ_11 | Ĝ_11 (no change; algebraic identity) |
| M13 | `G11 = G11 + G1j@P@Gj1` (loop) | Ĝ_11 | Ĝ_11 + 2 Σ_j T_j |
| M14 | init `G11 = solve(I − W, I)` | Ĝ_11 | Ĝ_11 − 2 Ĝ^(1)_11 (and Ĝ_N1 → −Ĝ_N1, Ĝ_1N → −Ĝ_1N) |
| M15 | M6 + M7 | Ĝ_N1 | as M6 |
| M16 | closing: GN1, G1N first; `G11 = G11 + G1N@Zi@GN1` | Ĝ_11 | Ĝ^(N−1)_11 + R_{N−1} Z⁻³ L_{N−1} |
| M11 | N = 1 boundary `2Δ_1 + Lm` | — | off-path for every N ≥ 2 (excluded) |

The oracle values are derived from the code changes by exact algebra, using the identities (verified at ≤ 1e-50) that the
correct recursion's intermediates equal the prefix-reference quantities. The oracle therefore uses **no recursion code**.

### 4.3 Independent reference for the mutants, and checks for shared conceptual errors

Three implementations of each mutant output:
1. the double mutant (`wz_mutants.py`);
2. the mpmath transcription;
3. the oracle (§4.2).

Required independence measures:
- **Clean-room transcription.** The mp transcription is written from §4.2's "exact code change" column and the correct-code
  listing only, without reading `wz_mutants.py`.
  - If the implementer has read `wz_mutants.py`, the transcription must instead be written by an independent implementer: a
    separate person or a fresh agent session given only this document.
  - Authorship is recorded.
- **Review gate R0 (mandatory, before any solver run, including the exploratory run).** An independent reviewer, who wrote
  neither the oracle nor the transcription, verifies the following and records a sign-off in
  `REVIEW_coupling_scaling_v2.md` (reviewer identity, date, reviewed commit):
  - (a) §4.2's code-change column against `wz_mutants.py`, line by line;
  - (b) each oracle formula against the code change, by an independent derivation;
  - (c) the oracle implementation against §4.2;
  - (d) the clean-room transcription against §4.2.
  - Running without a committed R0 sign-off → INVALID.
- **Three-way agreement** at every grid point (checks; violation → INVALID):
  - oracle vs transcription: ≤ 1e-50 relative to ‖Ĝ_full‖₂;
  - double vs oracle: |D_m − D*_oracle| ≤ 10·N_X(s) wherever SNR ≥ 10 (§4.5).
  - A shared error between (1) and (2), for example a common misreading of a mutant, shows up as disagreement with (3),
    which is derived independently from the algebra.
- **First-deviation check.** The first intermediate that differs in the transcription must equal the prefix-level closed
  form at the slice where the mutant first acts (e.g. M1: Ĝ^(1)_11 P_2 in place of P_2 Ĝ^(1)_11), to within 1e-50.
- **Residual risk (stated).** An error in §4.2's change column itself, shared by all three implementations, is detectable
  only by review step (a). No numerical check can detect it.

### 4.4 Fitted error quantity

- D_m(s) = ‖X_m^mut,double(s) − X_m^ref(s)‖_F (absolute, Ĝ units; no maximum over blocks; no s-dependent normalization).
- D*_m(s) is the same quantity with the oracle value (primary) or the transcription (check).
- All four blocks are recorded as secondary data.

### 4.5 Noise and SNR

- **Observed noise:** η_X(s) = ‖X^correct,double − X^ref‖_F (empirical).
- **Noise level:** N_X(s) = max(η_X, τ_abs, ε·‖X_ref‖_F), with τ_abs = τ_end·‖Ĝ_full‖₂ (heuristic; `PREDICTIONS_coupled_blocks.md` §5).
- **SNR\*:** D*_oracle / N_X.

| SNR* | meaning |
|---|---|
| ≥ 10 | detection-eligible |
| ≥ 1e3 | regression-usable |
| < 10 | **no claim; never interpreted as harmless** |

- **Detection** means D_m ≥ 10·N_X **and** |D_m − D*_oracle| ≤ 10·N_X.
- Including the heuristic τ in N_X makes SNR conservative. If the empirical η exceeds τ, N_X still contains η, so the SNR
  logic remains valid.

### 4.6 Classification of every error quantity and threshold

| quantity | nature | role in decisions |
|---|---|---|
| 1e-10 accuracy criterion (endpoint, prefix, W_j, P_j, ΔS) | **fixed acceptance criterion**, inherited unchanged | the only FAIL-W criterion |
| judgeability gate τ_end ≤ 1e-11 (inherited) | heuristic | a non-judgeable held-out case → DESIGN INSUFFICIENT (it cannot be skipped silently) |
| τ_end, τ_pre, τ_W, τ_P, τ_S and the U/Q/symmetry bounds | **heuristic** first-order block-elimination estimates (no theorem) | exceedance → anomaly → blocks PASS → DESIGN INSUFFICIENT; **never FAIL**; also enters N_X |
| reference accuracy: ‖A⁻¹ − X‖₂ ≤ ‖X‖₂·ρ/(1 − ρ), with ρ ≥ ‖I − AX‖₂ computed by the exact-residual method of §2.4.1 (X = the mp reference inverse rounded to dyadics; λ enclosed) | **validated** (exact residual; integer square-root enclosures) | ρ/(1 − ρ) > 1e-55 → escalation (dps 120), then INVALID |
| Neumann invertibility certificate, validated inverse-norm bound and continuous cond bound (§2.4.1) | **validated** (exact residual, exact rational coverage; conditional on the self-tests and the R0 review) | acceptance filter |
| Taylor tail majorant (§5.2 (ii)) | **validated** (from validated a, b) | INVALID check |
| Zero test of Taylor coefficients (§5.2) | **non-validated numerical classification** (precision scaling); not a proof that a coefficient is zero | order identification; unresolved → DESIGN INSUFFICIENT |
| dps-100 vs dps-60 agreement; `lu_solve` vs `inverse` | empirical cross-checks | INVALID on disagreement |
| η_X | **empirical** | noise estimate |
| SNR levels 10 and 1e3; detection factor 10; apparatus tolerance 10·N_X | design constants (the latter built on empirical η and heuristic τ) | detection, regression usability, INVALID |
| slope tolerances 0.01 (H1), 0.02 and 0.005 (H2), window adequacy 0.02 | derived from heuristic A1 plus the empirical noise level | FAIL-D (H1, H2); window adequacy → DESIGN INSUFFICIENT |
| zero-test thresholds of §5.2 | design constants justified by precision scaling | order identification |
| conditioning ≤ 1e3; δ_min = 0.02 | design constants | acceptance filter |

### 4.7 Represented-lattice and reference controls (violation → INVALID unless stated)

1. Off-diagonal Δ_fl(s) entries are bitwise −s·V_off(1)·dx². Diagonal rounding of the s-part is recorded.
2. e_glob(C(R_W) − C(L)) and e_glob(C(R_orig) − C(L)) are ≤ 10·cond₂(M)·ε.
3. |D*_oracle,R_W(s) − D*_oracle,L(s)| / D*_oracle,L(s) ≤ 1e-6 wherever D*_L ≥ 1e-30. Rationale: the lattices differ by O(ε)
   in commuting diagonal parts, which perturb C_1 by relative O(ε); the s-part rounding enters at O(s²).
4. Reference checks:
   - (a) validated a-posteriori bound ρ/(1 − ρ) ≤ 1e-55 (§2.4.1 exact-residual method, applied to the mp reference inverse)
     at every grid point, for the full and prefix systems of L and R_W.
     - R_W entries are exact dyadics except the lead terms (λ of δ_lead,fl, enclosed as in §2.4.1). R_orig entries are all
       doubles (σ_code), so ΔA = 0.
     - R_orig is validated at every grid point if the runtime allows. Otherwise, preregistered, it is validated at k ≡ 0
       (mod 4) only, and its remaining points are labelled empirical. B carries no decision weight, so this does not affect
       any decision. L and R_W validation is mandatory at every grid point.
     - The measured double errors are then bounded rigorously by measured value ± the validated reference error, which is
       ≤ 1e-55 relative and negligible against 1e-10.
   - (b) dps 100 agreement within 1e-55 at k ≡ 0 (mod 4), at s = 1, at 2⁻⁴⁰, and on the duplicate set;
   - (c) `mp.lu_solve` vs `mp.inverse` at s ∈ {1, 2⁻¹², 2⁻⁴⁰}.
   - Escalation to dps 120 on failure; INVALID if it persists.
5. Prefix and closing identities hold at ≤ 1e-50.
6. The instrumented copy is bit-identical to `wz_endpoint_blocks` at every execution.
7. Determinism: the first 5 configurations of X and of C, over the full grid with `OPENBLAS_NUM_THREADS=1`, are run twice
   and must be bit-identical.
8. Manifest hash and R0 sign-off: §8.

## 5. Scaling analysis (all procedures fixed before data collection)

### 5.1 Meaningfulness classification per (configuration, mutant); exact rules only

| category | rule |
|---|---|
| off_path | N < `MIN_N[m]`, or m = M11 |
| identical_by_algebra | M5 at N = 3, M8 at N = 2 (oracle change term exactly 0) |
| equivalent | M12 |
| class_D_order | a pure order mutant in class D |
| ordering-meaningful | a pure order mutant in class R or H, on path, not identical, with resolved finite order q ≥ 1 (§5.2) |
| ordering-zero | a pure order mutant in R/H whose coefficients C_0..C_K_max all test zero **and** whose oracle is ≤ the identity threshold at every grid point. A structural zero: excluded from detection claims and counted against coverage |
| order-unresolved | the zero test of §5.2 is ambiguous, or the order exceeds K_max |
| perturbation-control | M9 |
| structural | M4, M4L, M4R, M10, M13, M14, M16 |

Identical, equivalent, class_D_order and ordering-zero pairs are **never** counted as detections.

Expected ordering pairs (on-path, non-identical order mutants per N: N = 2 → 3, N = 3 → 7, N = 4 → 8, N = 6 → 8; 26 per
(N_y, class, draw) quadruple):

| set | R | H |
|---|---|---|
| C (held-out) | 312 | 312 |
| X (exploratory) | 156 | 156 |

### 5.2 Identification of the asymptotic order (replaces any restriction on q)

- **Taylor coefficients.** On L(s), each inverse needed by the oracle is expanded around s = 0 as the Neumann series
  A(s)⁻¹ = Σ_k (−s)^k (A0⁻¹A1)^k A0⁻¹, truncated at K_max = 12.
  - This covers every prefix matrix and the full M. Z⁻¹ = Ĝ_NN is a block of the full inverse; P_j = −Ĝ^(j)_jj and the prefix
    blocks come from the prefix inverses.
  - **Every oracle formula of §4.2 except M10 is a polynomial in these blocks.** For example Z⁻³ = (Ĝ_NN)³ and
    P_1 = −Ĝ^(1)_11; W_j does not appear in any oracle. So **no power-series inversion is needed**: only Neumann series and
    their Cauchy products.
  - M10 needs (Z − (Λ − I))⁻¹, which is not certified. Its series is not computed. M10 is structural and outside the scaling
    decision.
  - Evaluating the oracle formulas in truncated power-series arithmetic yields C_0..C_12 for the designated block of every
    other mutant.
  - The coefficients are computed in floating mp arithmetic at two precisions, so they are **approximate**. The zero test
    below is a numerical classification, not a validated proof.
  - Convergence of the series on |s| ≤ r_0 = 1/(2a_0 b) (r_0 = +∞ when b = 0) follows from the validated certificate at s = 0 (§2.4 item 3):
    ‖A0⁻¹A1‖₂ ≤ a_0 b, so every Neumann series converges absolutely there, and so do finite sums and products of them.
- **Zero test (precision scaling).** Compute everything at dps 100 and at dps 200. With z_k(dps) = ‖C_k‖_F / ‖Ĝ_full(0)‖₂:
  - C_k is **nonzero** if z_k(200) ≥ 1e-60 and |z_k(100) − z_k(200)| ≤ 1e-20·z_k(200).
  - C_k is **zero** if z_k(100) ≤ 1e-80 and z_k(200) ≤ 1e-170. An exactly zero coefficient computed in floating arithmetic
    scales with the working precision; a genuine coefficient does not.
  - Otherwise C_k is **unresolved**.
  - q = the smallest k with C_k nonzero, provided all C_j (j < k) test zero.
  - If any C_j (j ≤ q) is unresolved, or all of C_0..C_12 test zero but the oracle is not ≡ 0 on the grid, the pair is
    **order-unresolved** → DESIGN INSUFFICIENT for that pair. It is never INVALID merely because the order is large or the
    coefficients are tiny.
- **Consistency checks** (genuine contradictions → INVALID):
  - (i) **Pure order mutants and M9:** C_0 must test zero (F2, F6). If C_0 is nonzero → INVALID. If it is unresolved → this
    is precision-limited, not a contradiction → DESIGN INSUFFICIENT for that pair.
  - (ii) **Series vs direct evaluation:** at s = 2⁻⁶⁰, the truncated series Σ_{k≤12} C_k s^k must agree with the oracle
    evaluated directly from mp inverses of L(2⁻⁶⁰), to relative 1e-40.
    - **Validated tail majorant.** Each block series F_A(s) = Σ_k s^k F_{A,k} of a matrix A (in Ĝ units, factor t)
      satisfies ‖F_{A,k}‖₂ ≤ t·a_A,0·(a_A,0 b_A)^k, with a_A,0 and b_A the validated s = 0 bounds of §2.4.1.
      - b = 0 case (convention 0⁰ = 1): a constant matrix (b_A = 0) has only its k = 0 coefficient, and the bound gives 0
        for k ≥ 1, as required.
      - If b_A = 0 for every A, then β = 0, x = 0, r_0 = +∞ (the disk is all of ℂ) and the tail bound is exactly 0.
      - Let β = max_A a_A,0·b_A, x = |s|β < 1/2, and α = max(1, max_A t·a_A,0).
      - Write the discrepancy F_m (mutant formula minus the reference block) as a sum of monomials of degree ≤ m in such
        blocks, with T = the sum of the absolute values of the scalar coefficients. T and m are read off §4.2 per mutant
        and N; for example, M13 contributes coefficient 2 per loop term.
      - The coefficients of the oracle series are then majorized by those of T·α^m (1 − x)^(−m). Taking α ≥ 1 makes this
        majorant valid for every monomial of degree ≤ m.
      - The tail beyond order K = 12 is bounded by T·α^m·C(K + m, m − 1)·x^(K+1)/(1 − r·x), with r = (K + m + 1)/(K + 2).
        - Proof: (1 − x)^(−m) = Σ_k C(k + m − 1, m − 1) x^k. The ratio of consecutive terms is (k + m)/(k + 1)·x, which
          decreases in k, so for k ≥ K + 1 it is ≤ r·x. The tail is therefore at most its first term C(K + m, m − 1)x^(K+1)
          times a geometric factor 1/(1 − r·x).
        - The bound requires r·x < 1. This is checked exactly; if it fails, the check point s is halved until it holds.
      - The bound is evaluated in exact rationals.
    - The check uses s = 2⁻⁶⁰ if 2⁻⁶⁰ ≤ r_0, otherwise s = r_0·2⁻²⁰.
    - A disagreement larger than this validated tail bound plus the mp precision allowance (10^(−(dps−10)) relative) →
      INVALID (reference, oracle or series error). The precision allowance is the only non-validated term.
  - (iii) **Class D and M12:** every coefficient must test zero and the grid oracle must vanish (F4, F5). Otherwise →
    INVALID.
- **Finite-scale estimates.** Local slopes on the grid are **never** used to identify q, and a non-integer value at finite s
  is **never** INVALID. Their departure from q is the subject of §5.3 (window adequacy) and §5.4 (crossover).

### 5.3 Weak-coupling window and regression

**Window (fixed):** W_asym = {grid s ≤ 2⁻¹²} ∩ {SNR*(s) ≥ 1e3}.
- The edge follows from the first-order remainder. The local slope deviates from q by ≈ ρs/(1 + ρs). Assuming
  **A1: ρ ≤ 25** (heuristic), a model error ≤ 0.01 requires s ≤ 4e-4, giving the first dyadic grid point 2⁻¹².
- A1 is computable per pair from the Taylor coefficients: ρ_est = ‖C_{q+1}‖/‖C_q‖. ρ_est is **reported**, but it is **not**
  used to move the window. Window adequacy is checked directly, below.

**Regression:** unweighted OLS of ln D_m (double) and of ln D*_oracle (exact) on ln s, over all W_asym points.
- **Validity:** ≥ 8 points, ln-span ≥ ln 100, no removal except by the SNR rule.
- **Window adequacy** (exact data; otherwise DESIGN INSUFFICIENT for the pair, not a W/Z result): |p̂* − q| ≤ 0.02, and every
  exact centred local slope in W_asym is within 0.02 of q. A pair that fails is "insufficiently asymptotic". This is
  distinguished from a reference failure, which can only arise through §5.2 (ii)/(iii) or §4.7.
- **H2:** |p̂ − q| ≤ 0.02 and |p̂ − p̂*| ≤ 0.005.
- **Uncertainty.** Residuals are deterministic, so no probabilistic CI is claimed per pair. Reported per pair:
  - the deterministic interval p̂ ± Σ_i |w_i|·|ln(D_i/D*_i)| (w_i = OLS weights);
  - the OLS standard error, descriptive only.
- Across configurations: median, min and max per mutant and class, and a seeded bootstrap (10⁴ resamples) of the median.
- **Outliers:** none removed. Failed validity checks → INVALID, never dropped.

### 5.4 Local slopes, crossover and regimes (RQ1)

- Centred local slope: σ(s_k) = ln(D(s_{k+1})/D(s_{k−1})) / ln(s_{k+1}/s_{k−1}), exact and double. One-sided σ(1) uses
  s = 1 and 0.5.
- **H1:** |σ_double − σ_exact| ≤ 0.01 wherever both neighbours have SNR* ≥ 1e3, including s ∈ [2⁻⁴, 1].
- **Crossover** (exact):
  - r(s) = D*_oracle,R_W(s) / (‖C_q‖ s^q) − 1;
  - s_× = the largest grid s with |r(s')| ≤ 0.1 for all grid s' ≤ s;
  - secondary: the largest s with |σ_exact − q| ≤ 0.05 below it.
- **Regime models** (descriptive; fixed): ln D* = a + q ln s + ln|1 + Σ_{k=1}^m b_k s^k| for m = 0, 1, 2.
  - Select the smallest m with maximum |log residual| ≤ 0.01.
  - Otherwise label "strongly nonlinear in range".
  - No decision depends on the model choice.
- **RQ1 answer rule:** if H1 holds, any departure of σ from q near s = 1 is algebraic (a crossover), because the exact
  profile contains no rounding. If H1 fails, it is at least partly numerical.

### 5.5 Controls

- **class_D_order and M12:** the double D_m ≤ 10·N_X at every s (H4).
- **M9:** q reported (expected ≥ 1 in all classes, including D). This demonstrates that "D ∝ s" is not specific to ordering.
- **Structural mutants:** detected at every s with SNR* ≥ 10 (H4).
- **Noise-trend control:** for the correct code, the OLS slope of ln max(η_X, ε‖X_ref‖_F) on ln s over W_asym must lie in
  [−0.25, 0.25]. Otherwise the noise model is inadequate → DESIGN INSUFFICIENT.

## 6. Preregistered decision rules

Decisions use **set C (held-out) only**. Set X is reported and can never yield PASS.

Two outcomes are issued:
- **Physical-scope decision:** class R, with the class-D controls.
- **Algebraic stress-test outcome:** class H, with the same rules applied to its pairs.

The class-H outcome **never enters** the physical-scope decision and must not be cited as validation within the documented
scope of `wz_rgf.py`. INVALID conditions apply globally, because the apparatus is shared.

Precedence: **INVALID > FAIL-W > FAIL-D > DESIGN INSUFFICIENT > PASS**.

**INVALID** (global): any of
- a §4.7 control failure;
- a three-way mutant disagreement (§4.3);
- a missing R0 sign-off;
- a manifest hash mismatch;
- a §5.2 consistency contradiction (C_0 nonzero where a theorem requires zero; series vs direct disagreement; class D or
  M12 coefficients nonzero);
- a determinism mismatch.

**FAIL-W** (W/Z accuracy; per scope): any correct-code A endpoint, prefix, W_j or P_j error, or ΔS, above **1e-10** in a
judgeable case of that scope's classes.
- For the physical scope this is class R (and class D, where the correct code also runs).
- This is the only scientific W/Z failure criterion, and it is unchanged.

**FAIL-D** (a diagnostic or scaling hypothesis falsified; **not** evidence of a W/Z defect): any of
- H1 violated for a meaningful pair;
- H2 violated for a regression-valid, window-adequate meaningful pair;
- H4 violated (a class-D order mutant or M12 above 10·N_X, or a structural mutant undetected where SNR* ≥ 10).

**DESIGN INSUFFICIENT** (per scope): any of
- an anomaly (a heuristic bound exceeded while ≤ 1e-10);
- a non-judgeable case;
- an order-unresolved meaningful-candidate pair;
- a meaningful pair that is not regression-valid or not window-adequate;
- a failed noise-trend control;
- coverage below requirement: ordering-meaningful pairs fewer than 90 % of the expected count (R: 312), or fewer than 10 per
  pure order mutant;
- a C_0 that is unresolved for precision reasons (§5.2 (i)).

**PASS** (per scope): none of the above.

Physical-scope PASS claim:

> For 48 independently generated held-out class-R devices (N_y ≤ 4, N ≤ 6, equal lead widths, clean leads), conditional on
> the §2.4 acceptance filter, the correct W/Z recursion meets the 1e-10 criterion and its heuristic error estimates at
> 81 coupling values, and every meaningful ordering mutant's discrepancy follows its exactly identified asymptotic order
> down to the resolvable noise floor, with class-D controls and the M12 negative control behaving as required.

Class-H PASS claim: the same statement for class H, as an **algebraic stress test** only.

Interpretation:
- A pair with SNR* < 10 at some s carries **no** information about harm there.
- Class-D, M12, identical and ordering-zero non-detections are not detections.
- FAIL-D says nothing about W/Z correctness.

## 7. Scope

**Can establish** (conditional on §2.4):
- for small coupled devices, whether a strong-coupling crossover explains slope departures near s = 1;
- whether the exactly identified asymptotic order emerges in double precision;
- whether ordering errors are separable from rounding noise down to a quantified floor;
- whether the correct W/Z code meets the 1e-10 criterion across a dense coupling range.

Class R results refer to the documented physical scope. Class H results are algebraic stress tests only. Class D results
are controls.

**Cannot establish, and will not claim:**
- general W/Z validation;
- long-device stability (N ≤ 6);
- near-threshold behaviour (|δ_n| < 0.02 excluded);
- poorly conditioned, resonant or near-singular devices (excluded by the filter, §2.5);
- unequal lead widths;
- disordered systems (no disorder ensembles);
- anything about the B1 configuration F failure;
- any change to the previous experiment's DESIGN INSUFFICIENT decision;
- validation within the physical scope from class-H data.

## 8. Reproducibility and protocol

### 8.1 Planned files (none exist yet)

All under `experiments/wz_rgf_diagnostic/`:
- `coupling_scaling_v2_manifest.py` (generator, filter);
- `coupling_scaling_v2_certify.py` (validated exact-residual bounds, λ enclosures and the interval covering of §2.4.1, with
  self-tests T1–T4; also used for the validated reference bounds of §4.7);
- `coupling_scaling_v2.py` (driver);
- `wz_mutants_oracle.py` (§4.2 oracle and the truncated-series evaluation);
- `wz_mutants_mp.py` (clean-room transcription);
- `REVIEW_coupling_scaling_v2.md` (R0 sign-off);
- new directory `results_coupling_scaling_v2/`:
  - `manifest_X.json`, `manifest_C.json`: per configuration, the set, k, attempt, N_y, class, N, E, n_open, the integer arrays
    R_j and I_j, V_j(1) as exact rationals, the rejection log, the certificate (points, radii, the continuous cond bound,
    r_0), and the versions;
  - `manifest.sha256`;
  - `results_X.json`, `results_C.json`;
  - `points_X.csv`, `points_C.csv`: per (configuration, s), every A/B accuracy metric, conditioning, κ, the reference
    checks, and per mutant D_m, D*_oracle, D*_transcription for all four blocks (R_W and L), η_X, N_X, SNR*, and the category;
  - `pairs_X.csv`, `pairs_C.csv`: per (configuration, mutant), category, C_0..C_12 norms at dps 100 and 200, the zero-test
    results, q, ρ_est, p̂, p̂*, the deterministic interval, OLS SE, number of points, span, window adequacy, s_× (both),
    the regime model, σ(1) exact and double, and the H1/H2/H4 flags;
  - `coupling_scaling_v2_report.md`, with separate sections for class R, class H and class D;
  - `previous_files.sha256`.

### 8.2 Phases

1. This design is committed (the only action at this step).
2. **Implementation and manifest freeze:**
   - Write the code; generate the manifests (filter and certificate only); commit them with `manifest.sha256` **before any
     solver or reference run**.
   - Certificate self-tests T1–T4 must pass **before** the manifest is generated.
   - Complete review gate R0, which includes the certificate code, and commit its sign-off.
3. **Exploratory run (set X):**
   - Publish the report.
   - At most one amendment may follow, before the confirmatory run, limited to bug fixes, feasibility and added diagnostics.
     **No threshold, window, grid, SNR level, tolerance, coverage, K_max or decision rule may change**; any such change
     requires a new experiment ID.
   - An amendment touching the oracle or the transcription requires a new R0 sign-off.
4. **Confirmatory run (set C)**, then both decisions by §6.

### 8.3 Preserving the previous experiment

- Before any run, record the sha256 of every existing file under `experiments/wz_rgf_diagnostic/` (excluding
  `__pycache__`) in `previous_files.sha256`. Verify it after every phase; any change → INVALID.
- Code imports `wz_rgf.py` and `wz_mutants.py` read-only and writes only to `results_coupling_scaling_v2/`.

### 8.4 Independent reproduction

1. Fresh clone at the recorded commit.
2. Set `OPENBLAS_NUM_THREADS=1`.
3. Verify `manifest.sha256` (the frozen manifest is authoritative; numpy RNG streams may change across versions).
4. Run the driver.
5. Compare numeric outputs, excluding timestamps and paths: bitwise on the same platform and build. On another build,
   double values may differ at rounding level, and both decisions must be unchanged.
6. Record python, numpy, scipy and mpmath versions and the BLAS configuration.

**Cost** (uncertain; to be measured in the exploratory phase):
- 180 configurations × 81 points of mp references (≤ 24 × 24);
- oracle and transcription evaluations;
- series coefficients at two precisions per pair;
- certificate bisection.
- Estimated a few hours single-threaded. Parallelization is allowed only over configurations, with fixed output ordering.

## 9. Assumptions and remaining mathematical ambiguities

Assumptions:
- **A1 (heuristic):** ρ ≤ 25 for the window edge. Checked directly by window adequacy; ρ_est is reported.
- **A2 (heuristic):** |C_1| ≤ 1, used only to size the grid's lower edge.
- **A3 (heuristic):** the inherited τ bounds. Their exceedance is an anomaly, not FAIL.

Remaining ambiguities and open points:
1. **Zero-test thresholds** (1e-60, 1e-80, 1e-170, 1e-20) are design constants. A genuine coefficient smaller than 1e-60
   relative would be "unresolved" (→ DESIGN INSUFFICIENT), not misclassified. No theorem fixes these numbers.
2. **K_max = 12** is a finite search depth. A genuine order above 12 yields DESIGN INSUFFICIENT, not a claim.
3. **The certificate is on L(s) only.** R_W(s) is a rounded, non-affine perturbation of L(s). Its invertibility between grid
   points is irrelevant, because double data exist only at grid points. Its closeness to L is checked by §4.7 (2) and (3),
   but no continuous statement about R_W is made.
   - L(s) itself is defined with the double values of E and E_perp,n taken as exact. The certificate is a statement about
     that lattice, not about the mathematically exact E_perp,n = (nπ)²/2.
4. **The continuous bounds** (factor 2 from the half-radius covering) on invertibility, ‖A(s)⁻¹‖₂ and cond₂ are
   **validated** (§2.4.1), conditional on implementation correctness:
   - Python integer arithmetic and `math.isqrt` are assumed correct;
   - the code is covered by self-tests T1–T4 and review gate R0.
   A bug that the tests and review miss would invalidate the certificate. This is the only non-mathematical assumption of
   the certificate. The heuristic τ bounds evaluated at grid points are not extended continuously.
4a. **Non-validated numerical steps that remain:**
   - the zero test of Taylor coefficients (precision scaling);
   - the mp precision allowance in §5.2 (ii);
   - the dps-100 and `lu_solve` cross-checks;
   - the empirical grid conditioning filter;
   - every heuristic τ bound.
   None of them is called rigorous.
4b. **Validated statements available:**
   - **Nonzero coefficients:** a coefficient C_k could be proven nonzero by a validated enclosure (exact dyadic X0 plus the
     validated bound on ‖A0⁻¹ − X0‖₂ propagated through the products). This is **not** required and has no decision role;
     it may be recorded.
   - **Zero coefficients:** only C_0 = 0 for pure order mutants and M9 is established by theorem (F2, F6). Vanishing of
     C_1, …, C_{q−1} for q ≥ 2 rests on the numerical zero test alone.
5. **M10's matrix** is only checked pointwise. M10 is structural, outside the scaling decision, and its analyticity is not
   claimed.
6. **The oracle derivations** (§4.2) depend on the identity between recursion intermediates and prefix-reference
   quantities. That identity is verified numerically (≤ 1e-50) and algebraically (Schur complements), but the per-mutant
   algebra still requires the R0 review. An error in the code-change column itself is caught only by review.
7. **Strictness choices:**
   - coverage ≥ 90 %;
   - 100 % of meaningful pairs must be regression-valid and window-adequate;
   - anomaly → DESIGN INSUFFICIENT.
   These are deliberate preregistered choices and are not derived.
8. **v0 = 8** positions the crossover. A different strength would move s_× but not q.
9. **Runtime** may force a coarser dps-100 schedule. That would be a permitted feasibility amendment only if no decision
   threshold changes.
10. **Determinism across BLAS builds** is not guaranteed. Only same-build bitwise identity is required.
