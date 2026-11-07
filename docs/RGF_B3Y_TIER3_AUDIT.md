# B3Y Tier 3: mathematical and implementation audit

Date: 2026-11-07. Base: `origin/main` = `fdbff81ee891f8194b5923703f99c01efc10bd69` (merge of PR #7). Branch: `claude/b3y-tier3-audit`.
Files audited: `sparse_exact_lattice.py`, `tests/test_rgf_b3y.py`, `docs/RGF_B3Y_VALIDATION.md`, `docs/rgf_b3y_protocol.json`.

**Scope.** This is an audit of the Tier 3 infrastructure. **No B3Y numerical comparison was run, no B3Y result file was generated, and nothing was run on the real B3Y configurations.**
Every number below comes from synthetic problems (clean guide, y-uniform smooth barrier, a single-channel chain, small dense toys) at energies that are not B3Y points.
No protected solver, historical result, or predeclared threshold was changed (§5).
B3Y is a numerical-reliability campaign; it is not the project's research question (the effect of spatially correlated static disorder and of its strength on quantum wavepacket transport
in a two-dimensional hard-wall waveguide) and it makes no claim about Anderson localization, decoherence, or devices.

## 1. Summary of changes

| # | Finding | Change | File |
|---|---|---|---|
| F1 | `representation_dominance` tested `bulk ≥ 2·max(lead, solver)`; this does not exclude `lead + solver = bulk`, where the two other terms can cancel the bulk term completely | now `bulk ≥ 2·(lead + solver)` (dominance over the **sum**); strictly stronger; `term_resolution_factor` = 10 unchanged; unresolved terms stay UNEXPLAINED | `sparse_exact_lattice.py`, protocol JSON/doc |
| F2 | `error_decomposition` checked only the **two-term** identity (representation + solver); the three-term identity was never computed or checked | all three terms are computed explicitly; separate `identity_residual_three_term` and `identity_residual_two_term`; without `S_hybrid` the bulk and lead terms are `None`; documented as algebraic identities that do not validate the individual S-matrices | `sparse_exact_lattice.py` |
| F3 | Six extended-precision entry points had **no platform gate** (`apply_matrix_extended`, `build_sources`, `build_unit_sources`, `S_from_state`, `S_from_green_columns`, `iterative_refinement`); four of them ran with a simulated 53-bit longdouble and no refusal | gate added to all; registry `GATED_ENTRY_POINTS`/`EXEMPT_ENTRY_POINTS` with a test that every public callable is in exactly one and every gated one raises `NOT EVALUABLE` | `sparse_exact_lattice.py`, tests |
| F4 | The gate trusted `np.finfo` only | added a functional arithmetic probe `(1 + 2⁻⁶⁰) − 1 == 2⁻⁶⁰` (fails in float64); both must hold | `sparse_exact_lattice.py` |
| F5 | The resolution test used one e_ind (the exact solve's); the hybrid and rounded-matrix solves, which also enter the terms, had no indicator | `solve_matrix_extended` returns κ̂₁ and an error indicator; `combined_error_indicator` takes the maximum (NaN ⇒ unresolved) | `sparse_exact_lattice.py` |
| F6 | `e_ind` was reported as 0.0 when neither of its parts was available | now NaN | `sparse_exact_lattice.py` |
| F7 | The `x > 2` closed-channel branch and the retarded-root choice were not independently tested | tested against mpmath and against the E + iη root | tests |
| F8 | The hybrid (corner-swap) path was only exercised through a tautological identity | a non-tautological test attributes a known bulk perturbation to the bulk term | tests |

No formula was changed to make a test pass. One test tolerance (10⁻²⁰ for the retarded-root check, below the arithmetic's own resolution ε = 1.1×10⁻¹⁹) was a test bug and was corrected to 10⁻¹⁸.

## 2. Error decomposition

### 2.1 Definitions (max-abs norms)
* `S_exact`: solution of the **unmerged** operator in extended precision with exact (extended) lead data.
* `S_hybrid`: the method's **rounded bulk** operator with the lead terms replaced by the exact extended ones, solved to extended precision (the sources and the extraction use the exact lead data).
* `S_fl`: the method's **fully rounded matrix** (its double lead data in the matrix and in the extraction), solved to extended precision.
* `S_double`: the method's own double-precision result.

### 2.2 Identity
S_double − S_exact = (S_hybrid − S_exact) + (S_fl − S_hybrid) + (S_double − S_fl) = **bulk-operator rounding + lead data + linear-solver/extraction rounding**.
Without `S_hybrid` only S_double − S_exact = (S_fl − S_exact) + (S_double − S_fl) (representation + solver) exists.

**This is an algebraic telescoping sum.** Its residual measures only the floating-point evaluation of the differences. It does **not** validate the individual S-matrices, and the norms of the terms are
measured sizes obeying the triangle inequality, **not rigorous error bounds**. What *is* validated non-tautologically: the exact solution against mpmath (≤ 10⁻¹⁷ on the synthetic chains and ≤ 10⁻¹³ on the
dense toys), and the hybrid machinery (a known perturbation of every bulk diagonal entry is reproduced to 10⁻¹² by an independently built exact operator and attributed to the bulk term).

### 2.3 What each term contains
* **bulk-operator rounding**: rounding of the merged diagonal (E − 2t) − E⊥ − V (ulp of 2t: 5.8×10⁻¹¹ at dx = 1/512, 9.3×10⁻¹⁰ at 1/2048), the double evaluation of the V blocks, and the rounding of the
  merged corner entries (the hybrid subtracts the method's rounded t·Λ_d exactly, so the rounding of the sum base + t·Λ_d stays in this term).
* **lead data**: double-precision λ, Γ and incident amplitude in the matrix corners and in the source/extraction (measured Γ relative error 1.3×10⁻¹¹ at dx = 1/2048 for E = 88.5, channel 4).
* **solver/extraction**: double LU error, and the double-precision arithmetic of the source construction and the S extraction.

### 2.4 Which terms are available
| Method | Explicit matrix | Terms available | Claim allowed |
|---|---|---|---|
| M1 RGF | no (block recursion) | total error S_M1 − S_exact only | no mechanism claim; at most PRECISION-CONSISTENT (operational) |
| M2 Green-block | yes (merged double matrix, `smatrix_spec` lead data, Fisher–Lee on unit sources) | all three, if the driver solves the method's own matrix with unit sources and applies its extraction | representation dominance if supported by §2.5 |
| M3 wave-function | yes (merged double matrix, own lead data) | all three | same |
The B3Y driver (not yet written) needs each method's own rounded matrix, its double lead data (λ_d) and its extraction; this is a design requirement, not yet verified for M2 and M3.

### 2.5 Representation-dominance criterion (revision 2)
With b = bulk, l = lead, s = solver, e_ind = the largest convergence/error indicator among the extended solves used:
* **resolved**: b, l and s are each > 10 × e_ind (strict; NaN ⇒ unresolved);
* **dominant**: **b ≥ 2 × (l + s)**.
Both are required; otherwise the result is UNEXPLAINED (when a mechanism claim is attempted).

*Why the sum.* Because the identity is exact and the norm obeys the triangle inequality, |‖total‖ − b| ≤ l + s. Under b ≥ 2(l + s) this gives b/2 ≤ ‖total‖ ≤ 3b/2: the bulk term alone fixes the size of the total
error to within 50 %, **whatever the signs or cancellations** of the other two terms, and b is at least 2/3 of the sum of the three term norms. The old test b ≥ 2·max(l, s) allows l + s = b. A concrete case: b = 2, l = 1, s = 1
satisfies it, yet S_double = S_exact (total error 0): the "dominant" bulk term is completely cancelled (`test_dominance_case2…`). A property test over 200 random matrices confirms the 50 % statement under the new rule.
The criterion still concerns norms of measured terms; it does not say the bulk term *causes* the error beyond the definitions of §2.1.

## 3. Platform precision findings
* **Defect (F3).** Under the 53-bit simulation used by the tests, `lead_data_extended`, `potential_blocks_extended`, `UnmergedOperator`, `exact_lattice_S_from_blocks` and `solve_matrix_extended` refused, but `apply_matrix_extended`, `build_sources`,
  `build_unit_sources` and `iterative_refinement` ran without refusing, and `S_from_state`/`S_from_green_columns` had no gate (my dummy arguments failed for an unrelated reason). That is a silent-fallback path. It is closed.
* **Now.** 14 gated entry points (list in `GATED_ENTRY_POINTS`) and 17 exempt callables, each exempt one with a stated reason (mpmath-only references, pure classification/post-processing, the gate itself, exact power-of-two scalars, the double-precision preconditioner by design).
  A test asserts that these sets cover every public callable of the module, so a new function cannot be forgotten, and that **every gated entry point raises `ExtendedPrecisionUnavailable` whose message contains `NOT EVALUABLE`** under a simulated 53-bit longdouble.
* **Gate.** `np.longdouble` must provide ≥ 64 mantissa bits **and** pass a functional probe. The probe test shows that a platform whose `finfo` reports extended precision but whose arithmetic is double is refused. No fallback to float64 exists; none was added.
* **Reproducibility.** `platform_precision()` is deterministic (two calls give the same JSON-serialisable dict) and every `exact_lattice_S*` result carries it (`platform`). Measured here: x86-64 Linux, 64 mantissa bits, ε = 1.0842×10⁻¹⁹, 16-byte storage, Python 3.13.18, NumPy 2.5.3.
* **Platform dependence.** x86-64 Linux: 80-bit x87 extended. aarch64 Linux: IEEE quad in software (113 bits; accepted, tested by simulation only). Windows/MSVC and Apple arm64: longdouble = float64 → Tier 3 is NOT EVALUABLE there. Results from different extended types are not expected to be bitwise identical.

## 4. Mathematical audit (item by item)
Notation: slices j = 1…N, channels n = 1…N_y, t = 1/(2dx²), x_n = (E − E⊥ₙ)dx², Λ = diag(λₙ), P_j = diag(E − E⊥) − V_j.

| # | Item and expression implemented | Checks (signs, factors, indexing, dimensions, boundary terms) | Verification | Shared assumptions / independence |
|---|---|---|---|---|
| 1 | `UnmergedOperator.apply`: (Ac)_j = P_j c_j + t(c_{j−1} − 2c_j + c_{j+1}), ghosts c₀ = Λc₁, c_{N+1} = Λc_N | From M = E − H − Σ with H_jj = 2t + E⊥ + V, H_{j,j±1} = −t, Σ = −tΛ: (Mc)_j = (E−E⊥−V−2t)c_j + t c_{j−1} + t c_{j+1}, and −Σc = +tΛc = t·ghost ✓. `X[j]` ↔ slice j+1; arrays (N, N_y, nrhs); P real symmetric, applied to Re and Im separately ✓ | entry-by-entry mpmath M applied to a random vector: agreement ≤ 10⁻¹⁶·t (`test_unmerged_operator_equals_an_independently_built_mpmath_operator`); the unmerged operator differs from the merged double matrix by ≈ the ulp of 2t (test) | the operator is the specification §2.7.1/2.7.4 lattice problem: shared by M1–M4 (**mathematical** identity); M4's **implementation** is separate (extended, unmerged) |
| 2 | Lead data: open λ = (1−x) + i√(x(2−x)), Γ = 2t√(x(2−x)); x < 0: λ = 1/(c + √(y(2+y))), y = −x, c = 1−x; x > 2: λ = −1/((x−1) + √(x(x−2))) | 1−x never cancels (no c = 1 − x in double); Γ = v/dx ✓; retarded root has Im λ > 0 and |λ| < 1 (decaying) ✓; threshold rejection ✓ | vs mpmath (50 digits, **different formulas**: acos and polynomial roots): relative ≤ 1.1×10⁻¹⁹ for open and closed; the x > 2 branch now tested; the retarded choice equals the root with |λ| < 1 of the E + iη problem | M1/M2: `smatrix_spec.lead_modes` in double (acos(1−x): Γ error 8.9×10⁻¹² at dx = 1/2048): **same math, different and less accurate implementation**; M3: own closed form (double) |
| 3 | Source: a = Γ^{−1/2}, s = −t(1/λ − λ)·a e_n at slice 1 (left) or N (right); open channel: s = iΓa | Derived by splitting c₀ = Λ(c₁ − a e_n) + λ⁻¹ a e_n; sign: (Ac)₁ = −t(1/λ − λ)a e_n ✓; flux J = 2t Im(c_j* c_{j+1}) = |a|²Γ = 1 ✓ | flux = 1 to 10⁻¹⁷ (test); a flipped sign gives unitarity violation > 10⁻³ (test); S unitary to 10⁻¹⁵ on all gates | derivation shared by M3 and M4 (**math**); M1/M2 use Fisher–Lee; the two formulations are compared inside M4 (X5) |
| 4 | `S_from_state`: r = √Γ_m (c₁[m] − a δ_mn), t_RL = √Γ_m c_N[m]; right incidence mirrored; S = [[r_L, t_LR],[t_RL, r_R]] | outgoing amplitude b_m carries flux |b_m|²Γ_m ✓; only open rows/columns; `c1[o][:, :no]` ordering ascending ✓; dims (2n_o × 2n_o) ✓ | mpmath single-channel transfer matrix (r and t), clean guide (λⁿ⁻¹), mpmath dense solve; unitarity, reciprocity ≤ 10⁻¹⁵ | convention P, block order: shared by M1–M4 |
| 5 | `S_from_green_columns`: t_RL = iΓ_R^{1/2} G_N1 Γ_L^{1/2}, r_L = iΓ^{1/2}G₁₁Γ^{1/2} − I, … with G = A⁻¹e from separate unit-source solves | orientation: G_N1 = rows (slice N, open) × columns (unit source at slice 1) ✓; the −I term ✓ | X5: ≤ 5.9×10⁻¹⁷ vs `S_from_state` (separate solve) | spec §2.7.5 convention P; shared with M1/M2 (they implement it in double) |
| 6 | Refinement: r = b − A x (extended, unmerged), d = LU⁻¹ r (double), x ← x + d; η = ‖r‖∞/(‖b‖∞ + ‖M̃‖₁‖x‖∞); stop at η ≤ 10⁻¹⁷ after ≥ 2 iterations; contraction counted only while the previous η > 10⁻¹⁷ | η is the standard normwise relative backward error with the rounded matrix's 1-norm as proxy for ‖A‖₁; iterates stored in extended precision ✓; a small residual of the **rounded** matrix is not used as evidence ✓ | converges in 2 steps (contraction ≤ 1.6×10⁻⁶, η ≤ 4.9×10⁻²⁰) on all synthetic solves; a poor preconditioner (contraction 0.9) and a divergent one are flagged (test) | the LU is a preconditioner only; **implementation independent** of the other methods' solvers, **mathematically** the same linear system |
| 7 | Condition: κ̂₁ = ‖M̃‖₁·onenormest(M̃⁻¹) (block 1-norm estimator, t = 4, adjoint solves) of the merged double matrix | M̃ is complex symmetric ⇒ ‖·‖₁ = ‖·‖∞ and κ₁ = κ_∞ (tested), consistent with the ∞-norm η; the estimator is a lower bound, typically within a small factor | vs the exact `np.linalg.cond` on a small matrix (lower bound, within factor 4); e_ind = max(iterate difference, κ̂₁η) is a first-order **indicator, not a bound** | κ̂₁ is of the rounded proxy, not of the exact A; real B3Y devices (≈ 10⁷ extrapolated) untested |
| 8 | Exact lattice vs rounded matrix | S_exact is the problem defined by the **double** E, E⊥, V samples (exact binary rationals), dx = 2⁻ᵖ and exact irrational functions; fl(A) is a method-specific rounded matrix | decomposition tests vs mpmath; E and E⊥ bitwise identity is a B3Y criterion (P1) | S_exact is neither the continuum nor the "ideal-real" E⊥ problem; rounding of E⊥ and V relative to ideal reals is identical in all methods |
| 9 | Potential blocks: V_nm = (g_{|n−m|} − g_{n+m})/(N_q+1), g_k = Σ_l v_l cos(πkl/(N_q+1)) ⇔ ∫χₙVχₘ with χₙ = √(2/L) sin(nπy/L), y_l = lL/(N_q+1), weight L/(N_q+1) | L cancels; cos argument reduced by exact integer arithmetic (kl mod 2(N_q+1)); needs N_q ≥ 2N_y (enforced); the caller must supply samples at y_l = l/(N_q+1) | vs direct mpmath sine sums (10⁻¹⁷·|V|); constant potential → c·I; the analytic ∫ y χₙχₘ dy reproduced to 2×10⁻⁵ at N_q = 255 (a wrong point convention would give O(1)) | the quadrature definition and sample positions are shared by all methods (P1 verifies identical samples) |

**Mathematical vs implementation independence.** M1–M4 solve the *same mathematical problem* (§2.7 lattice with exact clean-lead elimination, convention P, same potential samples and E, E⊥ doubles). Independence is therefore of
*implementation and algorithm*: M1 (recursion) vs M2/M3 (SuperLU) vs M4 (unmerged operator, extended precision). Mathematical independence exists only in the S extraction (Fisher–Lee in M1/M2 vs scattering state in M3/M4; M4 computes both) and in the
lead formulas (acos-based double vs cancellation-free extended). An error in the shared lattice definition is not detectable by B3Y.

## 5. Tests and hash verification
Local runs (no CI).

| Command | Exit code | Result |
|---|---|---|
| `python -m pytest -q tests/test_rgf_b3y.py` | 0 | 50 passed, 0 failed, 0 xfailed (≈ 54 s) |
| `python -m pytest -q tests` | 0 | 302 passed, 0 failed, 2 xfailed (the two known expected failures); 4 warnings |

Warnings (all pre-existing, none from the B3Y tests): `tests/test_paired.py` NumPy "invalid value encountered in divide" (×2) and SciPy "Precision loss … catastrophic cancellation" (×2).
Before this audit the B3Y file had 34 tests and the full suite 286. The audit added 17 tests and replaced one (the old two-term identity test became the three-term test), net +16: the dominance boundary cases (3), the three-term identity (replacement) and the no-hybrid two-term case (1 new), the hybrid attribution test (1), the registry/gate/probe/metadata tests (3), and the mathematical-audit tests (8 + the rewritten classification test).

**Protected-file and result-tree hashes** (`docs/rgf_b3y_protocol.json`, recomputed by script and by `test_protected_solver_files_and_prior_artifacts_are_unchanged`):

| Check | Before the changes (at `fdbff81`) | After the changes |
|---|---|---|
| 18 protected files (solvers `stationary.py`, `smatrix_spec.py`, `sparse_reference.py`, `sparse_wavefunction.py`; validators B1–B3X; `potential.py`, `disorder.py`, `config.py`, `analytic.py`, `ode1d.py`; B1–B3X docs and protocol) | all match | all match |
| 4 historical result trees (B1, B2, B3, B3X; 39 files) | all match (`rgf_b1` `f4d84d54…`, `rgf_b2` `a7448587…`, `rgf_b3` `85c6116d…`, `rgf_b3x` `49fe3f51…`) | all match |
| 5 B3X artifact hashes | all match | all match |
| `git diff fdbff81` over these paths | — | empty |

Changed files: `sparse_exact_lattice.py`, `tests/test_rgf_b3y.py`, `docs/RGF_B3Y_VALIDATION.md`, `docs/rgf_b3y_protocol.json`, and this report. No B3Y result file exists.

## 6. Remaining limitations
1. All gates were exercised on synthetic problems only; the 63-channel disordered devices (κ̂₁ ≈ 10⁷ extrapolated), gate X8 and the B3Y-energy instances of X1–X3 are untested.
2. The terms of the decomposition are measured sizes, not bounds; e_ind is an indicator.
3. M1 (RGF) cannot be decomposed. M2/M3 decomposition needs each method's own rounded matrix, λ_d and extraction in the (unwritten) driver.
4. The bulk term includes the rounding of merged corner entries and of the double V blocks (§2.3); the lead term includes the double-precision extraction normalisation.
5. Extended precision is platform dependent (§3); results on other extended types may differ in the last digits.
6. The method independence is of implementation, not of the underlying lattice problem (§4).
7. One seed, one random draw per ζ (W = 200 is exactly twice W = 100); a finite set of grids cannot prove continuum convergence.

## 7. Issues requiring human approval
1. **Sum-form dominance** (`bulk ≥ 2·(lead + solver)`) in place of the max-form. It is stricter; the factor 2 is unchanged but remains one of the two flagged new numbers (`dominance_factor`, `term_resolution_factor` = 10).
2. **All three terms must be resolved** (> 10 × e_ind). A tiny (< 10 × e_ind) lead or solver term therefore makes the result UNEXPLAINED even if the bulk term is overwhelming (`test_dominance_case3…`). This follows the instruction to keep unresolved terms UNEXPLAINED; an alternative (use the unresolved terms' upper bounds l + e_ind, s + e_ind) would be less conservative and needs your decision.
3. **Combined indicator = maximum over the exact, hybrid and rounded-matrix solves**, with NaN making everything unresolved.
4. **Arithmetic probe as an additional gate condition** (strengthens the platform gate; cannot weaken it).
5. **Attribution of the corner-merge rounding to the bulk term** (§2.3) rather than to the lead term.
6. The exempt list (§3), in particular that `hopping` and the double-precision preconditioner are not gated.
7. The **new branch** `claude/b3y-tier3-audit` (as you requested) was pushed so the work is not lost; the earlier instruction set a different designated branch.
8. Still open from earlier reviews: A4 (buffer-control methods), A8, A9 (spread limit 2×10⁻⁶), A10 (gate tolerances 10⁻¹⁷/10⁻¹³).
9. **Authorization to write `validate_rgf_b3y.py` and to run any B3Y numerical comparison** is not given and was not assumed.
