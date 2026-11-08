# B3Y: independent fine-grid validation of the RGF (failed B3X configurations)

Phase 2, follow-up of [B3X](RGF_B3X_VALIDATION.md) (unchanged). Date: 2026-11-07.
Machine-readable protocol: [rgf_b3y_protocol.json](rgf_b3y_protocol.json).
Base commit: `5360a9f` (merge of PR #6).

**Status (reconciled after the campaign; four layers, kept separate):**

| layer | state | where |
|---|---|---|
| 1. Protocol definition (predeclared) | Frozen at revision 2 (commit `69c9539`); `docs/rgf_b3y_protocol.json` is **unchanged** and still carries its historical `status` text ("no B3Y numerical comparison has been run"), which was true when it was written. sha256 `489aea43f7d24aa913b59233d0f45ae78e4d347437ba0e405001f08a2ee8180a` (recorded in `run_metadata.json`; verified equal to the file). | Part I |
| 2. Implementation and synthetic testing | `sparse_exact_lattice.py` (M4, audited), `validate_rgf_b3y.py` (driver), `tests/test_rgf_b3y.py` and `tests/test_rgf_b3y_driver.py` (synthetic problems only). | §I.12, `RGF_B3Y_TIER3_AUDIT.md` |
| 3. Real B3Y numerical campaign | **Completed** (driver commit `1a5688b`, results commit `36055a7`); stopped_at = null. | `results/rgf_b3y_validation/`, Part II |
| 4. Campaign outcome | Overall **PASS** under the predeclared criteria; no exceedance of 10⁻¹⁰, nothing adjudicated; see Part II (II.0, II.8) and the post-campaign audit II.11. | Part II |

Statements in Part I, I-B and in `RGF_B3Y_TIER3_AUDIT.md` that say "no B3Y comparison has been run" describe the state **at the time they were written** and are marked *[historical]* where they occur. They are kept verbatim so the predeclared record is not rewritten.

Part I was written before any B3Y number existed; Part II (results) was appended after the runs. The only numerical work done before the campaign was the synthetic design-verification probe of §I.9
(a single-channel chain; it uses no B3Y point and no output of the RGF or of any B3X reference). *[historical]*

### Protocol revision 1 (clarifications; no existing criterion value changed)
Revision 0 is commit `5d39721`. Revision 1 clarifies five points and adds the Tier 3 implementation (synthetic tests only).
1. **Absolute S criterion preserved.** The criterion on S is max|ΔS| ≤ 10⁻¹⁰, absolute. ρ remains a diagnostic quantity and is never used to classify.
2. **The factor-10-per-halving rule is an operational precision-consistency classification, not a mathematical proof of round-off dominance.** It only says that an
   error that already met 10⁻¹⁰ at dx = 1/512 grew by at most 10× per halving; it says nothing about the mechanism.
3. **Representation dominance needs the decomposition.** Any claim that operator-representation error dominates must be supported by the Tier 3 error decomposition
   (§I.8.3 step 5). If the decomposition does not support it, or cannot be resolved, the result is classified **UNEXPLAINED**. M1 (the RGF) has no explicit matrix, so
   no mechanism claim is allowed for it.
4. **u_XP is renamed the convergence/error indicator e_ind.** It is a convergence indicator, not a rigorous uncertainty bound; its predeclared threshold 10⁻¹¹ is unchanged.
5. **Extended precision is mandatory.** `np.longdouble` must provide **at least 64 mantissa bits** (x86-64 80-bit extended). Where it does not (for example
   Windows/MSVC builds, or Apple arm64, where longdouble equals float64), Tier 3 stops and reports **NOT EVALUABLE**. Ordinary double precision is never substituted silently.

### Protocol revision 2 (Tier 3 audit; branch `claude/b3y-tier3-audit`; no numerical B3Y comparison)
Details and the full mathematical audit are in [RGF_B3Y_TIER3_AUDIT.md](RGF_B3Y_TIER3_AUDIT.md). Changes: (1) **representation dominance is now judged against the SUM** of the lead-data and solver terms,
`bulk_rounding ≥ 2 × (lead_data + solver)`, a strictly stronger criterion than the revision-1 `bulk_rounding ≥ 2 × max(lead_data, solver)`; `term_resolution_factor` = 10 and the rule that unresolved
terms are UNEXPLAINED are unchanged. (2) `error_decomposition` computes all three terms explicitly and reports a three-term and a two-term identity residual separately; these algebraic identities do
**not** validate the individual S-matrices. (3) Resolution is judged against the largest e_ind among the extended solves used. (4) The platform gate now covers **every** extended-precision entry point
(six had none) and adds a functional arithmetic probe. No threshold was weakened and no formula changed.

Two **new numbers** are introduced by item 3 and are flagged for review: `dominance_factor` = 2 and `term_resolution_factor` = 10 (§I.8.3). The B3Y energy-point instances of gates
X1–X3 and X8 remain for the later numerical stage; in this stage the gate code is exercised only on synthetic problems at energies that are not B3Y points.

## Part I: Protocol

### I.1 Scope, and what B3Y can and cannot establish
B3X found that the first-channel transmission T1 changes by up to 9.2×10⁻⁴ between dx = 1/512 and 1/1024 for the W = 200
configurations (limit 2×10⁻⁴). Its dx ladder (1/256 … 1/2048) used the RGF only, and the two sparse references were
compared with the RGF at dx = 1/512 only. Therefore **the dx dependence has not been shown to be solver-independent,
and it has not been shown to be a finite-difference discretization error.**

B3Y asks one question: *do three independent solvers, and a high-precision reference of the exact lattice problem,
reproduce the same complex S-matrix at dx = 1/1024 and 1/2048?*
* If yes, the dx dependence belongs to the lattice problem, not to a solver artifact (statement H1).
* B3Y cannot show that the lattice results converge to a continuum limit. That needs an independent discretization
  or a continuum reference. B3Y reports measured orders from three grids and flags every Richardson value as an estimate.
* All methods share the problem definition (§I.5), one random draw per ζ, and one seed.

### I.2 Authorization decisions and where they are implemented
| Decision | Implemented in |
|---|---|
| D5: keep absolute \|ΔS\| ≤ 10⁻¹⁰; ρ is a diagnostic only; an exceedance is precision-limited / inconclusive pending Tier 3, not an automatic failure | §I.8 (S1–S3, adjudication), §I.8.3 |
| Tier 3 separate high-precision exact-lattice reference, existing solver files untouched | §I.6, §I.7 |
| Tier 1: 15 W = 200 points + W100_z0.05 E = 371.36 control, dx 1/512 and 1/1024, 3 methods | §I.3, §I.8 |
| Tier 2: five failed points at dx = 1/2048, 3 methods + Tier 3 | §I.3 |
| Tier 3 at dx 1/1024 and 1/2048 for the five failed points | §I.3, §I.6 |
| Buffer control 205/512 vs 103/512 at dx = 1/1024, separate | §I.4, K1 |
| dx = 1/4096 not required | out of scope |
| Potential identity (direct call, coincident slices, hashes, clean-lead boundaries) | P1–P4 |
| Lead-tail guard 10⁻⁶ unchanged; failure = not evaluable | P4 |
| Diagnostics (condition numbers, residuals, lead-root and Γ errors, both unitarity products, reciprocity, T/R, off-diagonal T, channels) | §I.8, §I.10 |
| Convergence orders / Richardson only inside the predeclared window | O1 |
| Reproducibility record | §I.10 |

### I.3 Points (exact, from the B3X raw CSVs)
Seed 0, L_d = 1; energies E = (kx² + π²)/2 with kx = 20 + j·3.5355, j = −2…2 (4, 5, 6, 7, 8 open channels). The E values are
stored with full precision in the JSON; B3X artifact hashes are recorded there too. The **five failed points** are the
B3X rows of `refinement_convergence.csv` (name `dx`) with `passed = False`.

| Role | Configuration | E | Open | B3X \|ΔT1\| (dx 1/512→1/1024) |
|---|---|---|---|---|
| failed, Tier 2/3 | W200_z0.02 | 281.90 | 7 | 2.58×10⁻⁴ |
| failed, Tier 2/3 | W200_z0.05 | 140.47 | 5 | 5.97×10⁻⁴ |
| failed, Tier 2/3 | W200_z0.05 | 204.93 | 6 | 2.80×10⁻⁴ |
| failed, Tier 2/3 | W200_z0.05 | 371.36 | 8 | 2.80×10⁻⁴ |
| failed, Tier 2/3 | W200_z0.2 | 88.51 | 4 | 9.17×10⁻⁴ |
| control (passed in B3X, near the limit) | **W100**_z0.05 | 371.36 | 8 | 1.81×10⁻⁴ |
| context | the other 10 W = 200 energies (all three ζ) | — | 4–8 | < 2×10⁻⁴ |

The label "E = 371.4" occurs twice: W200_z0.05 (failed) and the W100_z0.05 control. Points are always identified by
configuration id **and** E.

| Tier | Points | dx | Methods |
|---|---|---|---|
| 1 | 16 (15 W = 200 + control) | 1/512, 1/1024 | M1 RGF, M2 Green-block, M3 wave-function |
| 2 | 5 failed | 1/2048 | M1, M2, M3, M4 (exact lattice) |
| 3 | 5 failed | 1/1024, 1/2048 | M4, compared with M1–M3 |
| 3 validation | all 16 | 1/512 | M4, compared with the B3X-era references M2, M3 (gate X8) |
| buffer control | 5 failed | 1/1024, buffer 205/512 | M1, M3 |

### I.4 Geometry, potential identity, lead tails
* N_y = 63, N_q = 127 (unchanged). The device is the physical interval [x_start, x_end] = [1.9 − 103/512, 2.9 + 103/512] at
  every dx (buffer 0.2012 on each side); slice counts 719, 1437, 2873. The lattices are nested: every 1/512 slice is also
  a 1/1024 and a 1/2048 slice. **No buffer change in the primary results.**
* **Buffer control (separate, K1):** 205/512 = 0.4004 at dx = 1/1024 (1845 slices): a physically longer device. The reference planes
  (first and last slice) move, so complex S differ by plane phases; only plane-independent quantities (|S_ij|², T1) are
  compared between the two buffers. Complex-S agreement of M1 vs M3 is evaluated within each geometry.
* **P1** bitwise identity of every method's sampled potential with a direct `Potential(cfg, disorder_from_config(cfg, 0)).numpy` call at each
  resolution; E and the E⊥ doubles bitwise identical across methods. **P2** bitwise equality at coincident physical slices across
  resolutions. **P3** sha256 of every potential and coordinate array recorded. **P4** |V| ≤ 10⁻⁶ on the first lead slices
  x_start − dx and x_end + dx of every device; a device above it is **NOT EVALUABLE** under the predeclared geometry
  (guard unchanged; geometry unchanged; no silent substitution).
* Pre-protocol facts from potential evaluation only: P2 holds for all three W = 200 configurations; the first-lead-slice tails are
  8.3, 9.1, 9.6 ×10⁻⁷ (W200_z0.02 at 1/512, 1/1024, 1/2048), 5.7–6.6×10⁻⁷ (z0.05) and 4.8–5.6×10⁻⁷ (z0.2). The smallest margin to the guard is 4.3 %.

### I.5 Methods; what is shared and what is independent
| Item | M1 RGF | M2 Green-block | M3 wave-function | M4 exact lattice (new) |
|---|---|---|---|---|
| Lattice Hamiltonian, DST-I definition, N_y truncation | shared | shared | shared | shared (the problem definition) |
| Potential callable / double samples | shared | shared | shared | shared |
| E and E⊥ doubles | shared | shared | shared | shared (used as exact inputs) |
| Lead roots, Γ | `smatrix_spec` | same function as M1 | own closed form (double) | own, cancellation-free, extended precision, checked against mpmath |
| V_j blocks | einsum | matrix product | cosine sums (double) | cosine sums, extended precision |
| Operator representation | blocks added per slice (double) | one merged sparse matrix (double) | one merged sparse matrix (double) | **unmerged**: (E − E⊥ − V) and t·(second difference) kept separate |
| Linear algebra | block recursion | SuperLU, COLAMD | SuperLU, MMD_AT_PLUS_A, 1 refinement step | SuperLU only as a preconditioner; **iterative refinement with extended-precision operator residuals** |
| S from | Fisher–Lee on G blocks | Fisher–Lee on G blocks (`smatrix_spec`) | scattering state | both: scattering state and Fisher–Lee, from the same solution |

### I.6 Tier 3: the exact-lattice reference (design; code to be written after approval)
New module `sparse_exact_lattice.py` (no edit of `sparse_wavefunction.py`, `sparse_reference.py`, `stationary.py`, `smatrix_spec.py`).

**Three objects are kept apart.**
1. *The exact lattice problem* A: the linear operator defined by the exact real numbers E, E⊥ₙ, the potential samples V(x_j, y_l) (all
   binary doubles, hence exact), dx = 2⁻ᵖ (so t = 1/(2dx²) is exact), the DST-I definition of V_j, and the lead elimination Σ = −tΛ.
   Only the irrational functions in the quadrature and the leads (sines, cosines, arcsine) are not exactly representable.
2. *Its floating-point representation* fl(A): the double-precision matrix a method assembles (rounded diagonal (E−2t)−E⊥−V,
   double V blocks, double Λ in Σ).
3. *The numerical solution*: the double-precision LU solve of fl(A)x = b.

**Unmerged operator.** With P_j = diag(E − E⊥) − V_j (extended precision) and second differences with lead ghosts c₀ = Λc₁, c_{N+1} = Λc_N,
the exact operator acts as (Ac)_j = P_j c_j + t(c_{j−1} − 2c_j + c_{j+1}). The large 2t is never added to the O(10²) physical part.
In double precision it is: the ulp of the diagonal entry is 5.8×10⁻¹¹ (dx = 1/512) and 9.3×10⁻¹⁰ (1/2048) in energy units.

**Extended precision.** `np.longdouble` (x86-64: 80-bit, 64-bit mantissa, ε = 1.08×10⁻¹⁹; measured 0.94 s per block matvec at N = 2873).
Precision is **platform dependent**: on x86-64 Linux it is 80-bit; on aarch64 Linux it is IEEE quad (software); on Apple/ARM and
Windows it equals double. **Gate X0 requires `np.finfo(np.longdouble)` to provide at least 64 mantissa bits and a functional arithmetic probe ((1 + 2⁻⁶⁰) − 1 = 2⁻⁶⁰); otherwise Tier 3 stops and reports NOT EVALUABLE on that platform**
(no silent fallback to double precision; the code raises `ExtendedPrecisionUnavailable`). The platform, `np.finfo(np.longdouble)` and the library versions are recorded.

**Lead data (cancellation-free).** x_n = (E − E⊥ₙ)dx² in extended precision. Open (x ∈ (0,2)): θ = 2·asin√(x/2), λ = e^{iθ}, Γ = 2t·sinθ,
sinθ = √(x(2−x)). Closed (x < 0): λ = 1/(c + √(c²−1)), c = 1 − x. This avoids the cancellation in c = 1 − x that gives Γ a relative
error of 1.3×10⁻¹¹ at dx = 1/2048 in the double formulas (measured against a 60-digit reference for E = 88.5, channel 4).
**Independent scalar check:** every channel's λ and Γ are recomputed with mpmath (50 digits) and compared (gate X1).

**Solution and refinement.** Preconditioner: SuperLU of the double matrix. Iterates x_k are stored in extended precision. Residual
r_k = b − A x_k with the *unmerged* operator in extended precision; correction δ = LU⁻¹ r_k; x_{k+1} = x_k + δ.
A small residual computed from the rounded matrix is **not** accepted as evidence (it measures only the solver). Convergence
must be demonstrated for every solve (gate X4): the relative operator residual η = ‖r‖∞/(‖b‖∞ + ‖A‖₁‖x‖∞) must reach ≤ 10⁻¹⁷ (≈ 100 ε_ld) within 10 iterations
(the loop stops at the first iterate k ≥ 1 with η ≤ 10⁻¹⁷, so an iterate difference exists), and the contraction ‖δ_k‖/‖δ_{k−1}‖ must be ≤ 0.1 for **every step taken
while the previous η still exceeded 10⁻¹⁷**. Below the target the correction stagnates at the evaluation-noise floor by construction, so contraction is not
meaningful there. (An earlier draft anchored the floor to a multiple of ε_ld·‖x‖; the implementation showed that noise floor is ≈ 150–200 ε_ld‖x‖, so the η-anchored definition is used.)

**Separating representation error from solver error.** With S(·) the S-matrix from a solution:
S_double − S_exact = [S(x_fl(A)) − S_exact] + [S_double − S(x_fl(A))], where x_fl(A) solves the rounded matrix **exactly** (to extended
precision, by refinement against the rounded matrix entries). The first bracket is representation error, the second is linear-solver
error. A further split of the first bracket (bulk-operator rounding vs lead data) uses a hybrid operator: rounded bulk diagonal, exact lead
terms. These identities need the explicit matrix, so they are available for M2 and M3. **M1 has no explicit matrix: for M1 only the total
error S_M1 − S_exact is measured.**

**Convergence / error indicator e_ind.** e_ind = max(max|S^(k) − S^(k−1)| over the last two iterates, κ̂₁·η_final), with κ̂₁ the 1-norm condition estimate of the
double matrix. It is an **indicator of convergence and of the size of the remaining error, not a rigorous uncertainty bound** (the bound κ̂₁·η is a first-order normwise estimate;
the iterate difference is a stagnation measure). Adjudication at a point requires e_ind ≤ 10⁻¹¹ (10× below the 10⁻¹⁰ target; threshold unchanged); otherwise the point is **unresolved**.

**Two S forms.** M4 extracts S from the scattering state and also by Fisher–Lee from the same Green columns; their difference (gate X5,
≤ 10⁻¹³) tests the equivalence of the two formulations free of round-off.

### I.7 Gates for M4 (run before M4 is used to interpret anything; failure stops the run order)
| Gate | Test | Tolerance | Why |
|---|---|---|---|
| X0 | `np.finfo(np.longdouble).nmant + 1 ≥ 64` | exact | extended precision available |
| X1 | λ, Γ of every channel/point/dx vs mpmath (50 digits) | relative 10⁻¹⁷ | ≈ 100 ε_ld; the probe reached 10⁻¹⁷ on t |
| X2 | clean guide (zero potential, N_y = 63): t_nn vs mpmath λⁿ⁻¹ at dx 1/512, 1/1024, 1/2048, at the five failed energies | 10⁻¹³ | probe error ≤ 2×10⁻¹⁷; 10⁻¹³ leaves room for 63 channels |
| X3 | y-uniform smooth barrier (the B2 barrier), channels decoupled: M4 vs mpmath single-channel transfer matrix for every open channel | 10⁻¹³ | same |
| X4 | refinement convergence: η ≤ 10⁻¹⁷ within 10 iterations; contraction ≤ 0.1 for every step taken while the previous η > 10⁻¹⁷ (§I.6) | as stated | §I.6 |
| X5 | scattering-state S vs Fisher–Lee S from the same solution | 10⁻¹³ | formulations equivalent; extended arithmetic |
| X6 | S stable under the preconditioner ordering (COLAMD vs MMD_AT_PLUS_A) | 10⁻¹³ | the result must not depend on the preconditioner |
| X7 | convergence/error indicator e_ind at every use | ≤ 10⁻¹¹ | §I.6 (an indicator, not a rigorous bound) |
| X8 | at dx = 1/512, all 16 points: max\|S_M4 − S_M2\|, max\|S_M4 − S_M3\|, max\|S_M4 − S_M1\| | 10⁻¹⁰ | validation against the existing B3X-era references before any finer-grid use |
| G0 | recomputed dx = 1/512 reproduces B3X: M1 T1 bitwise; M1–M3 max\|ΔS\| within 10⁻¹³ of the B3X value | as stated | the same code path gives the same numbers |

### I.8 Acceptance criteria and classification
**Tolerances are those of B3X (not modified):** 10⁻¹⁰ (S), 2×10⁻¹⁰ (T/R: |S_ij|² − |S′_ij|² ≤ 2 max|ΔS|), 10⁻¹⁰ (unitarity), 10⁻⁸ (reciprocity),
2×10⁻⁴ (buffer ΔT1), 10⁻⁶ (lead tail). New numbers are only the M4 gates (§I.7), the sparse backward error (10⁻¹⁴; B3X observed ≤ 8.8×10⁻¹⁶), the
method-spread limit of H1 (2×10⁻⁶ = 1 % of the 2×10⁻⁴ acceptance limit) and the growth factor of the adjudication rule (10).

| ID | Criterion | Limit | Applied to |
|---|---|---|---|
| P1–P4 | potential identity, cross-resolution identity, hashes, lead guard | bitwise / records / 10⁻⁶ | every device |
| N1 | open count equal in all methods and analytic; Γ_closed = 0; closed λ real, \|λ\| < 1; ≥ 1 closed channel | exact | all methods incl. M4 |
| S1 | max\|S_M1 − S_M2\| (complex, full S) | 10⁻¹⁰ | Tier 1 (512, 1024), Tier 2 (2048) |
| S2 | max\|S_M1 − S_M3\| | 10⁻¹⁰ | same |
| S3 | max\|S_M2 − S_M3\| | 10⁻¹⁰ | same |
| S4 | T, R, off-diagonal transmission agreement (M1 vs M2, M1 vs M3) | 2×10⁻¹⁰ | same |
| S5 | unitarity: max\|SᴴS − I\| and max\|SSᴴ − I\| reported separately | 10⁻¹⁰ | M1–M3 (M4 reported; gated by X5, X7) |
| S6 | reciprocity: full, blocks, probabilities | 10⁻⁸ | M1–M3 (M4 reported) |
| R1 | sparse backward error after refinement | 10⁻¹⁴ | M2, M3 |
| V1 | Tier 3 validated at dx = 1/512 (X0–X8) | gates | before M4 interprets finer grids |
| K1 | buffer control: \|ΔT1\| ≤ 2×10⁻⁴; complex-S agreement M1 vs M3 within each geometry ≤ 10⁻¹⁰ | as stated | 5 failed points, dx = 1/1024 |
| O1 | convergence orders (classification): T1 at 1/512, 1/1024, 1/2048, per method, 5 points; Richardson only if **both** orders ∈ [1.8, 2.2], always labelled an estimate | — | Tier 2 points |
| H1 | solver-independence statement for a point: S1–S4 hold at all its dx (or exceedances are adjudicated) and the T1 spread among methods at each dx ≤ 2×10⁻⁶ | statement | per point |

**Always recorded** (numerical diagnostics): condition estimates κ̂₁ of the sparse matrices, sparse backward errors, M4 operator-level residuals
and refinement histories, lead-root and Γ errors (double formulas vs M4/mpmath), S†S and SS† separately, reciprocity (all forms), T, R, off-diagonal
transmission, open/closed channel counts, closed-channel decay and amplitude content, and ρ.

#### I.8.1 Status vocabulary
**PASS**; **FAIL** (criterion not met and not attributable to precision limits, or an M4 gate failed); **PRECISION-LIMITED (inconclusive)** (S1–S3 absolute
discrepancy > 10⁻¹⁰); **NOT EVALUABLE** (lead guard exceeded, or `np.longdouble` has fewer than 64 mantissa bits); and, from Tier 3 adjudication, **PRECISION-CONSISTENT**
(operational), **PRECISION-LIMITED (measured; representation-dominated)** and **UNEXPLAINED** (§I.8.3).

#### I.8.2 The absolute 10⁻¹⁰ target (decision D5)
S1–S3 are judged by the absolute discrepancy only. A discrepancy > 10⁻¹⁰ is **not** automatically an implementation failure: it is
classified **PRECISION-LIMITED (inconclusive)** and stays so unless Tier 3 adjudicates it (§I.8.3). At the points without Tier 3 data (the 11 non-failed
points at 1/1024 and the 1/512 points other than the validation set) it stays inconclusive.

**ρ is a diagnostic and is never used to classify.** ρ = max|ΔS| / [(N−1)ε_mach/sin(k·dx)], with N the number of slices, ε_mach = 2.22×10⁻¹⁶, k the wave
number of the lowest-k open channel. *Why this is not a rigorous pairwise bound:* the B1 relation (|Δφ| ≤ 0.89·N·ε/sin(k dx)) was measured for one
method (the RGF) on a clean guide against an exact reference. For two independent implementations, |S_a − S_b| ≤ ε_a + ε_b with unrelated, partly systematic
(e.g. the common rounding of E − 2t) and partly random error components; the constant 0.89 has not been derived or tested for sparse-LU solves, for
disordered multichannel devices, or for differences between two methods. The B3X calibration (observed/bound at dx = 1/512: min 0.081, median 0.162, max 0.357) is used only
as a pre-run prior (§I.11).

#### I.8.3 Tier 3 adjudication (predeclared; revision 1)
Precondition: V1 passed and the convergence/error indicator e_ind ≤ 10⁻¹¹ at the point.
1. ε_m = max|S_m − S_M4| for M1, M2, M3 at the same point and dx.
2. ε_m ≤ 10⁻¹⁰ → m is consistent with the exact lattice solution at 10⁻¹⁰. ε_m > 10⁻¹⁰ → m **exceeds** 10⁻¹⁰ against the exact lattice (reported with its value; never upgraded to PASS).
3. For M2 and M3 the total error is split **exactly** (telescoping identity, §I.6) into a bulk-operator rounding term, a lead-data term and a linear-solver/extraction term (an algebraic identity: its residual tests the arithmetic, not the individual S-matrices; the terms are measured sizes, not rigorous bounds). For M1 only the
   total error is available, and **no mechanism claim is allowed for M1**.
4. **PRECISION-CONSISTENT** (an *operational* classification, **not a mathematical proof** of round-off dominance): ε_m(1/512) ≤ 10⁻¹⁰ and ε_m grows by at most a factor 10 per halving of dx.
   It makes no statement about the mechanism.
5. **Representation dominance** may be asserted only for M2/M3 and only if the decomposition supports it: the bulk-rounding term ≥ `dominance_factor` (2) × (lead-data term **+** solver term)
   (revision 2: dominance over the *sum*), **and** every term is resolved, i.e. > `term_resolution_factor` (10) × e_ind, with e_ind the largest indicator among the extended solves used.
   Meaning: by the exact identity and the triangle inequality, |‖total‖ − bulk| ≤ lead + solver ≤ bulk/2, so the bulk term alone fixes the size of the total error to within 50 % whatever the cancellations. A supported claim gives **PRECISION-LIMITED (measured; representation-dominated)**.
   If the claim is attempted and is not supported, or the terms are unresolved, the result is **UNEXPLAINED**.
6. A result that violates the growth rule of step 4 is **UNEXPLAINED** and is investigated before any conclusion.
7. **Overall:** B3Y is **PASS** only if every criterion is met with no exceedance. Any exceedance makes the overall status PRECISION-LIMITED / PRECISION-CONSISTENT (resolved by Tier 3), UNEXPLAINED or INCONCLUSIVE, never PASS.

### I.9 Design-verification probe (synthetic; not a B3Y result)
Script `probe_b3y_tier3_design.py` (imports no repository module). A single-channel chain (N_y = 1, E⊥ = threshold of channel 4, E = 88.5134459632352,
smooth two-bump potential), device lengths as in B3Y, exact answer from an mpmath transfer recursion at 60 digits using the same double inputs.

| dx | double LU solve \|t − t_exact\| | bulk-operator rounding alone | rounded matrix solved exactly (representation + lead data) | solver error \|t_dbl − t_fl(A)\| | M4-style extended refinement \|t − t_exact\| |
|---|---|---|---|---|---|
| 1/512 | 1.9×10⁻¹² | 1.2×10⁻¹² | 1.3×10⁻¹² | 6.4×10⁻¹³ | 6.9×10⁻¹⁹ |
| 1/1024 | 1.5×10⁻¹¹ | 1.8×10⁻¹¹ | 1.8×10⁻¹¹ | 3.4×10⁻¹² | 1.2×10⁻¹⁷ |
| 1/2048 | 7.9×10⁻¹¹ | 9.6×10⁻¹¹ | 9.5×10⁻¹¹ | 1.7×10⁻¹¹ | 1.8×10⁻¹⁷ |

* The refinement converges in two steps (update 8.8×10⁻¹³ → 3.6×10⁻¹⁹ at dx = 1/2048; contraction ≈ 4×10⁻⁷ per step), for both the exact and the rounded operator.
* Representation error dominates and grows 15× then 5× per halving; the solver error is smaller. A double-precision method can therefore be ≈ 8×10⁻¹¹ off the exact lattice
  answer at dx = 1/2048 even for this benign single-channel chain, close to the 10⁻¹⁰ target. This is the motivation for the precision-limited class.
* **What it does not show:** multichannel behaviour, the 63×63 dense blocks, the condition numbers of the real devices (κ̂₁ ≈ 10⁷ extrapolated), or the accuracy of M4 on real
  disorder. Those are exactly what gates X2–X8 test.

### I.10 Run order, outputs, reproducibility
Run order (a failed gate stops the sequence and is reported): G0 → X0–X7 on synthetic problems → X8/V1 → Tier 1 at 1/1024 → Tier 2 at 1/2048 with M4 → Tier 3 at 1/1024 → buffer control K1 → analysis.
Estimated cost (measured sparse solves, extrapolated to 2873 slices): about 10 s per sparse method per energy and 2.1 GB peak at 1/2048; M4 about two refinements of
≈ 10–20 s each. The whole matrix is expected to take on the order of 30 minutes.

Outputs in `results/rgf_b3y_validation/`: `per_point_comparison.csv`, `numerical_diagnostics.csv`, `tier3_exact_lattice.csv`, `tier3_error_decomposition.csv`, `xp_gates.csv`,
`convergence_orders.csv`, `buffer_control.csv`, `potential_identity.csv`, `summary.json`, `run_metadata.json`. Planned code: `sparse_exact_lattice.py`, `validate_rgf_b3y.py`,
`tests/test_rgf_b3y.py`, an additive `rgf_b3y` fingerprint kind in `provenance.py`. **Not modified:** `stationary.py`, `smatrix_spec.py`, `sparse_reference.py`, `sparse_wavefunction.py`,
`validate_rgf_b3x.py`, and every B1/B2/B3/B3X file and result (the B3X tests pin the hashes of the B3X solver files).
Metadata recorded: commit SHA and dirty flag, Python/NumPy/SciPy/mpmath versions, hardware (CPU model, cores, memory), precision details (`np.finfo` of float64 and longdouble, mpmath digits),
potential-array hashes, protocol hash and commit, exact energy points, solver settings (orderings, refinement counts, tolerances), B3X artifact hashes, test commands and results.

### I.11 Prediction register (recorded before any run; priors, not results and not bounds)
Basis: B3X observed RGF–wave-function ΔS divided by the B1-type bound at dx = 1/512, over 30 points (min 0.081, median 0.162, max 0.357), times the bound at the finer dx.
Predicted RGF–reference |ΔS| band at the priority points (full list in the JSON):

| Point | dx = 1/1024 | dx = 1/2048 |
|---|---|---|
| W200_z0.2, E = 88.5 | 6×10⁻¹² – 2.7×10⁻¹¹ | 2.4×10⁻¹¹ – 1.1×10⁻¹⁰ |
| W200_z0.05, E = 140.5 | 4.5×10⁻¹² – 2.0×10⁻¹¹ | 1.8×10⁻¹¹ – 8.0×10⁻¹¹ |
| W200_z0.05, E = 204.9 | 3.6×10⁻¹² – 1.6×10⁻¹¹ | 1.4×10⁻¹¹ – 6.3×10⁻¹¹ |
| W200_z0.02, E = 281.9 | 3.0×10⁻¹² – 1.3×10⁻¹¹ | 1.2×10⁻¹¹ – 5.2×10⁻¹¹ |
| W200_z0.05, E = 371.4 | 2.5×10⁻¹² – 1.1×10⁻¹¹ | 1.0×10⁻¹¹ – 4.4×10⁻¹¹ |

If these priors hold, no exceedance of 10⁻¹⁰ is expected at 1/1024 and a marginal one is possible at 1/2048 for W200_z0.2. The probe (§I.9) shows that a single method can sit at
8×10⁻¹¹ against the exact lattice answer, so exceedances at 1/2048 are plausible and are handled by §I.8.3.

### I.12 Tier 3 implementation and synthetic gate results (no B3Y data)
`sparse_exact_lattice.py` implements §I.6 and imports only the standard library, NumPy, SciPy and mpmath (an AST test forbids any import from, or reference to, `stationary`,
`smatrix_spec`, `sparse_reference`, `sparse_wavefunction` or the Fisher–Lee/solver function names). Main parts: platform check (`platform_precision`, `require_extended_precision`,
`ExtendedPrecisionUnavailable`), cancellation-free extended lead data with an mpmath scalar reference, extended-precision cosine-sum potential blocks, the unmerged operator
(`UnmergedOperator`) and an exact extended-precision product with a *rounded* double matrix (`apply_matrix_extended`), iterative refinement with the double SuperLU as a
preconditioner only, S from the scattering state and, from separate unit-source solves, by Fisher–Lee, the 1-norm condition estimate, the convergence/error indicator e_ind, the exact
solution of a rounded matrix (`solve_matrix_extended`, with an optional "hybrid" operator: rounded bulk, exact lead terms), the telescoping `error_decomposition`, the classification
functions of §I.8.3, and mpmath references (single-channel transfer matrix; dense-mixing small devices).

**Platform behaviour.** If `np.longdouble` has fewer than 64 mantissa bits every entry point raises `ExtendedPrecisionUnavailable("NOT EVALUABLE: …")`; there is no fallback to double precision. A test simulates a
53-bit longdouble (the Windows/MSVC or Apple arm64 case) and checks that every entry point refuses, and that a 113-bit type (aarch64 Linux quad) is accepted.
Measured here: x86-64 Linux, 64 mantissa bits, ε = 1.08×10⁻¹⁹. The extended precision is the hardware 80-bit x87 type; results are platform dependent in the sense of §I.6.

**Synthetic gate results** (clean guide and y-uniform B2 smooth barrier, N_y = 63, geometry as in B3Y, at synthetic energies E = 95.0, 150.25, 263.5 that are not B3Y energies;
single-channel and dense-mixing toys). Every number is measured; mpmath uses 50 digits.

| Gate | Problem | Result | Limit |
|---|---|---|---|
| X0 | this platform | 64 mantissa bits, evaluable | ≥ 64 |
| X1 | λ and Γ of all 63 channels, E = 95.0 and 263.5, dx 1/512, 1/1024, 1/2048 | relative 2.7×10⁻²⁰ – 1.1×10⁻¹⁹ | 10⁻¹⁷ |
| X2 | clean guide, t_nn vs mpmath λⁿ⁻¹, off-diagonal t, r (dx 1/512, 1/1024, 1/2048) | ≤ 3.7×10⁻¹⁸ (1/512), 6.2×10⁻¹⁸ (1/1024), 3.3×10⁻¹⁷ (1/2048) | 10⁻¹³ |
| X3 | smooth barrier, every open channel, t and r vs the mpmath single-channel transfer matrix | ≤ 6.2×10⁻¹⁸ (1/512), 1.3×10⁻¹⁷ (1/1024, 1/2048) | 10⁻¹³ |
| X4 | all 14 solves (state and Green) | 2 iterations; worst contraction 1.6×10⁻⁶; η_final ≤ 4.9×10⁻²⁰ | ≤ 0.1; ≤ 10⁻¹⁷ |
| X5 | scattering state vs Fisher–Lee (separate unit-source solve) | ≤ 5.9×10⁻¹⁷ | 10⁻¹³ |
| X6 | COLAMD vs MMD_AT_PLUS_A (barrier, 1/512) | 1.0×10⁻¹⁷ | 10⁻¹³ |
| X7 | convergence/error indicator e_ind | 2.7×10⁻¹⁵ (1/512) – 9.5×10⁻¹⁴ (1/2048); κ̂₁ up to 2.1×10⁶ | ≤ 10⁻¹¹ |
| dense mixing | 3-channel dense blocks, N = 24 (dx 1/512) and N = 30 (dx 1/2048), vs an exact mpmath dense solve | ≤ 10⁻¹³ (test limit); gates X4, X5, X7 met | 10⁻¹³ |
| dense toy | 6-channel dense device, dx = 0.05, vs RGF, Green-block and wave-function references | ≤ 10⁻¹² (test limit); unitarity, reciprocity ≤ 10⁻¹⁴ | 10⁻¹² |
| decomposition | single channel (E⊥ = threshold of channel 4), dx 1/512 and 1/2048: telescoping identity; M4 vs mpmath | identity ≤ 10⁻¹⁸; M4 error ≤ 10⁻¹³ | exact |

The e_ind values are dominated by κ̂₁·η; they are well below the 10⁻¹¹ threshold on these synthetic problems. Passing these gates on synthetic problems is **not** evidence about any B3Y point; gates X1–X3 at the B3Y energies and gate X8 are
run only in the numerical stage, after explicit authorization.

**Tests (actual commands and exit statuses; local runs, no CI).**

| Command | Exit status | Result |
|---|---|---|
| `python -m pytest -q tests/test_rgf_b3y.py` | 0 | 34 passed (≈ 60 s) |
| `python -m pytest -q tests` (complete suite) | 0 | 286 passed, 2 xfailed (the two known expected failures) |

The 34 tests cover: independence (AST) of the module, equality of every constant with the protocol, the revision-1 clarifications, protected-file and prior-artifact hashes, gate X0 including the simulated 53-bit
longdouble, X1 (six cases), the extended potential blocks, the unmerged-vs-merged operator, the exact extended product with a rounded matrix, X2 (four cases), X3 (three cases), X6, dense mixing
(mpmath dense solve and the three existing solvers on a toy), the error decomposition with a stand-in double method, gate X4 logic at its boundaries, e_ind, and the classification rules including their boundaries.
Audit results: no import of, or reference to, an existing solver module in the implementation; the 18 protected file hashes and the four prior result trees (B1, B2, B3, B3X) are unchanged; no B3Y result directory exists. *[historical: true at protocol revision 2; the directory `results/rgf_b3y_validation/` was created by the campaign, Part II]*

---

## Part I-B: Protocol self-audit

### I-B.1 Ambiguities found, and how each is resolved or flagged
| # | Ambiguity | Resolution | Needs your review |
|---|---|---|---|
| A1 | "exact failed energy points" | the five B3X rows with `passed = False`, identified by configuration id and full-precision E; the "371.4" label occurs twice | no |
| A2 | what happens when \|ΔS\| > 10⁻¹⁰ | §I.8.2–I.8.3: inconclusive until Tier 3 adjudicates. Revision 1: the factor-10 rule is an operational classification, not a proof | accepted in revision 1 |
| A3 | how reliable M4's own answer is | convergence/error indicator e_ind (§I.6), not a rigorous bound; adjudication requires e_ind ≤ 10⁻¹¹ | accepted in revision 1 (threshold unchanged) |
| A4 | the authorization does not name the methods of the buffer control | M1 and M3 (M4 excluded to keep it cheap) | yes |
| A5 | which k enters ρ | lowest-k open channel (largest open n); formula in the JSON | no (diagnostic only) |
| A6 | meaning of "full pass" | §I.8.3 item 5; "PASS" is unreachable with any exceedance | no |
| A7 | G0 cannot compare full S with B3X | B3X stored aggregated ΔS and T1 only; G0 checks T1 bitwise and ΔS within 10⁻¹³ | no |
| A8 | M4 vs the methods is not a pass/fail criterion at fine dx | it is the adjudication step (and gate X8 at 1/512); ε_m is always reported | yes |
| A9 | H1 spread limit 2×10⁻⁶ | 1 % of the 2×10⁻⁴ limit; a choice | yes |
| A10 | M4 gate tolerances (10⁻¹⁷, 10⁻¹³) are derived from the single-channel probe | if the 63-channel problem cannot meet them, the run stops (not silently relaxed) | yes |
| A11 | the buffer control compares different devices | plane-dependent complex S are not compared across buffers (§I.4) | no |
| A12 | M4 normalisation of the residual | η defined explicitly in §I.6 | no |

### I-B.2 Are the proposed comparisons mathematically independent?
| Pair | Independent | Shared |
|---|---|---|
| M1 vs M2 | solution algorithm (recursion vs sparse LU); V blocks; assembly | lead roots, Γ, Fisher–Lee, S assembly (`smatrix_spec`) |
| M1 vs M3 | algorithm, lead roots, V-block algebra, S extraction | problem definition |
| M2 vs M3 | lead data, V algebra, S extraction, sparse ordering | problem definition; both use SuperLU and one merged double matrix |
| M4 vs M1–M3 | arithmetic precision, operator form, lead formulas, Γ; M4 computes S in two ways | problem definition; the wave-function derivation is shared with M3 (M4's Fisher–Lee form, gate X5, covers the equivalence with M1/M2) |

**Common to all four:** the lattice Hamiltonian, the DST-I quadrature definition, the Dirichlet sine basis truncated at N_y = 63, the potential callable and its double samples,
the E and E⊥ doubles, the convention-P reference planes, the exact clean-lead elimination, and the neglect of the potential beyond the first lead slice (checked by P4 and
the buffer control). B3Y therefore cannot detect an error in those; they were examined in B1, B2 and the specification tests. The same libraries (NumPy/SciPy SuperLU) underlie M2–M4,
but M4's accuracy comes from the extended-precision residual, not from the LU.

### I-B.3 Can Tier 3 actually separate operator representation error from solver error?
**Argument.** The decomposition is an identity: S_double − S_exact = [S(x_fl(A)) − S_exact] + [S_double − S(x_fl(A))]. It is realized by two refinements, one against the unmerged exact
operator (giving S_exact) and one against the rounded matrix entries (giving S(x_fl(A))); both use extended-precision residuals and the same double LU as a preconditioner. It works if:
(i) both refinements converge (contraction ≈ κ̂ε_double ≈ 10⁻⁹ ≪ 1); (ii) the extended-precision error κ̂ε_ld ≲ 10⁻¹² is well below the quantities being separated (≈ 10⁻¹¹ – 10⁻¹⁰), i.e. e_ind ≤ 10⁻¹¹ (gate X7; an indicator, not a rigorous bound); (iii) the representation
term is large enough to resolve (predicted from the diagonal ulp to be ≈ 10⁻¹¹ – 10⁻¹⁰ at dx = 1/2048). A term is resolved only if it exceeds 10 × e_ind; otherwise it is unresolved and no dominance claim is allowed (UNEXPLAINED).
**Evidence so far (§I.9):** on a synthetic single-channel chain both refinements converge in two steps, M4-style error vs the mpmath exact lattice answer is ≤ 2×10⁻¹⁷ and the three terms are cleanly
separated (representation dominant, solver smaller).
**Limits:** (a) multichannel dense blocks and the real κ̂ are untested until gates X2–X8; (b) M1 has no explicit matrix, so its representation and solver errors cannot be separated, only its total;
(c) fl(A) in the decomposition is each method's own matrix, so the representation error is method-specific; (d) the exact inputs are the double E, E⊥, V samples, so rounding of E⊥ or V relative to the "ideal" real numbers is not part of the problem (it is identical in all methods).

### I-B.4 Residual risks
1. The 10⁻¹⁰ target may be exceeded by double-precision methods at dx = 1/2048 for round-off reasons (probe: 8×10⁻¹¹ for one benign channel); the protocol classifies this, it does not hide it.
2. M4's gates may fail on the 63-channel problem; then Tier 3 is reported as unresolved and nothing is inferred from it.
3. Extended precision depends on the platform (§I.6).
4. Only five points reach dx = 1/2048; orders from three grids are noisy; Richardson values are estimates.
5. Shared problem definition and one random draw (W = 200 is exactly twice W = 100).
6. The lead-tail margin at fine dx for W200_z0.02 is 4.3 %; a guard failure would make the point not evaluable, not pass.
7. A finite number of grids cannot prove continuum convergence; B3Y makes no such claim.

### I-B.5 Items I ask you to approve or change before any run
A2 and A3 were settled in revision 1. Still open from the first review: A4 (buffer control methods), A8, A9 (H1 spread limit), A10 (gate tolerances), plus the two new numbers of revision 1 (`dominance_factor` = 2, `term_resolution_factor` = 10).
`sparse_exact_lattice.py` is implemented with synthetic tests only; `validate_rgf_b3y.py` is not written, and nothing is run on B3Y points until the numerical stage is authorized. *[historical: the driver was written and the campaign authorized and run afterwards; see Part II]*

## Part II: Results

*Written after the run. Every number below is read from the files in `results/rgf_b3y_validation/` (generated by `validate_rgf_b3y.py` from the predeclared protocol, revision 2). Nothing in the protocol — points, realizations, thresholds, the 10⁻⁶ lead-tail guard, the absolute 10⁻¹⁰ S criterion — was changed to obtain these results.*

### II.0 Verdict and scope

**Under the predeclared criteria, B3Y has no failed criterion and no exceedance: overall status PASS** (`summary.json`: `overall_status`, `stopped_at = null`, 37 comparison rows, 0 NOT EVALUABLE). In words:

* The three double-precision implementations (M1 RGF, M2 Green-block, M3 wave-function) agree on the **complex S-matrix** to ≤ 9.1×10⁻¹¹ at every point and every resolution run (criterion 10⁻¹⁰, absolute). The closest approach to the limit is 9.09×10⁻¹¹ (M1–M2, `W200_z0.2@E88.51`, dx = 1/2048): **9 % below the limit**, see II.4.
* The independent extended-precision exact-lattice reference (M4) agrees with each of M1–M3 to ≤ 8.85×10⁻¹¹ at all 16 points (1/512) and at the five predeclared points at 1/1024 and 1/2048. Gates G0, X0–X8, V1 and P1–P4 all passed.
* The T1 changes between dx = 1/512, 1/1024, 1/2048 that made five W = 200 points fail the B3X dx criterion are **reproduced identically by all four methods** (spread ≤ 3.3×10⁻¹¹ ≪ 2×10⁻⁶, criterion H1). They are therefore a property of the **finite-difference lattice problem**, not of one solver. Their observed order is 1.996 – 2.006.
* Because no value exceeded 10⁻¹⁰, **no exceedance had to be adjudicated**: the PRECISION-CONSISTENT / PRECISION-LIMITED / UNEXPLAINED classification was not triggered anywhere (II.5 explains what the decomposition still shows).

What this does **not** establish: that the lattice T1 equals a continuum quantity (the lattice reference is the exact solution of the *discretized* problem, not of the continuum Schrödinger equation; the Richardson values are labelled estimates); anything about other realizations, seeds, energies or sizes (one realization, seed 0, 16 predeclared points); anything about Anderson localization, decoherence, ENAQT, or experimentally demonstrated devices. B3Y is a numerical-reliability campaign for the solver used in the physical question (how correlated static disorder and its strength affect wavepacket transport in a finite 2D hard-wall waveguide), not a physical result.

### II.1 What was run

| item | value |
|---|---|
| driver | `validate_rgf_b3y.py` (commit of the driver: `1a5688b`; protocol `docs/rgf_b3y_protocol.json`, sha256 in `run_metadata.json`) |
| command | `python validate_rgf_b3y.py --out results/rgf_b3y_validation --cache <scratch cache>` (first run: same command with a scratch output directory, 2351 s; the files in the repository were regenerated from the per-case caches in 117 s, content identical apart from timing columns — checked) |
| points | 16 predeclared (15 × W = 200 and the W = 100 control), exact energies from the protocol; dx = 1/512 (719 slices), 1/1024 (1437), 1/2048 (2873; five points); buffer control 205/512 (1845 slices) |
| methods | M1 `stationary.scattering_matrix`; M2 `sparse_reference` (COLAMD); M3 `sparse_wavefunction` (MMD_AT_PLUS_A + 1 refinement); M4 `sparse_exact_lattice` (x86-64 80-bit `np.longdouble`, 64 mantissa bits, arithmetic probe passed; mpmath 1.3.0, 50 digits) |
| platform | Linux x86-64, Intel Xeon 2.10 GHz, 4 cores, 16 GB; python 3.13.18, numpy 2.5.3 (`platform_precision.json`) |
| per-case cost | M1 0.4 / 0.8 / 1.6 s, M2 3.8 / 7.6 / 15 s, M3 2.3 / 4.7 / 9.3 s (mean, 1/512 / 1/1024 / 1/2048); M4 with Green solve 3 – 22 s |

Each method keeps its own implementation: M1–M3 are the unmodified protected solvers; the only object shared by all methods is the physical potential function (checked bitwise, P1). The decomposition for M2 and M3 re-runs the method (and checks that the recomputed double S equals the reported one bitwise: `M2_recomputed_bitwise`, `M3_recomputed_bitwise` true in all 30 rows).

### II.2 Gates (all passed)

| gate | rows | worst value | limit | passed |
|---|---|---|---|---|
| G0 | 16 | 0.00e+00 | 1e-13 | True |
| X0 | 1 | 6.40e+01 | 64 | True |
| X1 | 48 | 1.30e-19 | 1e-17 | True |
| X2 | 15 | 4.97e-17 | 1e-13 | True |
| X3 | 15 | 4.35e-17 | 1e-13 | True |
| X6 | 10 | 3.79e-17 | 1e-13 | True |
| X8 | 1 | 5.82e-12 | 1e-10 | True |

* **G0** (the recomputed 1/512 values reproduce the B3X artifacts): RGF T1 bitwise equal and the three reference differences equal to the B3X values (worst difference 0.0) at all 16 points.
* **X1** over all 16 energies × 3 dx; **X2** clean guide and **X3** y-uniform smooth barrier (B2 barrier shape, centred in the disorder interval) at the five predeclared energies × 3 dx, with the real geometries (up to 2873 slices); **X4/X5/X7** were evaluated for every M4 solve (26 solves): worst refinement contraction 1.7×10⁻⁵ (state), 3.3×10⁻⁵ (Green) vs 0.1; final η ≤ 4.1×10⁻²⁰ vs 10⁻¹⁷; Fisher–Lee vs state ≤ 5.8×10⁻¹⁷ vs 10⁻¹³; e_ind ≤ 5.7×10⁻¹³ vs 10⁻¹¹; κ̂₁ up to 1.8×10⁷.
* **X8/V1:** at 1/512 the largest |S_M4 − S_Mk| over all 16 points and k = 1,2,3 is 5.8×10⁻¹² (limit 10⁻¹⁰).
* **P1** potential identity: the V arrays used by the RGF device, M2 and M3 are bitwise identical to the direct `Potential(cfg, disorder_from_config(cfg, 0)).numpy` call in every one of the 42 device records, and the E_perp doubles of the three methods are bitwise identical. **P2:** potential at coincident physical slices is bitwise equal across 1/512, 1/1024, 1/2048 (6 comparisons). **P3:** hashes in `potential_identity.csv` and `run_metadata.json`. **P4** lead-tail guard: max |V| on the lead slices is 9.57×10⁻⁷ over all devices, **below the 10⁻⁶ guard but only by 4 %**; no device was NOT EVALUABLE.

### II.3 Tier 1 and Tier 2: M1 – M3 agreement, unitarity, reciprocity

| dx | points | max S1 (M1–M2) | max S2 (M1–M3) | max S3 (M2–M3) | max S4 | max unitarity | max reciprocity | max backward error | max ρ (diag.) |
|---|---|---|---|---|---|---|---|---|---|
| 1/512 | 16 | 5.01e-12 | 5.33e-12 | 1.08e-12 | 2.57e-12 | 9.85e-13 | 8.80e-13 | 8.58e-16 | 0.29 |
| 1/1024 | 16 | 1.55e-11 | 1.45e-11 | 4.10e-12 | 5.85e-12 | 2.31e-12 | 1.62e-12 | 1.64e-15 | 0.35 |
| 1/2048 | 5 | 9.09e-11 | 7.92e-11 | 1.95e-11 | 4.00e-11 | 1.39e-11 | 1.45e-11 | 3.16e-15 | 0.30 |

All of S1 (M1–M2), S2 (M1–M3), S3 (M2–M3), S4 (T, R and off-diagonal transmission probabilities, limit 2×10⁻¹⁰), S5 (unitarity, both products), S6 (reciprocity: full, blocks, probabilities, limit 10⁻⁸), R1 (sparse backward error ≤ 10⁻¹⁴) and N1 (channel counts equal to the analytic count, Γ_closed = 0, closed λ real with |λ| < 1, ≥ 1 closed channel) are PASS at every row (`status_*` columns). Open channel counts range 4 – 8 as predeclared; the ρ column is the diagnostic only (never used to classify).
The sparse solves are well conditioned at the backward-error level (≤ 3.2×10⁻¹⁵); the S discrepancy growth with dx is the round-off accumulation across a longer chain (II.4).

Reciprocity of the *physical* S (max |S − Sᵀ|) and unitarity are at the 10⁻¹² – 10⁻¹¹ level at the finest grid, consistent with the S differences. Closed-channel amplitude differences (report-only, columns `N1_closed_amp_*`; ≤ 2.8×10⁻¹⁴) are not part of any criterion.

### II.4 Tier 3: the exact lattice solution (M4)

M4 is the exact solution of the *lattice* problem (unmerged operator, cancellation-free lead data, extended-precision residuals). ε is max |S_Mk − S_M4| over all elements of the complex S.

| point | dx | ε(M1) | ε(M2) | ε(M3) | max\|T1(Mk)−T1(M4)\| | e_ind | κ̂₁ |
|---|---|---|---|---|---|---|---|
| W200_z0.02@E281.90 | 1/512 | 1.37e-12 | 1.49e-12 | 1.44e-12 | 1.20e-13 | 9.60e-15 | 3.1e+05 |
| W200_z0.02@E281.90 | 1/1024 | 3.65e-12 | 4.12e-12 | 4.15e-12 | 4.69e-13 | 3.21e-14 | 1.2e+06 |
| W200_z0.02@E281.90 | 1/2048 | 1.73e-11 | 1.46e-11 | 1.65e-11 | 1.67e-12 | 1.70e-13 | 4.7e+06 |
| W200_z0.05@E140.47 | 1/512 | 1.53e-12 | 1.03e-12 | 1.03e-12 | 7.05e-13 | 5.99e-15 | 2.3e+05 |
| W200_z0.05@E140.47 | 1/1024 | 4.22e-12 | 2.65e-12 | 3.97e-12 | 1.78e-12 | 3.70e-14 | 1.0e+06 |
| W200_z0.05@E140.47 | 1/2048 | 2.24e-11 | 1.90e-11 | 1.46e-11 | 7.36e-12 | 1.26e-13 | 4.6e+06 |
| W200_z0.05@E204.93 | 1/512 | 1.87e-12 | 5.16e-13 | 6.13e-13 | 2.52e-13 | 9.93e-15 | 3.4e+05 |
| W200_z0.05@E204.93 | 1/1024 | 4.10e-12 | 1.79e-12 | 1.56e-12 | 1.03e-12 | 3.50e-14 | 1.3e+06 |
| W200_z0.05@E204.93 | 1/2048 | 1.77e-11 | 7.96e-12 | 1.00e-11 | 1.51e-12 | 1.49e-13 | 4.5e+06 |
| W200_z0.05@E371.36 | 1/512 | 6.95e-13 | 8.34e-13 | 8.01e-13 | 1.51e-13 | 6.95e-15 | 2.0e+05 |
| W200_z0.05@E371.36 | 1/1024 | 5.68e-12 | 2.73e-12 | 3.19e-12 | 6.53e-13 | 2.54e-14 | 8.9e+05 |
| W200_z0.05@E371.36 | 1/2048 | 1.06e-11 | 1.64e-11 | 1.54e-11 | 3.31e-12 | 1.06e-13 | 2.8e+06 |
| W200_z0.2@E88.51 | 1/512 | 1.51e-12 | 5.50e-12 | 5.82e-12 | 2.16e-12 | 3.49e-14 | 1.2e+06 |
| W200_z0.2@E88.51 | 1/1024 | 1.82e-11 | 2.53e-11 | 2.56e-11 | 8.66e-12 | 1.10e-13 | 4.6e+06 |
| W200_z0.2@E88.51 | 1/2048 | 2.13e-11 | 8.85e-11 | 7.70e-11 | 3.27e-11 | 5.67e-13 | 1.8e+07 |

All 16 points at 1/512 (gate X8; the five points of the table above reappear here with the same numbers):

| point | ε(M1) | ε(M2) | ε(M3) |
|---|---|---|---|
| W100_z0.05@E371.36 | 7.16e-13 | 7.34e-13 | 7.73e-13 |
| W200_z0.02@E88.51 | 3.23e-12 | 1.25e-12 | 8.43e-13 |
| W200_z0.02@E140.47 | 2.72e-12 | 9.08e-13 | 5.99e-13 |
| W200_z0.02@E204.93 | 1.13e-12 | 4.86e-13 | 5.39e-13 |
| W200_z0.02@E281.90 | 1.37e-12 | 1.49e-12 | 1.44e-12 |
| W200_z0.02@E371.36 | 9.78e-13 | 6.72e-13 | 7.64e-13 |
| W200_z0.05@E88.51 | 1.19e-12 | 6.56e-13 | 6.39e-13 |
| W200_z0.05@E140.47 | 1.53e-12 | 1.03e-12 | 1.03e-12 |
| W200_z0.05@E204.93 | 1.87e-12 | 5.16e-13 | 6.13e-13 |
| W200_z0.05@E281.90 | 9.54e-13 | 1.60e-12 | 1.57e-12 |
| W200_z0.05@E371.36 | 6.95e-13 | 8.34e-13 | 8.01e-13 |
| W200_z0.2@E88.51 | 1.51e-12 | 5.50e-12 | 5.82e-12 |
| W200_z0.2@E140.47 | 2.29e-12 | 7.60e-13 | 8.58e-13 |
| W200_z0.2@E204.93 | 1.40e-12 | 3.69e-13 | 4.42e-13 |
| W200_z0.2@E281.90 | 6.21e-13 | 1.72e-12 | 1.65e-12 |
| W200_z0.2@E371.36 | 1.10e-12 | 8.65e-13 | 6.00e-13 |

* Every ε is below 10⁻¹⁰; the largest is **8.85×10⁻¹¹** (M2, `W200_z0.2@E88.51`, 1/2048), the next 7.7×10⁻¹¹ (M3) and 2.1×10⁻¹¹ (M1) at the same point.
* ε grows with every halving of dx (factors 1.2 – 12.1; the median is ≈ 4). For the three double methods at the point closest to the limit the factors from 1/512 → 1/1024 → 1/2048 are 12.1 → 1.2 (M1), 4.6 → 3.5 (M2), 4.4 → 3.0 (M3). One growth factor (M1, `W200_z0.2@E88.51`, 1/512 → 1/1024: 12.1) exceeds the operational factor 10 of rule §I.8.3; this rule is defined only for classifying a value that exceeds 10⁻¹⁰, so it was **not** applied (nothing exceeded), and it starts from a very small ε(1/512) = 1.5×10⁻¹². It is reported here as observed.
* **Margin.** At the worst point ε(1/2048) is 11 – 12 % below the limit and has been growing by ≈ 3 – 4× per halving. An exceedance at a hypothetical 1/4096 would therefore be expected, but dx = 1/4096 is not part of B3Y and was not run: that is an extrapolation, not a result.
* |T1(M_k) − T1(M4)| ≤ 3.3×10⁻¹¹ everywhere (T1 is much less sensitive than the complex S).
* The prediction register (§I.11, written before the runs) predicted RGF–reference |ΔS| bands at the five priority points at 1/1024 and 1/2048: **all 10 (point, dx) cases lie inside their predicted band, for both M1–M3 and M1–M2** (e.g. `W200_z0.2@E88.51`, 1/2048: predicted 2.4×10⁻¹¹ – 1.1×10⁻¹⁰, observed 7.9×10⁻¹¹ (M1–M3) and 9.1×10⁻¹¹ (M1–M2)). The priors were not bounds and could have failed; they did not.

### II.5 Error decomposition (M2 and M3 only)

Three-term decomposition S_double − S_exact = (bulk-operator rounding) + (lead data) + (solver/extraction), each method with *its own* rounded matrix (`sparse_reference.assemble_M` for M2, `sparse_wavefunction.assemble` for M3), double lead data, source construction and extraction (the mirrors reproduce the methods' own sources and S bitwise: `mirror_source_bitwise`, `mirror_extraction_bitwise` true in all 30 rows); the hybrid solve replaces only the two lead (corner) terms with exact extended-precision ones. M1 has no explicit matrix and receives no decomposition. The decomposition was computed at all five points at 1/512, 1/1024, 1/2048 although no exceedance made it necessary.

| point | dx | method | total | bulk rounding | lead data | solver/extraction | e_ind (combined) | resolved | bulk ≥ 2(lead+solver) | claim supported |
|---|---|---|---|---|---|---|---|---|---|---|
| W200_z0.02@E281.90 | 1/512 | M2 | 1.49e-12 | 1.42e-12 | 6.37e-14 | 3.26e-13 | 1.87e-14 | False | True | False |
| W200_z0.02@E281.90 | 1/512 | M3 | 1.44e-12 | 1.43e-12 | 5.08e-14 | 2.29e-13 | 1.76e-14 | False | True | False |
| W200_z0.05@E140.47 | 1/512 | M2 | 1.03e-12 | 9.44e-13 | 3.94e-14 | 5.17e-13 | 8.69e-15 | False | False | False |
| W200_z0.05@E140.47 | 1/512 | M3 | 1.03e-12 | 9.45e-13 | 4.13e-14 | 3.67e-13 | 1.06e-14 | False | True | False |
| W200_z0.05@E204.93 | 1/512 | M2 | 5.16e-13 | 5.33e-13 | 3.02e-14 | 3.00e-13 | 1.84e-14 | False | False | False |
| W200_z0.05@E204.93 | 1/512 | M3 | 6.13e-13 | 5.34e-13 | 2.53e-14 | 2.72e-13 | 1.98e-14 | False | False | False |
| W200_z0.05@E371.36 | 1/512 | M2 | 8.34e-13 | 5.83e-13 | 4.49e-14 | 2.68e-13 | 1.01e-14 | False | False | False |
| W200_z0.05@E371.36 | 1/512 | M3 | 8.01e-13 | 5.82e-13 | 5.36e-14 | 3.05e-13 | 1.29e-14 | False | False | False |
| W200_z0.2@E88.51 | 1/512 | M2 | 5.50e-12 | 6.33e-12 | 1.05e-13 | 9.48e-13 | 9.60e-14 | False | True | False |
| W200_z0.2@E88.51 | 1/512 | M3 | 5.82e-12 | 6.32e-12 | 2.36e-13 | 1.11e-12 | 8.37e-14 | False | True | False |
| W200_z0.02@E281.90 | 1/1024 | M2 | 4.12e-12 | 3.28e-12 | 2.57e-13 | 9.31e-13 | 7.41e-14 | False | True | False |
| W200_z0.02@E281.90 | 1/1024 | M3 | 4.15e-12 | 3.29e-12 | 2.83e-13 | 8.52e-13 | 8.24e-14 | False | True | False |
| W200_z0.05@E140.47 | 1/1024 | M2 | 2.65e-12 | 4.29e-12 | 9.23e-13 | 1.94e-12 | 5.22e-14 | True | False | False |
| W200_z0.05@E140.47 | 1/1024 | M3 | 3.97e-12 | 4.29e-12 | 7.90e-13 | 1.21e-12 | 6.48e-14 | True | True | True |
| W200_z0.05@E204.93 | 1/1024 | M2 | 1.79e-12 | 9.94e-13 | 2.40e-13 | 1.05e-12 | 8.36e-14 | False | False | False |
| W200_z0.05@E204.93 | 1/1024 | M3 | 1.56e-12 | 9.94e-13 | 3.44e-13 | 9.74e-13 | 7.38e-14 | False | False | False |
| W200_z0.05@E371.36 | 1/1024 | M2 | 2.73e-12 | 3.13e-12 | 1.16e-13 | 8.52e-13 | 4.57e-14 | False | True | False |
| W200_z0.05@E371.36 | 1/1024 | M3 | 3.19e-12 | 3.13e-12 | 1.26e-13 | 5.13e-13 | 5.49e-14 | False | True | False |
| W200_z0.2@E88.51 | 1/1024 | M2 | 2.53e-11 | 2.53e-11 | 1.47e-12 | 2.77e-12 | 3.85e-13 | False | True | False |
| W200_z0.2@E88.51 | 1/1024 | M3 | 2.56e-11 | 2.53e-11 | 8.49e-13 | 2.38e-12 | 2.95e-13 | False | True | False |
| W200_z0.02@E281.90 | 1/2048 | M2 | 1.46e-11 | 1.68e-11 | 5.04e-13 | 2.98e-12 | 2.97e-13 | False | True | False |
| W200_z0.02@E281.90 | 1/2048 | M3 | 1.65e-11 | 1.68e-11 | 9.75e-13 | 1.82e-12 | 2.84e-13 | False | True | False |
| W200_z0.05@E140.47 | 1/2048 | M2 | 1.90e-11 | 1.25e-11 | 2.96e-12 | 6.82e-12 | 1.77e-13 | True | False | False |
| W200_z0.05@E140.47 | 1/2048 | M3 | 1.46e-11 | 1.25e-11 | 2.33e-12 | 3.68e-12 | 2.08e-13 | True | True | True |
| W200_z0.05@E204.93 | 1/2048 | M2 | 7.96e-12 | 7.11e-12 | 1.93e-12 | 5.44e-12 | 2.60e-13 | False | False | False |
| W200_z0.05@E204.93 | 1/2048 | M3 | 1.00e-11 | 7.11e-12 | 2.60e-12 | 3.09e-12 | 2.66e-13 | False | False | False |
| W200_z0.05@E371.36 | 1/2048 | M2 | 1.64e-11 | 1.52e-11 | 4.52e-13 | 7.28e-12 | 2.14e-13 | False | False | False |
| W200_z0.05@E371.36 | 1/2048 | M3 | 1.54e-11 | 1.52e-11 | 3.12e-13 | 6.12e-12 | 1.87e-13 | False | True | False |
| W200_z0.2@E88.51 | 1/2048 | M2 | 8.85e-11 | 6.53e-11 | 4.53e-12 | 2.25e-11 | 1.55e-12 | False | True | False |
| W200_z0.2@E88.51 | 1/2048 | M3 | 7.70e-11 | 6.52e-11 | 6.67e-12 | 1.05e-11 | 1.68e-12 | False | True | False |

What the 30 rows show, and what they do not:

* The telescoping identity residuals are 0.0 in all rows (an arithmetic check of the differences only; it does **not** validate the individual S-matrices).
* The bulk-operator rounding term is the largest of the three in 29 of 30 rows and satisfies the predeclared sum-form criterion `bulk ≥ 2 × (lead + solver)` in 18 of 30 rows.
* The claim is **supported** only where the terms are also *resolved* (each term > 10 × e_ind, e_ind = the largest convergence indicator among the extended solves of the decomposition). That holds in 4 of 30 rows, and in 2 of those (**M3 at `W200_z0.05@E140.47`, 1/1024 and 1/2048**) the dominance criterion also holds. In the other 26 rows the lead-data term is ≤ 10 × e_ind (and in 6 of them also the solver term; e_ind is dominated by κ̂₁·η and is 0.9×10⁻¹⁴ – 1.7×10⁻¹²), so the decomposition is *unresolved* and, by the predeclared rule, supports **no** representation-dominance claim. No claim is made.
* M2 at the same point (`W200_z0.05@E140.47`) is resolved but not dominated by the bulk term (bulk 1.25×10⁻¹¹ vs lead + solver 9.8×10⁻¹²). The decomposition therefore does not give a uniform mechanism statement across methods or points; and the sizes (a few 10⁻¹² – 10⁻¹¹) are comparable with the resolution limit of the indicator.
* Because ε never exceeded 10⁻¹⁰, nothing needed to be classified; the decomposition is reported as a measurement of the sizes of the three contributions, not as an adjudication.

### II.6 Resolution convergence of T1 (the question B3X left open)

| point | method | T1(1/512) | T1(1/1024) | T1(1/2048) | observed order p | in [1.8, 2.2] | Richardson ESTIMATE |
|---|---|---|---|---|---|---|---|
| W200_z0.02@E281.90 | M1 | 0.5909559968 | 0.5906976826 | 0.5906329115 | 1.9957 | True | 0.59061124 |
| W200_z0.02@E281.90 | M2 | 0.5909559968 | 0.5906976826 | 0.5906329115 | 1.9957 | True | 0.59061124 |
| W200_z0.02@E281.90 | M3 | 0.5909559968 | 0.5906976826 | 0.5906329115 | 1.9957 | True | 0.59061124 |
| W200_z0.02@E281.90 | M4 | 0.5909559968 | 0.5906976826 | 0.5906329115 | 1.9957 | True | 0.59061124 |
| W200_z0.05@E140.47 | M1 | 0.3020774674 | 0.3014806098 | 0.3013316543 | 2.0025 | True | 0.30128212 |
| W200_z0.05@E140.47 | M2 | 0.3020774674 | 0.3014806098 | 0.3013316543 | 2.0025 | True | 0.30128212 |
| W200_z0.05@E140.47 | M3 | 0.3020774674 | 0.3014806098 | 0.3013316543 | 2.0025 | True | 0.30128212 |
| W200_z0.05@E140.47 | M4 | 0.3020774674 | 0.3014806098 | 0.3013316543 | 2.0025 | True | 0.30128212 |
| W200_z0.05@E204.93 | M1 | 0.1037248306 | 0.1040044352 | 0.1040744259 | 1.9982 | True | 0.10409780 |
| W200_z0.05@E204.93 | M2 | 0.1037248306 | 0.1040044352 | 0.1040744259 | 1.9982 | True | 0.10409780 |
| W200_z0.05@E204.93 | M3 | 0.1037248306 | 0.1040044352 | 0.1040744259 | 1.9982 | True | 0.10409780 |
| W200_z0.05@E204.93 | M4 | 0.1037248306 | 0.1040044352 | 0.1040744259 | 1.9982 | True | 0.10409780 |
| W200_z0.05@E371.36 | M1 | 0.3946762322 | 0.3943964766 | 0.3943268045 | 2.0055 | True | 0.39430370 |
| W200_z0.05@E371.36 | M2 | 0.3946762322 | 0.3943964766 | 0.3943268045 | 2.0055 | True | 0.39430370 |
| W200_z0.05@E371.36 | M3 | 0.3946762322 | 0.3943964766 | 0.3943268045 | 2.0055 | True | 0.39430370 |
| W200_z0.05@E371.36 | M4 | 0.3946762322 | 0.3943964766 | 0.3943268045 | 2.0055 | True | 0.39430370 |
| W200_z0.2@E88.51 | M1 | 0.3817580555 | 0.3826745561 | 0.3829029205 | 2.0048 | True | 0.38297871 |
| W200_z0.2@E88.51 | M2 | 0.3817580555 | 0.3826745561 | 0.3829029205 | 2.0048 | True | 0.38297871 |
| W200_z0.2@E88.51 | M3 | 0.3817580555 | 0.3826745561 | 0.3829029205 | 2.0048 | True | 0.38297871 |
| W200_z0.2@E88.51 | M4 | 0.3817580555 | 0.3826745561 | 0.3829029205 | 2.0048 | True | 0.38297871 |

* All four methods give T1 identical to the printed 10 digits at every dx (H1: spread ≤ 3.3×10⁻¹¹ ≪ 2×10⁻⁶, 15 of 15 cases). The B3X observation that T1 changes by up to 9.2×10⁻⁴ between dx = 1/512 and 1/1024 is therefore **solver-independent** and, being reproduced by the exact lattice solution M4, **is the dx dependence of the finite-difference problem itself**, not a defect of the RGF (or of M2, M3).
* The observed orders are 1.996 – 2.006, inside [1.8, 2.2] for all 20 (point, method) cases. **Protocol ambiguity (flagged):** O1 says "both orders in [1.8, 2.2]", but three grids give a single observed order. The window was applied to that one order. The Richardson column is an **estimate** (T₂₀₄₈ + (T₂₀₄₈ − T₁₀₂₄)/(2ᵖ − 1)), it differs from the dx = 1/2048 value by 2×10⁻⁵ (`W200_z0.02@E281.90`) to 7.6×10⁻⁵ (`W200_z0.2@E88.51`), and it is **not** a statement about the continuum (it assumes an asymptotic p ≈ 2 regime that three grids can only suggest).
* The five points keep failing the B3X dx criterion at 1/512 vs 1/1024 by construction (the predeclared values are reproduced: for all 16 points, |ΔT1| (1/512 → 1/1024) of M1 equals the B3X value to the printed precision, e.g. 9.165×10⁻⁴ at `W200_z0.2@E88.51`). B3Y does not turn this into a pass of the B3X criterion; it changes its interpretation from "possibly a solver defect" to "reproduced by the exact lattice solution".

### II.7 Buffer control K1 (205/512 vs 103/512 at dx = 1/1024)

| point | \|ΔT1\| M1 | \|ΔT1\| M3 | M1–M3 max\|ΔS\| (103/512) | M1–M3 max\|ΔS\| (205/512) | max\|Δ\|S_ij\|²\| M1 across buffers (diag.) |
|---|---|---|---|---|---|
| W200_z0.02@E281.90 | 4.60e-11 | 4.65e-11 | 3.38e-12 | 3.96e-12 | 8.94e-11 |
| W200_z0.05@E140.47 | 1.57e-11 | 1.57e-11 | 6.39e-12 | 9.63e-12 | 3.25e-10 |
| W200_z0.05@E204.93 | 2.51e-11 | 2.51e-11 | 5.00e-12 | 7.49e-12 | 1.49e-10 |
| W200_z0.05@E371.36 | 3.80e-12 | 3.95e-12 | 4.14e-12 | 7.40e-12 | 8.80e-11 |
| W200_z0.2@E88.51 | 2.03e-11 | 2.19e-11 | 7.49e-12 | 7.86e-12 | 2.61e-10 |

Both T1 differences are ≤ 4.7×10⁻¹¹ (limit 2×10⁻⁴; M1 and M3), and the complex-S agreement of M1 and M3 holds within each geometry (≤ 9.7×10⁻¹²). The reference planes move with the buffer, so only plane-independent quantities are compared across the two buffers; the largest element of |ΔS_ij|² (diagnostic only, not a criterion) is 3.2×10⁻¹⁰ — a longer device also contains the 10⁻⁶-level tails of the potential, which is the expected size of a buffer effect. Reported separately from the primary results.

### II.8 Classification summary

| category | cases |
|---|---|
| Successful validation (every predeclared criterion met) | S1–S6, R1, N1, P1, P2, P4, G0, X0–X8, V1, K1, H1 for all rows |
| PRECISION-CONSISTENT / PRECISION-LIMITED | none (no value exceeded 10⁻¹⁰, so the rules were not applied) |
| UNEXPLAINED | none |
| Inconclusive / NOT EVALUABLE | none. The dominance decomposition is *unresolved* in 26 of 30 rows (no claim made); this affects no criterion |
| Classification only | O1 orders (1.996 – 2.006; one order per case, ambiguity above) |

### II.9 Deviations, additions and limitations (all disclosed)

1. **Two edits to the unprotected M4 module `sparse_exact_lattice.py`** were made while writing the driver, before any B3Y number existed: (a) `apply_matrix_extended` now bounds its temporary by a 256 MiB budget (chunking changes no value: each row's sum is evaluated in the same order); (b) the 1-norm condition estimate `onenormest` is seeded (it uses the global NumPy RNG, which made κ̂₁ and therefore e_ind vary from run to run; found by the reproducibility test). The seed changes the *estimate* κ̂₁ (a lower bound on the condition number) and e_ind, never S. No protected file changed.
2. **File names.** The protocol lists `per_point_comparison.csv`, `convergence_orders.csv`, `tier3_error_decomposition.csv`; this stage's request named `per_energy_comparison.csv`, `refinement_convergence.csv`, `error_decomposition.csv`, `platform_precision.json`, which are the names produced. Also written: `solver_spread.csv` (H1), `xp_gates.csv`, `tier3_exact_lattice.csv`, `potential_identity.csv` (including the P2 rows), `buffer_control.csv`, `summary.json`, `run_metadata.json`. `numerical_diagnostics.csv` of the protocol is subsumed in the columns of `per_energy_comparison.csv`.
3. **Status of S4 exceedance and of S5/S6/R1:** the protocol's vocabulary names PRECISION-LIMITED only for S1–S3. The driver treats an S4 exceedance as inconclusive only if S1 or S2 also exceed (otherwise FAIL) and S5, S6, R1, N1 as FAIL when violated; this never mattered (no exceedance).
4. **Decomposition scope.** It was computed also where nothing exceeded; it uses M2/M3 solved with the same SuperLU orderings (COLAMD / MMD_AT_PLUS_A) as the methods, M2 on its open-channel unit sources only (the only columns its Fisher–Lee extraction uses).
5. **Scope limits.** One disorder realization (seed 0); 16 predeclared energies; dx ≥ 1/2048; fixed N_y = 63, N_q = 127 (N_y and N_q refinements were B3X, not repeated); the 5 W = 200 points at 1/2048 are the only ones with a third grid. ε is a maximum over matrix elements, not a bound for derived quantities. e_ind is an indicator, not an uncertainty bound. A pass of all criteria is a statement about these points and these lattices.
6. **No claim** is made about Anderson localization, decoherence, ENAQT or experimentally demonstrated devices, nor that the lattice reference approximates the continuum solution.

### II.10 Files and reproduction

`results/rgf_b3y_validation/`: `per_energy_comparison.csv` (37 rows), `tier3_exact_lattice.csv` (26), `error_decomposition.csv` (30), `refinement_convergence.csv` (20), `buffer_control.csv` (5), `potential_identity.csv` (48), `xp_gates.csv` (106), `solver_spread.csv` (15), `summary.json`, `run_metadata.json`, `platform_precision.json`.
Tests: `tests/test_rgf_b3y.py` (50, module) and `tests/test_rgf_b3y_driver.py` (20, driver; synthetic problem only): method-specific matrix construction, source and extraction conventions (bitwise mirrors), hybrid lead replacement, output schema, missing/unresolved terms, metadata reproducibility across two independent runs, protocol-threshold preservation, refusal to write into historical result directories. Full suite: 322 passed, 2 xfailed (the two pre-existing expected failures), exit status 0; `run_metadata.json` carries the recorded commands (`post_run_tests`).
```
python validate_rgf_b3y.py --synthetic --out <scratch> --cache <scratch>   # synthetic driver check (no B3Y point)
python validate_rgf_b3y.py --out results/rgf_b3y_validation --cache <cache dir>
python validate_rgf_b3y.py --record-tests --out results/rgf_b3y_validation
```

### II.11 Post-campaign documentation and provenance audit

Performed after the campaign (no new numerical result was computed; every number below is recomputed from the existing CSV files). Scope: documentation, metadata, consistency of the PASS with the criteria, and one precision-consistency observation.

**(a) Statements that implied "no B3Y comparison has been run", and their disposition**

| location | statement | disposition |
|---|---|---|
| `docs/RGF_B3Y_VALIDATION.md` header (lines 7–10) | status line and "only numerical work … synthetic probe" | replaced by the four-layer status table; remaining sentence marked *[historical]* |
| `docs/RGF_B3Y_VALIDATION.md` I-B.4/I.12 area | "no B3Y result directory exists" | kept, marked *[historical]* |
| `docs/RGF_B3Y_VALIDATION.md` I-B.5 | "`validate_rgf_b3y.py` is not written …" | kept, marked *[historical]* |
| `docs/RGF_B3Y_TIER3_AUDIT.md` | "no B3Y comparison was run", "driver not yet written", "no B3Y result file exists" | banner added at the top; text unchanged |
| `docs/rgf_b3y_protocol.json` `status`, `implementation_status` | "no B3Y numerical comparison has been run", "NO B3Y point has been computed" | **deliberately not edited** (predeclared record; its sha256 is recorded in `run_metadata.json` and is unchanged). Its historical meaning is stated in the status table |
| `results/rgf_b3y_validation/summary.json` `protocol_status` | copy of the protocol's historical text | the only metadata inconsistency in the results; corrected narrowly: key renamed `protocol_status_at_predeclaration`, new `campaign_status = "COMPLETED real B3Y numerical campaign"`, and a `metadata_corrections` entry. No numerical value changed; no other result file was touched |
| `validate_rgf_b3y.py` | wrote `protocol_status` | now writes `protocol_status_at_predeclaration` and `campaign_status` (so a re-run reproduces the corrected summary). `run_metadata.json` keeps the code fingerprint of the driver that generated the numbers (`1a5688b`); the later edit changes only these two summary keys |
| code docstrings / tests | "no B3Y point is used" | accurate (they describe the synthetic mode and tests); unchanged |

No other file (README, other reports) mentions B3Y as un-run.

**(b) Is the PASS consistent with the criteria and the adjudication logic?** Recomputed independently from `per_energy_comparison.csv`, `tier3_exact_lattice.csv`, `buffer_control.csv`, `solver_spread.csv`, `xp_gates.csv`, `potential_identity.csv` (not from `summary.json`):
S1–S3 ≤ 10⁻¹⁰ in 37/37 rows; S4 ≤ 2×10⁻¹⁰; unitarity ≤ 10⁻¹⁰ and reciprocity ≤ 10⁻⁸ for M1–M3; backward error ≤ 10⁻¹⁴; lead-tail ≤ 10⁻⁶; N1 in 37/37; all 26 M4 solves passed X4, X5, X7; every ε(Mk, M4) ≤ 10⁻¹⁰ (26 × 3); K1 5/5; H1 15/15; 16 points at 1/512 and 1/1024, 5 at 1/2048, as predeclared. Per the protocol's overall rule ("PASS only if every criterion is met with no exceedance") the PASS is consistent. Because there was **no exceedance**, the exceedance classification of §I.8.3 (PRECISION-CONSISTENT / -LIMITED / UNEXPLAINED) was never invoked; the `tier3_adjudication` entries say "CONSISTENT WITH THE EXACT LATTICE SOLUTION AT 1e-10", which is a statement about ε ≤ 10⁻¹⁰, not a precision classification.
Latent weaknesses of the driver's *overall-status logic* found during this audit (none of them was triggered in this campaign, because every gate row passed, all 26 M4 rows had `gates_passed = True` and all 30 decomposition rows were evaluable): (i) the hard-fail scan reads only the top-level `status` of each criterion, so a failed XP gate would propagate to the overall status only through V1/G0, not through the nested `XP_gates` entry; (ii) an M4 solve at 1/1024 or 1/2048 whose X4/X5/X7 failed would be dropped from the adjudication without turning the overall status into FAIL; (iii) `NOT EVALUABLE` decomposition rows do not enter the overall status. These should be fixed in any future campaign that reuses the driver; they do not affect this one and were not changed here (no retroactive change of the run logic).

**(c) The M1 growth factor ≈ 12.08 at `W200_z0.2@E88.51`.**
ε(M1, M4) = 1.509×10⁻¹², 1.823×10⁻¹¹, 2.133×10⁻¹¹ at 1/512, 1/1024, 1/2048: growth 12.08, then 1.17 (geometric mean over the two halvings 3.76).
* **Rule status.** The operational factor-10 rule (§I.8.3 steps 4–6) is a classification of a value that *exceeds* 10⁻¹⁰. M1 never exceeds 10⁻¹⁰ (maximum 2.13×10⁻¹¹), so the rule was not applicable and no classification depends on it. For completeness: if it had been applied, `precision_consistency` would return not-passed (12.08 > 10) and an exceedance would have been labelled UNEXPLAINED. It is the only growth factor > 10 among the 30 (point, method) factors (range 1.17 – 12.08, median 3.8).
* **What the data suggest (not prove).** Normalised by the B1-type scale (N−1)ε_mach/sin(k dx) used for the ρ diagnostic, M1's error is 0.081, 0.244, 0.071 at the three resolutions, whereas M2 gives 0.294, 0.339, 0.296 and M3 0.312, 0.343, 0.258. The B3X calibration range at 1/512 was 0.081 – 0.357, so M1's 1/512 value sits exactly at the **lower end** of the calibration range (and 1/2048 slightly below it), while the 1/1024 value is mid-range. The large factor therefore comes from an unusually *small* baseline at 1/512 for a statistic that is a maximum over matrix elements of accumulated round-off, not from a systematic growth of M1's error. M1's error at this point is the *smallest* of the three double methods at 1/2048 (2.1×10⁻¹¹ vs 8.9×10⁻¹¹, 7.7×10⁻¹¹).
* **Limits.** M1 has no explicit matrix, so no decomposition and no mechanism claim is available for it (protocol). The reading above is a statistical consistency argument based on the normalised scale; it is not a measured separation of error sources. Nothing in the data indicates an implementation defect: M1 agrees with M4 to 2.1×10⁻¹¹, with M2/M3 to ≤ 9.1×10⁻¹¹, and T1 agrees to 10⁻¹¹.
* **Recommendation (for a future protocol revision, not retroactive).** A growth-factor test between only two small, noisy numbers is fragile; a test on the normalised error or on a fit over ≥ 3 levels would be more robust. The predeclared rule and thresholds are left unchanged here.
