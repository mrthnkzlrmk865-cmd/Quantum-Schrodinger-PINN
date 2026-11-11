# Design and predictions: W/Z coupling-scaling experiment v2 (independent; preregistration document)

Status: design only. Nothing in this document has been implemented or run.
- No code exists for it.
- No matrix of it has been generated.
- No previous file is modified.

Relation to the previous experiment (`PREDICTIONS_coupled_blocks.md`, design `7d69260`, implementation `63e046d`, results
`e3c67a6`):
- That experiment's decision, **DESIGN INSUFFICIENT**, is final and is not re-scored, re-fitted or reinterpreted here.
- Its post-hoc 3-point slopes (0.985–0.989) are **not** evidence for anything in this design and are not used to choose any
  threshold.
- Its matrices are not part of this experiment's primary or confirmatory data.

Standing status:
- B1 configuration F still FAILS the unchanged 1e-10 criterion.
- W/Z is experimental and NOT validated.
- No production change is authorized.

## 0. Inspection of the previous experiment: verified definitions and reproducibility gaps

Verified from `coupled_blocks.py` and `wz_mutants.py` at `63e046d`:

| item | exact definition found |
|---|---|
| coupling parameter s | `V_j(s) = v0 · s · M_j` with v0 = 4 (`coupled_blocks.V`). **s multiplies the entire potential**, including the diagonal terms μ_j I, α_j A, … (not only the off-diagonal coupling). Values used: s ∈ {1, 2⁻⁴, 2⁻⁸, 2⁻¹²} |
| W/Z input | Δ_j = (E·I − diag(E_perp) − V_j(s))·dx² in double; leads from `lead_from_delta(fl(E − E_perp)·dx²)`; dx = 1/8, t = 32, L_y = 1 |
| mutant discrepancy D_mut | **max over the four endpoint blocks** of ‖X_mut − X_ref‖_F / ‖Ĝ_full(s)‖₂, against C(R_W) |
| slope | `np.polyfit(ln s, ln D_mut, 1)` over the 4 scales (unweighted OLS) |
| "D ≫ τ_end" | D ≥ 100·τ_end |
| mutants | `wz_mutants.MUTANTS`: M1–M16 plus M4L and M4R. Pure order mutants: M1, M2, M3, M5, M6, M7, M8, M15. `MIN_N` per mutant; M11 only on N = 1 |

Gaps in the previous design (documented, not repaired retroactively):
- G1. The fitting procedure (unweighted OLS over all 4 scales) and the meaning of "≫" existed only in code, not in the
  predictions document.
- G2. D_mut took a **maximum over blocks**. The block that attains the maximum can switch between s values, which makes
  ln D_mut piecewise and its slope not the slope of any single analytic quantity.
- G3. D_mut was normalized by ‖Ĝ_full(s)‖₂, which itself depends on s, so the fitted quantity carried an s-dependent
  normalization.
- G4. The scale grid (4 points, spanning s = 1 to 2⁻¹²) mixed a possible strong-coupling regime with the weak-coupling regime
  in one fit, with no rule for detecting a crossover.
- G5. "Coupling strength" was not separated into commuting (diagonal) and non-commuting (off-diagonal) parts. Both scale
  with s.
- G6. No exact (noise-free) mutant discrepancy was computed, so algebraic signal and rounding noise could not be separated.

v2 keeps the **same definition of s** (an overall multiplier of a fixed potential matrix, V_j(s) = s·V_j(1)) and the same W/Z
definitions. It addresses G1–G6 by explicit, preregistered definitions.

## 1. Research questions and hypotheses

| RQ | question | how it is answered |
|---|---|---|
| RQ1 | Is the deviation from the expected law near s = 1 a nonlinear/crossover effect rather than numerical noise? | compare the exact (noise-free) local slope with the double-precision local slope at every grid point; locate the crossover with a preregistered definition (§5.4) |
| RQ2 | Does the asymptotic law emerge at weak coupling? | the exact asymptotic order q from a high-precision asymptote check (§5.2), and the double-precision slope in a preregistered weak-coupling window (§5.3) |
| RQ3 | Is the scaling stable across independently generated matrices? | the held-out set of 120 freshly generated configurations (§2); every meaningful pair must satisfy the rules |
| RQ4 | Can genuine ordering errors be distinguished from numerical noise? | exact signal D* vs measured noise η; SNR rules; diagonal and equivalent controls; perturbation-driven non-ordering control (M9) |
| RQ5 | Can the scaling criterion be fixed without tuning to the previous slopes? | every threshold is derived in §5–§6 from precision bounds or algebra; none is derived from any previously observed slope |

**Mathematical facts** (provable; checked numerically as validity conditions, not as hypotheses):
- F1. Every mutant output and every reference block is a rational function of s, analytic on the grid interval as long as
  all inverted matrices stay invertible (enforced by the manifest conditioning rule, §2.4).
- F2. At s = 0 the device is clean, so every Δ_j and Λ is diagonal. Every pure order mutant therefore coincides exactly with
  the correct algorithm: D*(0) = 0.
- F3. From F1 and F2, for a pure order mutant the exact discrepancy of a block satisfies D*(s) = ‖s^q C_q + O(s^{q+1})‖ for
  some integer q ≥ 1 (or D* ≡ 0), so the exact local slope tends to q as s → 0.
  - At first order the commutators reduce to [W⁰, X¹] with W⁰ diagonal. Only the off-diagonal entries of the first-order
    terms enter, with weights (w⁰_a − w⁰_b), and w⁰_a ≠ w⁰_b because the E_perp,n are distinct. So the first-order
    coefficient C_1 depends only on the off-diagonal part of V(1).
  - **Whether C_1 ≠ 0 (q = 1) for a given configuration is not guaranteed** (symmetry can cancel it) and is determined
    numerically per pair.
- F4. In class D (diagonal V at all s), pure order mutants satisfy D*(s) ≡ 0 for all s.
- F5. M12 is algebraically equivalent to the correct code (commuting functions of one matrix), so D*(s) ≡ 0 for all s.
- F6. M9 (stale P) is not an ordering mutant, but D*(0) = 0: at s = 0 the clean device has W_j = w* (the lead fixed point),
  so W_2 = W_1. It therefore scales with s in every class, **including D**. This is the control showing that s-scaling alone
  does not identify an ordering error.
- F7. The remaining structural mutants (M4, M4L, M4R, M10, M13, M14, M16) have D*(0) ≠ 0 generically (q = 0). M11 is off-path
  for every N ≥ 2.

**Empirical hypotheses** (can fail; tested on the held-out set):
- H1 (attribution): at every grid point where the signal is resolvable (SNR ≥ 1e3), the double-precision local slope equals
  the exact local slope to within 0.01. Any deviation from q near s = 1 is therefore algebraic (crossover), not numerical.
- H2 (asymptotic emergence): in the preregistered weak-coupling window, the double-precision OLS slope agrees with the
  exact asymptotic order q to within 0.02, and with the exact-data slope over the same points to within 0.005.
- H3 (accuracy across s): the correct W/Z implementation reproduces the exact Green blocks at rounding level (within its own
  bounds and ≤ 1e-10) at every grid point of every held-out configuration.
- H4 (distinguishability): in class D and for M12, the double-precision discrepancy of every pure order mutant stays within
  the noise tolerance at every s. For every meaningful order mutant, the double result tracks the exact signal within the
  noise tolerance.

**Hypothesis for description only** (not a decision criterion): q = 1 for most ordering pairs in classes R and H. q = 2 is
possible where a symmetry cancels C_1. The distribution of q is reported.

## 2. Configurations: generation, classes, exploratory vs held-out

### 2.1 Fixed lattice (unchanged physics conventions)

- dx = 1/8 (t = 32), L_y = 1.
- E_perp,n from `stationary.transverse_basis(1.0, N_y)` (doubles).
- Clean leads.
- W/Z and production definitions exactly as in `PREDICTIONS_coupled_blocks.md` §1 and §4.

### 2.2 Dimensions and classes

- N_y ∈ {2, 3, 4}.
- N ∈ {2, 3, 4, 6}. N = 1 is excluded: no order mutant is on its code path. N = 6 gives 4 loop iterations; all assembled
  matrices are ≤ 24 × 24.
- Classes:
  - **R:** real symmetric V_j.
  - **H:** complex Hermitian V_j. This breaks time-reversal symmetry, is outside `wz_rgf.py`'s documented real-symmetric
    scope, and is algebra-only for this project.
  - **D:** diagonal V_j (control).

### 2.3 Random matrix generation (exact dyadic values)

Root entropies are fixed now, before any generation:

| set | role | root entropy |
|---|---|---|
| X | exploratory | 202611110001 |
| C | confirmatory / held-out | 202611110002 |

Per configuration index k (0-based, in the manifest order of §2.5):
- The generator is `numpy.random.Generator(numpy.random.PCG64(numpy.random.SeedSequence(entropy=ROOT, spawn_key=(k, attempt))))`.
- `attempt` = 0, 1, 2, … is incremented only by the deterministic rejection rule of §2.4.

The draw order is fixed and must not change:
1. `n_open = gen.integers(1, N_y + 1)`.
2. `u = gen.random()` (energy fraction; §2.4).
3. For j = 1..N, in order:
   - `R_j = gen.integers(-8, 9, size=(N_y, N_y))`;
   - `I_j = gen.integers(-8, 9, size=(N_y, N_y))` (drawn for every class, so the stream is class-independent; used only for
     class H).

Matrix construction (exact in double):
- B_j = (R_j + i·I_j·[class H])/8, so every entry is a multiple of 1/8 in [−1, 1].
- V_j(1) = v0 · (B_j + B_jᴴ)/2 with **v0 = 8**. Entries are multiples of 1/2 in [−8, 8].
- Class D: zero every off-diagonal entry of V_j(1) after symmetrization.
- Classes R and D: V_j(1) is real. Class H: Hermitian with a nonzero imaginary part.
- V_j(s) = s·V_j(1), with s from the dyadic grid (§3), so it is exact.

Why v0 = 8 (design rationale, not tuned to any result):
- At s = 1 it gives |V dx²| ≤ 1/8 entrywise and ‖V dx²‖₂ ≲ N_y/8 ≤ 0.5.
- That is comparable to the channel spacings of δ (≈ 0.05–1), so s = 1 lies in a strong-coupling regime. This is required
  to probe RQ1.

### 2.4 Energy rule and deterministic rejection

- Window: E ∈ (E_perp,n_open, E_upper), with E_upper = E_perp,n_open+1, or E_perp,1 + 64·(2 − δ_min) if n_open = N_y.
- E = round(64·(E_lo + u·(E_upper − E_lo)))/64 (a dyadic grid of 1/64).
- Accept only if every channel satisfies |δ_n| ≥ δ_min = 0.02 and δ_n ≤ 2 − δ_min, with δ_n = (E − E_perp,n)dx².
- **Conditioning rule (design-time, from direct inverses only, never from any recursion or solver).** On every grid point
  s, the following must all be ≤ 1e3:
  - cond₂(M(s)),
  - cond₂(M^(j)(s)) for j < N,
  - cond₂(I − W_j(s)), with W_j taken from the prefix inverses,
  - cond₂(Z_N(s)) = cond₂([(M/t)⁻¹]_NN⁻¹).
- On failure of any rule, increment `attempt` and redraw everything for that k.
- No selection on commutator size, slopes or mutant signals is allowed.

### 2.5 Configuration sets

Order of generation: for set X, then set C, the configuration index k enumerates N_y ∈ (2, 3, 4) → class ∈ (R, H, D) →
N ∈ (2, 3, 4, 6) → draw d.

| set | classes R, H | class D | configurations |
|---|---|---|---|
| X (exploratory) | 3 N_y × 2 classes × 4 N × 2 draws = 48 | 3 N_y × 4 N × 1 draw = 12 | 60 |
| C (held-out, confirmatory) | 3 × 2 × 4 × 4 draws = 96 | 3 × 4 × 2 draws = 24 | 120 |

There is also an optional legacy diagnostic set (G), excluded from every decision:
- the 4 previous scaling base configurations (e2b/R/N3, e2b/R/N4, e3b/R/N3, e3b/R/N4, with their original V_j(1) and
  energies) on the v2 grid;
- it serves RQ1 descriptively only, and **does not re-score or alter the previous experiment's decision**;
- it may be dropped at implementation time without consequence.

## 3. Coupling grid

- **Grid:** s ∈ {2⁻ᵏ : k = 0..40} ∪ {3·2⁻⁽ᵏ⁺²⁾ : k = 0..39}, i.e. 1, 0.75, 0.5, 0.375, 0.25, …, 2⁻⁴⁰.
  - That is 81 points, with alternating steps of 4/3 and 3/2 (about 6.6 points per decade).
  - Every point is dyadic, so s·V_j(1) is exact in double. This is why non-dyadic values such as 10⁻ᵏ or 2^(−k/2) are not
    used.
  - s = 1 is a diagnostic point.
- **Lower edge 2⁻⁴⁰ ≈ 9.1e-13 (justified by precision, not convenience):**
  - The double-precision noise level used for SNR is N_X(s) ≥ τ_abs ≥ 1e-13·‖Ĝ_full‖₂ (§4.3).
  - For an order-1 signal D* ≈ |C_1|·s with |C_1| ≤ 1 (the expected scale for dimensionless Ĝ with ‖Ĝ‖ = O(1) and
    ‖V dx²‖ ≤ 0.5; an assumption), SNR falls below 10 before s = 2⁻⁴⁰.
  - The grid therefore deliberately extends **past** the detectability limit, so that the transition into the noise floor
    is observed (RQ4), rather than ending where signals are still strong.
  - Going further down is pointless in double: also at 2⁻⁴⁰, the diagonal part of s·V dx² (~5e-14) is still represented
    in Δ_fl with relative rounding ~1e-3, which the represented lattice R_W accounts for exactly (§4.4).
- **Exact asymptote points (mp only, not double):** s ∈ {2⁻⁹⁰, 2⁻⁹⁵, 2⁻¹⁰⁰}, evaluated on the exact lattice L(s) at dps 120.
  - These give q and |C_q| without any double rounding.
  - At these s the double representation cannot carry the perturbation, so they are never used for double-precision data.
- **Crossover detection:** see §5.4. It uses the exact profile on the full grid. The grid density (≈ 6.6/decade) resolves a
  10 % deviation from first-order behaviour to within a factor 1.5 in s.
- **No grid choice uses previous slopes.** The previous grid points (1, 2⁻⁴, 2⁻⁸, 2⁻¹²) are contained in this grid only
  because the grid is dense.

## 4. Quantities, references and numerical controls

### 4.1 Solvers and exact computations per (configuration, s)

- **A:** `wz_rgf.wz_endpoint_blocks` and `wz_fisher_lee` (read-only).
  - Plus the instrumented copy `wz_mutants.wz_instrumented` (unchanged file from `63e046d`), which must be bit-identical
    at the four endpoint blocks.
  - Double mutants come from the same file, unmodified.
- **B:** production `stationary.scattering_matrix(E, dev, min_closed=0)` on a directly built `Device`, for an accuracy
  comparison only.
- **C (reference):**
  - The direct mpmath inverse of the full block-tridiagonal matrix, assembled from its definition. Σ_L and Σ_R are added
    as two separate operations at slices 1 and N.
  - Built for the represented lattices L(s), R_W(s) and R_orig(s) exactly as in `PREDICTIONS_coupled_blocks.md` §4, at
    dps 60.
  - Prefix references come from the same routine without Σ_R.
  - No recursion code is used for C.
- **D\* (exact mutant signal):** a **new mpmath transcription** of every mutant, run on R_W(s) at dps 60 (the same lattice
  the double code sees) and on L(s), compared with C of the same lattice.
  - Agreement of the transcription with the double mutants at SNR-adequate points is a validity condition (§4.5), which
    guards against transcription errors.
  - Residual risk: both are written by the same author (§9).

### 4.2 Fitted error quantity (fixes G2 and G3)

For each mutant m, a single designated block X_m is fixed in advance (the endpoint block it reaches first):

| mutants | X_m |
|---|---|
| M1, M3, M6, M15, M9 | Ĝ_N1 |
| M2, M7 | Ĝ_1N |
| M5, M8, M4, M4L, M4R, M13, M14, M16 | Ĝ_11 |
| M10 | Ĝ_NN |

- Fitted quantity: the **absolute** Frobenius discrepancy in dimensionless Ĝ units,
  D_m(s) = ‖X_m^mut,double(s) − X_m^ref(s)‖_F, and its exact counterpart D*_m(s) = ‖X_m^mut,mp(s) − X_m^ref(s)‖_F.
- Absolute rather than relative, because ‖X_ref(s)‖ depends on s and would add a spurious s-dependence. The asymptotic order
  is unaffected either way, since ‖X_ref(s)‖ → ‖X_ref(0)‖ > 0.
- No maximum over blocks enters any fit. All four blocks are recorded and reported as secondary data.

### 4.3 Noise estimation and SNR

- **Observed noise:** η_X(s) = ‖X^correct,double(s) − X^ref(s)‖_F (the correct W/Z code, same block).
- **Conservative noise level:** N_X(s) = max(η_X(s), τ_abs(s), ε·‖X_ref(s)‖_F), where:
  - τ_abs = τ_end·‖Ĝ_full‖₂;
  - τ_end is the endpoint bound of `PREDICTIONS_coupled_blocks.md` §5 (with c = 10, n = N·N_y, K_Schur and g from the
    references at that s).
- **Exact SNR:** SNR*(s) = D*_m(s) / N_X(s).
- Thresholds (fixed):

| SNR* | meaning |
|---|---|
| ≥ 10 | detection-eligible: a detection claim may be made |
| ≥ 1e3 | regression-usable: relative noise ≤ 1e-3, so the noise-induced error of a centred local slope is ≤ 2·1e-3/ln 2 ≈ 3e-3 |
| < 10 | **below the noise floor: no claim of any kind. This is never interpreted as "harmless" or "not a bug"** |

- **Detection** at s means D_m(s) ≥ 10·N_X(s) **and** |D_m(s) − D*_m(s)| ≤ 10·N_X(s).

### 4.4 Represented-lattice controls

1. Off-diagonal entries of Δ_fl(s) are bitwise equal to −s·V_off(1)·dx² (exact).
2. The diagonal representation error of the s-dependent part is recorded at every s.
3. e_glob(C(R_W) − C(L)) and e_glob(C(R_orig) − C(L)) are ≤ 10·cond₂(M)·ε, as before.
4. |D*_{R_W}(s) − D*_L(s)| / D*_L(s) ≤ 1e-6 wherever D*_L ≥ 1e-30.
   - Justification: the lattices differ by O(ε) in commuting diagonal parts, which perturb C_1 only by relative O(ε). The
     diagonal rounding of the s-part enters at O(s²) (§1, F3).
   - Violation → INVALID.
5. At the mp asymptote points the lattice is L(s), stated explicitly. No double data exist there.

### 4.5 Reference checks and precision escalation

- (a) ‖M·(Ĝ/t) − I‖_max ≤ 1e-55 at every grid point (full and prefix systems).
- (b) Agreement at dps 100 within 1e-55 at every 4th grid point (k ≡ 0 mod 4), at s = 1, at s = 2⁻⁴⁰, and at every point of
  the duplicate set.
- (c) Column-wise `mp.lu_solve` vs `mp.inverse` within 1e-55 at s ∈ {1, 2⁻¹², 2⁻⁴⁰}.
- (d) mp asymptote points at dps 120 are checked against dps 160 within 1e-90 (relative to D*).
- Escalation: if (a) or (b) fails, recompute that point at dps 120. If it still fails, the run is **INVALID**. No silent
  repair.
- Analytic-identity checks (any violation → INVALID):
  - D*(s) ≤ 10^(−(dps−10)) relative to ‖Ĝ‖ for every pure order mutant in class D, for M12 everywhere, and for M5 at N = 3
    and M8 at N = 2 (F4, F5, previous §6 structural identities);
  - D*(s) → 0 as s → 0 for pure order mutants and M9 (from the asymptote points: D*(2⁻¹⁰⁰) ≤ 2⁻⁹⁰·‖Ĝ‖);
  - the prefix and closing identities of the previous §1 hold at ≤ 1e-50.
- Apparatus consistency: at every point with SNR* ≥ 10, |D_m − D*_m| ≤ 10·N_X(s). Otherwise the double mutant and the mp
  transcription disagree → **INVALID**.
- The instrumented copy is bit-identical to `wz_endpoint_blocks` at every execution (else INVALID).
- **Determinism:** the first 5 configurations of each set (X and C) are run twice over the full grid, single-threaded BLAS
  (`OPENBLAS_NUM_THREADS=1`), and must be bit-identical (else INVALID).

## 5. Scaling analysis (all procedures fixed before data collection)

### 5.1 Meaningfulness classification per (configuration, mutant)

Classification is by exact rules only, never by double results:
- off_path: N < `wz_mutants.MIN_N[m]` (and M11 always).
- identical_by_algebra: M5 at N = 3, and M8 at N = 2.
- equivalent: M12.
- class_D_order: a pure order mutant in class D.
- **ordering-meaningful:** a pure order mutant in class R or H, on path, not identical, with D*(2⁻¹⁰⁰) > 0 (beyond the
  identity threshold of §4.5).
- perturbation-control: M9 (all classes).
- structural: M4, M4L, M4R, M10, M13, M14, M16.

Only ordering-meaningful pairs enter the scaling decision. Identical, equivalent and class_D_order pairs are **never** counted
as detections. Expected count of ordering-meaningful pairs (on-path, non-identical order mutants per N: N = 2 → 3, N = 3 → 7,
N = 4 → 8, N = 6 → 8): 26 per (N_y, class, draw) quadruple, i.e. **312 in set X and 624 in set C**, before exclusions.

### 5.2 Exact asymptotic order

- q̂ = ln(D*_L(2⁻⁹⁰)/D*_L(2⁻¹⁰⁰)) / ln(2¹⁰).
- It must lie within 1e-3 of an integer q ∈ {1, 2, 3}. Otherwise → **INVALID**: F3 requires an integer, so a non-integer
  indicates a reference or transcription error.
- |C_q| = D*_L(2⁻¹⁰⁰) / (2⁻¹⁰⁰)^q.
- For structural mutants, q = 0 is expected and |C_0| = lim D*.

### 5.3 Weak-coupling window and regression (fixes G1 and G4)

**Window (fixed):** W_asym = {grid s ≤ 2⁻¹²} ∩ {SNR*(s) ≥ 1e3}.
- The upper edge follows from the first-order remainder. With D* = |C_q| s^q (1 + ρ s + …), the local slope deviates from q
  by ≈ ρ s / (1 + ρ s).
- Assume ρ ≤ ρ_max = 25. This is a heuristic bound from ‖V(1) dx²‖₂ ≤ 0.5 times an amplification ≤ 50 consistent with the
  conditioning rule cond ≤ 1e3; it is not a theorem. Then the model error is ≤ 0.01 requires s ≤ 4e-4, and the first dyadic
  grid point below that is 2⁻¹² ≈ 2.4e-4.
- (That this edge coincides with the previous experiment's smallest scale is a consequence of the bound, not a choice made
  from its results.)
- The assumption is checked exactly per pair: window adequacy in the regression rules below.

**Regression (fixed):** unweighted OLS of y = ln D_m(s) on x = ln s over **all** points of W_asym, plus the same OLS on
D*_m over the same points:
- Validity requirements, all of which must hold, else the pair is not regression-valid:
  - ≥ 8 points;
  - x span ≥ ln(100) (two decades);
  - no point removed except by the preregistered SNR rule.
- Estimates: the double slope p̂ and the exact slope p̂*.
- **Window adequacy:** |p̂* − q| ≤ 0.02, and every exact centred local slope in W_asym is within 0.02 of q. Failure → pair
  "window-inadequate" (design insufficiency, not a W/Z result).
- **H2 criterion:** |p̂ − q| ≤ 0.02 and |p̂ − p̂*| ≤ 0.005.
  - Derivation: model error ≤ 0.01 (above), noise contribution ≤ 2·1e-3/ln(100) ≈ 4e-4, and a safety factor of about 2.
  - 0.005 ≫ 4e-4 bounds the pure noise effect.
- **Uncertainty.** The residuals are deterministic rounding, not random noise, so no probabilistic confidence interval is
  claimed for a single pair. Reported per pair:
  - (i) a deterministic worst-case interval p̂ ± δp_det, with δp_det = Σ_i |w_i|·|ln(D_i/D*_i)| and w_i the OLS weights;
  - (ii) the OLS standard error, as a descriptor only.
- Across configurations: the median, min and max of p̂ − q per mutant and class, and a nonparametric bootstrap over
  configurations (10⁴ resamples, seed 202611110003) for the median, as a descriptive 95 % percentile interval.
- **Outliers:** none are removed. A grid point is excluded from a fit only by the SNR rule. Any point failing a validity
  check makes the run INVALID; it is never silently dropped.

### 5.4 Local slopes, crossover and regime comparison (RQ1)

- Exact and double local slopes at interior grid points: σ(s_k) = ln(D(s_{k+1})/D(s_{k−1})) / ln(s_{k+1}/s_{k−1}).
- The local slope near s = 1 is the one-sided σ(1) = ln(D(1)/D(0.5))/ln 2 (exact and double).
- **H1 criterion:** |σ_double(s) − σ_exact(s)| ≤ 0.01 at every grid point where SNR* ≥ 1e3 for both neighbours, including
  s ∈ [2⁻⁴, 1].
- **Crossover** (primary, model-based, exact data): r(s) = D*_{R_W}(s) / (|C_q| s^q) − 1. The crossover scale is
  s_× = the largest grid s such that |r(s')| ≤ 0.1 for every grid s' ≤ s. Reported per pair; s_× = 1 means no crossover in
  range.
- **Secondary crossover:** the largest s with |σ_exact(s') − q| ≤ 0.05 for all s' ≤ s.
- **Regime model comparison** (descriptive; fixed now), on the exact profile over the full grid:
  - M_a, pure power law: ln D* = a + q ln s.
  - M_b, first-order-corrected: ln D* = a + q ln s + ln|1 + b s|.
  - M_c, second-order-corrected: ln D* = a + q ln s + ln|1 + b s + c s²|.
  - Fit by least squares in ln D*. Select the simplest model whose maximum absolute log residual is ≤ 0.01 (a 1 % relative
    description error).
  - If even M_c exceeds 0.01, label "strongly nonlinear in range".
  - No model choice affects any PASS/FAIL rule. It describes where the strong-coupling regime lies.
- **Answer rule for RQ1:**
  - If H1 holds, any deviation of σ from q near s = 1 is attributed to the algebra (a crossover), because the exact profile
    contains no rounding.
  - If H1 fails, the deviation is at least partly numerical.

### 5.5 Controls in the scaling analysis

| control | requirement |
|---|---|
| class_D_order and M12 | D*(s) ≡ 0 (validity); double D_m ≤ 10·N_X at every s (H4) |
| M9 (perturbation control) | q from the asymptote reported; expected q ≥ 1 in all classes including D. This demonstrates that "D ∝ s" is not specific to ordering errors |
| structural mutants | detected (SNR* ≥ 10 and detection) at every s (expected q = 0) |
| noise-trend control | OLS slope of ln max(η_X, ε‖X_ref‖_F) vs ln s over W_asym in [−0.25, 0.25] for the correct code, for every held-out configuration and designated block. Otherwise the noise model is inadequate → design insufficient |

## 6. Preregistered decision rules (precedence: INVALID > FAIL > DESIGN INSUFFICIENT > PASS)

The decision is made on **set C (held-out) only**. Set X results are reported but cannot produce PASS.

**INVALID** if any of the following holds, in any set:
- a reference check, escalation, analytic identity, represented-lattice control (§4.4–§4.5), apparatus-consistency check,
  instrumented-copy identity or determinism check fails;
- the manifest hash at run time differs from the hash committed before the run (§8);
- q̂ is non-integer (§5.2).

**FAIL** if, on set C, any of the following holds:
- **H3:**
  - a correct-code A endpoint, prefix, W_j or P_j error, or ΔS, exceeds 1e-10 (where τ ≤ 1e-11);
  - **or** exceeds its own §5 bound of `PREDICTIONS_coupled_blocks.md`, applied at each s. (Stricter than the previous
    design, where the latter was only "anomalous". This is a deliberate tightening, not a weakening.)
- **H1:** any ordering-meaningful pair violates the local-slope agreement.
- **H2:** any regression-valid, window-adequate ordering-meaningful pair violates |p̂ − q| ≤ 0.02 or |p̂ − p̂*| ≤ 0.005.
- **H4:** any class_D_order or M12 pair has D_m > 10·N_X at some s.
- A structural mutant is not detected at some s with SNR* ≥ 10.

**DESIGN INSUFFICIENT** if no INVALID or FAIL condition holds and, on set C, any of the following holds:
- any ordering-meaningful pair is not regression-valid (too few SNR-adequate points or too little span);
- any ordering-meaningful pair is window-inadequate;
- the noise-trend control fails;
- fewer than 90 % of the expected 624 ordering-meaningful pairs exist after the exact classification (coverage);
- any of the 8 pure order mutants has fewer than 20 ordering-meaningful held-out pairs.

**PASS** only if none of the above applies. The PASS claim is limited to: "for 120 independently generated held-out small
coupled devices (N_y ≤ 4, N ≤ 6, equal lead widths, clean leads, three potential classes), the correct W/Z recursion
reproduces the exact Green blocks at rounding level over 81 coupling values, and every meaningful ordering mutant's
discrepancy follows its exact algebraic scaling down to the resolvable noise floor, with the asymptotic order emerging in the
preregistered weak-coupling window."

Further interpretation rules:
- A mutant pair below the noise floor (SNR* < 10) at some s carries **no** information about harm at that s.
- Class-D and M12 non-detections are expected by algebra and are **not** detections.
- The exploratory phase never relaxes any threshold (§8.2).

## 7. Scope

This experiment **can** establish, for small deterministic coupled devices:
- whether the previous strong-coupling slope deviation is reproduced as an algebraic crossover on new matrices;
- whether the asymptotic order emerges;
- whether ordering errors are separable from rounding noise down to a quantified floor;
- whether the correct W/Z code is accurate across a dense coupling range.

It **cannot** establish, and will not claim:
- general validation of W/Z for arbitrary devices;
- long-device stability (N ≤ 6 only);
- near-threshold behaviour (|δ_n| ≥ 0.02);
- unequal lead widths;
- disordered-system validation (no disorder ensembles are included);
- any resolution of the B1 configuration F failure;
- any change to the previous experiment's DESIGN INSUFFICIENT decision.

## 8. Reproducibility and protocol

### 8.1 Files (planned; none exist yet)

All under `experiments/wz_rgf_diagnostic/`:
- `coupling_scaling_v2_manifest.py` (generator; §2.3–§2.4);
- `coupling_scaling_v2.py` (driver);
- `wz_mutants_mp.py` (mp transcription of the mutants);
- outputs in a new directory `results_coupling_scaling_v2/`:
  - `manifest_X.json`, `manifest_C.json`: per configuration, the set, k, attempt, N_y, class, N, E, n_open, the integer
    arrays R_j and I_j, V_j(1) as exact rationals, the rejection log, and the numpy/python versions;
  - `manifest.sha256`;
  - `results_X.json`, `results_C.json`;
  - `points_X.csv`, `points_C.csv`: one row per (configuration, s), with every accuracy metric of A and B, conditioning and
    κ values, the reference checks, and per mutant D_m, D*_m (R_W and L) for all four blocks, η_X, N_X, SNR* and the
    category;
  - `pairs_X.csv`, `pairs_C.csv`: one row per (configuration, mutant), with category, q̂, |C_q|, p̂, p̂*, δp_det, OLS SE,
    number of points, span, window adequacy, s_× (both definitions), the selected regime model, σ(1) exact and double, and
    the H1/H2/H4 flags;
  - `coupling_scaling_v2_report.md`;
  - `previous_files.sha256` (§8.3).

### 8.2 Phases

1. **This design is committed** (the only action of this step).
2. **Implementation and manifest freeze:**
   - Write the code, generate both manifests, and commit them together with `manifest.sha256` **before any solver or
     reference run**.
   - Manifest generation may compute only the conditioning rule of §2.4.
3. **Exploratory run (set X; optional legacy set G):**
   - Publish its report.
   - At most one amendment may be committed afterwards, before the confirmatory run. It may only:
     - fix implementation bugs;
     - address feasibility (runtime, memory);
     - add diagnostics.
   - It may **not** change any threshold, window, grid, SNR level, tolerance, coverage requirement or decision rule. Any
     such change requires a new experiment ID.
4. **Confirmatory run (set C), then the decision by §6.** No further amendment is allowed.

### 8.3 Preserving the previous experiment

- Before any run, record sha256 of every existing file under `experiments/wz_rgf_diagnostic/` (excluding `__pycache__`) in
  `previous_files.sha256`. Verify it after every phase.
- Any change → INVALID.
- New code imports `wz_rgf.py` and `wz_mutants.py` read-only and writes only to `results_coupling_scaling_v2/`.

### 8.4 Independent reproduction

1. Fresh clone at the recorded commit.
2. Set `OPENBLAS_NUM_THREADS=1`.
3. Run the generator and verify `manifest.sha256`. Because numpy `Generator.integers` streams are not guaranteed stable
   across numpy versions, the frozen manifest, not the RNG, is authoritative.
4. Run the driver.
5. Compare all numeric outputs, excluding timestamps and paths, against the committed results. Equality is expected bitwise
   on the same platform and numpy/BLAS build. On another build, double-precision values may differ at rounding level, and
   the decision must be unchanged.
6. Record python, numpy, scipy and mpmath versions and the BLAS configuration.

Estimated cost (uncertain; to be measured in the exploratory phase):
- 180 configurations × 81 points ≈ 14,600 point executions with mp references (≤ 24×24), plus mp mutant transcriptions.
- Plus 540 asymptote evaluations and about 810 duplicate executions.
- Roughly 2–6 h single-threaded.
- Parallelization is allowed only over configurations, with a fixed output order, single-threaded BLAS and the determinism
  check.

## 9. Assumptions and unresolved design questions

Assumptions, labelled:
- A1 (heuristic): ρ ≤ 25 for the window edge (§5.3). Checked exactly per pair; failure → design insufficient.
- A2 (heuristic): |C_1| ≤ 1, used only to argue that the grid extends past the noise floor. It is not used in any decision.
- A3 (heuristic): the τ bounds of the previous §5. They become hard FAIL criteria here, so a too-optimistic τ would produce
  FAIL rather than pass.

Open questions for review before implementation:
1. **Strictness:**
   - Should exceeding the heuristic τ bound (below 1e-10) be FAIL, as chosen here, or a separate "anomalous" outcome?
   - Should coverage require 100 % of ordering-meaningful pairs to be regression-valid? It currently does (any non-valid
     pair → DESIGN INSUFFICIENT). This makes DESIGN INSUFFICIENT more likely.
2. **mp transcription independence:** the transcription of the mutants is written by the same author as `wz_mutants.py`. A
   shared conceptual error would pass the apparatus check. An independent transcription, or a review, may be required.
3. **Designated blocks X_m:** the choice is structural. Alternatives (e.g. all four blocks fitted separately, with each
   required to pass) are stricter and would multiply the coverage requirements.
4. **v0 = 8 and the entry range [−8, 8]:** the strength is chosen so that s = 1 is strong coupling. A different strength
   changes where the crossover lies, but not the asymptotic order.
5. **Legacy set G:** include it (RQ1 diagnostic only) or drop it entirely, to avoid any appearance of revisiting the
   previous experiment.
6. **Runtime** may force a coarser dps-100 cross-check schedule. That would be a permitted feasibility amendment only if it
   does not touch any decision threshold.
7. **Production solver B** is included for accuracy only. Confirm whether its comparison should carry any decision weight.
   Currently it carries none.
8. **H3 at very small s:** the correct-code bound τ is evaluated per s. At s → 0 the device approaches the clean diagonal case
   tested before. No new risk is expected, but this is not proven.
