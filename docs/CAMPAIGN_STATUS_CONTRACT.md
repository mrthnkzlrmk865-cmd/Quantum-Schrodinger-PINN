# Fail-closed campaign status: validation contract notes (R6, R7)

Status: **validation-contract documentation, not a frozen protocol.** R7 records the reviewer's decisions (bootstrap RNG option A, tolerance, C0 selection procedure); the two numbers they require (root entropy, C0 seed) are still unrecorded. Nothing here changes the approved statistical procedure of
`docs/PHYSICAL_CAMPAIGN_PLAN.md` (68-test Holm family, B₁ = 10⁵, B₂ = 10⁶, Bonferroni Monte Carlo budget 0.01, Holm α = 0.05, sampling and
Monte Carlo interval definitions). The physical protocol is not frozen or hashed. No simulation, pilot or confirmatory run was made.
Implementation: `campaign_status.py` (aggregator), `verify_campaign_outputs.py` (independent verifier; imports neither the aggregator nor
the driver), `docs/physical_campaign_gates_draft.json` (DRAFT machine-readable transcription), tests in `tests/test_campaign_status.py`.

## 1. Bootstrap random-number streams — R7: option A APPROVED, root entropy not yet recorded

**Approved mechanism (rule `SEEDSEQ-PCG64-SPAWN-A`; machine-readable in `bootstrap_rng` of the draft):**

| item | specification |
|---|---|
| root entropy | one non-negative integer in decimal digits (no sign / exponent / decimal point), recorded **before any physical campaign data are generated**. **Not supplied: `root_entropy = null`.** |
| SeedSequence | `root = numpy.random.SeedSequence(entropy=root_entropy)` |
| spawn procedure | `children = root.spawn(2)`; Stage 1 = `children[0]` (`spawn_key == (0,)`), Stage 2 = `children[1]` (`spawn_key == (1,)`); equivalently `SeedSequence(entropy=root_entropy, spawn_key=(i,))` (verified identical by test). The root's own stream is not used. |
| bit generator / initialisation | PCG64: `rng_i = numpy.random.Generator(numpy.random.PCG64(children[i]))` |
| resampling API | `rng_i.integers(low=0, high=N, size=N, dtype=numpy.int64, endpoint=False)`, one call per replicate, replicates b = 0 … B−1 in order (B₁ = 10⁵, B₂ = 10⁶); replicate b of a stage is applied to all cells, energies and observables, i.e. to all 68 tests (paired). A batched implementation is admissible only if it reproduces exactly this index sequence. |
| numpy pin | `numpy==2.5.3` (from `requirements-lock.txt`) |
| stream identifier | `SEEDSEQ-PCG64-SPAWN-A\|entropy=<root_entropy>\|spawn_key=(<i>,)\|bitgen=PCG64\|numpy=2.5.3`, i = 0 (Stage 1), 1 (Stage 2) |
| recorded per test (`tests.csv`) | rule, bit generator, numpy version, resampling API (identical in all 68 rows); Stage-1 identifier, seed (= root entropy) and probe; for Stage-2 tests also the Stage-2 identifier, seed (= root entropy) and probe; Stage-2 fields empty for Stage-1 tests |
| probe | sha256 of the first replicate's index array of the stage (`<i8` bytes) — allows anyone with the pinned numpy to confirm that the recorded identifiers really produce the streams |

**Independent validation.** The aggregator (`_rule_a_validator`, child via `SeedSequence(entropy, spawn_key=(i,))`, identifier by
f-string) and the verifier (`_rule_a_v`, child via `SeedSequence(E).spawn(2)[i]`, identifier from the draft's `stream_id_format`
template) are separately coded. Each checks: rule, PCG64, numpy version and API strings equal the protocol; both seeds equal the root
entropy; identifiers carry spawn keys (0,) and (1,) in the declared format; one Stage-1 stream for all 68 tests and one Stage-2 stream for
all Stage-2 tests, different identifiers; probes recomputed with numpy. If the installed numpy differs from the pin the probes cannot be
reproduced and the gate is NOT EVALUABLE. A malformed protocol root entropy makes the registry unbuildable (overall FAIL).

**Wording.** The Stage-1 and Stage-2 streams are *independent by design under the declared RNG mechanism* (distinct SeedSequence spawn keys
of one root). This is not a mathematical proof of statistical independence; different integer seeds alone would prove nothing, and
distinct identifiers prove only that distinct spawn keys were declared.

**Gate state now:** `PHYS.statistics.bootstrap_rng` = NOT EVALUABLE ("root SeedSequence entropy … has not been recorded"), even for
structurally perfect metadata. Superseded R6 text: the RNG rule is no longer undecided; the R6 alternatives B–D were not adopted.

## 2. Numerical verification tolerance — R7: accepted as stated (rtol = 1e-9, atol = 0)

**Where.** Before R6: `campaign_status.statistics_table_evaluator` used `|a−b| ≤ 1e-9·max(|a|, |b|, 1e-300)` with an extra `+1e-300` on
lower endpoints; `verify_campaign_outputs._stats_ok` used `math.isclose(rel_tol=1e-9, abs_tol=1e-300)`. Both compare the *reported*
p̂, MC SE, Clopper–Pearson lower/upper endpoints, p_lo_final, p_hi_final, Holm threshold and Holm-adjusted p with values recomputed from the
integer counts (K, B, stage). The 1e-300 floors were numerically meaningless (equivalent to atol = 0).

**Measurements (this audit; planning arithmetic, not campaign data).**
* scipy Clopper–Pearson endpoints vs exact binomial tails (mpmath, 60 digits), α = 0.005/68, B ∈ {10⁵, 10⁶}, K ∈ {1, 2, 74, 735}:
  relative error ≤ 7.6×10⁻¹⁴; endpoints near one are accurate to ≈ 5×10⁻¹² absolute.
* MC SE `sqrt(p̂(1−p̂)/B)` vs the exact rational value: ≤ 10⁻¹⁶ in the bulk, 2.5×10⁻¹¹ relative at K = B−1, B = 10⁶ (cancellation in
  1 − p̂); exactly 0 at K = B.
* smallest nonzero compared quantity: the CP lower endpoint at K = 1, B = 10⁶ ≈ 3.7×10⁻¹¹. At K = 0 a naive `beta.ppf` returns NaN; the
  endpoint is defined as exactly 0 (and 1 at K = B).

**Rule (both implementations, separately coded):** a reported value a and a recomputed value b agree iff both are finite and
|a − b| ≤ atol + rtol·max(|a|, |b|) with **rtol = 10⁻⁹, atol = 0** (unchanged rtol; atol made explicit). rtol has a ≥ 40× margin over the
largest legitimate difference measured and is ≥ 10⁴ above the CP solver error. A positive atol is rejected: every compared quantity is
either structurally exact or ≥ 3.7×10⁻¹¹, and e.g. atol = 10⁻¹⁵ would let a K = 1 lower endpoint be wrong by 2.7×10⁻⁵ relative.
Added **structural exactness**: K = 0 ⇒ lower endpoint and p_lo_final exactly 0.0; K = B ⇒ upper endpoint and p_hi_final exactly 1.0;
ordering 0 ≤ p_lo_final ≤ p̂ ≤ p_hi_final ≤ 1. **The tolerance never affects a decision:** the Stage-2 selection, both Holm sets, the ranks
and the expected labels are recomputed from the integer counts (R6 change; previously the sets were formed from the reported floats).

## 3. B3Y decomposition rows at dx = 1/512 (compatibility) — R7: retained unchanged

The historical protocol `docs/rgf_b3y_protocol.json` (not modified) requires the M2/M3 Tier-3 error decomposition for the five
`tier3_xp_1024_2048` points at dx = 1/1024 and 1/2048 (20 rows). The B3Y driver additionally wrote the decomposition at dx = 1/512 for the
same five points (10 rows), where a Tier-3 solve exists because all five are also in `tier3_xp_512_validation`. The protocol neither
requires nor prohibits these rows.
* **protocol-required:** (p, 1024 | 2048, M2 | M3) for p in `tier3_xp_1024_2048` — each exactly once;
* **compatibility (permitted, not required):** (p, 512, M2 | M3) only if p is in both `tier3_xp_1024_2048` and `tier3_xp_512_validation`
  **and** `tier3_exact_lattice.csv` holds exactly one record (p, 512) with `gates_passed == True`; at most once each;
* **unexpected (FAIL):** anything else — another dx (256, 4096, …), method M1, a point outside `tier3_xp_1024_2048` (even with a 1/512
  Tier-3 record), a duplicate, or a 1/512 row whose Tier-3 record is missing or failed its gates.
The rule is implemented independently in `b3y_optional_decomposition_keys` (aggregator) and in `verify_b3y` (verifier). The real B3Y result
tree has exactly 20 required + 10 compatibility rows and still re-aggregates to PASS. Uncertainty stated: whether the 1/512 rows were
*intended* by the B3Y protocol is not documented; they are accepted as compatibility rows because they are evaluated (gates, mirrors,
finite terms) like the required rows and are supported by a validated Tier-3 solve, not because the protocol required them.

## 4. C2 solver cross-check and adjudication

* **Flagged C2 row:** a row of `solver_cross_check.csv` (key W, ζ, seed, energy index) whose `max_abs_dS` is not a finite number ≤ 10⁻¹⁰.
  `PHYS.C2.solver_cross_check` is FAIL whenever any row is flagged — with or without adjudication.
* **Association:** an adjudication row belongs to the flagged row with the identical key. Missing, duplicated, unexpected (no flagged
  counterpart) or mis-associated rows fail.
* **Two gates (R6 split of the former `PHYS.C2.adjudication_evaluable`):**
  `PHYS.C2.adjudication_recorded` — *administrative completion*: exactly one record per flagged row, status ADJUDICATED (with finite
  `eps_M1`, `eps_M3`) or NOT EVALUABLE; it says that an adjudication was recorded, **not** that anything was resolved.
  `PHYS.C2.adjudication_resolved` — *scientific resolution*: requires an approved numerical criterion. **None is approved**
  (`C2_adjudication_criteria = null`), so this gate is NOT EVALUABLE whenever anything is flagged (PASS if nothing is flagged).
* **Never a PASS by adjudication:** even with an approved criterion that is satisfied, `PHYS.C2.solver_cross_check` remains FAIL, in line
  with the plan ("at best COMPLETED WITH ADJUDICATED EXCEEDANCES, not PASS").
* **R7: no criterion approved; `C2_adjudication_criteria = null` retained.** Decisions still required from a reviewer (also in
  `C2_adjudication_open_decisions` of the draft): (1) the exact M4 comparison (which of ε_M1 / ε_M3 against the exact-lattice S, which
  norm); (2) the numerical acceptance limit; (3) the M4 validation gates the adjudicating solve must pass (X4 / X5 / X7, e_ind limit);
  (4) the recorded columns for those gates and comparisons; (5) the outcome vocabulary; (6) the reporting of adjudicated exceedances
  (the plan already fixes: never PASS, at best "COMPLETED WITH ADJUDICATED EXCEEDANCES").

## 5. C0 reproducibility subset — R7: selection procedure specified, seed not yet recorded

| item | specification (`reproducibility_selection` of the draft) |
|---|---|
| algorithm | `SEEDSEQ-PCG64-PERMUTATION-C0` |
| seed | one non-negative decimal integer **dedicated to C0** (must differ from the bootstrap root entropy; a separate SeedSequence, no spawn key). **Not supplied: `seed = null`.** |
| generator | `numpy.random.Generator(numpy.random.PCG64(numpy.random.SeedSequence(entropy=seed)))`, numpy 2.5.3 |
| units and ordering | all (W, ζ, seed) units of the confirmatory ensemble (16 cells × N), sorted ascending by W, then ζ, then seed (numeric); the keys are unique, so there are no ties |
| size | n = ⌈0.01 · 16 · N⌉ in exact rational arithmetic, i.e. (16N + 99) // 100 (aggregator: `Fraction`; verifier: `Decimal` ceiling; checked equal for every N in [100, 1000]) |
| API | `perm = rng.permutation(M)`, M = 16N; `selected = sorted(perm[:n])`; subset = `[units[i] for i in selected]` |
| recorded per subset row | `selection_algorithm`, `selection_seed`, `selection_numpy_version` (must equal the protocol values) |

The aggregator and the verifier compute the subset independently (sorted tuples vs nested ascending loops) and require the observed rows to
equal it exactly; a third implementation in the tests agrees. A malformed or non-dedicated protocol seed makes the registry unbuildable
(overall FAIL); a different installed numpy gives NOT EVALUABLE. **Gate state now:** `PHYS.C0.reproducibility` = NOT EVALUABLE (seed not
recorded); structural errors (duplicates, units outside the ensemble, wrong count) still FAIL.

## 6. Gates that are NOT EVALUABLE under the real draft (by design, until a protocol decision is made)

| gate | missing decision |
|---|---|
| `PHYS.C0.reproducibility` | the dedicated C0 selection seed (procedure approved in R7) |
| `PHYS.statistics.bootstrap_rng` | the bootstrap root SeedSequence entropy (option A approved in R7) |
| `PHYS.C2.adjudication_resolved` | numerical C2 resolution criterion (NOT EVALUABLE only if a C2 row is flagged) |

Consequently no run under the current draft can reach overall PASS.
