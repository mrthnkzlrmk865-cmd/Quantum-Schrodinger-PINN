# Fail-closed campaign status: validation contract notes (R6)

Status: **validation-contract documentation, not a protocol.** Nothing here changes the approved statistical procedure of
`docs/PHYSICAL_CAMPAIGN_PLAN.md` (68-test Holm family, B₁ = 10⁵, B₂ = 10⁶, Bonferroni Monte Carlo budget 0.01, Holm α = 0.05, sampling and
Monte Carlo interval definitions). The physical protocol is not frozen or hashed. No simulation, pilot or confirmatory run was made.
Implementation: `campaign_status.py` (aggregator), `verify_campaign_outputs.py` (independent verifier; imports neither the aggregator nor
the driver), `docs/physical_campaign_gates_draft.json` (DRAFT machine-readable transcription), tests in `tests/test_campaign_status.py`.

## 1. Bootstrap random-number streams

**Already specified by the plan** (section 5.2/5.3):
* one bootstrap replicate resamples the confirmatory seeds jointly for all cells, energies and both observables, so all 68 tests share the
  Stage-1 replicates — **one Stage-1 stream**;
* Stage 2 draws B₂ = 10⁶ new replicates "from an independent random stream (a separate, recorded child of the bootstrap SeedSequence, never
  overlapping Stage 1)" — **one Stage-2 stream, shared by all Stage-2 tests, derived as a child of the bootstrap SeedSequence**;
* "the bootstrap RNG seeds are fixed and recorded".

**Not specified (protocol decisions required; none is made here):**
1. the root entropy (seed) of the bootstrap SeedSequence and how it is fixed before data exist;
2. the bit generator (e.g. PCG64, Philox);
3. which SeedSequence node is the Stage-1 stream (the root, or a spawned child) and which spawn key is Stage 2;
4. the identifier format recording rule, entropy and spawn key;
5. the resampling call and the numpy version pin required for bit-for-bit reproducibility of the streams.

**Reporting contract (added schema; values are recorded, not prescribed).** Every row of `tests.csv` carries `bootstrap_rng_rule`,
`bootstrap_rng_stream_stage1`, `bootstrap_seed_stage1`, `bootstrap_rng_stream_stage2`, `bootstrap_seed_stage2`. Stage-1 tests carry
the Stage-1 stream and seed and **empty** Stage-2 fields; Stage-2 tests also carry the Stage-2 stream and seed. Gate
`PHYS.statistics.bootstrap_rng` (both implementations) requires: every field present and non-empty where required; seeds written as
non-negative decimal integers; exactly one rule, one Stage-1 (stream, seed) pair across all 68 rows, at most one Stage-2 pair across the
Stage-2 rows; Stage-2 stream identifier different from the Stage-1 identifier; if a rule is approved, the recorded rule must be it and the
rule's validator (if implemented) must confirm that the identifiers derive from the recorded seeds. Any violation is FAIL.
**While `bootstrap_rng` is null in the protocol the gate is NOT EVALUABLE** even when the metadata are structurally perfect, and an
approved rule without an implemented validator is also NOT EVALUABLE. The production registries of rule validators are empty.

**Seeds, streams, independence — three different claims.**
1. *Different recorded integer seeds* only say that two generator states were initialised from different numbers. Nothing follows about
   overlap or correlation of the sequences.
2. *Distinct streams* require a declared mechanism that guarantees distinct (ideally non-overlapping) subsequences: SeedSequence spawn keys,
   jump-ahead of a long-period generator, or disjoint counter/key spaces of a counter-based generator. The verifier can check only that
   the recorded identifiers are distinct and, once a rule and its validator exist, that they follow the rule.
3. *Statistical independence* of pseudo-random streams cannot be proved. The defensible statement is "independent by design of the
   declared mechanism" (plus, optionally, empirical test batteries). The gate therefore never reports independence as established;
   its PASS (only possible under an approved and implemented rule) means "metadata recorded, consistent and derived by the approved rule".

**Alternatives for review (none selected):** (A) root `SeedSequence(E)`; Stage 1 = child spawn key (0,), Stage 2 = child (1,); PCG64.
(B) root used directly for Stage 1; Stage 2 = `root.spawn(1)[0]`; PCG64 (literal reading of "child of the bootstrap SeedSequence").
(C) PCG64 for Stage 1 and the same generator `.jumped(1)` for Stage 2 (non-overlap by jump distance). (D) Philox with distinct key words
or counter blocks. For all: E fixed in the frozen protocol before data exist, as a recorded integer or derived from the frozen protocol
hash. The draft lists them in `bootstrap_rng_alternatives_for_review`.

## 2. Numerical verification tolerance

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

## 3. B3Y decomposition rows at dx = 1/512 (compatibility)

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
* **Decisions a reviewer must specify:** (1) the resolving comparison (which ε against M4, which limit); (2) the M4 gates the adjudicating
  solve must pass (X4/X5/X7, e_ind) and the columns that record them; (3) the outcome vocabulary; (4) how an adjudicated exceedance is
  reported in the campaign status beyond "not PASS".

## 5. Gates that are NOT EVALUABLE under the real draft (by design, until a protocol decision is made)

| gate | missing decision |
|---|---|
| `PHYS.C0.reproducibility` | selection algorithm and seed of the 1 % reproducibility subset |
| `PHYS.statistics.bootstrap_rng` | bootstrap RNG rule, root entropy, Stage-1/Stage-2 stream derivation, identifier format, numpy pin |
| `PHYS.C2.adjudication_resolved` | numerical C2 resolution criterion (NOT EVALUABLE only if a C2 row is flagged) |

Consequently no run under the current draft can reach overall PASS.
