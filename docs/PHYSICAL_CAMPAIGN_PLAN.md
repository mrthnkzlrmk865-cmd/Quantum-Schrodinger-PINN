# Physical campaign: proposal, revision 3 (NOT frozen, NOT hashed, NOT run)

Status: **revised proposal for review.** Revision 1 (commit `e51b973`) was approved in direction; revision 2 resolves the two statistical issues raised in review (primary/secondary observables; bootstrap resolution) and clarifies the control references. No campaign code exists, no simulation was run, no number below is a result, and this protocol has **not** been frozen or hashed. Cost figures are planning estimates from B3Y timings. Arithmetic shown (sample-size tables, χ² factors, grid counts) is pure calculation, not data.

## Changes in revision 3 (with respect to revision 2)

Approved in the review of revision 2 and not reopened: τ̄ primary / T̄1 secondary with the LOW PRECISION rules; the single 68-test Holm family; the two-stage 10⁵ / 10⁶ bootstrap with an independent Stage-2 stream and conservative handling of residual ambiguity; the per-control h definitions and TOST wording; separate control evaluation for τ̄ and T̄1.

| # | item | revision 2 | revision 3 | reason |
|---|---|---|---|---|
| R3-1 | Monte Carlo intervals used for the Holm decision | 99 % Clopper–Pearson interval **per test** (marginal) | **Bonferroni-adjusted Clopper–Pearson intervals with simultaneous coverage ≥ 99 % over everything used in the decision, across both stages**: per-test two-sided level α₁ = 0.005/68 in Stage 1 and α₂ = 0.005/68 in Stage 2 (total MC error budget 0.01) | marginal 99 % intervals give no 99 % coverage over 68 tests (review) |
| R3-2 | definitions | p-value formula only | explicit definitions of the exceedance count K, the plus-one estimator, the probability parameter π covered by the interval, the comparison of interval endpoints with Holm (monotonicity argument), Stage-2 selection and use, and the conservative final rule (§5.3) | review |
| R3-3 | Stage-2 counts | "recomputed with B₂ = 10⁶" | Stage-2 counts **replace** (are not pooled with) the Stage-1 counts for the selected tests; no Stage 3; the selection rule is fixed in advance | avoid data-dependent pooling/extension |

## Changes in revision 2 (with respect to revision 1)

| # | item | revision 1 | revision 2 | reason |
|---|---|---|---|---|
| R2-1 | observable hierarchy | τ̄ and T̄1 used side by side; δ defined for τ̄ only, Holm family over both | **τ̄ primary, T̄1 secondary** (§4.1). Sample size, δ and the gating C7 check are designed for τ̄ only. The 68-test Holm family is kept, with its justification and cost stated; T̄1 achieved precision reported separately, with predeclared interpretation rules for low precision | inconsistency between §3.1 and §5 |
| R2-2 | bootstrap count | 10⁴ resamples; claimed "adequate" for the smallest Holm threshold | **two-stage: B₁ = 10⁵ for every test and every CI; B₂ = 10⁶ for every test whose Holm decision is not Monte-Carlo stable at B₁**; Monte Carlo SE and a 99 % Clopper–Pearson interval reported for every p (marginal; superseded by R3-1); MC-ambiguous decisions after B₂ are *not resolved* (§5.2) | at 10⁴ the threshold 7.35×10⁻⁴ corresponds to ≈ 7 exceedances (MC SE of p ≈ 37 % relative); the "adequate" statement was wrong |
| R2-3 | control reference h | "h = 1.96 SE_main" without saying which estimate | **h is defined per control** (§6, table "reference estimate"): always the confirmatory estimate that is actually reported, at the reporting dx/basis/buffer, as 1.96 × its SE (not a bootstrap quantity, not a pilot quantity) | review item 6 |
| R2-4 | wording of the equivalence interval | "90 % CI … (TOST)" and "with 95 % confidence" | **TOST equivalence test at α = 0.05** (each one-sided test at 0.05), equivalently the 1 − 2α = 90 % two-sided interval inside the margin; not called a 95 % confidence interval | review item 6 |
| R2-5 | C1/C3/C4 for T̄1 | same gating as τ̄ | controls are evaluated for both observables; a failure on T̄1 blocks T̄1 claims for the affected cells **and** prevents overall PASS (fail-closed); it does not trigger new seeds | consistency with R2-1 and fail-closed design |

## Changes in revision 1 (with respect to revision 0; kept for the record)

| # | item | revision 0 | revision 1 | reason |
|---|---|---|---|---|
| 1 | grid | 4×4 grid with iso-line (reduced 3×3 offered) | **4×4 grid and W·ζ = 5 iso-line fixed**; 3×3 dropped | decision |
| 2 | sample-size rule | N = ⌈1.2(1.96 σ_max/δ)²⌉, clamp [100, 400]; "200 expected sufficient if σ̂ ≤ 0.25" | **N_req = ⌈(1.96 σ_UCB/δ)²⌉** with σ_UCB a χ²-based simultaneous upper confidence bound from a 40-seed pilot; clamp [100, 1000]; **200 is only a provisional planning figure**; behaviour above the cap predeclared (§3) | rev. 0 was inconsistent: σ = 0.25, δ = 0.02 gives 721 (601 without the 1.2 factor), not 200 |
| 3 | meaning of δ | undefined | **δ = target half-width of the 95 % CI of a cell mean (or confirmatory contrast) of the band-averaged τ̄** | needed to make the rule interpretable |
| 4 | pilot | 20 seeds, plain max of sample SDs | **40 seeds**, per-cell/contrast χ² upper bounds at one-sided level 0.05/50; pilot never pooled with the confirmatory data | 20 seeds give ≈ 16 % relative SE of s (and worse for skewed data); max of noisy estimates is not a bound |
| 5 | buffer | 205/512 proposed; C4 compared 205 vs 103 | **205/512 throughout**; no realization is ever discarded for amplitude or lead-tail; C4 now tests 205 vs 410/512 | decision; the 103/512 geometry cannot be used as a control without selecting on tail size |
| 6 | dense-energy module | inside "optional module P", 29 energies | **removed from this campaign**; Appendix A, separate follow-up, **86** points for Δk = 0.25 over ±3σ_k | 6σ_k/0.25 = 84.9 intervals ⇒ 86 points, not 29 |
| 7 | resolution acceptance | \|Δ\| ≤ 0.25 × SE_main (point estimate vs SE) | **equivalence rule**: 90 % CI of the paired mean difference inside ±0.25·h, h = 95 % CI half-width of the main estimate; power check; applied to cell means and to confirmatory contrasts (§6) | the old rule ignored the uncertainty of the control itself and compared a bias to a standard error rather than to the reported interval |
| 8 | controls | acceptance given, consequences partly implicit | every control has acceptance **and** an explicit fail-closed consequence (§6 table) | decision |
| 9 | status logic | not specified | **fail-closed driver design** with mutation tests for the three Phase 1 gaps, before any campaign code (§8) | decision |
| 10 | Holm family | "the contrast list" | explicit count: 34 contrasts × 2 observables = **68** | needed for Holm and for bootstrap p-value resolution |
| 11 | scope | two-level | **three-level separation** stationary / asymptotic energy-integrated packet / time-resolved (§9) | decision |

Unchanged in revision 1: 10⁴ bootstrap resamples (superseded by R2-2 in revision 2); Holm correction; paired resampling over common seeds; SD, SE and CIs reported separately; the exclusion of Anderson-localization and continuum claims.

## 0. Question and what B3Y provides

**Question.** How do the disorder strength W and the correlation length ζ affect the transmission of a quantum wavepacket through a finite 2D hard-wall waveguide (L_y = 1, slab [1.9, 2.9], L_d = 1)?

B3Y established, for **one** realization (seed 0), 16 predeclared points and dx ≥ 1/2048, that the RGF agrees with two independent sparse solvers and with an exact *lattice* solution to ≤ 10⁻¹⁰ in S, and that the dx dependence of T1 is a property of the lattice problem. It says nothing about how transmission depends on W or ζ in general.

## 1. Design principles

1. **Ensemble, not anecdote:** cells are characterised by distributions over independent realizations; conclusions are ensemble statistics with uncertainty.
2. **W and ζ separated by design:** full 4×4 factorial, common random numbers, and the W·ζ = 5 iso-line.
3. **M1 (validated RGF) produces the data;** M3 is an independent per-realization cross-check on a subset; M4 is available only to adjudicate flagged realizations.
4. **Lattice results are lattice results.** Statements are for the finite-difference model (N_y = 63, N_q = 127, M = 512 modes) at the stated dx; resolution controls bound the lattice error of the *ensemble statistics*; no continuum claim.
5. **Predeclaration:** grid, seeds, observables, contrasts, acceptance and stop rules are fixed before the first confirmatory run; later additions are exploratory and labelled.
6. **Fail closed:** a mandatory numerical gate that fails, is missing or was not run prevents an overall PASS (§8).

## 2. Disorder-model facts used (from `disorder.py`)

* V = W √(2/M) Σ_j cos(q_j·r + φ_j) S(x), q_j ~ N(0, I/ζ²), φ_j ~ U[0, 2π), M = 512; **no per-realization renormalisation** (sample variance fluctuates around W²). Fewer independent correlation cells fit in the slab at larger ζ, so realization-to-realization variability is expected to grow with ζ (an expectation to be checked by the pilot, not a result).
* `SeedSequence([seed, 20241105])`: integer seeds are statistically independent.
* For one seed V is exactly linear in W and the ζ-change is a coordinate dilation of the same draws: **common seeds in every cell give exactly paired W-contrasts and draw-paired ζ-contrasts.**
* Window tails at the lead slices scale with the realized amplitude (B3Y seed 0: 9.57×10⁻⁷ vs the 10⁻⁶ guard at buffer 103/512). **Buffer 205/512 is used for the entire main ensemble.** No realization is dropped or replaced because of amplitude or lead-tail value (selection on large amplitude would bias the means toward weak disorder). The lead-tail guard (≤ 10⁻⁶) is still a mandatory gate: a violation is **not** filtered; it stops the affected stage (§8) and is reported.

## 3. Parameter grid and sample size

| factor | values |
|---|---|
| W | 25, 50, 100, 200 |
| ζ | 0.025, 0.05, 0.10, 0.20 |
| clean control | W = 0 (deterministic) |

16 cells. **Iso-line W·ζ = 5:** (200, 0.025), (100, 0.05), (50, 0.10), (25, 0.20) — S(q→0) = 2πW²ζ² constant while ζ changes by 8×. This is one way of separating the parameters, not a unique "pure correlation effect" (the Born rate also depends on kζ); full-factorial contrasts remain primary.

**Energies.** kx = 20 + jσ_k, j = −2…2, σ_k = 1/(√2 σ_x) = 3.54 (E = 88.5 … 371.4; 4 – 8 open channels), as in B3/B3Y. Per-energy statistics are exploratory. The **primary scalar** is the within-realization 5-energy band average τ̄ (and T̄1), a *descriptive* quantity — not a packet transmission (§9).

### 3.1 Meaning of δ and the sample-size rule

* **δ** is the target **half-width of the 95 % confidence interval** of (a) each cell mean of τ̄ and (b) each confirmatory contrast of τ̄ (differences of cell means, paired). Setting δ = 0.02 means: with the planned N, the 95 % CI of a cell mean (or contrast) of τ̄ is ±0.02 or narrower. It is a precision target for the *estimate*; it is not a threshold on any physical effect and not a tolerance on numerical error.
* Pilot SD estimates s are of the **SD across realizations of τ̄** (cell) or of the **paired difference** (contrast).
* **Rule.** N_req = ⌈(1.96 σ_UCB / δ)²⌉ with σ_UCB = s · √((n_p − 1)/χ²_{α,n_p−1}) the one-sided upper confidence bound of the SD (valid for near-normal data; transmission is skewed, so it is a *planning* bound, not a guarantee), α = 0.05/50 (Bonferroni over the 16 cells and 34 contrasts of one observable), n_p = 40 ⇒ σ_UCB ≈ 1.50 s. N_req is the maximum over all cells and contrasts. The final N = clamp(N_req, 100, N_max) with **N_max = 1000**.
* **Table (planning arithmetic).** N_req = (147 s)²: s = 0.10 → 217; 0.15 → 487; 0.20 → 865; s ≈ 0.215 → 1000; 0.25 → 1351. (Rev. 0 formula for reference: σ = 0.25, δ = 0.02 ⇒ 601 without the 1.2 factor, 721 with it. The absolute worst case for a variable in [0, 1], s = 0.5, would be ≈ 5400 even before the pilot inflation; it is not expected.)
* **200 is a provisional planning estimate only.** With s = 0.25 it would give a CI half-width of 0.035, not 0.02. N = 200 may be used only if the pilot supports N_req ≤ 200.
* **If N_req > 1000:** the target δ is **not** silently relaxed and N is not silently capped. The pre-registered outcome is: *PRECISION TARGET NOT MET*; the campaign is paused for a decision by the reviewer (options: raise N_max, relax δ explicitly, or run at N_max and report achieved half-widths). If run at N_max, every statement is made with the **achieved** CIs and no claim is attached to the unmet target.
* **Pilot vs confirmatory.** Pilot: seeds 10001 … 10040, all 16 cells, 5 energies, dx = 1/512, same pipeline. Pilot data are used **only** to compute N_req; they are discarded from inference, never pooled with the confirmatory ensemble (seeds 1 … N), and reported separately. The computed N_req, N and the pilot s-values are recorded in an addendum to this protocol *before* it is frozen.
* **Is a 40-seed pilot adequate? (and is 20?)** The relative standard error of an SD estimate is ≈ 1/√(2(n−1)): 16 % for n = 20, 11 % for n = 40 (normal theory; skewed/heavy-tailed τ̄ makes it worse). With 16 cells + 34 contrasts, the *maximum* of noisy estimates is biased upward, which is conservative but not a bound. 20 seeds would need the χ² inflation factor 1.87 (α = 0.001) instead of 1.50 and still be fragile; therefore the pilot is increased to 40. Residual risk, stated: if the true SD of the binding cell is above its σ_UCB, the achieved CI half-width will exceed δ; this is detected afterwards (achieved half-widths are always reported) and has the same consequence as the N_req > N_max case for the affected cells.
* **The rule is designed for the primary observable τ̄ only.** N is not increased to meet any precision target for the secondary observable T̄1 or for exploratory quantities (mode-resolved P_mn, per-energy values). The pilot also records the SDs of T̄1 cell means and contrasts and the T̄1 half-widths *forecast* for the chosen N; these appear in the addendum for information only and do not enter N_req.

## 4. Observables (per realization r and energy E_j; S from M1)

Probabilities only (independent of reference-plane convention).

| symbol | definition |
|---|---|
| N_open(E) | open channels (deterministic, 4 – 8) |
| T_tot | Σ_mn \|t_mn\|² = Tr t†t |
| τ | T_tot / N_open (primary after band averaging: τ̄) |
| T1 | Σ_m \|t_m1\|² (incident mode n = 1; band average T̄1) |
| R_tot, R1 | Tr r†r, Σ_m \|r_m1\|²; unitarity T_tot + R_tot = N_open |
| P_mn = \|t_mn\|² | full mode-resolved transmission matrix, stored |
| q_mix | 1 − P_11/T1; ⟨m⟩₁ = Σ_m m P_m1 / T1 |
| ln τ | spread/typical value, **not** a localization diagnostic |

Energies of one realization are correlated: standard errors are always across realizations.

### 4.1 Primary and secondary observables

* **Primary: τ̄** — the within-realization 5-energy band average of τ = T_tot/N_open. It uses all open channels (the most information per realization), is bounded in [0, 1] and comparable across energies with different N_open. The sample-size rule (§3.1), δ = 0.02 and the gating precision check C7 refer to τ̄ only.
* **Secondary: T̄1** — the band average of T1 (incident transverse mode n = 1, the mode of the experiment's initial packet). It is scientifically important for the wavepacket question but is a single-channel quantity and is expected to vary more between realizations; **it is not guaranteed to meet δ**.
* **Exploratory:** per-energy values, P_mn, q_mix, ⟨m⟩₁, R, ln τ statistics.

**Interpretation of T̄1 when its achieved uncertainty is larger** (predeclared):
1. Achieved 95 % CI half-widths of every T̄1 cell mean and contrast are reported in their own columns and summarised separately from τ̄; the report never quotes δ for T̄1.
2. A T̄1 cell or contrast with achieved half-width > 2δ = 0.04 is labelled **LOW PRECISION (secondary)**.
3. A T̄1 contrast can still be *resolved* under the Holm rule of §5 (a large effect can be detected with a wide interval); for LOW PRECISION contrasts the report states the resolved sign and the interval, not a magnitude statement.
4. A non-rejected T̄1 contrast is reported as "not resolved at the achieved precision" — never as "no effect" or "equal".
5. No qualitative T̄1 trend statement ("T1 decreases with …") is made from unresolved or LOW PRECISION contrasts, and τ̄ conclusions are never extended to T̄1 (or vice versa) by analogy.
6. T̄1 precision does not gate the campaign (C7 is τ̄-only); T̄1 resolution controls (C1, C3, C4) are gating (R2-5).

## 5. Statistics and analysis plan

### 5.1 Estimates and intervals

* Per (cell, E, observable): mean, sample SD (ddof = 1), **SE = SD/√N**, 95 % percentile-bootstrap CI over realizations (B₁ = 10⁵ resamples), median, IQR, min, max, skewness. SD (spread of the population of realizations), SE (uncertainty of the mean) and the CI are separate columns. Primary (τ̄) and secondary (T̄1) results are in separate tables.
### 5.2 Confirmatory tests, Holm correction and bootstrap resolution

* **Contrasts** (defined identically for O = τ̄ and O = T̄1; definitions unchanged from revision 1): C-W: Ō(2W, ζ) − Ō(W, ζ), 3 steps × 4 ζ = 12; C-ζ: Ō(W, 2ζ) − Ō(W, ζ), 12; C-int: difference of differences over adjacent (W, ζ) steps, 3 × 3 = 9; C-iso: least-squares slope of Ō against ln ζ along the W·ζ = 5 line, 1. 34 per observable, **68 in total**.
* **Why one 68-test Holm family.** Both τ̄ (primary) and T̄1 (secondary, the incident mode of the packet) are to support confirmatory statements, so the family-wise error rate is controlled over all of them at α = 0.05. Cost, stated: the primary tests face stricter thresholds than in a τ̄-only family (smallest 0.05/68 ≈ 7.35×10⁻⁴ instead of 0.05/34 ≈ 1.47×10⁻³). A gatekeeping order (τ̄ first, T̄1 only after) was considered and not chosen, because it would make secondary claims conditional on primary outcomes. If the reviewer prefers to protect the primary family's power, the alternative is a 34-test Holm family for τ̄ with T̄1 tests reported as exploratory; the choice must be fixed before freezing.
* **Resampling:** paired over common seeds (a seed is resampled jointly across all cells, energies and observables); the bootstrap RNG seeds are fixed and recorded.
* **Monte Carlo resolution (why 10⁴ is not enough):** the MC standard error of an estimated p is √(p(1−p)/B). At the smallest Holm threshold p ≈ 7.35×10⁻⁴: B = 10⁴ → ≈ 7 expected exceedances, MC SE ≈ 2.7×10⁻⁴ (37 % relative); B = 10⁵ → ≈ 74, 8.6×10⁻⁵ (12 %); B = 10⁶ → ≈ 735, 2.7×10⁻⁵ (3.7 %).
* The p-value definitions, the Monte Carlo intervals and the two-stage Holm procedure are specified in §5.3.
* **Cost:** B₂ = 10⁶ with N ≤ 1000 seeds and 32 cell means per replicate is a few × 10¹⁰ elementary operations, i.e. minutes, chunked (estimate). It involves no new simulation.
* **A trend is called resolved** only if (i) it belongs to the final conservative rejection set R_claim of §5.3, and (ii) the resolution controls of §6 are accepted for the cells involved and the observable concerned. Everything else is exploratory and labelled so.
* Descriptive reference only (no inference): τ = 1 (ballistic) and the diffusive-slab form ℓ/(L_d + ℓ) with ℓ from the Born rate of `disorder.py` (weak-disorder approximation), used to place cells.
* No sentence of the form "transmission decreases with W/ζ" is written unless the corresponding contrasts are resolved over the range it covers.

### 5.3 Bootstrap p-values, simultaneous Monte Carlo intervals and the two-stage Holm decision

All quantities in this subsection are computed separately for each of the m = 68 confirmatory tests, indexed i = 1 … 68.

**(a) Bootstrap replicate and exceedance count.** A replicate b resamples N seeds with replacement from the confirmatory seeds; the same resampled seed list is applied to all 16 cells, all energies and both observables (paired design), and θ*_{i,b} is the contrast i recomputed from the resampled realizations. With θ̂_i the contrast on the original data, the **exceedance count** is
K_i = #{ b = 1 … B : |θ*_{i,b} − θ̂_i| ≥ |θ̂_i| }  (null-centred, two-sided; ties count as exceedances, which is conservative).

**(b) p-value estimator.** The reported p-value is the **plus-one** Monte Carlo estimator p̂_i = (K_i + 1)/(B + 1). It is never 0 (minimum 1/(B+1)) and is the standard valid Monte Carlo p-value (Davison & Hinkley 1997; Phipson & Smyth 2010). Its MC standard error is reported as √(p̂_i(1−p̂_i)/B).

**(c) Parameter covered by the interval.** Conditional on the observed data, the B replicates are independent draws from the bootstrap distribution, so K_i ~ Binomial(B, π_i) exactly, with
π_i = P*( |θ* − θ̂_i| ≥ |θ̂_i| ),
the **ideal (B → ∞) bootstrap p-value** of test i. The Clopper–Pearson interval [π_i^lo, π_i^hi] covers π_i; it quantifies **Monte Carlo error only**, not the sampling uncertainty of the data (which is what the bootstrap test itself addresses) and not the asymptotic accuracy of the bootstrap.

**(d) Simultaneous coverage.** Each interval is the exact two-sided Clopper–Pearson interval at per-test level α_MC (lower limit Beta(α_MC/2; K, B−K+1), upper limit Beta(1−α_MC/2; K+1, B−K); lower = 0 if K = 0, upper = 1 if K = B). Clopper–Pearson coverage is ≥ 1 − α_MC for every π (exact-conservative). By the Bonferroni inequality, the probability that *any* of a set of intervals misses its π is at most the sum of their α_MC, **whatever the dependence between them** — this matters because all 68 counts of one stage come from the same replicates and are strongly dependent. The total MC error budget is split equally between the stages:
* Stage 1 (B₁ = 10⁵, all 68 tests): α₁ = 0.005/68 ≈ 7.35×10⁻⁵ per test;
* Stage 2 (B₂ = 10⁶, only the selected tests, independent random stream): α₂ = 0.005/68 per test.
Probability that any interval used in the final decision misses: ≤ 68·α₁ + |A|·α₂ ≤ 0.005 + 0.005 = 0.01. Hence **simultaneous coverage ≥ 99 % over the whole final box**, including the data-dependent selection of A (Stage-2 draws are independent of Stage 1, so the Stage-2 intervals have their stated coverage conditional on the selection; the Stage-1 term is bounded unconditionally). Planning arithmetic at the smallest Holm threshold (π ≈ 7.35×10⁻⁴): Stage 1 interval ≈ [4.5, 11.4]×10⁻⁴ (−39 %/+55 %), Stage 2 ≈ [6.3, 8.5]×10⁻⁴ (−14 %/+15 %); for K = 0 the upper limits are ≈ 1.0×10⁻⁴ (Stage 1) and 1.0×10⁻⁵ (Stage 2).

**(e) Comparison with Holm.** Let H(p) denote the set of tests rejected by Holm's step-down procedure at family-wise α = 0.05 applied to a vector p of 68 p-values. H is **coordinate-wise monotone**: lowering any p-value can only enlarge the rejection set (Holm is a monotone step-down procedure). The endpoint vectors are compared with Holm directly on the p-value scale:
* lower vector p^lo with p_i^lo = π_i^lo;
* upper vector p^hi with p_i^hi = max(π_i^hi, p̂_i) (the max guarantees that the conservative decision is never more liberal than the plus-one point estimate).
On the event that all intervals cover (probability ≥ 0.99), p^lo ≤ π ≤ p^hi componentwise, so by monotonicity
H(p^hi) ⊆ H(π) ⊆ H(p^lo).
Thus every test in H(p^hi) is rejected by the ideal bootstrap Holm procedure, and every test outside H(p^lo) is not, simultaneously for all 68 tests with MC confidence ≥ 99 %.

**(f) Two stages.**
1. Stage 1: compute K_i, p̂_i and the α₁-intervals for all 68 tests with B₁ = 10⁵. Form R₁ = H(p^hi) and the **MC-ambiguous set** A = H(p^lo) \ H(p^hi).
2. If A is empty, Stage 1 is final.
3. Stage 2: draw B₂ = 10⁶ new replicates from an **independent random stream** (a separate, recorded child of the bootstrap SeedSequence, never overlapping Stage 1). For every test in A, the Stage-2 count **replaces** its Stage-1 count (no pooling of the two stages); its p̂ and its α₂-interval are recomputed. Tests outside A keep their Stage-1 values. The selection rule for A is fixed here, before any data exist; no test outside A is re-evaluated and there is no Stage 3.
4. Re-form the final vectors (Stage-1 values for tests outside A, Stage-2 values for tests in A) and compute R_claim = H(p^hi_final) and A_final = H(p^lo_final) \ H(p^hi_final).

**(g) Conservative decision if Stage 2 remains ambiguous.** Only tests in **R_claim** may be reported as rejected (and only R_claim enters "resolved", together with the control conditions of §5.2). Tests in A_final are reported as **"NOT RESOLVED — Monte Carlo ambiguous"** with their counts, p̂, intervals and the Holm thresholds they straddle; they are never counted as rejections and never as "no effect". Because R_claim is computed from the upper vector of *all* 68 tests, a residual ambiguity in one test cannot make another test's claim more liberal. No further resamples are drawn after Stage 2.

**(h) Reported per test:** θ̂_i; unadjusted 95 % percentile CI (Stage-1 replicates); K_i, B and the stage used; p̂_i and its MC SE; the simultaneous MC interval [π_i^lo, π_i^hi] with its per-test level; Holm rank and threshold; the decision (rejected / not rejected / MC-ambiguous). The bootstrap seeds of both stages and the total MC error budget (0.01) are recorded in the metadata.

**(i) What the 99 % refers to.** It is the probability, over the Monte Carlo resampling only, that the reported Holm decisions agree with those of the ideal (B → ∞) bootstrap procedure in the conservative direction. The family-wise error control of the *scientific* conclusions is the Holm α = 0.05 (asymptotic validity of the paired bootstrap test assumed); the two are separate and both are stated.

## 6. Controls: acceptance and fail-closed consequence

**Equivalence rule for resolution/basis/buffer controls.** For a control with n_k seeds, observable O ∈ {τ̄, T̄1} and cell c, let d_r = O_r(control setting) − O_r(reporting setting) for the same seed r, with mean d̄ and SE_d = s_d/√n_k. The equivalence margin is ε_c = 0.25 h_c, where **h_c is the half-width 1.96 × SE of the confirmatory estimate that is actually reported** for that cell and observable (computed from the confirmatory ensemble at the reporting dx, basis and buffer; not from the pilot, not from the control seeds, not a bootstrap quantity, so the margin carries no Monte Carlo noise). For each control the reference estimate is named in the table below.
The control is **accepted** iff the **TOST equivalence test at α = 0.05** rejects both H₀⁺: E[d] ≥ +ε_c and H₀⁻: E[d] ≤ −ε_c, each with a one-sided t-test at level 0.05 (n_k − 1 degrees of freedom). Equivalently, the two-sided (1 − 2α) = 90 % interval d̄ ± t_{0.95, n_k−1} SE_d lies inside (−ε_c, +ε_c). (This 90 % interval is the TOST device; it is not a 95 % confidence interval.) Meaning: the hypothesis that the setting changes the reported mean by at least a quarter of its stated sampling half-width is rejected at α = 0.05; if the bias is at most ε_c, adding it linearly to the half-width inflates the error bar by at most 25 %, in quadrature by at most 3.1 %. Because ε_c shrinks as N grows, a larger ensemble demands a smaller lattice error — intended.
**Under-powered controls:** if t_{0.95, n_k−1} SE_d > ε_c/2 the control is under-powered; n_k is extended with the next seeds (as listed per control) up to 100; if still under-powered the control is NOT DEMONSTRATED and treated as a failure.
**Contrasts:** the same test is applied to the paired control difference of each confirmatory contrast whose cells are all covered by the control, with ε = 0.25 × 1.96 × SE of that reported contrast.

| id | control | size | acceptance (reference estimate for h) | consequence of failure (fail-closed) |
|---|---|---|---|---|
| C0 | per-realization gates: unitarity ≤ 10⁻¹⁰ (both products), reciprocity ≤ 10⁻⁸, open-channel count = analytic, lead-tail ≤ 10⁻⁶, V-hash recorded, bitwise-reproducible rerun of a 1 % random subset | all | all pass | the stage stops; overall cannot be PASS; the realization is **kept** in the data table with its flag, never silently excluded; investigation reported |
| C1 | dx = 1/1024 vs 1/512 | first 20 seeds, all 16 cells, 5 E | TOST equivalence; d_r = O(1/1024) − O(1/512); **h = 1.96 SE of the confirmatory cell mean / contrast at dx = 1/512** (the reported estimate); extension seeds 21 … 100 | failing cell(s): main data for the affected cells (and every contrast containing them) are not accepted at 1/512; they are recomputed at 1/1024 with the full ensemble and the control is repeated with 1/2048; if that also fails, results for those cells are reported as lattice-dependent, no extrapolated or continuum statement is made |
| C1b | dx = 1/2048 vs 1/1024 | first 5 seeds, 4 corner cells | TOST equivalence; d_r = O(1/2048) − O(1/1024) on the same seeds; **h = 1.96 SE of the confirmatory corner-cell estimate at the reporting dx** (1/512 normally; 1/1024 if C1 moved that cell); extension seeds 6 … 100. Tests whether the 1/1024 level used by C1 is itself converged relative to the reported precision | same chain; overall PASS requires either C1b accepted for the corner cells or an explicit "lattice-dependent" label on the affected cells |
| C2 | M3 vs M1 | first 10 seeds, all cells, 5 E, dx = 1/512 | max \|ΔS\| ≤ 10⁻¹⁰ (B3Y criterion, unchanged) | the realization is flagged SOLVER-DISAGREEMENT, kept in the data; adjudicated with M4 (exact lattice) as in B3Y; overall cannot be PASS while any flag is unadjudicated, and an adjudicated exceedance yields at best "COMPLETED WITH ADJUDICATED EXCEEDANCES", not PASS |
| C3 | N_y = 127, N_q = 254 vs 63/127 | ζ = 0.025 cells (W = 25, 200), first 10 seeds | TOST equivalence; d_r = O(N_y 127/N_q 254) − O(63/127); **h = 1.96 SE of the confirmatory estimate of that cell at the reporting basis**; extension seeds 11 … 100 | refine the basis for **all** cells at ζ ≤ 0.05 and re-run those cells; repeat C3 |
| C4 | buffer 410/512 vs 205/512 | 4 corner cells, the 10 seeds of the main ensemble with the largest realized max \|V\| in that cell (an adverse, deliberately selected subset), dx = 1/512 | TOST equivalence; d_r = O(410/512) − O(205/512); **h = 1.96 SE of the confirmatory corner-cell estimate at buffer 205/512**. Because the subset is adverse, acceptance bounds the buffer effect for the largest-amplitude realizations (conservative for the ensemble mean); extension: next-largest seeds up to 100 | adopt the longer buffer for the whole ensemble and recompute |
| C5 | B3Y reproduction: seed 0 in (W = 200; ζ = 0.05, 0.20) at the B3Y energies | 10 solves | \|ΔT1\| ≤ 10⁻⁶ (buffers differ; B3Y K1 showed ≤ 5×10⁻¹¹) | stop: pipeline differs from the validated one; no campaign data accepted until explained |
| C6 | clean guide W = 0 | 5 solves | T_nn = 1 within 10⁻¹² | stop |
| C7 | pilot-sufficiency check, **primary observable only**: achieved 95 % CI half-widths of the confirmatory τ̄ cell means and contrasts ≤ δ | all | every τ̄ half-width ≤ δ (T̄1 half-widths reported, not gating; §4.1) | PRECISION TARGET NOT MET for the failing cells: they are reported with the achieved intervals and no resolved claim relies on the unmet target; overall cannot be PASS |

Controls C1, C1b, C3 and C4 are evaluated separately for τ̄ and T̄1. A failure for T̄1 blocks T̄1 claims for the affected cells and prevents overall PASS; it does not trigger additional seeds (N is fixed by τ̄). The τ̄ consequences in the table apply as written.

The mode count M = 512 defines the ensemble (it is not a numerical control) and is stated, not tested.

## 7. Cost (estimates, to be measured in the pilot)

≈ 4 s per (cell, seed), including the potential evaluation and 5 energies at dx = 1/512 (M1 ≈ 0.5 s per energy at 923 slices).

| stage | (cell, seed) units | est. time, 1 core (4 cores) |
|---|---|---|
| pilot, 40 seeds × 16 cells | 640 | ≈ 0.7 h (0.2 h) |
| main, N = 200 (provisional) | 3,200 | ≈ 3.6 h (0.9 h) |
| main, N = 1000 (cap) | 16,000 | ≈ 18 h (4.5 h) |
| controls C1–C6 (≈ 2,700 solves + 20 × 16 re-solves at 1/1024) | — | ≈ 1.5 h |

## 8. Fail-closed driver (to be done before any campaign code; design only here)

The three gaps found in the Phase 1 audit of the B3Y driver (nested gate failures reaching the overall status only indirectly; M4/XP rows with failed gates dropped silently; NOT EVALUABLE decomposition rows ignored) are addressed in the **new** campaign driver, by design, without touching `validate_rgf_b3y.py`, its results or its protocol:

1. **Gate registry.** Every mandatory gate (C0 items per realization, C1–C7, G-type reproduction, platform gate) is a registry entry with a state in {PASS, FAIL, NOT RUN, NOT EVALUABLE}. Entries are generated from the predeclared protocol, not from whatever the run happened to compute.
2. **Overall PASS = every mandatory registry entry is PASS.** Any FAIL, NOT RUN, NOT EVALUABLE, missing entry, unknown status string, exception, or stage stopped early ⇒ overall is not PASS (FAIL / INCOMPLETE / NOT EVALUABLE). Statuses are flattened recursively; nesting cannot hide a failure.
3. **No silent dropping.** Rows whose gates failed stay in the output tables with their flag; counts of expected vs completed rows per stage are checked (expected counts come from the protocol); a mismatch fails the stage.
4. **Independent verifier.** A separate small script recomputes the overall status from the raw CSVs alone and must agree with `summary.json`; disagreement ⇒ FAIL.
5. **Tests before campaign code runs** (mutation tests, synthetic problem): inject (a) a failing nested gate, (b) a failed extended-precision/solver-cross-check row, (c) a NOT EVALUABLE row, (d) a missing row, (e) a missing/unknown status, (f) a stage stopped early, (g) a control that is under-powered — each must produce an overall status ≠ PASS. A regression test re-aggregates the *existing* B3Y CSV files with the new aggregator and must reproduce the B3Y outcome (PASS) without modifying them.

## 9. Scientific scope: three levels that must not be mixed

| level | what it is | provided by this campaign? |
|---|---|---|
| **1. Stationary scattering results** | energy-resolved \|t_mn(E)\|², T_tot, T1, R at 5 fixed energies, per realization and ensemble statistics | **Yes — this is the campaign's product.** |
| **2. Asymptotic energy-integrated packet transmission** | for a single incident mode, P_trans(t→∞) = ∫dk w(k) T1(E(k)), w(k) ∝ exp(−(k−k_x0)²σ_x²), E = (k²+π²)/2, applied per realization and then averaged | **No.** Needs T1(E) on a dense grid (Appendix A). The 5-energy band average τ̄ is not this quantity. |
| **3. Time-resolved wavepacket dynamics** | arrival time, delay, spreading, dwell/trapping, density at intermediate times | **No, and not obtainable from \|t\|².** It needs the phase of S (Wigner delay) and the full state, or the time-dependent solver. |

Statements to be repeated in any report: stationary probabilities carry **no** time-resolved information; any packet-level statement is asymptotic and energy-integrated; results are for the lattice model; the time-dependent solver is a different (spectral, finite-window, CAP) discretization.

Not claimed: Anderson localization from reduced transmission alone (reduced τ at fixed L_d = 1 and N_open = 4 – 8 is also compatible with ballistic-to-diffusive attenuation, reflection by correlated barriers, channel mixing; separating them needs length scaling in L_d and the behaviour of ⟨ln τ⟩ and its variance, a possible follow-up); any continuum result; generalisation beyond this grid, energy window, L_d = 1, L_y = 1, M = 512, isotropic Gaussian-correlated fields; decoherence/ENAQT; devices.

## Appendix A. Follow-up validation study (NOT part of this campaign; separate proposal)

Purpose: test whether level-2 (asymptotic energy-integrated) packet transmission can be reconstructed from stationary data and how it compares with the time-dependent solver. Corrected energy grid: ±3σ_k = ±10.6, i.e. a span of 6σ_k = 21.2; with Δk = 0.25 that is 21.2/0.25 = 84.9 ⇒ **85 intervals, 86 energy points** (not 29). A coarser Δk = 0.5 gives 43 points; the spacing must be justified from the measured structure of T1(E) (resonance widths), which the 5-energy campaign cannot provide — so the follow-up begins with a spacing test on a few realizations. Provisional size: 86 energies × 4 corner cells × 50 seeds = 17,200 solves ≈ 2.4 h (1 core, estimate); plus time-dependent comparisons on a few seeds (cost to be measured), with the limitations listed in §9 and the non-equivalence of the DST and finite-difference discretizations stated. Separate approval, separate protocol.

## Decisions requested (revision 3)

Approved in earlier reviews (recorded, not reopened): 4×4 grid with the W·ζ = 5 iso-line; 40-seed pilot and χ² upper-bound rule with disjoint pilot/confirmatory seeds, N_max = 1000 and pause-for-decision; buffer 205/512; C4 and C7; fail-closed driver as prerequisite; equivalence framework (0.25 h, TOST at α = 0.05, under-powered rule); τ̄ primary / T̄1 secondary with LOW PRECISION rules; the single 68-test Holm family; the 10⁵ / 10⁶ two-stage bootstrap with an independent Stage-2 stream; per-control h; separate τ̄/T̄1 control evaluation.

Open for approval: the Monte Carlo procedure of §5.3 — Bonferroni-adjusted Clopper–Pearson intervals (α₁ = α₂ = 0.005/68, total budget 0.01, simultaneous coverage ≥ 99 % across both stages), the definitions (a)–(c), the endpoint rule p^hi = max(π^hi, p̂), the replacement (not pooling) of Stage-1 counts for the selected tests, and the conservative rule (g).

After approval: the next implementation task is **only** the fail-closed driver and its mutation tests; the protocol is frozen and hashed after the pilot addendum, before any confirmatory run.
