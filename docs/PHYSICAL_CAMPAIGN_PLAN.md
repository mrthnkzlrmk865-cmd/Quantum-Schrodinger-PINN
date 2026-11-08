# Physical campaign: proposal, revision 1 (NOT approved, NOT frozen, NOT run)

Status: **revised proposal for review.** No campaign code exists, no simulation was run, no number below is a result, and this protocol has **not** been frozen or hashed. Cost figures are planning estimates from B3Y timings. Arithmetic shown (sample-size tables, χ² factors, grid counts) is pure calculation, not data.

## Changes with respect to revision 0 (every changed criterion)

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

Unchanged: 10⁴ bootstrap resamples; Holm correction; paired resampling over common seeds; SD, SE and CIs reported separately; the exclusion of Anderson-localization and continuum claims.

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
* **Precision of the other observables** (T̄1, mode-resolved P_mn, per-energy values) is **not** targeted; their achieved intervals are reported, and no contrast on them is called resolved if its CI does not exclude 0.

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

## 5. Statistics and analysis plan

* Per (cell, E, observable): mean, sample SD (ddof = 1), **SE = SD/√N**, 95 % percentile-bootstrap CI over realizations (**10⁴ resamples**), median, IQR, min, max, skewness. SD (spread of the population of realizations), SE (uncertainty of the mean) and the CI are separate columns.
* **Confirmatory family** (observables τ̄ and T̄1, 34 contrasts each ⇒ **68 tests**): C-W: τ̄(2W, ζ) − τ̄(W, ζ), 3 steps × 4 ζ = 12; C-ζ: τ̄(W, 2ζ) − τ̄(W, ζ), 12; C-int: difference of differences over adjacent (W, ζ) steps, 3 × 3 = 9; C-iso: least-squares slope of τ̄ against ln ζ along the W·ζ = 5 line, 1. Resampling is **paired over common seeds** (a seed is resampled jointly across all cells, energies and observables). Two-sided bootstrap p-values (+1 rule, resolution 1/(10⁴+1) ≈ 10⁻⁴, adequate for the smallest Holm threshold 0.05/68 ≈ 7×10⁻⁴), **Holm** correction at family-wise α = 0.05. CIs are unadjusted 95 % percentile intervals, labelled as such.
* **A trend is called resolved** only if its Holm-adjusted p < 0.05 **and** the resolution controls of §6 are accepted for the cells involved (so the effect cannot be a lattice artefact within the demonstrated bound). Everything else — per-energy behaviour, mode-resolved structure, ln τ statistics — is exploratory and labelled so.
* Descriptive reference only (no inference): τ = 1 (ballistic) and the diffusive-slab form ℓ/(L_d + ℓ) with ℓ from the Born rate of `disorder.py` (weak-disorder approximation), used to place cells.
* No sentence of the form "transmission decreases with W/ζ" is written unless the corresponding contrasts are resolved over the range it covers.

## 6. Controls: acceptance and fail-closed consequence

**Equivalence rule for resolution/basis/buffer controls (replaces the 0.25 × SE point-estimate rule).** For a control with n_k seeds (the first n_k seeds of the main ensemble, so the control is paired with the main data), observable O ∈ {τ̄, T̄1}, cell c: d_r = O_r(fine) − O_r(main), d̄ and SE_d = s_d/√n_k. Let h_c = 1.96·SE_main,c be the half-width of the main 95 % CI (SE_main,c measured in the confirmatory ensemble) and ε_c = 0.25 h_c. The control is **accepted** iff the 90 % CI of the mean difference, d̄ ± t_{0.95, n_k−1} SE_d (two one-sided 5 % tests, TOST), lies inside (−ε_c, +ε_c). Interpretation: the lattice bias of the reported mean is, with 95 % confidence, at most a quarter of its stated sampling uncertainty; its worst-case linear addition to the error bar is ≤ 25 %, in quadrature ≤ 3.1 %. Because ε_c shrinks as N grows, a larger ensemble demands a smaller lattice error — intended, since the lattice bias must stay below the statistical resolution actually achieved. **Power check:** if t·SE_d > ε_c/2, the control is *under-powered*; n_k is extended with the next seeds (cheap) up to 100; if still under-powered the control is NOT DEMONSTRATED, treated as a failure. The same rule is applied to each confirmatory contrast (paired control difference of the contrast, ε = 0.25 × half-width of the main contrast CI); for the 1/2048 control only contrasts whose cells are all in the corner set are tested.

| id | control | size | acceptance | consequence of failure (fail-closed) |
|---|---|---|---|---|
| C0 | per-realization gates: unitarity ≤ 10⁻¹⁰ (both products), reciprocity ≤ 10⁻⁸, open-channel count = analytic, lead-tail ≤ 10⁻⁶, V-hash recorded, bitwise-reproducible rerun of a 1 % random subset | all | all pass | the stage stops; overall cannot be PASS; the realization is **kept** in the data table with its flag, never silently excluded; investigation reported |
| C1 | dx = 1/1024 vs 1/512 | first 20 seeds, all 16 cells, 5 E | equivalence rule (above), cell means and contrasts | failing cell(s): main data for the affected cells (and every contrast containing them) are not accepted at 1/512; they are recomputed at 1/1024 with the full ensemble and the control is repeated with 1/2048; if that also fails, results for those cells are reported as lattice-dependent, no extrapolated or continuum statement is made |
| C1b | dx = 1/2048 vs 1/1024 | first 5 seeds, 4 corner cells | equivalence rule against the 1/1024 data | same chain; overall PASS requires either C1b accepted for the corner cells or an explicit "lattice-dependent" label on the affected cells |
| C2 | M3 vs M1 | first 10 seeds, all cells, 5 E, dx = 1/512 | max \|ΔS\| ≤ 10⁻¹⁰ (B3Y criterion, unchanged) | the realization is flagged SOLVER-DISAGREEMENT, kept in the data; adjudicated with M4 (exact lattice) as in B3Y; overall cannot be PASS while any flag is unadjudicated, and an adjudicated exceedance yields at best "COMPLETED WITH ADJUDICATED EXCEEDANCES", not PASS |
| C3 | N_y = 127, N_q = 254 vs 63/127 | ζ = 0.025 cells (W = 25, 200), 10 seeds | equivalence rule | refine the basis for **all** cells at ζ ≤ 0.05 and re-run those cells; repeat C3 |
| C4 | buffer 410/512 vs 205/512 | 4 corner cells, 10 seeds with the largest realized amplitude, dx = 1/512 | equivalence rule | adopt the longer buffer for the whole ensemble and recompute |
| C5 | B3Y reproduction: seed 0 in (W = 200; ζ = 0.05, 0.20) at the B3Y energies | 10 solves | \|ΔT1\| ≤ 10⁻⁶ (buffers differ; B3Y K1 showed ≤ 5×10⁻¹¹) | stop: pipeline differs from the validated one; no campaign data accepted until explained |
| C6 | clean guide W = 0 | 5 solves | T_nn = 1 within 10⁻¹² | stop |
| C7 | pilot-sufficiency check: achieved CI half-widths of the confirmatory τ̄ cell means and contrasts ≤ δ | all | every half-width ≤ δ | PRECISION TARGET NOT MET for the failing cells: they are reported with the achieved intervals and no resolved claim relies on the unmet target; overall cannot be PASS |

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

## Decisions requested (revision 1)

1. Approve the revised sample-size rule (§3.1) including the 40-seed pilot, α = 0.05/50, N_max = 1000 and the "N_req > N_max ⇒ pause for decision" consequence; N = 200 stays provisional.
2. Approve the equivalence-rule acceptance for controls (0.25 h, 90 % TOST CI, power check) and the consequence table (§6) — the 0.25 factor and 90 % level are new numbers flagged for review.
3. Approve the control C4 redefinition (410/512 vs 205/512) and the new C7.
4. Approve the fail-closed driver design (§8) as a prerequisite work item (code + mutation tests, no simulation).
5. Confirm the 68-test Holm family.
