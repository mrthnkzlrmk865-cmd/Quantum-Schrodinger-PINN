# Physical campaign: proposal (NOT approved, NOT run)

Status: **proposal for approval.** No code has been written for it, no number in this document is a result, and nothing in `results/` was created or changed by it. Every cost figure is an estimate from the measured B3Y timings (`results/rgf_b3y_validation/per_energy_comparison.csv`). Once approved, this document becomes the predeclared protocol (it will be frozen and hashed before any run, as for B3Y).

## 0. Question and what B3Y does and does not give

**Question.** How do the disorder strength W and the spatial correlation length ζ affect the transmission of a quantum wavepacket through a finite 2D hard-wall waveguide (L_y = 1, disordered slab [1.9, 2.9], L_d = 1)?

B3Y established, for **one** realization (seed 0), 16 predeclared points and dx ≥ 1/2048, that the RGF agrees with two independent sparse solvers and with an exact *lattice* solution to ≤ 10⁻¹⁰ in S, and that the dx dependence of T1 is a property of the lattice problem (order ≈ 2.0). It says nothing about how transmission depends on W or ζ in general: a single realization is one sample of a random object, and the 5-energy values of seed 0 differ by up to 0.5 in T1 between neighbouring energies.

## 1. Design principles

1. **Ensemble, not anecdote.** Each (W, ζ) cell is characterised by a distribution over independent realizations; conclusions are about ensemble statistics with uncertainty, never about seed 0.
2. **Separate W from ζ by design, not by narrative.** Full factorial grid, common random numbers, plus a line along which one spectral quantity is held fixed (§3).
3. **Only the validated solver for the bulk of the data (M1, RGF)**, with M3 as an independent per-realization cross-check on a subset; no new solver.
4. **Lattice results are lattice results.** All statements are about the finite-difference lattice model (N_y = 63, N_q = 127, M = 512 field modes) at the stated dx; the resolution controls quantify the *lattice error of the ensemble statistics*; they do not turn the lattice into a continuum result.
5. **Predeclare** grid, seeds, observables, contrasts, acceptance rules and stop rules before the first run; no cell, observable or energy is added after seeing data; exploratory analyses are labelled as such.

## 2. Disorder model facts that the design relies on (from `disorder.py`)

* V(x,y) = W √(2/M) Σ_j cos(q_j·r + φ_j) S(x), q_j ~ N(0, I/ζ²), φ_j ~ U[0, 2π), M = 512, slab window S(x) (edge 0.02). The field is **not renormalised per realization**: the sample variance fluctuates around W² and the number of independent correlation cells in the slab is ≈ (L_d L_y)/(2πζ²)-ish (large for ζ = 0.025, only a few tens for ζ = 0.2), so **realization-to-realization fluctuations grow with ζ**.
* Seeding: `SeedSequence([seed, 20241105])`, so integer seeds 1, 2, … give statistically independent realizations.
* For one seed, **V is exactly linear in W** and the ζ-change is a coordinate dilation of the same random draws. Using the same seeds in every cell (common random numbers) therefore makes W-contrasts paired *exactly* and ζ-contrasts paired in the draws; paired contrasts have much smaller variance than unpaired ones.
* Window tails at the lead slices scale with the realized amplitude. B3Y's worst lead-tail was 9.57×10⁻⁷ against the 10⁻⁶ guard (seed 0, buffer 103/512): with many seeds some realizations *will* violate the guard. **Proposal: use buffer 205/512 (0.4004) for the whole ensemble** (tails ~e⁻⁴⁰ smaller). B3Y K1 showed T1 insensitive to this buffer (≤ 4.7×10⁻¹¹ at five points; verified here only for seed 0, hence control C4 below). Dropping violating realizations would select on large amplitude and bias the mean; every realization is therefore kept, and any guard violation is a reported stop condition, not a silent filter.

## 3. Parameter grid (minimal but defensible)

| factor | values | rationale |
|---|---|---|
| W | 25, 50, 100, 200 | factor-2 spacing; 100 and 200 are the B3/B3X/B3Y values; 25, 50 reach weaker scattering |
| ζ | 0.025, 0.05, 0.10, 0.20 | factor-2 spacing; 0.05 and 0.20 are shared with B3/B3X/B3Y (0.02 → 0.025 so that W·ζ closes on an exact line) |
| clean control | W = 0 | deterministic: T_nn = 1; checks the pipeline end to end, 1 realization |

16 cells (4 × 4). **Iso-line:** W·ζ = 5 passes through (200, 0.025), (100, 0.05), (50, 0.10), (25, 0.20). Along it the long-wavelength spectral weight S(q→0) = 2πW²ζ² is constant while the correlation length changes by a factor 8. This is the "ζ at fixed spectral power" contrast; it is **one** way of separating the two parameters and is not a unique definition of a pure correlation effect (the Born rate also depends on kζ; an optional second iso-line of fixed Born rate can be computed from `disorder.py` at planning time). Full-factorial contrasts (W at fixed ζ, ζ at fixed W, interaction) remain the primary description.

**Energies.** The B3/B3Y set kx = 20 + jσ_k, j = −2…2, σ_k = 1/(√2 σ_x) = 3.54 (E = 88.5 … 371.4; 4 – 8 open channels). Per-energy statistics are reported separately; the unweighted 5-energy mean is a *descriptive band average*, not a packet transmission (§7).

**Realizations.** Main ensemble: seeds 1 … N_real, N_real = **200** per cell (proposed default). Justification step (pilot, predeclared): 20 *other* seeds (10001 … 10020, discarded afterwards) in all 16 cells give σ̂ of τ (below); N_real = clamp(⌈1.2 (1.96 σ̂_max/δ)²⌉, 100, 400) with δ = 0.02 on τ; 200 is expected to be sufficient if σ̂_max ≤ 0.25. Seed 0 is added as an extra, separately labelled sample (not pooled) for cross-checking against B3Y values in the shared cells (W = 200, ζ = 0.05 and 0.20).

## 4. Observables (per realization r, energy E_j; S from M1)

Probabilities only, so they are independent of the reference-plane convention.

| symbol | definition | note |
|---|---|---|
| N_open(E) | number of open channels | deterministic (4 – 8) |
| T_tot | Σ_mn \|t_mn\|² = Tr t†t | total transmission (dimensionless conductance) |
| τ | T_tot / N_open | mean transmission per incident channel; primary scalar |
| T1 | Σ_m \|t_m1\|² | transmission for the experiment's incident mode n = 1 |
| R_tot, R1 | Tr r†r, Σ_m \|r_m1\|² | with unitarity: T_tot + R_tot = N_open |
| P_mn = \|t_mn\|² | mode-resolved transmission matrix | stored in full (≤ 8 × 8) |
| q_mix | 1 − P_11 / T1 | fraction of mode-1 transmission leaving mode 1 (intermode mixing) |
| ⟨m⟩_1 | Σ_m m P_m1 / T1 | mean outgoing transverse mode index for n = 1 |
| ln τ | | "typical" value; spread of ln τ is reported, **not** used as a localization diagnostic |

Aggregations: (i) per (cell, E); (ii) within-realization band average over the 5 energies, then ensemble statistics across realizations. Energies of one realization are correlated, so the standard error is **always across realizations**, never pooled over energies.

## 5. Statistics and analysis plan (predeclared)

* Per (cell, E, observable): mean, sample SD (ddof = 1), **SE = SD/√N**, 95 % CI by percentile bootstrap over realizations (10⁴ resamples), median, IQR, min, max, skewness. No normality is assumed (transmission is bounded and skewed). SD describes the population of realizations; SE/CI describe the uncertainty of the mean — they are reported in separate columns and never conflated.
* **Confirmatory contrasts** (the only inferential claims; Holm-adjusted across this list; paired bootstrap resampling *seeds jointly across cells*):
  C-W: τ̄(2W, ζ) − τ̄(W, ζ), for the 3 W-steps × 4 ζ and for T1;
  C-ζ: τ̄(W, 2ζ) − τ̄(W, ζ), 3 ζ-steps × 4 W;
  C-int: difference of differences (is the W effect ζ-dependent?);
  C-iso: slope of τ̄ along the W·ζ = 5 line (ζ = 0.025 → 0.20).
  A trend is called *resolved* only if its CI excludes 0 **and** exceeds the lattice bias of control C1. Everything else (energy dependence, mode-resolved structure, ln τ statistics) is exploratory and labelled so.
* Descriptive reference (not a fit target, not a claim): ballistic limit τ = 1, and the diffusive-slab expression τ ≈ ℓ/(L_d + ℓ) with ℓ from the Born rate of `disorder.py` (valid only for weak disorder; used only to place the cells, not to infer a regime).
* Results are reported cell by cell and as tables/figures with CIs; no sentence of the form "T decreases with W" is written unless the corresponding C-W contrasts are resolved at all energies it is meant to cover.

## 6. Controls (selected, not exhaustive)

| id | control | size | acceptance (predeclared) |
|---|---|---|---|
| C0 | numerical gates of B3Y per realization: unitarity ≤ 10⁻¹⁰, reciprocity ≤ 10⁻⁸, open-channel count = analytic, lead-tail ≤ 10⁻⁶, V-sample hash recorded, determinism of reruns | all | any failure stops the stage and is reported |
| C1 | dx = 1/1024 for the first 20 seeds in all 16 cells (paired with dx = 1/512); also dx = 1/2048 for the first 5 seeds in the 4 corner cells | 1,600 + 100 solves | cell-wise paired mean change \|Δτ̄\|, \|ΔT̄1\| ≤ 0.25 × SE of the main estimate; else the main dx is raised to 1/1024 for that cell and the study re-run. Richardson values only if the observed order ∈ [1.8, 2.2], always labelled estimates |
| C2 | M3 (wave-function) vs M1 on the first 10 seeds, all cells, 5 energies, dx = 1/512 | 800 solves | max \|ΔS\| ≤ 10⁻¹⁰ as in B3Y; a failure is reported, not hidden |
| C3 | N_y = 127, N_q = 254 on the ζ = 0.025 cells (W = 25, 200), 10 seeds | 100 solves | same acceptance as C1 (ζ = 0.025 has the finest structure; B3X refined N_y only at seed 0) |
| C4 | buffer 205/512 vs 103/512 for the 10 seeds with the smallest lead-tail margin, dx = 1/512 | 50 solves | \|ΔT1\| ≤ 2×10⁻⁴ (B3Y K1 limit) |
| C5 | seed 0 in the shared cells (W = 200; ζ = 0.05, 0.20) reproduces the B3Y values | 10 solves | \|ΔT1\| ≤ 10⁻⁶ (buffers differ) |
| C6 | W = 0 clean guide | 5 solves | T_nn = 1 within 10⁻¹² |

The M = 512 mode count is part of the *definition* of the ensemble (CLT-Gaussian one-point statistics approximately); changing M defines a different ensemble and is not a numerical control. This is stated, not tested.

## 7. Can stationary S-matrix data support a wavepacket-level interpretation?

**What it can support.** The experiment's initial state is a Gaussian packet in transverse mode 1 fully inside the clean left lead (x₀ = 1.2, σ_x = 0.2, k_x0 = 20, momentum density ∝ exp(−(k−k_x0)²σ_x²)). For a linear, time-independent Hamiltonian the transmitted probability at t → ∞ is P_trans = ∫ dk w(k) T1(E(k)) for a single incident mode, with the weight w(k) fixed by the initial condition and E(k) = (k² + π²)/2 — **provided** T1(E) is known on a grid fine enough for the integral, the packet components with other energies/channels are accounted for (single-mode packet: no cross-channel interference terms), and one compares with the same quantity (late-time transmitted probability) in the time-dependent solver. This is an *asymptotic, energy-integrated* statement.

**Limitations (to be repeated in any report).**
1. **Energy sampling.** The 5-energy set samples T(E) at spacing ΔE ≈ 50–90, while a slab of length 1 carries structure on a scale of order π v/L ~ 60 and finer resonances; the 5-point band average is not a packet transmission. A packet-level statement needs a dense grid (proposal: Δk = 0.25 over ±3σ_k, 29 energies).
2. **Time-resolved content** (arrival time, delay, packet spreading, dwell/trapping, intermediate-time spatial density) is not contained in |t|²; it requires the energy derivative of the phase of S (Wigner delay) and the full state, which is outside the minimal plan.
3. **Window and absorber.** The time-dependent solver has finite time T, finite domain and CAPs; its "transmitted weight at t_end" differs from the t → ∞ scattering quantity by what is still inside the slab/not yet arrived. The two are comparable only if t_end is long enough, which must be shown.
4. **Different discretizations.** The S-matrix is for the finite-difference lattice (dx = 1/512…), the time-dependent code is a DST (spectral) solver: lattice dispersion and potential projection differ; part of any discrepancy is discretization, not physics.
5. **Ensemble level.** The packet formula is applied realization by realization, then averaged; the mean of a packet transmission is not the packet transmission at the mean T(E).
6. **Disorder frozen** (static), single-particle, no interaction/decoherence: nothing here concerns dephasing or environment-assisted transport.

**Optional module P (separate approval; not part of the minimal plan):** dense-energy grid (29 energies) for 4 cells (the corners of the grid) × 50 seeds; plus a comparison of P_trans (S-matrix) with the transmitted weight of the existing time-dependent solver on 3 seeds in 2 cells, reporting the difference with the limitations above. Outcome is a *validation of the interpretation*, not a new physical result.

## 8. What will not be claimed

* No Anderson-localization claim from reduced transmission. Reduced τ with fixed L_d = 1 and N_open = 4 – 8 is consistent with ballistic-to-diffusive attenuation, reflection by correlated barriers, channel mixing, or localization; telling them apart needs length scaling at fixed (W, ζ) (e.g. L_d = 0.5, 1, 2, 4, behaviour of ⟨ln τ⟩ and its variance versus L_d, comparison with the localization length, scaling with N_open). That is a possible follow-up campaign, not this one.
* No continuum claim: results are for the lattice model; controls bound the *lattice error of the ensemble statistics*.
* No generalisation beyond the grid, the energy window, L_d = 1, L_y = 1, M = 512, Gaussian-correlated isotropic fields.
* No decoherence/ENAQT statement; no device claim.

## 9. Cost and outputs (estimates)

M1 ≈ 0.55 s per solve at 923 slices (dx = 1/512, buffer 205/512), ≈ 1.2 s at 1/1024.

| stage | solves | est. time (1 core) |
|---|---|---|
| pilot (20 seeds) | 1,600 | ~15 min |
| main (16 cells × 200 seeds × 5 E) | 16,000 | ~2.5 h (≈ 40 min on 4 cores) |
| C1 / C2 / C3 / C4 / C5 / C6 | ≈ 2,700 | ~1.2 h |
| optional module P | ≈ 5,800 + time-dependent runs | ~1.5 h + TD cost (to be measured) |

Outputs (new directory `results/transport_ensemble/`, protected results untouched): `realizations.csv` (one row per seed, cell, energy: all observables, N_open, V-hash, gates), `mode_resolved.npz` (P_mn), `cell_statistics.csv`, `contrasts.csv`, `controls.csv`, `summary.json`, `run_metadata.json` (commit, protocol sha256, versions, hardware, seeds, hashes, tests), plus a report appended to a new `docs/PHYSICAL_CAMPAIGN_REPORT.md`. A predeclared protocol JSON (hashed before the first run) and tests (synthetic problem, schema, metadata reproducibility, threshold preservation, refusal to write into historical directories) precede the run, as in B3Y.

**Stop rules.** A failed C0 gate, a lead-tail violation, a channel-count mismatch, or a failed C1 acceptance stops the affected stage; nothing is retuned to make a control pass.

## 10. Decisions requested

1. Grid W ∈ {25, 50, 100, 200} × ζ ∈ {0.025, 0.05, 0.10, 0.20} (16 cells) with the W·ζ = 5 iso-line — or the reduced 3 × 3 version (W ∈ {50, 100, 200}, ζ ∈ {0.025, 0.05, 0.20}, no exact iso-line).
2. N_real = 200 per cell, fixed by the pilot rule of §3 (min 100, max 400).
3. Buffer 205/512 for the ensemble (instead of the B3Y 103/512).
4. Dense-energy module P: include now, include later, or omit.
5. Confirmatory contrast list of §5 and the acceptance numbers of §6 (these are new numbers, flagged for review: 0.25 × SE for the dx control, δ = 0.02 for the pilot rule, 10⁴ bootstrap resamples, Holm correction).
