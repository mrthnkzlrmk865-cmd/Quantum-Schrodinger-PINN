# W/Z RGF diagnostic — worst B1 point only (experimental; NOT a production change, NOT a validation)

Reference commit `7422718`. Predictions committed before any run: `PREDICTIONS.md` (commit `3ce732d`,
sha256 `95d09070…a09522157`, recorded in the result JSON). Driver: `run_worst_point.py`; W/Z solver: `wz_rgf.py`
(numpy only, standalone); post-hoc check: `posthoc_amplitude_check.py` (written after the run, labelled as such).
Outputs: `results/worst_point.json`, `results/worst_point_comparisons.csv`, `results/worst_point_trajectory.csv`,
`results/posthoc_amplitude.json`. Reproduce: `python3 experiments/wz_rgf_diagnostic/run_worst_point.py` (~5 s).
The run was repeated after committing the driver (HEAD `ca07ab6`); all numbers were identical to the first run.
`git_dirty = true` in the JSON is caused only by the untracked `results/` directory the run itself writes.

**B1 status is unchanged: configuration F still FAILS the unchanged 1e-10 phase criterion** (`results/rgf_b1_validation/
summary.json`: `all_b1_passed = false`; `phase_audit_hp_summary.json` F: `b1_phase_criterion_pass = false`). Nothing in
this experiment modifies the production solver, B1 criteria, or any committed result.

## 1. Point

Config F: N_y = 63, dx = 1/1024 (t = 2^19), L_d = 3, N = 3483 slices; E = 125.69231487269683 (spec-grid double);
channel n = 5, E_perp = (5π)²/2 (double) = 123.37005501361698. δ = (E − E_perp)dx² = 2.21467958362565e-6,
θ = 2.1046046664843e-3, sin θ = 2.1046031128091e-3. The production run uses `stationary.scattering_matrix` on
`validate_rgf_b1.device(63, 1/1024, 3.0)` (imported read-only). The production Hamiltonian has zero off-diagonal mode coupling
(max 0.0), so channel 5 is an exactly decoupled scalar chain. Its exact references are given by a Chebyshev closed form at
110 digits, cross-checked by an mpmath recursion: the two agree to ≤ 8e-106. For L, R_dev and R_W the reference also agrees
with λ^(N−1) to ≤ 8e-106.

## 2. Represented lattices

| lattice | device diagonal | leads | c_device − c_lead | consistent |
|---|---|---|---|---|
| L | exact −2t c, c = 1 − δ (double E, E_perp exact) | exact root of c | 0 | yes |
| R_dev | stored a_dev = fl(E − fl(2t + E_perp)) | exact root of c_dev | 0 | yes |
| R_code | a_dev | production double σ, Γ (`lead_modes` via `scattering_matrix`) | 1.11e-16 | **no** |
| R_W | t(2Δ_fl − 2), Δ_fl = fl(E − E_perp)·dx² | exact root of c_W | 0 | yes |
| R_W_code_lead | t(2Δ_fl − 2) | production double σ, Γ | 1.37e-17 | **no** |

Further facts:
- c_dev − c = −9.7e-17 (one ulp of the diagonal).
- E − E_perp is exact here: Sterbenz applies, since E_perp/2 ≤ E ≤ 2E_perp, and dx² is a power of two. So
  **R_W is bit-for-bit L** (`R_W_identical_to_L = true`), as predicted.
- This property comes from the clean, near-threshold channel. It does not carry over to disordered devices.
- For every lattice, Γ equals 2t·Im λ exactly to working precision. The R_code inconsistency is purely the c mismatch.
  Full per-lattice data (diagonal, c, λ, σ, Γ) is in the JSON under `lattices`.

## 3. Comparisons (solver vs reference of the SAME represented lattice, or vs L for the total)

ε_φ = arg(t·conj(t_ref)) for t_RL; ε_A = |t|/|t_ref| − 1; |r| is channel-5 reflection; U = ‖S†S − I‖ and
Q = ‖S − Sᵀ‖ (solver's full open-channel S).

| id | solver | reference | ε_φ | ε_A | \|q\| | \|r\| solver | U | Q |
|---|---|---|---|---|---|---|---|---|
| C1 | original (production) | L | **+1.5425e-10** | −1.08e-15 | 1.54e-10 | 2.40e-11 | 1.14e-13 | 5.9e-15 |
| C2 | original (production) | R_code | −6.88e-12 | −1.08e-15 | 6.88e-12 | 2.40e-11 | 1.14e-13 | 5.9e-15 |
| C3 | W/Z, δ-based leads | L | **−1.19e-15** | **+4.82e-14** | 4.82e-14 | 1.87e-14 | 9.6e-14 | 0 |
| C4 | W/Z, δ-based leads | R_W (≡ L) | −1.19e-15 | +4.82e-14 | 4.82e-14 | 1.87e-14 | 9.6e-14 | 0 |
| C5 | W/Z, production leads | L | +3.18e-15 | +1.45e-13 | 1.45e-13 | 2.72e-12 | 2.9e-13 | 0 |
| C6 | W/Z, production leads | R_W_code_lead | −3.31e-15 | +1.45e-13 | 1.45e-13 | 2.72e-12 | 2.9e-13 | 0 |
| X1 | scalar replica of original | R_code | +3.13e-11 | +1.10e-13 | 3.13e-11 | 2.42e-11 | — | — |

**Reproduction gate (passed):** the production ε_φ reproduces the committed `phase_audit_hp.csv` values. The differences
are 4.5e-18 (L→R) and 4.3e-19 (R→L), against a gate of ≤ 1e-15.

The R→L direction (t_LR) is in the CSV. It differs from L→R by ≤ 2.3e-15 for the original and is identical for W/Z.

## 4. Exact telescoping decompositions (signed phase of t_RL; residual computed at 110 digits)

| formulation | device representation | lead mismatch | recurrence rounding | total | residual |
|---|---|---|---|---|---|
| original: [R_dev−L] + [R_code−R_dev] + [S−R_code] | **+1.6108e-10** | +5.3e-14 | −6.88e-12 | +1.5425e-10 | 0 |
| W/Z production leads: [R_W−L] + [R_W_code_lead−R_W] + [S−R_W_code_lead] | 0 (exact) | +6.5e-15 | −3.3e-15 | +3.2e-15 | 0 |
| W/Z δ-leads: [R_W−L] + [S−R_W] | 0 (exact) | — | −1.19e-15 (rounding incl. double lead data) | −1.19e-15 | 0 |

The same decomposition for r_L (magnitudes) is in the JSON. In the original, the reflection comes from the lead
mismatch (|r(R_code)| = 2.17e-11) and the recurrence (1.0e-11). R_dev itself is reflectionless.

**Interpretation for the original solver at this point.** About 104 % of the phase error is the stored-diagonal
representation (R_dev − L). Recurrence rounding offsets about 4 % with the opposite sign. The lead mismatch is
negligible in phase (5e-14), but it is the dominant source of |r|. This confirms the prior audit's root-cause attribution
at this point, now through an exact decomposition rather than a correlation.

## 5. Fixed-point trajectories (w_j − w*, every 50 slices in `worst_point_trajectory.csv`)

- **Original:** this is the scalar replica's trajectory, with D_j mapped exactly to W units, w_j = 1 + D_j/t, and compared with w* = 1 − 1/λ_dev. The deviation is
  5.3e-14 at j = 1 and stays at 5.3–7.1e-14 to j = N. It persists without decay (relative 3.4e-11), consistent with the
  marginal map |λ²| = 1.
- **W/Z, δ-leads:** this and the next item are compared with w* = 1 − 1/λ_W. Here w_j − w* = 1.25e-19 at every slice. W is reproduced to sub-ulp at each step and does not drift.
- **W/Z, production leads:** the offset is 6.5e-15 at j = 1 and 6.6e-15 at j = N. The mismatch offset persists, as
  predicted, without growth or decay.
- max |W_j| = 2.1046e-3 (= |w*|, no growth), and |Z_N| = 4.209e-3 for both lead variants.

## 6. Original vs W/Z — what changed

- **Phase:** at this point W/Z removes the device-representation error exactly, because R_W ≡ L. Its remaining phase
  error is about 1e-15, roughly 1.3 × 10⁵ times smaller than the original's 1.54e-10.
- **Reflection:** |r| = 1.9e-14 with δ-leads, against 2.4e-11 for the original.
- **Reciprocity:** Q = 0 exactly, against 5.9e-15. U is comparable (9.6e-14 vs 1.14e-13).
- **Amplitude: W/Z is worse.** ε_A = +4.82e-14 (δ-leads) and +1.45e-13 (production leads), against −1.08e-15 for the
  original. That is 44× and 134× larger.
  - Post-hoc mechanism (`posthoc_amplitude_check.py`; NOT pre-registered). In the clean channel W_j is the same double
    at every slice, so the per-slice propagator P = fl((I − W)⁻¹) is also one fixed double.
  - |1/(1 − W_double)| − 1 = −2.6e-22 exactly, but the rounded P has |P| − 1 = +1.37e-17 (rounding of Re P ≈ 1).
  - Ĝ_N1 ∝ P^(N−1) turns this into a coherent drift: (N − 1)(|P| − 1) = 4.76e-14. That is 98.7 % of the observed
    4.82e-14.
  - The error has moved from phase, where the original suffers 1/sin θ amplification, to amplitude, where it is
    O(N·ulp(1)) with no 1/sin θ factor. Its worst case for this N is about N·ε/2 ≈ 3.9e-13.
  - This mechanism relies on W being exactly constant, which is specific to a clean channel.

## 7. Predictions vs outcome (PREDICTIONS.md, unedited)

| prediction | outcome | verdict |
|---|---|---|
| R_dev − L ≈ +1.61e-10 | +1.6108e-10 | confirmed |
| lead mismatch phase ≲ 1e-11, sign unknown | +5.3e-14 | within bound, about 200× smaller than the heuristic |
| \|r\| of R_code ≲ 2.5e-11 | 2.17e-11 | confirmed |
| recurrence rounding ≲ 5e-11 | −6.88e-12 | within bound, smaller |
| total reproduces committed value to ≤ 1e-15 | 4.5e-18 | confirmed |
| R_W ≡ L exactly | exact | confirmed |
| W/Z δ-leads vs L phase about 1e-13 to 1e-12 (hypothesis) | 1.2e-15 | smaller than predicted |
| W/Z production leads vs L about 1e-11 (lead mismatch dominates) | phase 3.2e-15, \|q\| 1.45e-13 | **prediction wrong**: the mismatch appears in \|r\| (2.7e-12), not in phase |
| W/Z production leads vs R_W_code_lead about 1e-13 to 1e-12 | \|q\| 1.45e-13 | confirmed (lower end) |
| \|ε_A\| ≲ 1e-14 for double solvers | original 1.1e-15; W/Z 4.8e-14 and 1.45e-13 | **violated by W/Z** |
| unitarity defect ≲ 1e-13 | 1.14e-13; 9.6e-14; 2.9e-13 | marginal or violated (W/Z production leads) |
| trajectories: consistent leads stay at rounding level; mismatch offset persists without decay | as stated | confirmed |

## 8. Unexpected results and limitations

1. The amplitude error of W/Z exceeds the original's by more than 10× (section 6). This was not anticipated in the
   predictions.
2. The scalar replica (X1) does not reproduce the production solver bit-for-bit: max |difference| 3.8e-11, and its
   ε_φ = +3.1e-11 against −6.9e-12 for production. Production's per-slice internal state cannot be extracted without
   modifying it. The "original" trajectory in section 5 is therefore the replica's and is only indicative. The production
   numbers in C1 and C2 come from the real solver.
3. **One point, one N, one clean channel.**
   - The N-scaling hypotheses (original ∝ Nε/sin θ; W/Z ∝ Nε) are untested.
   - The exactness of R_W relies on Sterbenz plus a power-of-two dx² and a clean channel.
   - Constant W_j relies on a clean device. Disorder, mode mixing, closed channels, and configurations A–E were not run,
     by instruction.
4. The production-lead W/Z variant inherits the 1-ulp lead/device mismatch, which sets its |r| (2.7e-12).
5. The δ-based lead formula was exercised only for open channels at this point. The 62 other channels in the S used for
   U and Q include closed channels, so the closed-root branch is exercised in U and Q but not checked against an exact
   reference.

## 9. Conclusion (pre-registered decision rule applied as written)

The rule requires four things for "improves":
- the W/Z (δ-leads) phase error is ≥ 10× smaller than 1.54e-10;
- the telescoping residual is < 1e-30;
- amplitude, unitarity and reciprocity are not worse than the original's by more than 10×.

Outcome:
- Phase: met. 1.19e-15 is about 1.3 × 10⁵ times smaller.
- Telescoping residual: met (0).
- Unitarity: not worse (0.85×). Reciprocity: not worse (0 vs 5.9e-15).
- **Amplitude: NOT met.** |ε_A| = 4.82e-14 against 1.08e-15, so 44× worse.

The phase error is also not ≥ the original's.

**Conclusion: INCONCLUSIVE** (per the fixed rule).

At this single point, the W/Z formulation eliminates the device-representation phase error that makes B1 fail at F. In
exchange it introduces a coherent amplitude drift of about 5e-14, which is 44× the original's amplitude error, although
still far below any B1 tolerance.

This experiment does not validate W/Z. It does not change B1, and B1 still FAILS at configuration F under the unchanged
1e-10 criterion. Any extension (other energies or configurations, the 424-row sweep, N-scaling, disorder) requires
explicit approval.
