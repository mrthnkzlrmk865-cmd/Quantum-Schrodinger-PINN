# W/Z amplitude-roundoff investigation — single clean point (experimental; NOT a validation; NOT a production change)

**B1 is unchanged: configuration F still FAILS the unchanged 1e-10 phase criterion.** The production solver, its criteria and
every earlier result file are untouched. Nothing here validates W/Z.

## 0. Provenance

| item | value |
|---|---|
| production baseline | `74227189c524fb346247c92dc0a9602acf64730d` |
| predictions | `PREDICTIONS_amplitude_roundoff.md`, committed and pushed in `7d52af5` before any code of this investigation existed (sha256 `04badf0a…10827ee`, re-recorded in the JSON) |
| erratum for the prior report | `ERRATUM_worst_point.md` (same commit). `REPORT.md` and `results/` are unchanged |
| code | `wz_variants.py` (pluggable solve method; `solve` is checked bit-identical to `wz_rgf.wz_endpoint_blocks` on both paths), `amplitude_roundoff.py` (driver, committed before running in `3fe50fb`) |
| outputs | `results_amplitude_roundoff/amplitude_roundoff_results.json`, `results_amplitude_roundoff/amplitude_roundoff_comparisons.csv`; post-hoc `posthoc_N2_reflection.py` → `results_amplitude_roundoff/posthoc_N2_reflection.json` |
| reproduce | `python3 experiments/wz_rgf_diagnostic/amplitude_roundoff.py` (~6 min, mpmath at 110 digits) |

Run history, disclosed in full:
1. Run 1 completed every computation but crashed while writing the JSON (numpy int64 dict keys). Fixed in `d978d64`.
2. Run 2 completed. Its field `S_offdiag_max` was mis-defined: it measured the off-diagonal of the whole 10×10 S, which
   contains t itself (≈ 1).
3. Run 3 (`e274064`) replaced that field with a per-block off-diagonal measure and added the combined propagator term.
   Every field common to runs 2 and 3 is identical.

The published files are from run 3.

## 1. Setup

Configuration and inputs:
- Config F: N_y = 63, dx = 1/1024 (t = 2¹⁹), N = 3483, E = 125.69231487269683. Clean guide: V_j = 0 exactly.
- Delta_j = (E I − diag(E_perp) − V_j)·dx² is built as a full (3483, 63, 63) array. Every slice is identical and exactly
  diagonal.
- Channels: 5 open (n = 1..5) and 58 closed. Nearest closed is n = 6 (δ = −4.96e-5, |λ| = 0.9901); deepest is n = 63
  (δ = −1.86e-2, |λ| = 0.8250).
- All 63 lead roots come from `wz_rgf.lead_from_delta`.

References:
- Each W/Z run is compared with its own represented lattice R_W,n: δ_n = Delta_fl[n, n] exactly, with consistent exact
  leads. One lattice per channel, which is legitimate only because the clean device is exactly diagonal.
- Each reference uses the Chebyshev closed form at 110 digits. It agrees with the mp W/Z recursion and the mp original
  recursion to ≤ 7.7e-106 relative. The two mp recursions agree with each other to 5.5e-110. For open channels the
  reference also agrees with λ^(N−1) to ≤ 7.7e-106.
- Production is compared with R_code,n (stored a_dev,n with production σ, Γ) and with L_n (exact E − E_perp).

Paths and methods:
- "63x63" path: methods solve, inv, lu, diagdiv, refine.
- "1x1_ch5" path: the same five methods plus refine_exact.
- The methods are defined in `wz_variants.py`. No method normalizes |P|.

Metrics:
- ε_A = |t|/|t_ref| − 1 and ε_φ = arg(t·conj t_ref), for t_RL of each channel.
- U = entrywise max of |S†S − I| and |SS† − I| over the full open-channel S (10×10).
- Q = max|S − Sᵀ|.

## 2. MEASURED findings

### 2.1 Part B — full 63-channel matrix path (method `solve`)

| quantity | W/Z (63x63) | production original (same device) |
|---|---|---|
| off-diagonal entries of W_j, P_j, Ĝ blocks, Z_N, Z_N⁻¹ (all slices) | exactly 0 | — |
| channel mixing inside the S blocks | exactly 0 | exactly 0 |
| U (full open S) | 3.86e-13 (methods: 3.86e-13 … 5.24e-13) | 1.18e-13 |
| Q (full open S) | exactly 0 (all methods) | 5.9e-15 |
| max over open channels of \|ε_φ\| vs own lattice | 3.24e-15 | 6.9e-12 vs R_code; **1.54e-10 vs L** (ch5) |
| max over open channels of \|ε_A\| | 1.93e-13 | 5.7e-14 |
| channel-5 ε_A, ε_φ | +4.820e-14, −1.19e-15 | −1.08e-15, +1.54e-10 (vs L) |
| channel-5 1x1 vs 63x63 | **bit-identical** (W, P, Ĝ₁₁, Ĝ_N1; all 5 shared methods) | — |

Per open channel, method `solve`:

| n | W/Z ε_A | W/Z ε_φ | production ε_A | production ε_φ vs L | R_W − L phase |
|---|---|---|---|---|---|
| 1 | −1.43e-13 | +3.2e-15 | −3.5e-15 | −1.0e-11 | −3.9e-16 |
| 2 | −5.95e-15 | −2.0e-15 | +3.2e-15 | +2.8e-11 | +1.7e-15 |
| 3 | +1.93e-13 | −1.7e-17 | −5.4e-14 | +1.7e-11 | ~0 |
| 4 | +1.37e-13 | −1.7e-15 | +5.7e-14 | +3.3e-11 | ~0 |
| 5 | +4.82e-14 | −1.2e-15 | −1.1e-15 | +1.54e-10 | 0 (exact) |

So the "44×" amplitude ratio reported for channel 5 is channel-specific. Over all open channels the W/Z maximum |ε_A| is
3.4× the production maximum.

Closed channels and conditioning:

| quantity | measured |
|---|---|
| W/Z closed root, max rel. error of (λ − 1) vs exact root of its δ | 2.25e-16 (n = 6: 7.2e-17; n = 63: 1.9e-17); λ real, sin θ = 0 for all 58 |
| production closed root, rel. error of (λ − 1) vs exact root of its own c | max 4.67e-13 (n = 6); n = 63: 1.5e-15 |
| closed Ĝ₁₁ / Ĝ_NN max rel. error | 6.2e-15 / 1.0e-15 |
| closed Ĝ_N1 = Ĝ_1N max rel. error | 4.8e-13 (n = 6: 3.0e-13; n = 63: 1.0e-13) |
| max_j \|W_j\| (entrywise) | 0.21211 |
| cond₂(I − W_j) | 1.212111 at every slice |
| cond₂(Z_N); \|Z_N\| for ch5 | 91.97; 4.209e-3 |
| residual max\|I − (I − W_j)P_j\| (double), max over j | 2.2e-16 (solve, inv, lu); 1.1e-16 (diagdiv, refine) |
| residual \|I − Z_N Z_N⁻¹\| | 1.1e-16 |

### 2.2 Part C — stage-by-stage decomposition of ln(t_d / t_ref)

The decomposition is exact for every run and every open channel: the residual is ≤ 7.7e-106. The terms are defined in
`amplitude_roundoff.py::telescope`. Amplitude terms for channel 5 (identical on both paths):

| stage | solve / inv / lu / refine | diagdiv / refine_exact |
|---|---|---|
| P formation, Σ ln\|P_d (1 − W_d)\| | +4.757e-14 | +5.014e-14 |
| W recurrence, Σ ln\|(1 − W_d)⁻¹/λ\| (incl. W_1 initialization) | −9.2e-19 | −2.57e-15 |
| **combined propagator modulus Σ ln\|P_d/λ\|** | **+4.757e-14** | **+4.757e-14** |
| propagation multiplications (P @ G) | +3.9e-16 | +3.9e-16 |
| Z_N inputs + formation + inversion + GN1 product | −3.9e-17 | +2.7e-17 |
| Fisher–Lee (lead sin θ + product) | +2.9e-16 | +3.1e-16 |
| total ε_A | +4.820e-14 | +4.827e-14 |

Observations:
- **Every method, on both paths, returns the same double P at every slice in channel 5.** That P equals the correctly
  rounded λ_W, fl(λ_W) = 0.9999977853204164 + 0.0021046031128090737i.
- Its modulus defect is |fl(λ_W)| − 1 = +1.3661e-17. Multiplied by N − 1 this gives 4.757e-14, which is 98.7 % of ε_A.
- diagdiv and refine_exact give a different W sequence (3482 distinct values instead of 1) but the same P. The
  formation / W-recurrence split therefore shifts between the two rows while their sum does not.
- Per-step channel-5 propagator defect: amplitude +1.366e-17 (constant), phase −2.9e-20 (formation) plus +1.25e-19 (W).
- The channel-5 phase error is dominated by the propagation multiplications: −1.54e-15 out of −1.19e-15 total.

**Channels 1, 3 and 4 do discriminate between methods.** In these channels the methods produce different P sequences
(up to 517 distinct values) and clearly different ε_A:

| n | solve | inv | lu | diagdiv | refine | combined propagator term / ε_A (range over methods) |
|---|---|---|---|---|---|---|
| 1 | −1.43e-13 | −1.31e-13 | −1.43e-13 | +1.54e-13 | −2.34e-13 | 0.985 … 1.011 |
| 3 | +1.93e-13 | +2.07e-13 | +1.93e-13 | −2.62e-13 | +1.65e-13 | 0.983 … 1.023 |
| 4 | +1.37e-13 | +1.37e-13 | +1.37e-13 | −1.62e-13 | +1.37e-13 | 0.994 … 1.023 |
| 2 | −5.95e-15 (all) | | | | | 0.73 (ε_A is only 6e-15; the remaining 1.6e-15 is multiplication rounding) |

Iterative refinement (`refine`, residual in double):
- It halves the P residual (2.2e-16 → 1.1e-16 double; 2.4e-16 → 1.2e-16 exact).
- It does **not** reduce the amplitude error. Max |ε_A| goes from 1.93e-13 to 2.34e-13, and channel 1 gets worse.
- With the residual computed exactly (`refine_exact`, 1x1 only), it returns the correctly rounded P, which is what solve
  already returned. The drift is unchanged.

A smaller residual did not mean a smaller amplitude drift in any measured case.

### 2.3 Part D — N = 1 and N = 2, all 63 channels

| check | N = 1 | N = 2 |
|---|---|---|
| closed form vs mp W/Z vs mp original vs direct mp inverse of M (max rel. diff) | 1.7e-108 | 2.9e-108 |
| double W/Z vs exact, open channels (max rel., all methods) | 1.4e-16 | 3.2e-16 |
| double W/Z vs exact, closed near threshold n = 6 | 3.5e-17 | 3.5e-17 |
| double W/Z vs exact, closed far n = 63 | 1.6e-16 | 4.2e-17 |
| double W/Z vs exact, all closed (max) | 2.6e-16 | 5.7e-16 |
| U / Q of W/Z open S | ≤ 4.4e-16 / 0 | ≤ 6.7e-16 / 0 |
| production original (vs its own R_code lattice) | not constructible (`build_device` needs ≥ 1 interval) | t: rel. ≤ 4.5e-15; r: abs. ≤ 4.5e-15 (post-hoc); U 6.7e-16; Q 0 |

The production N = 2 relative error of r stored in the JSON (up to 1.5e3) is ill-posed: the exact r of R_code there is as
small as 2.9e-19. The post-hoc absolute errors in `posthoc_N2_reflection.json` are the meaningful measure.

## 3. Pre-registered hypothesis H_A (rule fixed in the predictions)

The rule is evaluated on channel 5. For all 11 runs:
- P formation explains ε_A within 10 % + 2e-15. The largest deviation is 1.9e-15 (diagdiv / refine_exact).
- The W-recurrence term is < 10 % of it: ≤ 5.1 % (diagdiv / refine_exact), 2e-5 % otherwise.

No other stage reaches 50 %. **Verdict by the pre-registered rule: SUPPORTED, but not independently discriminated in
channel 5**, because all methods produced bit-identical P there, as the rule anticipated.

Additional evidence, **not pre-registered** (the analysis of channels 1–4 was added after seeing the data, but every
number is a direct measurement):
- In channels 1, 3 and 4 the methods do produce distinct P sequences.
- ε_A changes sign and magnitude with the method.
- The combined propagator-modulus term Σ ln|P_d,j/λ| tracks ε_A to within 2.3 % in every run.

This is the method-discriminating test that channel 5 could not provide. Note that the pre-registered term (P formation
alone) does not explain channel 1 or 3 on its own. There W varies, and part of the modulus defect is attributed to "W
recurrence" relative to the intermediate double W. Only the combined term is robust to that attribution choice.

## 4. Predictions vs outcome (`PREDICTIONS_amplitude_roundoff.md`, unedited)

| prediction | outcome | verdict |
|---|---|---|
| P-formation term ≈ +4.76e-14 (solve, ch5) | +4.757e-14 | confirmed (a reproduction, not independent) |
| W recurrence contributes ≲ 1e-15 | −9.2e-19 (solve); −2.6e-15 (diagdiv, refine_exact) | confirmed for solve; **exceeded** when the W sequence differs (attribution shift; combined term unchanged) |
| propagation multiplications ≲ 6e-15 | ≤ 3.8e-15 (all channels and methods) | confirmed |
| Z_N / Fisher–Lee stages O(1e-16..1e-15) | ≤ 3.1e-16 (ch5) | confirmed |
| off-diagonals stay exactly 0 | exactly 0 | confirmed |
| ch5 bit-identical between 1x1 and 63x63 (uncertain) | bit-identical for all 5 shared methods | confirmed |
| methods give P within ~1 ulp of each other | ch5: identical. Ch1/3/4: distinct sequences; the ulp distance was **not measured** | partly confirmed |
| scalar division within 1 ulp of solve | ch5: identical P; ch1/3/4 P sequences differ (distance not measured) | partly confirmed |
| refine (double): ≤ 1 ulp change, no systematic reduction | no reduction (max \|ε_A\| increased 1.93e-13 → 2.34e-13) | confirmed |
| refine_exact ≈ correctly rounded P, drift not removed | P = correctly rounded (= solve's); drift unchanged | confirmed |
| no method brings ch5 \|ε_A\| below ~1e-15 | all 4.82e-14 … 4.83e-14 | confirmed |
| \|ε_φ\| ≤ 2e-15, all open channels and methods | ch5 ≤ 1.53e-15; **ch1 up to 3.24e-15** | **mildly violated** (channel 1) |
| U ~1e-13 (≤ 8e-13) for W/Z; production ~1e-13 | W/Z 3.9e-13 … 5.2e-13; production 1.18e-13 | within the bound, above the central estimate |
| W/Z Q exactly 0; production ~1e-14 | 0; 5.9e-15 | confirmed |
| W/Z closed root rel. error ≤ ~4 eps | 2.25e-16 (≈ 1 eps) | confirmed |
| production closed root up to ~eps/\|δ\| (~2e-12 at n = 6, ~1e-14 at n = 63) | 4.7e-13; 1.5e-15 | confirmed (within bound) |
| closed Ĝ₁₁/Ĝ_NN ≤ 1e-14; Ĝ_N1 ≤ 8e-13 | 6.2e-15 / 1.0e-15; 4.8e-13 | confirmed |
| max\|W\| ≈ 0.2121; cond(I − W) ≈ 1.212; cond(Z_N) ≈ 92 | 0.21211; 1.21211; 91.97 | confirmed |
| Part D: exact formulations agree ≤ 1e-90; double ≤ 1e-14; U ≤ 1e-14; Q = 0; production N = 2 ≤ 1e-13 | 2.9e-108; ≤ 5.7e-16; ≤ 6.7e-16; 0; ≤ 4.5e-15 | confirmed |

## 5. POST-HOC explanations (consistent with the data, not independently tested)

1. **The channel-5 drift is a representation floor of λ in double.** Every method returned P = fl(λ_W), so the defect is
   not a solve-method error. For P = cos θ + i sin θ with θ ≈ 2.1e-3, Re P lies just below 1, where ulp(1) = 1.1e-16.
   - Moving Re P by one ulp changes |P| by ~1.1e-16.
   - Moving Im P by one ulp (4.3e-19) changes |P| by only ~9e-22.
   - So no nearby double has both the correct phase and |P| = 1 to better than ~ulp(1)/2.
   - The nearest double has |P| − 1 = 1.37e-17. Applied N − 1 times to an unchanging P, this gives a coherent drift.
   - This argument was not checked by enumerating neighbouring doubles.
2. **Channels with varying W/P (1, 3) still drift coherently.** Their P sequences take a few hundred distinct values, but
   the per-step modulus defects evidently do not average out. The combined term tracks ε_A to ≤ 2.3 %. Why the sign and
   size depend on the method (e.g. ch1: −2.3e-13 … +1.5e-13) was not analysed further.
3. **Residual size is not the relevant metric.** Refinement lowers ‖I − AP‖ but not |P| − |λ|, because the target |λ| = 1
   lies between representable values. Normalizing |P| was not applied (instruction) and would be unjustified for coupled
   blocks.

## 6. HYPOTHESES (untested here)

- H1: amplitude drift grows ∝ N with a per-step bound of ~ulp(1)/2 ≈ 1.1e-16. The worst case is ~3.9e-13 at N = 3483.
  Observed max: 2.6e-13. N-scaling was not run.
- H2: in coupled (disordered) blocks W_j and P_j vary with j and are not diagonal. Whether the modulus defects then stay
  coherent is unknown.
- H3: the original formulation's amplitude errors (≤ 5.7e-14 here) come from a different mechanism. Its per-slice state
  was not instrumented (the production solver must not be modified).

## 7. UNRESOLVED issues

1. **Matrix multiplication order for coupled blocks is NOT validated.** The clean Delta_j is exactly diagonal, so every
   product here was elementwise. The left/right chains (Ĝ_j1 = P Ĝ, Ĝ_1j = Ĝ P) and the old-value use in the Ĝ₁₁ update
   were not exercised with non-commuting matrices. No disordered or synthetic coupled system was run, by instruction.
2. The ulp distance between the P sequences of different methods in channels 1/3/4 was not measured.
3. N-scaling of the drift (H1) and its behaviour under disorder (H2) are untested.
4. The production N = 1 case cannot be built with `build_device`. Only W/Z and the mp formulations were tested at N = 1.
5. Whether a reformulation that propagates a small quantity for P (e.g. P − 1) would reduce the drift is an open idea. It
   has not been tested and is not a recommendation.

## 8. Conclusion

At this single clean point, the W/Z amplitude error is explained, to within 2.3 % in every run for channels 1, 3, 4 and 5 (channel 2 has ε_A of only 6e-15), by the modulus defect
of the per-slice double propagator P accumulated over the device.
- In channel 5 that defect is the representation error of the correctly rounded λ (|fl(λ)| − 1 = 1.37e-17). Every solve
  method (numpy solve, inv, LU, scalar division, refinement with double or exact residual) reproduces it bit-for-bit, so
  no method tested removes it.
- In channels 1, 3 and 4 different methods produce different drifts of either sign (|ε_A| up to 2.6e-13), and the same
  propagator term tracks each of them.
- Iterative refinement lowers residuals but not the drift.
- The full 63-channel matrix path is bit-identical to the scalar path for channel 5. It handles all 58 closed channels to
  ≤ 4.8e-13 relative (Ĝ_N1) and ≤ 6.2e-15 (Ĝ₁₁). It is exactly reciprocal and has U ≤ 5.2e-13.
- The phase error stays ≤ 3.2e-15 in every open channel, against 1.0e-11 to 1.54e-10 for the original.

None of this validates W/Z. It is one energy and one device. The device is clean, so its blocks are diagonal, and the
coupled-block multiplication order is untested. B1 remains FAILED at configuration F under the unchanged 1e-10 criterion.
Further work (N-scaling, other points, the 424-row sweep, coupled blocks) needs explicit approval.
