# Predictions: W/Z RGF worst-point diagnostic (written BEFORE the experiment was executed)

Reference commit: `7422718`. Point: B1 configuration F (N_y = 63, dx = 1/1024, L_d = 3, N = 3483),
E = 125.69231487269683 (the spec-grid double), channel n = 5, E_perp = (5 pi)^2 / 2 (double), sin(k dx) = 2.1046e-3.
Committed facts used: signed phase error of the production RGF against the exact lattice = +1.5424928001550615e-10 (L->R),
+1.542470040583065e-10 (R->L) (`results/rgf_b1_validation/phase_audit_hp.csv`). B1 fails at F (criterion 1e-10, unchanged).

Only arithmetic on the stored inputs (no solve of the model) informs the numbers below; the device-representation value P1 was
computed in the preceding read-only audit. theta = 2 asin(sqrt(delta/2)) ~= 2.1046e-3, delta ~= 2.2147e-6, eps = 2.22e-16.

## Represented lattices (expected properties)
- L: exact lattice with the double E and E_perp taken as exact.
- R_dev: device diagonal a_dev = fl(E - fl(2t + E_perp)), consistent exact leads. Expected c_dev - c = O(1e-16).
- R_code: a_dev with the double lead data of smatrix_spec.lead_modes (via stationary's output). Expected |c_lead - c_dev| <= ~1.1e-16
  (mismatch) and Gamma not exactly 2t Im(lambda) (both rounded separately).
- R_W: Delta_fl = fl(E - E_perp) * dx^2. E_perp/2 <= E <= 2 E_perp holds here, so E - E_perp is exact (Sterbenz) and dx^2 is a
  power of two: **R_W is predicted to be IDENTICAL to L** (S(R_W) - S(L) = 0 exactly). This is specific to this clean, near-threshold
  channel and does not generalize to disordered devices.
- R_W_code_lead: Delta_fl with the arccos-based double lead data: lead/device mismatch of order 1e-16 in c.

## Original formulation: S_double - S(L) = [R_dev - L] + [R_code - R_dev] + [S_double - R_code]
| term | predicted signed phase contribution | basis |
|---|---|---|
| device representation R_dev - L | +1.61e-10 (P1 = (N-1)(arccos c_dev - arccos c)) | prior audit arithmetic |
| lead mismatch R_code - R_dev | magnitude <~ 1e-11 (heuristic eps/(4 theta^2) ~ 1.3e-11), sign unknown; boundary effect, not growing with N | heuristic, untested |
| recurrence rounding S_double - R_code | magnitude <~ 5e-11 | residual of the P1 fit (~30 %) |
| total | +1.5425e-10 (committed) | reproduction gate: the production run must reproduce the committed value to <= 1e-15 |
Amplitude errors: |eps_A| <~ 1e-14 for the double solver; reflection |r| of R_code (mismatch-induced) <~ eps/(2 theta^2) ~ 2.5e-11;
|r| = 0 exactly for L, R_dev, R_W. Unitarity defect of the double solvers <~ 1e-13.

## W/Z formulation
| comparison | prediction | status |
|---|---|---|
| R_W - L | exactly 0 | follows from Sterbenz + power-of-two dx^2 (checked in the run) |
| W/Z (delta-based leads) vs L | phase error of order N eps ~ 1e-13 to 1e-12 (hypothesis: no 1/sin theta amplification) | HYPOTHESIS |
| W/Z (code leads) vs L | dominated by the lead mismatch, order 1e-11 (heuristic eps/(4 theta^2)) | HYPOTHESIS |
| W/Z (code leads) vs R_W_code_lead | recurrence rounding only, order 1e-13 to 1e-12 | HYPOTHESIS |
Expected phase-error scaling (not testable in this single-N experiment): original ~ N eps / sin(theta); W/Z ~ N eps; lead mismatch
not growing with N. Fixed-point trajectory: with consistent leads w_1 equals the device fixed point up to rounding and the deviation
stays at rounding level; with mismatched leads the initial offset (~1e-16 relative in c, ~eps/theta in w) persists and rotates
(marginal map |lambda^2| = 1) without decaying.

## Decision rule for the conclusion (fixed in advance)
- "W/Z improves accuracy" at this point only if |phase error| vs L of W/Z with delta-based leads is at least 10x smaller than the
  original's (1.54e-10), the telescoping identities hold (residual < 1e-30), and amplitude/unitarity/reciprocity are not worse than
  the original's by more than 10x.
- "W/Z does not improve accuracy" if the W/Z phase error is >= the original's.
- Otherwise "inconclusive". In every case: ONE point only; no claim of validation; B1 still FAILS at F under the unchanged 1e-10.
