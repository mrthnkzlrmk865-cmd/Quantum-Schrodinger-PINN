# Erratum / review note for `REPORT.md` (worst-point experiment, commit `edf3881`)

`REPORT.md` and the files in `results/` are kept unchanged as historical records. This note corrects them. Each statement
below was checked against `run_worst_point.py` as committed in `ca07ab6`.

## E1. The W/Z run tested ONE scalar channel, not the 63-channel block implementation
`run_worst_point.py` line 186 builds `Delta = np.full((N, 1, 1), delta_fl)` and passes 1-element lead arrays
(`lead_from_delta(np.array([delta_fl]))`, `lead_from_code_data(np.array([lam_code]), ...)`).
Consequences:
- Every W/Z number in `REPORT.md` (C3–C6, trajectories, W_max, Z_N) belongs to the 1x1 scalar path for channel n = 5.
- `wz_endpoint_blocks` was run with 1x1 matrices only. Neither its 63x63 matrix path nor any off-diagonal handling was
  exercised.
- The production solver (C1, C2) did run the full 63-channel device. Only its channel-5 entries were compared.

## E2. The closed-channel branch of W/Z was NOT exercised
`REPORT.md` §8 item 5 says the closed-root branch "is exercised in U and Q". **This is incorrect.** The W/Z run supplied
only the open channel n = 5, so `lead_from_delta`'s closed branch (delta < 0) was never called. No closed channel entered
any W/Z computation.

## E3. U and Q in `REPORT.md` §3 are 2x2 channel-5 quantities, not full open-channel S quantities
`REPORT.md` §3 says U = ‖S†S − I‖ and Q = ‖S − Sᵀ‖ of "the solver's full open-channel S". **This is incorrect** for every
row, the production rows (C1, C2) included. The driver builds the 2x2 matrix [[r_L, t_LR], [t_RL, r_R]] from the channel-5
entries only (`S2`), and:
- U = max entrywise |S†S − I|, |SS† − I| of that 2x2 matrix;
- Q = |t_LR − t_RL| of channel 5.

For the scalar W/Z path, Q = 0 is therefore structural. The sentence in §6, "Reciprocity: Q = 0 exactly, against 5.9e-15",
compares channel-5 scalars only.

## E4. Scope of the conclusion
The "INCONCLUSIVE" conclusion of `REPORT.md` §9 stands as the result of the pre-registered rule applied to the scalar
channel-5 numbers. It says nothing about the 63-channel matrix implementation or about closed channels. Both are addressed
in `amplitude_roundoff_report.md`, where full open-channel U and Q are computed.

B1 is unaffected: configuration F still FAILS the unchanged 1e-10 phase criterion.
