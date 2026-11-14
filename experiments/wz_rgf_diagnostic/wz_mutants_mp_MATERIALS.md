# Clean-room transcription materials for `wz_mutants_mp.py` (exact text supplied to the implementer)

This file is the COMPLETE set of materials given to the clean-room implementation context (§4.3 of
`PREDICTIONS_coupling_scaling_v2.md`). It contains only: (1) the correct-code listing of §4.2, (2) the §4.2 change column
(mutant label -> code change; the designated-block and oracle columns are deliberately withheld), (3) the N = 1 rules, and
(4) the required interface and precision rules. The implementer is instructed not to read `wz_mutants.py`, `wz_rgf.py`,
`coupled_blocks.py`, any file under `cs2_t6_reference/` or any other repository source.

## 1. Variables (W/Z formulation, §2.1 / §4.2; dimensionless Green units)

- Inputs: Delta = [Delta_1, ..., Delta_N], each an N_y x N_y complex matrix (W/Z variable Delta_j); lam_minus_1 = the N_y
  diagonal entries of Lambda - I (lead roots minus one). Lm := diag(lam_minus_1) = Lambda - I. I = N_y x N_y identity.
- `solve(X, B)` means X^-1 B. `@` is the matrix product. `.T` is the plain (non-conjugate) transpose.

## 2. Correct-code listing (§4.2, verbatim)

Correct code (`wz_mutants.py`, `mutant=None`, identical to `wz_rgf.wz_endpoint_blocks`):
- init: `W = 2Δ_1 + Lm`; `G11 = solve(W − I, I)`; `Gj1 = G1j = G11`.
- loop: `W = 2Δ_j + solve(I − W, W)`; `P = solve(I − W, I)`; `G11 = G11 − G1j@P@Gj1`; `Gj1 = P@Gj1`; `G1j = G1j@P`.
- closing: `Z = 2Δ_N + Lm + solve(I − W, W)`; `Zi = solve(Z, I)`; `G11 = G11 + G1j@Zi@Gj1`; `GN1 = −Zi@Gj1`;
  `G1N = −G1j@Zi`; `GNN = Zi`.

N = 1 (no loop; one slice receives both leads): `Z = 2 Delta_1 + 2 Lm`; `Zi = solve(Z, I)`; all four endpoint blocks equal
`Zi` (G11 = GN1 = G1N = GNN = Zi).

For N >= 2 the loop runs over slices s = 2..N-1 (zero iterations when N = 2), using Delta_s.

## 3. Change column (§4.2; mutant label -> exact code change; nothing else)

| mutant | exact code change (§4.2 change column) |
|---|---|
| M1 | `Gj1 = Gj1@P` |
| M2 | `G1j = P@G1j` |
| M3 | M1 and M2 |
| M4 | chains updated first, then `G11 = G11 − G1j@P@Gj1` |
| M4L | `Gj1 = P@Gj1` before the G11 update; `G1j` after |
| M4R | `G1j = G1j@P` before the G11 update; `Gj1` after |
| M5 | `G11 = G11 − Gj1@P@G1j` |
| M6 | `GN1 = −(Gj1@Zi)` |
| M7 | `G1N = −(Zi@G1j)` |
| M8 | `G11 = G11 + Gj1@Zi@G1j` (closing) |
| M9 | `P = solve(I − W_old, I)` before the W update |
| M10 | `Z = 2Δ_N + solve(I − W, W)` (Lm dropped) |
| M12 | `solve((I−W).T, W.T).T` for (I − W)⁻¹W (loop and Z) |
| M13 | `G11 = G11 + G1j@P@Gj1` (loop) |
| M14 | init `G11 = solve(I − W, I)` |
| M15 | M6 + M7 |
| M16 | closing: GN1, G1N first; `G11 = G11 + G1N@Zi@GN1` |
| M11 | N = 1 boundary `2Δ_1 + Lm` |

Notes taken verbatim from §4.2 change-column semantics: each mutant changes exactly the stated expression(s) of the
correct code; everything else is identical. "M3: M1 and M2" means both changes. "M15: M6 + M7" means both changes.
"M12 ... (loop and Z)" means the replacement applies wherever (I - W)^-1 W is formed, in the loop and in the closing Z.

## 4. Required interface and precision rules

- Module `wz_mutants_mp.py`, function `mp_endpoint_blocks(Delta, lam_minus_1, mutant=None, dps=60)`.
  - `Delta`: list of N mpmath matrices (or nested lists convertible by `mpmath.matrix`); `lam_minus_1`: list of N_y values
    convertible by `mpmath.mpc`. Inputs are converted inside a private `mpmath.MPContext()` with `ctx.dps = dps`; the
    function must not modify the global `mpmath.mp` context.
  - `mutant`: None (correct code) or one of the labels of Section 3. Any other value raises ValueError.
  - Returns a dict with keys `Ghat11`, `GhatN1`, `Ghat1N`, `GhatNN` (mpmath matrices), `trace` (list of (name, matrix) in
    computation order, names: N = 1: "Z", "Zi", then "G11", "GN1", "G1N", "GNN"; N >= 2: "W_1", "G11^(1)", then for every
    loop slice s = 2..N-1: f"W_{s}", f"P_{s}", f"G11^({s})", f"G_{s}1^({s})", f"G_1{s}^({s})", then "Z", "Zi",
    "G11", "GN1", "G1N", "GNN"), `Z` and `Zi`.
  - All solves use `ctx.lu_solve` or `ctx.inverse` at the context precision; no numpy, no floating Python arithmetic on
    matrix entries.
- A module constant `MUTANTS` (tuple of the Section-3 labels plus the order in Section 3) and `PROVENANCE` (a string stating
  that the module was written only from this materials file).
- Do not import any repository module. Allowed imports: mpmath, the Python standard library.
