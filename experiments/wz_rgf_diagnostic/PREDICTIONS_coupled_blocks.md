# Predictions and design: non-commuting coupled-block validation of the W/Z recursion
(written and committed BEFORE any implementation or execution of this experiment)

Status:
- B1 configuration F still FAILS the unchanged 1e-10 criterion.
- W/Z is experimental and NOT validated.
- No production change is authorized.
- Earlier W/Z tests (`REPORT.md`, `amplitude_roundoff_report.md`) used exactly diagonal clean blocks, so every product was
  elementwise and **no non-commuting multiplication has been tested yet**.
- This document defines a small, deterministic experiment that exercises the block ordering. Nothing in it has been
  implemented or run.

Inspected for this design (read-only):
- `wz_rgf.py`, `wz_variants.py`, `amplitude_roundoff.py`, `amplitude_roundoff_report.md`, `ERRATUM_worst_point.md`;
- `stationary.py` (`Device`, `transverse_basis`, `device_hamiltonian_blocks`, `rgf_endpoint_blocks`, `scattering_matrix`);
- `smatrix_spec.py` (`hopping`, `lead_modes`, `fisher_lee`, `assemble_S`).

**Documentation gap found:** the `wz_rgf.py` docstring cites "REPORT.md section 1" for the audited equations, but
`REPORT.md` §1 describes the test point and contains no derivation. No committed W/Z derivation document exists. The
relations needed here are therefore derived in §1 below. Existing files are not edited.

## 1. Mathematical basis (to be tested, not assumed)

Device matrix, as in `stationary.rgf_endpoint_blocks`:
- t = 1/(2dx²).
- H_jj = 2t I + diag(E_perp) + V_j, and H_{j,j±1} = −t I.
- M = E − H_D − Σ_L ⊕ Σ_R, with Σ = diag(−t λ_n) on slice 1 and on slice N. The leads are clean (V = 0), so
  λ_n + 1/λ_n = 2c_n, c_n = 1 − δ_n, δ_n = (E − E_perp,n)dx², taking the retarded root.

W/Z variables:
- Δ_j = (E I − diag(E_perp) − V_j)dx².
- Since 2t dx² = 1, M/t has diagonal blocks 2Δ_j − 2I (+ Λ on slices 1 and N) and off-diagonal blocks I.
- Ĝ = tG.

Left-connected prefix systems. For j < N, let M^(j) be the leading j×j block submatrix of M/t, i.e. without the right lead.
The Schur recursion of the original RGF gives D_j/t = W_j − I, with:
- W_1 = 2Δ_1 + (Λ − I);
- W_j = 2Δ_j + (I − W_{j−1})⁻¹ W_{j−1};
- P_j = (I − W_j)⁻¹ = −[(M^(j))⁻¹]_jj (so W_j = I + ([(M^(j))⁻¹]_jj)⁻¹);
- Ĝ^(j)_j1 = P_j Ĝ^(j−1)_{j−1,1} (**left** multiplication);
- Ĝ^(j)_1j = Ĝ^(j−1)_{1,j−1} P_j (**right** multiplication);
- Ĝ^(j)_11 = Ĝ^(j−1)_11 − Ĝ^(j−1)_{1,j−1} P_j Ĝ^(j−1)_{j−1,1} (**old** values);
- Ĝ^(1)_11 = (W_1 − I)⁻¹.

Closing step:
- Z_N = 2Δ_N + (Λ − I) + (I − W_{N−1})⁻¹W_{N−1} = D_N/t.
- Ĝ_NN = Z_N⁻¹.
- Ĝ_N1 = −Z_N⁻¹ Ĝ^(N−1)_{N−1,1}.
- Ĝ_1N = −Ĝ^(N−1)_{1,N−1} Z_N⁻¹.
- Ĝ_11 = Ĝ^(N−1)_11 + Ĝ^(N−1)_{1,N−1} Z_N⁻¹ Ĝ^(N−1)_{N−1,1}.

Edge cases:
- N = 1: Z_1 = 2Δ_1 + 2(Λ − I), and all four blocks equal Z_1⁻¹.
- N = 2: the loop is empty. Z_2 = 2Δ_2 + Λ − 2I + (I − W_1)⁻¹, which is the Schur complement of [[W_1 − I, I], [I, 2Δ_2 − 2I + Λ]].

Fisher–Lee (open channels, s_n = sin θ_n = Γ_n/(2t)):
- t_RL = 2i √s Ĝ_N1 √s and t_LR = 2i √s Ĝ_1N √s;
- r_L = 2i √s Ĝ_11 √s − I and r_R = 2i √s Ĝ_NN √s − I.
These are elementwise in the channel indices, so they are valid for non-diagonal blocks.

Two consequences that drive the design:
- (i) If V_j is real symmetric (time-reversal symmetric, TRS), M is complex symmetric. Then Ĝ = Ĝᵀ, Ĝ_1N = Ĝ_N1ᵀ, and
  W_j, P_j and Ĝ_11 are symmetric.
- (ii) Any mutant that reverses a product at step j changes the result by a term that contains a commutator: [P_j, Ĝ_{j−1,1}]
  in the chains, [Z_N⁻¹, Ĝ_{N−1,1}] at the closing step. These commutators vanish exactly for diagonal blocks. **A
  diagonal device therefore cannot detect ordering bugs.** Non-commuting adjacent W_j are required.

## 2. Test devices (deterministic, small, Hermitian)

Lattice and leads:
- dx = 1/8 (t = 32, exact) and L_y = 1.
- E_perp,n = (nπ)²/2 from `stationary.transverse_basis(1.0, N_y)`: 4.934802200544679, 19.739208802178716,
  44.41321980490211.
- The leads are clean. Device slices may carry V_j ≠ 0, including the first and last slices.

Energies, with δ_n = (E − E_perp,n)dx² and every |c_n| ≥ 1.087 or ≤ 0.913, far from the thresholds:

| id | N_y | E | δ_n | open / closed |
|---|---|---|---|---|
| e2a | 2 | 12 | +0.110394, −0.120925 | 1 / 1 |
| e2b | 2 | 30 | +0.391644, +0.160325 | 2 / 0 |
| e3a | 3 | 12 | +0.110394, −0.120925, −0.506457 | 1 / 2 |
| e3b | 3 | 30 | +0.391644, +0.160325, −0.225207 | 2 / 1 |
| e3c | 3 | 50 | +0.704144, +0.472825, +0.087293 | 3 / 0 |

Potential matrices, for slices j = 1..4 (a device with N slices uses V_1..V_N). All coefficients are dyadic, so V_j is
exactly representable. The scale is v0 = 4 (|V_j dx²| ≲ 0.2).

N_y = 2, with A = diag(1, −1), B = [[0, 1], [1, 0]] and C = [[0, −i], [i, 0]]:

    V_j = v0 (μ_j I + α_j A + β_j B + γ_j C)
    μ = (0.5, −0.25, 1, 0)   α = (1, 2, −1.5, 0.5)   β = (3, −1, 2.5, 1.75)   γ = (1, 0.5, −2, 1.25)

N_y = 3:

    V_j = v0 (μ_j I + α_j diag(1,0,−1) + α'_j diag(1,−2,1) + β12_j S12 + β23_j S23 + β13_j S13 + γ12_j K12 + γ23_j K23)

- S_ab = E_ab + E_ba (real symmetric).
- K_ab = i(E_ba − E_ab) (imaginary Hermitian).
- Coefficients:
  - μ = (0.5, −0.25, 1, 0); α = (1, 2, −1.5, 0.5); α' = (0.25, −0.5, 0.75, 1);
  - β12 = (3, −1, 2.5, 1.75); β23 = (−2, 1.5, 0.5, −1); β13 = (1, 0, −1.25, 2);
  - γ12 = (1, 0.5, −2, 1.25); γ23 = (−0.5, 1, 0.75, −1.5).

Potential classes:
- **R** (real symmetric, TRS): γ = 0. This is the physically relevant class for this project.
- **H** (complex Hermitian, TRS broken): all coefficients. This is an algebraic stress test with no symmetry to hide
  transposition errors. It is not physical for this project, and S ≠ Sᵀ is expected even exactly.
- **D** (diagonal control): only μ, α, α'. Every block is diagonal and every product commutes. Ordering mutants must be
  invisible here.

Case matrix:
- Main: 5 energy sets × 3 classes × N ∈ {1, 2, 3, 4} = **60 cases**. N_y = 2 cases use e2a/e2b; N_y = 3 cases use
  e3a/e3b/e3c.
- Coupling-scaling sub-matrix (discriminates ordering error from rounding): class R, e2b and e3b, N ∈ {3, 4}, with
  v0 scaled by s ∈ {1, 2⁻⁴, 2⁻⁸, 2⁻¹²} (dyadic, exact). That is 4 base configurations × 4 scales = **16 cases**. The s = 1
  members duplicate main-matrix cases.

These are scaling tests of the coupling strength in N ≤ 4 devices. No N-scaling and no large devices are run.

## 3. Solvers and independent references

Solvers under test (double precision, unmodified code, read-only):
- **A. W/Z:**
  - Inputs: Δ_j = (E I − diag(E_perp) − V_j)·dx² in double, leads from `wz_rgf.lead_from_delta(fl(E − E_perp)·dx²)`.
  - Blocks from `wz_rgf.wz_endpoint_blocks`; S from `wz_rgf.wz_fisher_lee`.
  - Intermediate quantities (W_j, P_j, prefix blocks) come from an instrumented copy that must be bit-identical to
    `wz_endpoint_blocks` at the endpoints (checked; any mismatch invalidates the intermediate audit).
- **B. Original RGF (production):**
  - `stationary.scattering_matrix(E, dev, min_closed=0)` on a `stationary.Device` built directly with these V blocks
    (xs = dx·arange(N), V_lead_max = 0). Constructing the frozen dataclass is a read-only use.
  - This allows N = 1, which `build_device` cannot produce.
  - Leads come from `smatrix_spec.lead_modes`.

Represented lattices (all double inputs taken as exact):
- **L:** H_jj = 2tI + diag(E_perp) + V_j exactly, with exact leads from δ_L = (E − E_perp)dx².
- **R_W:** M/t blocks 2Δ_fl,j − 2I, with exact leads from δ_lead,fl. This is A's lattice.
- **R_orig:** D_j = fl(E·I − fl(2tI + diag(E_perp)) − V_j), as production forms it, with production double σ at the ends.
  This is B's lattice.

**C. References:**
- For each represented lattice, the full NN_y × NN_y matrix M (≤ 12 × 12) is assembled in mpmath directly from the
  block-tridiagonal definition and the lead self-energies. It is inverted with `mp.inverse` at dps = 60.
- Independent cross-checks of C:
  - (a) column-by-column `mp.lu_solve`;
  - (b) the same inverse at dps = 100;
  - (c) the residual ‖M·C − I‖_max at dps = 60.
- Prefix systems M^(j) (no right lead) are inverted the same way. They give the intermediate references for W_j, P_j,
  Ĝ^(j)_j1, Ĝ^(j)_1j and Ĝ^(j)_11.
- **Supplementary only, not used as references:** the W/Z equations and the original recursion transcribed in mpmath and
  run on their own lattices. Agreement between two recursions is not accepted as validation.

Primary comparisons:
- A vs C(R_W);
- B vs C(R_orig);
- representation differences C(R_W) − C(L) and C(R_orig) − C(L).

A vs B is reported but is not a validation.

## 4. Metrics (all predefined)

For a block X and its reference X_ref:
- normwise relative error e_rel = ‖X − X_ref‖_F / ‖X_ref‖_F;
- absolute error e_abs = ‖X − X_ref‖_max;
- entrywise error with a floor: max_ab |ΔX_ab| / max(|X_ref,ab|, 1e-6·‖X_ref‖_max). Entries below the floor are
  therefore measured relative to the block scale, which avoids division by near-zero references.

These are applied to:
- the four endpoint blocks Ĝ_11, Ĝ_N1, Ĝ_1N, Ĝ_NN;
- the intermediate W_j (absolute and normwise), P_j, Ĝ^(j)_j1, Ĝ^(j)_1j and Ĝ^(j)_11.

Scattering:
- ΔS_max = ‖S − S_ref‖_max over all open-channel blocks (absolute, since |S_ab| ≤ 1).
- Relative errors only for entries with |S_ref,ab| ≥ 1e-8.
- Unitarity U = max(‖S†S − I‖_max, ‖SS† − I‖_max).
- Reciprocity Q = ‖S − Sᵀ‖_max. For class H, Q_ref ≠ 0 is physical, and only |Q − Q_ref| and ΔS_max are judged.
- Symmetry residual ‖Ĝ_1N − Ĝ_N1ᵀ‖_max for classes R and D.

Residuals and conditioning:
- ‖I − (I − W_j)P_j‖_max, ‖I − (W_1 − I)Ĝ^(1)_11‖_max and ‖I − Z_N Z_N⁻¹‖_max in double.
- The mp residual of C.
- cond₂ of M (from C), I − W_j, W_1 − I and Z_N. Production min_rcond.

Non-commutativity:
- κ_W = max_j ‖[W_j, W_{j+1}]‖_F / (‖W_j‖_F ‖W_{j+1}‖_F) for j + 1 ≤ N − 1.
- κ_chain = max_j ‖[P_j, Ĝ_{j−1,1}]‖_F / (‖P_j‖_F ‖Ĝ_{j−1,1}‖_F).
- κ_Z = ‖[Z_N⁻¹, Ĝ_{N−1,1}]‖_F / (‖Z_N⁻¹‖_F ‖Ĝ_{N−1,1}‖_F).
- ‖[Δ_j, Λ]‖_F.

Reference agreement: C vs (a), (b), (c), for the full matrix and every prefix system.

## 5. Mutants: tests that must fail for ordering or stale-value bugs

Each mutant is a copy of `wz_endpoint_blocks` with exactly one change. Mutants live only in the diagnostic module and are
never used for results.

| id | change | exercised for | detectable when | expected in class D |
|---|---|---|---|---|
| M1 | left chain reversed: Ĝ_j1 ← Ĝ_j1 @ P | N ≥ 3 | κ_chain > 0 | invisible (exactly 0) |
| M2 | right chain reversed: Ĝ_1j ← P @ Ĝ_1j | N ≥ 3 | κ_chain > 0 | invisible |
| M3 | both chains reversed (M1 + M2) | N ≥ 3 | κ_chain > 0; **preserves Ĝ_1N = Ĝ_N1ᵀ in class R, so the reciprocity Q is blind to it** | invisible |
| M4 | Ĝ_11 update uses the NEW chain values (chains updated first) | N ≥ 3 | always (adds extra P factors) | detected |
| M5 | Ĝ_11 update order swapped: Ĝ_j1 P Ĝ_1j | N ≥ 3 | κ_chain > 0 | invisible |
| M6 | closing Ĝ_N1 reversed: −Ĝ_{N−1,1} Z⁻¹ | N ≥ 2 | κ_Z > 0 | invisible |
| M7 | closing Ĝ_1N reversed: −Z⁻¹ Ĝ_{1,N−1} | N ≥ 2 | κ_Z > 0 | invisible |
| M8 | closing Ĝ_11 order swapped: Ĝ_{N−1,1} Z⁻¹ Ĝ_{1,N−1} | N ≥ 2 | κ_Z > 0 | invisible |
| M9 | stale P: P = (I − W_{j−1})⁻¹ (before the W update) | N ≥ 3 | W_j ≠ W_{j−1} | detected (V_j differ) |
| M10 | Z_N missing the (Λ − I) term | N ≥ 2 | always | detected |
| M11 | N = 1 boundary 2Δ + (Λ − I) instead of 2(Λ − I) | N = 1 | always | detected |
| M12 | W update W (I − W)⁻¹ instead of (I − W)⁻¹ W | N ≥ 2 | **never: algebraically equivalent (commuting functions of one matrix)**. This is a negative control | not detected |
| M13 | sign error in the Ĝ_11 chain update (+ instead of −) | N ≥ 3 | always | detected |
| M14 | Ĝ_11 initialization (I − W_1)⁻¹ instead of (W_1 − I)⁻¹ | N ≥ 2 | always | detected |

The mutant discrepancy D_mut is the largest normwise relative error of the four endpoint blocks against C(R_W). In the
intermediate audit, the first slice j at which a prefix quantity deviates locates the bug.

How ordering bugs are distinguished from ordinary floating-point error (pre-registered):
1. **Magnitude gap.** The correct code is expected at ≤ ~1e-14, while ordering mutants are expected at
   ≥ ~1e-3 × (κ/1e-2). A discrepancy > 1e-10 is never rounding at these condition numbers.
2. **Commutator dependence.** For ordering mutants (M1–M3, M5–M8), D_mut must vanish in class D (κ = 0 exactly; diagonal
   products commute bitwise, so D_mut is expected to be exactly 0) and be large in classes R and H. Rounding errors do not
   switch off with κ.
3. **Coupling scaling.** In the s-sub-matrix, ordering-mutant D_mut must fall with s, roughly ∝ s or s² (log-log slope
   between 0.8 and 2.2). The correct implementation's error must show no such trend and stay at rounding level for every s.
4. **Locality.** The intermediate audit must show the first deviating prefix quantity at the slice and block the mutant
   changed (e.g. M1: Ĝ^(2)_21 at N ≥ 3).
5. **Structure.** M3 must be invisible to Q in class R but visible in ΔS_max and the block errors (and in class H through
   Q − Q_ref). This shows that reciprocity alone is not a sufficient ordering test.

## 6. Predictions

Algebraic equivalence (expected exact):
- **E1.** C via `mp.inverse`, via `mp.lu_solve` and at dps 100 agree to ≤ 1e-55 normwise. The mp residual is ≤ 1e-55.
- **E2.** The mp W/Z recursion on R_W and the mp original recursion on R_orig agree with C on their lattices to ≤ 1e-50,
  for every block and every case including N = 1 and N = 2 (supplementary).
- **E3.** Every prefix identity of §1 (W_j from [(M^(j))⁻¹]_jj, left/right chains, old-value Ĝ_11) holds in mp to ≤ 1e-50.

Numerical error scales (double):
- **N1.** A vs C(R_W): normwise errors of all four endpoint blocks are 1e-16 to 1e-14, and below
  τ = max(1e-13, 10·N·N_y·cond₂(M)·ε). ΔS_max ≤ τ.
  - U ≤ τ for every class.
  - Q ≤ τ for classes R and D; |Q − Q_ref| ≤ τ for class H.
  - Symmetry residual ≤ τ for classes R and D. It is not expected to be exactly 0, because the left and right chains
    round differently once the blocks are non-diagonal.
- **N2.** B vs C(R_orig): same scale and bound τ. Production was previously validated with disorder (B3), so no failure is
  expected.
- **N3.** Representation differences C(R_W) − C(L) and C(R_orig) − C(L): normwise ≤ 10·cond₂(M)·ε (≈ 1e-15 to 1e-14).
- **N4.** cond₂(M) is O(1) to O(1e2) for every case. Any case above 1e4 is flagged before interpretation.
  cond₂(I − W_j), cond₂(W_1 − I) and cond₂(Z_N) are O(1) to O(1e2).
- **N5.** Non-commutativity: in classes R and H, κ_chain ≥ 1e-2 (N ≥ 3) and κ_Z ≥ 1e-2 (N ≥ 2). In class D every κ is
  exactly 0. κ_W ≥ 1e-2 for N ≥ 3 in classes R and H (adjacent pairs (W_1, W_2) and, for N = 4, (W_2, W_3)).
- **N6.** Intermediate audit (A's instrumented copy vs mp prefix references): errors ≤ τ at every slice.
- **N7.** Mutants: each mutant marked "detectable" in §5 gives D_mut ≥ 1e-6 in every case where it is exercised and its
  condition holds. Ordering mutants in class D give D_mut = 0 (≤ 1e-15). M12 gives D_mut ≤ τ everywhere.
- **N8.** Coupling scaling: ordering-mutant D_mut decreases monotonically with s (slope 0.8–2.2). The correct-code error
  stays ≤ τ with no monotone trend.

## 7. Pre-registered decision rules

- **Experiment invalid** (no conclusion about W/Z) if any of the following holds:
  - E1 fails;
  - the instrumented copy is not bit-identical to `wz_endpoint_blocks` at the endpoints;
  - a class R or H case meant to test ordering has κ_chain < 1e-3 (N ≥ 3) or κ_Z < 1e-3 (N ≥ 2). That case is then
    reported as "does not test ordering".
- **W/Z ordering FAILS** if, in any case, any endpoint block, any intermediate prefix quantity, or ΔS_max of A vs C(R_W)
  exceeds **1e-10**.
- **Anomalous** if an error lies between τ and 1e-10. It is reported and investigated, not counted as a pass.
- **Test design insufficient** (not a W/Z failure) if any of the following holds:
  - a mutant predicted detectable has D_mut < 1e-6;
  - M12 is flagged (D_mut > 1e-10);
  - an ordering mutant is visible in class D (D_mut > 1e-12).
- **W/Z ordering PASSES at this scale** only if all of the following hold:
  - the experiment is valid;
  - no FAIL or anomaly;
  - the design is sufficient;
  - N1, N6 and N8 hold.

  Even then the claim is limited to "the W/Z recursion reproduces the exact Green blocks of these 76 small deterministic
  coupled devices". **This is not a validation of W/Z**: it covers no long devices, no near-threshold channels, no
  disorder ensembles and no B1 rows.
- **Production (B) vs C(R_orig)** above 1e-10 would be reported as a finding about the production solver. No change would
  be made.

## 8. Planned implementation (not yet written)

- `experiments/wz_rgf_diagnostic/coupled_blocks.py` (driver).
- `experiments/wz_rgf_diagnostic/wz_mutants.py` (instrumented faithful copy and the mutants M1–M14).
- Outputs under `experiments/wz_rgf_diagnostic/results_coupled_blocks/`: `coupled_blocks_results.json`,
  `coupled_blocks_comparisons.csv`, and the report `coupled_blocks_report.md`.

Read-only imports: `stationary`, `smatrix_spec`, `wz_rgf`. Expected cost: seconds.
