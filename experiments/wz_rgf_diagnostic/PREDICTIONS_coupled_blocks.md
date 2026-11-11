# Predictions and design: non-commuting coupled-block validation of the W/Z recursion
(Revision 1. Pre-implementation audit of the version committed in `c2e4d40`. Written and committed BEFORE any implementation
or execution of this experiment.)

Status:
- B1 configuration F still FAILS the unchanged 1e-10 criterion.
- W/Z is experimental and NOT validated.
- No production change is authorized.
- Earlier W/Z tests (`REPORT.md`, `amplitude_roundoff_report.md`) used exactly diagonal clean blocks, so **no non-commuting
  multiplication has been tested yet**.
- Nothing in this document has been implemented or run. The only computation is the design preflight of §3: exact
  quantities of the defined test matrices from direct mpmath inverses. It runs no W/Z code, no original-RGF code and no
  mutant.

## 0. Changes relative to `c2e4d40` (audit outcome)

| # | change | reason |
|---|---|---|
| R1 | Counts corrected: **76 executions, 72 unique configurations** (4 scaling executions at s = 1 duplicate main-matrix configurations) | count error |
| R2 | §3 preflight added: commutators and conditioning for all 72 unique configurations | the κ ≥ 1e-2 bounds were previously unsupported |
| R3 | N5 restated with the measured minima. The κ-validity rule now applies to the main matrix only (the scaling cases are small by design) | preflight |
| R4 | Error bounds split by quantity (§5): endpoint, prefix, W/P, S, residuals. cond₂(M) is no longer used as a universal bound | audit item 3 |
| R5 | **M5 needs N ≥ 4** (identical to the correct code at N = 3), and **M8 needs N ≥ 3** (identical at N = 2). Their detection conditions are corrected. M5, M8, M3 and the new M15 are reciprocity-blind in class R | audit item 4; algebra confirmed by the preflight |
| R6 | Mutants added: M4L and M4R (partial new-value reuse), M15 (M6 + M7), M16 (closing Ĝ₁₁ built from the new endpoint blocks). M12 implementation fixed (right division by `solve` on transposes) | audit item 4 |
| R7 | Reference assembly specified by one routine in physical units for all lattices. The production lattice is taken from production's own input formation. Lead terms are added separately at slices 1 and N (also for N = 1) | audit item 5 |
| R8 | The 4 duplicated s = 1 executions become a determinism check (bitwise identity predicted) | R1 |

**Documentation gap (unchanged):** the `wz_rgf.py` docstring cites "REPORT.md section 1" for its equations, but no
committed W/Z derivation document exists. §1 provides the derivation needed here. Existing files are not edited.

## 1. Mathematical basis (to be tested, not assumed)

Device matrix, as in `stationary.rgf_endpoint_blocks`:
- t = 1/(2dx²).
- H_jj = 2t I + diag(E_perp) + V_j, and H_{j,j±1} = −t I.
- M = E − H_D − Σ_L ⊕ Σ_R, with Σ = diag(−tλ_n) on slice 1 (Σ_L) and on slice N (Σ_R). For N = 1, both are on the same
  slice.
- The leads are clean: λ_n + 1/λ_n = 2(1 − δ_n), δ_n = (E − E_perp,n)dx², taking the retarded root.

W/Z variables:
- Δ_j = (E I − diag(E_perp) − V_j)dx².
- Since 2t dx² = 1, M/t has diagonal blocks 2Δ_j − 2I (+Λ at slices 1 and N) and off-diagonal blocks I.
- Ĝ = tG.

Prefix (left-connected) systems. For 1 ≤ j ≤ N − 1, M^(j) is the leading j-slice block submatrix of M/t with the right
lead omitted, i.e. only the left Λ on slice 1. Then:
- P_j = (I − W_j)⁻¹ = −[(M^(j))⁻¹]_jj and W_j = I + ([(M^(j))⁻¹]_jj)⁻¹;
- W_1 = 2Δ_1 + (Λ − I);
- W_j = 2Δ_j + (I − W_{j−1})⁻¹W_{j−1};
- Ĝ^(j)_j1 = P_j Ĝ^(j−1)_{j−1,1};
- Ĝ^(j)_1j = Ĝ^(j−1)_{1,j−1} P_j;
- Ĝ^(j)_11 = Ĝ^(j−1)_11 − Ĝ^(j−1)_{1,j−1} P_j Ĝ^(j−1)_{j−1,1};
- Ĝ^(1)_11 = (W_1 − I)⁻¹.

Closing step:
- Z_N = 2Δ_N + (Λ − I) + (I − W_{N−1})⁻¹W_{N−1}.
- Ĝ_NN = Z_N⁻¹, Ĝ_N1 = −Z_N⁻¹ Ĝ^(N−1)_{N−1,1} and Ĝ_1N = −Ĝ^(N−1)_{1,N−1} Z_N⁻¹.
- Ĝ_11 = Ĝ^(N−1)_11 + Ĝ^(N−1)_{1,N−1} Z_N⁻¹ Ĝ^(N−1)_{N−1,1}.

Edge cases:
- N = 1: Z_1 = 2Δ_1 + 2(Λ − I), i.e. both leads exactly once each, and all four blocks equal Z_1⁻¹.
- N = 2: the loop is empty, and Z_2 is the Schur complement of [[W_1 − I, I], [I, 2Δ_2 − 2I + Λ]].

The preflight (§3) confirms all of these identities at 60 digits for every configuration: the prefix chains, the old-value
Ĝ₁₁ update, the closing Z_N against [(M/t)⁻¹]_NN⁻¹, and the N = 1 form. The largest deviation is 3.5e-61.

Fisher–Lee (open channels, s_n = Γ_n/(2t)):
- t_RL = 2i √s Ĝ_N1 √s and t_LR = 2i √s Ĝ_1N √s;
- r_L = 2i √s Ĝ_11 √s − I and r_R = 2i √s Ĝ_NN √s − I.
These are elementwise, so they are valid for coupled blocks.

Structure used below:
- (i) For class R (real symmetric V_j), M is complex symmetric. Then Ĝ = Ĝᵀ, Ĝ^(j)_1j = (Ĝ^(j)_j1)ᵀ, and W_j, P_j, Ĝ^(j)_11
  and Z_N are symmetric.
- (ii) Reversing a product changes the result by a commutator term, which vanishes exactly for diagonal blocks.
- (iii) All blocks are N_y × N_y, because both leads have the same width. Swapped triple products such as those of M5 and
  M8 are therefore always dimensionally valid here, but they are index-inconsistent: they compose the paths 1 → j → j → 1
  in the wrong order. In a device with unequal lead widths they would not even be conformable. This design cannot test
  that.

## 2. Test devices

Lattice and leads:
- dx = 1/8 (t = 32, exact), L_y = 1, clean leads.
- E_perp,n = (nπ)²/2 (doubles from `stationary.transverse_basis(1.0, N_y)`): 4.934802200544679, 19.739208802178716,
  44.41321980490211.

Energies (all |c_n| far from 1):

| id | N_y | E | δ_n | open / closed |
|---|---|---|---|---|
| e2a | 2 | 12 | +0.110394, −0.120925 | 1 / 1 |
| e2b | 2 | 30 | +0.391644, +0.160325 | 2 / 0 |
| e3a | 3 | 12 | +0.110394, −0.120925, −0.506457 | 1 / 2 |
| e3b | 3 | 30 | +0.391644, +0.160325, −0.225207 | 2 / 1 |
| e3c | 3 | 50 | +0.704144, +0.472825, +0.087293 | 3 / 0 |

Potentials (unchanged from `c2e4d40`; dyadic, exactly representable; v0 = 4; slices j = 1..4, a device with N slices
uses V_1..V_N).

N_y = 2, with A = diag(1, −1), B = [[0, 1], [1, 0]] and C = [[0, −i], [i, 0]]:

    V_j = v0 s (μ_j I + α_j A + β_j B + γ_j C)
    μ = (0.5, −0.25, 1, 0)   α = (1, 2, −1.5, 0.5)   β = (3, −1, 2.5, 1.75)   γ = (1, 0.5, −2, 1.25)

N_y = 3:

    V_j = v0 s (μ_j I + α_j diag(1,0,−1) + α'_j diag(1,−2,1) + β12_j S12 + β23_j S23 + β13_j S13 + γ12_j K12 + γ23_j K23)

- S_ab = E_ab + E_ba and K_ab = i(E_ba − E_ab).
- Coefficients:
  - μ = (0.5, −0.25, 1, 0); α = (1, 2, −1.5, 0.5); α' = (0.25, −0.5, 0.75, 1);
  - β12 = (3, −1, 2.5, 1.75); β23 = (−2, 1.5, 0.5, −1); β13 = (1, 0, −1.25, 2);
  - γ12 = (1, 0.5, −2, 1.25); γ23 = (−0.5, 1, 0.75, −1.5).

Potential classes:
- **R** (γ = 0; real symmetric; time-reversal symmetric; the physically relevant class).
- **H** (all coefficients; complex Hermitian; time-reversal broken). An algebraic stress test, not physical for this
  project. `wz_rgf.py` documents real symmetric Δ, so class H tests the code beyond its documented scope.
- **D** (only μ, α, α'; diagonal control).

Case counts:

| matrix | configurations | executions |
|---|---|---|
| main: 5 energy sets × 3 classes × N ∈ {1, 2, 3, 4} (s = 1) | 60 | 60 |
| scaling: class R, {e2b, e3b} × N ∈ {3, 4} × s ∈ {1, 2⁻⁴, 2⁻⁸, 2⁻¹²} | 16, of which 4 (s = 1) duplicate main configurations | 16 |
| **total** | **72 unique** | **76** |

The 4 duplicated executions must reproduce their main-matrix counterparts bit-for-bit. This is a determinism check, not
additional coverage.

## 3. Design preflight: exact commutators and conditioning (computed; not an experiment run)

Method:
- mpmath at 60 digits, on the represented W/Z lattice R_W (double Δ_fl,j and double lead δ taken as exact, exact roots).
- Full matrix M/t and prefix matrices M^(j) are assembled from the block-tridiagonal definition and inverted directly.
- W_j, P_j and the prefix blocks are taken from the prefix inverses, Z_N from [(M/t)⁻¹]_NN⁻¹.
- κ_X = ‖[A, B]‖_F / (‖A‖_F ‖B‖_F):
  - κ_W uses (W_j, W_{j+1}) for j + 1 ≤ N − 1;
  - κ_chain uses (P_j, Ĝ^(j−1)_{j−1,1}) for 2 ≤ j ≤ N − 1;
  - κ_Z uses (Z_N⁻¹, Ĝ^(N−1)_{N−1,1}) for N ≥ 2;
  - κ_[Δ,Λ] uses (Δ_j, Λ), maximum over j.
- "M5 cond." = ‖X − Y‖_F/‖X‖_F with X = Ĝ_{1,j−1}P_jĜ_{j−1,1} and Y = Ĝ_{j−1,1}P_jĜ_{1,j−1} at the step j where M5 first
  differs (j = 3, so N = 4 only).
- "M8 cond." = the same quantity with Z_N⁻¹ at the closing step.
- cond₂ is computed from the double rounding of the exact matrices.
- Script: Appendix A (sha256 `74d2d54e…284172`). Output sha256 `59b844e1…bf0afa`. Neither file is committed; Appendix A
  reproduces the script.

| case | open | cond₂(M) | max cond₂(M^(j)) | max cond₂(I−W_j) | cond₂(Z_N) | max‖P_j‖₂ | κ_W (j=1,2) | κ_chain (j=2,3) | κ_Z | max κ_[Δ,Λ] | M5 cond. (j=3) | M8 cond. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| e2a/R/N1 | 1 | 1.81 | — | — | 1.81 | — | — | — | — | 4.3e-01 | — | — |
| e2a/R/N2 | 1 | 7.12 | 1.9 | 1.87 | 1.34 | 1.05 | — | — | 1.7e-01 | 4.4e-01 | — | 0 |
| e2a/R/N3 | 1 | 11.71 | 6.6 | 1.87 | 1.83 | 1.05 | 2.2e-01 | 2.9e-02 | 1.4e-01 | 4.4e-01 | — | 5.9e-02 |
| e2a/R/N4 | 1 | 17.69 | 14.9 | 2.17 | 1.50 | 1.14 | 2.2e-01, 1.3e-01 | 2.9e-02, 8.1e-02 | 1.3e-01 | 4.4e-01 | 4.2e-02 | 1.3e-01 |
| e2a/H/N1 | 1 | 1.87 | — | — | 1.87 | — | — | — | — | 4.3e-01 | — | — |
| e2a/H/N2 | 1 | 7.55 | 1.9 | 1.93 | 1.45 | 1.07 | — | — | 1.8e-01 | 4.5e-01 | — | 0 |
| e2a/H/N3 | 1 | 11.71 | 7.0 | 1.93 | 1.95 | 1.07 | 3.2e-01 | 4.6e-02 | 3.0e-01 | 4.5e-01 | — | 8.3e-02 |
| e2a/H/N4 | 1 | 17.52 | 15.5 | 2.28 | 1.46 | 1.19 | 3.2e-01, 5.9e-01 | 4.6e-02, 1.8e-01 | 2.0e-01 | 4.5e-01 | 6.0e-02 | 2.5e-01 |
| e2a/D/N1 | 1 | 1.02 | — | — | 1.02 | — | — | — | — | 0 | — | — |
| e2a/D/N2 | 1 | 5.48 | 1.3 | 1.34 | 1.21 | 0.85 | — | — | 0 | 0 | — | 0 |
| e2a/D/N3 | 1 | 10.64 | 5.0 | 1.34 | 1.76 | 0.85 | 0.0e+00 | 0.0e+00 | 0 | 0 | — | 0 |
| e2a/D/N4 | 1 | 17.78 | 11.4 | 1.85 | 1.48 | 1.03 | 0.0e+00, 0.0e+00 | 0.0e+00, 0.0e+00 | 0 | 0 | 0.0e+00 | 0 |
| e2b/R/N1 | 2 | 1.45 | — | — | 1.45 | — | — | — | — | 1.4e-01 | — | — |
| e2b/R/N2 | 2 | 3.17 | 1.8 | 1.76 | 1.73 | 1.29 | — | — | 7.1e-02 | 1.4e-01 | — | 0 |
| e2b/R/N3 | 2 | 7.88 | 6.0 | 1.93 | 1.10 | 1.38 | 3.6e-02 | 2.0e-02 | 1.2e-01 | 1.4e-01 | — | 2.8e-02 |
| e2b/R/N4 | 2 | 11.86 | 15.8 | 2.26 | 2.20 | 1.61 | 3.6e-02, 3.1e-01 | 2.0e-02, 2.1e-01 | 2.6e-01 | 1.4e-01 | 3.3e-02 | 5.1e-01 |
| e2b/H/N1 | 2 | 1.45 | — | — | 1.45 | — | — | — | — | 1.5e-01 | — | — |
| e2b/H/N2 | 2 | 3.21 | 1.8 | 1.80 | 1.77 | 1.31 | — | — | 7.9e-02 | 1.5e-01 | — | 0 |
| e2b/H/N3 | 2 | 7.62 | 6.0 | 1.92 | 1.18 | 1.35 | 8.6e-02 | 4.9e-02 | 1.9e-01 | 1.5e-01 | — | 5.7e-02 |
| e2b/H/N4 | 2 | 11.59 | 17.2 | 2.89 | 2.64 | 1.89 | 8.6e-02, 3.8e-01 | 4.9e-02, 2.8e-01 | 3.3e-01 | 1.5e-01 | 8.8e-02 | 6.6e-01 |
| e2b/D/N1 | 2 | 1.47 | — | — | 1.47 | — | — | — | — | 0 | — | — |
| e2b/D/N2 | 2 | 3.09 | 1.2 | 1.19 | 1.21 | 1.05 | — | — | 0 | 0 | — | 0 |
| e2b/D/N3 | 2 | 6.83 | 6.4 | 1.27 | 1.07 | 1.24 | 0.0e+00 | 0.0e+00 | 0 | 0 | — | 0 |
| e2b/D/N4 | 2 | 9.35 | 11.4 | 2.10 | 1.99 | 1.40 | 0.0e+00, 0.0e+00 | 0.0e+00, 0.0e+00 | 0 | 0 | 0.0e+00 | 0 |
| e3a/R/N1 | 1 | 3.31 | — | — | 3.31 | — | — | — | — | 2.1e-01 | — | — |
| e3a/R/N2 | 1 | 8.39 | 2.8 | 2.78 | 2.21 | 1.04 | — | — | 1.9e-01 | 2.1e-01 | — | 0 |
| e3a/R/N3 | 1 | 15.09 | 8.3 | 2.78 | 4.00 | 1.04 | 1.8e-01 | 7.6e-02 | 1.4e-01 | 2.1e-01 | — | 1.5e-01 |
| e3a/R/N4 | 1 | 24.70 | 17.7 | 3.38 | 3.91 | 1.13 | 1.8e-01, 1.4e-01 | 7.6e-02, 1.2e-01 | 1.3e-01 | 2.1e-01 | 1.2e-01 | 1.9e-01 |
| e3a/H/N1 | 1 | 3.36 | — | — | 3.36 | — | — | — | — | 2.2e-01 | — | — |
| e3a/H/N2 | 1 | 8.68 | 2.8 | 2.83 | 2.29 | 1.06 | — | — | 1.9e-01 | 2.2e-01 | — | 0 |
| e3a/H/N3 | 1 | 15.04 | 8.7 | 2.83 | 4.26 | 1.06 | 2.0e-01 | 8.2e-02 | 3.0e-01 | 2.2e-01 | — | 1.5e-01 |
| e3a/H/N4 | 1 | 24.64 | 18.3 | 3.58 | 3.93 | 1.19 | 2.0e-01, 1.6e-01 | 8.2e-02, 1.9e-01 | 1.4e-01 | 2.2e-01 | 1.2e-01 | 2.9e-01 |
| e3a/D/N1 | 1 | 2.37 | — | — | 2.37 | — | — | — | — | 0 | — | — |
| e3a/D/N2 | 1 | 7.47 | 2.2 | 2.17 | 2.32 | 0.83 | — | — | 0 | 0 | — | 0 |
| e3a/D/N3 | 1 | 14.03 | 6.6 | 2.17 | 3.79 | 0.83 | 0.0e+00 | 0.0e+00 | 0 | 0 | — | 0 |
| e3a/D/N4 | 1 | 23.41 | 14.0 | 2.88 | 3.50 | 0.97 | 0.0e+00, 0.0e+00 | 0.0e+00, 0.0e+00 | 0 | 0 | 0.0e+00 | 0 |
| e3b/R/N1 | 2 | 1.57 | — | — | 1.57 | — | — | — | — | 2.1e-01 | — | — |
| e3b/R/N2 | 2 | 4.84 | 2.6 | 2.60 | 1.95 | 1.32 | — | — | 1.2e-01 | 2.1e-01 | — | 0 |
| e3b/R/N3 | 2 | 10.20 | 6.5 | 2.60 | 1.88 | 1.32 | 2.6e-01 | 8.9e-02 | 2.0e-01 | 2.1e-01 | — | 1.9e-01 |
| e3b/R/N4 | 2 | 14.82 | 17.7 | 3.22 | 1.86 | 1.43 | 2.6e-01, 3.4e-01 | 8.9e-02, 2.0e-01 | 2.4e-01 | 2.1e-01 | 1.8e-01 | 4.9e-01 |
| e3b/H/N1 | 2 | 1.57 | — | — | 1.57 | — | — | — | — | 2.1e-01 | — | — |
| e3b/H/N2 | 2 | 4.83 | 2.6 | 2.64 | 2.03 | 1.34 | — | — | 1.3e-01 | 2.1e-01 | — | 0 |
| e3b/H/N3 | 2 | 9.84 | 6.5 | 2.64 | 1.83 | 1.34 | 2.9e-01 | 1.0e-01 | 2.6e-01 | 2.1e-01 | — | 2.0e-01 |
| e3b/H/N4 | 2 | 14.26 | 20.8 | 4.45 | 2.41 | 1.97 | 2.9e-01, 3.5e-01 | 1.0e-01, 3.2e-01 | 3.4e-01 | 2.1e-01 | 1.9e-01 | 7.3e-01 |
| e3b/D/N1 | 2 | 1.48 | — | — | 1.48 | — | — | — | — | 0 | — | — |
| e3b/D/N2 | 2 | 4.97 | 1.9 | 1.90 | 1.32 | 1.00 | — | — | 0 | 0 | — | 0 |
| e3b/D/N3 | 2 | 8.93 | 6.8 | 1.90 | 1.69 | 1.04 | 0.0e+00 | 0.0e+00 | 0 | 0 | — | 0 |
| e3b/D/N4 | 2 | 11.15 | 14.7 | 2.81 | 1.65 | 1.27 | 0.0e+00, 0.0e+00 | 0.0e+00, 0.0e+00 | 0 | 0 | 0.0e+00 | 0 |
| e3c/R/N1 | 3 | 2.32 | — | — | 2.32 | — | — | — | — | 9.9e-02 | — | — |
| e3c/R/N2 | 3 | 3.72 | 1.6 | 1.62 | 2.34 | 1.26 | — | — | 7.7e-02 | 9.9e-02 | — | 0 |
| e3c/R/N3 | 3 | 9.41 | 7.8 | 2.97 | 2.05 | 1.81 | 1.4e-01 | 9.2e-02 | 2.5e-01 | 9.9e-02 | — | 2.1e-01 |
| e3c/R/N4 | 3 | 13.35 | 14.4 | 2.97 | 2.70 | 1.81 | 1.4e-01, 2.8e-01 | 9.2e-02, 3.8e-01 | 5.0e-01 | 1.0e-01 | 2.9e-01 | 1.8e+00 |
| e3c/H/N1 | 3 | 2.33 | — | — | 2.33 | — | — | — | — | 1.0e-01 | — | — |
| e3c/H/N2 | 3 | 3.76 | 1.7 | 1.66 | 2.33 | 1.26 | — | — | 8.8e-02 | 1.0e-01 | — | 0 |
| e3c/H/N3 | 3 | 9.40 | 7.9 | 3.12 | 2.09 | 1.87 | 1.6e-01 | 1.0e-01 | 2.6e-01 | 1.0e-01 | — | 2.4e-01 |
| e3c/H/N4 | 3 | 12.94 | 13.9 | 3.12 | 2.60 | 1.87 | 1.6e-01, 2.9e-01 | 1.0e-01, 3.9e-01 | 5.1e-01 | 1.1e-01 | 3.1e-01 | 1.7e+00 |
| e3c/D/N1 | 3 | 2.35 | — | — | 2.35 | — | — | — | — | 0 | — | — |
| e3c/D/N2 | 3 | 3.58 | 1.1 | 1.12 | 1.93 | 1.03 | — | — | 0 | 0 | — | 0 |
| e3c/D/N3 | 3 | 10.06 | 7.4 | 1.51 | 1.55 | 1.43 | 0.0e+00 | 0.0e+00 | 0 | 0 | — | 0 |
| e3c/D/N4 | 3 | 14.74 | 15.7 | 1.61 | 1.98 | 1.43 | 0.0e+00, 0.0e+00 | 0.0e+00, 0.0e+00 | 0 | 0 | 0.0e+00 | 0 |
| e2b/R/N3/s=2^-4 | 2 | 7.07 | 5.7 | 1.04 | 1.43 | 1.03 | 1.5e-02 | 8.5e-03 | 8.7e-03 | 9.5e-03 | — | 8.8e-03 |
| e2b/R/N3/s=2^-8 | 2 | 7.06 | 5.6 | 1.00 | 1.46 | 1.00 | 1.0e-03 | 5.6e-04 | 5.6e-04 | 5.9e-04 | — | 5.9e-04 |
| e2b/R/N3/s=2^-12 | 2 | 7.06 | 5.6 | 1.00 | 1.46 | 1.00 | 6.4e-05 | 3.5e-05 | 3.5e-05 | 3.7e-05 | — | 3.7e-05 |
| e2b/R/N4/s=2^-4 | 2 | 10.11 | 11.9 | 1.04 | 1.49 | 1.03 | 1.5e-02, 1.2e-02 | 8.5e-03, 9.8e-03 | 1.3e-02 | 9.5e-03 | 1.1e-02 | 1.9e-02 |
| e2b/R/N4/s=2^-8 | 2 | 10.14 | 11.9 | 1.00 | 1.46 | 1.00 | 1.0e-03, 7.5e-04 | 5.6e-04, 6.1e-04 | 8.2e-04 | 5.9e-04 | 7.8e-04 | 1.2e-03 |
| e2b/R/N4/s=2^-12 | 2 | 10.14 | 11.9 | 1.00 | 1.46 | 1.00 | 6.4e-05, 4.7e-05 | 3.5e-05, 3.8e-05 | 5.1e-05 | 3.7e-05 | 4.9e-05 | 7.5e-05 |
| e3b/R/N3/s=2^-4 | 2 | 9.23 | 7.6 | 1.96 | 1.47 | 1.02 | 2.2e-02 | 9.0e-03 | 1.2e-02 | 1.4e-02 | — | 1.3e-02 |
| e3b/R/N3/s=2^-8 | 2 | 9.25 | 7.7 | 1.93 | 1.46 | 1.00 | 1.4e-03 | 5.7e-04 | 7.6e-04 | 8.7e-04 | — | 8.6e-04 |
| e3b/R/N3/s=2^-12 | 2 | 9.25 | 7.7 | 1.93 | 1.46 | 1.00 | 8.5e-05 | 3.6e-05 | 4.7e-05 | 5.4e-05 | — | 5.4e-05 |
| e3b/R/N4/s=2^-4 | 2 | 12.69 | 15.1 | 1.97 | 1.48 | 1.02 | 2.2e-02, 1.6e-02 | 9.0e-03, 1.2e-02 | 1.4e-02 | 1.4e-02 | 1.4e-02 | 2.5e-02 |
| e3b/R/N4/s=2^-8 | 2 | 12.78 | 15.1 | 1.94 | 1.46 | 1.00 | 1.4e-03, 1.0e-03 | 5.7e-04, 7.2e-04 | 8.7e-04 | 8.7e-04 | 9.2e-04 | 1.6e-03 |
| e3b/R/N4/s=2^-12 | 2 | 12.79 | 15.1 | 1.93 | 1.46 | 1.00 | 8.5e-05, 6.3e-05 | 3.6e-05, 4.5e-05 | 5.4e-05 | 5.4e-05 | 5.8e-05 | 9.9e-05 |

Preflight conclusions:
- **Ordering conditions are exercised in every main-matrix R and H case where the corresponding code path exists.**
  Minima over these cases:
  - κ_chain = 2.0e-2 (e2b/R/N3);
  - κ_W = 3.6e-2 (e2b/R);
  - κ_Z = 7.1e-2 (e2b/R/N2);
  - κ_[Δ,Λ] = 9.9e-2 (e3c/R);
  - M5 condition = 3.3e-2 (N = 4);
  - M8 condition (N ≥ 3) = 2.8e-2.
  - The previously claimed lower bound of 1e-2 is therefore **supported, but only with a margin of 2× for κ_chain**. The
    bound is kept and no matrix is changed.
- **Class D: every commutator and every M5/M8 condition is exactly 0** in 60-digit arithmetic, as required.
- **N = 2: the M8 condition is exactly 0 in every class**. At N = 2, Ĝ_{N−1,1} = Ĝ_{1,N−1} = Ĝ^(1)_11, so M8 is identical to
  the correct code (correction R5).
- Scaling cases: κ falls linearly with s. For κ_chain: 2.0e-2 → 8.5e-3 → 5.6e-4 → 3.5e-5 (e2b/N3), a factor 16 per factor
  2⁻⁴. At s = 2⁻⁸ and 2⁻¹² it lies below the main-matrix validity threshold of 1e-3, by design.
- Conditioning is mild everywhere:
  - cond₂(M) ≤ 24.7;
  - cond₂(M^(j)) ≤ 20.8;
  - cond₂(I − W_j) ≤ 4.45;
  - cond₂(Z_N) ≤ 4.26;
  - ‖P_j‖₂ ≤ 1.97.
- mp residuals of the full inverses are ≤ 3e-64.

## 4. Solvers and references

Solvers under test (double, unmodified code, read-only):
- **A. W/Z:**
  - Δ_j = (E·I − diag(E_perp) − V_j)·dx² in double; leads from `wz_rgf.lead_from_delta(fl(E − E_perp)·dx²)`.
  - Endpoint blocks from `wz_rgf.wz_endpoint_blocks`; S from `wz_rgf.wz_fisher_lee`.
  - Intermediate values come from an instrumented copy that must be bit-identical to `wz_endpoint_blocks` at all four
    endpoint blocks.
- **B. Original RGF:**
  - `stationary.scattering_matrix(E, dev, min_closed=0)` on `stationary.Device(xs=dx·arange(N), dx, basis, V, 0.0)`,
    built directly (read-only use; this allows N = 1).

References. One assembly routine is used for all lattices, in physical units:
- diagonal blocks A_j (= E − H_jj of the lattice) and off-diagonal blocks +tI;
- then Σ_L is subtracted on slice 1 and Σ_R on slice N, **as two separate operations**, so that for N = 1 the single
  block receives both terms exactly once;
- Ĝ_ref = t·M⁻¹ by `mp.inverse` at dps 60.

| lattice | A_j (taken exactly) | Σ_L = Σ_R | Fisher–Lee input | used for |
|---|---|---|---|---|
| L | E − (2t + E_perp + V_j) in exact arithmetic | −t λ_exact(δ_L), δ_L = (E − E_perp)dx² exact | exact s | representation differences |
| R_W | t(2Δ_fl,j − 2I), with Δ_fl,j the exact double array that A receives | −t λ_exact(δ_lead,fl), δ_lead,fl = fl(E − E_perp)·dx² (the exact input of `lead_from_delta`) | exact s from λ_exact | A |
| R_orig | the double array `E*eye − H[j]` with H from `stationary.device_hamiltonian_blocks(dev)`, i.e. production's own formation | production double `sigma` from `smatrix_spec.lead_modes(E, E_perp, dx)` | production double `gamma` | B |

Notes on the references:
- In R_orig the lead/device mismatch at the ulp level, and Γ ≠ 2t·Im λ, are part of the represented production lattice.
  So U_ref and Q_ref of R_orig are not exactly 0, and they are reported.
- In R_W the double lead data (λ − 1, sin θ) are counted as solver rounding, as in the earlier experiments.
- Production's in-loop rounding of D − Σ is solver arithmetic, not represented input.

Prefix references:
- M^(j) for j = 1..N − 1 is built by the same routine with **Σ_R omitted** (left lead on slice 1 only).
- M^(1) = A_1 − Σ_L. No prefix is formed for j = N.

Reference checks (an experiment is invalid if any fails):
- (a) `mp.lu_solve` column by column;
- (b) dps 100;
- (c) mp residual ‖M Ĝ_ref/t − I‖_max ≤ 1e-55;
- (d) the N = 1 assembled block equals A_1 − Σ_L − Σ_R entrywise (guards against double-counting or omitting a lead term);
- (e) C(R_W) − C(L) and C(R_orig) − C(L) are at representation level (§5).

The mp W/Z recursion and the mp original recursion are supplementary only.

## 5. Metrics and error bounds (separated by quantity)

Notation: ε = 2⁻⁵² ≈ 2.2e-16, n = N·N_y, c = 10 (safety constant).
- K_Schur = max(cond₂(W_1 − I), max_j cond₂(I − W_j), cond₂(Z_N)).
- g = max(1, max_j ‖P_j‖₂)².

All are computed from the references, not from the solver. The bounds below are first-order **heuristics for block
elimination without inter-block pivoting**, not theorems. Their only roles are (1) to flag anomalies and (2) to check that
the fixed FAIL threshold 1e-10 is far above legitimate rounding.

1. **Endpoint Green blocks** (Ĝ_11, Ĝ_N1, Ĝ_1N, Ĝ_NN).
   - Primary metric: the error relative to the whole inverse, e_glob = ‖ΔX‖_F / ‖Ĝ_full‖₂, which is the quantity a
     cond₂(M) argument controls.
   - Bound: τ_end = max(1e-13, c·n·ε·max(cond₂(M), K_Schur)·g).
   - The block-own relative error e_rel = ‖ΔX‖_F/‖X_ref‖_F is reported and judged against τ_end·‖Ĝ_full‖₂/‖X_ref‖_F. A
     small block (e.g. a strongly decayed closed-channel Ĝ_N1) is not expected to have a small own-relative error.
   - The entrywise error uses the floor max(|X_ref,ab|, 1e-6·‖X_ref‖_max).
   - Preflight maxima: τ_end ≤ max(1e-13, 10·12·ε·24.7·3.9) ≈ 2.6e-12.
2. **Prefix blocks** Ĝ^(j)_11, Ĝ^(j)_j1, Ĝ^(j)_1j.
   - Each is judged with **its own** conditioning: e_glob,j relative to ‖(M^(j))⁻¹‖₂ and
     τ_pre,j = max(1e-13, c·j·N_y·ε·max(cond₂(M^(j)), K_Schur,≤j)·g_≤j).
   - cond₂(M) of the full system is not used here: a prefix system lacks the right lead and can be conditioned
     differently. Preflight: cond₂(M^(j)) ≤ 20.8.
3. **W_j and P_j** (absolute errors; W_j is O(0.1–1)).
   - W errors propagate as ΔW_j ≈ P_{j−1} ΔW_{j−1} P_{j−1} + local rounding (since d[(I − W)⁻¹W] = P dW P).
   - Bound: τ_W,j = c·N_y·ε·j·g^(j−1)·(1 + K_Schur)·max(1, max_k ≤ j ‖W_k‖₂).
   - P bound: τ_P,j = ‖P_j‖₂²·τ_W,j + c·N_y·ε·cond₂(I − W_j)·‖P_j‖₂.
   - Neither uses cond₂(M). Preflight scale: ≲ 1e-12.
4. **Scattering matrix** (all open channels; |S_ab| ≤ 1, so absolute errors).
   - τ_S = max(1e-13, 2·max_n s_n·‖Ĝ_full‖₂·τ_end) + 4ε. The 4ε term covers the lead-data rounding of s.
   - Unitarity: U ≤ U_ref + 2·(2n_open)·τ_S.
   - Reciprocity: Q ≤ Q_ref + 2τ_S (classes R and D). For class H, Q_ref ≠ 0 is physical, and |Q − Q_ref| ≤ 2τ_S is required
     instead.
   - Relative errors are reported only for entries with |S_ref,ab| ≥ 1e-8.
   - Symmetry residual ‖Ĝ_1N − Ĝ_N1ᵀ‖_max ≤ 2τ_end·‖Ĝ_full‖₂ in classes R and D. It is not expected to be exactly 0.
5. **Residuals** of each double solve: ρ = ‖I − A X‖_F / (‖A‖_F ‖X‖_F) for (W_1 − I), (I − W_j) and Z_N.
   - For a backward-stable LU solve, ρ ≲ c·N_y·ε **independently of conditioning**.
   - ρ > 100·N_y·ε is an anomaly, not a FAIL. These are diagnostics of the linear solves only and say nothing about
     ordering.
6. **Non-commutativity:** κ_W, κ_chain, κ_Z and κ_[Δ,Λ] as in §3, recomputed in the experiment from the references.
7. **Reference agreement:** checks (a)–(e) of §4.
8. **Mutant discrepancy:** D_mut = max over the four endpoint blocks of e_glob against C(R_W), plus the first deviating
   prefix quantity in the intermediate audit.

FAIL threshold: unchanged at 1e-10 for e_glob (endpoint and prefix), for absolute W/P errors and for ΔS_max. Every τ in
the preflight is ≤ ~3e-12, so 1e-10 lies ≥ 30× above every legitimate rounding bound. If an experiment-time τ exceeds
1e-11 for some case, that case is reported as "insufficiently conditioned for the fixed threshold" instead of being judged.

## 6. Mutants (each a copy of `wz_endpoint_blocks` with exactly one change; diagnostic module only)

Correct code, for reference:
- Loop, slices 2..N−1: W ← 2Δ_j + solve(I − W, W); P ← solve(I − W, I); G11 ← G11 − G1j·P·Gj1; Gj1 ← P·Gj1; G1j ← G1j·P.
- Closing: Z ← 2Δ_N + Lm + solve(I − W, W); Zi ← solve(Z, I); G11 ← G11 + G1j·Zi·Gj1; GN1 ← −Zi·Gj1; G1N ← −G1j·Zi.
- Initialization: G11 ← solve(W_1 − I, I), Gj1 = G1j = copies of G11.

| id | exact change | min N | first differing quantity | detection condition (exact) | class D | Q-blind in class R? |
|---|---|---|---|---|---|---|
| M1 | Gj1 ← Gj1·P | 3 | Ĝ^(2)_21 | [P_2, Ĝ^(1)_11] ≠ 0 ⇔ [W_2, W_1] ≠ 0 (κ_chain) | invisible | no (breaks Ĝ_1N = Ĝ_N1ᵀ) |
| M2 | G1j ← P·G1j | 3 | Ĝ^(2)_12 | same as M1 | invisible | no |
| M3 | M1 + M2 | 3 | Ĝ^(2)_21 and Ĝ^(2)_12 (Ĝ^(2)_11 unchanged at j = 2) | κ_chain ≠ 0; at N = 3 it is seen only through the closing step | invisible | **yes**: maps each chain to the transpose of the other, which keeps S = Sᵀ |
| M4 | both chains updated before the G11 update: G11 − Ĝ_{1,j−1}P³Ĝ_{j−1,1} | 3 | Ĝ^(2)_11 | P_j(P_j² − I) ≠ 0 (P_j² ≠ I; true here) | detected | no |
| M4L | only Gj1 updated before G11: G11 − Ĝ_{1,j−1}P²Ĝ_{j−1,1} | 3 | Ĝ^(2)_11 | P_j(P_j − I) ≠ 0 (W_j ≠ 0) | detected | no (asymmetric) |
| M4R | only G1j updated before G11 (mirror of M4L) | 3 | Ĝ^(2)_11 | as M4L | detected | no |
| M5 | G11 ← G11 − Gj1·P·G1j | **4** | Ĝ^(3)_11 | Ĝ_{1,2}P_3Ĝ_{2,1} ≠ Ĝ_{2,1}P_3Ĝ_{1,2} ("M5 cond.", ≥ 3.3e-2). **At N = 3 it is identical to the correct code**: at j = 2, Gj1 = G1j = Ĝ^(1)_11 | invisible | **yes** (the swapped term is also symmetric) |
| M6 | GN1 ← −Gj1·Zi | 2 | Ĝ_N1 | [Z_N⁻¹, Ĝ_{N−1,1}] ≠ 0 (κ_Z). At N = 2 this needs Δ_1 ≠ Δ_2 (identical slices would commute) | invisible | no |
| M7 | G1N ← −Zi·G1j | 2 | Ĝ_1N | κ_Z ≠ 0 (as M6) | invisible | no |
| M8 | G11 ← G11 + Gj1·Zi·G1j | **3** | Ĝ_11 | Ĝ_{1,N−1}Z⁻¹Ĝ_{N−1,1} ≠ Ĝ_{N−1,1}Z⁻¹Ĝ_{1,N−1} ("M8 cond.", ≥ 2.8e-2 at N ≥ 3). **At N = 2 it is identical** (both chains equal Ĝ^(1)_11; preflight condition exactly 0). Dimensionally valid only because both leads have N_y channels | invisible | **yes** |
| M9 | stale P: P ← solve(I − W_old, I) computed before the W update | 3 | P_2, then Ĝ^(2)_* | W_2 ≠ W_1 (true in every class; V_j differ) | detected | no |
| M10 | Z ← 2Δ_N + solve(I − W, W) (Lm dropped) | 2 | Z_N, hence all four endpoint blocks | Λ ≠ I (always) | detected | no |
| M11 | N = 1: Z ← 2Δ_1 + Lm instead of 2Δ_1 + 2Lm | 1 (only) | Z_1 | Λ ≠ I (always) | detected | no |
| M12 | W ← 2Δ_j + W·(I − W)⁻¹, computed as `solve((I−W).T, W.T).T` (plain transpose; also in Z) | 2 | none: **algebraically equivalent** (see below) | never; negative control | not detected | — |
| M13 | G11 ← G11 + G1j·P·Gj1 (sign) | 3 | Ĝ^(2)_11 | Ĝ_{1,j−1}P_jĜ_{j−1,1} ≠ 0 (always) | detected | no |
| M14 | G11 init ← solve(I − W_1, I) (sign) | 2 | Ĝ^(1)_11 → endpoint blocks | always; at N = 2 it flips the sign of Ĝ_N1 and Ĝ_1N, so **\|t\| and transmission probabilities are unchanged**. Only amplitudes and r detect it | detected | no |
| M15 | M6 + M7 | 2 | Ĝ_N1, Ĝ_1N | κ_Z ≠ 0 | invisible | **yes** (keeps Ĝ_1N = Ĝ_N1ᵀ) |
| M16 | closing G11 from the NEW endpoint blocks: G11 + G1N·Zi·GN1 | 2 | Ĝ_11 | Z_N⁻² ≠ I (always here) | detected | no |

Verifications:
- **M12 is equivalent.** W and (I − W)⁻¹ commute: from (I − W)W = W(I − W), multiply on the left and on the right by
  (I − W)⁻¹ to get W(I − W)⁻¹ = (I − W)⁻¹W, provided I − W is invertible (cond₂ ≤ 4.45 here). The double results differ
  only by rounding, so D_mut ≤ τ_end is predicted, generally not exactly 0.
- **M5 and M8 algebra:**
  - The correct terms Ĝ_{1,j−1}P_jĜ_{j−1,1} and Ĝ_{1,N−1}Z⁻¹Ĝ_{N−1,1} follow from the Schur update of the (1,1) block:
    row 1 → slice j−1 → slice j → slice j−1 → column 1.
  - The swapped versions are products of the same N_y × N_y factors in an index-inconsistent order. They coincide with the
    correct terms whenever the two chain blocks are equal, which is exactly the case at the first opportunity (M5 at j = 2,
    M8 at N = 2), and whenever the chain blocks commute with the middle factor (class D).
  - In class R both the correct and the swapped terms are symmetric matrices, so neither M5 nor M8 breaks reciprocity.
- **Class D cannot expose pure order mutants** (M1, M2, M3, M5, M6, M7, M8, M15).
  - All blocks are diagonal, so in exact arithmetic every product commutes. The preflight shows every condition is exactly
    0.
  - In double, a product of diagonal matrices has a single nonzero term per entry, plus exact zeros. Complex
    multiplication is commutative in IEEE arithmetic, so bitwise equality is expected.
  - This relies on BLAS kernels not changing the rounding of a single-term product with operand order. Predicted:
    D_mut ≤ 1e-15, most likely exactly 0.
- **Reciprocity is not a sufficient ordering test.** M3, M5, M8 and M15 preserve S = Sᵀ in class R. They are detectable
  only through comparison with C, through the intermediate audit, and (M3/M15) through class H, where |Q − Q_ref| is
  judged. M14 is invisible to |t|.

Mutant coverage by N:

| N | mutants on the code path |
|---|---|
| 1 | M11 |
| 2 | M6, M7, M10, M12, M14, M15, M16 (M8 on the path but identical) |
| 3 | all except M5 (on the path but identical) and M11 |
| 4 | all except M11 |

## 7. Predictions

Algebraic (exact):
- **E1.** Reference checks (a)–(d) hold: agreement ≤ 1e-55; residual ≤ 1e-55; the N = 1 assembly is exact.
- **E2.** The mp W/Z recursion on R_W and the mp original recursion on R_orig agree with C on their lattices to ≤ 1e-50
  (supplementary).
- **E3.** The prefix identities and the closing identity hold to ≤ 1e-50. Already observed in the preflight: ≤ 3.5e-61.

Numerical (double):
- **N1.** A vs C(R_W):
  - e_glob of all four endpoint blocks between 1e-16 and 1e-14, and ≤ τ_end;
  - ΔS_max ≤ τ_S;
  - U, Q (or |Q − Q_ref| in class H) and the symmetry residual within the §5 bounds.
- **N2.** B vs C(R_orig): the same scales and bounds. No failure is expected: production was validated with disorder in B3.
- **N3.** C(R_W) − C(L) and C(R_orig) − C(L): e_glob ≤ 10·cond₂(M)·ε ≤ 6e-14.
- **N4 (revised from the preflight).** cond₂(M) is in [1.0, 24.7]; cond₂(M^(j)) ≤ 20.8; cond₂(I − W_j) ≤ 4.45;
  cond₂(Z_N) ≤ 4.26; ‖P_j‖₂ ≤ 1.97. The values recomputed at experiment time must match the preflight to ≤ 1e-6 relative
  (they come from the same definitions).
- **N5 (revised).** In main-matrix classes R and H:
  - κ_chain ≥ 2.0e-2 (N ≥ 3), κ_W ≥ 3.6e-2 (N ≥ 3), κ_Z ≥ 7.1e-2 (N ≥ 2), κ_[Δ,Λ] ≥ 9.9e-2;
  - M5 condition ≥ 3.3e-2 (N = 4);
  - M8 condition ≥ 2.8e-2 (N ≥ 3) and exactly 0 at N = 2.
  - In class D everything is exactly 0.
  - In the scaling matrix κ ∝ s, with log-log slope 1.00 ± 0.05.
  - (These values were observed in the preflight. The experiment must reproduce them, so they are consistency
    predictions, not independent ones.)
- **N6.** Intermediate audit (A's instrumented copy vs the mp prefix references): every W_j, P_j and prefix block is within
  its own §5 bound.
- **N7 (main matrix).**
  - Every mutant that is on the code path and whose condition is nonzero gives D_mut ≥ 1e-6. (Expected D_mut ~ condition
    magnitude × O(0.1–1), i.e. ≳ 1e-3; not proven.)
  - The first deviating prefix quantity is the one listed in §6.
  - Mutants that are identical by algebra give D_mut = 0 bitwise: M5 at N = 3, M8 at N = 2.
  - Pure order mutants in class D give D_mut ≤ 1e-15.
  - M12 gives D_mut ≤ τ_end.
- **N8 (scaling matrix).**
  - For ordering mutants, D_mut falls with s with log-log slope in [0.8, 2.2]. The expected slope is ≈ 1, because the
    mutant terms are linear in the commutators and those scale ∝ s. The pre-registered range is not narrowed after the
    preflight.
  - D_mut stays ≫ τ_end at every s. Expected ≳ 1e-6 even at s = 2⁻¹², where the conditions are ~3.5e-5.
  - The correct-code error stays ≤ τ_end with no monotone trend in s.
  - The detection threshold 1e-6 of N7 is not applied to the scaling cases.
- **N9.** The 4 duplicated s = 1 executions are bit-identical to their main-matrix counterparts.

## 8. Decision rules

- **Invalid** (no conclusion about W/Z) if any of the following holds:
  - a reference check (a)–(d) fails;
  - the instrumented copy is not bit-identical to `wz_endpoint_blocks`;
  - a main-matrix R/H case whose code path contains an ordering-sensitive mutant has its relevant condition below 1e-3.
    The preflight shows none does. The scaling cases are exempt.
- **FAIL** if any A-vs-C(R_W) endpoint or prefix e_glob, absolute W/P error, or ΔS_max exceeds 1e-10 in a case whose τ is
  ≤ 1e-11.
- **Anomalous** if an error lies between its own §5 bound and 1e-10. It is reported, not counted as a pass.
- **Design insufficient** (not a W/Z failure) if any of the following holds:
  - an N7 detection is missed;
  - an "identical by algebra" mutant gives D_mut ≠ 0;
  - M12 exceeds τ_end;
  - a pure order mutant is visible in class D (D_mut > 1e-12);
  - N8's slope range fails for an ordering mutant.
- **Pass** only if the experiment is valid and has no FAIL, no anomaly and no design insufficiency, and N1, N6 and N9 hold.
  - Claim, if it passes: "the W/Z recursion reproduces the exact Green blocks of these 72 small deterministic coupled
    configurations (76 executions)".
  - **This is not a validation of W/Z**: no long devices, near-threshold channels, disorder ensembles, unequal lead widths
    or B1 rows are covered.
- **B vs C(R_orig)** above 1e-10 would be reported as a production finding. No change would be made.

## 9. Remaining mathematical concerns (open before implementation)

1. **Thin margin.** The minimum κ_chain = 2.0e-2 (e2b/R/N3) is only 2× the stated lower bound. Detection of M1/M2/M3 in
   that case relies on the mutant term not being suppressed by the closing factors.
2. **D_mut magnitudes are not predicted quantitatively.** The preflight gives the exact detection conditions, not the
   endpoint discrepancies, because computing those would require running the mutants. N7's 1e-6 threshold relies on a
   ~2e4 margin below the smallest condition magnitude.
3. **Short chains.** The loop runs at most twice (N = 4), so chains of more than two P factors and repeated old-value
   accumulation beyond two steps are untested. Longer devices are excluded by scope.
4. **Unequal lead widths** cannot be represented: the code and the design use one N_y for both leads. Index-inconsistent
   swaps such as M5/M8 would be shape errors there, but are dimensionally legal here.
5. **Class H** goes beyond the documented scope of `wz_rgf.py` (real symmetric Δ) and beyond the project's physics.
   Results there test the algebra only.
6. **The τ formulas are heuristics.** Block elimination without inter-block pivoting has no simple a-priori bound. With the
   mild conditioning found (≤ 25), the fixed 1e-10 threshold is not sensitive to the exact form of τ.
7. **The production lattice R_orig** relies on reproducing production's double input formation (`E*eye − H[j]`, H from
   `device_hamiltonian_blocks`, σ/γ from `lead_modes`) by calling those read-only functions. A future change to them
   would change R_orig.
8. **The documentation gap** (no committed W/Z derivation document) remains. This document's §1 is the only written
   derivation.

## 10. Planned implementation (not yet written)

- `coupled_blocks.py` (driver).
- `wz_mutants.py` (instrumented faithful copy plus M1–M16 including M4L/M4R).
- `results_coupled_blocks/` containing `coupled_blocks_results.json`, `coupled_blocks_comparisons.csv` and
  `coupled_blocks_report.md`.
- The preflight script of Appendix A is to be committed with the implementation, unchanged.

Read-only imports: `stationary`, `smatrix_spec`, `wz_rgf`. Cost: seconds.

## Appendix A — preflight script (as run; not committed separately)

```python
"""Design preflight for PREDICTIONS_coupled_blocks.md (NOT the experiment).
Uses only the defined matrices/energies; exact quantities from direct mpmath inverses of the full and prefix matrices of the
represented W/Z lattice R_W (double Delta_fl and double lead delta taken as exact, exact retarded roots).
Does not import wz_rgf, wz_variants, stationary solvers or any mutant."""
import json
import sys

import mpmath as mp
import numpy as np

mp.mp.dps = 60
dx = 1 / 8
EP = [4.934802200544679, 19.739208802178716, 44.41321980490211]   # stationary.transverse_basis(1.0, 3).E_perp (doubles)
v0 = 4.0
C2 = dict(mu=(0.5, -0.25, 1, 0), alpha=(1, 2, -1.5, 0.5), beta=(3, -1, 2.5, 1.75), gamma=(1, 0.5, -2, 1.25))
C3 = dict(mu=(0.5, -0.25, 1, 0), a=(1, 2, -1.5, 0.5), a2=(0.25, -0.5, 0.75, 1), b12=(3, -1, 2.5, 1.75), b23=(-2, 1.5, 0.5, -1),
          b13=(1, 0, -1.25, 2), g12=(1, 0.5, -2, 1.25), g23=(-0.5, 1, 0.75, -1.5))
ENER = {"e2a": (2, 12.0), "e2b": (2, 30.0), "e3a": (3, 12.0), "e3b": (3, 30.0), "e3c": (3, 50.0)}


def V(ny, j, cls, s=1.0):
    k = j - 1
    if ny == 2:
        c = C2
        A = np.diag([1.0, -1.0]); B = np.array([[0, 1.0], [1.0, 0]]); C = np.array([[0, -1j], [1j, 0]])
        M = c["mu"][k] * np.eye(2) + c["alpha"][k] * A + 0j
        if cls in "RH":
            M = M + c["beta"][k] * B
        if cls == "H":
            M = M + c["gamma"][k] * C
    else:
        c = C3
        E = lambda a, b: np.eye(3)[:, [a]] @ np.eye(3)[[b], :]
        S = lambda a, b: E(a, b) + E(b, a)
        K = lambda a, b: 1j * (E(b, a) - E(a, b))
        M = c["mu"][k] * np.eye(3) + c["a"][k] * np.diag([1.0, 0, -1]) + c["a2"][k] * np.diag([1.0, -2, 1]) + 0j
        if cls in "RH":
            M = M + c["b12"][k] * S(0, 1) + c["b23"][k] * S(1, 2) + c["b13"][k] * S(0, 2)
        if cls == "H":
            M = M + c["g12"][k] * K(0, 1) + c["g23"][k] * K(1, 2)
    return v0 * s * M          # s and v0 dyadic: exact


def root(d):
    c = 1 - d
    return c + 1j * mp.sqrt(1 - c * c) if d > 0 else c - mp.sqrt(c * c - 1)


def tomp(A):
    return mp.matrix([[mp.mpc(complex(x).real, complex(x).imag) for x in row] for row in A])


def fro(A):
    return mp.sqrt(mp.fsum(abs(A[i, j]) ** 2 for i in range(A.rows) for j in range(A.cols)))


def comm_rel(A, B):
    return fro(A * B - B * A) / (fro(A) * fro(B))


def blk(G, a, b, ny):
    return G[a * ny:(a + 1) * ny, b * ny:(b + 1) * ny]


def cond2(A):
    a = np.array([[complex(A[i, j]) for j in range(A.cols)] for i in range(A.rows)])
    return float(np.linalg.cond(a))


def assemble(Dl, Lam, n, right_lead):
    """M/t of the first n slices: diag 2 Delta_j - 2I, off-diag I, + Lambda on slice 1, + Lambda on slice n if right_lead.
    The two lead terms are added as separate operations (for n = 1 both land on the single slice)."""
    ny = Lam.rows
    M = mp.matrix(n * ny, n * ny)
    for j in range(n):
        Bj = 2 * Dl[j] - 2 * mp.eye(ny)
        for a in range(ny):
            for b in range(ny):
                M[j * ny + a, j * ny + b] = Bj[a, b]
        if j + 1 < n:
            for a in range(ny):
                M[j * ny + a, (j + 1) * ny + a] = 1
                M[(j + 1) * ny + a, j * ny + a] = 1
    for a in range(ny):
        M[a, a] += Lam[a, a]                              # left lead on slice 1
    if right_lead:
        for a in range(ny):
            M[(n - 1) * ny + a, (n - 1) * ny + a] += Lam[a, a]   # right lead on slice n
    return M


def case(eid, cls, N, s=1.0):
    ny, E = ENER[eid]
    Ep = np.array(EP[:ny])
    I = mp.eye(ny)
    Dl_fl = [(E * np.eye(ny) - np.diag(Ep) - V(ny, j, cls, s)) * (dx * dx) for j in range(1, N + 1)]
    Dl = [tomp(D) for D in Dl_fl]
    dlead = [(E - e) * (dx * dx) for e in Ep]
    Lam = mp.diag([root(mp.mpf(float(d))) for d in dlead])
    M = assemble(Dl, Lam, N, True)
    G = M ** -1
    res_full = max(abs(x) for x in (M * G - mp.eye(N * ny)))
    out = dict(case=f"{eid}/{cls}/N{N}" + ("" if s == 1 else f"/s=2^{int(round(np.log2(s)))}"), energy=eid, cls=cls, N=N, s=s,
               n_open=int(sum(d > 0 for d in dlead)), ny=ny, cond_M=cond2(M), norm2_Ghat=float(np.linalg.norm(
                   np.array([[complex(G[i, j]) for j in range(G.cols)] for i in range(G.rows)]), 2)), mp_residual=float(res_full))
    out["kappa_DeltaLambda_max"] = float(max(comm_rel(D, Lam) if fro(D) > 0 else 0 for D in Dl))
    # prefixes j = 1..N-1 (no right lead)
    W, P, Gj1, G1j, G11, condPre, condIW = {}, {}, {}, {}, {}, {}, {}
    for j in range(1, N):
        Mj = assemble(Dl, Lam, j, False)
        Gp = Mj ** -1
        Gjj = blk(Gp, j - 1, j - 1, ny)
        W[j] = I + Gjj ** -1
        P[j] = -Gjj
        Gj1[j], G1j[j], G11[j] = blk(Gp, j - 1, 0, ny), blk(Gp, 0, j - 1, ny), blk(Gp, 0, 0, ny)
        condPre[j] = cond2(Mj)
        condIW[j] = cond2(I - W[j])
    # algebraic closing identity Z = 2 Delta_N + (Lambda - I) + (I - W_{N-1})^{-1} W_{N-1}  vs  Z = Ghat_NN^{-1}
    Zi = blk(G, N - 1, N - 1, ny)
    Z = Zi ** -1
    if N >= 2:
        Zf = 2 * Dl[N - 1] + (Lam - I) + (I - W[N - 1]) ** -1 * W[N - 1]
        out["closing_identity_err"] = float(max(abs(x) for x in (Zf - Z)))
        # prefix chain identities (left/right multiplication, old values) vs direct prefix inverses
        err = 0
        for j in range(2, N):
            err = max(err, max(abs(x) for x in (P[j] * Gj1[j - 1] - Gj1[j])), max(abs(x) for x in (G1j[j - 1] * P[j] - G1j[j])),
                      max(abs(x) for x in (G11[j - 1] - G1j[j - 1] * P[j] * Gj1[j - 1] - G11[j])))
        out["prefix_identity_err"] = float(err)
    else:
        out["N1_identity_err"] = float(max(abs(x) for x in (2 * Dl[0] + 2 * (Lam - I) - Z)))
    out["cond_Z"] = cond2(Z)
    out["cond_prefix_max"] = max(condPre.values()) if condPre else None
    out["cond_I_minus_W_max"] = max(condIW.values()) if condIW else None
    out["normP_max"] = max(cond2(P[j]) * 0 + float(np.linalg.norm(np.array([[complex(P[j][a, b]) for b in range(ny)] for a in range(ny)]), 2)) for j in P) if P else None
    kW = [float(comm_rel(W[j], W[j + 1])) for j in range(1, N - 1)]
    kC = [float(comm_rel(P[j], Gj1[j - 1])) for j in range(2, N)]
    out["kappa_W"] = kW
    out["kappa_chain"] = kC
    out["kappa_Z"] = float(comm_rel(Zi, Gj1[N - 1])) if N >= 2 else None
    # exact detection conditions for M5 (step j >= 3) and M8 (closing), relative to the correct term
    m5 = []
    for j in range(3, N):
        A = G1j[j - 1] * P[j] * Gj1[j - 1]; B = Gj1[j - 1] * P[j] * G1j[j - 1]
        m5.append(float(fro(A - B) / fro(A)))
    out["cond_M5"] = m5
    if N >= 2:
        A = G1j[N - 1] * Zi * Gj1[N - 1]; B = Gj1[N - 1] * Zi * G1j[N - 1]
        out["cond_M8"] = float(fro(A - B) / fro(A))
        out["chain_asym"] = float(fro(Gj1[N - 1] - G1j[N - 1]) / fro(Gj1[N - 1]))
    return out


rows = []
for eid in ENER:
    for cls in "RHD":
        for N in (1, 2, 3, 4):
            rows.append(case(eid, cls, N))
for eid in ("e2b", "e3b"):
    for N in (3, 4):
        for e in (-4, -8, -12):
            rows.append(case(eid, "R", N, 2.0 ** e))
json.dump(rows, open(sys.argv[1], "w"), indent=1)
print(len(rows), "unique configurations evaluated")
```
