# R0 review packet — W/Z coupling-scaling experiment v2 (implementation level)

**This is a review-preparation document, not the R0 sign-off.** It contains no PASS decision. The sign-off required by §4.3
is an independent reviewer's record in `REVIEW_coupling_scaling_v2.md` (identity, date, commit). That file does not exist and
must not be created by the implementer. **Implementation-level R0: PENDING.** Running any solver without a committed R0
sign-off is INVALID (§4.3).

Section and step numbers (§x.y) refer to `PREDICTIONS_coupling_scaling_v2.md`. Line numbers are those of the reviewed commit
(Part A) and must be re-checked by the reviewer against the commit actually reviewed.

---

## A. Scope and provenance

| item | value |
|---|---|
| frozen specification | `experiments/wz_rgf_diagnostic/PREDICTIONS_coupling_scaling_v2.md`, **Revision 3.7** (`0f502c9`) |
| implementation baseline | `baab2c5986e3bf40fc1cd0210f74b0de2bb85567` (covering, T1–T5, T6 review corrections) |
| this packet | first added in `25f48ee` (import-boundary hardening, runner docs); **updated in the main-experiment infrastructure phase** (see Part I). **The reviewer must record the exact commit SHA reviewed in Part H.** |
| `main` / branch reconciliation | `main` = `c7939df` (merge of PR #23), which contains every branch commit through `f1ecc60`; the branch was fast-forwarded to `c7939df` (no history rewrite) before the infrastructure phase |

**Source files under review**

| role | path (under `experiments/wz_rgf_diagnostic/`) |
|---|---|
| production certificate, §2.4.1 / §5.6, import boundary and exact assembly | `coupling_scaling_v2_certify.py` |
| production §2.4 item 3 covering and common Taylor radius | `coupling_scaling_v2_covering.py` |
| clean-room T6 reference (exact) and C4 numerical reference | `cs2_t6_reference/` (`qexact`, `rmat`, `rinv`, `rpsd`, `rlead`, `rlattice`, `rblocks`, `rmutants`, `rpoly`, `rcert`, `rpred`, `rnum`) |
| T6 harness, staged runner, mutation suite, import guard | `cs2_t6_harness/` |
| T1–T5 self-tests and their independent references | `cs2_t1t5/` |
| double-precision mutants (system whose change column is audited) | `wz_mutants.py` (unchanged since `63e046d`) |
| W/Z production solver (unchanged) | `wz_rgf.py` |
| repository wrapper tests | `tests/test_cs2_t6.py` |
| §4.2 oracle for every mutant (main run) | `wz_mutants_oracle.py` (Part I) |
| §4.3 clean-room mpmath transcription, its materials and provenance | `wz_mutants_mp.py`, `wz_mutants_mp_MATERIALS.md`, `wz_mutants_mp_PROVENANCE.md` (Part I) |
| §4.7 validated references | `coupling_scaling_v2_refs.py` (Part I) |
| manifest generator (§2.3, §2.4, §2.6, §8.1, §8.3) | `coupling_scaling_v2_manifest.py` (Part I) |
| §5 analysis procedures; driver (gated, not run) | `coupling_scaling_v2_analysis.py`, `coupling_scaling_v2.py` (Part I) |
| infrastructure tests | `tests/test_wz_mutants_oracle.py`, `tests/test_wz_mutants_mp.py`, `tests/test_cs2_refs.py`, `tests/test_cs2_manifest.py`, `tests/test_cs2_driver.py` |

**Existing design-level R0.** The independent design review of `f101224` (verdict **CONDITIONAL PASS**; corrections E1–E17
implemented in Revision 3) verified R0 items (a) and (b) **against the document**. Per §4.3 it does **not** replace the
implementation-level sign-off of items (c)–(e).

**What exists and what does not (important for items (c) and (d)).**
- Implemented: the production certificate (§2.4.1, §5.6), the covering (§2.4 item 3), T1–T6, the T6 clean-room reference,
  the C4 numerical cross-check, the mutation suite MU1–MU27.
- **Implemented in the infrastructure phase (Part I):** `wz_mutants_oracle.py` (the §4.2 oracle for *all* mutants),
  `wz_mutants_mp.py` (the §4.3 clean-room mpmath transcription), `coupling_scaling_v2_refs.py` (§4.7 validated references),
  `coupling_scaling_v2_manifest.py`, `coupling_scaling_v2_analysis.py` and the driver `coupling_scaling_v2.py`.
  - The driver is **gated and has not been run**: it refuses without the R0 sign-off, and while the §4.7 4(b)
    "duplicate set" is undefined (gap G-1).
- Parts C and D below describe the *T6-scope* oracle and the T6 transcription. Part I covers the main-run oracle and the
  clean-room transcription for R0 items (c) and (d).

**Authorship (to be recorded per §4.3).** All implementation code listed above — production, reference, harness, T1–T5, the
oracle, the §4.7 references, the manifest generator, the analysis and the driver — was written by the same AI coding agent
(Claude Code). The exception is `wz_mutants_mp.py`, which was written by a separately started AI subagent (same model) from a
materials file only (Part I.3). The T6 reference was written before the production module
(clean-room *ordering*), but **no part was written by an independent human implementer.** The reviewer should weigh this when
assessing shared-error risk (Parts C, D, F).

---

## B. Mutant derivation audit (R0 item (a) and (b))

### B.0 Correct-code semantics and the prefix identities

Notation (§N, §4.2): Ĝ = t·G = (M/t)⁻¹; Ĝ^(j) = (M^(j)/t)⁻¹ is the j-slice prefix inverse (left lead only). M/t has diagonal
blocks 2Δ_j − 2I (+Λ on slice 1 and on slice N; 2Λ when N = 1) and off-diagonal blocks I. Lm := Λ − I.

Correct code (`wz_mutants.py`, `mutant=None`; identical to `wz_rgf.wz_endpoint_blocks`):

| step | code (line) | value |
|---|---|---|
| W₁ | `W = 2.0 * Delta[0] + Lm` (L62) | W₁ = 2Δ₁ + Λ − I, and M^(1)/t = 2Δ₁ − 2I + Λ = W₁ − I |
| init | `G11 = solve(W - I, I)` (L64) | Ĝ^(1)_11 = (W₁ − I)⁻¹; P₁ := (I − W₁)⁻¹ = −Ĝ^(1)_11 |
| chains | `Gj1 = G1j = G11` (L66–67) | χ₁ = ψ₁ = Ĝ^(1)_11 |
| loop s = 2..N−1 | `W = 2Δ_s + solve(I - W, W)` then `P = solve(I - W, I)` (L74–75) | W_s = 2Δ_s + (I − W_{s−1})⁻¹W_{s−1}; P_s = (I − W_s)⁻¹ |
| | `G11 = G11 - G1j @ P @ Gj1` (L96) | Ĝ^(s)_11 = Ĝ^(s−1)_11 − ψ_{s−1}P_sχ_{s−1} (old chain values) |
| | `Gj1 = P @ Gj1`; `G1j = G1j @ P` (L97–98) | χ_s = P_sχ_{s−1}; ψ_s = ψ_{s−1}P_s |
| closing | `Z = 2Δ_N + Lm + solve(I - W, W)`; `Zi = solve(Z, I)` (L105–107) | Z = 2Δ_N + Λ − I + (I − W_{N−1})⁻¹W_{N−1}; Ĝ_NN = Z⁻¹ |
| | `G11 = G11 + G1j @ Zi @ Gj1` (L114) | Ĝ_11 = Ĝ^(N−1)_11 + ψ_{N−1}Z⁻¹χ_{N−1} |
| | `GN1 = -(Zi @ Gj1)`; `G1N = -(G1j @ Zi)` (L115–116) | Ĝ_N1 = −Z⁻¹χ_{N−1}; Ĝ_1N = −ψ_{N−1}Z⁻¹ |

**Prefix identities used below** (Schur complement of the last slice; verified exactly in the T6 reference, see D): for
j ≥ 2, with P_j := (I − W_j)⁻¹ = −Ĝ^(j)_jj:
(I1) χ_j := Ĝ^(j)_j1 = P_jχ_{j−1}; (I2) ψ_j := Ĝ^(j)_1j = ψ_{j−1}P_j; (I3) Ĝ^(j)_11 = Ĝ^(j−1)_11 − ψ_{j−1}P_jχ_{j−1}; (I4)
W_j = I + SC_j with SC_j = A_jj/t − Ĝ^(j−1)_{j−1,j−1}, i.e. I − W_j = −(Ĝ^(j)_jj)⁻¹. Hence χ_{N−1} = P_{N−1}⋯P₂Ĝ^(1)_11 and
ψ_{N−1} = Ĝ^(1)_11P₂⋯P_{N−1} (indices decreasing, resp. increasing, left to right). θ_j := ψ_{j−1}P_jχ_{j−1}. For N = 2 the
loop is empty and χ₁ = ψ₁ = Ĝ^(1)_11.

### B.1 Per-mutant derivations

For each mutant: the code change (line in `wz_mutants.py`), the designated block X_m, the derivation, and pitfalls.

**M1** — L97 `Gj1 = Gj1 @ P` (right instead of left). X = Ĝ_N1. On path N ≥ 3.
The chain becomes χ′_s = χ′_{s−1}P_s, so χ′_{N−1} = Ĝ^(1)_11P₂P₃⋯P_{N−1}. W, P and Z are untouched (they do not depend on the
chains). GN1 = −Z⁻¹χ′_{N−1} = **−Z⁻¹Ĝ^(1)_11P₂⋯P_{N−1}** (indices increasing). The G11 update also uses χ′ but G11 is not the
designated block. *Pitfalls:* the product order must be increasing; Z⁻¹ multiplies on the left (it is not the M6 change).

**M2** — L98 `G1j = P @ G1j`. X = Ĝ_1N. N ≥ 3. ψ′_s = P_sψ′_{s−1} ⇒ ψ′_{N−1} = P_{N−1}⋯P₂Ĝ^(1)_11; G1N = −ψ′_{N−1}Z⁻¹ =
**−(P_{N−1}⋯P₂Ĝ^(1)_11)Z⁻¹**. *Pitfall:* mirror of M1 — decreasing indices, Z⁻¹ on the right.

**M3** — L97 and L98 both mutated. X = Ĝ_N1 depends only on the Gj1 chain ⇒ **as M1**. *Pitfall:* G1N is also changed but
is not designated.

**M4** — L78–81: both chains updated *before* `G11 = G11 - G1j @ P @ Gj1`. X = Ĝ_11. N ≥ 3.
In iteration j the subtracted term is (ψ_{j−1}P_j)P_j(P_jχ_{j−1}) = ψ_{j−1}P_j³χ_{j−1} instead of ψ_{j−1}P_jχ_{j−1}. The chains
themselves are correct, so the closing term ψ_{N−1}Z⁻¹χ_{N−1} is unchanged, and G11 errors accumulate additively:
**Ĝ_11 − Σ_j ψ_{j−1}(P_j³ − P_j)χ_{j−1}**. *Pitfall:* P³ — rounding amplification; this is a structural (not order) mutant.

**M4L** — L82–85: Gj1 updated before, G1j after. Subtracted term ψ_{j−1}P_j(P_jχ_{j−1}) = ψ_{j−1}P_j²χ_{j−1} ⇒
**Ĝ_11 − Σ_j ψ_{j−1}(P_j² − P_j)χ_{j−1}**.

**M4R** — L86–89: G1j before, Gj1 after. Subtracted (ψ_{j−1}P_j)P_jχ_{j−1} = ψ_{j−1}P_j²χ_{j−1} ⇒ **as M4L**.

**M5** — L91–92 `G11 = G11 - Gj1 @ P @ G1j`. X = Ĝ_11. N ≥ 3 (nonzero from N = 4).
Subtracted χ_{j−1}P_jψ_{j−1} instead of θ_j ⇒ **Ĝ_11 − Σ_j (χ_{j−1}P_jψ_{j−1} − θ_j)**. For j = 2: χ₁ = ψ₁ = Ĝ^(1)_11, so
χ₁P₂ψ₁ = ψ₁P₂χ₁ and the term vanishes. *Pitfall:* index-inconsistent but conformable because both leads have width N_y
(§1 of `PREDICTIONS_coupled_blocks.md`, structure (iii)).

**M6** — L115 `GN1 = -(Gj1 @ Zi)`. X = Ĝ_N1. N ≥ 2. **−χ_{N−1}Z⁻¹** (factor order swapped).

**M7** — L116 `G1N = -(Zi @ G1j)`. X = Ĝ_1N. N ≥ 2. **−Z⁻¹ψ_{N−1}**.

**M8** — L114 `G11 = G11 + Gj1 @ Zi @ G1j` (closing). X = Ĝ_11. N ≥ 2. Added χZ⁻¹ψ instead of ψZ⁻¹χ ⇒
**Ĝ_11 + χ_{N−1}Z⁻¹ψ_{N−1} − ψ_{N−1}Z⁻¹χ_{N−1}**; zero change at N = 2 (χ₁ = ψ₁) — "not applicable" there.

**M9** — L70–72: `P = solve(I - W, I)` evaluated *before* the W update. X = Ĝ_N1. N ≥ 3.
In iteration s the stale P equals (I − W_{s−1})⁻¹ = P_{s−1} (P₁ at s = 2). The chain is χ′_s = P_{s−1}χ′_{s−1} ⇒
χ′_{N−1} = P_{N−2}⋯P₁Ĝ^(1)_11, while W and Z are updated correctly. GN1 = **−Z⁻¹(P_{N−2}⋯P₂P₁)Ĝ^(1)_11** (indices decreasing);
N = 3: −Z⁻¹P₁Ĝ^(1)_11. *Pitfall:* P₁ = −Ĝ^(1)_11 (sign).

**M10** — L102–103: `Z = 2Δ_N + solve(I - W, W)` (Lm dropped). X = Ĝ_NN. N ≥ 2. Z′ = Z − (Λ − I) ⇒ **(Z − (Λ − I))⁻¹**.
Unsupported by the §5.6 certificate (T6.8 P6).

**M11** — L53 (N = 1 branch) `Z = 2Δ₁ + Lm` instead of 2Δ₁ + 2Lm. Off-path for every N ≥ 2 (`on_path`, L26–29).

**M12** — L34–35 `_wpart`: `solve((I−W).T, W.T).T` = ((I − W)ᵀ)⁻¹Wᵀ)ᵀ = W(I − W)⁻¹, used in the loop and in Z. X = Ĝ_11.
W and (I − W)⁻¹ commute exactly (rational functions of the same matrix), so W(I − W)⁻¹ = (I − W)⁻¹W: **Ĝ_11 (algebraic
identity)**. *Pitfall:* plain transpose `.T`, **not** the conjugate transpose; with a conjugate transpose the identity
would fail for complex W.

**M13** — L93–94 `G11 = G11 + G1j @ P @ Gj1` (sign). X = Ĝ_11. N ≥ 3. Each iteration changes −θ_j into +θ_j ⇒
**Ĝ_11 + 2Σ_jθ_j**.

**M14** — L64 init `G11 = solve(I - W, I)` = −Ĝ^(1)_11. X = Ĝ_11. N ≥ 2.
G11, Gj1, G1j all start negated: χ′ = −χ, ψ′ = −ψ throughout. Loop: G11′ ← G11′ − ψ′Pχ′ = G11′ − θ_j, so after the loop
G11′ = −Ĝ^(1)_11 − Σθ_j = Ĝ^(N−1)_11 − 2Ĝ^(1)_11. Closing: + ψ′Z⁻¹χ′ = +ψZ⁻¹χ ⇒ **Ĝ_11 − 2Ĝ^(1)_11**; also Ĝ_N1 → −Ĝ_N1 and
Ĝ_1N → −Ĝ_1N.

**M15** — M6 + M7 (L115–116). X = Ĝ_N1 ⇒ **as M6**.

**M16** — L109–112: GN1 = −Zi Gj1 and G1N = −G1j Zi computed first, then `G11 = G11 + G1N @ Zi @ GN1`. X = Ĝ_11. N ≥ 2.
G11 = Ĝ^(N−1)_11 + (−ψZ⁻¹)Z⁻¹(−Z⁻¹χ) = **Ĝ^(N−1)_11 + ψ_{N−1}Z⁻³χ_{N−1}**. *Pitfall:* Z⁻³.

**Reviewer task (item (a)):** confirm each "code change" entry above against `wz_mutants.py` line by line, and each oracle
by an independent derivation (item (b)). The §4.3 residual risk — an error in the change column shared by all
implementations — is detectable only here.

---

## C. Oracle implementation audit (R0 item (c))

**Scope available now: the 8 pure order mutants** (M1, M2, M3, M5, M6, M7, M8, M15), as the §5.6 polynomial forms
F_m := (oracle value of X_m) − (designated reference block). The full §4.2 oracle for M4, M4L, M4R, M9, M10, M12, M13, M14, M16
on the experiment lattices is implemented in `wz_mutants_oracle.py` (Part I.2).

| mutant | §5.6 form | reference (`cs2_t6_reference/rpoly.py`, `monomials`) | production (`coupling_scaling_v2_certify.py`, `_monomials`, L232–259) |
|---|---|---|---|
| M1, M3 | −Ĝ_NN·Ĝ^(1)_11·P₂⋯P_{N−1} − Ĝ_N1 | `["G_NN", "G1_11"] + [("P", j) for j in range(2, N)]`, `["G_N1"]` | same keys, `# M1/M3 longest` / `single` |
| M2 | −P_{N−1}⋯P₂·Ĝ^(1)_11·Ĝ_NN − Ĝ_1N | `[("P", j) for j in range(N-1, 1, -1)] + ["G1_11", "G_NN"]`, `["G_1N"]` | same |
| M5 | −Σ_{j=2}^{N−1}(χ_{j−1}P_jψ_{j−1} − ψ_{j−1}P_jχ_{j−1}) | loop j = 2..N−1, two monomials each | same (`# M5 terms`) |
| M6, M15 | −χ_{N−1}Ĝ_NN − Ĝ_N1 | `[("chi", N-1), "G_NN"]`, `["G_N1"]` | same |
| M7 | −Ĝ_NNψ_{N−1} − Ĝ_1N | `["G_NN", ("psi", N-1)]`, `["G_1N"]` | same |
| M8 | χ_{N−1}Ĝ_NNψ_{N−1} − ψ_{N−1}Ĝ_NNχ_{N−1} | two length-3 monomials | same |

- **Multiplication order:** monomial key lists are multiplied left to right (`rpoly.evaluate`; production `F_poly`).
- **Prefix indexing:** blocks are 1-indexed slices; block (a, b) = rows (a−1)N_y..aN_y−1, columns (b−1)N_y..bN_y−1
  (`rmat.block`, production `_blk`). P_j = −(t·Y_j)_jj, χ_j = (t·Y_j)_j1, ψ_j = (t·Y_j)_1j from prefix matrix j; G1_11 from
  prefix 1; G_NN, G_N1, G_1N from the full matrix.
- **Boundary cases:** `on_path`: M1/M2/M3/M5 need N ≥ 3; M6/M7/M8/M15 N ≥ 2; M8 at N = 2 "not applicable".
  N = 1: no order mutant is on path (§2.2).
- **Conjugation/transpose:** none in the polynomial forms; all products are plain matrix products. Hermiticity is not assumed
  (class H lattices X2 exercise non-symmetric blocks).
- **Shared assumptions oracle ↔ transcription:** both production and reference take the forms from §5.6 (same source).
  They are checked against the *recursion* (`rmutants`, Part D) exactly on every T6.1 input (304 cases, stage 2 self-proof
  and T6.1), which is a different source (code-change column). A §5.6 form error would therefore be caught **unless** the
  transcription inherits the same misreading (Part D).

---

## D. Clean-room transcription audit (R0 item (d))

- **Implemented:** `cs2_t6_reference/rmutants.py`, `endpoint_blocks(Delta, lam_minus_1, mutant)`: exact-rational transcription
  of `wz_mutants.wz_instrumented` (inputs: the N_y × N_y blocks Δ_j and λ_n − 1 as exact rationals; outputs: G11, GN1, G1N,
  GNN in Ĝ units). `solve(X, B)` is X⁻¹B with an exact inverse proved by A·X = X·A = I (`rinv`). The docstring carries a
  line-correspondence table to `wz_mutants.py`.
- **Independence enforcement:** the AST import guard (stage 1) forbids any import of `wz_mutants` (and of production, harness,
  numpy, other repository modules) from the reference, and forbids dynamic imports (`__import__`, `import_module`, `exec`,
  `eval`); the guard is self-tested on synthetic sources.
- **Not clean-room in the §4.3 sense.** By the frozen T6 design (§2.4.3.4) this transcription was written **from the source
  text of `wz_mutants.py`** (read, never imported). §4.3's clean-room transcription for the main run must instead be written
  from §4.2's change column and the correct-code listing *without reading* `wz_mutants.py`, or by an independent implementer
  — that artefact is `wz_mutants_mp.py` (Part I.3; independence only partial).
- **Formula source:** the transcription follows the *code*; the oracle follows the *§5.6/§4.2 formulas*. They meet in the
  T6.1 equality (304 exact cases) and in the self-proof `rmutants(None) == rblocks` (direct inverses; every T6 lattice).
- **How R0 can detect a shared specification error:** (i) item (a) — re-derive each change-column entry from `wz_mutants.py`
  (Part B); (ii) item (b) — derive each oracle independently from the code semantics (B.0) rather than from §4.2; (iii) compare
  against the double mutants numerically (`cs2_t6_harness/anchors.py` does this only for the *correct* code, diagnostically).

---

## E. Certificate audit (R0 item (e)) — expression → code

All functions below are in `coupling_scaling_v2_certify.py` unless stated. Types: Q = `int`/`Fraction` (never `bool`);
QC = (Q, Q); QMat = tuple of row tuples of QC. Every non-exact input raises `TypeError` (`_q`, L35; `_qc`, L41; `_mat`, L47).

| expression (spec) | function (line) | exact? | failure conditions | tests |
|---|---|---|---|---|
| ū(f) = (isqrt(⌈f·4^ℓ⌉)+1)/2^ℓ, ū(0) = 0; l̄(f) = isqrt(⌊f·4^ℓ⌋)/2^ℓ (§2.4.1 step 4), ℓ = 256 | `_ub` L129, `_lb` L138 | yes (integers) | f < 0 → ValueError | T4 (ū KATs), T1–T6 indirectly |
| lead enclosure lo, hi, μ̃, e, λ̃ (§2.4.1) | `lead_enclosure` L146 (lo/hi L156–157, λ̃ L160) | yes | γ ≤ −1 or \|γ\| = 1 → ValueError | T3, T6.4, C4 Q0/Q1, MU6 |
| A_mid (W/Z Δ form) | `assemble_A_mid` L165 | yes | bad kind/prefix → ValueError | T6 (A_mid equality vs reference H form), T5 (R_W) |
| A_mid, A1 (L lattice, H form) | `assemble_L` L509, `assemble_A1` L544 | yes | float input → TypeError | T5, T1(f), T4 |
| R_W, R_orig assembly | `assemble_R_W` L564, `assemble_R_orig` L570 | yes | float input → TypeError | T5 |
| binary64 import | `import_binary64` L441, `import_complex128` L467, `import_matrix` L488 | yes (Fraction of a binary64) | non-binary64 format, NaN, ±inf | T5 |
| ε_A per matrix type (§2.4.1 lead perturbation) | `eps_A` L201 (prefix L207, N = 1 L209, N ≥ 2 L210); `eps_A_L` L558 | yes | — | T6.4 exact + targeted (G1–G3), MU1–MU5 |
| Ξ = I − A_mid·Y; f_Ξ, f_Y, f_mid (step 2–3) | `certify_point` L214 | yes | — | T1, T6.3, T6.4, C4 Q0 |
| ω = ū(f_Ξ) + ε_A·ū(f_Y) (step 5) | `certify_point` L220 | yes | ω ≥ 1 → "uncertified", a = None (L222) | T1(d,e), T6.3, T6.8 P1/P2/P4, MU14 |
| a = ū(f_Y)/(1 − ω) (step 6) | `certify_point` L226 | yes | only if ω < 1 | T1, T6.3, MU12 |
| h = ū(f_mid) + ε_A (step 7) | `certify_point` | yes | — | T1, C4 Q3 |
| b = ū(‖A1‖_F²) (step 7) | `b_bound` L598 | yes | — | T1(f), T4 |
| Y at dps 60, retry 120 (steps 1, 8) | `approx_inverse_Y` L398, `certify_with_retry` L603; exact mpf → Fraction `_mpf_to_fraction` L415 | Y is the exact dyadic value of the mpmath result | mp inverse failure → uncertified; NaN/±inf → ValueError | T1, T2, T6.7 (conversion), C4 Y_prod |
| r_κ = 1/(2a_κb), b = 0 ⇒ r = +∞ (never multiplied) (§2.4 item 3) | `coupling_scaling_v2_covering.py`: `half_radius` L38, `r_times_b` L44, `cover` L81 | yes | b = 0 → B0EvaluationError if evaluated | T2, T4 |
| coverage of [0, 1], bisection depth 30, 2000 points | `coverage_gaps` L62, `cover` L81 (constants L19–22) | yes | uncovered → configuration rejected | T2, T4 |
| r₀ = min_A r_0,A = 1/(2β) | `taylor_radius` L154, `common_taylor_radius` L164, `certify_configuration` L177 | yes | s = 0 uncertified → None (filter fails) | T2(f), T4 |
| η_A = t·ū(f_Y)·ω/(1 − ω) (§5.6) | `_eta_from_cert` L275, `eta_A_bound` L290 | yes | uncertified → ValueError | T6.2, C4 Q6 |
| G_lo = max(t/h_W, t·l̄(max\|Y_ab\|²) − η_full) (§5.6) | `_g_lo_from` L281 (L286), `G_lo_bound` L298 | yes | never clamped | T6.6, C4 Q7, MU13, MU27 |
| F_m polynomial (§5.6) | `_monomials` L232, `F_poly` L262 | yes | unsupported → ValueError | T6.1, MU7–MU9 |
| β̃_i = ū(‖B̃_i‖_F²); E_m = Σ\|c\|[Π(β̃ + η) − Πβ̃] | `pointwise_certificate` L314 (E_m L353) | yes | — | T6.2, MU10 |
| U_m = ū(‖F_m(B̃)‖_F²) + ū(N_y)·E_m | `pointwise_certificate` L355 | yes | U_m < 0 → INCONCLUSIVE (fail-closed) | T6.5, MU11 |
| U < THRESHOLD·G_lo, THRESHOLD = `Fraction(1, 10**10)` (L25) | `point_infeasible` L366 | yes | non-exact type; threshold ≠ literal; U < 0; G_lo < 0 | T6.6 iface a–i, T6.7, MU15–MU17, MU24–MU26 |
| every-point conclusion over exactly S_test (43 points, L30–31) | `conclude` L380 | yes | set ≠ S_test → INCONCLUSIVE | T6.8 P1–P7, MU19–MU21 |
| required prefixes; missing/uncertified/unsupported handling | `_required_prefixes` L306, `pointwise_certificate` | — | INCONCLUSIVE with reason | T6.8 P3, P6, MU22–MU23 |

---

## F. Validated-residual audit

Derivation (§2.4.1 steps 5–7; §5.6), with the norm used at each step:
1. A = A_mid + A_e with A_e diagonal and s-independent. |A_e,nn| = t·|λ_n − λ̃_n| = t·|√ξ_n − μ̃_n| ≤ t·e_n (both branches,
   since λ̃ differs from λ only through √ξ → μ̃). **Frobenius**: ‖A_e‖_F² ≤ Σ over lead entries of (c·t·e_n)².
   - full N ≥ 2: two entries per channel (slices 1 and N), c = 1 ⇒ ε_A = ū(2Σ(te_n)²);
   - full N = 1: one block receives both leads, c = 2 ⇒ ε_A = ū(Σ(2te_n)²) (*not* ū(2Σ(te_n)²));
   - prefix: slice 1 only ⇒ ε_A = ū(Σ(te_n)²).
   (`eps_A` L201–210; T6.4 targeted lattices prove sensitivity: G1 detects the N = 1 factor, G2/G3 the N ≥ 2 and prefix
   factors; MU1–MU5.)
2. R := I − A·Y = Ξ − A_e·Y. **2-norm bounded by Frobenius**: ‖R‖₂ ≤ ‖R‖_F ≤ ‖Ξ‖_F + ‖A_e‖_F·‖Y‖₂ ≤ ‖Ξ‖_F + ε_A‖Y‖_F
   ≤ ū(f_Ξ) + ε_A·ū(f_Y) = ω. Frobenius is only ever used as an *upper* bound of the 2-norm; there is no interchange in the
   other direction.
3. ω < 1 ⇒ A·Y = I − R invertible-Neumann; A⁻¹ = Y(I − R)⁻¹ ⇒ ‖A⁻¹‖₂ ≤ ‖Y‖₂/(1 − ω) ≤ ū(f_Y)/(1 − ω) = a.
4. h = ū(f_mid) + ε_A ≥ ‖A_mid‖_F + ‖A_e‖_F ≥ ‖A‖_F ≥ ‖A‖₂.
5. A⁻¹ − Y = Y(I − R)⁻¹R ⇒ ‖A⁻¹ − Y‖_F ≤ ‖Y‖_F·ω/(1 − ω), hence η_A = t·ū(f_Y)·ω/(1 − ω) bounds t‖A⁻¹ − Y‖ in **Frobenius**
   (stronger than the 2-norm statement of §5.6). Consequently E_m already bounds ‖F_m(B) − F_m(B̃)‖_F; the ū(N_y) factor in
   U_m is redundant but valid (documented in Revision 3.5).
6. G_lo: t/h_W ≤ t/‖A‖₂ ≤ t‖A⁻¹‖₂ = ‖Ĝ‖₂; and |Ĝ_ab| ≤ ‖Ĝ‖₂ with |tY_ab − Ĝ_ab| ≤ η, l̄ ≤ √ ⇒ t·l̄(max|Y_ab|²) − η ≤ ‖Ĝ‖₂.
   For every frozen input G_lo ≥ t/h_W > 0 (f_mid > 0 because A_mid has off-diagonal blocks t·I for N ≥ 2).

**Prefix vs full:** each prefix matrix M^(j) gets its own A_mid, ε_A (left lead only), Y, ω, a, h, η; its blocks feed P_j,
χ_j, ψ_j. **Exact vs approximate:** everything is exact rational except Y, which is the *exact* dyadic value of the mpmath
result — its accuracy affects only whether ω < 1, never the validity of an inequality. **Retry:** dps 60, then dps 120 if
ω ≥ 1 or if the mpmath inverse fails; still ω ≥ 1 ⇒ uncertified, no bound (`certify_with_retry`).

**Not yet implemented for the main run:** the §4.7 item 4(a) acceptance of a validated reference (ω/(1 − ω) ≤ 1e-55 with
escalation and INVALID), and the dps-100 / `lu_solve` cross-checks (driver not written).

---

## G. T1–T6 coverage matrix

Results at the reviewed commit come from `python -m cs2_t6_harness.run --out <dir>` (official mode); the reviewer should
re-run it. "Exact" = decided by exact rational predicates; C4 is a numerical cross-check only.

| test | mathematical claim | test code | independent reference | expected outcome | does NOT establish |
|---|---|---|---|---|---|
| T1 | ω < 1 ⇒ a ≥ ‖A⁻¹‖₂; h ≥ ‖A‖₂; b ≥ ‖A1‖₂; ω ≥ 1 ⇒ no bound | `cs2_t1t5/t1.py` | closed forms; exact 2×2 σ_max formula; exact PSD; mpmath dps-300 SVD with margin 10⁻²⁹⁰ for (f) | all bounds hold; (d) A⁻¹/2 for (b) and (e) Y = 0 uncertified | (f) is numerical-reference based; A1 := A in (a)–(c) is a test-file choice |
| T2 | no certified interval/disk contains a singular point; singular families rejected | `cs2_t1t5/t2.py` | exact singular points; squared distances for 1/2 ± i2⁻¹¹ | (a)–(e),(g) not covered and rejected; (f) disks exclude both points, r₀² < 1/4 + 2⁻²² | correctness of the covering for matrices far from the tested families |
| T3 | lo² ≤ ξ ≤ hi²; width; branch; e > 0; λ̃ ≠ λ | `cs2_t1t5/t3.py` | exact inequalities (mpmath dps 300 supplementary only) | all hold for the 7 deltas and the fixed N_y = 4 deltas | enclosure correctness for other γ |
| T4 | ū(0) = 0, ū(f)² > f; b = 0 single-point certificate with r = +∞; 1/(2ab), r·b never evaluated at b = 0; prefix-only b = 0 | `cs2_t1t5/t4.py` | exact; live instrumentation (`B0EvaluationError`) | no exception; flags correct; one point | that the full-matrix covering of the T4 lattice succeeds (reported only) |
| T5 | exact binary64 import; strict type boundary; L/R_W/R_orig assembly equals an independent Fraction-only routine; no float after import; A_mid(2⁻⁴⁰) − A_mid(0) = −2⁻⁴⁰ blockdiag(V(1)) | `cs2_t1t5/t5.py` | `cs2_t1t5/indep.py`; line tracer with liveness controls | all equal / rejected as specified | correctness of the *source* arrays (stationary, Δ_fl) themselves |
| T6.1 | F_m on exact blocks = recursion difference | `cs2_t6_harness/t6_1.py` | `rmutants` − `rblocks` | 304 exact equalities | the oracle for non-order mutants |
| T6.2 | η_A ≥ t·2⁻⁴⁰; ‖F(B) − F(B̃)‖₂ ≤ E_m | `t6_2.py` | exact PSD | holds | tightness |
| T6.3 | certified, ω < 1, a ≥ ‖A⁻¹‖₂ | `t6_3.py` | exact inverse + PSD | holds (incl. Y_half aux) | — |
| T6.4 | ε_A² ≥ ‖A_e‖_F²; ω ≥ ‖I − A_true Y‖₂ (rational-λ lattices); lead-factor sensitivity (G1–G3) | `t6_4.py`; witnesses | exact | holds; MU1–MU5 detected | irrational leads (→ C4) |
| C4 | Q1–Q8 bounds on the correct side of reference intervals, irrational leads | `cs2_t6_harness/c4.py`, `cs2_t6_reference/rnum.py` | mpmath dps 300/400/500 + exact Stage 2 under unverified assumption A1 | field CONSISTENT, non-blocking | **not a validation** ("numerical cross-check (not validated)") |
| T6.5 | U_m ≥ 0, U_m² ≥ D*² | `t6_5.py` | exact recursion difference | holds (incl. Y_adv) | tightness (deferred) |
| T6.6 | G_lo > 0 on frozen inputs; G_lo ≤ ‖Ĝ‖₂; interface cases a–i | `t6_6.py` | exact PSD; counted predicate | holds | — |
| T6.7 | THRESHOLD exact and typed; comparison cases | `t6_7.py` | literals | holds | — |
| T6.8 | INCONCLUSIVE paths P1–P7; P5 value agreement and status list; edge cases | `t6_8.py`, `p5.py` | `rcert` (N1–N5, N1′) | holds | — |
| mutations | MU1–MU27 (36 variants incl. supplementary MU18-N4) each detected | `cs2_t6_harness/mutations.py` | — | 36/36 detected | detection of conservative bugs (listed as not detectable in §2.4.3.10) |

---

## H. Reviewer checklist (to be completed by the independent reviewer only)

Reviewer identity: ____________________ Date: __________ Reviewed commit SHA: ________________________________________

Mark each item PASS / FAIL / NEEDS CLARIFICATION. Nothing below is pre-filled.

| # | item | PASS | FAIL | NEEDS CLARIFICATION | comments |
|---|---|---|---|---|---|
| a.1 | §4.2 change column matches `wz_mutants.py` for every mutant (Part B) | ☐ | ☐ | ☐ | |
| b.1 | every oracle formula independently re-derived (Part B) | ☐ | ☐ | ☐ | |
| c.1 | §5.6 polynomial forms implemented correctly in production and reference (Part C) | ☐ | ☐ | ☐ | |
| c.2 | full §4.2 oracle for the main run (`wz_mutants_oracle.py`, Part I.2), incl. off-path handling and Mag_m | ☐ | ☐ | ☐ | |
| c.3 | first-deviation closed forms (Part I.2, choice C-12) | ☐ | ☐ | ☐ | |
| d.1 | T6 exact transcription `rmutants.py` correct (line correspondence) | ☐ | ☐ | ☐ | |
| d.2 | §4.3 clean-room mpmath transcription (`wz_mutants_mp.py`); provenance and partial-independence disclosure (Part I.3) | ☐ | ☐ | ☐ | |
| e.1 | §2.4.1 certificate code (Parts E, F) | ☐ | ☐ | ☐ | |
| e.2 | §2.4 item 3 covering and Taylor radius | ☐ | ☐ | ☐ | |
| e.3 | validated-residual derivation and norms (Part F) | ☐ | ☐ | ☐ | |
| e.4 | §5.6 pointwise certificate code | ☐ | ☐ | ☐ | |
| e.5 | binary64 import boundary (strict binary64 only; lower precision rejected, not promoted) | ☐ | ☐ | ☐ | |
| e.6 | T1–T5 adequacy against §2.4.2 (including the test-file input choices) | ☐ | ☐ | ☐ | |
| e.7 | T6.1–T6.8, C4 and mutation suite adequacy against §2.4.3 | ☐ | ☐ | ☐ | |
| e.8 | official run re-executed by the reviewer at the reviewed SHA (record command and exit code) | ☐ | ☐ | ☐ | |
| e.9 | §4.7 validated-reference machinery (Part I.4) | ☐ | ☐ | ☐ | |
| e.10 | manifest generator, its choices C-1–C-4 and the committed manifests/hashes (Parts I.5, I.6) | ☐ | ☐ | ☐ | |
| e.11 | driver and §5/§6 implementation; gaps G-1–G-5, choices C-5–C-16 (Part I.5) | ☐ | ☐ | ☐ | |
| x.1 | independence of implementation (single-author disclosure, Part A; transcription provenance, Part I.3) acceptable | ☐ | ☐ | ☐ | |

**Open points the reviewer should adjudicate** (implementer's disclosures, not resolutions):
1. Main-experiment modules are now written (Part I). Open specification gaps G-1–G-5 and implementer choices C-1–C-16 are
   listed in Part I.5. G-1 (undefined §4.7 4(b) duplicate set) blocks the driver gate.
2. The T6 transcription is transcribed from `wz_mutants.py` source by design (§2.4.3.4), not clean-room in the §4.3 sense.
3. Test-file choices for T1–T5 (§2.4.2 delegates fixed integers): explicit Jr/Ji arrays and dyadic E in `cs2_t1t5/inputs.py`;
   A1 := A in T1(a)–(c); the T1(f) margin read as bound ≥ σ_ref·(1 + 10⁻²⁹⁰).
4. Import boundary: float16/float32/complex64 are **rejected, not promoted**; extended longdouble/clongdouble rejected; the
   decision is made from the actual format (`finfo`), so platforms with binary64 longdouble accept it.
5. Runner: `--no-c4` can still exit 0 when all other stages pass (documented in `run.py`); only the flag-free official run
   is a complete validation.
6. MU18 (frozen constant change) is detected by T6.7 and, inside P5, by production's literal threshold check raising; the
   N4 route is demonstrated by the supplementary MU18-N4.
7. Double-precision anchors use ‖X_double − X_exact‖_F/‖Ĝ‖₂ ≤ 10·cond₂·ε as the reading of §4.7 item 2 (diagnostic only).
8. File placement: T1–T5 live in `cs2_t1t5/` rather than inside the certificate module (§8.1 wording), to keep the
   references independent of production.

Signature (reviewer): ____________________ Decision recorded in `REVIEW_coupling_scaling_v2.md`: ☐ yes

---

**Status at preparation:** B1-F FAIL; coupled-block W/Z DESIGN INSUFFICIENT; W/Z NOT VALIDATED; main coupling-scaling
experiment NOT RUN; manifest NOT GENERATED; implementation-level R0 **PENDING**.
