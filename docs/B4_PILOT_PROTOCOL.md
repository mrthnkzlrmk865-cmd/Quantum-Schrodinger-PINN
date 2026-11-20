# B4 Limited Diagnostic Pilot: Protocol

**Status:** PROTOCOL ONLY. Nothing in this document has been run.

**Classification of every pilot output:** **DIAGNOSTIC / NOT B4**.

**Primary sources:**

- [B4_IMPLEMENTATION_READINESS_AUDIT.md](B4_IMPLEMENTATION_READINESS_AUDIT.md) (the audit)
- [PHASE2_METHOD_SELECTION.md](PHASE2_METHOD_SELECTION.md) §3 B4 (the B4 protocol)
- [RGF_B1_VALIDATION.md](RGF_B1_VALIDATION.md), [RGF_B2_VALIDATION.md](RGF_B2_VALIDATION.md),
  [RGF_B3_VALIDATION.md](RGF_B3_VALIDATION.md), [RGF_B3X_VALIDATION.md](RGF_B3X_VALIDATION.md) and
  [RGF_B3Y_VALIDATION.md](RGF_B3Y_VALIDATION.md)

---

## 0. Purpose, scope and prohibitions

**Purpose.** The pilot answers the four open readiness questions in the audit, plus one confirmation (P5):

| Question | What is unknown |
|---|---|
| **P1** RGF resolution | lattice resolution dependence outside the tested range E ≈ 88–371 |
| **P2** TD stopping | an operational P_in ≤ 10⁻⁴ stopping rule |
| **P3** CAP return | CAP return under disorder-induced mode conversion |
| **P4** threshold quadrature | quadrature behaviour at the threshold-moved energy points |
| **P5** same realization | confirmation, at run time, that the same realization is used |

**The pilot must not do any of the following:**

1. Compute ⟨T⟩_E or ⟨R⟩_E, the packet-weighted RGF average over the B4 window.
2. Compute any interval distance between an RGF quantity and a TD quantity.
3. Compute any other RGF-versus-TDSE comparison for the B4 realizations.
4. Tabulate RGF and TD quantities for the same realization side by side.
5. Change any B4 threshold. The following stay **unchanged**:
   - dist(⟨T⟩_E, [T_TD, T_TD + P_in]) ≤ 10⁻³, and the same for R;
   - P_in(t_end) ≤ 10⁻⁴;
   - the quadrature-halving change ≤ 2×10⁻⁴.
6. Change the B4 protocol. Any change found to be needed goes into a **dated B4 protocol amendment** (§6), written
   after the pilot and before any B4 run.
7. Modify any existing source file, validation report or result artifact.
8. Write to any official result directory (§5).

**Reuse of pilot data in B4 is prohibited.** The pilot's TD runs use the B4 realizations at production settings. Their
T(t) and R(t) are therefore **not** to be reused as B4 TD results. The future B4 driver re-runs TD under its own
protocol; the measured cost is seconds per run.

**What a successful pilot is and is not.** It is evidence that B4 *can be run and interpreted*. **It is not a B4
result, and it does not independently validate Phase 2 transport.**

**Note on numbers in this document.** The kₓ, E, channel and threshold values below are *design-time arithmetic* from
the grid definitions in `validate_rgf_b2.py` (`kx_window`, `adjust_energy`, `e_perp`). They contain no solver output.
The pilot must reproduce them exactly and stop if it does not (gate G0, §6).

### 0.1 Pilot-internal tolerances

The tolerances below are **implementation tolerances of the pilot, not B4 thresholds**. Each is set from an existing
measured rounding level or an existing convention. None is a derived error bound; that limitation is listed as OPEN-6
in §7.3.

| Tolerance | Value | Used for | Basis in the repository |
|---|---|---|---|
| control reproduction | ≤ 10⁻¹² | C1/C2 against the committed B3/B3X/B3Y T₁ | B3X reproduced B3 with difference 0.0 |
| TD bookkeeping drift \|ρ(t) − ρ(0)\| | ≤ 10⁻¹⁰ | probability bookkeeping (§2.2) | V6 measured max \|N + A_L + A_R − N(0)\| = 1.6×10⁻¹³ |
| TD separability check (A1) | ≤ 10⁻¹⁰ | mode-n vs mode-1 clean return (§3.2) | same rounding basis as ρ |
| order floor | 10⁻⁸ | changes below this give BELOW FLOOR or NOT RESOLVABLE (§1.4) | `validate_rgf_b2.ORDER_FLOOR` (B2 convention); ≥ 300× the B3Y method spread (≤ 3.3×10⁻¹¹) |
| solver agreement | ≤ 10⁻¹⁰ | M1–M2–M3 complex S | existing B3/B3Y criterion |

The **B4 thresholds** (10⁻³, 10⁻⁴, 2×10⁻⁴) and the existing B4 error-budget items (RGF dx ≈ 10⁻⁴, CAP ≈ 10⁻⁴,
quadrature ≤ 2×10⁻⁴) are used **only as stated in PHASE2_METHOD_SELECTION.md §3**. They are never relaxed.

---

## 1. Part 1 — RGF resolution pilot

### 1.1 Fixed inputs (unchanged from B3/B3X/B3Y)

| Input | Value |
|---|---|
| Realizations | seed 0, L_d = 1, disorder region (1.9, 2.9), edge 0.02, M = 512 |
| Configurations (B3 cases) | **A:** W = 100, ζ = 0.02; **B:** W = 100, ζ = 0.05; **C:** W = 100, ζ = 0.20; **D:** W = 200, ζ = 0.05 |
| Construction | `validate_rgf_b3.case_config` and `case_potential` |
| Geometry | `validate_rgf_b3.geometry` with the dx-independent buffer 103/512 |
| Slice counts | N = 719 / 1437 / 2873 at dx = 1/512 / 1/1024 / 1/2048 (nested grids) |
| Transverse basis | N_y = 63, N_q = 127 (base) |
| Packet convention | k₀ = 20, σₓ = 0.2, σ_k = 1/(√2σₓ) = 3.5355; E = (kₓ² + π²)/2 |
| Threshold rule | `adjust_energy` (rel 10⁻³, shift 2×10⁻³) |

### 1.2 Predeclared energy points

These are the same 14 points (12 diagnostic points and 2 controls) for all four configurations. The energies, open-channel counts and thresholds are lead
properties, so they do not depend on the configuration; the lattice T₁ and R₁ do.

- Every point except the controls is an exact node of the B4 coarse grid (`kx_window(20, 0.2)`, 1274 nodes,
  h = 0.024995919209…).
- **None of the 14 is moved by `adjust_energy`.** Moved nodes are treated in Part 4.

| ID | grid i | kₓ | E | open (closed at N_y = 63) | nearest threshold n (Eₙ) | E − Eₙ (relative) | (kₓ − k₀)/σ_k | Why informative |
|---|---|---|---|---|---|---|---|---|
| L1 | 0 | 4.090097 | 13.29925 | 1 (62) | 2 (19.7392) | −6.440 (4.8×10⁻¹) | −4.50 | bottom of the window; single channel; slowest incident wave; deepest classically forbidden regions relative to E |
| L2 | 116 | 6.989624 | 29.36222 | 2 (61) | 2 (19.7392) | +9.623 (3.3×10⁻¹) | −3.68 | lower tail, 2 channels; mid-gap between thresholds 2 and 3 |
| L3 | 256 | 10.489053 | 59.94492 | 3 (60) | 3 (44.4132) | +15.532 (2.6×10⁻¹) | −2.69 | lower tail just outside the B3 range (B3 lowest: E = 88.51) |
| T3− | 191 | 8.864318 | 44.22287 | 2 (61) | 3 (44.4132) | −0.190 (4.3×10⁻³) | −3.15 | nearest unmoved node **below** the n = 3 opening (low-E cusp) |
| T3+ | 193 | 8.914310 | 44.66726 | 3 (60) | 3 (44.4132) | +0.254 (5.7×10⁻³) | −3.14 | nearest unmoved node **above** the n = 3 opening; new channel with small longitudinal k |
| T10− | 1086 | 31.235666 | 492.76821 | 9 (54) | 10 (493.4802) | −0.712 (1.4×10⁻³) | +3.18 | just below the n = 10 opening, near the 10⁻³ rule boundary |
| T10+ | 1088 | 31.285658 | 494.33099 | 10 (53) | 10 (493.4802) | +0.851 (1.7×10⁻³) | +3.19 | just above the n = 10 opening; the n = 10 channel is very slow |
| T11− | 1212 | 34.385152 | 596.10412 | 10 (53) | 11 (597.1111) | −1.007 (1.7×10⁻³) | +4.07 | just below the last opening inside the window |
| T11+ | 1214 | 34.435143 | 597.82435 | 11 (52) | 11 (597.1111) | +0.713 (1.2×10⁻³) | +4.08 | just above it; 11 open channels |
| U1 | 1061 | 30.610768 | 473.44435 | 9 (54) | 10 (493.4802) | −20.036 (4.2×10⁻²) | +3.00 | upper tail, away from thresholds; 9 channels |
| U2 | 1157 | 33.010376 | 549.77726 | 10 (53) | 11 (597.1111) | −47.334 (8.6×10⁻²) | +3.68 | upper tail; 10 channels |
| U3 | 1273 | 35.909903 | 649.69535 | **11 (52)** | 11 (597.1111) | +52.584 (8.1×10⁻²) | +4.50 | top of the window; largest number of open channels; largest k, so the largest lattice-dispersion term (k dx)²/24 |
| C1 | — (B3 energy) | 16.464466 | 140.47412 | 5 (58) | 5 (123.3701) | +17.104 (1.2×10⁻¹) | −1.00 | **control**: reproduces the B3/B3X/B3Y values (config D's largest B3 dx change) |
| C2 | — (B3 energy) | 27.071068 | 371.35616 | 8 (55) | 9 (399.7190) | −28.363 (7.6×10⁻²) | +2.00 | **control**: the B3Y W100_z0.05 control point and a failed W200_z0.05 point |

**Coverage and limits of the point set:**

- The 12 non-control points cover the 4.6 % of packet weight outside ±2σ_k
  (2(1 − Φ(2)) = 0.0455; arithmetic from P(kₓ)), the near-threshold region at both ends, and up to 11 open channels.
- They do **not** cover the n = 2 opening (packet weight there ≈ 10⁻⁵ per grid node) or thresholds n = 4–9, which lie
  inside ±2σ_k between the B3 energies.
  - Thresholds 4–9 are addressed by Part 4's local refinement windows.
- **Missing information — not guessed:**
  - The positions of narrow resonances in the four realizations are unknown, for example quasi-bound states just below
    a threshold, or Fano features. They cannot be selected safely in advance.
  - The pilot therefore does not claim to sample them. Part 4's local fine scans are the only resonance-sensitive
    element.

### 1.3 Resolutions and methods per point

**Resolutions** (all points, all configurations): dx = 1/512, 1/1024, 1/2048, with nested slices and the same buffer.

**Methods:**

| ID | Implementation | Role | Applied at |
|---|---|---|---|
| M1 | RGF, `stationary.scattering_matrix` | primary | all points, all three dx |
| M2 | Green-block sparse reference, `sparse_reference.sparse_scattering_matrix` | independent implementation of the same lattice problem | all points, all three dx |
| M3 | wave-function reference, `sparse_wavefunction.wavefunction_scattering_matrix` | independent implementation | all points, all three dx |
| M4 | extended-precision exact lattice, `sparse_exact_lattice.exact_lattice_S` | adjudication only | L1, T3+, T10+, T11+, U3 at dx = 1/2048 (all four configurations), and at any point where an M1–M2/M3 difference exceeds 10⁻¹⁰ |

M4 is used under the B3Y platform gate: x86-64 80-bit long double with the arithmetic probe passed. **If the gate
fails, M4 is "NOT EVALUABLE"; this is not a pass.**

**Transverse checks:**

- **N_y.** RGF N_y 63 → 127 at fixed N_q = 254 (the formal B3 single-variable test) at U1, U2, U3, T11− and T11+, at
  dx = 1/1024, for all configurations.
- **N_q.** RGF N_q 127 → 254 at the same points.
- These are the only tests of transverse truncation above E = 371.

**Cost** — ESTIMATE from the measured per-energy costs (B3Y II.1: M1 0.4 / 0.8 / 1.6 s, M2 3.8 / 7.6 / 15 s, M3
2.3 / 4.7 / 9.3 s at the three dx; M4 3–22 s):

| Part | Estimated cost |
|---|---|
| 14 points × 4 configurations × 3 dx × (M1 + M2 + M3) | ≈ 50 min |
| N_y / N_q checks | ≈ 5 min |
| M4 (20 cases) | ≤ 10 min |

The total is about 1 h, single-threaded BLAS (`OPENBLAS_NUM_THREADS=1`).

### 1.4 Diagnostics recorded (per configuration × point × dx × method)

- **Solver outputs:** the full complex S (stored, `.npz`), T₁ = Σₘ|(t_RL)ₘ₁|² and R₁ = Σₘ|(r_L)ₘ₁|².
- **Checks:** 1 − T₁ − R₁; U = max|S†S − I| (both products); Q = max|S − Sᵀ| and the block-wise reciprocity
  residuals.
- **Channels:** n_open and n_closed.
- **M1 conditioning:** min rcond, and the slice at which it occurs.
- **M3 / M4 extras:** the solve residual; for M4, the error indicator and the refinement-gate status.
- **Geometry and potential:** V_lead_max on both sides (the 10⁻⁶ guard is unchanged; a failure makes the point "NOT
  EVALUABLE"), the potential-sample sha256, the slice count and x_start.
- **Bookkeeping:** wall time and the method version (code fingerprint).

**Derived quantities, computed per point from the stored rows by a separate summarizer:**

- the complex-S differences M1–M2, M1–M3 and M2–M3 at each dx;
- the consecutive changes ΔT₁ and ΔR₁ (1/512 → 1/1024, 1/1024 → 1/2048), per method;
- **Observed orders, separately for Q = T₁ and Q = R₁** (per method and point):
  - Δ₁(Q) = Q(1/1024) − Q(1/512) and Δ₂(Q) = Q(1/2048) − Q(1/1024).
  - **Status of each quantity:**

    | Status | Condition |
    |---|---|
    | **BELOW FLOOR** | \|Δ₁\| < 10⁻⁸ **and** \|Δ₂\| < 10⁻⁸: the dx dependence between these three grids is below the order floor (§0.1). **Establishes:** the lattice result changes by less than 10⁻⁸ between dx = 1/512 and 1/2048. **Does not establish:** an error bound, an asymptotic regime, or continuum convergence; with three grids, a pre-asymptotic coincidence or cancellation cannot be excluded |
    | **NOT RESOLVABLE** | exactly one of \|Δ₁\|, \|Δ₂\| is below 10⁻⁸; **or** sign(Δ₁) ≠ sign(Δ₂) (non-monotone: an order is undefined) |
    | **OUTSIDE WINDOW** | both changes ≥ 10⁻⁸ with the same sign, but p_Q = log₂(\|Δ₁\|/\|Δ₂\|) ∉ [1.8, 2.2] (the existing provisional window) |
    | **RESOLVED** | both changes ≥ 10⁻⁸, the same sign, and p_Q ∈ [1.8, 2.2] |

  - p_T and p_R are computed and classified **independently**. A RESOLVED T₁ does **not** make R₁ resolved, and vice
    versa.
  - For a single incident channel, unitarity gives T₁ + R₁ = 1 up to the unitarity residual (U is recorded at
    every point; the existing criterion is U ≤ 10⁻¹⁰). So Δ(R₁) ≈ −Δ(T₁) is *expected*, but it is **checked, not assumed**. If the two statuses differ,
    or |Δ(T₁) + Δ(R₁)| exceeds the recorded unitarity residuals, the point is flagged and adjudicated as a possible
    unitarity or extraction defect.
- **Richardson ESTIMATE** (only for a RESOLVED quantity, using the finest pair):
  - error of Q(1/2048) ≈ |Δ₂|/(2^{p_Q} − 1);
  - error of Q(1/1024) ≈ 2^{p_Q}|Δ₂|/(2^{p_Q} − 1).
  - For any other status it is **NOT APPLICABLE**.
  - These are estimates on the lattice. **They are not continuum statements.**
- the method spread of Δ₁ and Δ₂ for T₁ and R₁ (H1-type check, as in B3Y).
- **Differences from the B3Y convention, disclosed.** B3Y (`convergence_rows`) computed one order from T₁ only, without
  a floor and without a sign check. The pilot adds the floor, the sign check and the separate R₁ order. The B3Y
  results are unchanged.

**Controls C1 and C2.** M1 at dx = 1/512 and 1/1024 must reproduce the committed B3/B3X/B3Y T₁ values to ≤ 10⁻¹²
(B3X reproduced B3 with difference 0.0). A mismatch stops the pilot (gate G0).

### 1.5 How the pilot separates the three questions

1. **Agreement between independent implementations.**
   - Measured by M1 vs M2 vs M3 (and M4) on the **same lattice** at the same dx.
   - Criterion: max|ΔS| ≤ 10⁻¹⁰, the existing B3/B3Y value; exceedances are adjudicated by M4 as in B3Y.
   - This establishes only that the lattice problem is solved correctly.
2. **Resolution dependence of the discretized lattice problem.**
   - Measured by the change of T₁ and R₁ with dx for each method separately.
   - If all methods show the same changes (spread ≪ the changes), the dx dependence belongs to the finite-difference
     problem itself, which is B3Y's H1 interpretation.
   - The observed order and the size of the changes are lattice properties.
3. **Evidence about continuum convergence.**
   - The repository has **no** continuum reference for the disordered 2D problem; all four methods share one
     discretization.
   - The pilot can provide only *indirect* evidence: an observed order near 2 that is stable across three grids, and a
     Richardson ESTIMATE consistent with that order.
   - **These are not a demonstration of continuum convergence and must not be reported as one.**
   - The only independent discretization in the repository is the TD solver (spectral in x). Comparing with it is B4
     itself and is prohibited here.
   - A fourth grid (dx = 1/4096) would test the stability of the order. It is **not** part of the pilot, because the
     buffer and lead-tail guard at that dx and the B1-F-type rounding growth (N ≈ 5745 slices) would first need their
     own review. It may be added by a dated pilot amendment if gate G-RGF is inconclusive.

---

### 1.6 Deterministic resolution decision (input to gate G-RGF)

**Point set per configuration.**

- The 14 points of §1.2.
- **Plus** the three remaining B3 energies (E = 88.51, 204.93, 281.90), re-run with M1 at all three dx, so that T₁ and
  R₁ are both available on three grids. C1 and C2 already cover E = 140.47 and 371.36.
- The committed T₁ of these re-runs must be reproduced to ≤ 10⁻¹² (an extension of gate G0).
- Cost: 3 × 3 × 4 M1 evaluations, i.e. seconds.

**The algorithm** is applied per configuration c and per quantity Q ∈ {T₁, R₁} separately, using M1 values once
step 0 has passed.

**Step 0 — solver agreement.** The M1–M2–M3 agreement (≤ 10⁻¹⁰, or adjudicated as in B3Y) must hold at every point and
dx. Otherwise **STOP**: a solver defect blocks B4.

**Step 1 — status.** For every point x, compute s(x, Q) ∈ {BELOW FLOOR, NOT RESOLVABLE, OUTSIDE WINDOW, RESOLVED}
(§1.4).

**Step 2 — admissibility** of each candidate dx_c ∈ {1/512, 1/1024, 1/2048} for (x, Q):

- **The Richardson estimate** associated with dx_c always uses the finest pair Δ₂ and the single observed order p:
  e(dx_c) = 2^{k·p}·|Δ₂|/(2^p − 1), with k = 2, 1, 0 for dx_c = 1/512, 1/1024, 1/2048.
  - For k = 2 this equals 2^p·|Δ₁|/(2^p − 1), because p is defined from |Δ₁|/|Δ₂|.

| s(x, Q) | dx_c = 1/512 | dx_c = 1/1024 | dx_c = 1/2048 |
|---|---|---|---|
| RESOLVED | **ADMISSIBLE** iff \|Δ₁\| ≤ 2×10⁻⁴ (change to the next finer grid) **and** e(1/512) ≤ 1×10⁻⁴ | **ADMISSIBLE** iff \|Δ₂\| ≤ 2×10⁻⁴ (change to the next finer grid) **and** e(1/1024) ≤ 1×10⁻⁴ | the next-finer change (to 1/4096) is **not available**. **ADMISSIBLE-R** (Richardson-only) iff \|Δ₂\| ≤ 2×10⁻⁴ (the change *into* 1/2048, the only one available) **and** e(1/2048) ≤ 1×10⁻⁴ |
| BELOW FLOOR | **ADMISSIBLE-F** (below floor; not an error bound) | ADMISSIBLE-F | ADMISSIBLE-F |
| NOT RESOLVABLE or OUTSIDE WINDOW | **INADMISSIBLE-U** (unclassifiable) | INADMISSIBLE-U | INADMISSIBLE-U |

- The limits are the existing B3/C11 consecutive-change value (2×10⁻⁴) and the existing B4 budget item for the RGF dx
  error (≈ 10⁻⁴), used as a ≤ 1×10⁻⁴ test. **No new value is introduced.**
- **Can dx = 1/2048 satisfy the next-finer-grid condition?** No: that would need a 1/4096 run. The condition is
  therefore replaced for 1/2048 by the change into 1/2048 together with the Richardson estimate. Admissibility of
  1/2048 rests on the observed-order assumption from three grids, and that is why it carries the label "-R".

**Step 3 — selection.** dx_sel(c) is the **coarsest** candidate that is admissible (in any of its labelled forms) for
every point in the set and for both T₁ and R₁.

- A dated B4 amendment may choose a finer dx than dx_sel, never a coarser one.
- **If no candidate is admissible:**
  - **(a)** If any (x, Q) is INADMISSIBLE-U, the finest available pair cannot classify the order. **G-RGF is NOT
    EVALUABLE** for c. The remedy is a dated **pilot** amendment adding dx = 1/4096 after the separate review of §1.5.
  - **(b)** Otherwise, every (x, Q) is classified but a limit is exceeded even at 1/2048. **G-RGF FAILS** for c. The
    options are listed in §6.

**Step 4 — what is carried into the B4 amendment:**

- dx_sel(c);
- for every (x, Q): its status, its label (ADMISSIBLE, ADMISSIBLE-R or ADMISSIBLE-F), Δ₁, Δ₂, p and e(dx_sel);
- an explicit statement whenever dx_sel = 1/2048 (no next-finer change available) or any point is ADMISSIBLE-F (no
  error bound).

The Richardson values are estimates of the dx → 0 limit of the finite-difference problem at fixed N_y and N_q, under
the asymptotic assumption. **They are not demonstrated continuum convergence.**

## 2. Part 2 — TD stopping-rule pilot (TD only)

### 2.1 Runs

| Item | Value |
|---|---|
| Realizations | the four B4 realizations (seed 0, configurations A–D), built by `Potential(cfg, disorder_from_config(cfg, 0))` |
| Settings | **production** (`waveguide_config`): L_x = 3.9, nx = 249 (dx = 3.9/250), ny = 63, dt = 2×10⁻⁴, CAP width 0.6, strength 300, power 2, x₀ = 1.2, σₓ = 0.2, k_x0 = 20, mode 1 |
| Propagation | `SplitStepDST.run` with an `observables` callback; nothing written to `results/transport/` |
| Maximum time | **t_max = 6.0** (predeclared) |
| Frame spacing | **Δt_f = 0.005** (every 25 time steps), i.e. 1201 frames; scalars only, plus wave-function snapshots at t = 1, 2, 4, 6 |

**Why t_max = 6.0.** It is not taken from the audit's exponential fits (≈ 1.3–2.9). Those fits have serious
limitations:

- **(a)** They use 0.3 time units of one realization, sampled every 0.025.
- **(b)** The late-time P_in(t) is non-monotone in three of the four realizations, so a single exponential does not
  describe it.
- **(c)** The decay can have several time scales: resonance lifetimes, diffusive tails, and slow channels near a
  threshold.
- **(d)** CAP-returned waves can re-enter the slab, which no fit of the pre-return data captures.

The choice of 6.0 instead rests on:

- **(i)** ≥ 2× the largest audit estimate, so that a trailing window of 1.0 (§2.3) still fits after a late settling;
- **(ii)** cost: measured 5.8–6.2 s per single-seed production run to t = 1.5 (B2), so ≈ 25 s per run to t = 6;
- **(iii)** the bounded dt error. V3 measured 2×10⁻⁵ at t = 1, and its growth to t = 6 is measured by run TD-R (§2.4)
  rather than assumed.

If the rule (§2.3) is not satisfied by t = 6, the outcome is **NOT SETTLED**. t_max is **not** extended ad hoc; any
extension needs a dated pilot amendment.

### 2.2 Recorded quantities (every frame)

**Regions:** a grid point belongs to a region by its x coordinate, with strict inequalities as listed. Every grid point
belongs to exactly one region; this is a unit test.

| Region | Range | Note |
|---|---|---|
| left CAP | x < 0.6 | W_cap > 0 |
| left gap | 0.6 ≤ x < 1.9 | contains the initial packet |
| slab | 1.9 ≤ x ≤ 2.9 | the nominal disorder region, as used by `transport.py` |
| right gap | 2.9 < x ≤ 3.3 | |
| right CAP | x > 3.3 | |

**Per frame:**

- P_slab(t) = **direct** Σ|ψ|²dxdy over the slab. This is the pilot's P_in.
- P_gapL, P_gapR, P_capL, P_capR, computed directly in the same way.
- A_L(t) and A_R(t): the absorbed probability from the solver's exact bookkeeping (`absorbed_left`, `absorbed_right`).
- N(t): the total norm on the grid.
- Accounting residual ρ(t) = 1 − [N(t) + A_L(t) + A_R(t)], and the initial deviation 1 − N(0). These are reported
  separately.
- T(t) = A_R + P_gapR + P_capR and R(t) = A_L + P_gapL + P_capL. These are the existing `transport.py` definitions
  (x > b and x < a).
- P_in,res(t) = 1 − T − R (the existing residual definition). By the definitions, P_in,res = P_slab + ρ exactly, since
  1 − T − R = 1 − (A_L + A_R) − (N − P_slab). The check |P_in,res − P_slab − ρ| ≤ 10⁻¹⁰ is therefore an arithmetic
  consistency test of the implementation. The physically meaningful check is the bookkeeping drift
  |ρ(t) − ρ(0)| ≤ 10⁻¹⁰, with ρ(0) = 1 − N(0) reported separately (§0.1; the same checks are used in G-TD).
- Mode-resolved region probabilities: a DST-I in y of ψ restricted to each region, giving the probability per
  transverse mode n = 1…63. Used by Part 3.

### 2.3 Trailing-window stopping criterion (does not stop at the first crossing)

Let t_k be the frames, ε = 10⁻⁴ (the B4 value, unchanged) and τ_w = 1.0.

- **Settling time:** t_s = min{t_k : P_slab(t_j) ≤ ε for **every** j ≥ k with t_j ≤ t_max}. This is the time after the
  *last* up-crossing, not the first down-crossing.
- **SETTLED** if t_s exists and t_max − t_s ≥ τ_w; otherwise **NOT SETTLED**.
- Also reported: the first crossing time t_first; the number of up-crossings after t_first; max P_slab over
  [t_first, t_max]; and the contents of the gaps and CAPs at t_s and at t_max.
- **Unabsorbed exterior content** U(t) = P_gapL + P_gapR + P_capL + P_capR, reported at t_s and t_max.
  - B4's interval criterion covers only P_in. U is booked as T or R, but part of it can still return. Requiring
    U ≤ 10⁻⁴ at the candidate t_end is a **gate condition** of the pilot (G-TD), **not** a new B4 pass threshold.
- **τ_w = 1.0** is ≈ 1.5× the largest late-time decay time quoted in the protocol (0.19) and in the audit's fits
  (≤ 0.67). It also exceeds the round-trip time from the slab to the right CAP and back for kₓ ≥ 1 (≤ 0.8).

**Pilot output:** for each realization, SETTLED or NOT SETTLED, together with t_s, t_first and U.

**Candidate B4 evaluation time (deterministic):** t_end,cand = max over the four realizations of t_s. It is defined
only if all four are SETTLED.

- G-TD then requires t_max − t_end,cand ≥ τ_w.
- The B4 amendment adopts t_end,cand, or a later frame that satisfies the same trailing-window condition, as the time at
  which B4 reads T_TD, R_TD and P_in.

### 2.4 Robustness runs (TD only, informative)

- **TD-R1.** Configuration D (strongest disorder) to t_max at dt = 10⁻⁴. Purpose: dt error growth to t = 6, and whether
  the non-monotone P_slab structure is a time-step artifact.
- **TD-R2.** Configuration D at dx ≈ 1/128 (nx = round(3.9·128) − 1 = 498; dx = 3.9/499). Measured cost ≈ 46 s per
  1.5 time units, so ≈ 3 min. Purpose: whether the late-time P_slab structure is a spatial-resolution artifact.

**Frame alignment.** Production and both robustness runs record frames at the same times t_j = j·0.005, j = 0…1200:

- with dt = 2×10⁻⁴, frames are every 25 steps;
- with dt = 10⁻⁴, every 50 steps;
- TD-R2 keeps dt = 2×10⁻⁴.

The equality of the frame times is asserted by test. Frame-by-frame comparison is therefore defined.

**Predeclared robustness rule** (the only decision-relevant comparison). A robustness run is **ROBUST** if both hold:

- P_slab(t_j) ≤ ε = 10⁻⁴ for **every** frame with t_end,cand ≤ t_j ≤ t_max;
- t_max − t_end,cand ≥ τ_w.

Otherwise it is **NOT ROBUST**. This rule asks exactly the question that matters for B4: whether the stopping decision
taken on the production run would also hold under a refined dt or dx. It uses only the existing B4 value ε = 10⁻⁴ and
τ_w.

**Quantitative differences are report-only.** The following are reported but not used for any decision:

- Δt_s in frames;
- max_j |ΔP_slab(t_j)|;
- ΔT, ΔR and ΔU at t_end,cand.

No tolerance for them can be justified from the repository. V3 and V4 measured T and R changes of order 10⁻⁵ at t = 1,
which do not transfer to P_slab at the 10⁻⁴ level up to t = 6. Their quantitative comparison is therefore **NOT
EVALUABLE**, and no "exact t_s equality" or "maximum frame shift" criterion is imposed.

---

## 3. Part 3 — CAP sensitivity pilot (TD only; separate from Part 2)

### 3.1 Baseline and the question

**V5 baseline** (`results/validation/v5_cap_reflection.csv`): clean guide, **mode-1** packets with σₓ = 0.2
(σ_k ≈ 3.54):

| k_x0 | 10 | 12 | 16 | 20 | 25 | 30 |
|---|---|---|---|---|---|---|
| spurious return | 3.1×10⁻² | 1.3×10⁻² | 1.3×10⁻³ | 1.0×10⁻⁴ | 1.6×10⁻⁴ | 5.2×10⁻⁴ |

**The question.** Disorder converts flux into higher channels, whose longitudinal wave number kₓ,ₙ = √(2E − (nπ)²) is
small. At E = 205: n = 4, 5, 6 give kₓ ≈ 15.9, 12.8, 7.4 (arithmetic).

**Two assumptions are tested, not assumed:**

- **(A1)** For an x-only CAP in a clean guide, the transverse modes decouple, so the return depends only on the
  longitudinal kₓ distribution, *the same function for every mode*.
- **(A2)** V5's broadband values bound the return of narrowband low-kₓ flux. They may not: V5 averages over σ_k ≈ 3.5.

### 3.2 C-A: clean-guide, mode- and kₓ-resolved return (no disorder)

**Packets:**

- `packet = "waveguide_mode"`, `mode_n ∈ {1, 4, 6}` (all existing config fields).
- Narrowband σₓ = 1.0 (σ_k = 0.707).
- k_x0 ∈ {3, 5, 7, 9, 11, 13, 16, 20}.
- Mode 6 is open for E > 177.65; for a mode-n packet the relevant longitudinal energy is kₓ²/2.

**Domain:**

- left CAP [0, 0.6]; x₀ = 0.6 + 5.0;
- a right-CAP start at x₀ + 5.0;
- the production CAP (width 0.6, strength 300, power 2) on both ends;
- L_x = x₀ + 5.6, with nx chosen so that dx = 1/64 (or as close as `Lx/(nx+1)` allows; recorded);
- dt = 2×10⁻⁴;
- t_end = 5.0/max(k_x0 − 3σ_k, 0.5) + 1.0 (predeclared formula).

**Initial momentum distribution and free-propagation baseline** (derived from `initial_condition_numpy`).

- **Convention.** With `packet = "waveguide_mode"` the initial state is
  ψ₀ = (πσₓ²)^(−1/4) exp[−(x − x₀)²/(2σₓ²)] e^{ik_x0·x} · √(2/L_y) sin(nπy/L_y).
  - The Gaussian is in the **amplitude**, with width parameter σₓ, so |ψ₀|² has positional standard deviation σₓ/√2.
  - Its Fourier transform gives P(kₓ) = (σₓ/√π) exp[−(kₓ − k_x0)²σₓ²], a Gaussian with standard deviation
    σ_k = 1/(√2σₓ). This is the same P(kₓ) as in PHASE2_METHOD_SELECTION.md §3 B4.
- **Backward weight** of a run with parameters (k_x0, σₓ):
  W_back(k_x0, σₓ) = ∫_{−∞}^{0} P dkₓ = ½erfc(k_x0σₓ) = ½erfc(k_x0/(√2σ_k)).
  - **No extra √2 factor:** the argument k_x0σₓ already equals (k_x0/σ_k)/√2.
  - The value depends on each run's own (k_x0, σₓ). **It is never a single constant for all runs.**
- **Free-propagation baseline.** In free propagation on the infinite line, |ψ(x, t)|² is a Gaussian with mean
  x₀ + k_x0·t and variance s(t)² = σₓ²/2 + t²/(2σₓ²) (ħ = m = 1). The probability left of a plane x_c at time t is
  therefore
  B_free(x_c, t) = ½erfc[(x₀ + k_x0·t − x_c)/(√2·s(t))].
  - As t → ∞, B_free → W_back for any fixed x_c.
  - At finite t, B_free also contains **slow forward components** that have not yet passed x_c. These are not CAP
    returns.
  - **Assumptions:** infinite line, free propagation, no CAP and no walls. These are listed as OPEN-8.
  - The pilot also computes the grid state's backward weight W_back,grid, by a zero-padded FFT of the sampled ψ₀ along
    x, and reports |W_back,grid − W_back|.

**Sanity-check values for the C-A runs** (σₓ = 1.0, σ_k = 0.7071; t_end from the formula above; arbitrary-precision
erfc). These are arithmetic only:

| k_x0 | t_end | W_back = ½erfc(k_x0σₓ) | B_free(x₀, t_end) |
|---|---|---|---|
| 3 | 6.690 | 1.105×10⁻⁵ | 1.358×10⁻⁵ |
| 5 | 2.737 | 7.687×10⁻¹³ | 1.551×10⁻¹¹ |
| 7 | 2.025 | 2.092×10⁻²³ | 3.463×10⁻¹⁹ |
| 9 | 1.727 | 2.069×10⁻³⁷ | 1.627×10⁻²⁸ |
| 11 | 1.563 | 7.204×10⁻⁵⁵ | 1.557×10⁻³⁹ |
| 13 | 1.460 | 8.698×10⁻⁷⁶ | 2.937×10⁻⁵² |
| 16 | 1.360 | 1.164×10⁻¹¹³ | 1.463×10⁻⁷⁴ |
| 20 | 1.280 | 2.698×10⁻¹⁷⁶ | 2.500×10⁻¹¹⁰ |

At k_x0 = 3, subtracting only W_back would book ≈ 2.5×10⁻⁶ of not-yet-departed forward probability as CAP return;
B_free removes it.

**Spurious return:**

- **Primary (baseline-subtracted):** R_spur = A_L(t_end) + P(x < x₀, t_end) − B_free(x₀, t_end), using the run's own
  (x₀, k_x0, σₓ).
- **Report-only:**
  - the time-gated form [A_L(t_end) − A_L(t_g)] + P(x < x₀, t_end), where t_g is the time at which the trailing edge
    (x₀ − 5σₓ, moving at k_x0 − 3σ_k) passes x₀. It is not used for any decision, because P(x < x₀, t_end) still
    contains slow forward probability;
  - the ungated R_raw = A_L(t_end) + P(x < L_x − 0.6, t_end);
  - W_back, W_back,grid and B_free(x₀, t_end).

**Check of A1:** for each k_x0, |R_spur(n = 4) − R_spur(n = 1)| and |R_spur(n = 6) − R_spur(n = 1)|. Separability
predicts agreement to TD rounding. A difference larger than 10⁻¹⁰ (§0.1) falsifies A1 in this implementation and must
be reported.

**Output:** R_spur(k_x0, n), a narrowband CAP return curve for the production CAP, and for the CAP variants of §3.3 at
k_x0 ∈ {3, 7, 11}.

**C-A is not the calibration of the directional estimator.** Its packets are wider than its gap windows would allow.
The estimator is calibrated in §3.3.5 on clean (W = 0) runs in the C-B geometry.

### 3.3 C-B: the four B4 realizations with CAP variants

#### 3.3.1 Variants and the shifted frame

**Why integer-step shifts.** Extending the left CAP moves the domain origin, and the disorder field depends on the
absolute physical x. Every shift is therefore an **integer number m of production grid steps**, dx = 3.9/250 = 0.0156,
so that grid points of different variants coincide in physical coordinates.

| Variant | m | Shift s = m·dx | CAP width (each side) | Strength | Power | L_x | nx | Gaps (physical) |
|---|---|---|---|---|---|---|---|---|
| V0 (production) | 0 | 0 | 0.6 | 300 | 2 | 3.9 | 249 | left 1.3, right 0.4 |
| V1 (wide) | 39 | 0.6084 | 0.6 + s = 1.2084 | 300 | 2 | 5.1168 | 327 | unchanged |
| V2 (wide, weaker) | 39 | 0.6084 | 1.2084 | 150 | 2 | 5.1168 | 327 | unchanged |
| V3 (long gaps, production CAP) | 64 | 0.9984 | 0.6 | 300 | 2 | 5.8968 | 377 | left 2.2984, right 1.3984 |
| V4 (very wide) | 115 | 1.7940 | 0.6 + s = 2.3940 | 300 | 2 | 7.4880 | 479 | unchanged |

- In V1, V2 and V4 the inner CAP edges stay at the production physical positions (x = 0.6 and x = 3.3).
- **Design-time arithmetic** (IEEE double, this platform): for m = 39, 64 and 115, L_x/(nx + 1) equals the production
  dx **bitwise**. For m = 38 it does not, which is why m = 39 is used instead of a nominal 0.6.
- The pilot re-checks this equality by test. If it fails on the execution platform, the kinetic-operator dx differs
  at rounding level. That is recorded, and the potential identity below is unaffected.

**Construction (no source change).**

- **Physical coordinates.** x_phys,j = (j + 1 − m)·dx is built from integer indices and the production dx, using the
  same expression as the production grid, (i + 1)·dx with i = j − m. It is **not** computed as x_j − s, which would
  round differently.
- **Potential.** The potential array is V_j,l = Potential(cfg, disorder_from_config(cfg, 0)).numpy(x_phys,j, y_l), with
  the production y-grid. It is passed as `SplitStepDST(cfg_s, potentials=V)`.
- **Initial state.** ψ₀ is built by `initial_condition_numpy(cfg, x_phys, y)` with the **production** cfg (x₀ = 1.2) and
  passed as `run(psi0=...)`.
- **Region masks** (slab, gaps and windows) are defined in physical coordinates.

#### 3.3.2 Same realization: three levels, all required

Equal parameter hashes alone do **not** show that two differently sampled grids represent the same field.

| Level | Question | Check |
|---|---|---|
| **L1 parameter identity** | the same random draw and field definition | sha256 over (qx, qy, φ), together with (W, ζ, region, edge, M, seed), equal for all variants, the RGF devices and the B3 records |
| **L2 sampler fidelity** | each grid samples *that* field correctly | for each variant grid and each RGF device, the stored samples are **bitwise** equal to a direct `Potential(...).numpy` call at the same coordinates (as in B3) |
| **L3 common-coordinate identity** | the same physical values where grids share physical points | V_variant[j = i + m, :] equals V_production[i, :] **bitwise** for every production index i; ψ₀ likewise; for the RGF, the coincident slices across dx are bitwise equal (B3Y cross-resolution identity) |

- **Where L3 cannot apply.** The TD grid (x spacing 0.0156) and the RGF slices (1/512 and finer) share no common
  physical points. For RGF ↔ TD only L1 and L2 apply, and **no** claim of equal sampled values on different grids is
  made. That the two discretizations represent the same continuum field is then a consequence of L1 and L2 together.
- **Points outside the production domain.** Variant grid points outside [0, 3.9] (negative physical x, or beyond 3.9)
  have no production counterpart. Only L2 applies there. Their |V| is recorded; it is a tail of the tanh window.
- **Bitwise L3 is a test, not an assumption.** `DisorderRealization.numpy` reduces over chunks of 64 modes per point.
  Whether the result is independent of the array shape must be verified by the test.

#### 3.3.3 Runs

- All five variants for the four realizations.
- Production dx, dy and dt; t_max = 6.0.
- Scalar frames every Δt_f = 0.005 (as in Part 2).
- **Directional frames every Δt_d = 0.001** (every 5 time steps), computed by the `observables` callback.
- **Clean reference:** one W = 0 run per variant geometry, with the same packet and settings. It is used for the
  calibration in §3.3.5.

#### 3.3.4 Comparison time: separating CAP return from timing

**Why a common comparison time is needed.** The variants differ in CAP position, CAP gradient and gap length, so at a
fixed nominal time they can differ in how much probability has been *absorbed* versus how much is still in transit.

- By the existing definitions, T = A_R + P_gapR + P_capR and R = A_L + P_gapL + P_capL. Absorption timing therefore
  does **not** change T or R directly: probability in a gap or CAP is already booked to the side it is on.
- Timing enters only through probability whose final assignment is still open:
  - the slab content P_slab;
  - the unabsorbed exterior content U_L = P_gapL + P_capL and U_R = P_gapR + P_capR, which can still return.

**Procedure:**

1. **Per variant j.** Compute t_s,j and SETTLED_j by the §2.3 rule (ε = 10⁻⁴, τ_w = 1.0, t_max = 6.0).
2. **Common comparison time** for a set of compared variants: t_c = max_j t_s,j.
   - The comparison is **EVALUABLE** only if every compared variant is SETTLED and t_max − t_c ≥ τ_w.
   - Otherwise it is **NOT EVALUABLE**, and no difference is interpreted.
3. **Undetermined probability** of variant j at t_c: u_j = P_slab,j(t_c) + U_L,j(t_c) + U_R,j(t_c). The asymptotic
   values satisfy:
   - T_j(∞) ∈ [T_j(t_c) − U_R,j, T_j(t_c) + P_slab,j + U_L,j];
   - R_j(∞) ∈ [R_j(t_c) − U_L,j, R_j(t_c) + P_slab,j + U_R,j].

   These intervals follow from probability bookkeeping alone, with no amplitude argument:
   - absorbed probability is final, so T_j(∞) ≥ A_R,j(t_c) and R_j(∞) ≥ A_L,j(t_c);
   - region probabilities are non-negative;
   - T + R + P_slab is conserved, with P_slab(∞) = 0.

   Interference between components therefore cannot move a variant's asymptotic T or R outside its interval.
4. **Pairwise decision** for variants i, j, with Δ = |T_i(t_c) − T_j(t_c)| (the same for R) and u_ij = u_i + u_j:
   - **PASS** if Δ + u_ij ≤ 1×10⁻⁴;
   - **FAIL** if Δ − u_ij > 1×10⁻⁴;
   - otherwise **NOT EVALUABLE**.

   All undetermined probability is counted against the criterion, so a timing difference cannot be read as a CAP
   effect and cannot produce a false PASS. Because |T_i(∞) − T_j(∞)| ≤ Δ + u_ij and ≥ Δ − u_ij, PASS and FAIL are
   rigorous statements about the **difference** of the asymptotic T and R between two CAP variants.

   **What this measures.** It is the direct empirical test of the sensitivity of T and R to the CAP. It is **not** an
   absolute bound on the CAP-induced error of any one variant. A return common to all compared variants (common mode)
   would not appear in the differences (OPEN-13).
5. **Report-only.** The same quantities at t_max, and the full time series T_j(t), R_j(t), P_slab,j(t), U_L,j(t) and
   U_R,j(t).

#### 3.3.5 Direction-resolved probability (reproducible definition)

**Where it can be evaluated.** In a gap with no scatterer, left-moving probability in the right gap can only come from
the right CAP (or the wall behind it). Right-moving probability in a window to the left of the initial packet can only
come from the left CAP.

**Windows** (physical coordinates; a = 1.9, b = 2.9):

| Window | Variant | Support | Length L_w | Spectral resolution 2π/L_w (arithmetic) |
|---|---|---|---|---|
| right, W_R | V3 | [b + 0.10, x_capR − 0.02] = [3.0, 4.2784] | 1.2784 | ≈ 4.9 |
| right, W_R | V0 | [3.0, 3.28] | 0.28 | ≈ 22.4 |
| left, W_L | V3 | [x_capL + 0.02, x₀ − 5σₓ] = [−0.3784, 0.2] | 0.5784 | ≈ 10.9 |
| left, W_L | V0 | none: the initial packet extends into the left CAP | — | — |

**Why these window positions:**

- **0.10 from the slab edge** is 5 disorder-edge widths (edge 0.02), so the disorder tail factor is ≤ e⁻¹⁰ ≈ 4.5×10⁻⁵.
  The recorded max|V| in each window is reported.
- **0.02 from a CAP edge** places the window outside the CAP, where W_cap is exactly 0 by `cap_profile`.
- **The left window ends at x₀ − 5σₓ**, so the initial |ψ|² inside it is ≤ ½erfc(5) ≈ 7.7×10⁻¹³.
- Each window has cos² tapers of width 0.05 inside its support. w(x) is the taper function, and
  L_eff = Σ_x w(x)²·dx is its effective length.

**Spectral transform** (per directional frame):

- Mode amplitudes ψ_n(x, t) are obtained by a DST-I in y, orthonormal, so Σ_n|ψ_n|² = Σ_y|ψ|²dy.
- The window field φ_n = w·ψ_n is zero-padded to N_pad = the smallest power of two ≥ 8M (M window points) and
  FFT-transformed.
- k_m = 2πm/(N_pad·dx), signed. The normalization is chosen so that Σ_m|φ̂_n(k_m)|² = Σ_x|φ_n(x)|²dx (Parseval
  exact). The k = 0 bin is assigned to neither direction and reported as D⁰.

**Three quantities, kept distinct:**

| Quantity | Definition | Units | Role |
|---|---|---|---|
| **Instantaneous directional probability** | D⁻_R(t) = Σ_n Σ_{k_m<0} \|φ̂_n\|² in W_R; D⁺_L(t) the same with k_m > 0 in W_L | probability | diagnostic time series |
| **Maximum directional probability** | D⁻_R,max = max over directional frames in [0, t_max] of D⁻_R(t); likewise D⁺_L,max | probability | **diagnostic only** (see "Why D_max is not a lower bound" below) |
| **Time-integrated directional probability** | F⁻_R = ∫₀^{t_max} J⁻_R dt, with J⁻_R(t) = (1/L_eff) Σ_n Σ_{k_m<0} \|k_m\| \|φ̂_n\|² (density × velocity, ħ = m = 1); likewise F⁺_L | probability | the **estimator of the returned probability**; F_up is a **conditional upper bound on returned probability only** (conditional on SD and on the estimator assumptions). It is **not** a bound on the CAP-induced change of T or R (L-iv) |

- Only the time-integrated quantity has the meaning "probability returned by the CAP". Integrating D(t) over time
  would have units of probability × time and **cannot** support a 10⁻⁴ probability criterion.
- **Integration:** composite trapezoid over the directional frames (Δt_d = 0.001), over [0, t_max].
- **Sampling uncertainty:** δ_s = |F(Δt_d) − F(2Δt_d)|, where F(2Δt_d) uses every second frame.
- **Late returns:** δ_late = U_R(t_max) for F⁻_R and U_L(t_max) for F⁺_L. Unabsorbed content may still return after
  t_max.

**Leakage in both directions** (computed on the actual disordered frames). A finite window misassigns content both
ways:

- **(a1) Right-moving → k < 0 (spurious return).**
  - φ̃⁺ = w·IFFT[Θ(k > 0)·FFT(φ)] is the right-moving part, re-windowed.
  - D_leak⁺⁻(t) and F_leak⁺⁻ are the k < 0 measures of φ̃⁺, defined as above.
- **(a2) Left-moving → k > 0 (return missed).**
  - φ̃⁻ = w·IFFT[Θ(k < 0)·FFT(φ)] is the left-moving part, re-windowed.
  - D_leak⁻⁺(t) and F_leak⁻⁺ are its k > 0 measures.
  - Without this term, F⁻ could *underestimate* the return, so F_up would not be an upper bound.
- The left window W_L uses the mirrored definitions.
- **Both terms are model estimates of misassignment, not bounds.** Whether they bound the true misassignment is tested
  by SD (below).

**Clean-reference calibration** (one W = 0 run per geometry; the production packet, k₀ = 20, σₓ = 0.2, x₀ = 1.2
physical):

- **Direct clean return:** R_ret,clean = A_L(t_max) + P(x < x_capR, t_max) − B_free(x_capR, t_max), with B_free from
  §3.2 evaluated with **this run's** (x₀, k₀, σₓ).
  - For V3 (x_capR = 4.2984, t_max = 6): B_free = 1.79×10⁻⁸, of which W_back = ½erfc(4) = 7.71×10⁻⁹. This is
    arithmetic.
  - The previous form subtracted only W_back, and so booked the slow-forward part of B_free (≈ 1.0×10⁻⁸ here) as a
    return.
- **Estimator:** F⁻_R,clean, computed as above.
- **Calibration error:** ε_cal = |F⁻_R,clean − R_ret,clean|, reported and added as uncertainty.
- **Pre-return check.**
  - Before t_arr = (x_capR − x₀ − 5σₓ)/(k₀ + 5σ_k), no component within 5σ_k can have reached the right CAP
    (analytic). For V3, t_arr ≈ 0.0557.
  - Over [0, t_arr), the clean run's k < 0 content in W_R must be ≤ D_leak⁺⁻(t) frame by frame. If not, the leakage
    estimate (a1) is invalid and (iii) is **NOT EVALUABLE**.
- **Limitation:** the clean calibration contains only mode-1, broadband flux. It does not represent low-kₓ
  mode-converted flux.

**Why D_max is not a lower bound** (examined; not assumed). D⁻_R − D_leak⁺⁻ would be a lower bound on probability
returned by the right CAP only if all of the following held:

- **(L-i)** The window contains no scatterer. This holds up to the disorder tail (factor ≤ e⁻¹⁰), whose own
  reflection is not bounded here.
- **(L-ii)** The window contains no evanescent (closed-mode) content. Real-exponential tails have a spectrum symmetric
  about k = 0 and contribute to D⁻ without being left-moving. For a closed mode near its threshold, the decay length
  can exceed the 0.10 offset from the slab, and the repository gives no bound on it.
- **(L-iii)** D_leak⁺⁻(t) bounds the actual right→left misassignment at every frame. It is a model estimate (a1); no
  bound is established.
- **(L-iv)** Returned probability can be identified with the error the CAP causes in T or R. **It cannot, in either
  direction.**
  - **Not a lower bound.** A returned wave can re-scatter and end on the same side, so the effect on T or R can be far
    smaller than the returned amount.
  - **Not an upper bound.** The returned wave adds **coherently** to the amplitudes that are already there. If r is an
    outgoing reflection amplitude (energy- and channel-resolved) and δr the CAP-induced returned contribution in the
    same channel, then
    ΔR = |r + δr|² − |r|² = 2 Re(r*·δr) + |δr|².
    With |δr|² = p small, the interference term can be as large as 2|r|√p, which is ≫ p. The same holds for
    transmission amplitudes.
  - **Illustration only, not adopted.** Under an idealized linear decomposition, which the pilot does not establish
    (ψ = ψ_ideal + δψ, with all of the returned probability in δψ), Cauchy–Schwarz would give
    |ΔR| ≤ 2√(R·p) + p. At p = 10⁻⁴ this is of order 2×10⁻²·√R, far above 10⁻³. So even then a returned-probability
    bound at the 10⁻⁴ level would **not** support the 10⁻⁴ CAP budget item.

(L-ii) and (L-iii) cannot be established from repository evidence, and (L-iv) is false in general, in both
directions. **D_max is therefore a diagnostic only. It is not used to declare G-CAP FAIL, and no lower bound F_low is
used.** This is recorded as OPEN-9.

**Synthetic directional tests SD** (no physics; required before (iii) is evaluated). Each test uses analytic fields on
the actual window grids (V3 W_R and W_L; V0 W_R), per transverse mode, with known exact directional content and exact
flux through the window over time (free propagation):

| Test | Field | Checked |
|---|---|---|
| **SD1** one-directional outgoing | right-moving Gaussian packets, kₓ from the set below; the exact left-moving content is 0 | D⁻ ≤ D_leak⁺⁻ at every frame, and F⁻ ≤ F_leak⁺⁻ (validity of the spurious-content estimate used in the pre-return check) |
| **SD2** known incoming | left-moving packets of known probability p and known flux, same kₓ set | F⁻ + F_leak⁻⁺ ≥ exact returned flux (the upper-bound property needed for PASS) |
| **SD3** mixed direction | coherent superposition of SD1 and SD2 in the same mode, relative phases {0, π/2, π}, overlapping in time | both properties above under interference |
| **SD4** finite-window leakage | packets straddling the window edges, and kₓ below the resolution 2π/L_w | both properties above at and below resolution |
| **SD5** evanescent | real-exponential profiles e^{−κ(x − x_e)} with decay lengths 1/κ ∈ {1, 1/3, 1/10, 1/30} (test parameters spanning the window length and the 0.10 offset) | the size of the spurious directional content (reported; quantifies L-ii) |

- **kₓ set for SD1–SD4:** {0.25, 0.5, 1, 2, 4}·k_res and {7, 12.8, 15.9, 20} (the channel values of §3.1), where
  k_res = 2π/L_w is the window's resolution.
- **Outcome.** If SD2, SD3 or SD4 shows any case where F⁻ + F_leak⁻⁺ underestimates the exact returned flux, the
  upper-bound use of F is not established for that kₓ range, and (iii) is **NOT EVALUABLE** for flux in that range.

**Criterion (iii), evaluated on V3** (production CAP, resolvable windows):

- **Conditional upper bound on the returned probability** (conditional on SD and on the estimator assumptions). **It
  is not a bound on the CAP-induced change of T or R** (L-iv):
  F_up = F⁻_R + F_leak,R⁻⁺ + F⁺_L + F_leak,L⁺⁻ + δ_s,R + δ_s,L + δ_late,R + δ_late,L + ε_cal.
  - **Notation:** superscript "ab" means content moving in direction a, assigned by the window to direction b.
  - Only the **missed**-content terms are added: F_leak,R⁻⁺ (returned, left-moving content in W_R assigned to k > 0)
    and F_leak,L⁺⁻ (returned, right-moving content in W_L assigned to k < 0).
  - The spurious-content terms (F_leak,R⁺⁻, F_leak,L⁻⁺) already inflate the measured F and are conservatively **not**
    subtracted.
- **PASS** if SD passed for the relevant kₓ range **and** F_up ≤ 1×10⁻⁴; otherwise **NOT EVALUABLE**.
  - A PASS of (iii) means only that *the returned probability* is conditionally ≤ 10⁻⁴.
  - It does **not** mean that the CAP changes T or R by ≤ 10⁻⁴. By (L-iv), no such implication holds.
- **(iii) has no FAIL outcome.** A large F (or D_max) is reported as "CAP return not bounded below 10⁻⁴", which is
  diagnostic. A G-CAP FAIL can only come from (i) or (ii).
- **V0** (right window only, resolution ≈ 22.4) is **report-only**. A V0 value is interpreted only if its F_up itself
  is ≤ 10⁻⁴.
- The mode-resolved consistency estimate Σₙ∫wₙ(kₓ)R_spur(kₓ)dkₓ, using the C-A curves, is report-only.

#### 3.3.6 Predeclared CAP criterion (gate G-CAP)

For every realization, all of the following:

- **(i)** V0 vs each of V1–V4: PASS for T and for R by the §3.3.4 rule.
- **(ii)** Each pair among V1, V2 and V4: PASS for T and for R by the §3.3.4 rule. This shows that the enlarged CAPs
  agree with each other, so the V0 comparison is against a mutually consistent comparison set. This agreement does not
  show that any variant is free of CAP-induced error.
- **(iii)** PASS by §3.3.5. It requires SD to be passed and F_up ≤ 10⁻⁴; it has no FAIL outcome.
  - (iii) is a corroborating diagnostic of the returned probability under the assumptions of §3.3.5. It does not
    constrain a common-mode CAP-induced change in T or R, which may remain invisible to (i) and (ii).
  - It does **not** independently bound the CAP contribution to T or R.
- **A1:** separability holds (§3.2).

**Outcome.** Any FAIL means G-CAP FAILS. Any NOT EVALUABLE, without a FAIL, means G-CAP is **NOT EVALUABLE**, which
counts as not passed (§6).

The value 10⁻⁴ is the CAP item of the **existing** B4 error budget (≈ 10⁻⁴), i.e. 10 % of the 10⁻³ tolerance. It is
used as a gate on the evidence for the budget assumption. It is **not** a new B4 pass threshold, and it is not relaxed.

**What a G-CAP PASS establishes, and what it does not.**

- **It establishes:**
  - **(a)** rigorously (§3.3.4), that the asymptotic T and R of the production CAP differ from those of the enlarged
    CAP variants, and the enlarged variants from each other, by ≤ 10⁻⁴;
  - **(b)** conditionally (§3.3.5), that the probability returned by the production CAP is ≤ 10⁻⁴.
- **It does not establish** a rigorous upper bound on the absolute CAP-induced error of T or R, for two reasons:
  - **(a)** bounds only differences between variants (a common-mode return is invisible to it);
  - **(b)** bounds a probability, not an amplitude-level effect (L-iv).
- **Consequence for B4.** A PASS is evidence *consistent with* the budget item "CAP ≈ 10⁻⁴". It is **not a proof** of
  it, and the B4 amendment must say so (OPEN-13).

### 3.4 Code changes needed later (new files only; no existing file modified)

**Already possible with existing interfaces:**

- mode-resolved *region* probabilities and direction-resolved gap content, through the `observables` callback of
  `SplitStepDST.run`;
- shifted-frame potentials, through `potentials=array`.

**Mode-resolved *absorbed* probability:**

- This is not available: `SplitStepDST.run.absorb` sums over all modes.
- **Required:** a new module `b4_pilot_td.py` with a subclass whose `run` is a copy of `SplitStepDST.run`, in which
  `absorb` also accumulates Σ_y-DST|·|²·loss_frac per mode on the CAP slices.
- **Required test:** with mode resolution off, the subclass must reproduce `SplitStepDST.run` **bitwise**
  (T, R, norm, absorbed_left, absorbed_right at every frame), for a disordered and a clean configuration.
- `reference_solver.py` itself is **not** changed.

---

## 4. Part 4 — Threshold quadrature audit

### 4.1 Design-time enumeration (arithmetic; the pilot must reproduce it exactly)

**Coarse grid** (1274 nodes, h = 0.024996): 11 nodes are moved by `adjust_energy`.

| i | kₓ | kₓ' (moved) | Δkₓ (in units of h) | E → E' | n | side | side kept | trapezoid weight h·P(kₓ) |
|---|---|---|---|---|---|---|---|---|
| 54 | 5.439877 | 5.434138 | −0.230 | 19.7309 → 19.6997 | 2 | below | yes | 5.9×10⁻⁷ |
| 192 | 8.889314 | 8.895757 | +0.258 | 44.4448 → 44.5020 | 3 | above | yes | 2.0×10⁻⁵ |
| 323 | 12.163779 | 12.154351 | −0.377 | 78.9136 → 78.7989 | 4 | below | yes | 2.4×10⁻⁴ |
| 452 | 15.388253 | 15.374558 | −0.548 | 123.3340 → 123.1233 | 5 | below | yes | 1.2×10⁻³ |
| 580 | 18.587731 | 18.605020 | +0.692 | 177.6867 → 178.0082 | 6 | above | yes | 2.6×10⁻³ |
| 707 | 21.762212 | 21.743362 | −0.754 | 241.7317 → 241.3217 | 7 | below | yes | 2.5×10⁻³ |
| 834 | 24.936694 | 24.960937 | +0.970 | 315.8542 → 316.4590 | 8 | above | yes | 1.1×10⁻³ |
| 960 | 28.086180 | 28.070794 | −0.616 | 399.3516 → 398.9195 | 9 | below | yes | 2.1×10⁻⁴ |
| 961 | 28.111176 | 28.127695 | +0.661 | 400.0539 → 400.5184 | 9 | above | yes | 2.0×10⁻⁴ |
| 1087 | 31.260662 | 31.290010 | +1.174 | 493.5493 → 494.4672 | 10 | above | yes | 1.8×10⁻⁵ |
| 1213 | 34.410147 | 34.379705 | −1.218 | 596.9639 → 595.9168 | 11 | below | yes | 7.0×10⁻⁷ |

**Fine (halving) grid** (2547 nodes, h/2): 18 nodes are moved.

- Displacements reach **2.44 grid steps**.
- **Three pairs of distinct nodes are mapped to the same energy:**

  | Nodes | Shared energy |
  |---|---|
  | 1920, 1921 | E' = 398.91954 |
  | 2174, 2175 | E' = 494.46718 |
  | 2425, 2426 | E' = 595.91684 |

- All moved nodes keep their side of the threshold (the open-channel count is unchanged).

**Design-time observations** (from the definitions only; to be confirmed by the pilot):

- **(a) Fixed energy shift.** The rule moves a node to Eₙ(1 ± 2×10⁻³). This is a fixed *energy* distance, so at high n
  the displacement exceeds one grid step: the moved node lands beyond its neighbour.
- **(b) Weight assignment.** The existing helper `energy_average` assigns the trapezoid weight h·P(kₓ) (½h at the
  ends) at the **original** node, but evaluates the integrand at the **moved** energy. The quadrature is therefore not
  a trapezoid rule on any single grid near a threshold.
- **(c) Contaminated halving difference.** The coarse and fine grids move *different* sets of nodes, with different
  displacements, and the fine grid duplicates energies. The δkₓ/2 difference therefore mixes genuine quadrature error
  with differences between the two grids' threshold treatments. **It cannot be assumed to measure convergence near
  the thresholds.**
- **(d) Weight involved.** The total trapezoid weight on moved coarse nodes is ≈ 8.1×10⁻³ of the packet weight
  (arithmetic sum of the column above).

### 4.2 Diagnostics

1. **Enumeration (G0).** Reproduce §4.1 exactly from `kx_window`, `adjust_energy` and `e_perp`: indices, kₓ', Δkₓ,
   side, weights and shared energies. Any difference stops the pilot.
2. **Synthetic quadrature test** (no physics; exact reference). It characterizes the **quadrature rules** on
   idealized shapes. **It does not show that T₁ or R₁ has these shapes, and it does not validate any physical
   resonance structure.**

   - **Integrand classes.** Each class is tested separately, at every threshold n = 2…11, with E = E(kₓ):

     | Class | f(E) | Tests |
     |---|---|---|
     | **S0** smooth | s(E) = 0.5 + 0.3 sin(E/40) (an arbitrary smooth reference; any smooth choice serves) | baseline behaviour of halving on a smooth integrand, the case B4's halving check already relies on |
     | **S1** continuity with a kink | s(E) + max(E − Eₙ, 0)/Eₙ | continuous, with finite one-sided derivatives: separates "continuous but not smooth" from a cusp |
     | **S2** square-root opening cusp | s(E) + √(max(E − Eₙ, 0)/Eₙ) | the channel-opening behaviour stated in the protocol (√ cusp above threshold) |
     | **S3** two-sided square-root cusp | s(E) + [√(max(E − Eₙ, 0)) − √(max(Eₙ − E, 0))]/√Eₙ | cusp structure on both sides (closed-channel influence below threshold) |
     | **S4** narrow resonance | s(E) + Γ²/((E − E_r)² + Γ²), with E_r = Eₙ − 3Γ and Γ ∈ {0.1, 1, 10}·δEₙ, where δEₙ = kₓ(Eₙ)·h is the local coarse-grid energy spacing at that threshold | sensitivity to features narrower than, equal to and wider than the grid. The widths are defined relative to the grid, not to any physical width |

   - **Amplitude.** Unit amplitude. Every rule is linear in f, so errors and halving differences scale with the
     amplitude, and their ratio does not depend on it. No amplitude is chosen.
   - **Reference.** ∫P·f dkₓ by adaptive quadrature split at Eₙ (and at E_r for S4), with tolerance ≤ 10⁻¹³.
   - **Rules compared**, on the coarse and fine grids:
     - **R-a:** the existing rule (moved node, original weight);
     - **R-b:** plain trapezoid with no move;
     - **R-c:** a trapezoid split at the thresholds (piecewise grids ending exactly at kₓ(Eₙ)), as a possible amended
       method.
   - **Outputs per class, rule and threshold:**
     - the true error of the coarse result, e_c;
     - the halving difference H = |coarse − fine|;
     - the detection ratio r = H/|e_c|.
   - **Classification.** Halving is **DETECTING** for a class and rule if r ≥ r_S0, i.e. it detects at least the
     fraction of the error it detects for the smooth integrand, which B4 already relies on. It is **UNDER-DETECTED** if
     r < r_S0.
     - For S4, the smallest tested Γ at which halving is DETECTING is reported.
     - No numerical ratio other than r_S0 is introduced.
   - **What the classes separate.** S0 vs S1 shows the effect of a derivative jump; S1 vs S2/S3 the effect of an
     infinite derivative (cusp); S4 the effect of narrow structure. A rule can be DETECTING for S1–S3 and
     UNDER-DETECTED for narrow S4; these are reported separately.
3. **RGF local refinement windows** (per configuration; M1 only; at the resolution under consideration from Part 1, or
   dx = 1/1024 if undecided).
   - **Windows:** thresholds n ∈ {4, 6, 8, 10}, each ±8 coarse steps around kₓ(Eₙ).
   - **Grids:** h, h/2 and h/8 restricted to the window, with the rules R-a and R-c.
   - **Output:** only *differences* between rules and grids of the local window integral of P·T₁ and P·R₁, and the
     ratio as in item 2.
   - **Absolute window integrals are not reported**, so the pilot cannot be read as partial B4 averages.
   - **Resonance scan:** at the threshold of largest packet weight, n = 6, T₁ and R₁ at spacing h/32 over ±2 coarse
     steps.
     - Reported: the local window integral **difference** |I(h/2 subset of the scan) − I(h/32)|, and the narrowest
       feature resolved.
     - This real-RGF diagnostic is **separate** from the synthetic test, and it is the only resonance-sensitive element
       of the pilot.
     - Resonances at other thresholds, or narrower than h/32, remain **NOT VALIDATED**.
   - **Cost ESTIMATE:** about 4 windows × (17 + 33 + 129) points × 4 configurations ≈ 2900 RGF evaluations, at
     0.8 s ≈ 40 min.

**Not assumed.** The existing B2 halving change (2.3×10⁻¹⁰) concerned a single-channel ODE with no threshold. It does
not show that halving is reliable near the multichannel thresholds.

---

## 5. Part 5 — Outputs, separation and provenance

### 5.1 Output directory

```
results/b4_pilot/
  README.md                       "DIAGNOSTIC / NOT B4" banner, protocol hash, what is NOT computed
  rgf_resolution/   diagnostic_rgf_points.csv, diagnostic_rgf_S.npz, diagnostic_rgf_orders.csv, diagnostic_ny_nq.csv
  td_stopping/      diagnostic_td_frames_<cfg>.npz, diagnostic_td_settling.csv, diagnostic_td_robustness.csv
  cap_sensitivity/  diagnostic_cap_clean_scan.csv, diagnostic_cap_variants.csv, diagnostic_cap_direction.csv
  threshold_quadrature/ diagnostic_moved_nodes.csv, diagnostic_cusp_synthetic.csv, diagnostic_local_windows.csv
  gates.json                      gate evaluations (§6), each "DIAGNOSTIC / NOT B4"
  run_metadata.json
```

### 5.2 Labelling

- Every CSV has a first comment line `# DIAGNOSTIC / NOT B4 — docs/B4_PILOT_PROTOCOL.md`.
- Every JSON has `"classification": "DIAGNOSTIC / NOT B4"`, and every `.npz` has a `classification` array with that
  string.
- File names begin with `diagnostic_`.

### 5.3 Separation (fail-closed)

The pilot driver refuses to run if:

- the output path is not under `results/b4_pilot/`;
- the output path is non-empty, unless `--resume` is given with a cache whose protocol hash matches;
- any write would target `results/rgf_b*_validation/`, `results/stationary_validation/`, `results/rgf_b4_validation/`,
  `results/transport/` or `results/validation/`.

**Required tests:**

- the path guard;
- that no pilot code path calls `energy_average` on RGF data over the full B4 window;
- that no output contains both an RGF and a TD quantity for the same realization;
- that the committed B1–B3Y artifacts are unchanged (sha256 before and after the run, recorded in the metadata).

### 5.4 Metadata (`run_metadata.json`)

- the sha256 of this protocol file;
- the git commit and a clean-tree flag; **a dirty tree aborts the run**;
- `provenance.code_fingerprint` of every module used: `stationary`, `smatrix_spec`, `sparse_reference`,
  `sparse_wavefunction`, `sparse_exact_lattice`, `reference_solver`, `potential`, `disorder`, `config`,
  `initial_condition`, `validate_rgf_b2` (helpers), `validate_rgf_b3` (geometry), and the pilot modules;
- `provenance.environment_info`, `OPENBLAS_NUM_THREADS`, the CPU and core count, and the M4 platform-precision report;
- **per configuration:** the full config dict and the disorder-parameter sha256 over (qx, qy, φ) (L1);
- the potential-sample sha256 for every RGF device and every TD grid, including the shifted-frame variants. Each is
  recorded with the result of its bitwise direct-call check (L2);
- the result of every common-coordinate bitwise check (L3), listing the index offset m and the compared index range;
- for each shifted variant: m, s, L_x, nx, and whether L_x/(nx + 1) equals the production dx bitwise;
- the command line, wall times, the cache location and the cache-key scheme.

### 5.5 Reproducibility

- Deterministic inputs: no random number generator other than the disorder seed; seeded `onenormest` in M4, as in B3Y.
- A rerun from the cache must reproduce every CSV byte for byte, except the timing columns. This is checked once and
  recorded.
- One live re-run of a predeclared row (C1, config D, dx = 1/1024, M1) must equal the stored value bitwise.

---

## 6. Part 6 — Decision gates

**Gate G0 (preconditions) — blocking.** All of the following:

- the protocol hash matches;
- the clean tree;
- §1.2 and §4.1 reproduced exactly;
- the controls C1 and C2, and the three B3-energy re-runs of §1.6, reproduce the committed T₁ values (≤ 10⁻¹²);
- the committed artifacts are unchanged.

Any failure stops the pilot; there is no partial output.

**Amendment rule.** "Amendment" below means a **dated B4 protocol amendment**: a new document committed **before** any
B4 number. It records the pilot evidence and the change, and states that the B4 thresholds 10⁻³, 10⁻⁴ and 2×10⁻⁴ are
unchanged. **No gate may be passed by weakening any existing B4 threshold.**

| Gate | Required evidence | Pass condition | Fail condition | If it fails | Dated amendment required? |
|---|---|---|---|---|---|
| **G-REAL** (same realization) | the three levels of §3.3.2: **L1** parameter hashes and configuration equality across the RGF devices (all dx), the TD production grid, the CAP variants and the B3 records; **L2** bitwise sampler fidelity of every grid against direct `Potential.numpy` calls at its own coordinates (as in B3); **L3** bitwise equality at common physical coordinates (variant ↔ production TD grid; coincident RGF slices across dx), including ψ₀; the shifted-grid dx equality test | L1, L2 and L3 all true where applicable. RGF ↔ TD is judged on L1 + L2 only, and no equality of differently sampled grids is claimed | any mismatch at any level. **Equal hashes with an L2 or L3 mismatch is a FAIL**, not a pass | **blocking**: no B4 until the cause is found and fixed; any pilot output that depends on the mismatched item is void | no (unless the realization definition itself must change, which would be a new benchmark) |
| **G-RGF** (resolution) | Part 1 rows for all 14 points × 4 configurations × 3 dx × M1–M3 (M4 where specified), plus the three B3 energies re-run with M1 (§1.6) | (a) M1–M2–M3 agreement on S ≤ 10⁻¹⁰ (or adjudicated as in B3Y) at every point; **and** (b) the §1.6 algorithm yields a dx_sel for every configuration, i.e. a dx admissible (ADMISSIBLE, ADMISSIBLE-R or ADMISSIBLE-F) at every point of the §1.6 set, for T₁ and for R₁ separately. A resolved T₁ never substitutes for R₁; **and** (c) the N_y and N_q changes ≤ 2×10⁻⁴ at the top-window points; **and** (d) lead-tail guard met | (a) unexplained; or §1.6 step 3(b) for some configuration (all points classified, a limit exceeded even at 1/2048); or (c) not met. §1.6 step 3(a) (an unclassifiable point) gives **NOT EVALUABLE**, i.e. not passed | (a): treat as a solver defect, so B4 is blocked. (b): the B4 budget item cannot be met as designed for that configuration. Options, each needing a dated amendment: a finer dx after a separate dx = 1/4096 review (NOT EVALUABLE: a dated pilot amendment adding 1/4096); or a predeclared Richardson-extrapolated B4 quantity with the order gate; or recording that configuration's B4 outcome as discretization-limited. Never a relaxed 10⁻³. (c): a larger N_y or N_q in B4 | **yes, in every case**: the B4 text does not fix the RGF resolution, so even a pass requires the amendment to state the chosen dx, N_y and N_q |
| **G-TD** (stopping rule) | Part 2 frames for the four realizations; TD-R1 and TD-R2 | all four realizations SETTLED (§2.3); t_end,cand = max_r t_s,r with t_max − t_end,cand ≥ τ_w = 1.0; bookkeeping drift |ρ(t) − ρ(0)| ≤ 10⁻¹⁰ for all t and the consistency check |P_in,res − P_slab − ρ| ≤ 10⁻¹⁰ (§2.2, §0.1); ρ(0) reported; U(t_end,cand) ≤ 10⁻⁴ for every realization; TD-R1 and TD-R2 both ROBUST (§2.4) | any realization NOT SETTLED; or t_max − t_end,cand < τ_w; or U(t_end,cand) > 10⁻⁴; or an accounting check > 10⁻¹⁰; or TD-R1 or TD-R2 NOT ROBUST. The report-only differences of §2.4 never decide the gate | NOT SETTLED: B4 cannot meet its P_in ≤ 10⁻⁴ pass condition with the current setup. Options (dated pilot amendment, then a new pilot run): a larger t_max, or a geometry change (which would make it a different benchmark) | **yes**: the amendment must fix the B4 t_end (or rule), the frame spacing, and the frame at which T_TD, R_TD and P_in are read |
| **G-CAP** (CAP sensitivity) | C-A narrowband curves (modes 1, 4, 6); C-B variants V0–V4 for the four realizations, with settling times, the common comparison time t_c and the undetermined probabilities u_j (§3.3.4); the V3 direction-resolved quantities F, the two-way leakage terms, sampling, late-return and calibration terms (§3.3.5; D_max diagnostic only); the synthetic directional tests SD1–SD5; the clean-reference runs with B_free | §3.3.6 (i)–(iii) all PASS for every realization, by the conservative rules of §3.3.4–§3.3.5; A1 holds (≤ 10⁻¹⁰, §0.1). **A PASS means** that the T/R differences between CAP variants are ≤ 10⁻⁴ (rigorous) and that the returned probability is conditionally ≤ 10⁻⁴. It is **not** a proven bound on the absolute CAP-induced T/R error (§3.3.6; OPEN-13) | a FAIL of (i) or (ii) (a variant difference exceeding 10⁻⁴ even after subtracting all undetermined probability, §3.3.4); or A1 falsified. (iii) has **no FAIL outcome** (§3.3.5: returned probability does not bound the CAP effect from below). If there is no FAIL but any NOT EVALUABLE (an unsettled variant; t_max − t_c < τ_w; uncertainty terms too large to decide; SD not passed; or F_up > 10⁻⁴): **NOT EVALUABLE**, i.e. not passed | FAIL: the production CAP is not adequate for B4's budget. Options: B4 TD with a CAP variant that passes (i)–(iii) against the even larger variants, or the measured variant differences entered in the B4 budget, stating that they are not an absolute CAP error bound (OPEN-13). NOT EVALUABLE: a dated **pilot** amendment (e.g. longer t_max or longer gap windows), then a re-run; no interpretation of the differences. Never a relaxed 10⁻³ or 10⁻⁴ | **yes, if the B4 TD settings change**, because the B4 text says "production settings"; **no, if V0 passes** |
| **G-THR** (threshold quadrature) | §4.1 reproduced; the synthetic cusp test; the RGF local windows; the resonance scan | for the quadrature rule proposed for B4: halving DETECTING (r ≥ r_S0, §4.2) for S1, S2 and S3 at every threshold (S4 is reported, not gated); **and** in the real-RGF local windows plus the n = 6 resonance scan, the sum of the local differences \|rule(h) − rule(h/8)\| (and \|I(h/2) − I(h/32)\| for the scan), for T₁ and for R₁ separately, ≤ 2×10⁻⁴ (the existing quadrature item; a **necessary, not sufficient** condition, since the smooth part of the window is not included) | UNDER-DETECTED for S1–S3 with the proposed rule; or a local sum > 2×10⁻⁴ | adopt R-c (threshold-split quadrature) or local refinement for B4, re-verified by the same tests; if the scan shows structure that h/2 misses, local adaptive refinement near thresholds with a predeclared rule. A passing synthetic test never replaces the real-RGF diagnostic | **yes**, if the B4 quadrature changes from the protocol's uniform grid with moved nodes; **no**, if R-a passes |

**Proceeding to full B4 requires G0, G-REAL, G-RGF, G-TD, G-CAP and G-THR all PASS**, with the dated amendment
committed. A gate that is NOT EVALUABLE (for example the M4 platform gate failing) counts as **not passed**.

---

## 7. Part 7 — Final recommendation

### 7.1 Decision table

| Item | Already validated (existing, unchanged) | The pilot must establish | Remains unvalidated after a successful pilot |
|---|---|---|---|
| RGF solver on the lattice | B2 PASS; B3/B3X solver agreement (≤ 6×10⁻¹²); B3Y four-method agreement (≤ 9.1×10⁻¹¹) and H1 at its 16 points | agreement and lattice dx dependence at 14 points spanning the window (G-RGF) | the continuum limit; energies between the pilot points; other seeds |
| RGF resolution policy | B3 FAIL / B3X C11 FAIL (dx) at W = 200 — statuses unchanged | whether a dx (≤ 1/2048) meets the B4 budget item at all pilot points | continuum accuracy (Richardson values are estimates) |
| Transverse truncation | converged at E ≤ 371 (B3) | N_y and N_q at E 473–650 | N_y = 255 (not run) |
| TD solver | V1–V7; B2 TD vs ODE (barrier) PASS | settling of P_in below 10⁻⁴ with a trailing window; exterior content U | behaviour for other seeds; t > 6 |
| CAP | V5 (mode 1, broadband) | narrowband, mode-resolved return; the variant sensitivity at a common, settled comparison time; the V3 directional return estimate with leakage and calibration | CAP behaviour for other realizations; directional resolution for kₓ below the window resolution (OPEN-1); a rigorous bound on the absolute CAP-induced T/R error (OPEN-13) |
| Threshold quadrature | B2 halving (single channel, no thresholds) | the detection ratio of halving on synthetic kink, cusp and resonance shapes (rules only); real-RGF local window differences at n = 4, 6, 8, 10 and the n = 6 scan | the physical resonance structure in general; resonances narrower than h/32 or at unscanned thresholds |
| Same realization | identical disorder-parameter hashes for the RGF and TD paths (L1 only; audit §7) | L1 + L2 for every grid; L3 at every common physical coordinate (shifted TD variants, nested RGF slices) | equality of sampled values on grids without common points: not claimed (only L1 + L2) |
| B1-F | FAIL (provisional phase criterion), unchanged | rcond and U recorded at the near-threshold points at 1/2048 | phase accuracy on long fine grids (not used by B4) |
| **B4 itself** | **not performed** | **nothing; the pilot produces no B4 comparison** | **everything that B4 is designed to test** |

### 7.2 Prioritized implementation plan (none implemented)

1. **Commit this protocol.** Its sha256 becomes the G0 reference.
2. **Pilot tests** (new `tests/test_b4_pilot.py`): the G0 reproduction of §1.2 and §4.1; the region partition; the
   path guard and labels; the absence of full-window RGF averages; the bitwise equivalence of the subclass in §3.4;
   the trailing-window logic on synthetic P(t) curves (including non-monotone ones and late up-crossings); the
   L1/L2/L3 realization checks, including a mutation (shift by m ± 1) that L3 must detect; the dx-equality test of
   §3.3.1; the order-status classification of §1.4 on synthetic triples (BELOW FLOOR, one-sided floor, sign change,
   outside window, resolved), separately for T and R; the §1.6 decision algorithm on synthetic status tables
   (including the "-R", "-F" and "-U" cases and the coarsest-admissible rule); the frame-time equality of §2.4 and the
   ROBUST rule; the §3.3.4 PASS/FAIL/NOT EVALUABLE logic; the Parseval normalization, the two-way leakage terms and the
   directional tests SD1–SD5 of §3.3.5; B_free and W_back against the §3.2 table; and the synthetic-quadrature
   reference integrals (S0–S4).
3. **Pilot modules** (new): `b4_pilot_rgf.py` (Parts 1 and 4), `b4_pilot_td.py` (Parts 2 and 3, the subclass),
   `b4_pilot_summary.py` (gates). They reuse existing functions read-only.
4. **Run order:**
   - G0, then G-REAL;
   - Part 4 synthetic (seconds);
   - Part 1 (≈ 1 h ESTIMATE);
   - Part 2 (≈ 2 min, plus TD-R1/R2 ≈ 5 min);
   - Part 3 (C-A: 24 clean runs plus variant runs; ≈ 30–60 min ESTIMATE from the measured TD cost per time unit and
     the domain sizes);
   - Part 4 RGF windows (≈ 40 min).
5. **Pilot report** (`docs/B4_PILOT_REPORT.md`, a new file written after the run): every gate with its evidence,
   labelled DIAGNOSTIC / NOT B4.
6. **Dated B4 protocol amendment:** the RGF resolution, N_y and N_q, the TD t_end rule, the CAP settings, and the
   quadrature rule, with all B4 thresholds unchanged. It is committed before any B4 code is run.
7. **Only then:** the B4 driver and its tests, as in audit §10, P0.2–P0.3.

### 7.3 Open methodological issues (not resolvable from repository evidence)

- **OPEN-1 Directional resolution.**
  - The window lengths limit the spectral resolution to ≈ 4.9 (V3 right), ≈ 10.9 (V3 left) and ≈ 22.4 (V0 right).
    This is arithmetic, 2π/L_w.
  - Mode-converted flux with kₓ below these scales is poorly separated by direction. The operator-leakage term (a) and
    the calibration (b) quantify this. Whether they make criterion (iii) decidable cannot be known before the run, so
    NOT EVALUABLE is a possible outcome.
  - The clean calibration contains only mode-1 broadband flux.
- **OPEN-2 Comparison-time conservatism.** §3.3.4 counts all slab and exterior probability as undetermined. If u_j at
  t_c is not small compared with 10⁻⁴, the comparisons will be NOT EVALUABLE. No less conservative assignment is
  justified by the repository.
- **OPEN-3 Order floor.** The 10⁻⁸ floor is the B2 (ODE-context) convention. B3Y used no floor. Its adequacy for the
  RGF rests on the measured method spread (≤ 3.3×10⁻¹¹), not on a derived bound.
- **OPEN-4 Three grids.** One observed order per quantity cannot establish an asymptotic regime. dx = 1/4096 is not
  part of the pilot (see §1.5).
- **OPEN-5 Synthetic shapes.** The S1–S4 classes are idealized. Physical threshold behaviour, including Fano-type
  lineshapes, is not represented, and the synthetic test validates no physical structure.
- **OPEN-6 Pilot tolerances (§0.1).** These are set from measured rounding levels or existing conventions; they are not
  derived bounds.
- **OPEN-7 Bitwise L3.** Whether `DisorderRealization.numpy` gives shape-independent results bitwise, and whether the
  shifted-grid dx equality holds on the execution platform, must be shown by test. The design-time arithmetic of
  §3.3.1 is for this platform only.

- **OPEN-8 Free-propagation baseline.** B_free (§3.2) is exact only for free propagation on the infinite line. The
  box walls, the CAPs, the discrete DST representation and the truncation of ψ₀ at x = 0 are neglected. The pilot
  reports |W_back,grid − W_back| but cannot bound the remaining baseline error.
- **OPEN-9 Directional bounds.**
  - The assumptions (L-ii), no evanescent content in the window, and (L-iii), that the leakage estimate bounds the
    misassignment, cannot be established from the repository.
  - (L-iv) is false in general, in both directions: returned probability neither lower- nor upper-bounds the
    CAP-induced change of T or R, because returned amplitudes interfere coherently with existing ones.
  - D_max is diagnostic only, and criterion (iii) has no FAIL outcome.
  - F_up is a conditional upper bound on *returned probability* only (conditional on SD), never on ΔT or ΔR.
- **OPEN-10 BELOW FLOOR.** Two consecutive changes below 10⁻⁸ do not establish an error bound (§1.4). ADMISSIBLE-F
  points are carried into the amendment with that label.
- **OPEN-11 dx = 1/2048 selection.** No change to a finer grid is available. Selection of 1/2048 (ADMISSIBLE-R) rests on
  the observed-order assumption from three grids.
- **OPEN-12 TD robustness magnitudes.** No tolerance for Δt_s, ΔP_slab(t), ΔT or ΔR is justified by the repository, so
  their quantitative comparison is NOT EVALUABLE. Only the ROBUST / NOT ROBUST rule of §2.4 decides.
- **OPEN-13 No rigorous absolute CAP error bound.**
  - The pilot cannot establish a rigorous upper bound on the CAP-induced error of the production T and R.
  - The variant comparisons (§3.3.4) bound only *differences* between CAP variants; a common-mode return would be
    invisible to them.
  - The directional estimator (§3.3.5) bounds only returned *probability*, which does not bound the amplitude-level
    effect (L-iv).
  - A rigorous absolute bound would need a reference free of CAP return (for example an exact transparent boundary).
    No such reference exists in the repository.
  - The B4 budget item "CAP ≈ 10⁻⁴" therefore remains an assumption supported, at best, by consistent evidence. No
    replacement bound or tolerance is proposed.

All existing limitations and NOT VALIDATED statements of the audit and of B1–B3Y remain in force.

### 7.4 Statement of limits

- **A successful pilot does not constitute a B4 result.** It contains no comparison of stationary and time-dependent
  transmission or reflection for any disordered realization.
- **It does not independently validate Phase 2 transport.** Under PHASE2_METHOD_SELECTION.md §3, Phase 2 transport
  results are considered independently validated only if B1–B4 all pass. B1 (5/6, B1-F), B3 and B3X statuses are
  unchanged, and B4 has not been performed.
- **It establishes no continuum convergence, no physical transport statement, and no statement about Anderson
  localization.**
- **It establishes no rigorous bound on the CAP-induced error of T or R.** The CAP pilot bounds T/R differences
  between CAP variants and, conditionally, the returned probability. Neither is an absolute bound on the CAP
  contribution to the B4 error budget (§3.3.6, OPEN-13).
