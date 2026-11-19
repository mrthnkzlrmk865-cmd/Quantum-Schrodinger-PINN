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
- the observed order p = log₂(|T₁(1/512) − T₁(1/1024)| / |T₁(1/1024) − T₁(1/2048)|), reported with the B2 floor rule:
  a change below 10⁻⁸ makes the order "unresolvable";
- a Richardson **ESTIMATE** of the dx = 1/1024 and dx = 1/2048 errors, computed only if p ∈ [1.8, 2.2] (the existing
  provisional window). Otherwise it is "NOT APPLICABLE";
- the method spread of ΔT₁ (H1-type check, as in B3Y).

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
- P_in,res(t) = 1 − T − R (the existing residual definition). Check: |P_in,res − P_slab − ρ − (1 − N(0))| at rounding
  level.
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
- **τ_w = 1.0** is ≥ 1.5× the largest late-time decay time quoted in the protocol (0.19) and in the audit's fits
  (≤ 0.67). It also exceeds the round-trip time from the slab to the right CAP and back for kₓ ≥ 1 (≤ 0.8).

**Pilot output:** for each realization, SETTLED or NOT SETTLED, together with t_s, t_first and U.

**What the pilot recommends to the B4 amendment:** a fixed, predeclared t_end per realization, or a common t_end ≥ max
t_s that also satisfies the trailing window, at which B4 would evaluate T_TD, R_TD and P_in.

### 2.4 Robustness runs (TD only, informative)

- **TD-R1.** Configuration D (strongest disorder) to t_max at dt = 10⁻⁴. Purpose: dt error growth to t = 6, and whether
  the non-monotone P_slab structure is a time-step artifact.
- **TD-R2.** Configuration D at dx = 1/128 (nx = 3.9·128 − 1 = 498). Measured cost ≈ 46 s per 1.5 time units, so
  ≈ 3 min. Purpose: whether the late-time P_slab structure is a spatial-resolution artifact.

Both report only the TD-internal changes of P_slab(t), T(t), R(t) and t_s.

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

**Backward weight:** the analytic weight of kₓ < 0 is W_neg = ½erfc(k_x0σₓ) (1.1×10⁻⁵ at k_x0 = 3; < 10⁻¹² for
k_x0 ≥ 5). It is reported.

**Spurious return**, a time-gated definition:

- R_spur = [A_L(t_end) − A_L(t_g)] + P(x < x₀, t_end).
- t_g is the predeclared time at which the packet's trailing edge (x₀ − 5σₓ moving at k_x0 − 3σ_k) passes x₀. Backward
  components emitted at t = 0 have then left the gap, and are counted separately in A_L(t_g).
- The ungated R_raw = A_L(t_end) + P(x < L_x − 0.6, t_end) is also reported.

**Check of A1:** for each k_x0, |R_spur(n = 4) − R_spur(n = 1)| and |R_spur(n = 6) − R_spur(n = 1)|. Separability
predicts agreement to TD rounding. A difference larger than 10⁻⁸ falsifies A1 in this implementation and must be
reported.

**Output:** R_spur(k_x0, n), a narrowband CAP return curve for the production CAP, and for the CAP variants of §3.3 at
k_x0 ∈ {3, 7, 11}.

### 3.3 C-B: the four B4 realizations with CAP variants

**Same physical realization in a shifted frame.** Extending the left CAP requires moving the domain origin, but the
disorder field depends on absolute x.

- **Approach:** construct the potential as an array, V_s(X, Y) = Potential(cfg, disorder_from_config(cfg, 0)).numpy(X − s, Y).
- Pass it to `SplitStepDST(cfg_s, potentials=V_s)`, which already accepts an array; shift x₀ and all region boundaries
  by s.
- **The disorder parameters are unchanged**; the hash is checked under G-REAL.
- The sampled values differ from production only because the grid point positions differ, which is the same situation
  as any dx change.
- No source change is needed.

**Variants** (gap lengths unchanged except in V3):

| Variant | Left / right CAP width | Strength | Power | Shift s | Domain L_x |
|---|---|---|---|---|---|
| V0 (production) | 0.6 / 0.6 | 300 | 2 | 0 | 3.9 |
| V1 (wide) | 1.2 / 1.2 | 300 | 2 | 0.6 | 5.1 |
| V2 (wide, weaker) | 1.2 / 1.2 | 150 | 2 | 0.6 | 5.1 |
| V3 (long gaps) | 0.6 / 0.6 | 300 | 2 | 1.0 | 5.9 (left gap +1.0, right gap +1.0) |
| V4 (very wide) | 2.4 / 2.4 | 300 | 2 | 1.8 | 7.5 |

**Settings:** dx, dy and dt as production; t_max = 6.0 and the same frames and observables as Part 2.

**How CAP return is separated from genuine reflection by the slab:**

1. **Variant differences.** Genuine scattering by the slab does not depend on the CAP. Differences of T and R at the
   common evaluation time between V0 and V1–V4 estimate the CAP's contribution at production settings. V1, V2 and V4
   differ in CAP gradient; V3 changes only the time of flight. Agreement of V1, V2 and V4 with each other, but not with
   V0, isolates the production-CAP effect.
2. **Direction-resolved gap content.**
   - In the right gap there is no scatterer. Any **left-moving** probability there can come only from the right CAP
     or from the wall behind it.
   - Per frame and per transverse mode n, the right-gap field is multiplied by a predeclared smooth window (a cos²
     taper of 0.05 on each side) and Fourier-transformed in x. The probability with kₓ < 0 is the left-moving content.
   - The same is done for right-moving content in the left gap, after t_dep: the first frame at which the initial
     packet's content in the left gap is below 10⁻⁸ in the clean-guide reference run with the same settings.
   - These are approximate (windowing leakage). They are reported with the leakage measured on the clean-guide C-A
     runs.
3. **Mode-resolved outgoing content.** The probability per mode n in each gap (§2.2) gives the outgoing flux
   distribution over channels. Combined with the C-A curve R_spur(kₓ), it gives a *consistency estimate*
   Σₙ∫wₙ(kₓ)R_spur(kₓ)dkₓ. This is an estimate, not a measurement.

**Predeclared CAP criterion (gate G-CAP).** For every realization, at the candidate B4 t_end from Part 2:

- (i) |T_V0 − T_Vj| ≤ 1×10⁻⁴ and |R_V0 − R_Vj| ≤ 1×10⁻⁴ for j = 1…4; **and**
- (ii) max over the pairs among V1, V2, V4 of |ΔT| and |ΔR| ≤ 1×10⁻⁴, i.e. the enlarged CAPs agree with each other;
  **and**
- (iii) the time-integrated left-moving right-gap content plus right-moving left-gap content (method 2) is ≤ 1×10⁻⁴.

The value 10⁻⁴ is the CAP item of the **existing** B4 error budget (≈ 10⁻⁴), i.e. 10 % of the 10⁻³ tolerance. It is
used as a gate on the budget assumption. It is **not** a new B4 pass threshold.

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
2. **Synthetic cusp test (no physics; exact reference).**
   - Integrand f(kₓ) = s(E) + Σₙ aₙ·√(max(E − Eₙ, 0))/√Eₙ.
   - s is a smooth predeclared function: s(E) = 0.5 + 0.3 sin(E/40).
   - aₙ ∈ {0.05, 0.2} (two predeclared amplitudes), with E = E(kₓ).
   - Reference: ∫P f dkₓ by adaptive quadrature split at every threshold, with tolerance ≤ 10⁻¹³.
   - **Rules compared**, on the coarse and fine grids:
     - **R-a:** the existing rule (moved node, original weight);
     - **R-b:** plain trapezoid with no move (f is finite at the thresholds);
     - **R-c:** a trapezoid split at the thresholds (piecewise grids ending exactly at kₓ(Eₙ)), as the possible amended
       method.
   - **Output:** the true error of each rule and grid; the halving difference; and the ratio (halving difference) /
     (true coarse error).
   - **Reliability criterion.** Halving is "reliable" for a rule if that ratio is ≥ 0.5 for every amplitude, i.e. the
     halving difference does not underestimate the true error by more than 2×.
3. **RGF local refinement windows** (per configuration; M1 only; at the resolution under consideration from Part 1, or
   dx = 1/1024 if undecided).
   - **Windows:** thresholds n ∈ {4, 6, 8, 10}, each ±8 coarse steps around kₓ(Eₙ).
   - **Grids:** h, h/2 and h/8 restricted to the window, with the rules R-a and R-c.
   - **Output:** only *differences* between rules and grids of the local window integral of P·T₁ and P·R₁, and the
     ratio as in item 2.
   - **Absolute window integrals are not reported**, so the pilot cannot be read as partial B4 averages.
   - **Resonance scan:** at the threshold of largest packet weight, n = 6, T₁ and R₁ at spacing h/32 over ±2 coarse
     steps, to look for structure narrower than h/2. This is reported as found or not found, with the minimum width
     resolved.
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
- **per configuration:** the full config dict, the disorder-parameter sha256 over (qx, qy, φ), and the potential-sample
  sha256 for every RGF device and every TD grid (including the shifted-frame variants);
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
- the controls C1 and C2 reproduce the committed values (≤ 10⁻¹²);
- the committed artifacts are unchanged.

Any failure stops the pilot; there is no partial output.

**Amendment rule.** "Amendment" below means a **dated B4 protocol amendment**: a new document committed **before** any
B4 number. It records the pilot evidence and the change, and states that the B4 thresholds 10⁻³, 10⁻⁴ and 2×10⁻⁴ are
unchanged. **No gate may be passed by weakening any existing B4 threshold.**

| Gate | Required evidence | Pass condition | Fail condition | If it fails | Dated amendment required? |
|---|---|---|---|---|---|
| **G-REAL** (same realization) | disorder-parameter hashes and config equality across: the RGF devices (all dx), the TD production grid, the CAP variants (shifted frame), and the B3 records; bitwise RGF slice samples against direct `Potential.numpy` calls (as in B3) | all hashes equal; all bitwise checks true; the region, edge and M identical | any mismatch | **blocking**: no B4 until the cause is found and fixed; any pilot output that depends on the mismatched item is void | no (unless the realization definition itself must change, which would be a new benchmark) |
| **G-RGF** (resolution) | Part 1 rows for all 14 points × 4 configurations × 3 dx × M1–M3 (M4 where specified), plus the existing B3/B3X/B3Y points | (a) M1–M2–M3 agreement on S ≤ 10⁻¹⁰ (or adjudicated as in B3Y) at every point; **and** (b) for the resolution proposed for B4, at every point: the consecutive change to the next finer grid ≤ 2×10⁻⁴ (existing B3/C11 value), p ∈ [1.8, 2.2], and the Richardson ESTIMATE of the error of T₁ and R₁ ≤ 1×10⁻⁴ (the RGF dx item of the existing B4 budget); **and** (c) the N_y and N_q changes ≤ 2×10⁻⁴ at the top-window points; **and** (d) lead-tail guard met | (a) unexplained; or (b) not met at dx = 1/2048 for some configuration; or (c) not met | (a): treat as a solver defect, so B4 is blocked. (b): the B4 budget item cannot be met as designed for that configuration. Options, each needing a dated amendment: a finer dx after a separate dx = 1/4096 review; or a predeclared Richardson-extrapolated B4 quantity with the order gate; or recording that configuration's B4 outcome as discretization-limited. Never a relaxed 10⁻³. (c): a larger N_y or N_q in B4 | **yes, in every case**: the B4 text does not fix the RGF resolution, so even a pass requires the amendment to state the chosen dx, N_y and N_q |
| **G-TD** (stopping rule) | Part 2 frames for the four realizations; TD-R1 and TD-R2 | every realization SETTLED with t_max − t_s ≥ 1.0; ρ(t) ≤ 10⁻¹⁰ for all t; the residual and direct P_in consistent; U(t_end,candidate) ≤ 10⁻⁴; TD-R1/R2 do not move t_s past t_max − 1.0 | any realization NOT SETTLED; or U > 10⁻⁴; or ρ > 10⁻¹⁰; or the robustness runs change the settling | NOT SETTLED: B4 cannot meet its P_in ≤ 10⁻⁴ pass condition with the current setup. Options (dated pilot amendment, then a new pilot run): a larger t_max, or a geometry change (which would make it a different benchmark) | **yes**: the amendment must fix the B4 t_end (or rule), the frame spacing, and the frame at which T_TD, R_TD and P_in are read |
| **G-CAP** (CAP sensitivity) | C-A narrowband curves (modes 1, 4, 6); C-B variants V0–V4 for the four realizations; direction-resolved gap content | criterion §3.3 (i)–(iii) for every realization; the A1 separability check holds (≤ 10⁻⁸) | any of (i)–(iii) > 10⁻⁴; or A1 falsified | the production CAP is not adequate for B4's budget. Options: B4 TD with a CAP variant that passes (i)–(iii) against the even larger variants, or the CAP item entered in the B4 budget as measured. Never a relaxed 10⁻³ | **yes, if the B4 TD settings change**, because the B4 text says "production settings"; **no, if V0 passes** |
| **G-THR** (threshold quadrature) | §4.1 reproduced; the synthetic cusp test; the RGF local windows; the resonance scan | for the quadrature rule proposed for B4: halving reliable (ratio ≥ 0.5) in the synthetic test; and in the local RGF windows, the sum over thresholds of \|rule(h) − rule(h/8)\| ≤ 1×10⁻⁴ (half of the existing 2×10⁻⁴ halving limit, leaving room for the smooth part); and no unresolved structure narrower than h/2 found in the scan | halving unreliable for R-a; or the local error > 10⁻⁴; or narrow structure found | adopt R-c (threshold-split quadrature) or local refinement for B4, re-verified by the same tests; if narrow resonances exist, local adaptive refinement near thresholds with a predeclared rule | **yes**, if the B4 quadrature changes from the protocol's uniform grid with moved nodes; **no**, if R-a passes |

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
| CAP | V5 (mode 1, broadband) | narrowband, mode-resolved return; the variant sensitivity for the four realizations | CAP behaviour for other realizations |
| Threshold quadrature | B2 halving (single channel, no thresholds) | reliability of halving at multichannel thresholds; a quadrature rule for B4 | resonances narrower than the scan resolution |
| Same realization | identical disorder hashes for the RGF and TD paths (audit §7) | run-time assertion across all pilot devices and variants | — |
| B1-F | FAIL (provisional phase criterion), unchanged | rcond and U recorded at the near-threshold points at 1/2048 | phase accuracy on long fine grids (not used by B4) |
| **B4 itself** | **not performed** | **nothing; the pilot produces no B4 comparison** | **everything that B4 is designed to test** |

### 7.2 Prioritized implementation plan (none implemented)

1. **Commit this protocol.** Its sha256 becomes the G0 reference.
2. **Pilot tests** (new `tests/test_b4_pilot.py`): the G0 reproduction of §1.2 and §4.1; the region partition; the
   path guard and labels; the absence of full-window RGF averages; the bitwise equivalence of the subclass in §3.4;
   the trailing-window logic on synthetic P(t) curves (including non-monotone ones and late up-crossings); the
   shifted-frame potential hash; and the synthetic-cusp reference integral.
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

### 7.3 Statement of limits

- **A successful pilot does not constitute a B4 result.** It contains no comparison of stationary and time-dependent
  transmission or reflection for any disordered realization.
- **It does not independently validate Phase 2 transport.** Under PHASE2_METHOD_SELECTION.md §3, Phase 2 transport
  results are considered independently validated only if B1–B4 all pass. B1 (5/6, B1-F), B3 and B3X statuses are
  unchanged, and B4 has not been performed.
- **It establishes no continuum convergence, no physical transport statement, and no statement about Anderson
  localization.**
