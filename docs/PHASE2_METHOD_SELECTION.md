# Phase 2 preparation: reproducibility closure, method selection, validation and statistics plan

Date: 2026-11-05. Starting point: `origin/main` at `3c42f7a` (merge of PR #2,
Phase 1). Main scientific question of Phase 2:

> How does the spatial correlation length ζ of a Gaussian-correlated random
> potential affect the transmission of a quantum wave packet through a
> finite disordered section of a 2D hard-wall waveguide?

Nothing in this document changes the Hamiltonian, the physical model, the
committed baseline results or the conclusions in the README. No new
parameter sweeps were run. No independent solver was implemented (see
§2.6). Anderson localization is not claimed. No new physical discovery is
claimed.

Contents
1. Task A: PINN thread count and environment compatibility (implemented)
2. Task B: independent fixed-energy validation method (selection), with
   §2.7 the mathematical specification (contract for the solver)
3. Task C: minimal validation benchmark with pre-defined criteria
4. Task D: paired statistical analysis of the correlation-length sweep (performed)
5. Implementation plan
6. Required tests and validation criteria (checklist)
7. Remaining scientific limitations
8. References

---

## 1. Task A: PINN cache thread count and environment compatibility

### 1.1 Problem
Phase 1 showed that PINN training depends on the PyTorch intra-op thread
count (2 threads: final rel. L2 0.48850; 4 threads: 0.48987, max |Δw| = 0.011),
but `torch_threads` was only recorded as metadata. A model trained with 2
threads could have been silently reused for a 4-thread request.

### 1.2 Changes (implemented and tested)
* `experiment_pinn.pinn_spec(cfg, torch_threads)`: the spec now contains
  `execution.torch_threads`, so the thread count is part of the cache
  fingerprint.
* `run_one(..., threads=None)` applies `--threads` **before** the cache
  spec is built and validated. The spec uses the *effective* count
  `torch.get_num_threads()`, so a run without `--threads` is also pinned to
  the count actually used. A guard checks that the count did not change
  during training.
* `provenance.py` (schema 2): an **environment fingerprint** for each
  experiment kind, built from the library versions that determine the
  floating-point results, the Python version and the CPU architecture:

  | kind | fields |
  |---|---|
  | transport, localization, validation | python, numpy, scipy, machine |
  | pinn | python, numpy, scipy, torch, machine |

  `check_cache(..., strict_env=True)` refuses a cache from a different
  environment. By default it warns, like a source-code change.
* New flags on all experiment drivers: `--strict-env`, and `--strict`
  (= `--strict-code` + `--strict-env`). **Runs whose results enter a report
  must use `--strict`.**
* Schema-1 (Phase 1) records have no environment fingerprint and no
  `machine` field. The check derives both from the versions and the
  platform string that were recorded, so Phase 1 results stay verifiable.
* Existing results were not modified. The two Phase 1 PINN reproduction
  runs have no thread count in their spec, so they can no longer match any
  request. This is intended (a test asserts it), and they stay on disk as
  records.

### 1.3 Environment review: what is and is not in the fingerprint
| Candidate | Decision | Evidence / reason |
|---|---|---|
| NumPy, SciPy versions | **included** | they implement all reference-solver arithmetic (FFT, exp, sums) |
| PyTorch version (PINN only) | **included** | kernels and reduction order differ between versions; the reference solver does not use it (a test checks that a torch change does not affect transport caches) |
| Python version | **included** | floating-point formatting and `math` functions; cheap to include |
| CPU architecture | **included** | x86_64 vs ARM changes SIMD paths and libm |
| FFT worker / CPU count (reference solver) | **excluded** | measured: wavefunctions are **bit-identical** for 1, 2 and 4 FFT workers (new regression test) |
| PyTorch thread count | **in the PINN spec** (not the environment) | measured to change training results (Phase 1) |
| BLAS/LAPACK vendor, CPU micro-architecture (AVX level), OS kernel | **recorded only** (`platform` string) | cannot be read portably and reliably; listed as a residual risk (§7) |
| GPU / device | **excluded, recorded** | all results are CPU. A future GPU run should get its own spec field |

Verified outcomes: `python -m pytest -q tests` → **78 passed, 1 xfailed**
(the xfail is the pre-existing, documented PINN 20 % accuracy target).
There are 17 new tests, including
`test_pinn_cache_trained_with_2_threads_not_reused_with_4`,
`test_pinn_thread_count_applied_before_cache_validation`,
`test_pinn_default_threads_use_effective_count`,
`test_environment_change_warns_or_fails_in_strict_mode`,
`test_reference_solver_bit_identical_across_fft_workers`, and the 6 tests in
`tests/test_paired.py`.

Observation: the code fingerprint hashes whole driver files. Adding the
`--strict-env` option to `experiment_disorder.py` therefore changed the
transport code fingerprint, so the Phase 1 results now get a "source code
changed" warning, and `--strict-code` refuses them. This is conservative by
design, but it means that **non-numerical edits to a driver also invalidate
strict reuse**. See §7.

---

## 2. Task B: independent fixed-energy validation method

### 2.1 The target problem (stationary form)
Same Hamiltonian as the time-dependent (TD) solver, ħ = m = 1:

  H = −½(∂²ₓ + ∂²ᵧ) + V(x, y),  0 < y < L_y,  ψ(x, 0) = ψ(x, L_y) = 0,

with V = V_dis·S(x) (+ optional barrier). V is negligible outside
a − 0.2 < x < b + 0.2: the window tail is ≈ 2e^{−2d/δ} ≤ 4×10⁻⁹ of W at
d = 0.2 with δ = 0.02. Outside that region the guide is clean, which
defines the asymptotic **leads**. Transverse modes are
χₙ(y) = √(2/L_y) sin(nπy/L_y) with energies Eₙ^⊥ = (nπ/L_y)²/2. At total
energy E the open channels are those with E > Eₙ^⊥, with
kₙ = √(2(E − Eₙ^⊥)). For the baseline (k_x0 = 20) there are N = 6 open
channels. The quantity to compute is the total transmission for incidence
from the left in mode 1, T₁(E) = Σₘ T^{L→R}_{m←1}(E) (sum over the open
outgoing channels of the right lead), with R₁(E) = Σₘ T^{L→L}_{m←1}(E) and
T₁ + R₁ = 1 (unitarity). The precise definitions, normalizations and
conventions are in §2.7, which is the contract for the solver.

In the TD solver the CAPs replace the semi-infinite leads. A stationary
method treats the leads exactly (self-energies or mode matching), so it
tests the CAP approximation as well as the propagation.

### 2.2 Candidates compared

| Criterion | (G) Direct Green's function / sparse solve of the whole region | (RGF) Recursive Green's function, slice by slice | (MM) Mode matching / coupled-channel stationary scattering | (TMM) Transfer matrix |
|---|---|---|---|---|
| Same Hamiltonian and Dirichlet walls | yes (any discretization) | yes | yes (transverse sine basis = exact Dirichlet eigenfunctions) | yes |
| Leads / open channels | exact lead self-energies Σ_L, Σ_R; transmission via Fisher–Lee | same (Σ only at the two ends) | exact: plane waves e^{±ikₙx} per channel | exact at the ends |
| Evanescent (closed) channels | included automatically, stable | included, **stable** (only bounded quantities are propagated) | must be kept for convergence; plain shooting is unstable; needs a log-derivative/R-matrix or S-matrix propagator | **unstable** without stabilization: closed channels grow like e^{κₙL}. For n = 20, κ ≈ 59 gives e^{59} ≈ 10²⁵ over L_d = 1, which destroys all precision. Needs QR or S-matrix composition, which is effectively RGF |
| Cost per energy (L_d = 1, dx = 1/512 ≈ 614 slices, N_y = 63 modes) | sparse LU of ≈ 4×10⁴ unknowns with dense 63×63 blocks: about the same as RGF; more memory | ≈ 614 × O(N_y³) ≈ 10⁹ flop ≈ 0.1–1 s; memory ≈ 40 MB | similar to RGF when a stable propagator is used | ≈ RGF when stabilized |
| Numerical stability | excellent | excellent (standard for disordered wires) | good only with a stable propagator | poor unless stabilized |
| Uses exactly the same disorder realization | yes: evaluate `Potential(cfg, disorder_from_config(cfg, seed))` at the grid | yes, same | yes, via V_nm(x) = ∫χₙVχₘ dy | yes |
| Independence from the split-step TD solver | high: energy domain, linear solve, exact leads, no CAP, different x-discretization | high, same | high | high |
| Implementation risk | low (scipy.sparse) but memory scales with the whole system | low–moderate; textbook algorithm | moderate (propagator) | high (stabilization) |
| Literature status for disordered wires | standard (e.g. Kwant) | standard since Lee & Fisher 1981, MacKinnon 1985 | standard for waveguides (Szafer & Stone 1989; Usuki et al. 1995) | standard, with stabilization (Kramer & MacKinnon 1993) |

### 2.3 Recommendation: mode-space recursive Green's function (primary)

**Primary method: RGF in a hybrid basis.** The transverse direction uses
the sine basis χₙ(y) (n = 1…N_y, the same exact Dirichlet eigenfunctions
as the DST solver). The propagation direction x uses second-order finite
differences, and both ends get exact discrete-lead self-energies.
Transmission follows from the discrete Fisher–Lee relation (§2.7.5), with
t_RL = i Γ_R^{1/2} G^r_{N1} Γ_L^{1/2} restricted to the open channels, where
Γ_α = diag(vₙ/dx) contains the *discrete* lead group velocities, and from
the Caroli trace formula T = Tr[Γ_R G^r_{N1} Γ_L G^a_{1N}] (§2.7.6). (An
earlier draft wrote "t = i√v_R G_{N1} √v_L". That omits the 1/dx factor of
the lattice broadening and does not specify the velocities or the phase
convention. It is superseded by §2.7.)

Reasons:
1. **Exact Hamiltonian and boundary conditions.** The sine basis satisfies
   the Dirichlet walls exactly and gives exact transverse energies, so the
   channel thresholds Eₙ^⊥ are not shifted. A pure finite-difference
   transverse grid would shift them by ≈ (nπ dy)²/12, about 0.7 % for
   n = 6 at dy = 1/64. That matters for near-threshold channels.
2. **Exact treatment of the open and closed channels** through the
   mode-diagonal lead self-energies Σₙ = −t_x λₙ, t_x = 1/(2dx²), where λₙ
   is the retarded lead root (§2.7.2–2.7.3): λₙ = e^{ikₙdx} with
   cos(kₙdx) = 1 − (E − Eₙ^⊥)dx² for open channels, and the real root with
   |λₙ| < 1 for closed channels. Leads and scattering region share the same
   discretization, so a clean guide transmits exactly (benchmark B1).
3. **Stable** for the long evanescent decay that kills a plain transfer
   matrix.
4. **Cheap**: one realization and one energy cost ≲ 1 s, so the ≈ 1300
   energies needed for energy averaging (§3, B4) take ≈ 20 min.
5. **Same realization**: V(x_j, y) is evaluated from the identical
   `Potential` object. V_nm(x_j) is obtained by sine-transforming
   V(x_j, y_l) on a y-grid finer than the TD grid, so the stationary result
   samples the continuum potential rather than copying the TD solver's
   collocation.
6. **Independent of the split-step solver** in every step that could hide
   an error: energy domain instead of time domain, linear solve instead of
   exponential splitting, finite differences in x instead of spectral, exact
   leads instead of CAPs. What is shared on purpose is only the physical
   definition (H, V, transverse basis).

**Internal cross-check (not an independent method):** a direct sparse solve
(G) of the same discretized system must agree with the RGF to 10⁻¹⁰. This
catches coding errors in the recursion, but it does not test the physics or
the discretization.

Rejected as primary:
* **Transfer matrix**: unstable with the closed channels that are needed
  for convergence at small ζ (§2.4). The stabilized forms add complexity
  and gain nothing over RGF.
* **Continuum mode matching / log-derivative propagator**: a good
  method, but it needs a carefully implemented stable propagator (Johnson
  1973), which is more new code to validate than RGF.
* **Whole-system sparse Green's function**: kept as the internal coding
  cross-check.

Optional, later: if the package `kwant` can be installed, a pure
tight-binding Kwant calculation of the same system is a third,
externally written implementation. Its transverse discretization differs
(finite differences), so it must agree only within that discretization
error. This is not a dependency.

### 2.4 Discretization requirements (estimated; to be verified by B3)
* **Longitudinal dispersion error.** The FD lead dispersion is
  E = Eₙ^⊥ + (1 − cos k dx)/dx² instead of k²/2, a relative error in k of
  about (k dx)²/24. That is ≈ 6×10⁻⁵ at dx = 1/512 for k ≈ 20 (35 at the
  top of the packet's spectrum). Plan: dx = 1/512, check against 1/1024,
  and Richardson-extrapolate.
* **Number of transverse modes N_y.** The disorder couples modes n and m
  through wave numbers |q_y| ≈ |n − m|π/L_y up to a few × 1/ζ. For
  ζ = 0.02 that means q ≈ 150, i.e. coupling to n ≈ 50. Plan: N_y = 63 (as
  in the TD grid), checked against 127; ζ = 0.02 is the stringent case.
* **Scattering-region buffer:** 0.2 on each side of the slab (window tail
  ≤ 4×10⁻⁹ W), checked against 0.4.

### 2.5 How transmission is obtained independently
For each energy E: build the lead self-energies, run the RGF recursion
through the slices, and compute t_{m1} from G_{N1}. T₁(E), R₁(E) and the
full matrices come out directly, with no time propagation, no CAP, and no
wave packet. The wave-packet comparison (§3, B4) is a separate step that
weights T₁(E) with the TD packet's known momentum distribution (§3, B4).

### 2.6 Feasibility verdict on implementing the full solver now
The method and the benchmark are well defined, and the cost is small.
I did **not** implement it in this task, because doing it safely requires:
(i) settling the analytic reference formulas against the original sources.
The rectangular-barrier formula is now verified (§2.7.10). The smooth-step
formula in Landau & Lifshitz could not be checked against the book text and
is not used (B2).
(ii) Fixing the mathematical conventions (§2.7) and writing the benchmark
criteria down and committing them **before** the code exists. This
document does that.
Implementing the solver is the first step of Phase 2 (§5).

### 2.7 Mathematical specification of the stationary problem (contract for the solver)

This section is the contract that `stationary.py` must satisfy. Its formulas
are implemented, without any Hamiltonian assembly or recursion, in
`smatrix_spec.py`, and verified on toy lattices in
`tests/test_smatrix_spec.py` (§2.7.10). The derivation below is done for
this specific discretization; formulas from other discretizations
(continuum leads, other normalizations) are not used.

Units: ħ = m = 1. Indexing: every amplitude matrix is indexed
**[outgoing, incoming]**.

#### 2.7.1 Discretization, basis and matrix dimensions
* Slices x_j = x_ref + j·dx. **Device** slices j = 1…N; **left lead**
  j ≤ 0; **right lead** j ≥ N + 1. The leads are clean (V = 0). The device
  includes the buffers around the slab (§2.4).
* On each slice the wavefunction is expanded in the transverse Dirichlet
  modes, ψ(x_j, y) = Σₙ c_{j,n} χₙ(y) with χₙ(y) = √(2/L_y) sin(nπy/L_y),
  n = 1…N_y. All N_y modes (open **and** closed) are kept in the device
  basis. c_j ∈ ℂ^{N_y}.
* Hamiltonian (second-order differences in x, exact transverse operator),
  with t_x = 1/(2dx²):

  (Hc)_j = −t_x c_{j−1} + (2t_x I + E^⊥ + V_j) c_j − t_x c_{j+1},

  E^⊥ = diag(Eₙ^⊥), Eₙ^⊥ = (nπ/L_y)²/2. V_j ∈ ℝ^{N_y×N_y} is defined by
  (V_j)_{nm} = ∫ χₙ V(x_j, y) χₘ dy, computed with the DST-I quadrature
  y_l = l·L_y/(N_q + 1), l = 1…N_q, N_q ≥ 2N_y:
  V_j = Xᵀ diag(V(x_j, y_l)) X · L_y/(N_q + 1), with X_{ln} = χₙ(y_l). This
  form is **real symmetric by construction**. The discrete sines are exactly
  orthogonal (Σ_l χₙ(y_l)χₘ(y_l) L_y/(N_q+1) = δ_{nm} for n, m ≤ N_q), so a
  y-independent potential gives exactly V_j = V(x_j) I.
* Inner product Σ_j c_j†c_j (no dx factor). With this choice G^r has units
  of 1/energy, Γ has units of energy, and the amplitudes t, r are
  dimensionless. Including a dx factor would rescale G and Γ inversely and
  leave t, r unchanged.
* Dimensions: device Hamiltonian H_D ∈ ℝ^{NN_y×NN_y}, block tridiagonal,
  **real symmetric**. Lead self-energies Σ_L, Σ_R ∈ ℂ^{N_y×N_y} act on the
  blocks (1,1) and (N,N). Γ_L, Γ_R ∈ ℝ^{N_y×N_y}, diagonal, with rank
  N_open. Green-function blocks G_{ab} ∈ ℂ^{N_y×N_y}. Scattering matrix
  S ∈ ℂ^{2N_open×2N_open}.

#### 2.7.2 Lead solutions, dispersion and retarded root
In a clean lead the channels decouple. For channel n,
−t_x c_{j−1} + (2t_x + Eₙ^⊥ − E) c_j − t_x c_{j+1} = 0. With c_j = λ^j:

  λ + λ⁻¹ = 2cₙ,  cₙ = 1 − (E − Eₙ^⊥)/(2t_x) = 1 − (E − Eₙ^⊥) dx².

The two roots satisfy λ₊λ₋ = 1.
* **Open (propagating)** channel: |cₙ| < 1, i.e. 0 < E − Eₙ^⊥ < 2/dx².
  Then λ = e^{±ikₙdx} with kₙdx = arccos cₙ ∈ (0, π). **Discrete dispersion:**
  E = Eₙ^⊥ + (1 − cos kₙdx)/dx². **Discrete group velocity:**
  vₙ = dE/dk = sin(kₙdx)/dx (> 0), which tends to √(2(E − Eₙ^⊥)) as dx → 0.
* **Closed (evanescent)** channel: |cₙ| > 1; both roots are real. For
  E < Eₙ^⊥ (cₙ > 1), λ = cₙ − √(cₙ² − 1) ∈ (0, 1). For
  E − Eₙ^⊥ > 2/dx² (cₙ < −1, a lattice artifact; it never occurs here
  because 2/dx² ≈ 5×10⁵ ≫ E), λ = cₙ + √(cₙ² − 1) ∈ (−1, 0).
* **Retarded root convention:** λₙ is the root with |λₙ| < 1 for E → E + i0⁺.
  For an open channel this is λₙ = e^{+ikₙdx} with kₙ > 0, the wave that
  moves away from the device with positive group velocity. For a closed
  channel it is the root that decays away from the device. One rule, used
  for both leads: in the left lead the outgoing/decaying solution is
  λₙ^{−(j−1)}, in the right lead λₙ^{j−N}.
* At a threshold, |cₙ| = 1 (E = Eₙ^⊥), the velocity vanishes and the
  channel count changes. Such energies are excluded (`lead_modes` raises an
  error).

#### 2.7.3 Surface Green function, self-energy and broadening
For channel n of a semi-infinite clean lead the surface Green function
satisfies the Dyson equation gₙ = [E − (2t_x + Eₙ^⊥) − t_x² gₙ]⁻¹. Its two
solutions are −λ/t_x and −1/(t_xλ). The retarded one is

  g^r_n = −λₙ / t_x   (|λₙ| < 1; Im g^r ≤ 0).

The lead–device coupling is τ_L = H_{0,1} = −t_x I and τ_R = H_{N+1,N} = −t_x I. Hence

  Σ^r_α = τ_α† g^r_α τ_α = t_x² g^r_α = −t_x Λ,  Λ = diag(λₙ)  (identical leads: Σ_L = Σ_R),

  Γ_α = i(Σ^r_α − Σ^r_α†) = 2 t_x diag(Im λₙ).

* Open channel: Γₙ = 2t_x sin(kₙdx) = sin(kₙdx)/dx² = **vₙ / dx** (discrete group velocity).
* Closed channel: λₙ is real, so Σₙ is real and **Γₙ = 0 exactly**. Evanescent
  modes stay in the device Green-function basis, where they mediate coupling
  between open channels, but they carry no asymptotic flux and do not
  appear in S.
* Flux: on the lattice the probability current of c_j = a λ^j e_n between
  neighbouring slices is J = 2t_x Im(c_j* c_{j+1}) = |a|² 2t_x sin(kₙdx)
  = |a|² Γₙ. A unit-flux wave therefore has amplitude |a| = Γₙ^{−1/2}. The
  velocities in the normalization must be the **discrete** vₙ.
  `tests/test_smatrix_spec.py::test_continuum_velocities_break_unitarity`
  shows that continuum velocities violate unitarity at O((k dx)²).

#### 2.7.4 Device Green function
G^r(E) = [E·I − H_D − Σ^r_L ⊕ Σ^r_R]⁻¹ (Σ_L on block (1,1), Σ_R on block
(N,N)), G^a = (G^r)†. Required blocks: G^r_{11}, G^r_{N1}, G^r_{1N},
G^r_{NN}. The RGF recursion computes exactly these. The direct sparse solve
(internal cross-check) computes them from the same matrix.

#### 2.7.5 Phase convention and channel-resolved Fisher–Lee amplitudes
**Convention P (reference planes at the first and last device slices).**
For incidence in channel n from the left with amplitude a:

  left lead (j ≤ 1): c_j = a λₙ^{j−1} eₙ + Σₘ bₘ λₘ^{−(j−1)} eₘ,  right lead (j ≥ N): c_j = Σₘ dₘ λₘ^{j−N} eₘ.

Here the lead solution's formal continuation onto the adjacent device slice
defines the reference amplitude, so c_1 = a eₙ + b and c_N = d. Eliminating
the leads gives (E − H_D − Σ_L − Σ_R)c = s with source s_1 = iΓ_L a eₙ, hence
c = G^r s. With unit-flux normalization (amplitudes scaled by Γ^{1/2}), the
flux-normalized amplitudes (open channels only) are

  t_RL = i Γ_R^{1/2} G^r_{N1} Γ_L^{1/2},  t_LR = i Γ_L^{1/2} G^r_{1N} Γ_R^{1/2},
  r_L = i Γ_L^{1/2} G^r_{11} Γ_L^{1/2} − I,  r_R = i Γ_R^{1/2} G^r_{NN} Γ_R^{1/2} − I,

where Γ^{1/2} = diag(Γₙ^{1/2}) over the open channels, and the rows and
columns of G are restricted to the open channels. Notation:
* t_RL = "to R from L", the transmission for incidence from the left;
* t_LR = transmission for incidence from the right;
* channel-resolved probabilities T^{L→R}_{m←n} = |(t_RL)_{mn}|² = Γ_{R,m} Γ_{L,n} |(G^r_{N1})_{mn}|²
  and T^{R→L}_{n←m} = |(t_LR)_{nm}|².

The −I in r is the incident wave's own contribution at the reference plane.

**Phase convention dependence.** Moving the reference planes, e.g. to the
lead surface slices 0 and N+1, multiplies each amplitude by unimodular
channel phases (λ factors, |λₙ| = 1 for open channels): S → D_out S D_in with
diagonal unitary D. Probabilities |S_{ij}|² and unitarity are invariant.
S = Sᵀ is preserved if each channel's incoming and outgoing waves are
referred to the same plane (D_out = D_in); otherwise it is not (test
`test_phase_convention_invariance`). Convention P satisfies this, so the
reciprocity test of §2.7.8 applies to it as stated. The same convention is
used throughout this document and the solver.

#### 2.7.6 Caroli (trace) formula
  T(E) = Tr[Γ_R G^r_{N1} Γ_L G^a_{1N}],  G^a_{1N} = (G^r_{N1})†,

which equals Tr(t_RL† t_RL) = Σ_{m,n} T^{L→R}_{m←n}. It is independent of
the channel phase convention. The total transmission for incidence in
channel 1 alone is T₁(E) = Σₘ T^{L→R}_{m←1} = [t_RL† t_RL]_{11}.

#### 2.7.7 Scattering matrix, block ordering and unitarity
Incoming channel order: [left open channels (ascending n), right open
channels (ascending n)]. Outgoing channel order: the same. Then

  S = [[ r_L , t_LR ],
       [ t_RL, r_R  ]]   (rows: outgoing L; R; columns: incoming L; R).

Incidence from left channel n into right channel m sits at S[N_open + m, n].
**Unitarity** (flux conservation, open channels only): S†S = SS† = I_{2N_open}.
It implies, for every incoming left channel n,
Σₘ T^{L→R}_{m←n} + Σₘ |(r_L)_{mn}|² = 1, and Tr(t_RL†t_RL) = Tr(t_LR†t_LR)
(equal total transmission in both directions). Unitarity relies on the
discrete velocities in Γ and on Γ_closed = 0.

#### 2.7.8 Reciprocity (time-reversal symmetry, no magnetic field)
H_D is real symmetric (real V, real χₙ, real hopping), and Σ_L, Σ_R are
diagonal and therefore complex symmetric. So M = E − H_D − Σ is complex
symmetric (Mᵀ = M), and **G^r = (G^r)ᵀ** exactly: the inverse of a complex
symmetric matrix is complex symmetric. With convention P:

  **S = Sᵀ**, i.e. r_L = r_Lᵀ, r_R = r_Rᵀ, t_LR = t_RLᵀ.

At the level of probabilities this compares **opposite** propagation
directions:

  T^{L→R}_{m←n} = T^{R→L}_{n←m}  (in: left n, out: right m  ⟷  in: right m, out: left n).

The same-direction relation T^{L→R}_{m←n} = T^{L→R}_{n←m} does **not**
follow from time-reversal symmetry. It holds only in special cases, e.g. a
y-uniform potential (t_RL diagonal), or with an additional spatial symmetry.
A generic toy device violates it, while S = Sᵀ holds to rounding
(`test_same_direction_symmetry_is_not_a_valid_criterion`). The earlier B3
criterion "|T_{m←n} − T_{n←m}| ≤ 10⁻⁸" was this incorrect same-direction
form and has been replaced (§3, B3). Reciprocity applies to the open
channels of S. Closed channels are not part of S.

Because the RGF recursion does not impose Mᵀ = M explicitly, the residual
‖S − Sᵀ‖ measures accumulated rounding and coding errors. Typical coding
errors it catches: a transposed coupling block, an asymmetric V quadrature,
or swapped index order.

#### 2.7.9 Relation to the literature
The trace formula is due to Caroli et al. (1971), in the general Γ form of
Meir & Wingreen (1992). The conductance–transmission relation is Fisher & Lee
(1981, continuum derivation). The equivalence of Green-function and
mode-matching amplitudes for lattice leads is treated in Datta (1995) and
Khomyakov et al. (2005). The symmetric S for time-reversal-invariant systems
without spin is the orthogonal symmetry class (Beenakker 1997). The discrete
formulas above were **derived for this discretization and verified
numerically** (§2.7.10). They are not copied from a source with a different
normalization.

#### 2.7.10 Numerical verification of this specification (toy lattices, no physics)
`tests/test_smatrix_spec.py` (19 tests, all passing):
* **Independent check of Fisher–Lee under convention P, phases included:**
  a 1D chain computed with a Green-function-free transfer-matrix wave
  function agrees with t_RL and r_L to ≤ 10⁻¹⁰ (complex values) at three
  energies. The toy is non-trivial: |t|² = 0.16, 0.56, 0.60.
* Lattice continuum limit vs the analytic rectangular-barrier transmission.
  The analytic formula is itself checked against exact plane-wave matching
  to 10⁻¹⁰. The test covers dx = 1/256…1/1024. An additional ad-hoc run
  (not part of the suite) extended this to 1/2048. Errors from 1/256 to
  1/2048 are 7.9×10⁻³ → 5.0×10⁻⁴ (ε = 150) and 4.2×10⁻³ → 2.3×10⁻⁴
  (ε = 320). The **observed orders are
  irregular** (1.97, 0.00, 2.01 and 2.55, −0.25, 1.86), because the sharp
  barrier edges fall differently on each grid. This is why convergence
  orders are measured, not assumed (B2).
* Random multichannel toys (5 modes, 3 open, real symmetric V):
  ‖S†S − I‖ < 10⁻¹², ‖S − Sᵀ‖ < 10⁻¹², r_L, r_R symmetric, t_LR = t_RLᵀ,
  Caroli = Tr(t†t), and equal total transmission in both directions.
* Same-direction symmetry is violated (> 10⁻³) by a generic device:
  counterexample to the old criterion.
* A y-uniform potential gives no channel mixing; a clean guide gives
  |t_nn| = 1, r = 0.
* Continuum velocities in Γ violate unitarity (> 10⁻³ at dx = 0.05).
* Phase-convention test (D S D vs D₁ S D₂), block ordering and index
  mapping, retarded root selection (E + i·10⁻⁹), discrete dispersion,
  Γ = v/dx, Γ_closed = 0, rejection of threshold energies.

---

## 3. Task C: minimal validation benchmark (criteria fixed in advance)

All thresholds below are fixed now, before any stationary-solver result
exists. If a criterion fails, the cause is investigated and reported. The
threshold is not changed and parameters are not tuned to pass.

Notation follows §2.7: S = [[r_L, t_LR], [t_RL, r_R]] over the open channels,
T^{L→R}_{m←n} = |(t_RL)_{mn}|², R_n = Σₘ |(r_L)_{mn}|², and the residuals
U = max|S†S − I| and Q = max|S − Sᵀ|. "Both incidence directions" means
that the left and right incident columns of S are both computed and checked.
**Status of each threshold:** *kept* means unchanged from the previous
version of this document; *new* means it was added in this audit, before
any solver exists; *provisional* means not fully justified a priori (the
reason is given), so a failure must be investigated and reported, never
fixed by changing the threshold.

### B1. Clean, matched waveguide (analytic: t_RL diagonal, |t_nn| = 1, r = 0)
* Energies: 50 values spanning the 1st to 6th channel windows, avoiding
  |E − Eₙ^⊥| < 10⁻³E (thresholds are singular, §2.7.2).
* Pass, for **every** open channel n and both incidence directions:
  * |T^{L→R}_{n←n} − 1| ≤ 10⁻¹²: unit transmission in the incident channel (*kept*; it was stated for n = 1 only, now for all n)
  * Σ_{m≠n} T^{L→R}_{m←n} ≤ 10⁻¹²: no spurious channel mixing (*kept*)
  * R_n ≤ 10⁻¹²: zero total reflection (*kept*)
  * U ≤ 10⁻¹⁰: unitarity of the complete S (*kept*, stated for the full matrix)
  * phase check under convention P, |(t_RL)_{nn} − λₙ^{N−1}| ≤ 10⁻¹⁰: free propagation from slice 1 to slice N (*new*, provisional)
* Provisional because the round-off of the RGF recursion grows with the
  number of slices (≈ 600–1200); about 10⁻¹³ is expected, but this is not
  derived. The 10⁻¹² thresholds are kept. If one is exceeded while staying
  ≤ 10⁻¹⁰, the conditioning is to be analysed and reported.
* TD side (already known): the clean-guide TD packet has T = 0.99989 at
  k_x0 = 20 (V5/V6). The deviation from 1 is the CAP spurious return
  (≈ 10⁻⁴), a known TD error.

### B2. Smooth rectangular barrier, uniform in y
**Same potential in both calculations.** The barrier is the existing
`Potential` barrier (tanh edges, edge width e = 0.01, width 0.1,
V_b ∈ {150, 250}). The 1D ODE evaluates exactly the callable
`Potential(cfg).barrier_numpy(x)` at its integration points. The RGF
samples the same function at the slice positions x_j. Test: the sampled
values are bitwise identical to `barrier_numpy(x_j)`. Because V does not
depend on y, V_j = V_b(x_j)·I exactly (§2.7.1), the channels decouple, and
T^{L→R}_{1←1}(E) = T_1D(E − E₁^⊥) for −½φ″ + V_b(x)φ = εφ.

**Independent 1D reference (ODE), validated first:**
* (a) Sharp rectangular barrier, integrated piecewise (the integration is
  split at the discontinuities), against the closed form
  `analytic.rectangular_barrier_transmission`. That formula is verified
  against exact plane-wave matching to 10⁻¹⁰
  (`test_rectangular_barrier_formula_vs_matching`).
  Pass: |T_ODE − T_formula| ≤ 10⁻⁸ (*new*).
* (b) Smooth barrier: ODE self-convergence. Pass: changing rtol from 10⁻¹⁰
  to 10⁻¹² changes T by ≤ 10⁻⁹, and doubling the clean integration
  range changes it by ≤ 10⁻¹⁰ (*new*).
* The Landau–Lifshitz smooth-step formula is **not used**. Its exact
  parameter convention and problem number could not be verified against
  the text (§2.7.9, §9). The fallback (a) + (b) is the validation.

**RGF vs ODE** (60 energies across the packet spectrum):
* |ΔT| ≤ 2×10⁻⁴ at dx = 1/512 (*kept*, binding).
* Convergence order: **measured** from dx = 1/256, 1/512, 1/1024 and
  reported for every energy.
* The previous criteria "order between 1.8 and 2.2" and
  "Richardson-extrapolated |ΔT| ≤ 2×10⁻⁵" are **provisional**. Second-order
  convergence holds only in the asymptotic regime, and the tanh edge
  (e = 0.01) is resolved by only ≈ 5 points at dx = 1/512. The toy
  sharp-barrier study showed irregular orders (§2.7.10). If the measured
  order is outside [1.8, 2.2], Richardson extrapolation is not applied and
  only the binding 2×10⁻⁴ criterion is evaluated. The outcome is reported
  either way.
* Off-diagonal T^{L→R}_{m≠1←1} ≤ 10⁻¹² (*kept*); U ≤ 10⁻¹⁰ (*kept*);
  Q ≤ 10⁻⁸ (*new*).
* TD vs energy-averaged ODE (exactly as in B4): |ΔT| ≤ 1×10⁻³ for both barrier heights (*kept*).

### B3. One fixed disorder realization: RGF consistency, unitarity, reciprocity, convergence
Realizations: seed 0 at W = 100, L_d = 1, ζ ∈ {0.02, 0.05, 0.20}, plus seed 0
at W = 200, ζ = 0.05 (kℓ_tr ≈ 2.6). 5 test energies each.

* **Same physical realization:** V(x_j, y_l) is computed by
  `Potential(cfg, disorder_from_config(cfg, seed)).numpy` at the RGF
  quadrature points. Test: bitwise equality with a direct call at the same
  points, and the seed recorded in the run metadata.
* **RGF vs direct sparse Green function** (same discretized M, base
  resolution dx = 1/512, N_y = 63): max|S_RGF − S_sparse| ≤ 10⁻¹⁰
  (*kept*, now on the full S).
* **Unitarity, complete S, both incidence directions:** U ≤ 10⁻¹⁰
  (*kept*; it replaces the weaker |T₁ + R₁ − 1| ≤ 10⁻¹⁰, which is implied).
* **Reciprocity (replaces the incorrect same-direction criterion):** with
  the full S from both incidence directions,
  * Q = max|S − Sᵀ| ≤ 10⁻⁸;
  * block-wise: max|r_L − r_Lᵀ| ≤ 10⁻⁸, max|r_R − r_Rᵀ| ≤ 10⁻⁸,
    max|t_LR − t_RLᵀ| ≤ 10⁻⁸;
  * probability form: max_{m,n} |T^{L→R}_{m←n} − T^{R→L}_{n←m}| ≤ 10⁻⁸.

  The numerical value 10⁻⁸ is *kept* from the previous criterion. Only
  the compared quantities are corrected. The same-direction quantity
  max|T^{L→R}_{m←n} − T^{L→R}_{n←m}| is **reported only**; it is not a
  criterion.
* **Closed-channel (transverse) convergence:** |ΔT₁(E)| ≤ 2×10⁻⁴ for
  N_y 63 → 127 (*kept*). The number of closed channels retained,
  N_y − N_open(E), is reported. For ζ = 0.02, the most demanding case
  (§2.4), if the criterion fails at 127, N_y = 255 is run and reported. The
  threshold is not changed.
* Other refinements, each |ΔT₁(E)| ≤ 2×10⁻⁴ (*kept*): dx 1/512 → 1/1024;
  y-quadrature N_q ×2; buffer 0.2 → 0.4.

### B4. Fixed realization: stationary vs time-dependent, with the energy width accounted for
Same realizations as B3, TD solver at the production settings.

**The energy-width correction is mandatory.** The TD packet is
φ(x) ∝ e^{−(x−x₀)²/2σₓ²}e^{ik₀x} in mode 1, with momentum distribution

  P(kₓ) = (σₓ/√π) exp[−(kₓ − k₀)²σₓ²],  E(kₓ) = (kₓ² + π²/L_y²)/2.

The TD kinetic operator is spectral, so this E(kₓ) mapping is exact on the
TD side. The stationary solver is evaluated at the same energy E, not at
a lattice wave number. For σₓ = 0.2 the standard deviation of kₓ is
1/(√2σₓ) ≈ 3.5, so the packet covers E ≈ 205 ± 70 (1σ). For linear
scattering, with the packet initially entirely in the left lead and moving
right, the transmitted probability at t → ∞ in the problem with ideal
leads is the energy average

  T_pkt = ∫ P(kₓ) T₁(E(kₓ)) dkₓ,  T₁(E) = Σₘ T^{L→R}_{m←1}(E)  (sum over the open channels at that E)

(standard wave-packet scattering theory, e.g. Taylor 1972). The
contributions of different energies do not interfere asymptotically.
Conditions, all checked for the baseline: initial overlap with the slab
≈ ½erfc(3.5) ≈ 4×10⁻⁷; weight of kₓ < 0 ≈ ½erfc(4) ≈ 8×10⁻⁹; weight of
kₓ < 5 (slow, poorly absorbed by the CAP) ≈ 10⁻⁵.

The stationary result is compared as ⟨T⟩_E = ∫P T₁ dkₓ (and ⟨R⟩_E
likewise), **never as T₁(E₀) at the central energy alone**.

* Energy grid: T₁(E) of a single realization fluctuates on the
  correlation-energy scale E_c ~ ħ/τ, where τ is the longest relevant dwell
  time. Measured from the committed strength sweep (∫P_in dt and the
  late-time exponential decay of ⟨P_in(t)⟩): mean dwell 0.078 (W = 100) and
  0.075 (W = 200); late-time decay time 0.17 (W = 100) and 0.19 (W = 200).
  So E_c ≈ 5–6. Plan: uniform in kₓ over k₀ ± 4.5/(√2σₓ) with
  δkₓ = 0.025 (≈ 1270 energies; δE = kₓδkₓ ≤ 0.9 ≈ E_c/6 over the whole
  window).
  * The number of open channels changes across the window (from 1 to about
    8). T₁(E) is continuous at a channel threshold, but its derivative is
    singular (√ cusp). Grid points within 10⁻³E of a threshold are moved off
    it, and the quadrature error there is covered by the halving check.
  * Convergence check: halving δkₓ changes ⟨T⟩_E by ≤ 2×10⁻⁴ (*kept*).
  * Estimated cost ≈ 20 min per realization. (An earlier draft estimate of
    E_c ~ 10, from the mean dwell time alone, was too optimistic.)
* **Remaining probability in the slab.** The TD run (same realization) is
  extended until P_in(t_end) ≤ 10⁻⁴ (*kept*, now a pass condition). The
  remaining probability is undetermined between T and R, so the TD result
  is the interval 𝒯 = [T_TD, T_TD + P_in], and likewise
  ℛ = [R_TD, R_TD + P_in], with T_TD + R_TD + P_in = 1.
* Error budget (estimated from the Phase 1 validation; RGF items to be
  confirmed by B2/B3): TD Δt 2×10⁻⁵, dx 1.4×10⁻⁵, CAP ≈ 10⁻⁴,
  P_in ≤ 10⁻⁴; RGF dx ≈ 10⁻⁴, energy quadrature ≤ 2×10⁻⁴. Total ≲ 5×10⁻⁴.
* **Pass** (*kept* thresholds; P_in handling made explicit), for every B4
  realization:
  P_in(t_end) ≤ 10⁻⁴ **and** dist(⟨T⟩_E, 𝒯) ≤ 1×10⁻³ **and**
  dist(⟨R⟩_E, ℛ) ≤ 1×10⁻³, where dist(x, [a, b]) = 0 if a ≤ x ≤ b and
  otherwise the distance to the nearest endpoint. This is ≈ 2× the error
  budget and ≈ 10× below the ensemble standard error of ⟨T⟩ (≈ 10⁻²).
* B4b (narrow packet, informative, not pass/fail): a TD run with σₓ = 0.8
  (σ_k ≈ 0.88) in a longer domain. It is compared (i) with the
  energy-averaged stationary result using the *same* criterion, and (ii)
  with the monochromatic T₁(E₀), only to **quantify** how far a naive
  monochromatic comparison would be off.

Phase 2 transport results are considered independently validated only if
B1–B4 all pass.

---

## 4. Task D: paired statistical analysis of the correlation-length sweep (post-hoc re-analysis, performed)

### 4.1 Design
Realization i of every ζ point uses disorder seed i. For a given seed,
`make_disorder` uses the same Gaussian draws and phases at every ζ
(q = g/ζ), so the seeds are **paired**. Verified on the data: identical
seed arrays (0–47) at all 7 ζ, and specs that differ only in
`disorder_corr_length` (`paired_analysis.load_sweep` refuses anything
else). Analysis unit: the per-seed difference dᵢ = Σⱼ wⱼ Tᵢ(ζⱼ).

* **Bootstrap:** B = 20 000 resamples of the seed index (with replacement),
  shared by all contrasts (bootstrap RNG seed 20261105). Percentile CIs:
  95 % and Bonferroni 98.33 % for the 3 primary contrasts.
* **Bootstrap p-value (corrected in this audit).** Null-centering procedure:
  each bootstrap mean Δ*_b is centered at the observed mean Δ̂, so
  Δ*_b − Δ̂ approximates the sampling distribution of (Δ̂ − Δ) and hence the
  null distribution of Δ̂ when the true Δ = 0. A replicate counts as at
  least as extreme as the observation if |Δ*_b − Δ̂| ≥ |Δ̂| (two-sided). With
  b such replicates,

    p_boot = (b + 1)/(B + 1).

  This finite-resampling (+1) rule (North, Curtis & Sham 2002; Phipson &
  Smyth 2010; Davison & Hinkley 1997) never reports zero. Its smallest value
  is 1/(B + 1) = 5.0×10⁻⁵. It is a valid, slightly conservative test p-value,
  **not an unbiased estimate of the true p** (North et al.'s 2003 reply;
  Ewens 2003; Broman & Caffo 2003). The cited papers treat Monte Carlo and
  permutation tests; here the rule is applied to the centered-bootstrap test
  as a conservative finite-resolution floor. The previous version reported
  p_boot = 0 when b = 0; that was not a valid p-value.
* **Other p-values:** paired t and Wilcoxon signed-rank. All three kinds are
  Holm-adjusted over the primary family {C1, C2, C3}. The smallest
  Holm-adjusted bootstrap p-value is therefore 3/(B + 1) = 1.5×10⁻⁴.
* **Effect sizes:** d_z = mean(d)/sd(d), with a bootstrap CI, and
  d_av = mean(d)/mean(sd).

Primary contrasts:
* **C1** T(0.05) − T(0.02);
* **C2** T(0.05) − T(0.20);
* **C3** curvature K = [T(0.03) + T(0.07)]/2 − T(0.05).

Secondary (not multiplicity-controlled): the two sides of C3; the bootstrap
distribution of the grid argmin; the vertex of a parabola in ln ζ through
0.03, 0.05, 0.07.

**Status: post-hoc re-analysis.** The contrasts were chosen after the
*unpaired* results had been seen (README §9.3). This analysis re-examines
the same 48 realizations with the correct (paired) dependence structure. It
is **not** a confirmatory test and **not** an independent replication. An
independent confirmation needs new seeds, fixed in advance (§5, step 5).

Input: `results/reproduction/transport/corr_length/` (fingerprinted,
bit-identical to the committed baseline). Command: `python paired_analysis.py`.
Outputs: `results/analysis/paired_corr_length.{json,csv}`,
`figures/paired_corr_length.png`. Re-running after the p-value correction
changed only the p_boot fields. CIs, effect sizes, the argmin distribution
and the vertex are bit-identical, and the figure is unchanged.

### 4.2 Results (n = 48 seeds, W = 100, L_d = 1, k ≈ 20.25)

| Contrast | mean ΔT | paired 95 % CI | Bonferroni 98.3 % CI | unpaired Welch 95 % CI | r (pairs) | d_z [95 % CI] | p_boot (b/B) | Holm p: boot / t / Wilcoxon |
|---|---|---|---|---|---|---|---|---|
| C1 T(.05)−T(.02) | −0.176 | [−0.198, −0.154] | [−0.203, −0.149] | [−0.200, −0.152] | +0.14 | −2.25 [−2.89, −1.86] | 5.0×10⁻⁵ (0/20000) | 1.5×10⁻⁴ / 1×10⁻¹⁹ / 2×10⁻¹⁴ |
| C2 T(.05)−T(.20) | −0.240 | [−0.284, −0.193] | [−0.294, −0.182] | [−0.290, −0.190] | +0.18 | −1.49 [−2.19, −1.07] | 5.0×10⁻⁵ (0/20000) | 1.5×10⁻⁴ / 3×10⁻¹³ / 6×10⁻¹¹ |
| C3 curvature K | +0.046 | [0.023, 0.068] | [0.018, 0.073] | – | – | +0.58 [0.27, 1.00] | 1.0×10⁻⁴ (1/20000) | 1.5×10⁻⁴ / 2×10⁻⁴ / 7×10⁻⁵ |
| C3a T(.03)−T(.05) | +0.051 | [0.029, 0.073] | [0.024, 0.078] | [0.024, 0.078] | +0.31 | +0.64 | 5.0×10⁻⁵ (0/20000) | secondary, unadjusted |
| C3b T(.07)−T(.05) | +0.041 | [0.006, 0.074] | **[−0.002, 0.081]** | [0.006, 0.075] | −0.01 | +0.34 | 1.9×10⁻² (382/20000) | secondary, unadjusted |

A p_boot of 5.0×10⁻⁵ with b = 0 means only that **p ≤ ~5×10⁻⁵ at this
resolution**. Smaller values cannot be resolved with B = 20 000. The
parametric t / Wilcoxon values are much smaller, but they rest on their
own distributional assumptions.

* Grid argmin of ⟨T⟩: ζ = 0.05 in 98.9 % of bootstrap replicates, ζ = 0.07 in 1.1 %.
  This locates the minimum only to the resolution of the 7-point grid.
* Parabola vertex through three grid points in ln ζ: kζ* = 0.95, bootstrap
  95 % CI [0.88, 1.12]. This number depends on the parabolic model and on
  three coarse points. It **does not precisely determine the physical
  minimum**; it only indicates where a minimum lies for that model. For
  comparison, the minimum of the Born transport mean free path is at
  kζ = 0.889, and that of the energy-averaged Born + diffusion estimate of T
  is at kζ = 0.883.

### 4.3 Interpretation (limited to what the analysis supports)
* In this post-hoc paired re-analysis, with Holm correction over three
  primary contrasts, ⟨T⟩ at ζ = 0.05 is lower than at ζ = 0.02 and at
  ζ = 0.20, with large effect sizes, and ζ = 0.05 lies below the mean of its
  two grid neighbours (C3).
* **The ζ = 0.07 side is not established after multiplicity correction.**
  The secondary contrast T(0.07) − T(0.05) has a small effect (d_z = 0.34),
  and its Bonferroni interval [−0.002, 0.081] includes zero. The data
  therefore locate the minimum only between ζ ≈ 0.03 and 0.07
  (kζ ≈ 0.6–1.4). A finer, pre-registered ζ grid with new seeds is needed for
  a sharper statement.
* Pairing helps little here: the per-seed correlation of T across ζ is only
  r ≈ −0.01 to +0.31, because the same draw rescaled in ζ gives a quite
  different scattering landscape. The paired and unpaired intervals are
  therefore similar. The paired analysis is still the correct one.
* These results describe this geometry (N = 6 channels, L_d = 1,
  W/E₀ ≈ 0.49, energy-averaged packet) and this solver. Their numerical
  correctness is what the B1–B4 benchmark is designed to confirm
  independently.

---

## 5. Implementation plan (Phase 2, in order)

1. **Reference formulas.** Done for the rectangular barrier
   (`analytic.rectangular_barrier_transmission`, verified against exact
   plane-wave matching). Remaining: `ode1d.py` (independent 1D integrator,
   piecewise for sharp barriers), tests vs the formula to 10⁻⁸, and
   self-convergence for the smooth barrier. The smooth-step formula is used
   only if verified against the book.
2. **`stationary.py`: mode-space RGF.** It builds H_D (§2.7.1, with
   symmetric V_j), runs the recursion for G_{11}, G_{N1}, G_{1N}, G_{NN},
   and obtains lead modes, Σ, Γ, Fisher–Lee amplitudes, S and the checks
   **only through `smatrix_spec.py`** (§2.7), which is already specified
   and tested. A direct sparse-solve variant of the same M serves as the
   internal cross-check. It uses the `Potential`/`disorder` objects
   unchanged.
3. **`validate_stationary.py`**: runs B1–B4 with `provenance` metadata and
   `--strict`, writes `results/stationary_validation/` and a pass/fail table
   against the criteria in §3. Estimated cost: B1–B3 a few minutes. B4 is
   4 realizations × (≈ 20 min RGF energy grid, ≈ 40 min with the halving
   check, + ≈ 10 min TD with extended t_end) ≈ 3–3.5 h. B4b ≈ 40 min.
4. **Report** in `docs/PHASE2_VALIDATION.md`: every criterion with its
   measured value, including any failure.
5. Only after B1–B4 pass: a **finer ζ grid** around the minimum
   (e.g. ζ ∈ {0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.07}, 48–96 paired
   seeds, `--strict`). Its pre-registered contrasts are the location of
   the minimum (vertex CI) and C3-type curvature. Optionally, a
   stationary-solver ζ sweep at fixed energy to separate the
   energy-averaging effect.
6. Phase-2 science items from README §9.11 (longer slabs, channel-number
   dependence) only after step 5.

## 6. Required tests and validation criteria (checklist)

Implemented now (passing):
- [x] S-matrix specification (`tests/test_smatrix_spec.py`, 19 tests): lead modes (open/closed, retarded root, discrete dispersion, Γ = v/dx, Γ_closed = 0, threshold rejection); 1D Fisher–Lee vs independent transfer matrix (amplitudes and phases); continuum limit vs the rectangular barrier; analytic barrier vs exact matching; multichannel unitarity, S = Sᵀ, symmetric reflection blocks, t_LR = t_RLᵀ, Caroli = Tr t†t; opposite-direction mode reciprocity; same-direction counterexample; y-uniform decoupling; clean guide incl. phase λ^{N−1}; continuum-velocity unitarity violation; phase-convention invariance; block ordering and index mapping
- [x] Bootstrap p-value (+1) correction, finite resolution, Holm floor, determinism on the real data, preserved analysis settings
- [x] PINN cache rejects a different torch thread count; the thread count is applied before validation; the effective count is used by default
- [x] Environment fingerprint; strict refusal; warning by default; schema-1 derivation; torch irrelevant for transport
- [x] Reference solver bit-identical across FFT worker counts
- [x] Paired analysis: CI coverage, pairing benefit on correlated data, curvature/vertex, Holm, refusal of unpaired or inconsistent inputs, loads the reproduced data

Required for the stationary solver (to be written with it):
- [ ] `test_solver_uses_spec_conventions`: the solver's S comes from `smatrix_spec.fisher_lee`/`assemble_S` (no private copy of the formulas)
- [ ] `test_V_blocks_symmetric`: V_j = V_jᵀ exactly (quadrature form of §2.7.1)
- [ ] `test_clean_guide_unit_transmission` (B1 thresholds incl. the phase λ^{N−1}, both incidence directions)
- [ ] `test_unitarity_and_reciprocity`: full S from both incidence directions, U ≤ 10⁻¹⁰, Q = max|S − Sᵀ| ≤ 10⁻⁸, reflection blocks symmetric, t_LR = t_RLᵀ, T^{L→R}_{m←n} = T^{R→L}_{n←m}
- [ ] `test_rgf_equals_direct_sparse` (full S, 10⁻¹⁰)
- [ ] `test_y_uniform_barrier_decouples_channels` (off-diagonal ≤ 10⁻¹²)
- [ ] `test_ode1d_vs_rectangular_formula` (10⁻⁸; `test_ode1d_vs_smooth_step` only if the source is verified)
- [ ] `test_rgf_barrier_vs_ode1d` (2×10⁻⁴ binding; measured order reported)
- [ ] `test_same_realization` (V on the RGF slices equals `Potential.numpy` at the same points, bitwise)
- [ ] `test_energy_average_weights` (P(kₓ) normalized; reproduces T_pkt = 1 for the clean guide; threshold energies moved off the grid)
- [ ] `test_b4_interval_criterion` (dist(⟨T⟩_E, [T_TD, T_TD + P_in]) logic)
- [ ] B3/B4 as slow, marked integration tests (run by `validate_stationary.py`, not by the default suite)

## 7. Remaining scientific limitations

1. **No independent validation yet.** All transport results still come from
   one solver family (split-step + CAP), validated by convergence,
   analytic free/eigenmode tests and conservation. B1–B4 are designed to
   close this gap.
2. **Energy averaging.** Every TD result is an average over a broad band
   (E ≈ 205 ± 70). Statements about "transmission at kζ ≈ 1" refer to the
   central k. The stationary solver will make fixed-energy statements
   possible.
3. **Coarse ζ grid at the minimum.** The minimum is located only to
   kζ ≈ 0.6–1.4. The ζ = 0.07 side is not separated after multiplicity
   correction.
4. **Re-analysis, not replication.** The paired contrasts were chosen
   after the unpaired results were seen. A finer, pre-registered sweep with
   new seeds (e.g. 1000–1095) would be an independent confirmation.
5. **Single geometry.** N = 6 channels, L_d = 1, one disorder strength for
   the ζ sweep. Born-approximation comparisons are order-of-magnitude
   (unbounded 2D, no channel structure).
6. **Finite observation time and CAP.** P_in(t_end) up to 10⁻³; CAP return
   ≈ 10⁻⁴ at k ≈ 20 and 1.3 % at k_x0 = 12.
7. **No localization claim.** Nothing here addresses localization. The
   ζ-dependence concerns ⟨T⟩ in a diffusive/quasi-ballistic regime
   (kℓ_tr ≈ 10–50 for W = 100).
8. **Specification residuals.** The B1 10⁻¹² thresholds and the B2
   order window / Richardson criterion are provisional (§3). The RGF
   round-off behaviour over ≈ 10³ slices is not derived, and the asymptotic
   convergence regime for the tanh-edged barrier is not established.
9. **Statistical resolution.** Bootstrap p-values are limited to
   ≥ 1/(B + 1) = 5×10⁻⁵ (Holm ≥ 1.5×10⁻⁴). The (+1) rule is a valid,
   slightly conservative test p-value, not an unbiased p estimate.
10. **Reproducibility residuals.** BLAS vendor and CPU micro-architecture are
   recorded but not fingerprinted. Bit-identity has been shown on one
   machine only. The code fingerprint is per file, so non-numerical driver
   edits trigger it (conservative, sometimes inconvenient). The PINN study
   (17 of 18 runs) has not been re-run with a recorded, fixed thread
   count.

## 8. References (new in this document; others in docs/LITERATURE.md)

* P. A. Lee, D. S. Fisher, "Anderson localization in two dimensions," *Phys. Rev. Lett.* **47**, 882 (1981). doi:10.1103/PhysRevLett.47.882
* D. S. Fisher, P. A. Lee, "Relation between conductivity and transmission matrix," *Phys. Rev. B* **23**, 6851 (1981). doi:10.1103/PhysRevB.23.6851
* A. MacKinnon, "The calculation of transport properties and density of states of disordered solids," *Z. Phys. B* **59**, 385 (1985). doi:10.1007/BF01328846
* S. Datta, *Electronic Transport in Mesoscopic Systems*, Cambridge University Press (1995).
* A. Szafer, A. D. Stone, "Theory of quantum conduction through a constriction," *Phys. Rev. Lett.* **62**, 300 (1989). doi:10.1103/PhysRevLett.62.300
* Usuki et al., "Numerical analysis of ballistic-electron transport in magnetic fields by using a quantum point contact and a quantum wire," *Phys. Rev. B* **52**, 8244 (1995). doi:10.1103/PhysRevB.52.8244 (author order differs between indexes; check on the APS page before citing)
* B. R. Johnson, "The multichannel log-derivative method for scattering calculations," *J. Comput. Phys.* **13**, 445 (1973). doi:10.1016/0021-9991(73)90049-1
* C. W. Groth, M. Wimmer, A. R. Akhmerov, X. Waintal, "Kwant: a software package for quantum transport," *New J. Phys.* **16**, 063065 (2014). doi:10.1088/1367-2630/16/6/063065
* L. D. Landau, E. M. Lifshitz, *Quantum Mechanics: Non-Relativistic Theory* (Course of Theoretical Physics Vol. 3), §25 and its problems. The smooth-step reflection coefficient is attributed there by citing papers, but they disagree on the problem number and the parametrization. **Not verified against the text; not used.**
* J. R. Taylor, *Scattering Theory: The Quantum Theory of Nonrelativistic Collisions*, Wiley (1972): wave-packet description of scattering.
* B. Kramer, A. MacKinnon, *Rep. Prog. Phys.* **56**, 1469 (1993) (transfer-matrix methods; already in docs/LITERATURE.md).

### 8.1 Added in the specification audit (verified bibliographic data)
* C. Caroli, R. Combescot, P. Nozières, D. Saint-James, "Direct calculation of the tunnelling current," *J. Phys. C: Solid State Phys.* **4**, 916 (1971). Volume and page confirmed through citations in Parts II and III of the same series (doi:10.1088/0022-3719/4/16/025 and /026). **No DOI was verified for Part I**, so none is given.
* Y. Meir, N. S. Wingreen, "Landauer formula for the current through an interacting electron region," *Phys. Rev. Lett.* **68**, 2512 (1992). doi:10.1103/PhysRevLett.68.2512
* P. A. Khomyakov, G. Brocks, V. Karpan, M. Zwierzycki, P. J. Kelly, "Conductance calculations for quantum wires and interfaces: Mode matching and Green's functions," *Phys. Rev. B* **72**, 035450 (2005). doi:10.1103/PhysRevB.72.035450
* C. W. J. Beenakker, "Random-matrix theory of quantum transport," *Rev. Mod. Phys.* **69**, 731 (1997). doi:10.1103/RevModPhys.69.731 (symmetry of S in the orthogonal class; already in docs/LITERATURE.md)
* B. V. North, D. Curtis, P. C. Sham, "A note on the calculation of empirical P values from Monte Carlo procedures," *Am. J. Hum. Genet.* **71**, 439 (2002) (PMID 12111669), and the 2003 exchange in *Am. J. Hum. Genet.* **72**: W. J. Ewens (496–498), K. W. Broman & B. S. Caffo (496), North, Curtis & Sham (reply, 498–499). The reply states that (r+1)/(n+1) gives the correct type-I error rate but is not an unbiased estimate of p.
* B. Phipson, G. K. Smyth, "Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn," *Stat. Appl. Genet. Mol. Biol.* **9**(1), Article 39 (2010). doi:10.2202/1544-6115.1585
* A. C. Davison, D. V. Hinkley, *Bootstrap Methods and their Application*, Cambridge University Press (1997). (Cited for bootstrap tests in general; no chapter or equation is claimed.)

### 8.2 Sources by topic
| Topic | Source(s) | How it is used here |
|---|---|---|
| Fisher–Lee relation | Fisher & Lee 1981 (continuum, conductance = Tr t†t); Datta 1995; Khomyakov et al. 2005 (lattice leads, Green function ↔ mode matching) | the discrete, phase-explicit formulas of §2.7.5 are **derived here** and verified against an independent transfer-matrix calculation |
| Caroli transmission formula | Caroli et al. 1971; Meir & Wingreen 1992 (general Γ form) | §2.7.6, verified as Tr(t†t) on toys |
| Time-reversal reciprocity of S | Beenakker 1997 (orthogonal class, S symmetric); follows algebraically from Mᵀ = M (§2.7.8) | derivation given here; verified on toys |
| Recursive Green function | Lee & Fisher 1981; MacKinnon 1985; Datta 1995; Kwant (Groth et al. 2014) | method for G_{11}, G_{N1}, G_{1N}, G_{NN} |
| Rectangular-barrier transmission | standard textbook result (no specific textbook location is claimed) | `analytic.rectangular_barrier_transmission`, **verified here** against exact plane-wave matching to 10⁻¹⁰ |
| Bootstrap (+1) p-value | North et al. 2002 (+2003 exchange); Phipson & Smyth 2010; Davison & Hinkley 1997 | §4.1 |


---

## 9. Specification audit (mathematical and methodological; no solver implemented)

Audit of commit `1c72630`. Scope: mathematics, conventions, criteria and
statistics only. No change to the Hamiltonian, the physical model, the
committed simulation data or the physical conclusions. No new simulations
were run. `stationary.py` was not started.

### 9.1 Mathematical corrections
1. **Reciprocity (B3 was wrong).** The old criterion
   |T_{m←n} − T_{n←m}| ≤ 10⁻⁸ compared the same propagation direction,
   which time-reversal symmetry does not imply. A generic toy device
   violates it by > 10⁻³ while satisfying the correct relation to
   10⁻¹². It is replaced by S = Sᵀ for the full flux-normalized S (both
   incidence directions): r_L = r_Lᵀ, r_R = r_Rᵀ, t_LR = t_RLᵀ, equivalently
   T^{L→R}_{m←n} = T^{R→L}_{n←m} (§2.7.8). The threshold 10⁻⁸ is unchanged.
2. **Fisher–Lee normalization (§2.3 was underspecified and dimensionally
   incomplete).** The draft "t = i√v_R G_{N1} √v_L" lacked the 1/dx of the
   lattice broadening and did not define the velocities or the phase
   convention. Now: Γ = diag(vₙ/dx) with **discrete** velocities
   vₙ = sin(kₙdx)/dx, t_RL = iΓ_R^{1/2}G^r_{N1}Γ_L^{1/2}, the reflection
   formula with its −I term, and the explicit phase convention P (§2.7.5).
3. **Lead solutions and conventions made explicit:** discrete dispersion,
   open/closed criteria, retarded root |λ| < 1, Σ = −t_xΛ, g^r = −Λ/t_x,
   Γ_closed = 0, matrix dimensions, [outgoing, incoming] indexing, S block
   ordering, unitarity of the complete S (§2.7.1–2.7.7).
4. **Benchmark criteria (§3):** B1 now applies to every open channel and both
   incidence directions, with an added phase check. B2 requires the
   identical barrier callable in the ODE and the RGF, validates the ODE
   against the (now verified) rectangular formula, and measures the
   convergence order instead of assuming it. B3 adds the full-S unitarity,
   the corrected reciprocity, closed-channel reporting and the bitwise
   same-realization check. B4 makes the remaining slab probability part of
   the pass condition (interval criterion) and handles channel thresholds
   inside the energy window.
5. **Thresholds:** all numerical thresholds were kept. Two are labelled
   **provisional** with reasons, and no threshold was relaxed: the B1
   10⁻¹² values (RGF round-off over ≈ 10³ slices is not derived), and the
   B2 order window [1.8, 2.2] with its Richardson criterion (asymptotic
   regime not established; toy sharp-barrier orders were irregular).
   New thresholds added before any solver exists: the B1 phase check
   10⁻¹⁰, the B2 ODE checks 10⁻⁸, 10⁻⁹ and 10⁻¹⁰, and the B2 reciprocity
   check 10⁻⁸.
6. **Bootstrap p-values (§4).** When no replicate was as extreme as the
   observation, the old code reported p = 0. Now p = (b + 1)/(B + 1) with an
   explicit centered-null procedure. Holm adjustment is applied to the
   bootstrap p-values as well. Re-running changed only the p_boot fields:
   for C1, C2 and C3a, 0 → 5.0×10⁻⁵ (b = 0, the resolution limit); for C3,
   5.0×10⁻⁵ → 1.0×10⁻⁴ (b = 1); for C3b, 0.0191 → 0.01915. CIs, effect
   sizes, argmin and vertex are bit-identical, and the figure is unchanged.
   The wording was corrected from "confirmatory" and "pre-specified" to
   **post-hoc re-analysis** (also in README §9.3 and the
   `paired_analysis.py` docstring).

### 9.2 Verified formulas and references
* Verified **numerically in this repository** (toy lattices, no physics;
  `tests/test_smatrix_spec.py`):
  * Fisher–Lee t and r under convention P, against a Green-function-free
    transfer-matrix wave function: complex amplitudes agree to 10⁻¹⁰;
  * Caroli = Tr(t†t);
  * unitarity and S = Sᵀ to 10⁻¹²;
  * Γ = v/dx and the retarded root;
  * the clean-guide phase λ^{N−1};
  * the rectangular-barrier formula against exact matching (10⁻¹⁰);
  * the lattice continuum limit against that formula.
* Bibliographic data verified for: Fisher & Lee 1981; Lee & Fisher 1981;
  MacKinnon 1985; Meir & Wingreen 1992; Khomyakov et al. 2005; Kwant 2014;
  Johnson 1973; Phipson & Smyth 2010; North et al. 2002 and the 2003
  exchange. Caroli et al. 1971 Part I: journal, volume and page verified
  through citing papers; **no DOI verified, none given**.
* **Not verified, not used:** the Landau–Lifshitz smooth-step formula. The
  fallback is the rectangular formula plus ODE self-convergence.
* No page numbers, equation numbers or quotations are claimed for any
  textbook.

### 9.3 Tests executed (actual outcomes)
`python -m pytest -q tests` → **104 passed, 1 xfailed, 4 warnings** (21.9 s).
* The xfail is the pre-existing, documented PINN 20 % accuracy target.
* All 4 warnings come from deliberately degenerate synthetic data in
  `tests/test_paired.py` (a constant series gives an undefined correlation
  coefficient and SciPy precision-loss warnings). None comes from the
  real-data analysis.
* New in this audit: 19 tests in `tests/test_smatrix_spec.py` and 7
  bootstrap tests in `tests/test_paired.py`.
* Re-run of the statistical analysis only: `python paired_analysis.py` on
  `results/reproduction/transport/corr_length/` (fingerprinted Phase 1 data).
  No physical simulation was run.

### 9.4 Remaining uncertainties
* The specification has been verified on small toy matrices: up to 5
  transverse modes, 10–12 slices for the multichannel tests, and up to 110
  slices for the 1D continuum-limit test (213 in the ad-hoc dx = 1/2048 run
  quoted in §2.7.10). The production problem (N_y = 63–127,
  600–1200 slices) may show round-off or conditioning effects that the
  toys cannot reveal. This is why the B1 thresholds are provisional.
* The 1D transfer-matrix check is independent of the Green functions, but
  it is single-channel. Multichannel amplitudes are checked only through
  consistency relations (unitarity, reciprocity, Caroli, decoupling), not
  against an independent multichannel wave-function calculation. Kwant
  (§2.3) would provide such a check.
* The second-order convergence regime for the tanh-edged barrier at
  dx ≈ 1/512 is not established a priori.
* The energy-correlation estimate E_c ≈ 5–6 comes from ensemble-averaged
  dwell times. Single realizations may have narrower resonances, which only
  the δkₓ-halving check can detect.
* The (+1) bootstrap p-value is conservative and resolution-limited. The
  analysis remains post-hoc.

### 9.5 Is the specification ready for the RGF implementation?
**Yes, with two conditions.**
1. The solver must obtain lead modes, Σ, Γ, the amplitudes, S and the checks
   exclusively from `smatrix_spec.py`, which fixes the conventions under test.
2. The provisional thresholds (B1 10⁻¹², B2 order window and Richardson) are
   to be evaluated and reported as provisional. Any failure is investigated,
   not resolved by changing a threshold.

Everything the solver needs is defined: discretization and basis,
dimensions, retarded root, self-energies and broadenings, the phase
convention, the amplitude formulas, S ordering, unitarity, reciprocity and
the B1–B4 pass criteria.
