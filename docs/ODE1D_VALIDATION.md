# Validation of the independent 1D ODE scattering reference (`ode1d.py`)

Phase 2, implementation stage 1. Date: 2026-11-05. Base commit: `d51697f`.
Scope: the 1D reference required by benchmark B2 of
[PHASE2_METHOD_SELECTION.md](PHASE2_METHOD_SELECTION.md). **This is not a
validation of the stationary (RGF) solver. That solver does not exist
yet, and nothing here tests it.** No simulation data were changed and no
disorder sweeps were run.

Reproduce: `python validate_ode1d.py` (≈ 10 s) writes
`results/ode1d_validation/{sharp.csv, smooth.csv, rtol_convergence.csv, summary.json, run_metadata.json}`.
Tests: `python -m pytest -q tests/test_ode1d.py tests/test_analytic.py`.

## 1. Mathematical formulation
Stationary 1D Schrödinger equation (ħ = m = 1, energy ε > 0):

  −½ φ″(x) + V(x) φ(x) = ε φ(x),  i.e.  φ″ = 2(V(x) − ε) φ.

V is real and vanishes in the clean leads x ≤ x_L and x ≥ x_R.
k = √(2ε).

In B2 the 1D problem is the y-uniform barrier of the waveguide in the
incident transverse mode 1. There T^{L→R}_{1←1}(E) = T_1D(ε) with
ε = E − E₁^⊥, and T_1D is computed here.

## 2. Boundary-condition convention
Incidence from the left. Amplitudes are referred to the planes x_L
(incident and reflected waves) and x_R (transmitted wave):

  x ≤ x_L: φ = e^{ik(x−x_L)} + r e^{−ik(x−x_L)},  x ≥ x_R: φ = t e^{ik(x−x_R)}.

Both leads have V = 0, so the flux-normalized probabilities are T = |t|² and
R = |r|², and R + T = 1 for real V. Moving a reference plane multiplies t or
r by a unimodular phase (`shift_reference`):

  t′ = t·e^{ik[(x_R′ − x_R) + (x_L − x_L′)]},  r′ = r·e^{2ik(x_L − x_L′)}.

T and R are unchanged by such a move. The plane-wave matching reference
(`analytic.rectangular_barrier_matching`) uses the same convention with
x_L = 0 and x_R = w.

## 3. Numerical method
Backward shooting with `scipy.integrate.solve_ivp`. It uses no lattice and
no Green functions, and shares no code with the time-dependent or the
planned RGF solver.
1. At x_R, start from the pure outgoing wave φ = 1, φ′ = ik (a provisional
   t′ = 1).
2. Integrate the first-order system (φ, φ′) from x_R to x_L.
   * Method: DOP853 (explicit Runge–Kutta of order 8 with adaptive steps),
     with rtol = atol = 10⁻¹². These defaults come from the specification
     and were not tuned per case.
   * Sharp potentials: the integration is split at the discontinuities.
     φ and φ′ are carried continuously across each breakpoint. Inside each
     segment V is evaluated at x clipped into the segment interior, so the
     adaptive stepper never samples V on the wrong side of a step.
3. At x_L, decompose into plane waves: A = (φ + φ′/(ik))/2 (incident) and
   B = (φ − φ′/(ik))/2 (reflected). Then t = 1/A and r = B/A. This rescaling
   is exact because the ODE is linear.

**Smooth barrier.** The callable `Potential(cfg).barrier_numpy` of the
existing code is used, with cfg = `waveguide_config(barrier_height=V_b)`:
centre 2.4, width 0.1, tanh edge width 0.01. This is bitwise the same
function the RGF must sample (test
`test_smooth_barrier_uses_the_potential_callable`). Its tails never vanish
exactly. `clean_interval` widens [x₁, x₂] in steps of 0.01 until
|V| ≤ 10⁻¹² at the interval ends, which gives [2.18, 2.62] for V_b = 150 and
for V_b = 250. The neglected tail beyond these points is below 10⁻¹² in
potential (relative ≤ 7×10⁻¹⁵). Its effect is bounded by the
range-doubling test.

**Threshold ε = V₀.** The ODE is regular there (φ″ = 0 inside a sharp
barrier), so it is evaluated directly. The closed form
`analytic.rectangular_barrier_transmission` is now written as
T = [1 + V₀²w²s²/(2ε)]⁻¹, with s = sin(qw)/(qw), sinh(κw)/(κw) or 1. This
is algebraically identical to the textbook form, has no singular
intermediate expression, and gives the continuous limit
T(V₀) = 1/(1 + V₀w²/2) exactly at ε = V₀. Away from the threshold it agrees
with the textbook form to ≤ 2.2×10⁻¹⁶. The exact-matching reference uses
the linear interior solution C + Dx at ε = V₀.

## 4. Test cases
* **Sharp rectangular barriers:** (V₀, w) ∈ {(150, 0.1), (250, 0.1),
  (250, 0.3)}. Each uses 60 energies uniformly spanning the longitudinal
  energies of the B4 packet window, ε = kₓ²/2 with kₓ ∈ [4.09, 35.91], i.e.
  ε ∈ [8.36, 644.8]. To these are added 5 near-threshold energies
  V₀(1 ± 10⁻³), V₀(1 ± 10⁻⁶) and V₀ exactly. This covers energies below,
  at and above the barrier, and deep tunnelling (w = 0.3: T down to
  9.7×10⁻⁷). In total 195 cases.
* **Smooth barriers (existing Potential):** V_b ∈ {150, 250}, the same 65
  energies each, 130 cases.
* **Convergence study:** rtol = atol ∈ {10⁻⁶, 10⁻⁸, 10⁻¹⁰, 10⁻¹², 10⁻¹³}
  at ε ∈ {60, 249, 400}, for the sharp barrier (V₀ = 250, w = 0.1) and the
  smooth barrier (V_b = 250).
* **Other tests:** a clean lead is transparent (T = 1, t = e^{ikL}); input
  validation; the smooth barrier tends to the sharp one as its edge width
  shrinks (e = 0.01, 0.003, 0.001; the difference decreases monotonically,
  informative only).

## 5. Criteria and measured results
Criteria (a) and (b) are from PHASE2_METHOD_SELECTION.md B2, unchanged.
The flux and amplitude criteria were added and fixed before the first
measurement. All values are maxima over the cases in §4, at the default
settings.

| Quantity | Criterion | Measured (max) | Pass |
|---|---|---|---|
| (a) sharp: \|T_ODE − T_formula\| | ≤ 1×10⁻⁸ | **1.5×10⁻¹³** | yes |
| sharp: \|t_ODE − t_match\| (complex) | ≤ 1×10⁻⁸ | **6.6×10⁻¹³** | yes |
| sharp: \|r_ODE − r_match\| (complex) | ≤ 1×10⁻⁸ | **2.3×10⁻¹³** | yes |
| flux \|R + T − 1\| (sharp / smooth) | ≤ 1×10⁻⁹ | **8.5×10⁻¹⁴ / 2.3×10⁻¹³** | yes |
| sharp: lead range +0.5 on each side, \|ΔT\| / \|Δt\| after reference shift | – (test: ≤ 10⁻¹⁰ / 10⁻⁹) | 3.9×10⁻¹³ / 2.9×10⁻¹² | yes |
| (b) smooth: rtol 10⁻¹⁰ → 10⁻¹², \|ΔT\| | ≤ 1×10⁻⁹ | **1.4×10⁻¹⁰** | yes |
| (b) smooth: clean distance doubled, \|ΔT\| | ≤ 1×10⁻¹⁰ | **6.9×10⁻¹³** | yes |
| exactly at ε = V₀ vs analytic limit | ≤ 1×10⁻⁸ (test) | ≤ 3×10⁻¹⁶ (spot values) | yes |

Per sharp configuration, max |T_ODE − T_formula| is 4.8×10⁻¹⁴ (150, 0.1),
4.2×10⁻¹⁴ (250, 0.1) and 1.5×10⁻¹³ (250, 0.3). In deep tunnelling
(w = 0.3, ε = 5–60, T < 10⁻³) the relative error is ≤ 10⁻⁶ (test).

The smallest margin is criterion (b)(rtol): 1.4×10⁻¹⁰ against 10⁻⁹,
about 7×. This is still an order of magnitude below the threshold.

## 6. Convergence results (`rtol_convergence.csv`)
Error of the sharp barrier against the closed form, and change of the
smooth barrier relative to rtol = 10⁻¹³:

| ε | rtol | sharp \|T − T_f\| | sharp nfev | smooth \|ΔT\| | smooth nfev |
|---|---|---|---|---|---|
| 60 | 10⁻⁶ | 1.2×10⁻⁸ | 38 | 2.3×10⁻⁸ | 302 |
| 60 | 10⁻⁸ | 2.0×10⁻¹⁰ | 50 | 9.7×10⁻¹¹ | 506 |
| 60 | 10⁻¹⁰ | 1.9×10⁻¹² | 74 | 2.0×10⁻¹¹ | 806 |
| 60 | 10⁻¹² | 1.9×10⁻¹⁴ | 134 | 4.9×10⁻¹⁴ | 1322 |
| 249 | 10⁻⁶ … 10⁻¹³ | ≤ 1.9×10⁻¹⁵ | 26 | 1.2×10⁻⁷ … 9.3×10⁻¹⁴ | 326 … 1766 |
| 400 | 10⁻⁶ | 1.6×10⁻⁸ | 26 | 2.9×10⁻⁷ | 314 |
| 400 | 10⁻⁸ | 2.9×10⁻¹⁰ | 50 | 1.8×10⁻⁹ | 506 |
| 400 | 10⁻¹⁰ | 1.5×10⁻¹² | 74 | 1.2×10⁻¹¹ | 854 |
| 400 | 10⁻¹² | 6.2×10⁻¹⁵ | 122 | 1.3×10⁻¹³ | 1442 |

* The error decreases roughly in proportion to the tolerance (about ×100
  per ×100 in rtol). This is the behaviour expected from an adaptive
  integrator. It is not an order in a fixed step size, so none is quoted.
* At ε = 249, just below V₀ = 250, the sharp-barrier error is at rounding
  level for every tolerance. The interior solution is then nearly linear,
  and DOP853 integrates it almost exactly in a few steps.
* The smooth barrier needs ≈ 10× more function evaluations than the sharp
  one, because of its 0.01-wide tanh edges.
* The flux error follows the same pattern: from ≈ 10⁻⁸ at rtol = 10⁻⁶ to
  ≈ 10⁻¹⁴ at 10⁻¹².

## 7. Remaining limitations
* Only single-channel 1D scattering. It validates the y-uniform barrier
  case of B2. It says nothing about channel mixing, disorder, or the RGF
  implementation.
* **No RGF validation has been performed.** The RGF-vs-ODE comparison of B2
  (|ΔT| ≤ 2×10⁻⁴, measured convergence order) can only be run once
  `stationary.py` exists.
* The smooth barrier has no closed form here. Its accuracy rests on
  self-convergence (tolerance and lead range), flux conservation, and the
  monotone approach to the sharp barrier as the edge shrinks. The
  Landau–Lifshitz smooth-step formula is not used (it is not verified;
  PHASE2_METHOD_SELECTION.md §8).
* The smooth-barrier tails are truncated at |V| ≤ 10⁻¹². The range-doubling
  test bounds the effect at 7×10⁻¹³ for these parameters, but it is not a
  proof for arbitrary potentials.
* Backward shooting can lose relative accuracy for very small T (thick or
  high barriers), because the solution grows exponentially across the
  barrier. It is verified here down to T ≈ 10⁻⁶ (relative error ≤ 10⁻⁶).
  Much smaller T would need a check before use.
* Errors were measured on one machine and one library set
  (`requirements-lock.txt`; the environment is recorded in
  `run_metadata.json`).
