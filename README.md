# Quantum Wave Dynamics in Disordered Potential Landscapes

A computational quantum-physics project on how **static, spatially
correlated disorder** affects the propagation, scattering, transmission and
reflection of a quantum wave packet in two dimensions. It also looks at
whether localization can be detected at the accessible system sizes. Two
independent solvers are used for the time-dependent Schrödinger equation:

* a **classical spectral split-step solver** (DST-based, validated against
  exact solutions and by convergence tests). All physics experiments use it.
* a **physics-informed neural network (PINN)**. It is evaluated on the same
  disordered potentials (identical realizations) to find out where PINNs are
  reliable and where they fail.

The main contribution is the **physical investigation** and an honest
**method assessment**. No new neural-network architecture or new
mathematical result is claimed. See [docs/LITERATURE.md](docs/LITERATURE.md)
for the literature, the research gap, and the line between established
knowledge and this project's own results.

> Project history: an earlier version described the work as "inspired by
> Environment-Assisted Quantum Transport (ENAQT)". The model has no
> environment, dephasing or open-system dynamics, so this framing has been
> removed. An audit of the earlier version, which found that its main
> results were not supported, is in [docs/AUDIT.md](docs/AUDIT.md).

---

## 1. Research question

> **How do the strength W and the spatial correlation length ζ of a
> Gaussian-correlated random potential control the transmission and
> reflection of a quantum wave packet through a finite disordered region
> of a 2D waveguide? Are the observed trends consistent with
> perturbative (Born-approximation) transport theory, and can
> localization be distinguished from finite-size effects at the accessible
> system sizes?**

Secondary (methodological) question: *how accurate is a PINN solution of
the same disordered 2D TDSE, compared with a validated classical solver,
and what limits it?*

## 2. Scientific motivation

Disorder is everywhere: impurities in conductors, roughness in waveguides,
random optical (speckle) potentials for cold atoms. Waves that scatter many
times off a random potential interfere. Depending on the disorder strength,
its correlation length and the dimensionality, transport can be
**ballistic**, **diffusive**, or suppressed by **Anderson localization**
(Anderson 1958; Abrahams et al. 1979; Lee & Ramakrishnan 1985). Correlated
disorder matters in practice: speckle potentials have a finite grain size,
and the ratio of grain size to wavelength (kζ) controls how strongly and in
which directions a wave is scattered (Kuhn et al. 2007).

Simulating the time-dependent dynamics directly lets us watch these
processes and measure transmission statistics over many disorder
realizations. The theory used for comparison is established. The
numerical experiment tests where its perturbative predictions hold for a
specific, finite geometry.

## 3. Physical model

Dimensionless units with ħ = m = 1. The wavefunction ψ(x, y, t) obeys

$$ i\,\partial_t\psi = -\tfrac12\left(\partial_x^2+\partial_y^2\right)\psi + V(x,y)\,\psi $$

on a rectangle [0, L_x] × [0, L_y] with **hard walls** (ψ = 0 on the
boundary). Lengths are in units of the waveguide width.

Two geometries are used (`config.py`):

| Geometry | Used for | Domain | Initial state | Boundaries |
|---|---|---|---|---|
| **Open waveguide** (`waveguide_config`) | transport experiments | L_y = 1, L_x = 2.9 + L_d | Gaussian packet in transverse mode 1, k_x0 = 20, σ_x = 0.2 | hard walls in y; complex absorbing potentials (CAPs) at both x-ends |
| **Closed box** (`closed_box_config`) | spreading / localization diagnostics | L × L, L ∈ {1, 2} | 2D Gaussian, k_x0 = 20, σ = 0.08 | hard walls |
| **PINN problem** (`pinn_config`) | PINN vs reference | 1 × 1, t ∈ [0, 0.02] | 2D Gaussian at (0.35, 0.5), k_x0 = 15, σ = 0.07 | hard walls |

Waveguide numbers: k = (k_x0² + π²)^{1/2} ≈ 20.25, E₀ = k²/2 ≈ 205,
wavelength ≈ 0.31, **N = 6 open transverse channels**. The disorder occupies
the slab 1.9 < x < 1.9 + L_d (default L_d = 1).

### Real/imaginary form (used by the PINN)

With ψ = u + iv, the equation splits into two coupled real PDEs:

$$ f_u \equiv \partial_t u + \tfrac12\nabla^2 v - V v = 0, \qquad f_v \equiv \partial_t v - \tfrac12\nabla^2 u + V u = 0. $$

### Initial conditions (`initial_condition.py`)

* 2D Gaussian: ψ₀ = (π σ_x σ_y)^{-1/2} exp[−(x−x₀)²/2σ_x² − (y−y₀)²/2σ_y²] e^{i(k_x0 x + k_y0 y)}
* Waveguide mode: ψ₀ = (π σ_x²)^{-1/4} e^{−(x−x₀)²/2σ_x²} e^{i k_x0 x} · (2/L_y)^{1/2} sin(π y / L_y)

Both are normalized to 1 (checked on the box to < 10⁻⁶ by the tests). The
momentum spread is σ_k = 1/(√2 σ_x).

### Conservation law

Without absorbing layers the dynamics is unitary: ∫|ψ|² dx dy = 1 for all t.
With CAPs, the probability removed at each end is recorded exactly, so
**T + R + P_in = 1** holds to rounding error (validation V6).

## 4. Disorder generation (`disorder.py`)

The disorder is a **continuum** random field defined by a seed, the RMS
strength W and the correlation length ζ:

$$ V_{\rm dis}(\mathbf r) = W\sqrt{\tfrac{2}{M}}\sum_{j=1}^{M}\cos(\mathbf q_j\cdot\mathbf r+\varphi_j)\;S(x),\qquad \mathbf q_j\sim\mathcal N(0,\zeta^{-2}\mathbb 1_2),\ \ \varphi_j\sim U[0,2\pi) $$

with M = 512 modes and an optional smooth slab window
S(x) = ½[tanh((x−a)/δ) − tanh((x−b)/δ)], δ = 0.02. This is the
spectral-representation method of Shinozuka & Deodatis (1991). Averaged
over the ensemble, for **any** M:

$$ \langle V_{\rm dis}\rangle = 0,\qquad \langle V_{\rm dis}(\mathbf r)V_{\rm dis}(\mathbf r')\rangle = W^2 e^{-|\mathbf r-\mathbf r'|^2/2\zeta^2}. $$

The power spectrum is S(q) = 2πW²ζ² e^{−q²ζ²/2}. Individual realizations
are **not** re-normalized, so the ensemble statistics stay exactly as
stated.

Properties (all tested in `tests/test_potential.py`):
* the **same seed gives the same physical field at any grid resolution**;
* NumPy (reference solver) and PyTorch (PINN) evaluations **agree
  pointwise** to < 10⁻⁹ relative in float64;
* the field is differentiable, so it can be used inside the PINN residual;
* the ensemble mean, variance and correlation function match the
  definition within sampling error (3000 realizations).

Optional deterministic barrier (`potential.py`): a strip of height V_b and
width w with tanh edges, or a sharp step when the edge width is zero.

### Perturbative benchmark (Born approximation)

For a particle with wave number k in unbounded 2D Gaussian-correlated
disorder, Fermi's golden rule gives (`disorder.born_rates`, with
a = k²ζ²):

$$ \frac{1}{\tau_s}=2\pi W^2\zeta^2 e^{-a}I_0(a),\qquad \frac{1}{\tau_{tr}}=2\pi W^2\zeta^2e^{-a}\left[I_0(a)-I_1(a)\right],\qquad \ell_{s,tr}=k\,\tau_{s,tr}. $$

These hold only for weak disorder (kℓ ≫ 1) and ignore the waveguide's
discrete channels, so they serve as **order-of-magnitude benchmarks**. A
simple diffusive estimate of the slab transmission is
T ≈ ℓ/(L_d + ℓ) with ℓ = (π/2)ℓ_tr. This is an order-of-magnitude form,
not a fit.

## 5. Classical reference solver (`reference_solver.py`)

Strang split-step propagation (Strang 1968; Feit, Fleck & Steiger 1982)
using a **2D type-I discrete sine transform**. Sine modes sin(nπx/L) are
the eigenfunctions of the Dirichlet Laplacian, so the kinetic step is exact
for the sine-series representation and ψ = 0 on the walls holds by
construction. The scheme is second order in Δt and spectrally accurate in
space. A CAP W(x) = W₀ s² (s = depth into a 0.6-wide layer, W₀ = 300)
absorbs outgoing waves. Several disorder realizations are propagated as a
batch.

The reference solver is **not assumed to be exact**. Its accuracy is
measured in Section 7.

## 6. PINN (`model.py`, `pinn_core.py`, `train.py`)

* MLP with 4 hidden layers of 64 tanh units. Inputs are normalized to
  [0, 1] (x/L_x, y/L_y, t/T), followed by an optional Fourier-feature
  embedding (6 linearly spaced frequencies from 1 to 6 cycles per unit;
  Tancik et al. 2020).
* **Hard constraints**: u = u₀ + A·B(x,y)·τ·N_u and v = v₀ + A·B·τ·N_v,
  with B = 16 x̂(1−x̂)ŷ(1−ŷ) and τ = t/T. The initial condition holds
  exactly. The Dirichlet condition holds up to |ψ₀| on the walls
  (< 4×10⁻⁵).
* Loss: mean squared **scaled** residual (T/A)·(f_u, f_v), plus a
  normalization penalty (∫|ψ|² − 1)² at random times.
* Collocation: a fresh uniform batch of 4000 points **every Adam step**, or
  residual-adaptive sampling (RAD; Wu et al. 2023). Optional causal time
  weighting (Wang, Sankaran & Perdikaris 2024).
* Optimizers: Adam (3000 steps, lr 10⁻³ halved every 1000 steps, PDE-weight
  warm-up), then L-BFGS on a fixed set of 16000 points.
* The potential is the same `Potential` object (same disorder realization)
  that the reference solver uses.

RESULTS_PINN_METHOD_NOTE

## 7. Numerical validation (`validate_solver.py` → `results/validation/`)

All numbers below come from `results/validation/*.csv`
(figure: `figures/validation_solver.png`).

| Check | Setting | Result |
|---|---|---|
| V1 exact eigenmode (constant V = 50) | 63², Δt = 10⁻⁴, t ≤ 0.1 | max \|error\| = 4.6×10⁻¹⁴ (exact to rounding, as expected) |
| V2 analytic free Gaussian | 127² and 255², Δt = 10⁻⁴ and 2×10⁻⁵, t ≤ 0.01 | relative L2 error ≤ 9.4×10⁻¹⁰ |
| V3 Δt convergence, PINN problem with W = E₀ | 127², error vs Δt = 1.25×10⁻⁶ | observed order 2.0; **Δt = 10⁻⁵ (production): 5.3×10⁻⁷** |
| V3 Δt convergence, waveguide with W = 200 | dx = 1/64, t = 1 | order 2; **Δt = 2×10⁻⁴ (production): \|ΔT\| = 2.0×10⁻⁵, \|ΔR\| = 2.3×10⁻⁵** |
| V4 spatial convergence, waveguide W = 200, ζ = 0.02 and 0.05 | same continuum disorder on dx = 1/32, 1/64, 1/128 | **dx = 1/64 (production): \|ΔT\| ≤ 1.4×10⁻⁵**; dx = 1/32: ≤ 2×10⁻³ |
| V4 spatial convergence, PINN problem W = E₀ | 31² … 255² | **127² (production): 1.9×10⁻⁶** relative L2 |
| V5 CAP spurious return (clean waveguide) | vs incident k_x0 | 1.0×10⁻⁴ (k = 20), 1.6×10⁻⁴ (25), 5×10⁻⁴ (30), 1.3×10⁻³ (16), **1.3×10⁻² (12), 3.1×10⁻² (10)** |
| V6 conservation | closed box, W = E₀ / open waveguide | \|N(t) − 1\| ≤ 1.5×10⁻¹² / bookkeeping error ≤ 1.6×10⁻¹³ |
| V7 number of disorder modes M | ⟨T⟩ over 24 realizations, W = 100 | M = 64: 0.488 ± 0.016; 256: 0.515 ± 0.018; 512: 0.518 ± 0.015; 2048: 0.501 ± 0.014 (no significant dependence; M = 512 used) |

**Conclusion.** At the production settings the discretization error of the
transport observables (~10⁻⁵) is about three orders of magnitude below the
statistical uncertainty from finite disorder ensembles (SE ~10⁻²). The one
non-negligible systematic error is the **CAP for slow components**: at
k_x0 = 12 about 1.3 % of the probability is spuriously returned, so
reflection at the lowest momentum in the momentum sweep is biased upward by
about 0.01. The clean-waveguide values in `transport_vs_momentum.png` show
this bias directly.

**Contrast with the original configuration** (docs/AUDIT.md): with the
original settings (k₀ = 40, T = 1, sharp barrier, 128², Δt = 2.5×10⁻⁴) the
reference wavefunction at t = T changed by 51 % when Δt was divided by
four, and by 177 % between grids. It was not a valid reference.

## 8. Experimental methodology

**Transport experiments** (`experiment_disorder.py`, reference solver only):

| Sweep | Varied | Fixed | Realizations per point |
|---|---|---|---|
| `strength` | W ∈ {0, 25, 50, 75, 100, 150, 200, 300} (W/E₀ = 0 … 1.46) | ζ = 0.05, L_d = 1, k_x0 = 20 | 48 |
| `corr_length` | ζ ∈ {0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2} (kζ = 0.4 … 4) | W = 100 | 48 |
| `momentum` | k_x0 ∈ {12, 16, 20, 25, 30} | W = 100, ζ = 0.05 | 32 |
| `length` | L_d ∈ {0.25, 0.5, 1, 2, 3} | W ∈ {100, 200}, ζ = 0.05 | 24 |
| `barrier` | barrier height V_b ∈ {150, 250}, W ∈ {0, 50, 100, 200} | width 0.1, ζ = 0.05 | 24 |
| `convergence` | number of realizations (running mean) | W = 100, ζ = 0.05 | 192 |

* **Observables** (`transport.py`) at the end time t_end (1, or longer for
  slow packets and long slabs): transmission T = absorbed right + P(x > b),
  reflection R = absorbed left + P(x < a), and the probability still inside
  the slab P_in = 1 − T − R. P_in(t_end) is always reported. When it is
  small, T and R are the asymptotic scattering probabilities, averaged
  over the packet's energy spread.
* **Statistics** (`stats_utils.py`): mean, sample standard deviation,
  standard error, Student-t and percentile-bootstrap 95 % confidence
  intervals, median, and the distribution over realizations. Realization
  *i* always uses disorder seed *i*, so different parameter points share
  the same random draws (rescaled in W and ζ). This is a paired design,
  which reduces the noise of differences between neighbouring points. As a
  consequence, the points of one sweep are not statistically independent
  of each other.
* **Number of realizations**: chosen from the `convergence` sweep (SE
  ∝ 1/√n, see `figures/realization_convergence.png`) and the CPU budget.
  48 realizations give SE(⟨T⟩) ≈ 0.01–0.02, well above the
  discretization error (~10⁻⁵).
* **Raw data**: `results/transport/<sweep>/<point>.npz` stores T(t), R(t),
  P_in(t) for every realization, the seeds and the full configuration (JSON).
  `summary.csv` holds the statistics.

**Localization diagnostics** (`experiment_localization.py`): closed boxes
of size L = 1 and L = 2 with disorder everywhere, W ∈ {0, 100, 200, 300},
16 realizations. Mean-square displacement and participation ratio vs time.
The late-time PR/area is compared between box sizes (finite-size test).

**PINN study** (`experiment_pinn.py`): 18 runs, each compared with the
reference solution of the identical configuration (`evaluate.py`). Metrics:
relative L2 error of ψ and of |ψ|², norm drift, centre-of-mass error, error
in P(x > 0.5), training time and reference-solver time.

## 9. Results

RESULTS_SECTION

## 10. Limitations

LIMITATIONS_SECTION

## 11. Reproducing everything

```bash
pip install -r requirements.txt          # numpy, scipy, matplotlib, torch, pytest (+ plotly for visualize.py)

python -m pytest -q tests                # unit/validation tests (~1 min)
python validate_solver.py                # solver validation -> results/validation (~15 min, 4 cores)
python experiment_disorder.py --sweep all          # transport sweeps -> results/transport (~2-3 h, 4 cores)
python experiment_disorder.py --sweep strength --n-real 48   # a single sweep
python experiment_localization.py --n-real 16      # closed-box diagnostics -> results/localization
python experiment_pinn.py                # 18 PINN runs + evaluation -> results/pinn (~4-5 h CPU)
python experiment_pinn.py --only clean_fourier_uniform_s0    # one run
python make_figures.py                   # all figures + figures/CAPTIONS.md from saved data
python visualize.py results/pinn/W56_fourier_uniform_s0      # optional interactive 3D animation
```

Every result directory contains the configuration (JSON) and the seeds
used. Re-running a sweep reuses existing `.npz` files; delete them to
recompute. Timings are for a 4-core CPU without a GPU. The numbers in this
README were produced with Python 3.13, NumPy 2.x, SciPy 1.18, PyTorch 2.x
(CPU).

## 12. Repository structure

```text
config.py                Config dataclass + presets (pinn / waveguide / closed box / legacy)
disorder.py              continuum Gaussian-correlated disorder (NumPy + PyTorch), Born rates
potential.py             V = barrier + disorder, one definition for both solvers
initial_condition.py     Gaussian and waveguide-mode packets (NumPy + PyTorch)
analytic.py              exact solutions for validation
reference_solver.py      DST-I Strang split-step solver, batched, with CAPs
observables.py           norm, region probabilities, MSD, participation ratio
transport.py             T, R, P_in for the open waveguide
stats_utils.py           mean/SE/t-CI/bootstrap-CI, running mean
embeddings.py, model.py  PINN (Fourier features, hard IC/BC constraints)
pinn_core.py             residual, sampling (uniform / RAD), causal weighting, losses
train.py, evaluate.py    PINN training and evaluation against the reference
validate_solver.py       solver validation study
experiment_disorder.py   transport sweeps
experiment_localization.py  closed-box spreading / finite-size test
experiment_pinn.py       PINN performance study
make_figures.py          all figures + captions
visualize.py             interactive Plotly animation (optional)
tests/                   pytest suite
results/                 raw data (CSV/NPZ/JSON) of every experiment
figures/                 figures + CAPTIONS.md
docs/AUDIT.md            audit of the original version
docs/LITERATURE.md       literature review, research gap, verified references
```

## 13. Debugging history (engineering record)

The first PINN version silently collapsed to ψ ≈ 0. Four compounding
issues were found and fixed: (1) IC points too sparse around the packet;
(2) spectral bias (fixed with Fourier features); (3) PDE-loss blow-up from
high Fourier frequencies (narrower frequency range, gradient clipping,
warm-up); (4) an **unnormalized initial condition** that contradicted the
norm = 1 target (fixed by normalizing ψ₀ and imposing the IC and BC as hard
constraints). The audit ([docs/AUDIT.md](docs/AUDIT.md)) then showed that
the resulting model still had relative L2 error ≈ 1.3, because the
original problem (k₀ = 40, T = 1, ~40 wall reflections, ~127 temporal
phase cycles) was far beyond what the network could represent. A fifth
issue, overfitting of a fixed collocation set, is described in Section 9.

## 14. References

The full, verified list is in [docs/LITERATURE.md](docs/LITERATURE.md#references).

## License

MIT (see `LICENSE`).
