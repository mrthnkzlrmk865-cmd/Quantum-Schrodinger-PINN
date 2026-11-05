"""
Central configuration for the project.

Every solver, experiment and test receives an explicit ``Config`` object, so
the PINN and the classical reference solver are always run on the same
physical setup. Experiments make modified copies with ``cfg.replace(...)``
instead of mutating a global object. Each experiment saves
``cfg.to_dict()`` next to its outputs.

Units: dimensionless, hbar = m = 1. The equation solved everywhere is

    i dpsi/dt = -1/2 (d^2/dx^2 + d^2/dy^2) psi + V(x, y) psi

on the rectangle [0, Lx] x [0, Ly] with psi = 0 on the walls (Dirichlet).
An optional complex absorbing potential (CAP) near x = 0 and x = Lx turns the
box into an open scattering geometry (reference solver only; see
reference_solver.py).

Presets
-------
legacy_config()     the original configuration (unit box, k0 = 40, T = 1,
                    sharp barrier). Kept only so the audit is reproducible.
pinn_config()       short-time, small-box problem used for the PINN-vs-
                    reference study (a regime where a PINN can be trained on
                    a CPU).
waveguide_config()  open disordered waveguide used for the transport
                    (transmission / reflection) experiments.
closed_box_config() closed square box with disorder everywhere, used for the
                    spreading / localization diagnostics.
"""

from __future__ import annotations

import dataclasses
import json
import math
from dataclasses import dataclass, field
from typing import Optional, Tuple


@dataclass
class Config:
    # ------------------------------------------------------------------
    # Domain: [0, Lx] x [0, Ly], time in [0, T]
    # ------------------------------------------------------------------
    Lx: float = 1.0
    Ly: float = 1.0
    T: float = 0.02

    # ------------------------------------------------------------------
    # Initial wave packet
    #   "gaussian":       2D Gaussian envelope x plane wave exp(i k.r)
    #   "waveguide_mode": Gaussian in x x plane wave in x, times the
    #                     transverse Dirichlet mode sqrt(2/Ly) sin(n pi y/Ly)
    # ------------------------------------------------------------------
    packet: str = "gaussian"
    x0: float = 0.35
    y0: float = 0.5
    sigma_x: float = 0.07
    sigma_y: float = 0.07
    kx0: float = 15.0
    ky0: float = 0.0
    mode_n: int = 1

    # ------------------------------------------------------------------
    # Deterministic barrier (height 0 = no barrier). The barrier is the
    # strip |x - barrier_center| < barrier_width/2. barrier_edge > 0
    # smooths its edges with tanh profiles of that width; 0 = sharp step.
    # ------------------------------------------------------------------
    barrier_height: float = 0.0
    barrier_center: float = 0.55
    barrier_width: float = 0.06
    barrier_edge: float = 0.0

    # ------------------------------------------------------------------
    # Correlated random disorder (see disorder.py for the definition).
    #   disorder_strength      W: RMS amplitude of V_dis (ensemble value)
    #   disorder_corr_length   zeta: <V(r)V(r')> = W^2 exp(-|r-r'|^2/(2 zeta^2))
    #   disorder_seed          one seed -> one physical field, at any
    #                          resolution, for both solvers
    #   disorder_n_modes       number of Fourier modes M in the series
    #   disorder_region        (x_start, x_end) slab containing disorder,
    #                          None = whole domain
    #   disorder_edge          tanh smoothing width of the slab edges
    # ------------------------------------------------------------------
    disorder_strength: float = 0.0
    disorder_corr_length: float = 0.08
    disorder_seed: int = 0
    disorder_n_modes: int = 512
    disorder_region: Optional[Tuple[float, float]] = None
    disorder_edge: float = 0.02

    # ------------------------------------------------------------------
    # Complex absorbing potential (reference solver only). cap_width = 0
    # disables it. -i * cap_strength * s^cap_power with s in [0, 1] the
    # depth into the absorbing layer, on both x-ends.
    # ------------------------------------------------------------------
    cap_width: float = 0.0
    cap_strength: float = 0.0
    cap_power: int = 2

    # ------------------------------------------------------------------
    # Reference solver discretization. nx/ny = number of interior grid
    # points (DST-I), dx = Lx/(nx+1).
    # ------------------------------------------------------------------
    nx: int = 127
    ny: int = 127
    dt: float = 1e-5
    n_frames: int = 21

    # ------------------------------------------------------------------
    # PINN architecture and training
    # ------------------------------------------------------------------
    hidden_layers: int = 4
    hidden_units: int = 64
    use_fourier: bool = True
    fourier_num_freq: int = 6
    fourier_max_freq: float = 6.0

    n_interior: int = 4000
    n_norm_side: int = 48       # quadrature grid side for the norm loss
    n_norm_times: int = 4
    sampling: str = "uniform"   # "uniform" or "rad" (residual-adaptive)
    rad_pool_factor: int = 8    # candidate pool size for "rad"

    epochs_adam: int = 2000
    lr_adam: float = 1e-3
    lr_decay_every: int = 1000
    lr_decay_gamma: float = 0.5
    use_lbfgs: bool = True
    epochs_lbfgs: int = 200
    loss_weights: dict = field(default_factory=lambda: dict(pde=1.0, norm=1.0))
    pde_warmup_epochs: int = 300
    grad_clip_norm: float = 5.0
    resample_every: int = 1         # Adam: fresh collocation batch every step
    rad_resample_every: int = 100   # "rad" pool re-scoring interval
    lbfgs_points_factor: int = 4    # L-BFGS uses a fixed set of factor * n_interior
    causal_eps: float = 0.0         # >0 enables causal weighting
    causal_bins: int = 16
    log_every: int = 250

    seed: int = 42
    device: str = "cpu"

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------
    def replace(self, **changes) -> "Config":
        return dataclasses.replace(self, **changes)

    def to_dict(self) -> dict:
        d = dataclasses.asdict(self)
        if d["disorder_region"] is not None:
            d["disorder_region"] = list(d["disorder_region"])
        return d

    @classmethod
    def from_dict(cls, d: dict) -> "Config":
        d = dict(d)
        if d.get("disorder_region") is not None:
            d["disorder_region"] = tuple(d["disorder_region"])
        names = {f.name for f in dataclasses.fields(cls)}
        return cls(**{k: v for k, v in d.items() if k in names})

    def save_json(self, path: str) -> None:
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, indent=2)

    @classmethod
    def load_json(cls, path: str) -> "Config":
        with open(path) as f:
            return cls.from_dict(json.load(f))

    # Physical derived quantities -----------------------------------------
    @property
    def k_transverse(self) -> float:
        return self.mode_n * math.pi / self.Ly if self.packet == "waveguide_mode" else self.ky0

    @property
    def mean_energy(self) -> float:
        """Kinetic energy of the packet's central wave vector (not <H>)."""
        return 0.5 * (self.kx0**2 + self.k_transverse**2)

    @property
    def dx(self) -> float:
        return self.Lx / (self.nx + 1)

    @property
    def dy(self) -> float:
        return self.Ly / (self.ny + 1)


# ----------------------------------------------------------------------
# Presets
# ----------------------------------------------------------------------
def legacy_config() -> Config:
    """Original project setup (commit 42cea06). Documented as NOT converged
    and physically ill-posed for transmission (see docs/AUDIT.md)."""
    return Config(
        Lx=1.0, Ly=1.0, T=1.0, packet="gaussian", x0=0.3, y0=0.5,
        sigma_x=0.06, sigma_y=0.06, kx0=40.0, ky0=0.0,
        barrier_height=600.0, barrier_center=0.55, barrier_width=0.06,
        barrier_edge=0.0, nx=128, ny=128, dt=1.0 / 4000, n_frames=60,
        fourier_max_freq=10.0, epochs_adam=4000, epochs_lbfgs=400,
    )


def pinn_config(**overrides) -> Config:
    """PINN-vs-reference problem: unit box, Gaussian packet (sigma = 0.07,
    5 sigma from the nearest wall, so |psi0| < 4e-5 on the walls) with
    kx0 = 15 (E0 = 112.5) and T = 0.02. The packet centre moves 0.3 and the
    packet width grows from 0.07 to ~0.29 (free dispersion), so it
    starts to touch the right wall only at the very end. Disorder (if
    enabled) covers the whole box. Reference grid 127^2, dt = 1e-5
    (convergence shown in results/validation)."""
    cfg = Config()
    return cfg.replace(**overrides)


def waveguide_config(**overrides) -> Config:
    """Open quasi-1D waveguide for transport experiments.

    Width Ly = 1 with hard walls, incident transverse mode n = 1 and
    kx0 = 20 (E0 ~ 205, wavelength ~0.31, 6 open transverse channels).
    Layout along x (lengths in units of Ly):
        [0, 0.6]                left CAP (absorbs reflected waves)
        packet at x0 = 1.2, sigma_x = 0.2
        [1.9, 1.9 + Ld]         disordered slab (default Ld = 1)
        gap 0.4, then 0.6 right CAP (absorbs transmitted waves)
    Grid: dx = dy = 1/64. dt = 2e-4. Convergence of all of these is
    checked in results/validation.
    """
    Ld = overrides.pop("slab_length", 1.0)
    x_s = 1.9
    Lx = x_s + Ld + 0.4 + 0.6
    dx_target = 1.0 / 64
    nx = int(round(Lx / dx_target)) - 1
    cfg = Config(
        Lx=Lx, Ly=1.0, T=1.0, packet="waveguide_mode", x0=1.2, y0=0.5,
        sigma_x=0.2, sigma_y=0.2, kx0=20.0, ky0=0.0, mode_n=1,
        barrier_height=0.0, barrier_center=x_s + Ld / 2, barrier_width=0.1,
        barrier_edge=0.01,
        disorder_strength=0.0, disorder_corr_length=0.05,
        disorder_region=(x_s, x_s + Ld), disorder_edge=0.02,
        cap_width=0.6, cap_strength=300.0, cap_power=2,
        nx=nx, ny=63, dt=2e-4, n_frames=51,
    )
    return cfg.replace(**overrides)


def closed_box_config(L: float = 1.0, **overrides) -> Config:
    """Closed L x L box, disorder everywhere, packet launched from the
    centre with kx0 = 20. Used to measure spreading (MSD) and the
    participation ratio, and their dependence on the box size L."""
    cfg = Config(
        Lx=L, Ly=L, T=0.5, packet="gaussian", x0=L / 2, y0=L / 2,
        sigma_x=0.08, sigma_y=0.08, kx0=20.0, ky0=0.0,
        disorder_strength=0.0, disorder_corr_length=0.05,
        disorder_region=None, nx=int(round(64 * L)) - 1, ny=int(round(64 * L)) - 1,
        dt=2e-4, n_frames=51,
    )
    return cfg.replace(**overrides)


# Backwards-compatible default instance (PINN problem).
cfg = pinn_config()
