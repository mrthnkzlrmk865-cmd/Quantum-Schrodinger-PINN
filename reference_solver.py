"""
Classical (non-neural) reference solver for the 2D TDSE

    i dpsi/dt = [-1/2 (d_xx + d_yy) + V(x,y) - i W_cap(x)] psi

on [0, Lx] x [0, Ly] with psi = 0 on the walls.

Method: Strang split-step propagation in a 2D Discrete Sine Transform
(DST-I) basis (Feit, Fleck & Steiger 1982 split-operator method with sine
instead of Fourier modes, Strang 1968 splitting):

    psi(t+dt) = e^{-(iV + W) dt/2} DST^-1[ e^{-i |k|^2 dt/2} DST[ e^{-(iV + W) dt/2} psi ] ]

DST-I basis functions sin(n pi x / Lx) are the Dirichlet-Laplacian
eigenfunctions, with eigenvalues (n pi / Lx)^2, so the kinetic step is
exact for the sine-series representation and psi = 0 on the walls by
construction. Grid: x_i = (i+1) dx, dx = Lx/(nx+1), i = 0..nx-1 (the wall
points are implicit). Same in y.

Error: the scheme is second order in dt (splitting error from the
commutator [T, V]) and spectrally accurate in space for smooth V and psi.
Without the CAP the scheme is exactly unitary (norm conserved to rounding),
which does NOT by itself prove accuracy; the accuracy is demonstrated
separately in validate_solver.py (exact eigenmode phase, analytic free
Gaussian packet, dt- and dx-refinement).

Complex absorbing potential (optional, cfg.cap_width > 0):
    W_cap(x) = cap_strength * s(x)^cap_power,  s = depth into the layer / width,
on both x-ends. It removes outgoing waves so that the box behaves like an
open waveguide. The probability removed in each half-step,
sum |psi|^2 (1 - e^{-W dt}) dx dy, is accumulated separately for the left
and right layers. The kinetic step is unitary, so this bookkeeping is exact:
    norm(t) + absorbed_left(t) + absorbed_right(t) = 1.
absorbed_right(t -> inf) is the transmission probability and
absorbed_left the reflection probability. Spurious CAP reflection is
measured in validate_solver.py.

Multiple disorder realizations can be propagated as a batch (leading axis),
which is faster than running them one by one.
"""

from __future__ import annotations

import numpy as np
from scipy.fft import dstn

from initial_condition import initial_condition_numpy
from potential import Potential


def build_grid(cfg, nx=None, ny=None):
    nx = nx or cfg.nx
    ny = ny or cfg.ny
    dx = cfg.Lx / (nx + 1)
    dy = cfg.Ly / (ny + 1)
    xs = np.arange(1, nx + 1) * dx
    ys = np.arange(1, ny + 1) * dy
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return X, Y, dx, dy


def cap_profile(cfg, xs):
    """W_cap(x) >= 0 on the 1D x-grid."""
    W = np.zeros_like(xs)
    if cfg.cap_width <= 0:
        return W
    w = cfg.cap_width
    sl = np.clip((w - xs) / w, 0, None)
    sr = np.clip((xs - (cfg.Lx - w)) / w, 0, None)
    W = cfg.cap_strength * (sl**cfg.cap_power + sr**cfg.cap_power)
    return W


class SplitStepDST:
    def __init__(self, cfg, potentials=None, nx=None, ny=None, workers=-1):
        """potentials: None (use Potential(cfg)), a Potential, a list of
        Potentials (batch), or an array (B, nx, ny) / (nx, ny) of values."""
        self.cfg = cfg
        self.X, self.Y, self.dx, self.dy = build_grid(cfg, nx, ny)
        self.nx, self.ny = self.X.shape
        self.workers = workers
        if potentials is None:
            potentials = [Potential(cfg)]
        if isinstance(potentials, Potential):
            potentials = [potentials]
        if isinstance(potentials, np.ndarray):
            V = potentials if potentials.ndim == 3 else potentials[None]
        else:
            V = np.stack([p.numpy(self.X, self.Y) for p in potentials])
        self.V = V.astype(np.float64)
        self.B = self.V.shape[0]

        xs = self.X[:, 0]
        self.W_cap = cap_profile(cfg, xs)[None, :, None]  # (1, nx, 1)
        # index ranges of the left / right absorbing layers (contiguous)
        idx = np.nonzero(cap_profile(cfg, xs) > 0)[0]
        self.cap_left = slice(0, int(idx[idx < self.nx // 2].max()) + 1) if np.any(idx < self.nx // 2) else slice(0, 0)
        self.cap_right = slice(int(idx[idx >= self.nx // 2].min()), self.nx) if np.any(idx >= self.nx // 2) else slice(0, 0)

        kx = np.arange(1, self.nx + 1) * np.pi / cfg.Lx
        ky = np.arange(1, self.ny + 1) * np.pi / cfg.Ly
        KX, KY = np.meshgrid(kx, ky, indexing="ij")
        self.kinetic_energy = 0.5 * (KX**2 + KY**2)

    # ------------------------------------------------------------------
    def dst(self, psi):
        return dstn(psi, type=1, norm="ortho", axes=(-2, -1), workers=self.workers)

    # DST-I with orthonormal scaling is its own inverse.
    idst = dst

    def norm(self, psi):
        return np.sum(np.abs(psi) ** 2, axis=(-2, -1)) * self.dx * self.dy

    def initial_state(self):
        psi0 = initial_condition_numpy(self.cfg, self.X, self.Y)
        return np.broadcast_to(psi0, (self.B,) + psi0.shape).copy()

    # ------------------------------------------------------------------
    def run(self, psi0=None, t_end=None, dt=None, n_frames=None,
            store_frames=True, observables=None, normalize_initial=False):
        """Propagate to t_end. Returns a dict with
            times (F,), norm (F, B), absorbed_left (F, B), absorbed_right (F, B),
            frames (F, B, nx, ny) complex  [if store_frames],
            and one (F, B) array per observable returned by the callback
            ``observables(psi, solver) -> dict``.
        dt is adjusted down slightly so that t_end / dt is an integer."""
        cfg = self.cfg
        t_end = cfg.T if t_end is None else t_end
        dt = cfg.dt if dt is None else dt
        n_frames = cfg.n_frames if n_frames is None else n_frames
        n_steps = max(1, int(np.ceil(t_end / dt - 1e-9)))
        dt = t_end / n_steps
        record = set(np.unique(np.round(np.linspace(0, n_steps, n_frames)).astype(int)).tolist())

        psi = self.initial_state() if psi0 is None else np.array(psi0, dtype=complex, copy=True)
        if psi.ndim == 2:
            psi = np.broadcast_to(psi, (self.B,) + psi.shape).copy()
        if normalize_initial:
            psi /= np.sqrt(self.norm(psi))[:, None, None]

        half_V = np.exp(-1j * self.V * dt / 2.0)
        has_cap = np.any(self.W_cap > 0)
        if has_cap:
            damp = np.exp(-self.W_cap * dt / 2.0)
            loss_frac = 1.0 - damp**2
            half = half_V * damp
        else:
            half = half_V
        kin = np.exp(-1j * self.kinetic_energy * dt)

        abs_l = np.zeros(self.B)
        abs_r = np.zeros(self.B)
        out = dict(times=[], norm=[], absorbed_left=[], absorbed_right=[])
        frames = []
        extra = {}

        dA = self.dx * self.dy
        sl, sr = self.cap_left, self.cap_right

        def absorb(p):
            # exact probability removed by the next half-step damping,
            # evaluated only on the absorbing layers
            if has_cap:
                pl, pr = p[:, sl, :], p[:, sr, :]
                abs_l[:] += np.sum((pl.real**2 + pl.imag**2) * loss_frac[:, sl, :], axis=(-2, -1)) * dA
                abs_r[:] += np.sum((pr.real**2 + pr.imag**2) * loss_frac[:, sr, :], axis=(-2, -1)) * dA

        for step in range(n_steps + 1):
            if step in record:
                out["times"].append(step * dt)
                out["norm"].append(self.norm(psi))
                out["absorbed_left"].append(abs_l.copy())
                out["absorbed_right"].append(abs_r.copy())
                if store_frames:
                    frames.append(psi.copy())
                if observables is not None:
                    for k, v in observables(psi, self).items():
                        extra.setdefault(k, []).append(np.asarray(v))
            if step == n_steps:
                break
            absorb(psi)
            psi = half * psi
            psi = self.idst(kin * self.dst(psi))
            absorb(psi)
            psi = half * psi

        res = {k: np.array(v) for k, v in out.items()}
        res.update({k: np.array(v) for k, v in extra.items()})
        if store_frames:
            res["frames"] = np.array(frames)
        res["dt"] = dt
        res["n_steps"] = n_steps
        return res


# ----------------------------------------------------------------------
# Backwards-compatible helpers
# ----------------------------------------------------------------------
def propagate(cfg, nx=None, ny=None, dt=None, n_frames=None, t_end=None, potential=None):
    """Single-realization propagation. Returns X, Y, frame_times, frames."""
    solver = SplitStepDST(cfg, potential, nx=nx, ny=ny)
    res = solver.run(dt=dt, n_frames=n_frames, t_end=t_end)
    return solver.X, solver.Y, res["times"], res["frames"][:, 0]


def total_probability(psi, dx, dy=None):
    dy = dx if dy is None else dy
    return float(np.sum(np.abs(psi) ** 2) * dx * dy)


if __name__ == "__main__":
    import argparse
    import os

    from config import pinn_config

    ap = argparse.ArgumentParser(description="Reference solution for the PINN problem")
    ap.add_argument("--strength", type=float, default=0.0, help="disorder RMS W")
    ap.add_argument("--zeta", type=float, default=0.08)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="results/reference_solution.npz")
    a = ap.parse_args()
    cfg = pinn_config(disorder_strength=a.strength, disorder_corr_length=a.zeta, disorder_seed=a.seed)
    X, Y, t, frames = propagate(cfg)
    dx, dy = cfg.dx, cfg.dy
    probs = [total_probability(p, dx, dy) for p in frames]
    print(f"norm: start={probs[0]:.10f} end={probs[-1]:.10f}")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    np.savez_compressed(a.out, X=X, Y=Y, frame_times=t, psi_real=frames.real,
                        psi_imag=frames.imag, config=str(cfg.to_dict()))
    print(f"Saved -> {a.out}")
