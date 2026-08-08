"""
Classical (non-neural) reference solver for the 2D TDSE, used as ground
truth to check the PINN's accuracy.

Method: Strang-split split-step propagation, but using a 2D Discrete Sine
Transform (DST-I) instead of the usual FFT.

Why DST instead of the textbook FFT split-step method
-------------------------------------------------------
The standard split-step Fourier method assumes *periodic* boundary
conditions, because the FFT diagonalizes the Laplacian on a ring. Our PINN
is trained with *Dirichlet* boundaries (Psi = 0 at the walls of an infinite
well) -- a periodic reference would therefore be solving a different
physical problem and any comparison would be unfair.

The DST-I basis functions sin(n*pi*x/L) are exactly the eigenfunctions of
the Dirichlet Laplacian on [0, L], with eigenvalues (n*pi/L)^2. Using DST-I
in place of FFT gives a spectral method that enforces Psi = 0 at both walls
exactly, on both spatial axes, matching the PINN's boundary condition
by construction rather than by penalty.

Grid convention: DST-I operates on N interior points, with the (implicit,
not stored) boundary points at index -1 and N fixed at zero. So
x_i = (i + 1) * dx, dx = L / (N + 1), i = 0 .. N-1.
"""

import numpy as np
from scipy.fft import dstn, idstn

from config import cfg
from potential import potential_numpy


def _dst2(real_field: np.ndarray) -> np.ndarray:
    return dstn(real_field, type=1, norm="ortho")


def _idst2(coeffs: np.ndarray) -> np.ndarray:
    return idstn(coeffs, type=1, norm="ortho")


def _complex_dst2(psi: np.ndarray) -> np.ndarray:
    return _dst2(psi.real) + 1j * _dst2(psi.imag)


def _complex_idst2(coeffs: np.ndarray) -> np.ndarray:
    return _idst2(coeffs.real) + 1j * _idst2(coeffs.imag)


def build_grid(n_grid: int = None):
    n_grid = n_grid or cfg.eval_n_grid
    dx = cfg.L / (n_grid + 1)
    xs = (np.arange(1, n_grid + 1)) * dx
    ys = xs.copy()
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return X, Y, dx


def _kinetic_phase(n_grid: int, dt: float) -> np.ndarray:
    """Diagonal kinetic propagator exp(-i * 0.5 * k^2 * dt) in DST space."""
    n = np.arange(1, n_grid + 1)
    k = n * np.pi / cfg.L
    KX, KY = np.meshgrid(k, k, indexing="ij")
    energy = 0.5 * (KX**2 + KY**2)
    return np.exp(-1j * energy * dt)


def initial_wavepacket(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    envelope = np.exp(-((X - cfg.x0) ** 2 + (Y - cfg.y0) ** 2) / (2 * cfg.sigma**2))
    phase = cfg.kx0 * X + cfg.ky0 * Y
    psi = envelope * (np.cos(phase) + 1j * np.sin(phase))
    # Normalize on the grid so integral(|psi|^2) dx dy = 1 exactly, matching
    # the (approximately normalized) analytic Gaussian used as the PINN's IC.
    dx = X[1, 0] - X[0, 0]
    norm = np.sqrt(np.sum(np.abs(psi) ** 2) * dx * dx)
    return psi / norm


def propagate(n_grid: int = None, n_steps: int = None, n_frames: int = None):
    """
    Propagate the initial wave packet from t=0 to t=T with Strang splitting.

    Returns
    -------
    X, Y : spatial grids
    frame_times : array of times at which frames were recorded
    frames : list of complex psi arrays (n_grid x n_grid) at each frame time
    """
    n_grid = n_grid or cfg.eval_n_grid
    n_steps = n_steps or cfg.eval_n_steps
    n_frames = n_frames or cfg.eval_n_frames

    X, Y, dx = build_grid(n_grid)
    V = potential_numpy(X, Y)
    dt = cfg.T / n_steps

    half_potential_phase = np.exp(-1j * V * dt / 2.0)
    kinetic_phase = _kinetic_phase(n_grid, dt)

    psi = initial_wavepacket(X, Y)

    frame_stride = max(1, n_steps // n_frames)
    frame_times = []
    frames = []

    for step in range(n_steps + 1):
        t = step * dt
        if step % frame_stride == 0 or step == n_steps:
            frame_times.append(t)
            frames.append(psi.copy())
        if step == n_steps:
            break
        psi = half_potential_phase * psi
        coeffs = _complex_dst2(psi)
        coeffs = kinetic_phase * coeffs
        psi = _complex_idst2(coeffs)
        psi = half_potential_phase * psi

    return X, Y, np.array(frame_times), frames


def total_probability(psi: np.ndarray, dx: float) -> float:
    return float(np.sum(np.abs(psi) ** 2) * dx * dx)


if __name__ == "__main__":
    print(f"Running reference split-step (DST) solver | scenario={cfg.SCENARIO}")
    X, Y, frame_times, frames = propagate()
    dx = X[1, 0] - X[0, 0]
    probs = [total_probability(p, dx) for p in frames]
    print(f"Total probability: start={probs[0]:.6f}  end={probs[-1]:.6f}  "
          f"drift={abs(probs[-1] - probs[0]):.2e}")
    np.savez(
        "reference_solution.npz",
        X=X, Y=Y, frame_times=frame_times,
        psi_real=np.array([p.real for p in frames]),
        psi_imag=np.array([p.imag for p in frames]),
    )
    print("Saved -> reference_solution.npz")
