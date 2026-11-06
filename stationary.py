r"""
Stationary scattering solver: mode-space recursive Green's function (RGF).

Implements the contract of docs/PHASE2_METHOD_SELECTION.md §2.7. All lead
quantities (lead roots, self-energies, broadenings), the Fisher-Lee
amplitudes, the S-matrix assembly and the unitarity/reciprocity checks
come from smatrix_spec.py, which is the tested specification. This module
only builds the device Hamiltonian and computes the endpoint
Green-function blocks.

Model (hbar = m = 1):  E psi = [-1/2 (d_xx + d_yy) + V(x, y)] psi on the strip
0 < y < L_y with psi = 0 at y = 0, L_y; clean leads (V = 0) for x outside
the device.

Discretization (spec §2.7.1):
  * slices x_j = x_start + (j - 1) dx, j = 1..N (device; it includes buffers);
    left lead j <= 0, right lead j >= N + 1, both clean;
  * transverse basis chi_n(y) = sqrt(2/L_y) sin(n pi y / L_y), n = 1..N_y,
    exact Dirichlet eigenfunctions, E_perp_n = (n pi / L_y)^2 / 2;
  * H_jj = 2 t_x I + diag(E_perp) + V_j,  H_{j,j+1} = H_{j+1,j} = -t_x I,
    t_x = 1/(2 dx^2);
  * (V_j)_nm = int chi_n V(x_j, y) chi_m dy by the DST-I quadrature
    y_l = l L_y/(N_q + 1), l = 1..N_q (N_q >= 2 N_y):
    V_j = X^T diag(V(x_j, y_l)) X L_y/(N_q + 1), X_ln = chi_n(y_l).

Retarded device Green function G = M^{-1}, M = E I - H_D - Sigma_L (+) Sigma_R,
with Sigma_L = Sigma_R = diag(-t_x lambda_n) from smatrix_spec.lead_modes
for ALL N_y channels: closed channels keep their real self-energy (they
shift the energy of the evanescent modes) although their broadening is 0.

Recursion (block Gaussian elimination from the left; M_{j,j+-1} = +t_x I):
    A_j  = E I - H_jj - Sigma_L delta_{j1} - Sigma_R delta_{jN}
    D_1  = A_1,  D_j = A_j - t_x^2 g_{j-1}           (Schur complement)
    g_j  = D_j^{-1}                                   (left-connected G of slices 1..j)
    G_{j1} = -t_x g_j G_{j-1,1}                       (G_11 = g_1 at j = 1)
    G_{1j} = -t_x G_{1,j-1} g_j
    G_11  <- G_11 + t_x^2 G_{1,j-1} g_j G_{j-1,1}
  at j = N: G_NN = g_N and G_N1, G_1N, G_11 of the full device.
Each D_j is LU-factorized once; g_j and g_j G_{j-1,1} come from one
lu_solve with stacked right-hand sides (no explicit inverse of the full
device matrix, and no explicit inverse of D_j other than the solve with I
that defines g_j). G_1N is computed by its own recursion, NOT as G_N1^T,
so that the reciprocity check S = S^T is a genuine numerical test.

The reciprocal 1-norm condition number of every D_j is estimated
(LAPACK gecon) and reported. A singular D_j or a non-finite result raises
an error. S and the Green-function outputs are never clipped or
renormalized.

Performance note (measured). With a multithreaded OpenBLAS (the default,
4 threads on the 4-core validation machine), each 63x63 complex BLAS call
in this loop suffered quantized ~4-16 ms thread wake-up latency; one
718-slice energy took 16-22 s. With OPENBLAS_NUM_THREADS=1 (or 2) it took
0.4-0.6 s. RGF scripts therefore set OPENBLAS_NUM_THREADS=1 before NumPy
is imported (validate_rgf_b1.py, tests/conftest.py). The effect of the
thread count on the results is measured in docs/RGF_B1_VALIDATION.md.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.linalg as sla

import smatrix_spec as sm


# ----------------------------------------------------------------------
# Transverse basis and potential blocks
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class TransverseBasis:
    Ly: float
    n_modes: int          # N_y
    n_quad: int           # N_q
    y: np.ndarray         # (N_q,) quadrature points
    X: np.ndarray         # (N_q, N_y) chi_n(y_l)
    weight: float         # L_y / (N_q + 1)
    E_perp: np.ndarray    # (N_y,)


def transverse_basis(Ly: float, n_modes: int, n_quad: int | None = None) -> TransverseBasis:
    if not (Ly > 0 and np.isfinite(Ly)):
        raise ValueError(f"Ly must be finite and > 0, got {Ly}")
    if int(n_modes) != n_modes or n_modes < 1:
        raise ValueError(f"n_modes must be a positive integer, got {n_modes}")
    n_modes = int(n_modes)
    n_quad = 2 * n_modes + 1 if n_quad is None else n_quad
    if int(n_quad) != n_quad or n_quad < 2 * n_modes:
        raise ValueError(f"n_quad must be an integer >= 2 n_modes = {2 * n_modes} (spec §2.7.1), got {n_quad}")
    n_quad = int(n_quad)
    l = np.arange(1, n_quad + 1)
    n = np.arange(1, n_modes + 1)
    y = l * Ly / (n_quad + 1)
    X = np.sqrt(2.0 / Ly) * np.sin(np.outer(l, n) * np.pi / (n_quad + 1))
    E_perp = 0.5 * (n * np.pi / Ly) ** 2
    return TransverseBasis(Ly=float(Ly), n_modes=n_modes, n_quad=n_quad, y=y, X=X,
                           weight=Ly / (n_quad + 1), E_perp=E_perp)


def potential_blocks(V_values: np.ndarray, basis: TransverseBasis) -> np.ndarray:
    """V_j = X^T diag(V(x_j, y_l)) X * weight for every slice.
    V_values: (N, N_q) real samples. Returns (N, N_y, N_y), exactly
    symmetric (the product is symmetric up to rounding; it is symmetrized,
    which changes it by at most the rounding of the product)."""
    V_values = np.asarray(V_values, dtype=float)
    if V_values.ndim != 2 or V_values.shape[1] != basis.n_quad:
        raise ValueError(f"V_values must have shape (N, {basis.n_quad}), got {V_values.shape}")
    if not np.all(np.isfinite(V_values)):
        raise ValueError("non-finite potential values")
    X = basis.X
    B = np.einsum("ln,jl,lm->jnm", X, V_values, X, optimize=True) * basis.weight
    return 0.5 * (B + np.transpose(B, (0, 2, 1)))


# ----------------------------------------------------------------------
# Device
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Device:
    xs: np.ndarray          # (N,) slice positions
    dx: float
    basis: TransverseBasis
    V: np.ndarray           # (N, N_y, N_y) potential blocks
    V_lead_max: float       # max |V| sampled on the first lead slice on each side

    @property
    def n_slices(self) -> int:
        return self.xs.size

    @property
    def t_x(self) -> float:
        return sm.hopping(self.dx)


def build_device(potential_numpy, x_start: float, x_end: float, dx: float, Ly: float, n_modes: int,
                 n_quad: int | None = None, lead_tol: float = 1e-6) -> Device:
    """Device from a potential callable V(X, Y) (vectorized, e.g.
    potential.Potential(cfg).numpy) on x_start..x_end with spacing dx
    (x_end - x_start must be an integer multiple of dx).
    The leads are clean by assumption: |V| on the slices x_start - dx and
    x_end + dx (first lead slices) must be <= lead_tol, else ValueError."""
    if not (dx > 0 and np.isfinite(dx)):
        raise ValueError(f"dx must be finite and > 0, got {dx}")
    span = (x_end - x_start) / dx
    n_int = int(round(span))
    if n_int < 1 or abs(span - n_int) > 1e-9 * max(1.0, span):
        raise ValueError(f"(x_end - x_start)/dx = {span} must be a positive integer")
    basis = transverse_basis(Ly, n_modes, n_quad)
    xs = x_start + dx * np.arange(n_int + 1)
    Xg, Yg = np.meshgrid(xs, basis.y, indexing="ij")
    Vv = np.asarray(potential_numpy(Xg, Yg), dtype=float)
    lead_x = np.array([x_start - dx, x_end + dx])
    Xl, Yl = np.meshgrid(lead_x, basis.y, indexing="ij")
    V_lead_max = float(np.max(np.abs(np.asarray(potential_numpy(Xl, Yl), dtype=float))))
    if V_lead_max > lead_tol:
        raise ValueError(f"potential on the first lead slices is {V_lead_max:g} > lead_tol = {lead_tol:g}; "
                         f"the clean-lead assumption is violated (enlarge the device)")
    return Device(xs=xs, dx=float(dx), basis=basis, V=potential_blocks(Vv, basis), V_lead_max=V_lead_max)


def device_hamiltonian_blocks(dev: Device):
    """Diagonal blocks H_jj (N, N_y, N_y) and the scalar hopping t_x
    (off-diagonal blocks are -t_x I)."""
    t = dev.t_x
    ny = dev.basis.n_modes
    diag = 2.0 * t * np.eye(ny) + np.diag(dev.basis.E_perp)
    return diag[None, :, :] + dev.V, t


# ----------------------------------------------------------------------
# RGF
# ----------------------------------------------------------------------
def _rcond(lu_piv, anorm):
    lu, piv = lu_piv
    gecon, = sla.get_lapack_funcs(("gecon",), (lu,))
    rc, info = gecon(lu, anorm, norm="1")
    if info != 0:
        raise np.linalg.LinAlgError(f"gecon failed (info={info})")
    return float(rc)


def rgf_endpoint_blocks(E: float, dev: Device, sigma: np.ndarray):
    """Endpoint blocks of G = [E - H_D - Sigma_L (+) Sigma_R]^{-1} by the
    left-to-right recursion. sigma: (N_y,) mode-diagonal lead self-energy
    (identical leads, all channels). Returns dict(G11, GN1, G1N, GNN,
    min_rcond, argmin_rcond)."""
    H, t = device_hamiltonian_blocks(dev)
    N, ny = dev.n_slices, dev.basis.n_modes
    sigma = np.asarray(sigma)
    if sigma.shape != (ny,):
        raise ValueError(f"sigma must have shape ({ny},), got {sigma.shape}")
    eye = np.eye(ny, dtype=complex)
    Sig = np.diag(sigma.astype(complex))
    min_rc, arg_rc = np.inf, -1
    g = G_j1 = G_1j = G11 = None
    for j in range(N):
        D = E * eye - H[j]
        if j == 0:
            D = D - Sig
        if j == N - 1:
            D = D - Sig
        if j > 0:
            D = D - (t * t) * g
        try:
            lu_piv = sla.lu_factor(D, check_finite=True)
        except (ValueError, np.linalg.LinAlgError) as exc:
            raise np.linalg.LinAlgError(f"RGF: LU of slice {j + 1} failed at E = {E}: {exc}") from exc
        rc = _rcond(lu_piv, np.linalg.norm(D, 1))
        if rc < min_rc:
            min_rc, arg_rc = rc, j + 1
        if rc == 0.0 or not np.isfinite(rc):
            raise np.linalg.LinAlgError(f"RGF: singular Schur complement at slice {j + 1}, E = {E}")
        if j == 0:
            g = sla.lu_solve(lu_piv, eye)
            G_j1, G_1j, G11 = g.copy(), g.copy(), g.copy()
            continue
        sol = sla.lu_solve(lu_piv, np.hstack([eye, G_j1]))        # [g_j, g_j G_{j-1,1}]
        g_new, gG = sol[:, :ny], sol[:, ny:]
        G11 = G11 + (t * t) * (G_1j @ gG)                          # + t^2 G_{1,j-1} g_j G_{j-1,1}
        G_1j = -t * (G_1j @ g_new)                                 # G_{1j} = -t G_{1,j-1} g_j
        G_j1 = -t * gG                                             # G_{j1} = -t g_j G_{j-1,1}
        g = g_new
    out = dict(G11=G11, GN1=G_j1, G1N=G_1j, GNN=g, min_rcond=min_rc, argmin_rcond=arg_rc)
    for k in ("G11", "GN1", "G1N", "GNN"):
        if not np.all(np.isfinite(out[k])):
            raise FloatingPointError(f"RGF: non-finite {k} at E = {E}")
    return out


def scattering_matrix(E: float, dev: Device, min_closed: int = 1):
    """Flux-normalized S = [[r_L, t_LR], [t_RL, r_R]] over the open channels
    (spec §2.7.5-2.7.7, via smatrix_spec). Requires at least `min_closed`
    closed channels in the basis (otherwise the transverse truncation
    would cut open channels). Returns dict with S, the four blocks,
    the lead modes, the Green blocks and diagnostics."""
    modes = sm.lead_modes(E, dev.basis.E_perp, dev.dx)
    n_open = int(np.sum(modes["open"]))
    if n_open == 0:
        raise ValueError(f"no open channel at E = {E}")
    if dev.basis.n_modes - n_open < min_closed:
        raise ValueError(f"only {dev.basis.n_modes - n_open} closed channels in the basis at E = {E} "
                         f"(need >= {min_closed}); increase n_modes")
    G = rgf_endpoint_blocks(E, dev, modes["sigma"])
    amp = sm.fisher_lee(G["G11"], G["GN1"], G["G1N"], G["GNN"], modes["gamma"], modes["gamma"],
                        modes["open"], modes["open"])
    S = sm.assemble_S(amp["r_L"], amp["t_LR"], amp["t_RL"], amp["r_R"])
    return dict(S=S, n_open=n_open, modes=modes, G=G, **amp,
                caroli=sm.caroli_transmission(G["GN1"], modes["gamma"], modes["gamma"]),
                min_rcond=G["min_rcond"], argmin_rcond=G["argmin_rcond"])
