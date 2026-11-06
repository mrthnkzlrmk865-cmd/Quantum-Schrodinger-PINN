r"""
Validation-only direct sparse Green-function reference for benchmark B3
(docs/PHASE2_METHOD_SELECTION.md §3, docs/RGF_B3_VALIDATION.md).

It computes the same quantity as stationary.py (the S-matrix of the
mode-space lattice problem of spec §2.7) by a different algorithm:

  * the complete device matrix  M = E I - H_D - Sigma_L (+) Sigma_R  of size
    (N N_y) x (N N_y) is assembled as ONE scipy.sparse matrix;
  * M is factorized once (SuperLU, scipy.sparse.linalg.splu, COLAMD column
    ordering) and the boundary Green-function blocks are obtained by direct
    triangular solves  M X = B  with the 2 N_y right-hand sides
    B = [e_slice1 (x) I, e_sliceN (x) I]:
        G_{a1} = X[block a, :N_y],   G_{aN} = X[block a, N_y:];
  * no recursion over slices, no Schur complements, no code from
    stationary.py (this module does not import it), no dense inverse.

Shared with the RGF by construction (the specified problem, spec §2.7):
  the Hamiltonian discretization (second-order differences, t = 1/(2 dx^2),
  Dirichlet sine modes, DST-I quadrature of the potential), the clean-lead
  roots, self-energy and broadening (smatrix_spec.lead_modes), the
  Fisher-Lee amplitudes and the S assembly (smatrix_spec.fisher_lee,
  assemble_S), and the channel ordering. The potential blocks are
  RE-DERIVED here from the potential callable (own coordinates, own
  transverse basis, own contraction), not taken from the RGF device.
So agreement validates the Green-function algebra and its numerical
implementation, not the physical model or the Fisher-Lee convention.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

import smatrix_spec as sm


def transverse_data(Ly: float, n_modes: int, n_quad: int):
    """Quadrature points y_l, sine matrix X_ln = chi_n(y_l), weight and E_perp
    (spec §2.7.1), computed here independently of stationary.transverse_basis."""
    if n_quad < 2 * n_modes:
        raise ValueError(f"n_quad = {n_quad} must be >= 2 n_modes = {2 * n_modes}")
    l = np.arange(1, n_quad + 1)
    n = np.arange(1, n_modes + 1)
    y = l * Ly / (n_quad + 1)
    X = np.sqrt(2.0 / Ly) * np.sin(np.pi * np.outer(y, n) / Ly)
    return y, X, Ly / (n_quad + 1), 0.5 * (n * np.pi / Ly) ** 2


def sample_potential(potential_numpy, x_start: float, dx: float, n_slices: int, y: np.ndarray):
    """V(x_j, y_l) on the slices x_j = x_start + j dx (j = 0..N-1) via the
    potential callable (vectorized V(X, Y))."""
    xs = x_start + dx * np.arange(n_slices)
    Xg, Yg = np.meshgrid(xs, y, indexing="ij")
    return xs, np.asarray(potential_numpy(Xg, Yg), dtype=float)


def potential_blocks(V_samples: np.ndarray, X: np.ndarray, weight: float) -> np.ndarray:
    """V_j = X^T diag(V(x_j, .)) X weight, one slice at a time (plain matmul)."""
    N = V_samples.shape[0]
    out = np.empty((N, X.shape[1], X.shape[1]))
    for j in range(N):
        out[j] = ((X.T * V_samples[j][None, :]) @ X) * weight
    return out


def assemble_M(E: float, V_blocks: np.ndarray, E_perp: np.ndarray, dx: float, sigma: np.ndarray):
    """Sparse M = E I - H_D - Sigma_L (+) Sigma_R (spec §2.7.4), complex CSC.
    H_jj = 2t I + diag(E_perp) + V_j, H_{j,j+-1} = -t I, so M_{j,j+-1} = +t I."""
    N, ny = V_blocks.shape[0], V_blocks.shape[1]
    t = sm.hopping(dx)
    eye = np.eye(ny)
    Sig = np.diag(np.asarray(sigma, dtype=complex))
    blocks = []
    for j in range(N):
        A = (E - 2.0 * t) * eye - np.diag(E_perp) - V_blocks[j]
        A = A.astype(complex)
        if j == 0:
            A = A - Sig
        if j == N - 1:
            A = A - Sig
        blocks.append(sp.csr_matrix(A))
    D = sp.block_diag(blocks, format="csr")
    hop = sp.diags([np.ones(N - 1), np.ones(N - 1)], [1, -1], shape=(N, N))
    return (D + t * sp.kron(hop, sp.identity(ny), format="csr")).tocsc()


def boundary_green_blocks(M: sp.csc_matrix, n_slices: int, n_modes: int):
    """G_{11}, G_{N1}, G_{1N}, G_{NN} by one sparse LU factorization and direct
    solves for the 2 N_y boundary right-hand sides. Returns the blocks and
    diagnostics (relative residual of M X = B, LU fill)."""
    N, ny = n_slices, n_modes
    lu = spla.splu(M)
    B = np.zeros((N * ny, 2 * ny), dtype=complex)
    B[:ny, :ny] = np.eye(ny)
    B[(N - 1) * ny:, ny:] = np.eye(ny)
    X = lu.solve(B)
    res = M @ X - B
    resid = float(np.max(np.abs(res)) / (spla.norm(M, 1) * np.max(np.abs(X)) + np.max(np.abs(B))))
    G = dict(G11=X[:ny, :ny], GN1=X[(N - 1) * ny:, :ny], G1N=X[:ny, ny:], GNN=X[(N - 1) * ny:, ny:])
    if not all(np.all(np.isfinite(g)) for g in G.values()):
        raise FloatingPointError("sparse reference: non-finite Green-function block")
    return G, dict(relative_residual=resid, lu_nnz=int(lu.L.nnz + lu.U.nnz), M_nnz=int(M.nnz))


def sparse_S_from_blocks(E: float, V_blocks: np.ndarray, E_perp: np.ndarray, dx: float):
    """Complete S for given potential blocks (used directly by the synthetic tests)."""
    N, ny = V_blocks.shape[0], V_blocks.shape[1]
    modes = sm.lead_modes(E, E_perp, dx)
    M = assemble_M(E, V_blocks, E_perp, dx, modes["sigma"])
    G, diag = boundary_green_blocks(M, N, ny)
    amp = sm.fisher_lee(G["G11"], G["GN1"], G["G1N"], G["GNN"], modes["gamma"], modes["gamma"],
                        modes["open"], modes["open"])
    S = sm.assemble_S(amp["r_L"], amp["t_LR"], amp["t_RL"], amp["r_R"])
    return dict(S=S, n_open=int(np.sum(modes["open"])), modes=modes, G=G, **amp, **diag)


def sparse_scattering_matrix(E: float, potential_numpy, x_start: float, dx: float, n_slices: int,
                             Ly: float, n_modes: int, n_quad: int):
    """S-matrix of the device x_j = x_start + j dx, j = 0..n_slices-1, for the
    potential callable, by the direct sparse method. Also returns the potential
    samples V(x_j, y_l) it evaluated (for the same-realization check)."""
    y, X, w, E_perp = transverse_data(Ly, n_modes, n_quad)
    xs, V = sample_potential(potential_numpy, x_start, dx, n_slices, y)
    Vb = potential_blocks(V, X, w)
    out = sparse_S_from_blocks(E, Vb, E_perp, dx)
    out.update(xs=xs, V_samples=V, V_blocks=Vb)
    return out
