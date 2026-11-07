r"""
B3X reference solver: S-matrix from the scattering WAVE FUNCTION by one sparse
linear solve (docs/RGF_B3X_VALIDATION.md, Part I.3).

Independent of the RGF and of the other validation reference:
  * does not import stationary.py, smatrix_spec.py or sparse_reference.py;
  * lead roots lambda, broadening Gamma: own closed form (the retarded choice
    is verified against the root of the E + i eta problem, see retarded_root_eta);
  * potential blocks V_j: cosine-sum form of the DST-I quadrature
        V_nm = (g_{|n-m|} - g_{n+m}) / (N_q + 1),  g_k = sum_l V(x_j, y_l) cos(pi k l / (N_q + 1)),
    instead of X^T diag(V) X;
  * linear algebra: SuperLU with the MMD_AT_PLUS_A column ordering (the other
    reference uses COLAMD) plus one step of iterative refinement;
  * S: from the scattering state. For incidence a lambda_n^{j-1} e_n in the left lead
    the elimination of the leads gives  M c = s  with s_1 = -t (1/lambda_n - lambda_n) a e_n,
    M = E - H_D + t Lambda on the two end blocks. Unit flux: a = Gamma_n^{-1/2}.
    r = Gamma^{1/2} (c_1 - a e_n), t = Gamma^{1/2} c_N (open channels), and the
    mirror problem for incidence from the right. No Green-function block and no
    Fisher-Lee formula is used.
Shared with the RGF (the problem definition, not an algorithm): the lattice
Hamiltonian (second-order differences, t = 1/(2 dx^2)), the Dirichlet sine basis,
the DST-I quadrature points, the potential callable, the channel ordering
(open channels ascending; S = [[r_L, t_LR], [t_RL, r_R]], rows outgoing).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

THRESHOLD_TOL = 1e-12


def hopping(dx: float) -> float:
    return 1.0 / (2.0 * dx * dx)


def lead_data(E: float, E_perp: np.ndarray, dx: float):
    """Retarded root lambda (|lambda| <= 1), open mask, Gamma = 2 t sin(k dx) (0 if closed).
    c_n = 1 - (E - E_perp_n) dx^2;  open: lambda = c + i sqrt(1 - c^2) = e^{+ i k dx};
    closed: the real root with |lambda| < 1."""
    E_perp = np.asarray(E_perp, float)
    t = hopping(dx)
    c = 1.0 - (E - E_perp) * dx * dx
    if np.any(np.abs(np.abs(c) - 1.0) < THRESHOLD_TOL):
        raise ValueError("energy at a channel threshold or band edge")
    is_open = np.abs(c) < 1.0
    lam = np.empty(c.shape, dtype=complex)
    s = np.sqrt(1.0 - c[is_open] ** 2)
    lam[is_open] = c[is_open] + 1j * s
    cc = c[~is_open]
    r = np.sqrt(cc * cc - 1.0)
    lam[~is_open] = np.where(cc > 0, cc - r, cc + r)
    gamma = np.zeros(c.shape)
    gamma[is_open] = 2.0 * t * s
    return dict(lam=lam, open=is_open, gamma=gamma, c=c)


def retarded_root_eta(E: float, E_perp: np.ndarray, dx: float, eta: float = 1e-9) -> np.ndarray:
    """Independent determination of the retarded root: the root with |lambda| < 1 of
    lambda^2 - 2 c lambda + 1 = 0 for the complex energy E + i eta."""
    c = 1.0 - (E + 1j * eta - np.asarray(E_perp, float)) * dx * dx
    out = np.empty(c.shape, dtype=complex)
    for i, ci in enumerate(c):
        roots = np.roots([1.0, -2.0 * ci, 1.0])
        out[i] = roots[np.argmin(np.abs(roots))]
    return out


def transverse_points(Ly: float, n_modes: int, n_quad: int):
    if n_quad < 2 * n_modes:
        raise ValueError("n_quad must be >= 2 n_modes")
    return np.arange(1, n_quad + 1) * Ly / (n_quad + 1), 0.5 * (np.arange(1, n_modes + 1) * np.pi / Ly) ** 2


def potential_blocks_cosine(V_samples: np.ndarray, n_modes: int, n_quad: int) -> np.ndarray:
    """V_j in the sine basis from cosine sums (identical quadrature, different algebra)."""
    th = np.pi / (n_quad + 1)
    l = np.arange(1, n_quad + 1)
    k = np.arange(0, 2 * n_modes + 1)
    g = V_samples @ np.cos(th * np.outer(k, l)).T                    # (N, 2 N_y + 1)
    n = np.arange(1, n_modes + 1)
    return (g[:, np.abs(n[:, None] - n[None, :])] - g[:, n[:, None] + n[None, :]]) / (n_quad + 1)


def assemble(E: float, V_blocks: np.ndarray, E_perp: np.ndarray, dx: float, lam: np.ndarray) -> sp.csc_matrix:
    """M = E - H_D + t Lambda on the first and last diagonal blocks (Sigma = -t Lambda); off-diagonal +t."""
    N, ny = V_blocks.shape[0], V_blocks.shape[1]
    t = hopping(dx)
    base = (E - 2.0 * t) * np.eye(ny) - np.diag(E_perp)
    blocks = [(base - V_blocks[j]).astype(complex) for j in range(N)]
    blocks[0] = blocks[0] + t * np.diag(lam)
    blocks[-1] = blocks[-1] + t * np.diag(lam)
    D = sp.block_diag([sp.csr_matrix(b) for b in blocks], format="csr")
    hop = sp.diags([np.ones(N - 1), np.ones(N - 1)], [-1, 1], shape=(N, N))
    return (D + t * sp.kron(hop, sp.identity(ny), format="csr")).tocsc()


def wavefunction_S_from_blocks(E: float, V_blocks: np.ndarray, E_perp: np.ndarray, dx: float,
                               _source_sign: float = 1.0):
    """S-matrix for given potential blocks (see the module docstring). `_source_sign` exists only so a
    test can show that a wrong source sign is caught by the unitarity check."""
    N, ny = V_blocks.shape[0], V_blocks.shape[1]
    t = hopping(dx)
    ld = lead_data(E, E_perp, dx)
    lam, op, gam = ld["lam"], ld["open"], ld["gamma"]
    o = np.where(op)[0]
    no = o.size
    if no == 0:
        raise ValueError("no open channel")
    M = assemble(E, V_blocks, E_perp, dx, lam)
    a = gam[o] ** -0.5                                               # unit incident flux
    B = np.zeros((N * ny, 2 * no), dtype=complex)
    src = _source_sign * (-t) * (1.0 / lam[o] - lam[o]) * a          # s = -t (1/lambda - lambda) a
    B[o, np.arange(no)] = src                                        # left incidence: slice 1
    B[(N - 1) * ny + o, no + np.arange(no)] = src                    # right incidence: slice N
    lu = spla.splu(M, permc_spec="MMD_AT_PLUS_A")
    X = lu.solve(B)
    res0 = B - M @ X
    X = X + lu.solve(res0)                                           # one step of iterative refinement
    res1 = B - M @ X
    scale = spla.norm(M, 1) * np.max(np.abs(X)) + np.max(np.abs(B))
    c1, cN = X[:ny, :], X[(N - 1) * ny:, :]
    sg = np.sqrt(gam[o])[:, None]
    r_L = sg * (c1[o][:, :no] - np.diag(a))
    t_RL = sg * cN[o][:, :no]
    t_LR = sg * c1[o][:, no:]
    r_R = sg * (cN[o][:, no:] - np.diag(a))
    S = np.block([[r_L, t_LR], [t_RL, r_R]])
    cl = np.where(~op)[0]
    closed_w = float(np.max(np.linalg.norm(c1[cl][:, :no], axis=0) / np.linalg.norm(c1[o][:, :no], axis=0))) if cl.size else 0.0
    return dict(S=S, r_L=r_L, t_RL=t_RL, t_LR=t_LR, r_R=r_R, n_open=int(no), lam=lam, open=op, gamma=gam,
                c1=c1, cN=cN, source=src, a=a,
                residual_before_refinement=float(np.max(np.abs(res0)) / scale),
                residual_after_refinement=float(np.max(np.abs(res1)) / scale),
                closed_to_open_weight_slice1=closed_w, lu_nnz=int(lu.L.nnz + lu.U.nnz), M_nnz=int(M.nnz))


def wavefunction_scattering_matrix(E: float, potential_numpy, x_start: float, dx: float, n_slices: int,
                                   Ly: float, n_modes: int, n_quad: int):
    """S-matrix of the device x_j = x_start + j dx (j = 0..n_slices-1) for a potential callable V(X, Y)."""
    y, E_perp = transverse_points(Ly, n_modes, n_quad)
    xs = x_start + dx * np.arange(n_slices)
    Xg, Yg = np.meshgrid(xs, y, indexing="ij")
    V = np.asarray(potential_numpy(Xg, Yg), dtype=float)
    Vb = potential_blocks_cosine(V, n_modes, n_quad)
    out = wavefunction_S_from_blocks(E, Vb, E_perp, dx)
    out.update(xs=xs, V_samples=V, V_blocks=Vb, E_perp=E_perp)
    return out
