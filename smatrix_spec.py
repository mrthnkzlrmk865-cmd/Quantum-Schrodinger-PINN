r"""
Reference implementation of the scattering-matrix CONVENTIONS for the
planned stationary (RGF) solver. See docs/PHASE2_METHOD_SELECTION.md,
section 2.3 "Mathematical specification".

This module is the executable form of the specification. It contains no
Hamiltonian assembly and no Green-function recursion: it implements only
(i) the discrete lead solutions, (ii) the lead self-energies and
broadenings, (iii) the Fisher-Lee / Caroli formulas that turn the endpoint
Green-function blocks into the flux-normalized scattering matrix, and
(iv) the unitarity and reciprocity checks. The future solver must call
these functions, so its conventions are fixed and tested independently.

Discretization (longitudinal finite differences, transverse sine modes):
  slices j along x with spacing dx; at each slice a vector c_j in C^{N_y}
  of transverse-mode amplitudes. With t_x = 1/(2 dx^2) (hbar = m = 1):
    (H c)_j = -t_x c_{j-1} + (2 t_x I + E_perp + V_j) c_j - t_x c_{j+1},
  where E_perp = diag((n pi / L_y)^2 / 2) and V_j is the real symmetric
  matrix V_nm(x_j). Inner product: sum_j c_j^dagger c_j (no dx factor).
  Leads: V_j = 0 for j <= 0 (left) and j >= N+1 (right); device j = 1..N.

Phase convention "P" (reference planes at the first and last device slices):
  left lead (j <= 1):  c_j = a_n lambda_n^{j-1} e_n + sum_m b_m lambda_m^{-(j-1)} e_m
  right lead (j >= N): c_j = sum_m d_m lambda_m^{j-N} e_m
  (the formal continuation of each lead solution to the adjacent device
  slice defines the reference amplitude there).

Index conventions: matrices are indexed [outgoing, incoming]. Channel
order inside a lead: ascending transverse index n, open channels only.
Full S-matrix (incoming order [L open..., R open...], outgoing order
[L open..., R open...]):

    S = [[ r_L , t_LR ],
         [ t_RL, r_R  ]]

  r_L  : L -> L (reflection, incidence from the left)
  t_RL : L -> R (transmission, incidence from the left), "to R from L"
  t_LR : R -> L (transmission, incidence from the right), "to L from R"
  r_R  : R -> R
"""

from __future__ import annotations

import numpy as np

THRESHOLD_TOL = 1e-12


def hopping(dx: float) -> float:
    return 1.0 / (2.0 * dx * dx)


def lead_modes(E: float, E_perp, dx: float):
    """Retarded lead solutions for each transverse channel.

    Channel n obeys  -t c_{j-1} + (2t + E_perp_n - E) c_j - t c_{j+1} = 0.
    With c_j = lambda^j:  lambda + 1/lambda = 2 c_n,
        c_n = 1 - (E - E_perp_n)/(2t) = 1 - (E - E_perp_n) dx^2.
    Discrete dispersion (open channel, lambda = e^{i k dx}):
        E = E_perp_n + (1 - cos(k_n dx)) / dx^2.
    Open:   |c_n| < 1, k_n dx = arccos(c_n) in (0, pi), lambda = e^{+i k_n dx}
            (the retarded root, |lambda| -> 1^- for E + i0),
            discrete group velocity v_n = dE/dk = sin(k_n dx)/dx.
    Closed: |c_n| > 1, lambda real with |lambda| < 1 (decays away from the
            device); v_n = 0.
    Returns dict(lam, open, k, v, sigma, gamma) with
        sigma_n = -t lambda_n     (retarded self-energy, mode-diagonal)
        gamma_n = i(sigma_n - conj(sigma_n)) = 2 t Im(lambda_n) = v_n / dx.
    """
    E_perp = np.asarray(E_perp, dtype=float)
    t = hopping(dx)
    c = 1.0 - (E - E_perp) / (2.0 * t)
    if np.any(np.abs(np.abs(c) - 1.0) < THRESHOLD_TOL):
        raise ValueError("energy at a channel threshold or band edge (|c_n| = 1): undefined")
    is_open = np.abs(c) < 1.0
    lam = np.empty(c.shape, dtype=complex)
    kdx = np.zeros(c.shape)
    kdx[is_open] = np.arccos(c[is_open])
    lam[is_open] = np.exp(1j * kdx[is_open])
    below = (~is_open) & (c > 1.0)
    above = (~is_open) & (c < -1.0)
    lam[below] = c[below] - np.sqrt(c[below] ** 2 - 1.0)      # in (0, 1)
    lam[above] = c[above] + np.sqrt(c[above] ** 2 - 1.0)      # in (-1, 0); lattice artifact
    v = np.where(is_open, np.sin(kdx) / dx, 0.0)
    sigma = -t * lam
    gamma = np.where(is_open, 2.0 * t * np.sin(kdx), 0.0)     # = v / dx, exactly 0 for closed
    return dict(lam=lam, open=is_open, k=kdx / dx, v=v, sigma=sigma, gamma=gamma)


def fisher_lee(G11, GN1, G1N, GNN, gamma_L, gamma_R, open_L, open_R):
    """Flux-normalized scattering amplitudes (convention P), open channels only.

        t_RL = i Gamma_R^{1/2} G_{N1} Gamma_L^{1/2}
        t_LR = i Gamma_L^{1/2} G_{1N} Gamma_R^{1/2}
        r_L  = i Gamma_L^{1/2} G_{11} Gamma_L^{1/2} - I
        r_R  = i Gamma_R^{1/2} G_{NN} Gamma_R^{1/2} - I

    G_ab are N_y x N_y blocks of the retarded device Green function
    G^r = [E - H_D - Sigma_L - Sigma_R]^{-1} (blocks of slices a, b);
    gamma_* are the diagonal broadenings; only open rows and columns are kept.
    """
    sL = np.sqrt(gamma_L[open_L])
    sR = np.sqrt(gamma_R[open_R])
    oL, oR = np.where(open_L)[0], np.where(open_R)[0]
    t_RL = 1j * sR[:, None] * GN1[np.ix_(oR, oL)] * sL[None, :]
    t_LR = 1j * sL[:, None] * G1N[np.ix_(oL, oR)] * sR[None, :]
    r_L = 1j * sL[:, None] * G11[np.ix_(oL, oL)] * sL[None, :] - np.eye(len(oL))
    r_R = 1j * sR[:, None] * GNN[np.ix_(oR, oR)] * sR[None, :] - np.eye(len(oR))
    return dict(r_L=r_L, t_RL=t_RL, t_LR=t_LR, r_R=r_R)


def assemble_S(r_L, t_LR, t_RL, r_R):
    """S = [[r_L, t_LR], [t_RL, r_R]]; rows outgoing [L; R], columns incoming [L; R]."""
    nL, nR = r_L.shape[0], r_R.shape[0]
    if t_LR.shape != (nL, nR) or t_RL.shape != (nR, nL):
        raise ValueError("block shapes inconsistent with [[r_L, t_LR], [t_RL, r_R]]")
    return np.block([[r_L, t_LR], [t_RL, r_R]])


def split_S(S, n_L):
    return dict(r_L=S[:n_L, :n_L], t_LR=S[:n_L, n_L:], t_RL=S[n_L:, :n_L], r_R=S[n_L:, n_L:])


def caroli_transmission(GN1, gamma_L, gamma_R):
    """T(E) = Tr[Gamma_R G^r_{N1} Gamma_L G^a_{1N}],  G^a_{1N} = (G^r_{N1})^dagger.
    Independent of the channel phase convention; equals Tr(t_RL^dagger t_RL)."""
    GL = np.diag(gamma_L)
    GR = np.diag(gamma_R)
    return float(np.real(np.trace(GR @ GN1 @ GL @ GN1.conj().T)))


def unitarity_residual(S):
    """max |S^dagger S - I| (flux conservation over all open channels, both leads)."""
    return float(np.max(np.abs(S.conj().T @ S - np.eye(S.shape[0]))))


def reciprocity_residual(S):
    """max |S - S^T| (time-reversal symmetry, real symmetric H, convention P)."""
    return float(np.max(np.abs(S - S.T)))


def mode_reciprocity_residual(t_RL, t_LR):
    """max_{m,n} |T^{L->R}_{m<-n} - T^{R->L}_{n<-m}|, with
    T^{L->R}_{m<-n} = |t_RL[m, n]|^2 (in: left channel n, out: right channel m) and
    T^{R->L}_{n<-m} = |t_LR[n, m]|^2 (in: right channel m, out: left channel n).
    This is the correct probability-level reciprocity relation. It compares
    OPPOSITE propagation directions. Same-direction symmetry
    |t_RL[m,n]|^2 = |t_RL[n,m]|^2 does NOT follow from time-reversal symmetry."""
    return float(np.max(np.abs(np.abs(t_RL) ** 2 - (np.abs(t_LR) ** 2).T)))
