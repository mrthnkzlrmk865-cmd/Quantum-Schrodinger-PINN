r"""
EXPERIMENTAL: W/Z endpoint-block recursion of wz_rgf.wz_endpoint_blocks with a pluggable linear-solve method and per-slice
recording (diagnostic only; NOT a production solver). Imports numpy, scipy.linalg, mpmath only.

Same equations and the same operation order as wz_rgf.wz_endpoint_blocks; with method="solve" the output must be bit-identical
to it (checked by the driver). Every solve X = A^{-1} B in the recursion goes through `solve_with(method)`:
    W_1-I inverse (Ghat_11), (I - W_{j-1})^{-1} W_{j-1} (W update), (I - W_j)^{-1} (P), Z_N^{-1}.
Methods (no method normalizes |P|):
    solve        numpy.linalg.solve(A, B)
    inv          numpy.linalg.inv(A) @ B
    lu           scipy.linalg.lu_solve(lu_factor(A), B)
    diagdiv      B / diag(A)[:, None]  -- scalar complex division; ONLY valid when A is exactly diagonal (checked, else error)
    refine       X0 = solve(A, B); R = B - A @ X0 (double); X = X0 + solve(A, R)
    refine_exact X0 = solve(A, B); R = B - A X0 evaluated exactly (mpmath) and rounded once to double; X = X0 + solve(A, R)
                 (1x1 only)
"""

from __future__ import annotations

import mpmath as mp
import numpy as np
import scipy.linalg as sla

METHODS = ("solve", "inv", "lu", "diagdiv", "refine", "refine_exact")


def _offdiag_max(A):
    if A.shape[0] == 1:
        return 0.0
    return float(np.max(np.abs(A - np.diag(np.diag(A)))))


def solve_with(method):
    if method == "solve":
        return lambda A, B: np.linalg.solve(A, B)
    if method == "inv":
        return lambda A, B: np.linalg.inv(A) @ B
    if method == "lu":
        return lambda A, B: sla.lu_solve(sla.lu_factor(A), B)
    if method == "diagdiv":
        def f(A, B):
            if _offdiag_max(A) != 0.0:
                raise ValueError("diagdiv requires an exactly diagonal matrix")
            return B / np.diag(A)[:, None]
        return f
    if method == "refine":
        def f(A, B):
            X0 = np.linalg.solve(A, B)
            R = B - A @ X0
            return X0 + np.linalg.solve(A, R)
        return f
    if method == "refine_exact":
        def f(A, B):
            if A.shape != (1, 1):
                raise ValueError("refine_exact is implemented for 1x1 only")
            X0 = np.linalg.solve(A, B)
            a, b, x = (mp.mpc(complex(M[0, 0]).real, complex(M[0, 0]).imag) for M in (A, B, X0))
            r = b - a * x                                  # exact at the working mpmath precision
            R = np.array([[complex(float(mp.re(r)), float(mp.im(r)))]])
            return X0 + np.linalg.solve(A, R)
        return f
    raise ValueError(method)


def wz_blocks_variant(Delta, lam_minus_1, method="solve", record=True, cond=False):
    """Delta: (N, ny, ny); lam_minus_1: (ny,). Returns Ghat blocks plus per-slice diagnostics:
       W_diag (N-1, ny) = diag(W_j) for j = 1..N-1;  P_diag (N-2, ny) = diag(P_j) for j = 2..N-1;
       G11_init = Ghat_11 after slice 1; G_Nm1_1 = Ghat_{N-1,1} before the Z_N step; Z, Zi;
       res_P (N-2,) = max|I - (I - W_j) P_j| in double; offdiag maxima; cond_IW (N-1,) if cond."""
    Delta = np.asarray(Delta)
    N, ny = Delta.shape[0], Delta.shape[1]
    S = solve_with(method)
    I = np.eye(ny, dtype=complex)
    Lm = np.diag(np.asarray(lam_minus_1, dtype=complex))
    rec = dict(method=method, N=N, ny=ny)
    if N == 1:
        Z = 2.0 * Delta[0] + 2.0 * Lm
        Zi = S(Z, I)
        rec.update(Ghat11=Zi, GhatN1=Zi.copy(), Ghat1N=Zi.copy(), GhatNN=Zi.copy(), Z=Z, Zi=Zi, W_max=0.0,
                   offdiag_max=_offdiag_max(Zi), res_Z=float(np.max(np.abs(I - Z @ Zi))), cond_Z=float(np.linalg.cond(Z)))
        return rec
    W = 2.0 * Delta[0] + Lm
    Wd, Pd, resP, condIW = [np.diag(W).copy()], [], [], []
    off = _offdiag_max(W)
    w_max = float(np.max(np.abs(W)))
    if cond:
        condIW.append(float(np.linalg.cond(I - W)))
    G11 = S(W - I, I)
    rec["G11_init"] = np.diag(G11).copy()
    Gj1 = G11.copy()
    G1j = G11.copy()
    for j in range(1, N - 1):
        W = 2.0 * Delta[j] + S(I - W, W)
        P = S(I - W, I)
        G11 = G11 - G1j @ P @ Gj1
        Gj1 = P @ Gj1
        G1j = G1j @ P
        if record:
            Wd.append(np.diag(W).copy())
            Pd.append(np.diag(P).copy())
            resP.append(float(np.max(np.abs(I - (I - W) @ P))))
            off = max(off, _offdiag_max(W), _offdiag_max(P), _offdiag_max(Gj1), _offdiag_max(G1j), _offdiag_max(G11))
        if cond:
            condIW.append(float(np.linalg.cond(I - W)))
        w_max = max(w_max, float(np.max(np.abs(W))))
    rec["G_Nm1_1"] = np.diag(Gj1).copy()
    rec["G_1_Nm1"] = np.diag(G1j).copy()
    Z = 2.0 * Delta[N - 1] + Lm + S(I - W, W)
    Zi = S(Z, I)
    G11 = G11 + G1j @ Zi @ Gj1
    GN1 = -(Zi @ Gj1)
    G1N = -(G1j @ Zi)
    off = max(off, _offdiag_max(Z), _offdiag_max(Zi), _offdiag_max(G11), _offdiag_max(GN1), _offdiag_max(G1N))
    rec.update(Ghat11=G11, GhatN1=GN1, Ghat1N=G1N, GhatNN=Zi, Z=Z, Zi=Zi, W_max=w_max, offdiag_max=off,
               W_diag=np.array(Wd), P_diag=np.array(Pd) if Pd else np.zeros((0, ny), complex), res_P=np.array(resP),
               res_Z=float(np.max(np.abs(I - Z @ Zi))), cond_Z=float(np.linalg.cond(Z)), cond_IW=np.array(condIW),
               W_last=np.diag(W).copy())
    return rec
