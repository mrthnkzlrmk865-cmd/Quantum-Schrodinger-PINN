r"""
EXPERIMENTAL small-quantity W/Z formulation of the mode-space RGF (diagnostic only; NOT a production solver).

Standalone: imports only numpy. It does not import or modify stationary.py / smatrix_spec.py.

Same discrete problem as stationary.rgf_endpoint_blocks (block-tridiagonal M = E - H_D - Sigma_L (+) Sigma_R,
M_{j,j+-1} = +t I, t = 1/(2 dx^2), Sigma = diag(-t lambda), convention P), rewritten in the variables
    Delta_j = (E I - diag(E_perp) - V_j) dx^2           (so A_j = t (2 Delta_j - 2 I) exactly, since 2 t dx^2 = 1)
    W_j     = D_j / t + I       (1 <= j < N),           Z_N = D_N / t,
    Ghat    = t G.
Audited equations (docs: experiments/wz_rgf_diagnostic/REPORT.md section 1):
    N >= 2:  W_1 = 2 Delta_1 + (Lambda - I)
             W_j = 2 Delta_j + (I - W_{j-1})^{-1} W_{j-1}                      1 < j < N
             Z_N = 2 Delta_N + (Lambda - I) + (I - W_{N-1})^{-1} W_{N-1}
    N == 1:  Z_1 = 2 Delta_1 + 2 (Lambda - I)
    Ghat_11 <- (W_1 - I)^{-1};  Ghat_{j1} <- (I - W_j)^{-1} Ghat_{j-1,1}  (LEFT);  Ghat_{1j} <- Ghat_{1,j-1} (I - W_j)^{-1}  (RIGHT);
             Ghat_11 <- Ghat_11 + Ghat_{1,j-1} (W_j - I)^{-1} Ghat_{j-1,1}  (old values)
    j = N:   Ghat_N1 = -Z^{-1} Ghat_{N-1,1};  Ghat_1N = -Ghat_{1,N-1} Z^{-1};  Ghat_NN = Z^{-1};  Ghat_11 += Ghat_{1,N-1} Z^{-1} Ghat_{N-1,1}
    Fisher-Lee (open channels, s = sin theta = Gamma/(2t)):  t_RL = 2i sqrt(s_m s_n) Ghat_N1, r_L = 2i sqrt(s s) Ghat_11 - I, ...
"""

from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------
# lead data from delta = (E - E_perp) dx^2, cancellation-free (open and closed below threshold)
# ----------------------------------------------------------------------
def lead_from_delta(delta):
    """Returns dict(lam_minus_1, sin_theta, open) for each channel, computed directly from delta:
       open   (0 < delta < 2):  lambda - 1 = -delta + i sqrt(delta (2 - delta)),   sin theta = sqrt(delta (2 - delta))
       closed (delta < 0):      lambda - 1 = -2|delta| / (|delta| + sqrt(|delta| (2 + |delta|))),  sin theta = 0
       delta >= 2 (c <= -1, band-top artifact) and delta == 0 (threshold) are rejected."""
    d = np.asarray(delta, dtype=float)
    if np.any(d >= 2.0) or np.any(d == 0.0):
        raise ValueError("delta outside the supported range (threshold or c <= -1)")
    lm1 = np.empty(d.shape, dtype=complex)
    s = np.zeros(d.shape)
    op = d > 0
    so = np.sqrt(d[op] * (2.0 - d[op]))
    lm1[op] = -d[op] + 1j * so
    s[op] = so
    ad = -d[~op]
    lm1[~op] = -2.0 * ad / (ad + np.sqrt(ad * (2.0 + ad)))
    return dict(lam_minus_1=lm1, sin_theta=s, open=op)


def lead_from_code_data(lam, gamma, t):
    """Lead data supplied as numbers (e.g. the existing double lead data): Lambda - I by one double subtraction,
    sin theta = Gamma / (2t) (exact division by a power of two when dx = 2^-p)."""
    lam = np.asarray(lam, dtype=complex)
    g = np.asarray(gamma, dtype=float)
    return dict(lam_minus_1=lam - 1.0, sin_theta=g / (2.0 * t), open=g > 0)


# ----------------------------------------------------------------------
# W/Z recursion
# ----------------------------------------------------------------------
def wz_endpoint_blocks(Delta, lam_minus_1, record_w=False):
    """Delta: (N, ny, ny) real symmetric; lam_minus_1: (ny,) complex (Lambda - I, diagonal).
    Returns dict(Ghat11, GhatN1, Ghat1N, GhatNN, W_max, Z, w_trace)."""
    Delta = np.asarray(Delta)
    N, ny = Delta.shape[0], Delta.shape[1]
    I = np.eye(ny, dtype=complex)
    Lm = np.diag(np.asarray(lam_minus_1, dtype=complex))
    trace = []
    if N == 1:
        Z = 2.0 * Delta[0] + 2.0 * Lm
        Zi = np.linalg.solve(Z, I)
        return dict(Ghat11=Zi, GhatN1=Zi.copy(), Ghat1N=Zi.copy(), GhatNN=Zi.copy(), W_max=0.0, Z=Z, w_trace=trace)
    W = 2.0 * Delta[0] + Lm                                        # W_1
    if record_w:
        trace.append(W.copy())
    w_max = float(np.max(np.abs(W)))
    G11 = np.linalg.solve(W - I, I)                                # Ghat_11 = (W_1 - I)^{-1}
    Gj1 = G11.copy()                                               # Ghat_{1,1} as the start of both chains
    G1j = G11.copy()
    for j in range(1, N - 1):                                      # slices 2 .. N-1
        W = 2.0 * Delta[j] + np.linalg.solve(I - W, W)             # W_j = 2 Delta_j + (I - W_{j-1})^{-1} W_{j-1}
        if record_w:
            trace.append(W.copy())
        w_max = max(w_max, float(np.max(np.abs(W))))
        P = np.linalg.solve(I - W, I)                              # (I - W_j)^{-1};  (W_j - I)^{-1} = -P
        G11 = G11 - G1j @ P @ Gj1                                  # += Ghat_{1,j-1} (W_j - I)^{-1} Ghat_{j-1,1}  (old values)
        Gj1 = P @ Gj1                                              # LEFT multiplication
        G1j = G1j @ P                                              # RIGHT multiplication
    Z = 2.0 * Delta[N - 1] + Lm + np.linalg.solve(I - W, W)        # Z_N
    Zi = np.linalg.solve(Z, I)
    G11 = G11 + G1j @ Zi @ Gj1
    GN1 = -(Zi @ Gj1)
    G1N = -(G1j @ Zi)
    return dict(Ghat11=G11, GhatN1=GN1, Ghat1N=G1N, GhatNN=Zi, W_max=w_max, Z=Z, w_trace=trace)


def wz_fisher_lee(G, sin_theta, is_open):
    """Flux-normalized S = [[r_L, t_LR], [t_RL, r_R]] over open channels from the scaled Green blocks."""
    o = np.where(is_open)[0]
    sq = np.sqrt(np.asarray(sin_theta)[o])
    F = 2j * sq[:, None] * sq[None, :]
    no = o.size
    I = np.eye(no)
    ix = np.ix_(o, o)
    t_RL = F * G["GhatN1"][ix]
    t_LR = F * G["Ghat1N"][ix]
    r_L = F * G["Ghat11"][ix] - I
    r_R = F * G["GhatNN"][ix] - I
    return dict(S=np.block([[r_L, t_LR], [t_RL, r_R]]), t_RL=t_RL, t_LR=t_LR, r_L=r_L, r_R=r_R, n_open=no)


def original_scalar_recurrence(a, sigma, t, N, record_d=False):
    """Scalar replica of the ORIGINAL double recurrence of stationary.rgf_endpoint_blocks for one decoupled channel:
       D = a (= fl(E - H_jj));  D -= sigma at j = 1 and j = N;  D -= (t*t) * g  for j > 1;  g = 1 / D.
    Used only for the fixed-point trajectory; the production solver itself is run separately."""
    d_trace = []
    g = G_j1 = G_1j = G11 = None
    for j in range(N):
        D = complex(a)
        if j == 0:
            D = D - sigma
        if j == N - 1:
            D = D - sigma
        if j > 0:
            D = D - (t * t) * g
        if record_d:
            d_trace.append(D)
        gn = 1.0 / D
        if j == 0:
            g = gn
            G_j1 = G_1j = G11 = gn
            continue
        gG = gn * G_j1
        G11 = G11 + (t * t) * (G_1j * gG)
        G_1j = -t * (G_1j * gn)
        G_j1 = -t * gG
        g = gn
    return dict(G11=G11, GN1=G_j1, G1N=G_1j, GNN=g, d_trace=d_trace)
