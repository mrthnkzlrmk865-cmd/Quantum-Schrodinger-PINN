r"""
EXPERIMENTAL (diagnostic only; NOT a production solver): instrumented copy of wz_rgf.wz_endpoint_blocks and the mutants
M1-M16 (incl. M4L, M4R) of PREDICTIONS_coupled_blocks.md section 6. numpy only.

With mutant=None the operations, their order and their operands are those of wz_rgf.wz_endpoint_blocks, so the four endpoint
blocks must be bit-identical (checked by the driver for every execution). Each mutant changes exactly one expression.

Trace names, in computation order:
  N = 1:  Z, Zi, then the endpoint blocks G11, GN1, G1N, GNN
  N >= 2: W_1, G11^(1); for every loop slice s = 2..N-1: W_s, P_s, G11^(s), G_s1^(s), G_1s^(s); then Z, Zi, G11, GN1, G1N, GNN
"""

from __future__ import annotations

import numpy as np

MUTANTS = ("M1", "M2", "M3", "M4", "M4L", "M4R", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14", "M15", "M16")
ORDER_MUTANTS = ("M1", "M2", "M3", "M5", "M6", "M7", "M8", "M15")          # pure multiplication-order mutants
MIN_N = {"M1": 3, "M2": 3, "M3": 3, "M4": 3, "M4L": 3, "M4R": 3, "M5": 3, "M9": 3, "M13": 3,
         "M6": 2, "M7": 2, "M8": 2, "M10": 2, "M12": 2, "M14": 2, "M15": 2, "M16": 2}
PREDICTED_FIRST = {"M1": "G_21^(2)", "M2": "G_12^(2)", "M3": "G_21^(2)", "M4": "G11^(2)", "M4L": "G11^(2)", "M4R": "G11^(2)",
                   "M5": "G11^(3)", "M6": "GN1", "M7": "G1N", "M8": "G11", "M9": "P_2", "M10": "Z", "M11": "Z", "M12": None,
                   "M13": "G11^(2)", "M14": "G11^(1)", "M15": "GN1", "M16": "G11"}


def on_path(m, N):
    if m == "M11":
        return N == 1
    return N >= MIN_N[m]


def _wpart(I, W, m):
    """(I - W)^{-1} W; M12: W (I - W)^{-1} by right division (plain transposes)."""
    if m == "M12":
        return np.linalg.solve((I - W).T, W.T).T
    return np.linalg.solve(I - W, W)


def wz_instrumented(Delta, lam_minus_1, mutant=None):
    m = mutant
    if m is not None and m not in MUTANTS:
        raise ValueError(m)
    Delta = np.asarray(Delta)
    N, ny = Delta.shape[0], Delta.shape[1]
    I = np.eye(ny, dtype=complex)
    Lm = np.diag(np.asarray(lam_minus_1, dtype=complex))
    trace = []

    def rec(name, X):
        trace.append((name, X.copy()))

    if N == 1:
        Z = 2.0 * Delta[0] + Lm if m == "M11" else 2.0 * Delta[0] + 2.0 * Lm
        rec("Z", Z)
        Zi = np.linalg.solve(Z, I)
        rec("Zi", Zi)
        out = dict(Ghat11=Zi, GhatN1=Zi.copy(), Ghat1N=Zi.copy(), GhatNN=Zi.copy())
        for nm, k in (("G11", "Ghat11"), ("GN1", "GhatN1"), ("G1N", "Ghat1N"), ("GNN", "GhatNN")):
            rec(nm, out[k])
        out.update(trace=trace, Z=Z, Zi=Zi)
        return out
    W = 2.0 * Delta[0] + Lm
    rec("W_1", W)
    G11 = np.linalg.solve(I - W, I) if m == "M14" else np.linalg.solve(W - I, I)
    rec("G11^(1)", G11)
    Gj1 = G11.copy()
    G1j = G11.copy()
    for j in range(1, N - 1):
        s = j + 1
        if m == "M9":
            P = np.linalg.solve(I - W, I)                          # stale: from W_{s-1}
            W = 2.0 * Delta[j] + _wpart(I, W, m)
        else:
            W = 2.0 * Delta[j] + _wpart(I, W, m)
            P = np.linalg.solve(I - W, I)
        rec(f"W_{s}", W)
        rec(f"P_{s}", P)
        if m == "M4":
            Gj1 = P @ Gj1
            G1j = G1j @ P
            G11 = G11 - G1j @ P @ Gj1
        elif m == "M4L":
            Gj1 = P @ Gj1
            G11 = G11 - G1j @ P @ Gj1
            G1j = G1j @ P
        elif m == "M4R":
            G1j = G1j @ P
            G11 = G11 - G1j @ P @ Gj1
            Gj1 = P @ Gj1
        else:
            if m == "M5":
                G11 = G11 - Gj1 @ P @ G1j
            elif m == "M13":
                G11 = G11 + G1j @ P @ Gj1
            else:
                G11 = G11 - G1j @ P @ Gj1
            Gj1 = Gj1 @ P if m in ("M1", "M3") else P @ Gj1
            G1j = P @ G1j if m in ("M2", "M3") else G1j @ P
        rec(f"G11^({s})", G11)
        rec(f"G_{s}1^({s})", Gj1)
        rec(f"G_1{s}^({s})", G1j)
    if m == "M10":
        Z = 2.0 * Delta[N - 1] + _wpart(I, W, m)
    else:
        Z = 2.0 * Delta[N - 1] + Lm + _wpart(I, W, m)
    rec("Z", Z)
    Zi = np.linalg.solve(Z, I)
    rec("Zi", Zi)
    if m == "M16":
        GN1 = -(Zi @ Gj1)
        G1N = -(G1j @ Zi)
        G11 = G11 + G1N @ Zi @ GN1
    else:
        G11 = G11 + (Gj1 @ Zi @ G1j if m == "M8" else G1j @ Zi @ Gj1)
        GN1 = -(Gj1 @ Zi) if m in ("M6", "M15") else -(Zi @ Gj1)
        G1N = -(Zi @ G1j) if m in ("M7", "M15") else -(G1j @ Zi)
    out = dict(Ghat11=G11, GhatN1=GN1, Ghat1N=G1N, GhatNN=Zi)
    for nm, k in (("G11", "Ghat11"), ("GN1", "GhatN1"), ("G1N", "Ghat1N"), ("GNN", "GhatNN")):
        rec(nm, out[k])
    out.update(trace=trace, Z=Z, Zi=Zi)
    return out
