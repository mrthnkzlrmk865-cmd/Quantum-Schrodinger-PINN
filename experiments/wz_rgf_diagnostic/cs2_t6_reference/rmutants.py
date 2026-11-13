"""Exact transcription of wz_mutants.wz_instrumented (source text read, never imported) in exact rationals.

Line correspondence (wz_mutants.py at 63e046d; reviewed in R0):
  wz_mutants L32-36  _wpart                  -> _wpart (M12: plain transposes, W (I-W)^-1)
  L52-61   N == 1 branch (M11)               -> N == 1 branch
  L62      W = 2 Delta_0 + Lm                -> W_1
  L64      G11 init (M14: solve(I - W, I))   -> G11 init
  L68-75   loop W/P update (M9 stale P)      -> loop W/P
  L78-89   M4, M4L, M4R chain/G11 orderings  -> same
  L91-98   M5, M13, default G11; M1/M2/M3 chains -> same
  L102-105 closing Z (M10 drops Lm)          -> closing Z
  L109-116 M16, M8, M6/M15, M7/M15 closing    -> same
solve(X, B) is X^-1 B with an exact, proved inverse.
"""
from .qexact import RC
from . import rmat, rinv

MUTANTS = ("M1", "M2", "M3", "M4", "M4L", "M4R", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M12", "M13", "M14",
           "M15", "M16")


def _solve(X, B):
    return rmat.matmul(rinv.inverse(X), B)


def _wpart(I, W, m):
    if m == "M12":  # solve((I - W).T, W.T).T  ==  W (I - W)^-1
        return rmat.transpose(_solve(rmat.transpose(rmat.sub(I, W)), rmat.transpose(W)))
    return _solve(rmat.sub(I, W), W)


def endpoint_blocks(Delta, lam_minus_1, mutant=None):
    m = mutant
    if m is not None and m not in MUTANTS:
        raise ValueError(m)
    N, ny = len(Delta), len(Delta[0])
    I = rmat.eye(ny)
    Lm = rmat.diag([RC.of(x) for x in lam_minus_1])
    two = RC(2)
    if N == 1:
        Z = rmat.add(rmat.scale(two, Delta[0]), Lm) if m == "M11" else rmat.add(rmat.scale(two, Delta[0]),
                                                                                rmat.scale(two, Lm))
        Zi = _solve(Z, I)
        return dict(G11=Zi, GN1=Zi, G1N=Zi, GNN=Zi)
    W = rmat.add(rmat.scale(two, Delta[0]), Lm)
    G11 = _solve(rmat.sub(I, W), I) if m == "M14" else _solve(rmat.sub(W, I), I)
    Gj1 = G11
    G1j = G11
    for j in range(1, N - 1):
        if m == "M9":
            P = _solve(rmat.sub(I, W), I)
            W = rmat.add(rmat.scale(two, Delta[j]), _wpart(I, W, m))
        else:
            W = rmat.add(rmat.scale(two, Delta[j]), _wpart(I, W, m))
            P = _solve(rmat.sub(I, W), I)
        mm = rmat.matmul
        if m == "M4":
            Gj1 = mm(P, Gj1)
            G1j = mm(G1j, P)
            G11 = rmat.sub(G11, mm(mm(G1j, P), Gj1))
        elif m == "M4L":
            Gj1 = mm(P, Gj1)
            G11 = rmat.sub(G11, mm(mm(G1j, P), Gj1))
            G1j = mm(G1j, P)
        elif m == "M4R":
            G1j = mm(G1j, P)
            G11 = rmat.sub(G11, mm(mm(G1j, P), Gj1))
            Gj1 = mm(P, Gj1)
        else:
            if m == "M5":
                G11 = rmat.sub(G11, mm(mm(Gj1, P), G1j))
            elif m == "M13":
                G11 = rmat.add(G11, mm(mm(G1j, P), Gj1))
            else:
                G11 = rmat.sub(G11, mm(mm(G1j, P), Gj1))
            Gj1 = mm(Gj1, P) if m in ("M1", "M3") else mm(P, Gj1)
            G1j = mm(P, G1j) if m in ("M2", "M3") else mm(G1j, P)
    if m == "M10":
        Z = rmat.add(rmat.scale(two, Delta[N - 1]), _wpart(I, W, m))
    else:
        Z = rmat.add(rmat.add(rmat.scale(two, Delta[N - 1]), Lm), _wpart(I, W, m))
    Zi = _solve(Z, I)
    mm = rmat.matmul
    if m == "M16":
        GN1 = rmat.neg(mm(Zi, Gj1))
        G1N = rmat.neg(mm(G1j, Zi))
        G11 = rmat.add(G11, mm(mm(G1N, Zi), GN1))
    else:
        G11 = rmat.add(G11, mm(mm(Gj1, Zi), G1j) if m == "M8" else mm(mm(G1j, Zi), Gj1))
        GN1 = rmat.neg(mm(Gj1, Zi)) if m in ("M6", "M15") else rmat.neg(mm(Zi, Gj1))
        G1N = rmat.neg(mm(Zi, G1j)) if m in ("M7", "M15") else rmat.neg(mm(G1j, Zi))
    return dict(G11=G11, GN1=GN1, G1N=G1N, GNN=Zi)
