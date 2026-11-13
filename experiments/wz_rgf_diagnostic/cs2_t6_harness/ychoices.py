"""Frozen Y matrices of §2.4.3.1 (built from reference data only)."""
from fractions import Fraction as F
from cs2_t6_reference import rmat
from cs2_t6_reference.qexact import RC
from .common import ytrue, yt6, ref_blocks
from . import inputs as I


def Y_T6(name, N, s, ell=256):
    return yt6(name, N, s, ell)


def Y_exact(name, N, s):
    return ytrue(name, N, s)


def Y_pert(name, N, s):
    out = {}
    for k, Y in ytrue(name, N, s).items():
        n = len(Y)
        out[k] = tuple(tuple(Y[i][j] + (RC(F(1, 2**40)) if (i, j) == (0, 0) else RC(0)) for j in range(n)) for i in range(n))
    return out


def Y_half(name, N, s):
    return {k: rmat.scale(RC(F(1, 2)), Y) for k, Y in ytrue(name, N, s).items()}


def Y_zero(name, N, s, ell=256):
    return {k: rmat.zeros(len(Y)) for k, Y in yt6(name, N, s, ell).items()}


def Y_adv():
    """X1, N = 2, s = 2^-19, M6: block (2,1) of the full Y replaced so that F_M6(B~) = 0 exactly."""
    name, N, s = "X1", 2, F(1, 2**19)
    ny = 2
    yt = ytrue(name, N, s)
    Yp = yt[("prefix", 1)]
    Yf = [list(r) for r in yt["full"]]
    repl = rmat.scale(RC(-I.T), rmat.matmul(Yp, rmat.block(yt["full"], 2, 2, ny)))
    for a in range(ny):
        for b in range(ny):
            Yf[ny + a][b] = repl[a][b]
    return {"full": tuple(tuple(r) for r in Yf), ("prefix", 1): Yp}
