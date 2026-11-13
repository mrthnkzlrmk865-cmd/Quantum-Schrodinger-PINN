"""Exact reference blocks from direct inverses: G_hat = t*A^-1 (full) and the prefix inverses."""
from .qexact import RC
from . import rmat, rinv, rlattice


def blocks(lat, N, s, lam):
    ny = len(lat["gammas"])
    A = rlattice.M_phys(lat, N, s, "full", lam=lam)
    Ainv = rinv.inverse(A)
    Gh = rmat.scale(RC(rlattice.T), Ainv)
    out = {"A": A, "Ainv": Ainv, "Ghat": Gh,
           "G_NN": rmat.block(Gh, N, N, ny), "G_N1": rmat.block(Gh, N, 1, ny),
           "G_1N": rmat.block(Gh, 1, N, ny), "G_11": rmat.block(Gh, 1, 1, ny), "prefix": {}}
    for j in range(1, N):
        Ap = rlattice.M_phys(lat, N, s, "prefix", j=j, lam=lam)
        Apinv = rinv.inverse(Ap)
        Gp = rmat.scale(RC(rlattice.T), Apinv)
        out["prefix"][j] = {"A": Ap, "Ainv": Apinv, "Ghat": Gp}
        out[("P", j)] = rmat.neg(rmat.block(Gp, j, j, ny))
        out[("chi", j)] = rmat.block(Gp, j, 1, ny)
        out[("psi", j)] = rmat.block(Gp, 1, j, ny)
        if j == 1:
            out["G1_11"] = rmat.block(Gp, 1, 1, ny)
    return out
