"""Physical lattice assembly M = E*I - H - (Sigma_L (+) Sigma_R), built from H explicitly (never from the Delta form).

A lattice spec is a dict: gammas (list of exact rationals, channel order), V (list of N_y x N_y matrices V_j(1), entries
exact rationals or (re, im) pairs). E := 0 and E_perp,n := -delta_n/dx^2 with delta_n = 1 - gamma_n, so that the clean-lead
condition delta_n = (E - E_perp,n) dx^2 holds exactly (§2.4.3.1).
"""
from fractions import Fraction
from .qexact import RC, exact_rational
from . import rmat

T = Fraction(32)
DX2 = Fraction(1, 64)


def H_blocks(lat, n_slices, s):
    s = exact_rational(s, "s")
    ny = len(lat["gammas"])
    E_perp = [-(1 - exact_rational(g)) / DX2 for g in lat["gammas"]]
    diag_blocks = []
    for j in range(n_slices):
        V = rmat.mat(lat["V"][j])
        Hjj = tuple(tuple((RC(2 * T + E_perp[a]) if a == b else RC(0)) + RC(s) * V[a][b] for b in range(ny))
                    for a in range(ny))
        diag_blocks.append(Hjj)
    return diag_blocks


def M_phys(lat, N, s, kind="full", j=None, lam=None):
    """Full matrix (kind='full', N slices, both leads) or prefix j (kind='prefix', slices 1..j, left lead only).
    lam: list of RC lead roots per channel (exact, midpoint, or C4 centre)."""
    ny = len(lat["gammas"])
    n = N if kind == "full" else j
    if kind not in ("full", "prefix") or n is None or n < 1 or (kind == "prefix" and not 1 <= j <= N - 1):
        raise ValueError("bad kind/j")
    E = Fraction(0)
    Hd = H_blocks(lat, n, s)
    size = n * ny
    M = [[RC(0)] * size for _ in range(size)]
    for b in range(n):
        for a1 in range(ny):
            for a2 in range(ny):
                M[b * ny + a1][b * ny + a2] = (RC(E) if a1 == a2 else RC(0)) - Hd[b][a1][a2]
        if b + 1 < n:
            for a in range(ny):  # H_{j,j+1} = -t I  =>  -H = +t I
                M[b * ny + a][(b + 1) * ny + a] = RC(T)
                M[(b + 1) * ny + a][b * ny + a] = RC(T)
    sigma = [-(RC(T) * RC.of(l)) for l in lam]  # Sigma = -t lambda
    for a in range(ny):
        M[a][a] = M[a][a] - sigma[a]                       # Sigma_L on slice 1
    if kind == "full":
        for a in range(ny):
            k = (n - 1) * ny + a
            M[k][k] = M[k][k] - sigma[a]                   # Sigma_R on slice N (same block when N = 1)
    return tuple(tuple(r) for r in M)


def delta_blocks(lat, N, s):
    """Delta_j(s) = diag(1 - gamma_n) - s V_j(1) dx^2 (W/Z variables; used only to drive the transcribed recursion)."""
    s = exact_rational(s, "s")
    ny = len(lat["gammas"])
    out = []
    for j in range(N):
        V = rmat.mat(lat["V"][j])
        out.append(tuple(tuple((RC(1 - exact_rational(lat["gammas"][a])) if a == b else RC(0)) - RC(s * DX2) * V[a][b]
                               for b in range(ny)) for a in range(ny)))
    return out
