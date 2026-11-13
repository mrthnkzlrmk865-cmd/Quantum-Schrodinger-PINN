"""Exact PSD / PD decisions for Hermitian matrices (§2.4.3.9), and exact 2-norm comparisons."""
from .qexact import RC, exact_rational
from . import rmat


class NotHermitian(ValueError):
    pass


def _check_hermitian(H):
    n, m = rmat.shape(H)
    if n != m or any(H[i][j] != H[j][i].conj() for i in range(n) for j in range(n)):
        raise NotHermitian("matrix is not Hermitian")


def is_psd(H):
    """Diagonal pivoting. Negative diagonal -> not PSD; all-zero diagonal -> PSD iff remaining block is 0;
    otherwise pivot on the largest positive diagonal entry (ties: smallest index) and recurse on the Schur complement."""
    _check_hermitian(H)
    S = [list(r) for r in H]
    while S:
        n = len(S)
        d = [S[i][i].re for i in range(n)]
        if any(x < 0 for x in d):
            return False
        p = max(range(n), key=lambda i: (d[i], -i))
        if d[p] == 0:
            return all(x.is_zero() for r in S for x in r)
        hp = S[p][p]
        rest = [i for i in range(n) if i != p]
        S = [[S[i][j] - S[i][p] * S[p][j] / hp for j in rest] for i in rest]
    return True


def is_pd(H):
    """Sylvester-type: H > 0 iff H11 > 0 and the Schur complement is PD."""
    _check_hermitian(H)
    S = [list(r) for r in H]
    while S:
        h = S[0][0]
        if h.re <= 0:
            return False
        S = [[S[i][j] - S[i][0] * S[0][j] / h for j in range(1, len(S))] for i in range(1, len(S))]
    return True


def _gram_gap(X, c):
    c = exact_rational(c, "norm bound c")  # strict: int (not bool) or Fraction only
    if c < 0:
        raise ValueError("c must be >= 0")
    G = rmat.matmul(rmat.ctranspose(X), X)
    n = len(G)
    return tuple(tuple((RC(c * c) if i == j else RC(0)) - G[i][j] for j in range(n)) for i in range(n))


def norm2_le(X, c):
    """||X||_2 <= c  <=>  c^2 I - X^H X is PSD (c >= 0)."""
    return is_psd(_gram_gap(X, c))


def norm2_lt(X, c):
    """||X||_2 < c  <=>  c^2 I - X^H X is PD (c >= 0)."""
    return is_pd(_gram_gap(X, c))
