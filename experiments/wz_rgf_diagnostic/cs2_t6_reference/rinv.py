"""Exact inverse by Gauss-Jordan elimination on [A | I] with first-nonzero pivoting, returned only with a proof."""
from .qexact import ZERO
from . import rmat


class SingularMatrix(ArithmeticError):
    pass


def inverse(A):
    n, m = rmat.shape(A)
    if n != m:
        raise ValueError("square matrix required")
    I = rmat.eye(n)
    M = [list(A[i]) + list(I[i]) for i in range(n)]
    for c in range(n):
        p = next((r for r in range(c, n) if not M[r][c].is_zero()), None)
        if p is None:
            raise SingularMatrix("no nonzero pivot")
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and not M[r][c].is_zero():
                f = M[r][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    X = tuple(tuple(row[n:]) for row in M)
    if rmat.matmul(A, X) != I or rmat.matmul(X, A) != I:
        raise SingularMatrix("inverse proof failed (A.X != I or X.A != I)")
    return X
