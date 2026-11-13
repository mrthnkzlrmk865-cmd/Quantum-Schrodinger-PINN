"""Exact matrices: tuples of row tuples of RC."""
from fractions import Fraction
from .qexact import RC, ZERO, ONE


def mat(rows):
    rows = tuple(tuple(RC.of(x) for x in r) for r in rows)
    if not rows or any(len(r) != len(rows[0]) for r in rows):
        raise ValueError("ragged or empty matrix")
    return rows


def shape(A):
    return len(A), len(A[0])


def eye(n):
    return tuple(tuple(ONE if i == j else ZERO for j in range(n)) for i in range(n))


def zeros(n, m=None):
    m = n if m is None else m
    return tuple(tuple(ZERO for _ in range(m)) for _ in range(n))


def matmul(A, B):
    if len(A[0]) != len(B):
        raise ValueError("shape mismatch")
    Bt = tuple(zip(*B))
    out = []
    for r in A:
        row = []
        for c in Bt:
            re = Fraction(0)
            im = Fraction(0)
            for x, y in zip(r, c):
                re += x.re * y.re - x.im * y.im
                im += x.re * y.im + x.im * y.re
            row.append(RC(re, im))
        out.append(tuple(row))
    return tuple(out)


def add(A, B):
    return tuple(tuple(x + y for x, y in zip(ra, rb)) for ra, rb in zip(A, B))


def sub(A, B):
    return tuple(tuple(x - y for x, y in zip(ra, rb)) for ra, rb in zip(A, B))


def scale(c, A):
    c = RC.of(c)
    return tuple(tuple(c * x for x in r) for r in A)


def neg(A):
    return tuple(tuple(-x for x in r) for r in A)


def ctranspose(A):
    return tuple(tuple(A[i][j].conj() for i in range(len(A))) for j in range(len(A[0])))


def transpose(A):
    return tuple(tuple(A[i][j] for i in range(len(A))) for j in range(len(A[0])))


def block(X, a, b, ny):
    """Block (a, b), 1-indexed slices: rows (a-1)ny..a*ny-1, columns (b-1)ny..b*ny-1."""
    return tuple(tuple(r[(b - 1) * ny:b * ny]) for r in X[(a - 1) * ny:a * ny])


def fro2(X):
    return sum((x.abs2() for r in X for x in r), Fraction(0))


def is_zero(X):
    return all(x.is_zero() for r in X for x in r)


def diag(vals):
    n = len(vals)
    return tuple(tuple(RC.of(vals[i]) if i == j else ZERO for j in range(n)) for i in range(n))


def max_abs2(X):
    return max(x.abs2() for r in X for x in r)
