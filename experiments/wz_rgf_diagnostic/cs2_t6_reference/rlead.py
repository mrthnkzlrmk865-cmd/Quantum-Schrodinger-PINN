"""Lead roots, §2.4.1 isqrt enclosures (reimplemented), rigorous square roots and ties-to-even dyadic rounding."""
from fractions import Fraction
from math import isqrt
from .qexact import RC, exact_rational

ELL = 256
ELL_Y = 512


def _floor(q):
    return q.numerator // q.denominator


def _ceil(q):
    return -((-q.numerator) // q.denominator)


def ub(f, ell=ELL):
    """Rigorous upper square root (§2.4.1 step 4): ub(0) = 0; ub(f) = (isqrt(ceil(f 4^ell)) + 1)/2^ell."""
    f = exact_rational(f, "ub argument")
    if f < 0:
        raise ValueError("ub of a negative number")
    if f == 0:
        return Fraction(0)
    return Fraction(isqrt(_ceil(f * 4 ** ell)) + 1, 2 ** ell)


def lb(f, ell=ELL):
    """Rigorous lower square root (§2.4.1 step 4): lb(f) = isqrt(floor(f 4^ell))/2^ell."""
    f = exact_rational(f, "lb argument")
    if f < 0:
        raise ValueError("lb of a negative number")
    return Fraction(isqrt(_floor(f * 4 ** ell)), 2 ** ell)


def branch_and_xi(gamma):
    g = exact_rational(gamma, "gamma")
    if abs(g) < 1:
        return "open", 1 - g * g
    if g > 1:
        return "closed", g * g - 1
    raise ValueError("gamma <= -1 or |gamma| = 1 is excluded (§2.4.1)")


def _rational_sqrt(q):
    a, b = q.numerator, q.denominator
    ra, rb = isqrt(a), isqrt(b)
    if ra * ra != a or rb * rb != b:
        return None
    return Fraction(ra, rb)


def is_rational_square(q):
    return _rational_sqrt(Fraction(q)) is not None


def exact_root(gamma):
    """Retarded lead root for rational sqrt(xi) only; proves lambda + 1/lambda = 2 gamma, the branch and Im > 0."""
    g = exact_rational(gamma, "gamma")
    br, xi = branch_and_xi(g)
    r = _rational_sqrt(xi)
    if r is None:
        raise ValueError(f"sqrt(xi) is irrational for gamma = {g}")
    lam = RC(g - r) if br == "closed" else RC(g, r)
    if lam + RC(1) / lam != RC(2 * g):
        raise AssertionError("lambda + 1/lambda != 2 gamma")
    if br == "closed" and not (0 < lam.re < 1):
        raise AssertionError("closed root not in (0, 1)")
    if br == "open" and not (lam.abs2() == 1 and lam.im > 0):
        raise AssertionError("open root not on the unit circle with Im > 0")
    return lam


def enclosure(gamma, ell=ELL):
    """§2.4.1: lo = isqrt(floor(xi 4^ell))/2^ell, hi = (isqrt(ceil(xi 4^ell)) + 1)/2^ell, mu = (lo+hi)/2, e = (hi-lo)/2;
    the midpoint root is gamma - mu (closed) or gamma + i mu (open)."""
    g = exact_rational(gamma, "gamma")
    br, xi = branch_and_xi(g)
    a = xi * 4 ** ell
    lo = Fraction(isqrt(_floor(a)), 2 ** ell)
    hi = Fraction(isqrt(_ceil(a)) + 1, 2 ** ell)
    mu, e = (lo + hi) / 2, (hi - lo) / 2
    lam_mid = RC(g - mu) if br == "closed" else RC(g, mu)
    return dict(lo=lo, hi=hi, mu=mu, e=e, branch=br, xi=xi, lam_mid=lam_mid)


def round_dyadic(x, ell_y=ELL_Y):
    """Nearest, ties-to-even, on the grid 2^-ell_y; integer-only."""
    x = exact_rational(x, "x")
    v = x * (1 << ell_y)
    p, q = v.numerator, v.denominator
    n0 = p // q
    rem = p - n0 * q
    if 2 * rem < q:
        n = n0
    elif 2 * rem > q:
        n = n0 + 1
    else:
        n = n0 if n0 % 2 == 0 else n0 + 1
    return Fraction(n, 1 << ell_y)


def round_matrix(X, ell_y=ELL_Y):
    return tuple(tuple(RC(round_dyadic(z.re, ell_y), round_dyadic(z.im, ell_y)) for z in r) for r in X)
