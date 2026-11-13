"""Exact predicates: THRESHOLD_REF literal, the §3 grid, S_TEST_REF, bound checks, and the point_infeasible contract."""
from fractions import Fraction
from .qexact import exact_rational

THRESHOLD_REF = Fraction(1, 10**10)


class InvalidLowerBound(ValueError):
    code = "INVALID_LOWER_BOUND"


class InvalidUpperBound(ValueError):
    code = "INVALID_UPPER_BOUND"


def grid():
    """§3: {2^-k : k = 0..40} U {3 * 2^-(k+2) : k = 0..39} (81 points)."""
    g = {Fraction(1, 2**k) for k in range(41)} | {Fraction(3, 2**(k + 2)) for k in range(40)}
    assert len(g) == 81
    return frozenset(g)


S_TEST_REF = frozenset(s for s in grid() if s <= Fraction(1, 2**19))
_EXPLICIT = frozenset({Fraction(1, 2**k) for k in range(19, 41)} | {Fraction(3, 2**k) for k in range(21, 42)})
assert S_TEST_REF == _EXPLICIT and len(S_TEST_REF) == 43


def check_lower_bound(G):
    G = exact_rational(G, "G_lo")
    if G < 0:
        raise InvalidLowerBound("INVALID_LOWER_BOUND")
    return "ZERO" if G == 0 else "POSITIVE"


def check_upper_bound(U):
    U = exact_rational(U, "U")
    if U < 0:
        raise InvalidUpperBound("INVALID_UPPER_BOUND")
    return "ZERO" if U == 0 else "POSITIVE"


def expected_point_infeasible(U, G, threshold):
    """The five-step contract of §2.4.3.3, with the reference's own checks. Returns True/False or raises."""
    for x, w in ((U, "U"), (G, "G_lo"), (threshold, "threshold")):
        exact_rational(x, w)
    if Fraction(threshold) != THRESHOLD_REF:
        raise ValueError("threshold must equal Fraction(1, 10**10)")
    check_upper_bound(U)
    check_lower_bound(G)
    return Fraction(U) < Fraction(threshold) * Fraction(G)
