"""T6.6: G_lo precondition / frozen-input invariant (0 < G_lo) and G_lo <= ||G_hat||_2; aux Y_half; interface cases a-i."""
from fractions import Fraction as F
import numpy
import mpmath
from cs2_t6_reference import rmat, rpsd, rpred
from cs2_t6_reference.qexact import RC
from . import inputs as I
from .common import Results, ref_blocks, ref_certified, required_keys, expect_status
from .ychoices import Y_T6
from .adapters import check_mid, prod_refs, InterfaceFail
from .ychoices import Y_half
from . import pointwise_cache


class Counted:
    """Counts calls to the reference matrix predicate (case c requires zero calls)."""
    calls = 0

    @classmethod
    def is_pd(cls, H):
        cls.calls += 1
        return rpsd.is_pd(H)


def _gap(g, Gh):
    GG = rmat.matmul(rmat.ctranspose(Gh), Gh)
    n = len(GG)
    return tuple(tuple((RC(g * g) if i == j else RC(0)) - GG[i][j] for j in range(n)) for i in range(n))


def check_glo(R, G_lo, Gh, tag):
    """Order: exact type -> sign (frozen invariant requires POSITIVE) -> matrix predicate."""
    if isinstance(G_lo, bool) or not isinstance(G_lo, (int, F)):
        R.add(f"T6.6 G_lo exact type {tag}", False, "INTERFACE-FAIL", type(G_lo).__name__)
        return
    try:
        kind = rpred.check_lower_bound(G_lo)
    except rpred.InvalidLowerBound:
        R.add(f"T6.6 G_lo >= 0 {tag}", False, "CERT-VIOLATION", "INVALID_LOWER_BOUND")
        return
    if kind == "ZERO":
        R.add(f"T6.6 G_lo > 0 on frozen input {tag}", False, "CERT-VIOLATION", "ZERO_LOWER_BOUND")
        return
    R.add(f"T6.6 G_lo <= ||G_hat||_2 {tag}", not Counted.is_pd(_gap(F(G_lo), Gh)), "CERT-VIOLATION")


def interface_cases(R, prod):
    K = rmat.mat(I.K)
    G = I.G_IFACE
    TH = rpred.THRESHOLD_REF

    def prod_call(U, Gl):
        try:
            return ("value", prod.point_infeasible(U, Gl, TH))
        except TypeError:
            return ("TypeError", None)
        except ValueError as e:
            return ("ValueError", str(e))
        except Exception as e:  # noqa: BLE001
            return (type(e).__name__, str(e))

    def ref_call(U, Gl):
        try:
            return ("value", rpred.expected_point_infeasible(U, Gl, TH))
        except TypeError:
            return ("TypeError", None)
        except (rpred.InvalidUpperBound, rpred.InvalidLowerBound) as e:
            return ("ValueError", e.code)
        except ValueError as e:
            return ("ValueError", str(e))

    def agree(tag, U, Gl, expected):
        r, p = ref_call(U, Gl), prod_call(U, Gl)
        R.add(f"T6.6 iface {tag} reference == expected", r == expected, "KAT-FAIL", f"{r} vs {expected}")
        R.add(f"T6.6 iface {tag} production == expected", p == expected, "CERT-VIOLATION", f"{p} vs {expected}")

    # a: valid positive G_lo (matrix predicate with K; 3 and 5 hold, control 51/10 flags a bound above ||K||_2)
    for g, holds in ((F(3), True), (F(5), True), (F(51, 10), False)):
        R.add(f"T6.6 iface a: check_lower_bound({g}) = POSITIVE", rpred.check_lower_bound(g) == "POSITIVE", "KAT-FAIL")
        R.add(f"T6.6 iface a: predicate G_lo={g} holds={holds}", (not rpsd.is_pd(_gap(g, K))) is holds, "KAT-FAIL")
    agree("a (1, G)", 1, G, ("value", True))
    agree("a (3, G)", 3, G, ("value", False))
    # b: zero G_lo (generic contract: admissible, never infeasible)
    R.add("T6.6 iface b: check_lower_bound(0) = ZERO", rpred.check_lower_bound(0) == "ZERO", "KAT-FAIL")
    for U in (0, 1, 10**20):
        agree(f"b ({U}, 0)", U, 0, ("value", False))
    # c: negative G_lo, rejected before the matrix predicate (zero calls)
    before = Counted.calls
    for g in (-1, -6):
        ok = False
        try:
            check_glo_probe = rpred.check_lower_bound(g)  # noqa: F841
        except rpred.InvalidLowerBound:
            ok = True
        R.add(f"T6.6 iface c: G_lo={g} rejected (INVALID_LOWER_BOUND)", ok, "KAT-FAIL")
        agree(f"c (1, {g})", 1, g, ("ValueError", "INVALID_LOWER_BOUND"))
    R.add("T6.6 iface c: zero calls to the matrix predicate", Counted.calls == before, "KAT-FAIL")
    # d: non-exact G_lo
    for bad, w in ((3.0, "float"), (numpy.float64(3), "numpy.float64"), (mpmath.mpf(3), "mpmath.mpf"), (True, "bool")):
        agree(f"d (1, {w})", 1, bad, ("TypeError", None))
    # e/f: U = 0 and U > 0
    agree("e (0, G)", 0, G, ("value", True))
    for U, exp in ((1, True), (3, False), (4, False)):
        agree(f"f ({U}, G)", U, G, ("value", exp))
    # g: U < 0 (checked before G_lo)
    agree("g (-1, G)", -1, G, ("ValueError", "INVALID_UPPER_BOUND"))
    agree("g (-1, 0)", -1, 0, ("ValueError", "INVALID_UPPER_BOUND"))
    # h/i: float and bool U
    for bad, w in ((1.0, "float"), (numpy.float64(1), "numpy.float64"), (mpmath.mpf(1), "mpmath.mpf"), (True, "True"),
                   (False, "False")):
        agree(f"h/i ({w}, G)", bad, G, ("TypeError", None))


def run(prod, R=None, cache=None):
    R = R or Results("T6.6")
    cache = cache if cache is not None else pointwise_cache.compute(R, prod)
    for (name, N, s, m), res in cache.items():
        if not expect_status(R, res, ref_certified(name, N, s, required_keys(m, N), Y_T6(name, N, s)),
                             f"T6.6 [{name} N={N} s={s} {m}]", accept_bounds=True):
            continue
        check_glo(R, res["G_lo"], ref_blocks(name, N, s)["Ghat"], f"[{name} N={N} s={s} {m}]")
    s = F(1, 2)
    base = check_mid(R, prod, "X0", 2, s)
    res = prod.pointwise_certificate("M6", 2, 1, prod_refs(base, Y_half("X0", 2, s)))
    if expect_status(R, res, ref_certified("X0", 2, s, required_keys("M6", 2), Y_half("X0", 2, s)),
                     "T6.6 aux [X0 N=2 s=1/2 M6 Y_half]", accept_bounds=True):
        check_glo(R, res["G_lo"], ref_blocks("X0", 2, s)["Ghat"], "[aux Y_half X0 N=2 s=1/2 M6]")
    interface_cases(R, prod)
    return R
