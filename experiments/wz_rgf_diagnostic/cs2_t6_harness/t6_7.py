"""T6.7: THRESHOLD value/type, operand rejections, comparison cases (a)-(c)."""
from fractions import Fraction as F
import numpy
from cs2_t6_reference import rpred
from .common import Results


def run(prod, R=None):
    R = R or Results("T6.7")
    TH = getattr(prod, "THRESHOLD", None)
    R.add("T6.7 THRESHOLD == Fraction(1, 10**10)", TH == F(1, 10**10), "CERT-VIOLATION", repr(TH))
    R.add("T6.7 THRESHOLD != Fraction(1e-10)", TH != F(1e-10), "CERT-VIOLATION")
    R.add("T6.7 type(THRESHOLD) is Fraction", type(TH) is F, "CERT-VIOLATION", type(TH).__name__)
    lit = rpred.THRESHOLD_REF
    for U, G, exp, tag in ((1, 10**20, True, "a"), (1, 10**5, False, "b"), (1, 10**10, False, "c (equality)")):
        try:
            got = prod.point_infeasible(U, G, lit)
        except Exception as e:  # noqa: BLE001
            got = f"raised {type(e).__name__}"
        R.add(f"T6.7 case {tag}: point_infeasible({U}, {G}) = {exp}", got is exp, "CERT-VIOLATION", got)
    for thr, w in ((1e-10, "float"), (numpy.float64(1e-10), "numpy.float64"), (F(1e-10), "Fraction(1e-10)"),
                   (F(1, 10**9), "Fraction(1, 10**9)")):
        try:
            prod.point_infeasible(1, 10**20, thr)
            ok = False
        except (TypeError, ValueError):
            ok = True
        R.add(f"T6.7 threshold {w} rejected", ok, "INTERFACE-FAIL" if w in ("float", "numpy.float64") else "CERT-VIOLATION")
    try:
        prod.point_infeasible(0.5, 10**20, lit)
        ok = False
    except TypeError:
        ok = True
    R.add("T6.7 float U rejected", ok, "INTERFACE-FAIL")
    return R
