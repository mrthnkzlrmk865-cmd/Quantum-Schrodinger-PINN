"""T6.5: U_m >= 0 and U_m^2 >= D*_true^2 (Y_T6), plus aux Y_adv (X1, N = 2, s = 2^-19, M6)."""
from fractions import Fraction as F
from cs2_t6_reference import rmat
from .common import Results, recursion_diff, ref_certified, required_keys, expect_status
from .ychoices import Y_T6
from .adapters import check_mid, prod_refs, exact_out
from .ychoices import Y_adv
from . import pointwise_cache


def _check(R, res, D2, tag, ref_ok=True):
    if not expect_status(R, res, ref_ok, f"T6.5 {tag}", accept_bounds=True):
        return
    U = exact_out(res["U_m"], "U_m")
    R.add(f"T6.5 U_m >= 0 {tag}", U >= 0, "CERT-VIOLATION")
    R.add(f"T6.5 U_m^2 >= D*^2 {tag}", U * U >= D2, "CERT-VIOLATION")


def run(prod, R=None, cache=None):
    R = R or Results("T6.5")
    cache = cache if cache is not None else pointwise_cache.compute(R, prod)
    for (name, N, s, m), res in cache.items():
        _check(R, res, rmat.fro2(recursion_diff(name, N, s, m)), f"[{name} N={N} s={s} {m}]",
               ref_certified(name, N, s, required_keys(m, N), Y_T6(name, N, s)))
    s = F(1, 2**19)
    base = check_mid(R, prod, "X1", 2, s)
    res = prod.pointwise_certificate("M6", 2, 2, prod_refs(base, Y_adv()))
    D2 = rmat.fro2(recursion_diff("X1", 2, s, "M6"))
    R.add("T6.5 aux witness D*_true > 0 [X1 N=2 s=2^-19 M6 Y_adv]", D2 > 0, "PRECONDITION-UNMET")
    _check(R, res, D2, "[aux Y_adv X1 N=2 s=2^-19 M6]", ref_certified("X1", 2, s, required_keys("M6", 2), Y_adv()))
    return R
