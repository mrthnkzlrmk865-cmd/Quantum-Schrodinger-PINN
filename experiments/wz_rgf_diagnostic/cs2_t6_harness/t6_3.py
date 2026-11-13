"""T6.3: certified, omega < 1, ||A_true^-1||_2 <= a, for every full and prefix matrix (Y_T6), plus aux Y_half."""
from fractions import Fraction as F
from cs2_t6_reference import rpsd
from .common import Results, ytrue, all_cert_cases, ref_certified
from .adapters import check_mid, exact_out, to_prod
from .ychoices import Y_T6, Y_half


def _one(R, prod, base, Ys, name, N, s, tag):
    inv = ytrue(name, N, s)
    for key, (Ap, e) in base.items():
        c = prod.certify_point(Ap, to_prod(Ys[key]), e)
        t = f"{tag} [{name} N={N} s={s} {key}]"
        ref_ok = ref_certified(name, N, s, [key], Ys)
        if not R.add(f"T6.3 certified {t}", c["status"] == "certified",
                     "CERT-VIOLATION" if ref_ok else "PRECONDITION-UNMET", f"reference certified={ref_ok}"):
            continue
        om, a = exact_out(c["omega"], "omega"), exact_out(c["a"], "a")
        R.add(f"T6.3 omega < 1 {t}", om < 1, "CERT-VIOLATION")
        R.add(f"T6.3 ||A^-1||_2 <= a {t}", rpsd.norm2_le(inv[key], a), "CERT-VIOLATION")


def run(prod, R=None, cases=None):
    R = R or Results("T6.3")
    for name, N, s in (cases or list(all_cert_cases())):
        base = check_mid(R, prod, name, N, s)
        _one(R, prod, base, Y_T6(name, N, s), name, N, s, "Y_T6")
    base = check_mid(R, prod, "X0", 2, F(1, 2))
    _one(R, prod, base, Y_half("X0", 2, F(1, 2)), "X0", 2, F(1, 2), "aux Y_half")
    return R
