"""P5 reference computations (§2.4.3.6): N1-N5, N1', reference point statuses, and U_ref/G_lo_ref self-checks.
Everything here uses the clean-room reference only."""
from fractions import Fraction as F
from functools import lru_cache

from cs2_t6_reference import rmat, rpsd, rcert, rpred
from cs2_t6_reference.qexact import RC
from . import inputs as I
from .common import ref_blocks, recursion_diff, mids, yt6, eps

NAME, N, M, NY = "X1n", 2, "M6", 2
TH = rpred.THRESHOLD_REF


def refs_ref(s):
    A, Y = mids(NAME, N, s), yt6(NAME, N, s)
    return {"full": (A["full"], Y["full"], eps(NAME, N, "full")),
            ("prefix", 1): (A[("prefix", 1)], Y[("prefix", 1)], eps(NAME, N, "prefix"))}


def _gram(s):
    Gh = ref_blocks(NAME, N, s)["Ghat"]
    return rmat.matmul(rmat.ctranspose(Gh), Gh)


def _rho_ge(s, c):
    """rho(s) >= c  <=>  (D*^2/(c T)^2) I - G^H G is PSD."""
    D2 = rmat.fro2(recursion_diff(NAME, N, s, M))
    GG = _gram(s)
    k = D2 / ((c * TH) ** 2)
    n = len(GG)
    return rpsd.is_psd(tuple(tuple((RC(k) if i == j else RC(0)) - GG[i][j] for j in range(n)) for i in range(n)))


@lru_cache(maxsize=None)
def data():
    out = {}
    for s in I.S_TEST:
        r = rcert.pointwise_ref(M, N, NY, refs_ref(s))
        out[s] = dict(U=r["U"], G_lo=r["G_lo"], D2=rmat.fro2(recursion_diff(NAME, N, s, M)),
                      status="infeasible" if rpred.expected_point_infeasible(r["U"], r["G_lo"], TH) else "not_infeasible")
    return out


def selfproofs():
    d = data()
    s19 = F(1, 2**19)
    items = [("P5 N1': U_ref(2^-19) >= THRESHOLD*G_lo_ref(2^-19)", d[s19]["U"] >= TH * d[s19]["G_lo"],
              "REF-SELFCHECK-FAIL", "")]
    for s in I.S_TEST:
        items.append((f"P5 U_ref^2 >= D*^2 [s={s}]", d[s]["U"] ** 2 >= d[s]["D2"], "REF-SELFCHECK-FAIL", ""))
        GG = _gram(s)
        g2 = d[s]["G_lo"] ** 2
        H = tuple(tuple((RC(g2) if i == j else RC(0)) - GG[i][j] for j in range(len(GG))) for i in range(len(GG)))
        items.append((f"P5 G_lo_ref <= G_true [s={s}]", not rpsd.is_pd(H), "REF-SELFCHECK-FAIL", ""))
    return items


def preconditions():
    d = data()
    s19, s40 = F(1, 2**19), F(1, 2**40)
    items = [("P5 N1: rho(2^-19) >= 1", _rho_ge(s19, F(1)), "PRECONDITION-UNMET", "")]
    items.append(("P5 N2: rho(s) < 5/2 for all s in S_test", all(not _rho_ge(s, F(5, 2)) for s in I.S_TEST),
                  "PRECONDITION-UNMET", ""))
    items.append(("P5 N3: rho(2^-40) < 1/5", not _rho_ge(s40, F(1, 5)), "PRECONDITION-UNMET", ""))
    n4 = [s for s in I.S_TEST if not d[s]["U"] < 10 * TH * d[s]["G_lo"]]
    items.append(("P5 N4: U_ref < 10*THRESHOLD*G_lo_ref at all 43 points", not n4, "PRECONDITION-UNMET", f"violations {n4}"))
    items.append(("P5 N5: U_ref(2^-40) < THRESHOLD*G_lo_ref(2^-40)", d[s40]["U"] < TH * d[s40]["G_lo"],
                  "PRECONDITION-UNMET", ""))
    ni = frozenset(s for s in I.S_TEST if d[s]["status"] == "not_infeasible")
    items.append(("P5 reference status list == recorded list (5 not_infeasible / 38 infeasible)",
                  ni == I.P5_NOT_INFEASIBLE, "PRECONDITION-UNMET", f"reference not_infeasible = {sorted(ni)}"))
    mx = max(d[s]["U"] / (TH * d[s]["G_lo"]) for s in I.S_TEST)
    items.append(("P5 info: max U_ref/(THRESHOLD*G_lo_ref)", True, None, f"{float(mx):.6f}"))
    return items
