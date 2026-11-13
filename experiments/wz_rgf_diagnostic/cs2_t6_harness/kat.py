"""Stage 2: reference known-answer tests and self-proofs (§2.4.3.7), including N1'. Must pass before production is read."""
from fractions import Fraction as F
import numpy  # noqa: F401  (only to build the numpy rejection case; the reference itself never imports numpy)
import mpmath

from cs2_t6_reference import rmat, rinv, rpsd, rlead, rpoly, rpred
from cs2_t6_reference.qexact import RC
from cs2_t6_reference.rinv import SingularMatrix
from . import inputs as I
from .common import Results, ref_blocks, recursion, recursion_diff, t6_cases
from . import p5


def _raises(fn, exc):
    try:
        fn()
    except exc:
        return True
    except Exception:
        return False
    return False


def run():
    R = Results("2 reference KATs and self-proofs")
    m = rmat.mat
    # RC
    R.add("RC product", RC(F(1, 2), F(1, 3)) * RC(2, -1) == RC(F(4, 3), F(1, 6)), "KAT-FAIL")
    R.add("RC quotient", RC(F(4, 3), F(1, 6)) / RC(2, -1) == RC(F(1, 2), F(1, 3)), "KAT-FAIL")
    R.add("|3/5+4i/5|^2 = 1", RC(F(3, 5), F(4, 5)).abs2() == 1, "KAT-FAIL")
    for bad, w in ((0.5, "float"), (1 + 0j, "complex"), (numpy.float64(1), "numpy.float64"), (True, "bool"),
                   (mpmath.mpf(1), "mpmath")):
        R.add(f"RC rejects {w}", _raises(lambda: RC(bad), TypeError), "KAT-FAIL")
    # rinv
    R.add("inv [[2,1],[1,1]]", rinv.inverse(m([[2, 1], [1, 1]])) == m([[1, -1], [-1, 2]]), "KAT-FAIL")
    R.add("inv permutation (pivoting)", rinv.inverse(m([[0, 1], [1, 0]])) == m([[0, 1], [1, 0]]), "KAT-FAIL")
    R.add("inv complex upper-triangular",
          rinv.inverse(m([[(1, 1), F(1, 2)], [0, (1, -1)]])) == m([[(F(1, 2), F(-1, 2)), F(-1, 4)], [0, (F(1, 2), F(1, 2))]]),
          "KAT-FAIL")
    R.add("inv singular raises", _raises(lambda: rinv.inverse(m([[1, 2], [2, 4]])), SingularMatrix), "KAT-FAIL")
    # rpsd
    for H, exp in (([[1, 0], [0, 0]], True), ([[0, 0], [0, 0]], True), ([[1, 1], [1, 1]], True),
                   ([[2, (0, 1)], [(0, -1), 2]], True), ([[1, (1, 1)], [(1, -1), 2]], True),
                   ([[0, 1], [1, 0]], False), ([[1, 2], [2, 1]], False), ([[1, (1, 1)], [(1, -1), 1]], False),
                   ([[-1]], False)):
        R.add(f"is_psd {H} = {exp}", rpsd.is_psd(m(H)) is exp, "KAT-FAIL")
    R.add("is_psd non-Hermitian raises", _raises(lambda: rpsd.is_psd(m([[1, 1], [0, 1]])), rpsd.NotHermitian), "KAT-FAIL")
    for H, exp in (([[1, 1], [1, 1]], False), ([[2, (0, 1)], [(0, -1), 2]], True), ([[1, 0], [0, 0]], False)):
        R.add(f"is_pd {H} = {exp}", rpsd.is_pd(m(H)) is exp, "KAT-FAIL")
    for X, c, exp in (([[3, 0], [0, 4]], F(4), True), ([[3, 0], [0, 4]], F(399, 100), False),
                      ([[3, 0], [4, 0]], F(5), True), ([[3, 0], [4, 0]], F(49, 10), False)):
        R.add(f"norm2_le {X} {c} = {exp}", rpsd.norm2_le(m(X), c) is exp, "KAT-FAIL")
    for bad, w in ((4.0, "float"), (numpy.float64(4), "numpy.float64"), (True, "bool"), (mpmath.mpf(4), "mpmath.mpf")):
        R.add(f"norm2_le rejects c of type {w}", _raises(lambda: rpsd.norm2_le(m([[3, 0], [0, 4]]), bad), TypeError), "KAT-FAIL")
        R.add(f"norm2_lt rejects c of type {w}", _raises(lambda: rpsd.norm2_lt(m([[3, 0], [0, 4]]), bad), TypeError), "KAT-FAIL")
    R.add("norm2_le accepts int and Fraction c", rpsd.norm2_le(m([[3, 0], [0, 4]]), 4) is True
          and rpsd.norm2_le(m([[3, 0], [0, 4]]), F(399, 100)) is False and rpsd.norm2_lt(m([[3, 0], [0, 4]]), F(401, 100)), "KAT-FAIL")
    R.add("fro2 KAT = 15", rmat.fro2(m([[(1, 1), 2], [0, (0, -3)]])) == 15, "KAT-FAIL")
    # rlead
    for g, lam in ((F(5, 4), RC(F(1, 2))), (F(5, 3), RC(F(1, 3))), (F(29, 20), RC(F(2, 5))), (F(3, 5), RC(F(3, 5), F(4, 5)))):
        R.add(f"exact_root({g})", rlead.exact_root(g) == lam, "KAT-FAIL")
    R.add("exact_root irrational raises", _raises(lambda: rlead.exact_root(F(1, 2)), ValueError), "KAT-FAIL")
    enc = rlead.enclosure(F(5, 4), 256)  # xi = 9/16
    R.add("enclosure(xi=9/16) = (3/4, 3/4 + 2^-256)", (enc["lo"], enc["hi"]) == (F(3, 4), F(3, 4) + F(1, 2**256)), "KAT-FAIL")
    R.add("round_dyadic KAT", [rlead.round_dyadic(x, 2) for x in (F(1, 8), F(3, 8), F(5, 8), F(-3, 8), F(1, 3), F(7, 8))]
          == [0, F(1, 2), F(1, 2), F(-1, 2), F(1, 4), 1], "KAT-FAIL")
    R.add("ub(0)=0, ub(f)^2 > f", rlead.ub(0) == 0 and all(rlead.ub(f) ** 2 > f for f in (F(1), F(2), F(1, 4), F(1, 2**600))),
          "KAT-FAIL")
    R.add("lb(f)^2 <= f", all(rlead.lb(f) ** 2 <= f for f in (F(0), F(1), F(2), F(1, 4), F(7, 3))), "KAT-FAIL")
    # rnum exact mpf conversion keeps the sign (mpmath 1.3.0 mpf.man_exp drops it; found during C4 debugging)
    from cs2_t6_reference import rnum
    ctx = mpmath.MPContext()
    ctx.dps = 60
    R.add("rnum exact mpf conversion: -3, -3/4, 0, 5, 2^-1000", [rnum._exact_mpf(ctx, ctx.mpf(v)) for v in (-3, -0.75, 0, 5)]
          == [F(-3), F(-3, 4), F(0), F(5)] and rnum._exact_mpf(ctx, ctx.mpf(2) ** -1000) == F(1, 2**1000), "KAT-FAIL")
    for sp, w in ((ctx.nan, "NaN"), (ctx.inf, "+inf"), (-ctx.inf, "-inf")):
        R.add(f"rnum exact mpf conversion rejects {w}", _raises(lambda: rnum._exact_mpf(ctx, sp), ValueError), "KAT-FAIL")
    # rpred
    R.add("S_TEST_REF = 43 explicit points = grid <= 2^-19", len(rpred.S_TEST_REF) == 43
          and rpred.S_TEST_REF == frozenset(I.S_TEST), "KAT-FAIL")
    R.add("THRESHOLD_REF literal", rpred.THRESHOLD_REF == F(1, 10**10) and type(rpred.THRESHOLD_REF) is F, "KAT-FAIL")
    # self-proofs on every rational-lead T6 lattice
    for name in I.MAIN + I.TARGETED + ("X1n",):
        for N in I.NS[name]:
            for s in (I.S6 if name != "X1n" else [F(1, 2**19), F(1, 2**40)]):
                b, rec = ref_blocks(name, N, s), recursion(name, N, s, None)
                ok = (rec["G11"] == b["G_11"] and rec["GN1"] == b["G_N1"] and rec["G1N"] == b["G_1N"]
                      and rec["GNN"] == b["G_NN"])
                R.add(f"rmutants(None) == rblocks [{name} N={N} s={s}]", ok, "REF-SELFCHECK-FAIL")
    for name, N, s, mu in t6_cases():
        b = ref_blocks(name, N, s)
        R.add(f"F_ref == recursion diff [{name} N={N} s={s} {mu}]",
              rpoly.F_ref(mu, N, b) == recursion_diff(name, N, s, mu), "REF-SELFCHECK-FAIL")
    # P5 self-proofs (N1', U_ref^2 >= D*^2, G_lo_ref <= G_true at every point)
    for item in p5.selfproofs():
        R.add(*item)
    return R
