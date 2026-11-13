"""T6.8: INCONCLUSIVE paths P1-P7 over S_test, including P5 value agreement and statuses; edge cases."""
from fractions import Fraction as F
from cs2_t6_reference import rmat, rpoly, rpred
from . import inputs as I
from .common import Results, ref_blocks, recursion_diff
from .adapters import check_mid, prod_refs, exact_out, from_prod
from .ychoices import Y_T6, Y_half, Y_zero
from . import p5
from .t6_1 import prod_blocks


def _run_path(R, prod, name, N, mutant, ychoice, tag, ell=256, drop=None, sets=None):
    ny = len(I.LATTICES[name]["gammas"])
    per = {}
    for s in (sets if sets is not None else I.S_TEST):
        base = check_mid(R, prod, name, N, s, ell)
        Ys = ychoice(name, N, s)
        refs = prod_refs(base, Ys)
        if drop is not None:
            refs.pop(drop, None)
        try:
            per[s] = prod.pointwise_certificate(mutant, N, ny, refs)
        except Exception as e:  # noqa: BLE001
            R.add(f"{tag} pointwise raised at s={s}", False, "PATH-FAIL", f"{type(e).__name__}: {e}")
            return None, None
    try:
        out = prod.conclude(per)
    except Exception as e:  # noqa: BLE001
        R.add(f"{tag} conclude raised", False, "PATH-FAIL", f"{type(e).__name__}: {e}")
        return per, None
    return per, out


def _expect(R, tag, out, outcomes, reason):
    if out is None:
        return
    R.add(f"{tag} outcome {outcomes} reason {reason}", out["outcome"] in outcomes and out["reason"] == reason, "PATH-FAIL",
          out)
    R.add(f"{tag} never 'infeasible'", out["outcome"] != "infeasible", "PATH-FAIL", out)


def run(prod, R=None, only=None):
    R = R or Results("T6.8")
    want = (lambda p: only is None or p in only)
    if want("P1"):
        per, out = _run_path(R, prod, "X0", 2, "M6", Y_zero, "T6.8 P1 (Y = 0)")
        if per:
            R.add("T6.8 P1 every point uncertified_reference",
                  all(r["reason"] == "uncertified_reference" for r in per.values()), "PATH-FAIL")
        _expect(R, "T6.8 P1", out, ("INCONCLUSIVE",), "uncertified_reference")
    if want("P2"):
        per, out = _run_path(R, prod, "X1", 2, "M6", Y_half, "T6.8 P2 (Y = A^-1/2, n_tot = 4)")
        if per:
            R.add("T6.8 P2 full omega >= 1 at every point",
                  all(exact_out(r.get("omega", {}).get("full", 0), "omega") >= 1 for r in per.values()), "PATH-FAIL")
        _expect(R, "T6.8 P2", out, ("INCONCLUSIVE",), "uncertified_reference")
    if want("P3"):
        per, out = _run_path(R, prod, "X1", 3, "M1", Y_T6, "T6.8 P3 (missing prefix 2)", drop=("prefix", 2))
        _expect(R, "T6.8 P3", out, ("INCONCLUSIVE",), "missing_reference")
    if want("P4"):
        y4 = (lambda name, N, s: Y_T6(name, N, s, 4))
        per, out = _run_path(R, prod, "C8", 2, "M6", y4, "T6.8 P4 (ell = 4)", ell=4)
        if per:
            R.add("T6.8 P4 production full omega >= 1",
                  all(exact_out(r.get("omega", {}).get("full", 0), "omega") >= 1 for r in per.values()), "PATH-FAIL")
        _expect(R, "T6.8 P4", out, ("INCONCLUSIVE",), "uncertified_reference")
    if want("P5"):
        per, out = _run_path(R, prod, "X1n", 2, "M6", Y_T6, "T6.8 P5 (near-threshold)")
        if per:
            d = p5.data()
            mism = [s for s in I.S_TEST if per[s].get("status") != "computed"
                    or per[s].get("U_m") != d[s]["U"] or per[s].get("G_lo") != d[s]["G_lo"]]
            R.add("T6.8 P5 value agreement U_m, G_lo == reference at all 43 points", not mism, "CERT-VIOLATION",
                  f"P5_VALUE_MISMATCH at {mism[:5]}")
            bad = [s for s in I.S_TEST if per[s].get("point") != d[s]["status"]]
            R.add("T6.8 P5 point statuses == reference prediction (5 not_infeasible / 38 infeasible)", not bad,
                  "PATH-FAIL", f"{bad[:5]}")
            R.add("T6.8 P5 s* = 2^-19 not_infeasible", per[F(1, 2**19)].get("point") == "not_infeasible", "PATH-FAIL")
        _expect(R, "T6.8 P5", out, ("INCONCLUSIVE",), "bound_not_below_threshold")
    if want("P6"):
        for m in I.UNSUPPORTED_P6:
            per, out = _run_path(R, prod, "X0", 2, m, Y_T6, f"T6.8 P6 ({m})")
            if per:
                R.add(f"T6.8 P6 ({m}) every point reason unsupported_mutant",
                      all(r["status"] in ("INCONCLUSIVE", "not applicable") and r["reason"] == "unsupported_mutant"
                          for r in per.values()), "PATH-FAIL")
            _expect(R, f"T6.8 P6 ({m})", out, ("INCONCLUSIVE", "not applicable"), "unsupported_mutant")
    if want("P7"):
        for sname, sset in I.P7_SETS.items():
            per, out = _run_path(R, prod, "X0", 2, "M6", Y_T6, f"T6.8 P7 ({sname})", sets=sset)
            _expect(R, f"T6.8 P7 ({sname})", out, ("INCONCLUSIVE",), "evaluated_set_mismatch")
    if want("edge"):
        for name in ("X0", "XD"):
            for N in I.NS[name]:
                for s in I.S6:
                    b = ref_blocks(name, N, s)
                    for m in I.ORDER_MUTANTS:
                        if I.on_path(m, N):
                            z = from_prod(prod.F_poly(m, N, prod_blocks(b, N)))
                            R.add(f"edge F_poly(B) = 0 [{name} N={N} s={s} {m}]", rmat.is_zero(z), "CERT-VIOLATION")
        for name in I.MAIN:
            for s in I.S6:
                b = ref_blocks(name, 3, s)
                R.add(f"edge M5 j=2 term vanishes [{name} N=3 s={s}]",
                      rmat.is_zero(from_prod(prod.F_poly("M5", 3, prod_blocks(b, 3)))), "CERT-VIOLATION")
        base = check_mid(R, prod, "X1", 2, F(1, 2))
        r = prod.pointwise_certificate("M8", 2, 2, prod_refs(base, Y_T6("X1", 2, F(1, 2))))
        R.add("edge M8 at N = 2 is 'not applicable'", r["status"] == "not applicable", "PATH-FAIL", r)
        R.add("edge production S_TEST == S_TEST_REF (43 points)",
              getattr(prod, "S_TEST", None) == rpred.S_TEST_REF and len(rpred.S_TEST_REF) == 43, "CERT-VIOLATION")
    return R
