"""C4 supplementary irrational-lead case: "numerical cross-check (not validated)" (§2.4.2 T6.4; inputs §2.4.3.8).
Its outcome never contributes to an exact predicate; only its blocking outcomes gate the manifest freeze."""
from fractions import Fraction as F

from cs2_t6_reference import rmat, rinv, rlead, rlattice, rcert, rnum
from cs2_t6_reference.qexact import RC
from . import inputs as I
from .common import Results
from .adapters import to_prod, from_prod, exact_out

BLOCKING = {"CONTRADICTED", "PRECONDITION-FAILED", "NONCONVERGENT", "REFERENCE-INVALID"}


def _prod_case(prod, lat, N, s, mid_ref):
    gam = lat["gammas"]
    Dp = [to_prod(D) for D in rlattice.delta_blocks(lat, N, s)]
    A, eps = {}, {}
    q0 = []
    for k in mid_ref:
        kind, j = ("full", None) if k == "full" else ("prefix", k[1])
        A[k] = from_prod(prod.assemble_A_mid(Dp, gam, kind, j), "A_mid")
        eps[k] = exact_out(prod.eps_A(gam, N, kind), "eps_A")
        q0.append((f"Q0 A_mid production == reference [{k}]", A[k] == mid_ref[k]))
        q0.append((f"Q0 eps_A production == reference [{k}]", eps[k] == rcert.eps_A_ref(lat, N, kind)))
    lo, hi = [], []
    for g in gam:
        pe, re_ = prod.lead_enclosure(g), rlead.enclosure(g)
        q0.append((f"Q0 mu_n, e_n production == reference [gamma={g}]", pe["mu"] == re_["mu"] and pe["e"] == re_["e"]))
        lo.append(exact_out(pe["mu"]) - exact_out(pe["e"]))
        hi.append(exact_out(pe["mu"]) + exact_out(pe["e"]))
    return A, eps, lo, hi, q0


def run(prod, R=None, quick=False):
    R = R or Results("C4 numerical cross-check (not validated)")
    bad = rnum.selection_check(I.C4_OPEN, I.C4_CLOSED, I.C4_EDGE)
    field = []
    if bad:
        field.append("REFERENCE-INVALID")
    R.add("C4 selection rules re-evaluated exactly == frozen selection", not bad, "REFERENCE-INVALID", bad)
    cases = []
    for lname, lat, Ns in I.c4_lattices():
        for N in Ns:
            for s in (I.C4_S if not quick else I.C4_S[:1]):
                cases.append((lname, lat, N, s))
    rows_all = []
    for lname, lat, N, s in cases:
        ny = len(lat["gammas"])
        mid_ref = rcert.mid_matrices(lat, N, s)
        A, eps, lo, hi, q0 = _prod_case(prod, lat, N, s, mid_ref)
        tag = f"[{lname} N={N} s={s}]"
        mutants = [m for m in I.ORDER_MUTANTS if I.on_path(m, N)] if N >= 2 else []
        st1_cache = {}
        for yname in ("Y_fine", "Y_prod"):
            if yname == "Y_fine":
                Ys = {k: rlead.round_matrix(rinv.inverse(Am), rlead.ELL_Y) for k, Am in mid_ref.items()}
            else:
                Ys = {}
                for k, Am in mid_ref.items():
                    Yq = None
                    for dps in (60, 120):
                        Yq = from_prod(prod.approx_inverse_Y(to_prod(Am), dps), "Y_prod")
                        if prod.certify_point(to_prod(Am), to_prod(Yq), eps[k])["status"] == "certified":
                            break
                    Ys[k] = Yq
            certs = {k: prod.certify_point(to_prod(A[k]), to_prod(Ys[k]), eps[k]) for k in A}
            q0y = list(q0)
            for k, c in certs.items():
                rc = rcert.cert_point_ref(mid_ref[k], Ys[k], eps[k])
                q0y.append((f"Q0 f_Xi, f_Y, f_mid production == reference [{k}]",
                            (c["f_Xi"], c["f_Y"], c["f_mid"]) == (rc["f_Xi"], rc["f_Y"], rc["f_mid"])))
            q0_fail = [n for n, ok in q0y if not ok]
            if q0_fail:
                field.append("PRECONDITION-FAILED")
                R.add(f"C4 Q0 {yname} {tag}", False, "PRECONDITION-FAILED", q0_fail)
                continue
            if any(c["status"] != "certified" for c in certs.values()):
                field.append("NONCONVERGENT")
                R.add(f"C4 {yname} certified {tag}", False, "NONCONVERGENT", "production certificate uncertified")
                continue
            prodvals = dict(A_mid=A, eps=eps, lo=lo, hi=hi, h={k: exact_out(c["h"]) for k, c in certs.items()},
                            omega={k: exact_out(c["omega"]) for k, c in certs.items()},
                            a={k: exact_out(c["a"]) for k, c in certs.items()}, eta={}, U={}, G_lo=None)
            # eta_A for every matrix and G_lo for the full matrix, independent of mutant availability (also N = 1)
            for k in A:
                prodvals["eta"][k] = exact_out(prod.eta_A_bound(to_prod(A[k]), to_prod(Ys[k]), eps[k]), "eta_A")
            prodvals["G_lo"] = exact_out(prod.G_lo_bound(to_prod(A["full"]), to_prod(Ys["full"]), eps["full"]), "G_lo")
            for m in mutants:
                refs = {k: (to_prod(A[k]), to_prod(Ys[k]), eps[k]) for k in A}
                r = prod.pointwise_certificate(m, N, ny, refs)
                if r["status"] == "computed":
                    prodvals["U"][m] = exact_out(r["U_m"], "U_m")
                    R.add(f"C4 pointwise G_lo == G_lo_bound {yname} {tag} {m}", r["G_lo"] == prodvals["G_lo"],
                          "PRECONDITION-FAILED")
                    for k, v in r["eta"].items():
                        R.add(f"C4 pointwise eta == eta_A_bound {yname} {tag} {m} {k}", v == prodvals["eta"][k],
                              "PRECONDITION-FAILED")
            outcome, rows, notes = rnum.classify_case(lat, N, s, Ys, prodvals, mutants, st1_cache)
            field.append(outcome)
            qcount = {}
            for r_ in rows:
                qcount[r_["Q"]] = qcount.get(r_["Q"], 0) + 1
            expectQ = {"Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"} | ({"Q8"} if mutants else set())
            R.add(f"C4 coverage {yname} {tag}: Q present = {sorted(qcount)}", expectQ <= set(qcount), "PRECONDITION-FAILED",
                  f"missing {sorted(expectQ - set(qcount))}")
            rows_all.append(dict(case=f"{yname} {tag}", N=N, outcome=outcome, n_rows=len(rows), qcount=qcount, notes=notes[:3],
                                 contradicted=[r_ for r_ in rows if r_["cls"] == "CONTRADICTED"][:3],
                                 unresolved=sum(1 for r_ in rows if r_["cls"] == "UNRESOLVED")))
            R.add(f"C4 {yname} {tag} not blocking", outcome not in BLOCKING, outcome, f"{outcome}; {notes[:2]}")
    fo = rnum.worst(field)
    cov = {}
    for c in rows_all:
        for q, k in c["qcount"].items():
            cov.setdefault(q, [0, 0])
            cov[q][0] += k
            cov[q][1] += k if c["N"] == 1 else 0
    return R, dict(field_outcome=fo, blocking=fo in BLOCKING, n_cases=len(rows_all), cases=rows_all,
                   coverage={q: dict(rows=v[0], rows_at_N1=v[1]) for q, v in sorted(cov.items())},
                   label="numerical cross-check (not validated)")
