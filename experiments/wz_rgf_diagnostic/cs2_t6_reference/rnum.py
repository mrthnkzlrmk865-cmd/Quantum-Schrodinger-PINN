"""C4 supplementary irrational-lead numerical cross-check (§2.4.2 T6.4, protocol of Revision 3.5; inputs §2.4.3.8).

Label: "numerical cross-check (not validated)". Stage 1 uses mpmath at dps 300/400/500 (the only mpmath use in the
reference); Stage 2 is exact. Assumption A1 (faithful rounding of mp.sqrt) is NOT verified, which is why nothing here can
contribute to a validated claim. Own ub/lb at ell_ref(p) = 2*prec(p).
"""
from fractions import Fraction
from math import isqrt

import mpmath

from .qexact import RC
from . import rmat, rlead, rlattice, rpoly

T = rlattice.T
DPS = (300, 400, 500)
EXPECTED_PREC = {300: 1000, 400: 1332, 500: 1664}
ORDER = ("CONTRADICTED", "PRECONDITION-FAILED", "NONCONVERGENT", "REFERENCE-INVALID", "UNRESOLVED", "CONSISTENT")


def worst(outcomes):
    return min(outcomes, key=ORDER.index) if outcomes else "CONSISTENT"


def _exact_mpf(ctx, x):
    """Exact dyadic value of an mpf from its raw (sign, mantissa, exponent, bitcount) tuple. NaN/inf are rejected.
    (mpf.man_exp is NOT used: in mpmath 1.3.0 it returns |mantissa| and drops the sign.)"""
    sign, man, exp, _ = ctx.mpf(x)._mpf_
    if man == 0 and exp != 0:
        raise ValueError("non-finite mpmath value")
    mag = Fraction(man) * Fraction(2) ** exp if exp >= 0 else Fraction(man, 2 ** (-exp))
    return -mag if sign else mag


def _to_mp(ctx, X):
    M = ctx.matrix(len(X), len(X[0]))
    for i, r in enumerate(X):
        for j, z in enumerate(r):
            M[i, j] = ctx.mpc(ctx.mpf(z.re.numerator) / z.re.denominator, ctx.mpf(z.im.numerator) / z.im.denominator)
    return M


def _from_mp(ctx, M, n, m):
    return tuple(tuple(RC(_exact_mpf(ctx, ctx.re(M[i, j])), _exact_mpf(ctx, ctx.im(M[i, j]))) for j in range(m))
                 for i in range(n))


def lead_lams(gammas, s_hat):
    out = []
    for g, sh in zip(gammas, s_hat):
        br, _ = rlead.branch_and_xi(g)
        out.append(RC(g - sh) if br == "closed" else RC(g, sh))
    return out


def stage1(lat, N, s, p):
    """mpmath at dps p: s_hat_n, X_hat and X_hat' (lu_solve) for every matrix, and the svd test vector of each X_hat."""
    ctx = mpmath.MPContext()
    ctx.dps = p
    prec = ctx.prec
    info = dict(dps=p, prec=prec, invalid=[])
    if prec != EXPECTED_PREC[p]:
        info["invalid"].append(f"prec {prec} != expected {EXPECTED_PREC[p]}")
    s_hat = []
    for g in lat["gammas"]:
        _, xi = rlead.branch_and_xi(g)
        xm = ctx.mpf(xi.numerator) / xi.denominator
        if _exact_mpf(ctx, xm) != xi:
            info["invalid"].append(f"xi import not exact for gamma={g}")
        s_hat.append(_exact_mpf(ctx, ctx.sqrt(xm)))
    lam = lead_lams(lat["gammas"], s_hat)
    mats = {"full": rlattice.M_phys(lat, N, s, "full", lam=lam)}
    for j in range(1, N):
        mats[("prefix", j)] = rlattice.M_phys(lat, N, s, "prefix", j=j, lam=lam)
    X, X2, v = {}, {}, {}
    for k, A in mats.items():
        n = len(A)
        Am = _to_mp(ctx, A)
        X[k] = _from_mp(ctx, ctx.inverse(Am), n, n)
        cols = [ctx.lu_solve(Am, ctx.matrix([[1 if i == c else 0] for i in range(n)])) for c in range(n)]
        Mx = ctx.matrix(n, n)
        for c in range(n):
            for i in range(n):
                Mx[i, c] = cols[c][i]
        X2[k] = _from_mp(ctx, Mx, n, n)
        try:
            Um, S, V = ctx.svd_c(_to_mp(ctx, X[k]))
            kk = max(range(n), key=lambda i: S[i])
            v[k] = tuple(RC(_exact_mpf(ctx, ctx.re(V[kk, i])), -_exact_mpf(ctx, ctx.im(V[kk, i]))) for i in range(n))
        except Exception:  # noqa: BLE001  (fallback: basis vector of the column with the largest norm)
            col = max(range(n), key=lambda c: sum(X[k][i][c].abs2() for i in range(n)))
            v[k] = tuple(RC(1 if i == col else 0) for i in range(n))
        if all(z.is_zero() for z in v[k]):
            col = max(range(n), key=lambda c: sum(X[k][i][c].abs2() for i in range(n)))
            v[k] = tuple(RC(1 if i == col else 0) for i in range(n))
    return dict(info=info, s_hat=s_hat, A_c=mats, X=X, X2=X2, v=v)


def _lead_entries(key, N, ny):
    if key == "full":
        if N == 1:
            return [(p, 2) for p in range(ny)]
        return [(p, 1) for p in range(ny)] * 2
    return [(p, 1) for p in range(ny)]


def _iv(lo, hi):
    return (max(Fraction(0), lo), hi)


def stage2(lat, N, s, st1, Ys, prod, mutants):
    """Exact Stage 2 at one precision: reference intervals I_p(Q) for every Q and matrix (both inverse methods)."""
    p = st1["info"]["dps"]
    prec = st1["info"]["prec"]
    ell = 2 * prec
    u = Fraction(1, 2 ** (prec - 1))
    ny = len(lat["gammas"])
    ub = (lambda f: rlead.ub(f, ell))
    lb = (lambda f: rlead.lb(f, ell))
    r_n = [u * sh / (1 - u) for sh in st1["s_hat"]]
    out = {}
    for n_, g in enumerate(lat["gammas"]):
        out[("Q1", n_)] = {"A": (st1["s_hat"][n_] - r_n[n_], st1["s_hat"][n_] + r_n[n_])}
    rho, rX = {}, {}
    for k, Ac in st1["A_c"].items():
        rho[k] = ub(sum((Fraction(c * c) * T * T * r_n[q] ** 2 for q, c in _lead_entries(k, N, ny)), Fraction(0)))
        Amid = prod["A_mid"][k]
        Aec = rmat.sub(Ac, Amid)
        out[("Q2", k)] = {"A": _iv(lb(rmat.fro2(Aec)) - rho[k], ub(rmat.fro2(Aec)) + rho[k])}
        out[("Q3", k)] = {"A": _iv(lb(rmat.fro2(Ac)) - rho[k], ub(rmat.fro2(Ac)) + rho[k])}
        Y = Ys[k]
        n = len(Ac)
        Rc = rmat.sub(rmat.eye(n), rmat.matmul(Ac, Y))
        fY = ub(rmat.fro2(Y))
        out[("Q4", k)] = {"A": _iv(lb(rmat.fro2(Rc)) - rho[k] * fY, ub(rmat.fro2(Rc)) + rho[k] * fY)}
        for meth, Xs in (("A", st1["X"]), ("B", st1["X2"])):
            Xh = Xs[k]
            Xi2 = rmat.sub(rmat.eye(n), rmat.matmul(Ac, Xh))
            fX = ub(rmat.fro2(Xh))
            om = ub(rmat.fro2(Xi2)) + rho[k] * fX
            if om >= 1:
                out.setdefault("nonconvergent", []).append(f"omega' >= 1 [{k} method {meth} dps {p}]")
                continue
            r = fX * om / (1 - om)
            rX[(k, meth)] = r
            out.setdefault(("Q5", k), {})[meth] = _iv(lb(rmat.fro2(Xh)) - r, fX + r)
            D = rmat.sub(Xh, Y)
            out.setdefault(("Q6", k), {})[meth] = _iv(T * (lb(rmat.fro2(D)) - r), T * (ub(rmat.fro2(D)) + r))
            if k == "full":
                vv = st1["v"][k]
                Xv = [sum((Xh[i][c] * vv[c] for c in range(n)), RC(0)) for i in range(n)]
                q = sum((z.abs2() for z in Xv), Fraction(0)) / sum((z.abs2() for z in vv), Fraction(0))
                lo7 = T * max(lb(rmat.max_abs2(Xh)), lb(q)) - T * r
                out.setdefault(("Q7", k), {})[meth] = (lo7, T * (fX + r))
    for m in mutants:
        for meth, Xs in (("A", st1["X"]), ("B", st1["X2"])):
            need = ["full"] + [("prefix", j) for j in range(1, N)]
            if any((k, meth) not in rX for k in need):
                continue
            tX = {k: rmat.scale(RC(T), Xs[k]) for k in need}
            B, owner = {}, {}
            for key, (a, b) in (("G_NN", (N, N)), ("G_N1", (N, 1)), ("G_1N", (1, N))):
                B[key] = rmat.block(tX["full"], a, b, ny)
                owner[key] = "full"
            for j in range(1, N):
                B[("P", j)] = rmat.neg(rmat.block(tX[("prefix", j)], j, j, ny))
                B[("chi", j)] = rmat.block(tX[("prefix", j)], j, 1, ny)
                B[("psi", j)] = rmat.block(tX[("prefix", j)], 1, j, ny)
                owner[("P", j)] = owner[("chi", j)] = owner[("psi", j)] = ("prefix", j)
            if N >= 2:
                B["G1_11"] = rmat.block(tX[("prefix", 1)], 1, 1, ny)
                owner["G1_11"] = ("prefix", 1)
            E = Fraction(0)
            for c, keys in rpoly.monomials(m, N):
                up, pr = Fraction(1), Fraction(1)
                for kk in keys:
                    beta = ub(rmat.fro2(B[kk]))
                    up *= beta + T * rX[(owner[kk], meth)]
                    pr *= beta
                E += abs(c) * (up - pr)
            Fm = rpoly.F_ref(m, N, B)
            out.setdefault(("Q8", m), {})[meth] = _iv(lb(rmat.fro2(Fm)) - E, ub(rmat.fro2(Fm)) + E)
    return out


def _intersect(ivs):
    lo = max(i[0] for i in ivs)
    hi = min(i[1] for i in ivs)
    return (lo, hi) if lo <= hi else None


def classify_case(lat, N, s, Ys, prod, mutants, st1_cache=None):
    """Run Stage 1 + 2 at all three precisions; convergence; classification. Returns (outcome, rows, notes)."""
    notes, outcomes, rows = [], [], []
    st1s, st2s = {}, {}
    for p in DPS:
        try:
            if st1_cache is not None and p in st1_cache:
                st1s[p] = st1_cache[p]
            else:
                st1s[p] = stage1(lat, N, s, p)
                if st1_cache is not None:
                    st1_cache[p] = st1s[p]
        except Exception as e:  # noqa: BLE001
            return "NONCONVERGENT", rows, [f"mpmath exception at dps {p}: {type(e).__name__}: {e}"]
        if st1s[p]["info"]["invalid"]:
            outcomes.append("REFERENCE-INVALID")
            notes += st1s[p]["info"]["invalid"]
        st2s[p] = stage2(lat, N, s, st1s[p], Ys, prod, mutants)
        if st2s[p].get("nonconvergent"):
            outcomes.append("NONCONVERGENT")
            notes += st2s[p]["nonconvergent"]
    ell_B = 2 * st1s[500]["info"]["prec"]
    qkeys = sorted({k for p in DPS for k in st2s[p] if isinstance(k, tuple)}, key=str)
    for qk in qkeys:
        ivs = [iv for p in DPS for iv in st2s[p].get(qk, {}).values()]
        if not ivs:
            continue
        Istar = _intersect(ivs)
        if Istar is None:
            outcomes.append("NONCONVERGENT")
            notes.append(f"empty intersection for {qk}")
            continue
        Q, idx = qk
        if Q == "Q1":
            g = lat["gammas"][idx]
            encB = rlead.enclosure(g, ell_B)
            if _intersect([Istar, (encB["lo"], encB["hi"])]) is None:
                outcomes.append("NONCONVERGENT")
                notes.append(f"Q1 isqrt-enclosure disagreement for gamma={g}")
                continue
            lo, hi = prod["lo"][idx], prod["hi"][idx]
            cls = ("CONSISTENT" if lo <= Istar[0] and Istar[1] <= hi else
                   "CONTRADICTED" if Istar[1] < lo or Istar[0] > hi else "UNRESOLVED")
            claim = (lo, hi)
        elif Q == "Q7":
            b = prod["G_lo"]
            if b is None:
                continue
            cls = "CONSISTENT" if b < Istar[0] else "CONTRADICTED" if b > Istar[1] else "UNRESOLVED"
            claim = b
        else:
            b = {"Q2": prod["eps"], "Q3": prod["h"], "Q4": prod["omega"], "Q5": prod["a"], "Q6": prod["eta"],
                 "Q8": prod["U"]}[Q].get(idx)
            if b is None:
                continue
            cls = "CONSISTENT" if Istar[1] < b else "CONTRADICTED" if Istar[0] > b else "UNRESOLVED"
            claim = b
        outcomes.append(cls)
        rows.append(dict(Q=Q, item=str(idx), cls=cls, claim=claim, inf=Istar[0], sup=Istar[1]))
    return worst(outcomes), rows, notes


def selection_check(open_sel, closed_sel, edges, ell=rlead.ELL):
    """Re-evaluate the Rev3.5 eligibility/selection rules exactly; return list of violations (REFERENCE-INVALID)."""
    def eligible(g):
        _, xi = rlead.branch_and_xi(g)
        return not rlead.is_rational_square(xi)

    def sensitive(g):
        enc = rlead.enclosure(g, ell)
        return abs(enc["xi"] - enc["mu"] ** 2) >= Fraction(7, 8) * enc["e"] * (enc["hi"] + enc["mu"])

    bad = []
    o = [Fraction(k, 1024) for k in range(1, 1024) if eligible(Fraction(k, 1024)) and sensitive(Fraction(k, 1024))][:2]
    c = [1 + Fraction(k, 1024) for k in range(1, 1024) if eligible(1 + Fraction(k, 1024)) and sensitive(1 + Fraction(k, 1024))][:2]
    if o != list(open_sel):
        bad.append(f"open selection {o} != frozen {list(open_sel)}")
    if c != list(closed_sel):
        bad.append(f"closed selection {c} != frozen {list(closed_sel)}")
    for g in edges:
        if not eligible(g):
            bad.append(f"band edge {g} not eligible")
    return bad
