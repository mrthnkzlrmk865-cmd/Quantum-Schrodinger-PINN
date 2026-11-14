r"""§4.7 item 4 validated references (PREDICTIONS_coupling_scaling_v2.md, Revision 3.7), plus the §4.7 cross-checks.

Every reference inverse is the exact dyadic Y of §2.4.1 step 1 (`approx_inverse_Y`, no further rounding; R3.4), certified
by §2.4.1 steps 2-6 (`certify_point`) with the lattice's own enclosure bound ε_A. Acceptance is an exact rational
comparison:

    ω < 1  and  ω/(1 − ω) ≤ 10^-55          (dps 60; on failure dps 120; otherwise INVALID)

Validated error bounds (exact rationals; Ĝ = t·A⁻¹ units):
  * err_Ĝ := t · ū(‖Y‖_F²) · ω/(1 − ω)  ≥ ‖Ĝ_true − t·Y‖₂, which bounds every entry and every block;
  * P_j = −Ĝ^(j)_jj:  err_P,j := err_Ĝ^(j);
  * W_j = I + SC_j, SC_1 = A_11/t (slice 1 with the left lead), SC_j = A_jj/t − Ĝ^(j−1)_{j−1,j−1}:
    err_W,1 := ε_A(prefix 1)/t,  err_W,j := err_Ĝ^(j−1) (A_jj, j ≥ 2, is exact);
  * S_ref,nn′ = 2i·ū(μ_n μ_n′)·Ĝ^Y_nn′ − δ_nn′ on open channels (exact rationals; ū(μμ′) ∈ [l̄(lo lo′), ū(hi hi′)]), with
    |ΔS_ref,nn′| ≤ 2·ū(hi_n hi_n′)·err_Ĝ + 2·|Ĝ^Y_nn′|^up·(ū(hi_n hi_n′) − l̄(lo_n lo_n′)), |z|^up := ū(re² + im²).
A bound that is not computed is reported as None and the quantity is then *not reliably measured* (§6).

Cross-checks (empirical; INVALID on disagreement after escalation):
  * (b) dps-100 agreement: ‖Ĝ_100 − t·Y‖_F / ‖Ĝ_100‖_F ≤ 1e-55 (relative Frobenius, the inherited convention of
    coupled_blocks.reference);
  * (c) mp.lu_solve vs mp.inverse: ‖t·X_lu − t·X_inv‖_F / ‖t·X_inv‖_F ≤ 1e-55 at dps 60 (tolerance not stated in §4.7;
    the inherited coupled_blocks threshold — implementer choice C-6, flagged for R0).
This module imports only the production certificate module and mpmath; it contains no recursion code.
"""
from __future__ import annotations

from fractions import Fraction

import coupling_scaling_v2_certify as C

ACCEPT = Fraction(1, 10**55)
DPS_LEVELS = (60, 120)
CHECK_TOL = "1e-55"


class ReferenceInvalid(RuntimeError):
    """A §4.7 validated-reference control failed after escalation (INVALID)."""


def validated_inverse(A_mid, eps_A, levels=DPS_LEVELS):
    """Validated reference of A⁻¹ (physical units). Returns dict(status, dps, Y, omega, rel, normY_ub, err_inv, attempts).
    status 'validated' or 'INVALID' (never raises for a numerical failure; the caller decides)."""
    attempts = []
    for dps in levels:
        try:
            Y = C.approx_inverse_Y(A_mid, dps)
        except (ZeroDivisionError, ValueError, ArithmeticError) as exc:
            attempts.append(dict(dps=dps, status="inverse-failed", reason=repr(exc)))
            continue
        c = C.certify_point(A_mid, Y, eps_A)
        om = c["omega"]
        if om < 1:
            rel = om / (1 - om)
            ok = rel <= ACCEPT
        else:
            rel, ok = None, False
        attempts.append(dict(dps=dps, status="accepted" if ok else "rejected", omega=om, rel=rel))
        if ok:
            nY = C._ub(c["f_Y"])
            return dict(status="validated", dps=dps, Y=Y, omega=om, rel=rel, normY_ub=nY, err_inv=nY * rel,
                        attempts=attempts)
    return dict(status="INVALID", dps=None, Y=None, omega=None, rel=None, normY_ub=None, err_inv=None,
                attempts=attempts)


def err_ghat(ref):
    """t · ‖Y‖₂ · ω/(1−ω) (Ĝ units), or None if the reference is not validated."""
    if ref["status"] != "validated":
        return None
    return C.T_HOP * ref["err_inv"]


def ghat_from_Y(Y):
    """t·Y as an exact matrix (Ĝ units)."""
    t = C.T_HOP
    return tuple(tuple((t * re, t * im) for re, im in row) for row in Y)


def block(X, a, b, ny):
    """1-indexed ny × ny block (a, b) of an exact matrix."""
    return tuple(tuple(X[(a - 1) * ny + p][(b - 1) * ny + r] for r in range(ny)) for p in range(ny))


def w_p_bounds(prefix_refs, eps_prefix1, N):
    """Validated error bounds of the W_j and P_j references, j = 1..N−1.
    prefix_refs: {j: validated_inverse record of M^(j)}; eps_prefix1: ε_A of M^(1) (physical units)."""
    out = {}
    for j in range(1, N):
        eP = err_ghat(prefix_refs[j]) if j in prefix_refs else None
        if j == 1:
            eW = Fraction(eps_prefix1) / C.T_HOP
        else:
            eW = err_ghat(prefix_refs[j - 1]) if (j - 1) in prefix_refs else None
        out[j] = dict(err_P=eP, err_W=eW)
    return out


def w_reference(A_full_mid, prefix_refs, N, ny):
    """W_j^ref = I + SC_j (exact, from A_mid and the validated prefix Y), j = 1..N−1, as exact ny × ny matrices."""
    t = C.T_HOP
    out = {}
    for j in range(1, N):
        Ajj = block(A_full_mid, j, j, ny)
        SC = [[(z[0] / t, z[1] / t) for z in row] for row in Ajj]
        if j >= 2:
            Gp = ghat_from_Y(prefix_refs[j - 1]["Y"])
            B = block(Gp, j - 1, j - 1, ny)
            SC = [[(SC[p][r][0] - B[p][r][0], SC[p][r][1] - B[p][r][1]) for r in range(ny)] for p in range(ny)]
        out[j] = tuple(tuple(((SC[p][r][0] + (1 if p == r else 0)), SC[p][r][1]) for r in range(ny)) for p in range(ny))
    return out


def _abs_up(z):
    return C._ub(z[0] * z[0] + z[1] * z[1])


def s_reference(full_ref, gammas, N, ny):
    """Validated S reference on open channels: returns dict(open=[n], S={(a,b): (re, im)}, bound={(a,b): Fraction|None}).
    Index layout: a, b ∈ 0..2·n_open−1, first n_open = left lead (block 1), next n_open = right lead (block N)."""
    enc = [C.lead_enclosure(g) for g in gammas]
    opn = [n for n in range(ny) if enc[n]["branch"] == "open"]
    no = len(opn)
    validated = full_ref["status"] == "validated"
    G = ghat_from_Y(full_ref["Y"]) if full_ref.get("Y") is not None else None
    eG = err_ghat(full_ref)
    S, bnd = {}, {}
    blocks = {(0, 0): (1, 1), (0, 1): (1, N), (1, 0): (N, 1), (1, 1): (N, N)}
    for (sa, sb), (ba, bb) in blocks.items():
        for i, n in enumerate(opn):
            for k, m in enumerate(opn):
                a, b = sa * no + i, sb * no + k
                if G is None:
                    S[(a, b)], bnd[(a, b)] = None, None
                    continue
                g = G[(ba - 1) * ny + n][(bb - 1) * ny + m]
                pref = C._ub(enc[n]["mu"] * enc[m]["mu"])
                # 2i·pref·g − δ
                re = -2 * pref * g[1] - (1 if a == b else 0)
                im = 2 * pref * g[0]
                S[(a, b)] = (re, im)
                if validated:
                    up = C._ub(enc[n]["hi"] * enc[m]["hi"])
                    lo = C._lb(enc[n]["lo"] * enc[m]["lo"])
                    bnd[(a, b)] = 2 * up * eG + 2 * _abs_up(g) * (up - lo)
                else:
                    bnd[(a, b)] = None
    return dict(open=opn, S=S, bound=bnd)


# ------------------------------------------------------------------------------------------------ mp cross-checks
def _mp_matrix(ctx, A):
    n = len(A)
    M = ctx.matrix(n, len(A[0]))
    for i, r in enumerate(A):
        for k, (re, im) in enumerate(r):
            M[i, k] = ctx.mpc(ctx.mpf(re.numerator) / re.denominator, ctx.mpf(im.numerator) / im.denominator)
    return M


def _fro(ctx, X):
    return ctx.sqrt(ctx.fsum(abs(X[i, k]) ** 2 for i in range(X.rows) for k in range(X.cols)))


def check_dps100(A_mid, Y):
    """(b): relative Frobenius distance between the dps-100 inverse and Y. Returns (value as mpf string, ok)."""
    import mpmath
    ctx = mpmath.MPContext()
    ctx.dps = 100
    Xh = ctx.inverse(_mp_matrix(ctx, A_mid))
    Ym = _mp_matrix(ctx, Y)
    v = _fro(ctx, Xh - Ym) / _fro(ctx, Xh)
    return ctx.nstr(v, 8), bool(v <= ctx.mpf(CHECK_TOL))


def check_lu_vs_inverse(A_mid, dps=60):
    """(c): mp.lu_solve (column by column) vs mp.inverse at the given dps."""
    import mpmath
    ctx = mpmath.MPContext()
    ctx.dps = int(dps)
    M = _mp_matrix(ctx, A_mid)
    n = M.rows
    Xi = ctx.inverse(M)
    Xl = ctx.matrix(n, n)
    for k in range(n):
        col = ctx.lu_solve(M, ctx.matrix([[1 if i == k else 0] for i in range(n)]))
        for i in range(n):
            Xl[i, k] = col[i]
    v = _fro(ctx, Xl - Xi) / _fro(ctx, Xi)
    return ctx.nstr(v, 8), bool(v <= ctx.mpf(CHECK_TOL))


def cross_checks(A_mid, Y, do_b, do_c, lu_dps=60):
    out = {}
    if do_b:
        v, ok = check_dps100(A_mid, Y)
        out["b"] = dict(value=v, ok=ok)
    if do_c:
        v, ok = check_lu_vs_inverse(A_mid, lu_dps)
        out["c"] = dict(value=v, ok=ok, dps=lu_dps)
    out["ok"] = all(x["ok"] for k, x in out.items() if k in ("b", "c"))
    return out


def reference(A_mid, eps_A, do_b=False, do_c=False):
    """Validated reference plus requested cross-checks (§4.7 4(a)-(c)). On a (b)/(c) failure at dps 60 the reference is
    rebuilt at dps 120 (escalation) and the checks repeated; INVALID if it persists. record['valid'] is the flag."""
    r = validated_inverse(A_mid, eps_A)
    r["eps"] = Fraction(eps_A)
    if r["status"] != "validated":
        r["checks"], r["valid"], r["escalated"] = {}, False, False
        return r
    chk = cross_checks(A_mid, r["Y"], do_b, do_c, lu_dps=60)
    r["escalated"] = False
    if not chk["ok"] and r["dps"] == 60:
        first = chk
        r2 = validated_inverse(A_mid, eps_A, levels=(120,))
        r2["eps"] = Fraction(eps_A)
        r2["attempts"] = r["attempts"] + r2["attempts"]
        r = r2
        r["escalated"] = True
        if r["status"] != "validated":
            r["checks"], r["valid"] = dict(first=first, ok=False), False
            return r
        chk = cross_checks(A_mid, r["Y"], do_b, do_c, lu_dps=120)
        chk["first"] = first
    r["checks"] = chk
    r["valid"] = chk["ok"]
    return r
