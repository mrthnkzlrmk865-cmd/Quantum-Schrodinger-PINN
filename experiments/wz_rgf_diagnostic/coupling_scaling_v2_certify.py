r"""Production certificate module for the W/Z coupling-scaling experiment v2 (PREDICTIONS_coupling_scaling_v2.md).

Scope of this file (Revision 3.7):
  * the production black-box API of §2.4.3.3 used by T6 (lead enclosure, A_mid assembly, eps_A, the point certificate of
    §2.4.1 steps 2-7, the §5.6 pointwise infeasibility certificate, point_infeasible, conclude), and
  * approx_inverse_Y, the §2.4.1 step 1 mpmath approximate inverse (the only place mpmath is used).
NOT in this file yet: the T1-T5 self-tests and the §2.4.1 interval covering / Taylor radius r_0 machinery. Those are
specified in §2.4.2 but are not implemented, so the manifest-freeze condition "T1-T5 pass" (§2.4.3.11) cannot be met.

Rules (§2.4.1, §2.4.3.3, §2.4.3.5):
  * Every number accepted or returned is int or Fraction (never bool); complex values are (re, im) pairs; matrices are
    tuples of row tuples of such pairs. Anything else raises TypeError. No floating-point arithmetic after import.
  * This module must not import the T6 reference (cs2_t6_reference) or the T6 harness (cs2_t6_harness).
  * U_m and G_lo are never clamped.
"""
from __future__ import annotations

from fractions import Fraction
from math import isqrt

ELL = 256
T_HOP = Fraction(32)              # t = 1/(2 dx^2), dx = 1/8
THRESHOLD = Fraction(1, 10**10)
SUPPORTED_MUTANTS = ("M1", "M2", "M3", "M5", "M6", "M7", "M8", "M15")
_MIN_N = {"M1": 3, "M2": 3, "M3": 3, "M5": 3, "M6": 2, "M7": 2, "M8": 2, "M15": 2}

# §3 design grid and S_test (all preregistered grid points with s <= 2^-19)
_GRID = frozenset({Fraction(1, 2**k) for k in range(41)} | {Fraction(3, 2**(k + 2)) for k in range(40)})
S_TEST = frozenset(s for s in _GRID if s <= Fraction(1, 2**19))


# ----------------------------------------------------------------------------------------------------- exact types
def _q(x, what="value"):
    if isinstance(x, bool) or not isinstance(x, (int, Fraction)):
        raise TypeError(f"{what}: exact int or Fraction required, got {type(x).__name__}")
    return Fraction(x)


def _qc(z, what="entry"):
    if not (isinstance(z, tuple) and len(z) == 2):
        raise TypeError(f"{what}: complex value must be an (re, im) pair of exact rationals")
    return (_q(z[0], what), _q(z[1], what))


def _mat(A, what="matrix"):
    if not isinstance(A, tuple) or not A or not all(isinstance(r, tuple) for r in A):
        raise TypeError(f"{what}: tuple of row tuples required")
    w = len(A[0])
    if any(len(r) != w for r in A):
        raise ValueError(f"{what}: ragged")
    return tuple(tuple(_qc(z, what) for z in r) for r in A)


def _cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def _csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def _abs2(a):
    return a[0] * a[0] + a[1] * a[1]


_Z = (Fraction(0), Fraction(0))
_O = (Fraction(1), Fraction(0))


def _eye(n):
    return tuple(tuple(_O if i == j else _Z for j in range(n)) for i in range(n))


def _zeros(n):
    return tuple(tuple(_Z for _ in range(n)) for _ in range(n))


def _mm(A, B):
    cols = list(zip(*B))
    out = []
    for r in A:
        row = []
        for c in cols:
            re = Fraction(0)
            im = Fraction(0)
            for x, y in zip(r, c):
                re += x[0] * y[0] - x[1] * y[1]
                im += x[0] * y[1] + x[1] * y[0]
            row.append((re, im))
        out.append(tuple(row))
    return tuple(out)


def _msub(A, B):
    return tuple(tuple(_csub(x, y) for x, y in zip(a, b)) for a, b in zip(A, B))


def _madd(A, B):
    return tuple(tuple(_cadd(x, y) for x, y in zip(a, b)) for a, b in zip(A, B))


def _mscale(c, A):
    return tuple(tuple(_cmul(c, x) for x in r) for r in A)


def _fro2(A):
    return sum((_abs2(x) for r in A for x in r), Fraction(0))


def _blk(A, a, b, ny):
    return tuple(r[(b - 1) * ny:b * ny] for r in A[(a - 1) * ny:a * ny])


# ------------------------------------------------------------------------------------- §2.4.1 rigorous square roots
def _ceil_q(q):
    return -((-q.numerator) // q.denominator)


def _floor_q(q):
    return q.numerator // q.denominator


def _ub(f, ell=ELL):
    """Upper square root: 0 for f = 0, else (isqrt(ceil(f*4^ell)) + 1) / 2^ell >= sqrt(f)."""
    if f < 0:
        raise ValueError("_ub of a negative number")
    if f == 0:
        return Fraction(0)
    return Fraction(isqrt(_ceil_q(Fraction(f) * 4**ell)) + 1, 2**ell)


def _lb(f, ell=ELL):
    """Lower square root: isqrt(floor(f*4^ell)) / 2^ell <= sqrt(f)."""
    if f < 0:
        raise ValueError("_lb of a negative number")
    return Fraction(isqrt(_floor_q(Fraction(f) * 4**ell)), 2**ell)


# ---------------------------------------------------------------------------------------------- lead enclosure
def lead_enclosure(gamma, ell=ELL):
    """§2.4.1 enclosure of sqrt(xi). `ell` is the enclosure resolution only (T6.8 P4 uses 4); _ub/_lb keep ELL."""
    g = _q(gamma, "gamma")
    ell = int(ell)
    if abs(g) < 1:
        branch, xi = "open", 1 - g * g
    elif g > 1:
        branch, xi = "closed", g * g - 1
    else:
        raise ValueError("gamma <= -1 and |gamma| = 1 are excluded")
    lo = Fraction(isqrt(_floor_q(xi * 4**ell)), 2**ell)
    hi = Fraction(isqrt(_ceil_q(xi * 4**ell)) + 1, 2**ell)
    mu = (lo + hi) / 2
    e = (hi - lo) / 2
    lam_mid = (g - mu, Fraction(0)) if branch == "closed" else (g, mu)
    return dict(lo=lo, hi=hi, mu=mu, e=e, branch=branch, lam_mid=lam_mid)


# ------------------------------------------------------------------------------------ A_mid assembly (Delta form)
def assemble_A_mid(Delta_blocks, gammas, kind, j=None, ell=ELL):
    """Physical units: A = t * (M/t), with M/t diagonal blocks 2 Delta_j - 2I (+ Lambda~ on slice 1, and on slice N for
    'full'; 2 Lambda~ when N = 1) and off-diagonal blocks I. Delta_blocks: the N slices' Delta_j (N_y x N_y)."""
    D = [_mat(B, "Delta block") for B in Delta_blocks]
    lam = [lead_enclosure(g, ell)["lam_mid"] for g in gammas]
    ny = len(lam)
    if kind == "full":
        n = len(D)
    elif kind == "prefix":
        n = int(j)
        if not 1 <= n <= len(D) - 1:
            raise ValueError("prefix index out of range")
    else:
        raise ValueError("kind must be 'full' or 'prefix'")
    size = n * ny
    A = [[_Z] * size for _ in range(size)]
    t = (T_HOP, Fraction(0))
    for b in range(n):
        for p in range(ny):
            for r in range(ny):
                v = _cmul((Fraction(2), Fraction(0)), D[b][p][r])
                if p == r:
                    v = _csub(v, (Fraction(2), Fraction(0)))
                A[b * ny + p][b * ny + r] = _cmul(t, v)
        if b + 1 < n:
            for p in range(ny):
                A[b * ny + p][(b + 1) * ny + p] = t
                A[(b + 1) * ny + p][b * ny + p] = t
    for p in range(ny):
        A[p][p] = _cadd(A[p][p], _cmul(t, lam[p]))                            # left lead on slice 1
        if kind == "full":
            k = (n - 1) * ny + p
            A[k][k] = _cadd(A[k][k], _cmul(t, lam[p]))                        # right lead on slice N
    return tuple(tuple(r) for r in A)


def eps_A(gammas, N, kind, ell=ELL):
    """Validated bound of ||A_e||_F (§2.4.1 'Lead perturbation, per matrix type')."""
    root = _ub  # eps_A rounding
    s_te2 = sum(((T_HOP * lead_enclosure(g, ell)["e"]) ** 2 for g in gammas), Fraction(0))
    s_2te2 = sum(((2 * T_HOP * lead_enclosure(g, ell)["e"]) ** 2 for g in gammas), Fraction(0))
    if kind == "prefix":
        return root(s_te2)                        # prefix: left lead only
    if N == 1:
        return root(s_2te2)                       # full N = 1: one block receives both leads
    return root(2 * s_te2)                        # full N >= 2: slices 1 and N


# --------------------------------------------------------------------------------- §2.4.1 steps 2-7 point certificate
def certify_point(A_mid, Y, eps_A_value):
    A = _mat(A_mid, "A_mid")
    Yq = _mat(Y, "Y")
    eps = _q(eps_A_value, "eps_A")
    Xi = _msub(_eye(len(A)), _mm(A, Yq))
    f_Xi, f_Y, f_mid = _fro2(Xi), _fro2(Yq), _fro2(A)
    omega = _ub(f_Xi) + eps * _ub(f_Y)
    h = _ub(f_mid) + eps
    certified = omega < 1
    out = dict(status="certified" if certified else "uncertified", f_Xi=f_Xi, f_Y=f_Y, f_mid=f_mid, omega=omega, h=h,
               a=None)
    if certified:
        a = _ub(f_Y) / (1 - omega)
        out["a"] = a
    return out


# --------------------------------------------------------------------------------------------- §5.6 polynomials
def _monomials(mutant, N):
    if mutant in ("M1", "M3"):
        longest = ["G_NN", "G1_11"] + [("P", j) for j in range(2, N)]  # M1/M3 longest
        single = (-1, ["G_N1"])  # M1/M3 single
        return [(-1, longest), single]
    if mutant == "M2":
        longest = [("P", j) for j in range(N - 1, 1, -1)] + ["G1_11", "G_NN"]  # M2 longest
        single = (-1, ["G_1N"])  # M2 single
        return [(-1, longest), single]
    if mutant == "M5":
        terms = []
        for j in range(2, N):  # M5 terms
            first = [("chi", j - 1), ("P", j), ("psi", j - 1)]  # M5 longest
            terms.append((-1, first))
            terms.append((+1, [("psi", j - 1), ("P", j), ("chi", j - 1)]))
        return terms
    if mutant in ("M6", "M15"):
        longest = [("chi", N - 1), "G_NN"]  # M6/M15 longest
        single = (-1, ["G_N1"])  # M6/M15 single
        return [(-1, longest), single]
    if mutant == "M7":
        longest = ["G_NN", ("psi", N - 1)]  # M7 longest
        single = (-1, ["G_1N"])  # M7 single
        return [(-1, longest), single]
    if mutant == "M8":
        longest = [("chi", N - 1), "G_NN", ("psi", N - 1)]  # M8 longest
        return [(+1, longest), (-1, [("psi", N - 1), "G_NN", ("chi", N - 1)])]
    raise ValueError(f"unsupported mutant {mutant}")


def F_poly(mutant, N, blocks):
    """F_m = (§4.2 oracle) - (designated reference block), evaluated on the given N_y x N_y blocks."""
    B = {k: _mat(v, f"block {k}") for k, v in blocks.items()}
    ny = len(B["G_NN"])
    F = _zeros(ny)
    for c, keys in _monomials(mutant, N):
        X = B[keys[0]]
        for k in keys[1:]:
            X = _mm(X, B[k])
        F = _madd(F, _mscale((Fraction(c), Fraction(0)), X))
    return F


def _required_prefixes(mutant, N):
    return list(range(1, N)) if mutant in ("M1", "M2", "M3", "M5") else [N - 1]


def _inconclusive(reason, **extra):
    return dict(status="INCONCLUSIVE", reason=reason, point="not_infeasible", **extra)


def pointwise_certificate(mutant, N, N_y, refs):
    """§5.6 pointwise certificate. refs: "full" and ("prefix", j) -> (A_mid, Y, eps_A)."""
    if mutant not in SUPPORTED_MUTANTS:
        return _inconclusive("unsupported_mutant")
    if N < _MIN_N[mutant] or (mutant == "M8" and N == 2):
        return dict(status="not applicable", reason=None, point="not_infeasible")
    keys = ["full"] + [("prefix", j) for j in _required_prefixes(mutant, N)]
    for key in keys:
        if key not in refs:
            return _inconclusive("missing_reference")
    certs = {k: certify_point(*refs[k]) for k in keys}
    omegas = {k: c["omega"] for k, c in certs.items()}
    if any(c["status"] != "certified" for c in certs.values()):
        return _inconclusive("uncertified_reference", omega=omegas)
    t = (T_HOP, Fraction(0))
    eta = {}
    for k, c in certs.items():
        omega = c["omega"]
        eta[k] = T_HOP * _ub(c["f_Y"]) * omega / (1 - omega)
    blocks, owner = {}, {}
    tYf = _mscale(t, _mat(refs["full"][1]))
    blocks["G_NN"], blocks["G_N1"], blocks["G_1N"] = _blk(tYf, N, N, N_y), _blk(tYf, N, 1, N_y), _blk(tYf, 1, N, N_y)
    owner.update({"G_NN": "full", "G_N1": "full", "G_1N": "full"})
    for j in _required_prefixes(mutant, N):
        tYp = _mscale(t, _mat(refs[("prefix", j)][1]))
        blocks[("P", j)] = _mscale((Fraction(-1), Fraction(0)), _blk(tYp, j, j, N_y))
        blocks[("chi", j)] = _blk(tYp, j, 1, N_y)
        blocks[("psi", j)] = _blk(tYp, 1, j, N_y)
        owner[("P", j)] = owner[("chi", j)] = owner[("psi", j)] = ("prefix", j)
        if j == 1:
            blocks["G1_11"] = _blk(tYp, 1, 1, N_y)
            owner["G1_11"] = ("prefix", 1)
    beta = {}
    E_m = Fraction(0)
    for c, mono in _monomials(mutant, N):
        prod_up = Fraction(1)
        prod = Fraction(1)
        for k in mono:
            if k not in beta:
                beta[k] = _ub(_fro2(blocks[k]))
            prod_up *= beta[k] + eta[owner[k]]
            prod *= beta[k]
        E_m += abs(c) * (prod_up - prod)
    F = F_poly(mutant, N, blocks)
    U_m = _ub(_fro2(F)) + _ub(Fraction(N_y)) * E_m
    h_W = certs["full"]["h"]
    maxabs2 = max(_abs2(z) for r in _mat(refs["full"][1]) for z in r)
    eta_full = eta["full"]
    G_lo = max(T_HOP / h_W, T_HOP * _lb(maxabs2) - eta_full)
    base = dict(eta=eta, beta=beta, E_m=E_m, U_m=U_m, G_lo=G_lo, omega=omegas)
    if U_m < 0:
        return _inconclusive("invalid_upper_bound", **base)
    if G_lo < 0:
        return _inconclusive("invalid_lower_bound", **base)
    point = "infeasible" if point_infeasible(U_m, G_lo, THRESHOLD) else "not_infeasible"
    return dict(status="computed", reason=None, point=point, **base)


def point_infeasible(U, G_lo, threshold):
    """Generic contract (§2.4.3.3): types, threshold, U >= 0, G_lo >= 0, then U < threshold * G_lo."""
    U = _q(U, "U")
    G_lo = _q(G_lo, "G_lo")
    threshold = _q(threshold, "threshold")
    if threshold != Fraction(1, 10**10):
        raise ValueError("threshold must equal Fraction(1, 10**10)")
    if U < 0:  # generic upper-bound check
        raise ValueError("INVALID_UPPER_BOUND")
    if G_lo < 0:  # generic lower-bound check
        raise ValueError("INVALID_LOWER_BOUND")
    return U < threshold * G_lo


def conclude(per_point):
    """Every-point rule over exactly S_test (§5.6 conclusion rule)."""
    if frozenset(per_point) != S_TEST:
        return dict(outcome="INCONCLUSIVE", reason="evaluated_set_mismatch")
    rs = [per_point[s] for s in sorted(per_point, reverse=True)]
    if all(r["status"] == "not applicable" for r in rs):
        return dict(outcome="not applicable", reason=None)
    for r in rs:
        if r["status"] == "INCONCLUSIVE":
            return dict(outcome="INCONCLUSIVE", reason=r["reason"])
        if r["status"] == "not applicable":
            return dict(outcome="not applicable", reason=None)
    if all(r["point"] == "infeasible" for r in rs):
        return dict(outcome="infeasible", reason=None)
    return dict(outcome="INCONCLUSIVE", reason="bound_not_below_threshold")


# -------------------------------------------------------------------- §2.4.1 step 1: mpmath approximate inverse Y
def approx_inverse_Y(A_mid, dps=60):
    """Y = exact dyadic value of the mpmath inverse of A_mid at `dps` (60; 120 on retry). The only mpmath use."""
    import mpmath

    A = _mat(A_mid, "A_mid")
    ctx = mpmath.MPContext()
    ctx.dps = int(dps)
    M = ctx.matrix(len(A), len(A))
    for i, r in enumerate(A):
        for k, (re, im) in enumerate(r):
            M[i, k] = ctx.mpc(ctx.mpf(re.numerator) / re.denominator, ctx.mpf(im.numerator) / im.denominator)
    X = ctx.inverse(M)

    def exact(x):
        # mpf._mpf_ = (sign, man, exp, bc); value = (-1)^sign * man * 2^exp. (.man_exp drops the sign.)
        sign, man, e, _bc = ctx.mpf(x)._mpf_
        if man == 0 and e != 0:
            raise ValueError("non-finite mpmath value")
        v = Fraction(man) * (Fraction(2) ** e) if e >= 0 else Fraction(man, 2 ** (-e))
        return -v if sign else v

    return tuple(tuple((exact(ctx.re(X[i, k])), exact(ctx.im(X[i, k]))) for k in range(len(A))) for i in range(len(A)))
