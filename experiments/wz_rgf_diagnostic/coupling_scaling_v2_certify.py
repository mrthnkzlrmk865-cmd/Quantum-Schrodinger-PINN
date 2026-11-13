r"""Production certificate module for the W/Z coupling-scaling experiment v2 (PREDICTIONS_coupling_scaling_v2.md).

Scope of this file (Revision 3.7):
  * the production black-box API of §2.4.3.3 used by T6 (lead enclosure, A_mid assembly, eps_A, the point certificate of
    §2.4.1 steps 2-7, eta_A and G_lo, the §5.6 pointwise infeasibility certificate, point_infeasible, conclude);
  * approx_inverse_Y / certify_with_retry, the §2.4.1 step 1 and step 8 mpmath approximate inverse with retry (the only
    place mpmath is used);
  * the §2.4.1 binary64 import boundary and the exact L / R_W / R_orig assembly, A1 and b (T5).
The §2.4 item 3 covering and the common Taylor radius live in coupling_scaling_v2_covering.py. The T1-T5 self-tests are in
the cs2_t1t5 package and the T6 harness in cs2_t6_harness (independent of this module).

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


def _eta_from_cert(cert):
    """eta_A := t * ub(f_Y) * omega / (1 - omega) (§5.6), from a certified point certificate."""
    omega = cert["omega"]
    return T_HOP * _ub(cert["f_Y"]) * omega / (1 - omega)


def _g_lo_from(cert_full, Y_full):
    """G_lo := max(t/h_W, t*lb(max_ab |Y_full,ab|^2) - eta_full) (§5.6); never clamped."""
    h_W = cert_full["h"]
    maxabs2 = max(_abs2(z) for r in Y_full for z in r)
    eta_full = _eta_from_cert(cert_full)
    G_lo = max(T_HOP / h_W, T_HOP * _lb(maxabs2) - eta_full)
    return G_lo


def eta_A_bound(A_mid, Y, eps_A_value):
    """Per-matrix reference-error bound eta_A (§5.6) for any certified matrix, independent of mutants."""
    c = certify_point(A_mid, Y, eps_A_value)
    if c["status"] != "certified":
        raise ValueError("eta_A is defined only for a certified matrix (omega < 1)")
    return _eta_from_cert(c)


def G_lo_bound(A_mid_full, Y_full, eps_A_full):
    """Full-matrix lower bound G_lo (§5.6), independent of mutants."""
    c = certify_point(A_mid_full, Y_full, eps_A_full)
    if c["status"] != "certified":
        raise ValueError("G_lo is defined only for a certified full matrix (omega < 1)")
    return _g_lo_from(c, _mat(Y_full, "Y_full"))


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
    eta = {k: _eta_from_cert(c) for k, c in certs.items()}
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
    G_lo = _g_lo_from(certs["full"], _mat(refs["full"][1]))
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

    return tuple(tuple((_mpf_to_fraction(ctx, ctx.re(X[i, k])), _mpf_to_fraction(ctx, ctx.im(X[i, k])))
                       for k in range(len(A))) for i in range(len(A)))


def _mpf_to_fraction(ctx, x):
    """Exact dyadic value of an mpf: mpf._mpf_ = (sign, man, exp, bc), value = (-1)^sign * man * 2^exp.
    NaN and +/-inf (man == 0 with a special exponent) raise ValueError; nothing is approximated or mapped to zero.
    (mpf.man_exp is NOT used: in mpmath 1.3.0 it drops the sign.)"""
    sign, man, e, _bc = ctx.mpf(x)._mpf_
    if man == 0 and e != 0:
        raise ValueError("non-finite mpmath value (NaN or infinity)")
    v = Fraction(man) * (Fraction(2) ** e) if e >= 0 else Fraction(man, 2 ** (-e))
    return -v if sign else v


# ================================================================================================================
# §2.4.1 "Exact representation of the inputs; the binary64 import boundary" and the exact assembly paths (T5)
# ================================================================================================================
DX2 = Fraction(1, 64)


def import_binary64(x):
    """The ONLY entry point for binary64 values: exact dyadic rational of a finite binary64 (Python float or a
    numpy.floating). NaN and +/-inf raise ValueError; bool/int/str/complex raise TypeError."""
    if isinstance(x, bool) or not isinstance(x, float):
        try:
            import numpy as _np  # only for the numpy.floating type check; no arithmetic
            if isinstance(x, _np.floating) and not isinstance(x, _np.longdouble):
                x = float(x)
            else:
                raise TypeError
        except (ImportError, TypeError):
            raise TypeError(f"import_binary64: binary64 value required, got {type(x).__name__}") from None
    if x != x or x in (float("inf"), float("-inf")):
        raise ValueError("import_binary64: NaN and infinities are rejected")
    return Fraction(x)  # exact for every finite binary64


def import_complex128(z):
    """Exact (re, im) pair of a complex128 value; both parts go through import_binary64."""
    if isinstance(z, bool) or not hasattr(z, "real") or not hasattr(z, "imag") or isinstance(z, (int, Fraction)):
        raise TypeError(f"import_complex128: complex binary64 value required, got {type(z).__name__}")
    return (import_binary64(float(z.real)), import_binary64(float(z.imag)))


def import_matrix(Mx):
    """Exact QMat of a 2-D array of binary64 (real or complex) values."""
    out = []
    for row in Mx:
        r = []
        for z in row:
            if isinstance(z, complex) or (hasattr(z, "imag") and not isinstance(z, float) and z.imag != 0) or \
                    type(z).__name__.startswith("complex"):
                r.append(import_complex128(z))
            else:
                r.append((import_binary64(z), Fraction(0)))
        out.append(tuple(r))
    return tuple(out)


def _check_exact_matrix_inputs(blocks, what):
    return [_mat(B, what) for B in blocks]


def assemble_L(E, E_perp, V1_blocks, s, kind, j=None, ell=ELL):
    """L lattice, physical units: M = E*I - H(s) - Sigma~_L (+) Sigma~_R with H_jj = 2t I + diag(E_perp) + s V_j(1),
    H_{j,j+1} = -t I and Sigma~ = -t lambda~ from the enclosure midpoints of gamma_n = 1 - (E - E_perp,n) dx^2.
    Exact inputs only (imported rationals); a float raises TypeError."""
    E = _q(E, "E")
    Ep = [_q(e, "E_perp") for e in E_perp]
    s = _q(s, "s")
    V = _check_exact_matrix_inputs(V1_blocks, "V(1) block")
    ny, N = len(Ep), len(V)
    n = N if kind == "full" else int(j)
    if kind not in ("full", "prefix") or (kind == "prefix" and not 1 <= n <= N - 1):
        raise ValueError("bad kind / prefix index")
    lam = [lead_enclosure(1 - (E - e) * DX2, ell)["lam_mid"] for e in Ep]
    t = (T_HOP, Fraction(0))
    size = n * ny
    A = [[_Z] * size for _ in range(size)]
    for b in range(n):
        for p in range(ny):
            for r in range(ny):
                h = _cmul((s, Fraction(0)), V[b][p][r])
                if p == r:
                    h = _cadd(h, (2 * T_HOP + Ep[p], Fraction(0)))
                A[b * ny + p][b * ny + r] = _csub((E, Fraction(0)) if p == r else _Z, h)
        if b + 1 < n:
            for p in range(ny):
                A[b * ny + p][(b + 1) * ny + p] = t
                A[(b + 1) * ny + p][b * ny + p] = t
    for p in range(ny):
        A[p][p] = _cadd(A[p][p], _cmul(t, lam[p]))                     # -Sigma_L = t lambda~
        if kind == "full":
            k = (n - 1) * ny + p
            A[k][k] = _cadd(A[k][k], _cmul(t, lam[p]))                 # -Sigma_R (same block when N = 1)
    return tuple(tuple(r) for r in A)


def assemble_A1(V1_blocks, kind, j=None):
    """A1 = dA/ds = -blockdiag(V_j(1)) over the slices of the matrix (Sigma does not depend on s)."""
    V = _check_exact_matrix_inputs(V1_blocks, "V(1) block")
    ny, N = len(V[0]), len(V)
    n = N if kind == "full" else int(j)
    size = n * ny
    A = [[_Z] * size for _ in range(size)]
    for b in range(n):
        for p in range(ny):
            for r in range(ny):
                A[b * ny + p][b * ny + r] = _csub(_Z, V[b][p][r])
    return tuple(tuple(r) for r in A)


def eps_A_L(E, E_perp, N, kind, ell=ELL):
    """eps_A of the L lattice (lead parameters gamma_n = 1 - (E - E_perp,n) dx^2)."""
    E = _q(E, "E")
    return eps_A([1 - (E - _q(e, "E_perp")) * DX2 for e in E_perp], N, kind, ell)


def assemble_R_W(Delta_fl_blocks, dlead_fl, kind, j=None, ell=ELL):
    """R_W lattice from the imported Delta_fl blocks and the imported lead deltas delta_lead,fl (gamma = 1 - delta)."""
    gam = [1 - _q(d, "delta_lead") for d in dlead_fl]
    return assemble_A_mid(Delta_fl_blocks, gam, kind, j, ell)


def assemble_R_orig(A_orig_blocks, sigma_code, kind, j=None):
    """R_orig from the imported production arrays E*eye - H[j] and sigma_code (no enclosure: eps_A = 0).
    M = blocktridiag(A_orig_j; +t I) - Sigma_L on slice 1 - Sigma_R on slice N."""
    B = _check_exact_matrix_inputs(A_orig_blocks, "A_orig block")
    sig = [_qc(z, "sigma_code") for z in sigma_code]
    ny, N = len(sig), len(B)
    n = N if kind == "full" else int(j)
    if kind not in ("full", "prefix") or (kind == "prefix" and not 1 <= n <= N - 1):
        raise ValueError("bad kind / prefix index")
    t = (T_HOP, Fraction(0))
    size = n * ny
    A = [[_Z] * size for _ in range(size)]
    for b in range(n):
        for p in range(ny):
            for r in range(ny):
                A[b * ny + p][b * ny + r] = B[b][p][r]
        if b + 1 < n:
            for p in range(ny):
                A[b * ny + p][(b + 1) * ny + p] = t
                A[(b + 1) * ny + p][b * ny + p] = t
    for p in range(ny):
        A[p][p] = _csub(A[p][p], sig[p])
        if kind == "full":
            k = (n - 1) * ny + p
            A[k][k] = _csub(A[k][k], sig[p])
    return tuple(tuple(r) for r in A)


def b_bound(A1):
    """b := ub(||A1||_F^2) >= ||A1||_2 (§2.4.1 step 7); b = 0 iff A1 = 0."""
    return _ub(_fro2(_mat(A1, "A1")))


def certify_with_retry(A_mid, eps_A_value):
    """§2.4.1 steps 1-8: Y at dps 60; if omega >= 1 (or the mpmath inverse fails), retry once at dps 120; if still not
    certified the point is 'uncertified' (no bound is returned)."""
    last = None
    for dps in (60, 120):
        try:
            Y = approx_inverse_Y(A_mid, dps)
        except (ZeroDivisionError, ValueError, ArithmeticError):
            last = dict(status="uncertified", reason=f"mpmath inverse failed at dps {dps}", a=None, dps=dps)
            continue
        c = certify_point(A_mid, Y, eps_A_value)
        c["dps"] = dps
        if c["status"] == "certified":
            return c
        last = c
    return last
