r"""§5 analysis procedures of the coupling-scaling experiment v2 (PREDICTIONS_coupling_scaling_v2.md, Revision 3.7).

Pure functions (no solver, no I/O), testable on fixtures:
  * SeriesBackend: truncated matrix power series (Cauchy products), so the §4.2 oracle runs unchanged on Taylor series;
  * neumann_series: Ĝ(s) = t·A(s)⁻¹ = t·Σ_k (−s)^k (A0⁻¹A1)^k A0⁻¹ up to K_max (§5.2);
  * zero_test / identify_order (§5.2 zero test, thresholds unchanged);
  * s_star, b_tail (validated majorant, exact rationals), b_prec (§5.2 (ii));
  * ols / interval / regression_validity / window_adequacy (§5.3);
  * local_slopes, crossover, regime_models (§5.4); noise_trend (§5.5); bootstrap_median (§5.3).
Implementer choices (flagged for R0, R0_REVIEW_PACKET.md Part I): C-8 regime-model fit (relative linear least squares
of D*/s^q on [1, s, …, s^mreg] over all grid points with D* > 0); C-9 bootstrap spawn key (k,) per reported group in
the fixed report order.
"""
from __future__ import annotations

import math
from fractions import Fraction

K_MAX = 12
Z_NONZERO_MIN = "1e-60"
Z_REL_AGREE = "1e-20"
Z_ZERO_100 = "1e-80"
Z_ZERO_200 = "1e-170"
S_STAR_DEFAULT = Fraction(1, 2**60)
WINDOW_EDGE = Fraction(1, 2**12)
SNR_REG = 1e3
SNR_DET = 10.0
BOOT_N = 10**4


# --------------------------------------------------------------------------------------------- series backend
class SeriesBackend:
    """Matrices are lists [C_0, …, C_K] of ctx matrices; products are Cauchy products truncated at K."""

    def __init__(self, ctx, K=K_MAX):
        self.ctx, self.K = ctx, K

    def mul(self, A, B):
        out = []
        for k in range(self.K + 1):
            acc = None
            for i in range(k + 1):
                term = A[i] * B[k - i]
                acc = term if acc is None else acc + term
            out.append(acc)
        return out

    def add(self, A, B):
        return [a + b for a, b in zip(A, B)]

    def sub(self, A, B):
        return [a - b for a, b in zip(A, B)]

    def neg(self, A):
        return [-a for a in A]

    def scal(self, c, A):
        return [self.ctx.mpf(c) * a for a in A]

    def inv(self, A):
        raise NotImplementedError("series inverse is not used (M10 has no series, §5.2)")

    def block(self, X, a, b, ny):
        return [x[(a - 1) * ny:a * ny, (b - 1) * ny:b * ny] for x in X]

    def diag(self, vals):
        raise NotImplementedError("M10 has no series (§5.2)")


def neumann_series(ctx, A0, A1, t, K=K_MAX):
    """Coefficients of Ĝ(s) = t·(A0 + s A1)⁻¹ (ctx matrices): C_k = t·(−1)^k (A0⁻¹A1)^k A0⁻¹."""
    A0i = ctx.inverse(A0)
    T = A0i * A1
    out, cur = [], A0i
    for k in range(K + 1):
        out.append((ctx.mpf(t) * (-1) ** k) * cur)
        cur = T * cur
    return out


# ----------------------------------------------------------------------------------------------- zero test, order
def zero_test(z100, z200):
    """§5.2: 'nonzero', 'zero' or 'unresolved' (mpf inputs, any context)."""
    import mpmath
    if z200 >= mpmath.mpf(Z_NONZERO_MIN) and abs(z100 - z200) <= mpmath.mpf(Z_REL_AGREE) * z200:
        return "nonzero"
    if z100 <= mpmath.mpf(Z_ZERO_100) and z200 <= mpmath.mpf(Z_ZERO_200):
        return "zero"
    return "unresolved"


def identify_order(tests):
    """tests: zero-test results for C_0..C_K. Returns ('q', q) | ('all_zero', None) | ('unresolved', j)."""
    for k, r in enumerate(tests):
        if r == "nonzero":
            return ("q", k)
        if r == "unresolved":
            return ("unresolved", k)
    return ("all_zero", None)


# ------------------------------------------------------------------------------------------------ §5.2 (ii)
def s_star(r0):
    """2⁻⁶⁰ if 2⁻⁶⁰ ≤ r_0, otherwise r_0·2⁻²⁰ (r0: Fraction or the string '+inf')."""
    if r0 == "+inf":
        return S_STAR_DEFAULT
    r0 = Fraction(r0)
    return S_STAR_DEFAULT if S_STAR_DEFAULT <= r0 else r0 / 2**20


def b_tail(T_m, deg_m, alpha, beta, s, K=K_MAX):
    """Validated tail majorant, exact rationals. Halves s until ζx < 1. Returns (B_tail, s_used, x, ζ).
    B_tail = T_m·α^deg·C(K+deg, deg−1)·x^(K+1)/(1 − ζx), x = |s|β, ζ = (K+deg+1)/(K+2)."""
    T_m, alpha, beta, s = Fraction(T_m), Fraction(alpha), Fraction(beta), Fraction(s)
    d = int(deg_m)
    zeta = Fraction(K + d + 1, K + 2)
    while True:
        x = abs(s) * beta
        if zeta * x < 1:
            break
        s = s / 2
    if beta == 0:
        return Fraction(0), s, Fraction(0), zeta
    if x > Fraction(1, 2):
        raise ValueError("x = |s*| beta > 1/2 outside the Taylor disk")
    B = T_m * alpha**d * math.comb(K + d, d - 1) * x ** (K + 1) / (1 - zeta * x)
    return B, s, x, zeta


def b_prec(dps, mag):
    """Non-validated precision allowance 10^−(dps−10)·Mag_m(s*)."""
    import mpmath
    return mpmath.mpf(10) ** (-(dps - 10)) * mag


def alpha_from(a0_list):
    """α = max(1, max_A t·a_A,0)."""
    from coupling_scaling_v2_certify import T_HOP
    return max([Fraction(1)] + [T_HOP * Fraction(a) for a in a0_list])


# ------------------------------------------------------------------------------------------------------- §5.3
def ols(x, y):
    """Unweighted OLS slope, intercept, weights ϑ_i (slope = Σ ϑ_i y_i) and the descriptive SE."""
    n = len(x)
    if n < 2:
        return None
    xm = sum(x) / n
    sxx = sum((xi - xm) ** 2 for xi in x)
    if sxx == 0:
        return None
    th = [(xi - xm) / sxx for xi in x]
    p = sum(t * yi for t, yi in zip(th, y))
    ym = sum(y) / n
    a = ym - p * xm
    res = [yi - (a + p * xi) for xi, yi in zip(x, y)]
    se = math.sqrt(sum(r * r for r in res) / (n - 2) / sxx) if n > 2 else float("nan")
    return dict(slope=p, intercept=a, weights=th, se=se, n=n)


def window(points):
    """Ω_asym: grid s ≤ 2⁻¹² with SNR*_m ≥ 1e3. points: [(s, D, Dstar, snr)]. Returns the filtered list."""
    return [p for p in points if Fraction(p[0]) <= WINDOW_EDGE and p[3] is not None and p[3] >= SNR_REG]


def regression(points, q):
    """§5.3 on Ω_asym points [(s, D, Dstar, snr)]. Returns the record with validity, adequacy-ready slopes and Hyp2."""
    W = sorted(window(points), key=lambda p: Fraction(p[0]))
    rec = dict(n_points=len(W), valid=False, span=None, p_hat=None, p_hat_star=None, interval=None, se=None, hyp2=None)
    if len(W) < 2 or any(p[1] <= 0 or p[2] <= 0 for p in W):
        rec["reason"] = "fewer than 2 positive window points"
        return rec
    x = [math.log(float(p[0])) for p in W]
    span = max(x) - min(x)
    rec["span"] = span
    o = ols(x, [math.log(p[1]) for p in W])
    os_ = ols(x, [math.log(p[2]) for p in W])
    rec.update(p_hat=o["slope"], p_hat_star=os_["slope"], se=o["se"],
               interval=sum(abs(t) * abs(math.log(p[1] / p[2])) for t, p in zip(o["weights"], W)))
    rec["valid"] = len(W) >= 8 and span >= math.log(100)
    if q is not None:
        rec["hyp2"] = abs(o["slope"] - q) <= 0.02 and abs(o["slope"] - os_["slope"]) <= 0.005
    return rec


# ------------------------------------------------------------------------------------------------------- §5.4
def local_slopes(svals, D):
    """Centred φ(s_k) = ln(D_{k+1}/D_{k−1})/ln(s_{k+1}/s_{k−1}) on the sorted grid; one-sided φ(1) from s = 1 and 1/2.
    svals increasing. Returns {s: φ} (only where both values are > 0)."""
    out = {}
    for k in range(1, len(svals) - 1):
        a, b = D[k - 1], D[k + 1]
        if a > 0 and b > 0:
            out[svals[k]] = math.log(b / a) / math.log(float(svals[k + 1]) / float(svals[k - 1]))
    if svals[-1] == 1 and Fraction(1, 2) in svals:
        i = svals.index(Fraction(1, 2))
        if D[-1] > 0 and D[i] > 0:
            out[Fraction(1)] = math.log(D[-1] / D[i]) / math.log(2.0)
    return out


def window_adequacy(reg, phi_exact, q):
    """|p̂* − q| ≤ 0.02 and every exact centred φ in Ω_asym within 0.02 of q."""
    if reg["p_hat_star"] is None:
        return False
    return abs(reg["p_hat_star"] - q) <= 0.02 and all(abs(v - q) <= 0.02 for v in phi_exact)


def crossover(svals, Dstar, Cq, q, phi_exact):
    """s_× = the largest grid s with |dev(s′)| ≤ 0.1 for all grid s′ ≤ s, dev = D*/(‖C_q‖ s^q) − 1; secondary: the largest
    s with |φ_exact − q| ≤ 0.05 for every computed φ at and below it."""
    sx = None
    for s, d in zip(svals, Dstar):
        if Cq <= 0:
            break
        if abs(d / (Cq * float(s) ** q) - 1) <= 0.1:
            sx = s
        else:
            break
    sec = None
    for s in svals:
        if s in phi_exact:
            if abs(phi_exact[s] - q) <= 0.05:
                sec = s
            else:
                break
    return sx, sec


def regime_models(svals, Dstar, q, max_m=2):
    """Descriptive (choice C-8): fit D*/s^q ≈ c0 (1 + Σ b_i s^i) by relative linear least squares; select the smallest
    mreg with max |log residual| ≤ 0.01; otherwise 'strongly nonlinear in range'."""
    import numpy as np
    pts = [(float(s), d) for s, d in zip(svals, Dstar) if d > 0]
    if len(pts) < max_m + 2:
        return dict(selected=None, note="too few points")
    s = np.array([p[0] for p in pts])
    y = np.array([p[1] for p in pts]) / s**q
    for mreg in range(max_m + 1):
        X = np.vstack([s**i for i in range(mreg + 1)]).T / y[:, None]
        c, *_ = np.linalg.lstsq(X, np.ones_like(y), rcond=None)
        model = X @ c * y
        if np.any(model <= 0):
            continue
        r = np.max(np.abs(np.log(y) - np.log(model)))
        if r <= 0.01:
            return dict(selected=mreg, max_log_residual=float(r), a=float(np.log(abs(c[0]))),
                        b=[float(ci / c[0]) for ci in c[1:]])
    return dict(selected=None, note="strongly nonlinear in range")


# ------------------------------------------------------------------------------------------------------- §5.5
def noise_trend(svals, eta, eps_ref):
    """OLS slope of ln max(η_X, ε‖X_ref‖_F) on ln s over all grid s ≤ 2⁻¹² (57 points); pass iff in [−0.25, 0.25]."""
    pts = [(math.log(float(s)), math.log(max(e, f))) for s, e, f in zip(svals, eta, eps_ref) if Fraction(s) <= WINDOW_EDGE]
    o = ols([p[0] for p in pts], [p[1] for p in pts])
    if o is None:
        return dict(slope=None, n=len(pts), ok=False)
    return dict(slope=o["slope"], n=len(pts), ok=-0.25 <= o["slope"] <= 0.25)


def bootstrap_median(values, key, root=202611110003, n=BOOT_N):
    """Seeded bootstrap of the median (choice C-9: SeedSequence(entropy=root, spawn_key=(key,)))."""
    import numpy as np
    if not values:
        return None
    g = np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy=root, spawn_key=(int(key),))))
    v = np.asarray(values, dtype=float)
    meds = np.median(v[g.integers(0, len(v), size=(n, len(v)))], axis=1)
    return dict(median=float(np.median(v)), lo=float(np.percentile(meds, 2.5)), hi=float(np.percentile(meds, 97.5)), n=n)
