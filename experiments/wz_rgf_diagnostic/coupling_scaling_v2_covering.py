r"""§2.4 item 3 continuous-interval invertibility certificate on [0, 1], and the common Taylor radius r_0 (Fact 3, §5.2).

Production module (imports only the production certificate module). For an affine matrix A(s) = A0 + s*A1 with a
validated lead-enclosure bound eps_A:
  * every certification point s_k is evaluated by §2.4.1 (certify_with_retry on A_mid(s_k) = A0 + s_k*A1); uncertified
    points are discarded and never receive a radius;
  * a certified point gets the half radius r_k = 1/(2 a_k b) (b > 0), or r = +inf (b = 0, symbolic only);
  * coverage of [0, 1] by the certified points, exactly as in §2.4 item 3, with dyadic-midpoint bisection of every
    uncovered gap, a depth limit of 30 per gap and at most 2000 certification points per matrix;
  * for b = 0 a single certified point certifies all of [0, 1]; 1/(2 a b) and r*b are never evaluated.
b = 0 instrumentation: `half_radius` and `r_times_b` raise B0EvaluationError when called with b = 0 (T4).
"""
from __future__ import annotations

from fractions import Fraction

import coupling_scaling_v2_certify as C

MAX_DEPTH = 30
MAX_POINTS = 2000
GRID = sorted(C._GRID)                                      # the 81 preregistered §3 points
INITIAL_POINTS = sorted({Fraction(0)} | set(GRID) | {Fraction(1)})


class B0EvaluationError(ArithmeticError):
    """Raised if 1/(2 a b) or r*b is evaluated with b = 0 (forbidden by §2.4 item 3)."""


class _Infinity:
    """Symbolic +inf for radii; supports only the coverage comparisons, never arithmetic with b."""
    def __repr__(self):
        return "+inf"


INF = _Infinity()


def half_radius(a, b):
    if b == 0:
        raise B0EvaluationError("1/(2 a b) evaluated with b = 0")
    return Fraction(1) / (2 * a * b)


def r_times_b(r, b):
    if b == 0:
        raise B0EvaluationError("r*b evaluated with b = 0")
    return r * b


def _lo(s, r):          # s - r  (or -inf)
    return None if r is INF else s - r


def _hi(s, r):          # s + r  (or +inf)
    return None if r is INF else s + r


def _affine(A0, A1, s):
    return tuple(tuple((x[0] + s * y[0], x[1] + s * y[1]) for x, y in zip(ra, rb)) for ra, rb in zip(A0, A1))


def coverage_gaps(cert):
    """Uncovered gaps of [0, 1] given certified points [(s, r)] sorted by s (exact; +inf symbolic)."""
    if not cert:
        return [(Fraction(0), Fraction(1))]
    gaps = []
    s1, r1 = cert[0]
    if r1 is not INF and s1 - r1 > 0:
        gaps.append((Fraction(0), s1))
    for (sa, ra), (sb, rb) in zip(cert, cert[1:]):
        if ra is INF or rb is INF:
            continue
        if sb - sa > ra + rb:
            gaps.append((sa, sb))
    sK, rK = cert[-1]
    if rK is not INF and sK + rK < 1:
        gaps.append((sK, Fraction(1)))
    return gaps


def cover(A0, A1, eps_A, initial_points=None, max_depth=MAX_DEPTH, max_points=MAX_POINTS):
    """Continuous certificate of A(s) = A0 + s A1 on [0, 1]. Returns a record; 'covered' is the acceptance flag."""
    A0 = C._mat(A0, "A0")
    A1 = C._mat(A1, "A1")
    eps = C._q(eps_A, "eps_A")
    b = C.b_bound(A1)
    b_zero = (b == 0)
    pts = list(INITIAL_POINTS if initial_points is None else initial_points)
    evaluated = {}            # s -> record
    depth = {}

    def evaluate(s, d):
        c = C.certify_with_retry(_affine(A0, A1, s), eps)
        rec = dict(s=s, depth=d, status=c["status"], dps=c.get("dps"))
        if c["status"] == "certified":
            rec.update(a=c["a"], h=c["h"], omega=c["omega"])
            rec["r"] = INF if b_zero else half_radius(c["a"], b)
        evaluated[s] = rec
        depth[s] = d
        return rec

    if b_zero:
        # A(s) == A0 for all s: one certified point certifies [0, 1]; stop at the first certified point.
        for s in pts:
            if len(evaluated) >= max_points:
                break
            if evaluate(s, 0)["status"] == "certified":
                break
    else:
        for s in pts:
            evaluate(s, 0)
        while True:
            cert = sorted((s, r["r"]) for s, r in evaluated.items() if r["status"] == "certified")
            gaps = coverage_gaps(cert)
            if not gaps or len(evaluated) >= max_points:
                break
            new = 0
            for (l, u) in gaps:
                inside = sorted(x for x in evaluated if l <= x <= u)
                if inside[0] != l:
                    inside = [l] + inside
                if inside[-1] != u:
                    inside = inside + [u]
                for p, q in zip(inside, inside[1:]):
                    if len(evaluated) >= max_points:
                        break
                    m = (p + q) / 2
                    d = max(depth.get(p, 0), depth.get(q, 0)) + 1
                    if d > max_depth or m in evaluated:
                        continue
                    evaluate(m, d)
                    new += 1
            if new == 0:
                break
    cert = sorted((s, r["r"]) for s, r in evaluated.items() if r["status"] == "certified")
    covered = bool(cert) and not coverage_gaps(cert)
    certified = [evaluated[s] for s, _ in cert]
    cond = None
    min_r = None
    if certified:
        if b_zero:
            cond = max(r["a"] * r["h"] for r in certified)          # A(s) constant: cond <= a * h
            min_r = INF
        else:
            cond = max(2 * r["a"] * (r["h"] + r_times_b(r["r"], b)) for r in certified)
            min_r = min(r["r"] for r in certified)
    zero = evaluated.get(Fraction(0))
    return dict(covered=covered, b=b, b_zero=b_zero, n_points=len(evaluated), certified=certified,
                uncertified=[evaluated[s] for s in sorted(evaluated) if evaluated[s]["status"] != "certified"],
                cond_bound=cond, min_r=min_r, zero=zero,
                budget_exhausted=len(evaluated) >= max_points)


def taylor_radius(record):
    """r_0,A from the s = 0 certificate: 1/(2 a_0 b) (b > 0) or +inf (b = 0); None if s = 0 is not certified (R3.20)."""
    z = record["zero"]
    if z is None or z["status"] != "certified":
        return None
    if record["b_zero"]:
        return INF
    return half_radius(z["a"], record["b"])


def common_taylor_radius(records):
    """r_0 := min_A r_0,A = 1/(2 beta), beta := max_A a_A,0 b_A (r_0 = +inf iff beta = 0). None if any s = 0 is
    uncertified (the configuration then fails the filter)."""
    beta = Fraction(0)
    for rec in records:
        z = rec["zero"]
        if z is None or z["status"] != "certified":
            return None, None
        if not rec["b_zero"]:
            beta = max(beta, z["a"] * rec["b"])
    return (INF if beta == 0 else Fraction(1) / (2 * beta)), beta


def certify_configuration(matrices):
    """matrices: list of (name, A0, A1, eps_A). Accept iff [0, 1] is covered for every matrix and every s = 0 is
    certified. Records the validated cond bound, the minimum r_k and the common Taylor radius."""
    recs = {name: cover(A0, A1, eps) for name, A0, A1, eps in matrices}
    r0, beta = common_taylor_radius(list(recs.values()))
    accepted = all(r["covered"] for r in recs.values()) and r0 is not None
    return dict(accepted=accepted, records=recs, r0=r0, beta=beta)
