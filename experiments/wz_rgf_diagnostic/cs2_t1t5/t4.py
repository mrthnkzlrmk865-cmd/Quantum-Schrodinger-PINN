"""T4 - the b = 0 path (§2.4.2), with the production covering's live b = 0 instrumentation."""
from fractions import Fraction as F

from . import inputs as I
from .lattice import matrices
from .result import Results


def _M(rows):
    return tuple(tuple((F(x), F(0)) for x in r) for r in rows)


def run(prod, cov):
    R = Results("T4 b = 0 path")
    for f in (F(0), F(1), F(2), F(1, 4), F(1, 2**600)):
        u = prod._ub(f)
        if f == 0:
            R.add("T4 ub(0) = 0 exactly", u == 0 and isinstance(u, F), "CERT-VIOLATION", u)
        else:
            R.add(f"T4 ub({f if f.denominator < 2**20 else '2^-600'})^2 > f", u * u > f, "CERT-VIOLATION")
    # the instrumentation is live: evaluating 1/(2 a b) or r*b with b = 0 raises
    for fn, args in ((cov.half_radius, (F(1), F(0))), (cov.r_times_b, (F(1), F(0)))):
        try:
            fn(*args)
            live = False
        except cov.B0EvaluationError:
            live = True
        R.add(f"T4 instrumentation live: {fn.__name__}(b = 0) raises B0EvaluationError", live, "CERT-VIOLATION")
    # A(s) = diag(2, 3), A1 = 0: certified on [0, 1] from a single point with r = +inf, no instrumented exception
    try:
        rec = cov.cover(_M([[2, 0], [0, 3]]), _M([[0, 0], [0, 0]]), 0)
        exc = None
    except cov.B0EvaluationError as e:
        rec, exc = None, e
    R.add("T4 diag(2,3), A1 = 0: no instrumented exception", exc is None, "CERT-VIOLATION", exc)
    if rec:
        R.add("T4 diag(2,3): b = 0 flag set", rec["b_zero"] and rec["b"] == 0, "CERT-VIOLATION")
        R.add("T4 diag(2,3): [0, 1] covered", rec["covered"], "CERT-VIOLATION")
        R.add("T4 diag(2,3): exactly one certification point used", rec["n_points"] == 1 and len(rec["certified"]) == 1,
              "CERT-VIOLATION", rec["n_points"])
        R.add("T4 diag(2,3): r = +inf", rec["certified"] and rec["certified"][0]["r"] is cov.INF, "CERT-VIOLATION")
        r0 = cov.taylor_radius(rec)
        R.add("T4 diag(2,3): Taylor radius r_0 = +inf", r0 is cov.INF, "CERT-VIOLATION")
    # N = 2, N_y = 2 L-lattice with V_1(1) = 0, V_2(1) != 0: prefix 1 has b = 0, the full matrix has b != 0
    mats, V1, Ep = matrices(prod, I.T4_CONFIG, F(0))
    for key, m in mats.items():
        A0 = m["A"]
        try:
            rec = cov.cover(A0, m["A1"], m["eps"])
            exc = None
        except cov.B0EvaluationError as e:
            rec, exc = None, e
        tag = f"T4 L N=2 N_y=2 {key}"
        R.add(f"{tag}: no instrumented exception", exc is None, "CERT-VIOLATION", exc)
        if rec is None:
            continue
        want_b0 = key != "full"
        R.add(f"{tag}: b = 0 flag = {want_b0}", rec["b_zero"] == want_b0, "CERT-VIOLATION", rec["b"])
        if want_b0:
            R.add(f"{tag}: constant prefix certified from one point", rec["covered"] and rec["n_points"] == 1,
                  "CERT-VIOLATION", rec["n_points"])
        else:
            R.add(f"{tag}: full matrix (b > 0) covering reported: covered={rec['covered']}, points={rec['n_points']}", True)
    return R
