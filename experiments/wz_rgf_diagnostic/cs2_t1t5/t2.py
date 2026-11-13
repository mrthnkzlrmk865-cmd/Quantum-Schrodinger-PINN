"""T2 - no certification across singularities (§2.4.2), on the production covering of §2.4 item 3."""
from fractions import Fraction as F

from .result import Results


def _M(rows):
    return tuple(tuple((F(x[0]), F(x[1])) if isinstance(x, tuple) else (F(x), F(0)) for x in r) for r in rows)


def _diag(*v):
    n = len(v)
    return _M([[v[i] if i == j else 0 for j in range(n)] for i in range(n)])


def _e24(c):
    return _M([[(c if (i, j) == (0, 0) else 0) for j in range(24)] for i in range(24)])


def cases():
    I2 = _diag(1, 1)
    I24 = _M([[1 if i == j else 0 for j in range(24)] for i in range(24)])
    e = F(1, 2**10)
    return [  # (id, A0, A1, real singular points, complex singular points (re, im^2), coverage must fail)
        ("a: diag(1 - 2s, 1), singular at 1/2 (dyadic, bisection midpoint)", I2, _diag(-2, 0), [F(1, 2)], [], True),
        ("b: diag(1 - 3s, 1), singular at 1/3 (non-dyadic)", I2, _diag(-3, 0), [F(1, 3)], [], True),
        ("c: diag(s, 1), singular at 0 (endpoint)", _diag(0, 1), _diag(1, 0), [F(0)], [], True),
        ("d: diag(1 - s, 1), singular at 1 (endpoint)", I2, _diag(-1, 0), [F(1)], [], True),
        ("e: 24x24 I - 3s e1e1^T, singular at 1/3", I24, _e24(-3), [F(1, 3)], [], True),
        ("f: [[1-2s, 2^-10], [-2^-10, 1-2s]], singular at 1/2 +- i 2^-11", _M([[1, e], [-e, 1]]), _diag(-2, -2), [],
         [(F(1, 2), F(1, 2**22))], False),
        ("g: constant singular diag(0, 1), b = 0", _diag(0, 1), _diag(0, 0), "everywhere", [], True),
    ]


def run(prod, cov):
    R = Results("T2 no certification across singularities")
    for cid, A0, A1, real_sing, cplx_sing, must_fail in cases():
        rec = cov.cover(A0, A1, 0)
        tag = f"T2({cid})"
        if must_fail:
            R.add(f"{tag}: coverage of [0, 1] fails", not rec["covered"], "CERT-VIOLATION")
            acc = cov.certify_configuration([("A", A0, A1, 0)])
            R.add(f"{tag}: configuration rejected by the filter", not acc["accepted"], "CERT-VIOLATION")
        if real_sing == "everywhere":
            R.add(f"{tag}: no point certified for a constant singular matrix", not rec["certified"], "CERT-VIOLATION",
                  len(rec["certified"]))
            continue
        bad = []
        for c in rec["certified"]:
            r = c["r"]
            for sg in real_sing:          # certified interval [s - r, s + r] must exclude sg:  r < |s - sg|
                if r is cov.INF or not r < abs(c["s"] - sg):
                    bad.append((c["s"], sg))
            for (re_, im2) in cplx_sing:  # certified disk must exclude both s = re +- i sqrt(im2):  r^2 < (s - re)^2 + im2
                if r is cov.INF or not r * r < (c["s"] - re_) ** 2 + im2:
                    bad.append((c["s"], (re_, im2)))
        R.add(f"{tag}: no certified interval/disk contains a singular point ({len(rec['certified'])} certified)", not bad,
              "CERT-VIOLATION", bad[:3])
        if cplx_sing:
            r0 = cov.taylor_radius(rec)
            re_, im2 = cplx_sing[0]
            R.add(f"{tag}: s = 0 certified (Taylor disk exists)", r0 is not None, "CERT-VIOLATION")
            if r0 is not None:
                R.add(f"{tag}: Taylor radius r_0 < |s_sing| (r_0^2 < 1/4 + 2^-22, exact)",
                      r0 is not cov.INF and r0 * r0 < re_ * re_ + im2, "CERT-VIOLATION", r0)
            R.add(f"{tag}: real-interval coverage reported (may succeed): covered={rec['covered']}", True)
    return R
