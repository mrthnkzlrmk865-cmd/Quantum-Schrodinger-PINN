"""T1 - validated inverse bound (§2.4.2)."""
from fractions import Fraction as F
import mpmath

from cs2_t6_reference import rmat, rpsd, rinv
from . import inputs as I, indep
from .lattice import matrices
from .result import Results


def _pairs(rows):
    return tuple(tuple((F(x[0]), F(x[1])) if isinstance(x, tuple) else (F(x), F(0)) for x in r) for r in rows)


def _exact_types(c):
    return all(isinstance(c[k], F) and not isinstance(c[k], bool) for k in ("omega", "h") if c.get(k) is not None) and \
        (c.get("a") is None or isinstance(c["a"], F))


def run(prod):
    R = Results("T1 validated inverse bound")
    Ad = _pairs([[2, 0, 0], [0, -3, 0], [0, 0, F(1, 2)]])
    Ap = _pairs([[0, 2, 0, 0], [0, 0, 2, 0], [0, 0, 0, 2], [2, 0, 0, 0]])         # 2 * cyclic permutation
    Ac = _pairs([[(1, 1), F(1, 2)], [0, (1, -1)]])
    inv = {"a": _pairs([[F(1, 2), 0, 0], [0, F(-1, 3), 0], [0, 0, 2]]),
           "b": _pairs([[0, 0, 0, F(1, 2)], [F(1, 2), 0, 0, 0], [0, F(1, 2), 0, 0], [0, 0, F(1, 2), 0]]),
           "c": _pairs([[(F(1, 2), F(-1, 2)), F(-1, 4)], [0, (F(1, 2), F(1, 2))]])}
    closed = {"a": (F(2), F(3)), "b": (F(1, 2), F(2))}                               # ||A^-1||_2, ||A||_2
    mats = {"a": Ad, "b": Ap, "c": Ac}
    for k, A in mats.items():
        Ainv = rmat.mat(inv[k])
        R.add(f"T1({k}) reference inverse is exact (A.X = X.A = I)", rinv.inverse(rmat.mat(A)) == Ainv, "KAT-FAIL")
        for yname, Y in (("Y_prod (dps 60/120)", None), ("Y = A^-1 exact (Xi = 0)", tuple(tuple((z.re, z.im) for z in r) for r in Ainv))):
            c = prod.certify_with_retry(A, 0) if Y is None else prod.certify_point(A, Y, 0)
            tag = f"T1({k}) {yname}"
            if not R.add(f"{tag} certified", c["status"] == "certified", "CERT-VIOLATION"):
                continue
            R.add(f"{tag} exact rational omega, a, h", _exact_types(c), "INTERFACE-FAIL")
            a, h = c["a"], c["h"]
            b = prod.b_bound(A)                                                        # A1 := A (test-file choice)
            if k in closed:
                ninv, nA = closed[k]
                R.add(f"{tag} a >= ||A^-1||_2 = {ninv} (closed form)", a >= ninv, "CERT-VIOLATION")
                R.add(f"{tag} h >= ||A||_2 = {nA} (closed form)", h >= nA, "CERT-VIOLATION")
                R.add(f"{tag} b >= ||A1||_2 = {nA} (A1 := A, closed form)", b >= nA, "CERT-VIOLATION")
            else:
                R.add(f"{tag} a >= sigma_max(A^-1) (exact 2x2 formula)", indep.sigma_max_2x2_le(inv[k], a), "CERT-VIOLATION")
                R.add(f"{tag} h >= sigma_max(A) (exact 2x2 formula)", indep.sigma_max_2x2_le(A, h), "CERT-VIOLATION")
                R.add(f"{tag} b >= sigma_max(A1) (A1 := A, exact 2x2 formula)", indep.sigma_max_2x2_le(A, b), "CERT-VIOLATION")
            R.add(f"{tag} PSD cross-check a^2 I - A^-H A^-1 >= 0", rpsd.norm2_le(Ainv, a), "CERT-VIOLATION")
            if Y is not None:
                R.add(f"{tag} edge: Xi = 0 so omega = eps_A * ub(f_Y) = 0", c["omega"] == 0 and c["f_Xi"] == 0, "CERT-VIOLATION")
    # (d) poor Y
    Yd = tuple(tuple((z.re / 2, z.im / 2) for z in r) for r in rmat.mat(inv["a"]))
    c = prod.certify_point(Ad, Yd, 0)
    R.add("T1(d) Y = A^-1/2 for (a): ||Xi||_F = sqrt(3)/2 < 1 -> certified", c["status"] == "certified", "CERT-VIOLATION",
          c.get("omega"))
    if c["status"] == "certified":
        R.add("T1(d) loose bound still valid: a >= ||A^-1||_2 = 2", c["a"] >= 2, "CERT-VIOLATION")
    Yb = tuple(tuple((z.re / 2, z.im / 2) for z in r) for r in rmat.mat(inv["b"]))
    c = prod.certify_point(Ap, Yb, 0)
    R.add("T1(d) Y = A^-1/2 for (b): ||Xi||_F = 1 -> uncertified, no bound", c["status"] == "uncertified" and c["a"] is None,
          "CERT-VIOLATION", c.get("omega"))
    c = prod.certify_point(Ad, tuple(tuple((F(0), F(0)) for _ in range(3)) for _ in range(3)), 0)
    R.add("T1(e) Y = 0 -> uncertified, no bound", c["status"] == "uncertified" and c["a"] is None, "CERT-VIOLATION")
    R.add("T1 edge: eps_A = 0 (no leads) used in (a)-(e)", True)
    # (f) lead-enclosure paths at s = 0: mpmath dps-300 SVD reference, margin 10^-290 * sigma_ref
    ctx = mpmath.MPContext()
    ctx.dps = 300
    margin = F(1) + F(1, 10**290)

    def to_mp(A):
        M = ctx.matrix(len(A), len(A[0]))
        for i, r in enumerate(A):
            for j, z in enumerate(r):
                M[i, j] = ctx.mpc(ctx.mpf(z[0].numerator) / z[0].denominator, ctx.mpf(z[1].numerator) / z[1].denominator)
        return M

    def mpq(x):  # exact rational value of a (finite, signed) mpf
        sgn, man, e, _ = ctx.mpf(x)._mpf_
        v = F(man) * F(2) ** e if e >= 0 else F(man, 2 ** (-e))
        return -v if sgn else v

    for cfg in I.T1F_CONFIGS:
        name, ny, N = cfg[0], cfg[1], cfg[2]
        mats_, V1, Ep = matrices(prod, cfg, F(0))
        E = cfg[4]
        for key, m in mats_.items():
            tag = f"T1(f) {name} {key}"
            c = prod.certify_with_retry(m["A"], m["eps"])
            if not R.add(f"{tag} certified", c["status"] == "certified", "CERT-VIOLATION"):
                continue
            R.add(f"{tag} exact rational omega, a, h", _exact_types(c), "INTERFACE-FAIL")
            # true lattice: replace the midpoint lead roots by the true roots (mp.sqrt at dps 300)
            At = to_mp(m["A"])
            for p in range(ny):
                g = 1 - (E - Ep[p]) * I.DX2
                enc = indep.enclosure(g)
                xi = ctx.mpf(enc["xi"].numerator) / enc["xi"].denominator
                lam_t = (ctx.mpf(g.numerator) / g.denominator - ctx.sqrt(xi)) if enc["branch"] == "closed" else \
                    ctx.mpc(ctx.mpf(g.numerator) / g.denominator, ctx.sqrt(xi))
                lam_m = ctx.mpc(ctx.mpf(enc["lam"][0].numerator) / enc["lam"][0].denominator,
                                ctx.mpf(enc["lam"][1].numerator) / enc["lam"][1].denominator)
                d = 32 * (lam_t - lam_m)
                At[p, p] += d
                if m["full"]:
                    k = (m["n"] - 1) * ny + p
                    At[k, k] += d
            S = ctx.svd_c(At, compute_uv=False)
            smax, smin = max(S), min(S)
            S1 = ctx.svd_c(to_mp(m["A1"]), compute_uv=False)
            b = prod.b_bound(m["A1"])
            R.add(f"{tag} a >= ||A_true^-1||_2 (1 + 10^-290) [dps-300 SVD]", c["a"] >= mpq(1 / smin) * margin, "CERT-VIOLATION")
            R.add(f"{tag} h >= ||A_true||_2 (1 + 10^-290) [dps-300 SVD]", c["h"] >= mpq(smax) * margin, "CERT-VIOLATION")
            R.add(f"{tag} b >= ||A1||_2 (1 + 10^-290) [dps-300 SVD]", b >= mpq(max(S1)) * margin, "CERT-VIOLATION")
    return R
