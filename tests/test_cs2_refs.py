"""§4.7 item 4 validated-reference machinery (coupling_scaling_v2_refs): acceptance, rejection, escalation, validated
W/P/S bounds against exact inverses, the (b)/(c) cross-checks, and module independence. Fixtures only (no manifest)."""
import ast
import pathlib
import sys
from fractions import Fraction as F

import pytest

HERE = pathlib.Path(__file__).resolve().parent.parent / "experiments" / "wz_rgf_diagnostic"
sys.path.insert(0, str(HERE))
mpmath = pytest.importorskip("mpmath")
import coupling_scaling_v2_certify as C  # noqa: E402
import coupling_scaling_v2_refs as R  # noqa: E402

T = C.T_HOP


def cm(rows):
    return tuple(tuple((F(a), F(b)) for a, b in r) for r in rows)


def exact_inv(A):
    """Exact Gauss-Jordan inverse over Q(i) (pairs)."""
    n = len(A)
    M = [[complex_q(A[i][k]) for k in range(n)] + [((F(1) if i == k else F(0)), F(0)) for k in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != (0, 0))
        M[c], M[p] = M[p], M[c]
        inv = cinv(M[c][c])
        M[c] = [cmul(inv, x) for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != (0, 0):
                f = M[r][c]
                M[r] = [csub(x, cmul(f, y)) for x, y in zip(M[r], M[c])]
    return tuple(tuple(r[n:]) for r in M)


def complex_q(z):
    return (F(z[0]), F(z[1]))


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def csub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def cinv(a):
    d = a[0] ** 2 + a[1] ** 2
    return (a[0] / d, -a[1] / d)


def well(n, shift=5):
    return cm([[(shift * 32 if i == k else (i - k) * 3, (i + 2 * k) % 3) for k in range(n)] for i in range(n)])


def test_accepts_well_conditioned_with_exact_comparison_and_bounds_hold():
    A = well(4)
    r = R.reference(A, F(0), do_b=True, do_c=True)
    assert r["status"] == "validated" and r["valid"] and r["dps"] == 60
    assert isinstance(r["rel"], F) and r["rel"] <= R.ACCEPT
    X = exact_inv(A)
    e = R.err_ghat(r)
    G = R.ghat_from_Y(r["Y"])
    for i in range(4):
        for k in range(4):
            d = csub((T * X[i][k][0], T * X[i][k][1]), G[i][k])
            assert d[0] ** 2 + d[1] ** 2 <= e * e
    assert r["checks"]["b"]["ok"] and r["checks"]["c"]["ok"]


def test_rejects_when_enclosure_too_large_then_invalid_after_escalation():
    A = well(3)
    r = R.reference(A, F(1, 10**40))          # eps * ||Y|| >> 1e-55
    assert r["status"] == "INVALID" and not r["valid"]
    assert [a["dps"] for a in r["attempts"]] == [60, 120]
    assert all(a["status"] == "rejected" for a in r["attempts"])


def test_escalates_to_dps120_for_ill_conditioned_matrix():
    # cond ~ 1e40: dps-60 residual ~ 1e-20 > 1e-55; dps 120 accepts
    big = F(10**20)
    A = cm([[(big, 0), (big - 1, 0)], [(big + 1, 0), (big, 0)]])     # det = 1
    r = R.validated_inverse(A, F(0))
    assert r["status"] == "validated" and r["dps"] == 120
    assert [a["status"] for a in r["attempts"]] == ["rejected", "accepted"]


def test_singular_is_invalid_never_silent():
    A = cm([[(1, 0), (2, 0)], [(2, 0), (4, 0)]])
    r = R.reference(A, F(0), do_b=True, do_c=True)
    assert r["status"] == "INVALID" and not r["valid"] and R.err_ghat(r) is None


def test_w_p_bounds_against_exact_prefix_inverses():
    ny, N = 2, 3
    gam = [F(1, 2), F(3, 2)]                            # one open (|gamma| < 1), one closed (gamma > 1) channel
    D = [cm([[(F(1, 3), 0), (F(1, 8), 0)], [(F(1, 8), 0), (F(-1, 5), 0)]]),
         cm([[(F(1, 7), 0), (F(-1, 16), 0)], [(F(-1, 16), 0), (F(2, 9), 0)]]),
         cm([[(F(1, 4), 0), (F(1, 32), 0)], [(F(1, 32), 0), (F(1, 6), 0)]])]
    full = C.assemble_A_mid(D, gam, "full")
    pre = {j: R.reference(C.assemble_A_mid(D, gam, "prefix", j), C.eps_A(gam, N, "prefix")) for j in (1, 2)}
    assert all(p["valid"] for p in pre.values())
    eps1 = C.eps_A(gam, N, "prefix")
    bnd = R.w_p_bounds(pre, eps1, N)
    Wref = R.w_reference(full, pre, N, ny)
    # exact W_j with the midpoint lattice: W_j = I + (Ĝ^(j)_jj)^-1 (algebraically identical to I + SC_j)
    for j in (1, 2):
        Gx = exact_inv(C.assemble_A_mid(D, gam, "prefix", j))
        Gjj = tuple(tuple((T * Gx[(j - 1) * ny + p][(j - 1) * ny + r][0], T * Gx[(j - 1) * ny + p][(j - 1) * ny + r][1])
                          for r in range(ny)) for p in range(ny))
        Wx = exact_inv(Gjj)
        for p in range(ny):
            for r in range(ny):
                w = (Wx[p][r][0] + (1 if p == r else 0), Wx[p][r][1])
                d = csub(w, Wref[j][p][r])
                assert d[0] ** 2 + d[1] ** 2 <= (bnd[j]["err_W"] + F(1, 10**70)) ** 2
                Pd = csub((-Gjj[p][r][0], -Gjj[p][r][1]), (-T * pre[j]["Y"][(j - 1) * ny + p][(j - 1) * ny + r][0],
                                                          -T * pre[j]["Y"][(j - 1) * ny + p][(j - 1) * ny + r][1]))
                assert Pd[0] ** 2 + Pd[1] ** 2 <= bnd[j]["err_P"] ** 2
    assert bnd[1]["err_W"] == eps1 / T and bnd[2]["err_W"] == R.err_ghat(pre[1])


def test_s_reference_bound_contains_true_s():
    ny, N = 2, 2
    gam = [F(1, 2), F(-1, 4)]                           # two open channels
    D = [cm([[(F(1, 3), 0), (F(1, 8), 0)], [(F(1, 8), 0), (F(-1, 5), 0)]]),
         cm([[(F(1, 7), 0), (F(-1, 16), 0)], [(F(-1, 16), 0), (F(2, 9), 0)]])]
    A = C.assemble_A_mid(D, gam, "full")
    ref = R.reference(A, C.eps_A(gam, N, "full"))
    out = R.s_reference(ref, gam, N, ny)
    assert out["open"] == [0, 1] and all(v is not None for v in out["bound"].values())
    X = exact_inv(A)                                      # midpoint lattice inverse (|A_e| is inside the bound)
    with mpmath.workdps(120):
        s = [mpmath.sqrt(1 - (mpmath.mpf(g.numerator) / g.denominator) ** 2) for g in gam]
        blocks = {(0, 0): (1, 1), (0, 1): (1, N), (1, 0): (N, 1), (1, 1): (N, N)}
        for (sa, sb), (ba, bb) in blocks.items():
            for i in range(2):
                for k in range(2):
                    a, b = sa * 2 + i, sb * 2 + k
                    g = X[(ba - 1) * ny + i][(bb - 1) * ny + k]
                    G = mpmath.mpc(mpmath.mpf(g[0].numerator) / g[0].denominator, mpmath.mpf(g[1].numerator) / g[1].denominator) * 32
                    St = 2j * mpmath.sqrt(s[i] * s[k]) * G - (1 if a == b else 0)
                    Sr = out["S"][(a, b)]
                    Srm = mpmath.mpc(mpmath.mpf(Sr[0].numerator) / Sr[0].denominator, mpmath.mpf(Sr[1].numerator) / Sr[1].denominator)
                    bd = out["bound"][(a, b)]
                    assert abs(St - Srm) <= mpmath.mpf(bd.numerator) / bd.denominator


def test_s_bound_none_when_reference_not_validated():
    gam = [F(1, 2), F(3, 2)]
    D = [cm([[(F(1, 3), 0), (0, 0)], [(0, 0), (F(-1, 5), 0)]])] * 2
    A = C.assemble_A_mid(D, gam, "full")
    bad = R.reference(A, F(1, 10**30))
    assert bad["status"] == "INVALID"
    out = R.s_reference(bad, gam, 2, 2)
    assert all(v is None for v in out["bound"].values())


def test_refs_module_imports_no_solver_or_test_reference():
    tree = ast.parse((HERE / "coupling_scaling_v2_refs.py").read_text())
    mods = {(n.names[0].name if isinstance(n, ast.Import) else (n.module or "")).split(".")[0]
            for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))}
    assert mods <= {"__future__", "fractions", "coupling_scaling_v2_certify", "mpmath"}, mods
