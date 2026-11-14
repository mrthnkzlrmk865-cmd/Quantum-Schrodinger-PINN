"""Tests of wz_mutants_oracle.py: every §4.2 oracle formula against the independent exact recursion
(cs2_t6_reference/rmutants.py, a transcription of the code-change column) on exact rational-lead lattices, the M14
secondary blocks, off-path handling, and the mpmath backend."""
import sys
import pathlib
from fractions import Fraction as F

import pytest

HERE = pathlib.Path(__file__).resolve().parent.parent / "experiments" / "wz_rgf_diagnostic"
sys.path.insert(0, str(HERE))
mpmath = pytest.importorskip("mpmath")

import wz_mutants_oracle as O  # noqa: E402
from cs2_t6_reference import rmat, rinv  # noqa: E402
from cs2_t6_reference.qexact import RC  # noqa: E402
from cs2_t6_harness import inputs as I  # noqa: E402
from cs2_t6_harness.common import ref_blocks, recursion, true_lams  # noqa: E402


class ExactBackend:
    def mul(self, A, B):
        return rmat.matmul(A, B)

    def add(self, A, B):
        return rmat.add(A, B)

    def sub(self, A, B):
        return rmat.sub(A, B)

    def neg(self, A):
        return rmat.neg(A)

    def scal(self, c, A):
        return rmat.scale(RC(c), A)

    def inv(self, A):
        return rinv.inverse(A)

    def block(self, X, a, b, ny):
        return rmat.block(X, a, b, ny)

    def diag(self, vals):
        return rmat.diag(vals)


EX = ExactBackend()
KEY = {"G_11": "G11", "G_N1": "GN1", "G_1N": "G1N", "G_NN": "GNN"}
CASES = [(name, N, s) for name in I.MAIN for N in I.NS[name] for s in (F(1), F(1, 2))]


def _inputs(name, N, s):
    b = ref_blocks(name, N, s)
    return b["Ghat"], {j: b["prefix"][j]["Ghat"] for j in range(1, N)}, [l - RC(1) for l in true_lams(name)]


@pytest.mark.parametrize("name,N,s", CASES)
def test_every_mutant_oracle_equals_exact_recursion(name, N, s):
    ny = len(I.LATTICES[name]["gammas"])
    Gf, Gp, lm1 = _inputs(name, N, s)
    checked = 0
    for m in O.MUTANTS:
        if not O.on_path(m, N) or m == "M11":
            with pytest.raises(O.OffPath):
                O.oracle(m, N, ny, Gf, Gp, lm1, EX)
            continue
        res = O.oracle(m, N, ny, Gf, Gp, lm1, EX)
        rec = recursion(name, N, s, m)
        assert res["X"] == rec[KEY[res["designated"]]], (m, N, name, s)
        for k, v in res["secondary"].items():
            assert v == rec[KEY[k]], (m, k)
        checked += 1
    assert checked == sum(1 for m in O.MUTANTS if O.on_path(m, N) and m != "M11")


def test_identities_named_in_section_4_2():
    """M12 equals the reference; M5 at N = 3 and M8 at N = 2 equal the reference (identical by algebra)."""
    for name in I.MAIN:
        ny = len(I.LATTICES[name]["gammas"])
        for N, m in ((3, "M5"), (2, "M8"), (2, "M12"), (3, "M12")):
            Gf, Gp, lm1 = _inputs(name, N, F(1, 2))
            assert O.oracle(m, N, ny, Gf, Gp, lm1, EX)["X"] == O.reference_block(O.DESIGNATED[m], N, ny, Gf, EX)


def test_off_path_is_never_a_zero_discrepancy():
    Gf, Gp, lm1 = _inputs("X1", 2, F(1, 2))
    for m in ("M1", "M2", "M3", "M4", "M5", "M9", "M13", "M11"):
        with pytest.raises(O.OffPath):
            O.oracle(m, 2, 2, Gf, Gp, lm1, EX)
    with pytest.raises(ValueError):
        O.oracle("M99", 3, 2, Gf, Gp, lm1, EX)
    with pytest.raises(ValueError):
        O.oracle("M10", 2, 2, Gf, Gp, None, EX)


def test_mp_backend_matches_exact_backend():
    ctx = mpmath.MPContext()
    ctx.dps = 60
    be = O.MpBackend(ctx)

    def to_mp(X):
        M = ctx.matrix(len(X), len(X[0]))
        for i, r in enumerate(X):
            for j, z in enumerate(r):
                M[i, j] = ctx.mpc(ctx.mpf(z.re.numerator) / z.re.denominator, ctx.mpf(z.im.numerator) / z.im.denominator)
        return M

    for name, N in (("X1", 4), ("X2", 4)):
        ny = 2
        Gf, Gp, lm1 = _inputs(name, N, F(1, 2))
        b = ref_blocks(name, N, F(1, 2))
        Gf_mp = ctx.mpf(32) * ctx.inverse(to_mp(b["A"]))
        Gp_mp = {j: ctx.mpf(32) * ctx.inverse(to_mp(b["prefix"][j]["A"])) for j in Gp}
        lm_mp = [ctx.mpc(ctx.mpf(z.re.numerator) / z.re.denominator, ctx.mpf(z.im.numerator) / z.im.denominator) for z in lm1]
        for m in O.MUTANTS:
            if not O.on_path(m, N) or m == "M11":
                continue
            Xe = to_mp(O.oracle(m, N, ny, Gf, Gp, lm1, EX)["X"])
            Xm = O.oracle(m, N, ny, Gf_mp, Gp_mp, lm_mp, be)["X"]
            assert ctx.mnorm(Xm - Xe, "f") <= ctx.mpf(10) ** -50 * max(1, ctx.mnorm(Gf_mp, "f")), (m, name)


def test_monomial_counts_match_section_4_2_examples():
    assert sum(c for c, _ in O.monomials("M13", 4)) == 1 + 2 * 2 + 1   # Ĝ_11 + 2 theta_2 + 2 theta_3, plus reference
    assert max(d for _, d in O.monomials("M16", 3)) == 5               # psi Ĝ_NN^3 chi
    assert max(d for _, d in O.monomials("M1", 6)) == 6
    with pytest.raises(O.OffPath):
        O.monomials("M10", 3)
