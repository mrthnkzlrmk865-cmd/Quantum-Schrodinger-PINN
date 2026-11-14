"""Tests of the clean-room mpmath transcription wz_mutants_mp.py: correct code = direct inverses; every mutant's designated
block agrees with the §4.2 oracle (§4.3 oracle-vs-transcription check, 1e-50*||G||) on exact rational-lead lattices; the
§4.3 first-deviation example (M1); determinism; private context; import independence."""
import ast
import sys
import pathlib
from fractions import Fraction as F

import pytest

HERE = pathlib.Path(__file__).resolve().parent.parent / "experiments" / "wz_rgf_diagnostic"
sys.path.insert(0, str(HERE))
mpmath = pytest.importorskip("mpmath")

import wz_mutants_mp as T  # noqa: E402
import wz_mutants_oracle as O  # noqa: E402
from cs2_t6_reference import rlattice  # noqa: E402
from cs2_t6_reference.qexact import RC  # noqa: E402
from cs2_t6_harness import inputs as I  # noqa: E402
from cs2_t6_harness.common import ref_blocks, true_lams  # noqa: E402
from test_wz_mutants_oracle import EX  # noqa: E402

KEY = {"G_11": "Ghat11", "G_N1": "GhatN1", "G_1N": "Ghat1N", "G_NN": "GhatNN"}
CASES = [(name, N, s) for name in I.MAIN for N in I.NS[name] for s in (F(1), F(1, 2))]


@pytest.fixture(autouse=True)
def _workdps():
    """Test-side conversions and comparisons at 80 digits (the transcription itself uses its own private context)."""
    with mpmath.workdps(80):
        yield


def _q(x):
    return mpmath.mpf(x.numerator) / x.denominator


def _c(z):
    return mpmath.mpc(_q(z.re), _q(z.im))


def _inputs(name, N, s):
    lat = I.LATTICES[name]
    D = [[[_c(z) for z in r] for r in Dj] for Dj in rlattice.delta_blocks(lat, N, s)]
    lm1 = [_c(l - RC(1)) for l in true_lams(name)]
    return D, lm1


def _maxdiff(A_mp, B_exact):
    return max(abs(A_mp[i, j] - _c(B_exact[i][j])) for i in range(len(B_exact)) for j in range(len(B_exact[0])))


def _fro_exact(X):
    return mpmath.sqrt(sum(_q(z.abs2()) for r in X for z in r))


@pytest.mark.parametrize("name,N,s", CASES)
def test_transcription_correct_code_equals_direct_inverse(name, N, s):
    D, lm1 = _inputs(name, N, s)
    out = T.mp_endpoint_blocks(D, lm1, None, dps=60)
    b = ref_blocks(name, N, s)
    tol = mpmath.mpf(10) ** -50 * max(1, _fro_exact(b["Ghat"]))
    for k, rk in (("Ghat11", "G_11"), ("GhatN1", "G_N1"), ("Ghat1N", "G_1N"), ("GhatNN", "G_NN")):
        assert _maxdiff(out[k], b[rk]) <= tol, (k, name, N, s)


@pytest.mark.parametrize("name,N,s", CASES)
def test_oracle_vs_transcription_every_mutant(name, N, s):
    """§4.3: ||X_m^transcription - X_m^oracle|| <= 1e-50 ||G_full|| (here entrywise max, a stronger per-entry statement)."""
    ny = len(I.LATTICES[name]["gammas"])
    D, lm1 = _inputs(name, N, s)
    b = ref_blocks(name, N, s)
    Gp = {j: b["prefix"][j]["Ghat"] for j in range(1, N)}
    lm1_exact = [l - RC(1) for l in true_lams(name)]
    tol = mpmath.mpf(10) ** -50 * max(1, _fro_exact(b["Ghat"]))
    n = 0
    for m in O.MUTANTS:
        if not O.on_path(m, N) or m == "M11":
            continue
        orc = O.oracle(m, N, ny, b["Ghat"], Gp, lm1_exact, EX)
        tr = T.mp_endpoint_blocks(D, lm1, m, dps=60)
        assert _maxdiff(tr[KEY[orc["designated"]]], orc["X"]) <= tol, (m, name, N, s)
        for k, v in orc["secondary"].items():
            assert _maxdiff(tr[KEY[k]], v) <= tol, (m, k)
        n += 1
    assert n > 0


def test_first_deviation_M1_matches_prefix_closed_form():
    """§4.3 example: for M1 the first deviating intermediate is G_21^(2), equal to Ĝ^(1)_11 P_2 instead of P_2 Ĝ^(1)_11."""
    name, N, s = "X1", 3, F(1, 2)
    D, lm1 = _inputs(name, N, s)
    ok = dict(T.mp_endpoint_blocks(D, lm1, None, dps=60)["trace"])
    mu = dict(T.mp_endpoint_blocks(D, lm1, "M1", dps=60)["trace"])
    names = [nm for nm, _ in T.mp_endpoint_blocks(D, lm1, None, dps=60)["trace"]]
    def dev(A, B):   # entrywise (the matrices live in the transcription's private context)
        if (A.rows, A.cols) != (B.rows, B.cols):
            return mpmath.inf
        return max(abs(mpmath.mpc(A[i, j]) - mpmath.mpc(B[i, j])) for i in range(A.rows) for j in range(A.cols))
    assert list(ok) == list(mu) == names
    first = next(nm for nm in names if dev(mu[nm], ok[nm]) > mpmath.mpf(10) ** -50)
    assert first == "G_21^(2)"
    b = ref_blocks(name, N, s)
    G1 = b["G1_11"]
    P2 = b[("P", 2)]
    from cs2_t6_reference import rmat
    closed = rmat.matmul(G1, P2)
    assert _maxdiff(mu["G_21^(2)"], closed) <= mpmath.mpf(10) ** -50


def test_deterministic_private_context_and_invalid_labels():
    D, lm1 = _inputs("X2", 4, F(1, 2))
    dps_before = mpmath.mp.dps                        # 80 inside the fixture; the transcription must not change it
    a = T.mp_endpoint_blocks(D, lm1, "M5", dps=60)
    b = T.mp_endpoint_blocks(D, lm1, "M5", dps=60)
    assert mpmath.mp.dps == dps_before
    assert all(a[k][i, j] == b[k][i, j] for k in KEY.values() for i in range(2) for j in range(2))
    for bad in ("M0", "m1", 3, "", ("M1",)):
        with pytest.raises(ValueError):
            T.mp_endpoint_blocks(D, lm1, bad, dps=60)


def test_transcription_and_oracle_import_independence():
    for fname, forbidden in (("wz_mutants_mp.py", {"wz_mutants", "wz_rgf", "coupled_blocks", "numpy", "wz_mutants_oracle",
                                                   "cs2_t6_reference", "cs2_t6_harness"}),
                             ("wz_mutants_oracle.py", {"wz_mutants", "wz_rgf", "coupled_blocks", "wz_mutants_mp",
                                                       "cs2_t6_reference", "cs2_t6_harness"})):
        tree = ast.parse((HERE / fname).read_text())
        mods = {(n.names[0].name if isinstance(n, ast.Import) else (n.module or "")).split(".")[0]
                for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))}
        assert not (mods & forbidden), (fname, mods & forbidden)
        calls = {getattr(n.func, "id", getattr(n.func, "attr", None)) for n in ast.walk(tree) if isinstance(n, ast.Call)}
        assert not (calls & {"__import__", "import_module", "exec", "eval"}), fname
