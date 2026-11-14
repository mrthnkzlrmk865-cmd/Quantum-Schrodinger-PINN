"""Driver (coupling_scaling_v2) and §5 analysis (coupling_scaling_v2_analysis): gate refusal, fixture pipeline and output
schema, first-deviation logic, decision precedence, analysis procedures. FIXTURES ONLY: configurations come from a
non-frozen test root; the XP/HO manifests are never read and no main-experiment run is performed."""
import csv
import json
import math
import pathlib
import sys
from fractions import Fraction as F

import pytest

HERE = pathlib.Path(__file__).resolve().parent.parent / "experiments" / "wz_rgf_diagnostic"
sys.path.insert(0, str(HERE))
mpmath = pytest.importorskip("mpmath")
import coupling_scaling_v2_analysis as AN  # noqa: E402
import wz_mutants_oracle as O  # noqa: E402

TEST_ROOT = 12345


@pytest.fixture(scope="module")
def fixture_cfg():
    import coupling_scaling_v2_manifest as M
    M.ROOTS["TEST"] = TEST_ROOT
    try:   # no solver_guard here: the driver is already imported in this test process (the guard is tested separately)
        rec = M.generate_configuration("TEST", 0, 2, "R", 3, 0, max_attempts=200)
    finally:
        M.ROOTS.pop("TEST", None)
    import coupling_scaling_v2 as D
    return D.config_from_manifest(json.loads(M.dumps(rec)))


@pytest.fixture(scope="module")
def D():
    import coupling_scaling_v2 as D
    return D


# ------------------------------------------------------------------------------------------------------------ gate
def test_gate_refuses_in_current_state(D):
    with pytest.raises(D.GateRefused) as e:
        D.gate()
    msg = str(e.value)
    assert "duplicate set" in msg and "REVIEW_coupling_scaling_v2.md" in msg


def test_main_refuses_before_any_solver_import(D, monkeypatch, tmp_path):
    called = []
    monkeypatch.setattr(D, "_solvers", lambda: called.append(1))
    with pytest.raises(D.GateRefused):
        D.main(["--set", "XP", "--out", str(tmp_path)])
    assert not called and not list(tmp_path.iterdir())


def _fake_out(tmp_path):
    (tmp_path / "manifest_XP.json").write_text("{}\n")
    import hashlib
    h = hashlib.sha256(b"{}\n").hexdigest()
    (tmp_path / "manifest.sha256").write_text(f"{h}  manifest_XP.json\n")
    (tmp_path / "previous_files.sha256").write_text("")
    return h


def test_gate_detects_hash_mismatch_review_and_gap(D, tmp_path, monkeypatch):
    h = _fake_out(tmp_path)
    rv = tmp_path / "REVIEW.md"
    rv.write_text(f"R0 SIGN-OFF\n{h}\n")
    monkeypatch.setattr(D, "DUPLICATE_SET", frozenset())
    assert D.gate(str(tmp_path), str(rv))
    (tmp_path / "manifest_XP.json").write_text("{ }\n")
    with pytest.raises(D.GateRefused, match="hash mismatch"):
        D.gate(str(tmp_path), str(rv))
    _fake_out(tmp_path)
    rv.write_text("R0 SIGN-OFF\n")
    with pytest.raises(D.GateRefused, match="does not quote"):
        D.gate(str(tmp_path), str(rv))
    rv.write_text(f"{h}\n")
    with pytest.raises(D.GateRefused, match="SIGN-OFF"):
        D.gate(str(tmp_path), str(rv))
    monkeypatch.setattr(D, "DUPLICATE_SET", None)
    rv.write_text(f"R0 SIGN-OFF\n{h}\n")
    with pytest.raises(D.GateRefused, match="duplicate set"):
        D.gate(str(tmp_path), str(rv))


def test_gate_detects_changed_previous_file(D, tmp_path, monkeypatch):
    h = _fake_out(tmp_path)
    (tmp_path / "previous_files.sha256").write_text("0" * 64 + "  wz_rgf.py\n")
    rv = tmp_path / "REVIEW.md"
    rv.write_text(f"R0 SIGN-OFF\n{h}\n")
    monkeypatch.setattr(D, "DUPLICATE_SET", frozenset())
    with pytest.raises(D.GateRefused, match="previous file changed"):
        D.gate(str(tmp_path), str(rv))


# --------------------------------------------------------------------------------------------- fixture pipeline
@pytest.fixture(scope="module")
def fixture_run(D, fixture_cfg):
    grid = [F(1, 2**40), F(1, 2**12), F(3, 8), F(1, 2), F(1)]
    return D.run_configuration(fixture_cfg, grid)


def test_fixture_pipeline_valid_and_schema(D, fixture_run):
    r = fixture_run
    assert r["invalid"] == []
    p = r["points"][-1]
    assert p["faithful"] and p["item1_ok"] and not p["aborted"]
    assert set(p["A"]["blocks"]) == set(D.BLOCKS) and p["A"]["dS_ref_bound"] is not None
    for m in O.MUTANTS:
        mr = p["mutants"][m]
        if m == "M11":
            assert mr["status"] == "off_path"
            continue
        assert mr["status"] == "computed"
        assert mr["agree_oracle_transcription"]["R_W"] <= 1e-50 * p["normG"]["R_W"]
        assert mr["first_deviation"]["ok"]
        assert set(mr["blocks"]) == set(D.BLOCKS)
    for lat in ("L", "R_W", "R_orig"):
        assert all(v["dps"] in (60, 120) for v in p["ref_meta"][lat].values())
    s1 = p["ref_meta"]["R_W"]["full"]["checks"]
    assert s1["b"]["ok"] and s1["c"]["ok"]                    # s = 1 is in both check sets


def test_fixture_pair_categories_and_orders(fixture_run):
    pr = fixture_run["pairs"]
    for m in ("M1", "M2", "M3", "M6", "M7", "M8", "M15"):
        assert pr[m]["category"] == "ordering-meaningful" and pr[m]["q"] == 1
        assert pr[m]["series_check"]["status"] in ("pass", "pass_after_escalation")
    assert pr["M5"]["category"] == "identical_by_algebra"      # N = 3
    assert pr["M9"]["category"] == "perturbation-control" and pr["M9"]["q"] >= 1   # Fact 6: C_0 = 0
    assert pr["M12"]["category"] == "equivalent" and pr["M12"]["order"][0] == "all_zero"
    assert pr["M11"]["category"] == "off_path"
    for m in ("M4", "M4L", "M4R", "M10", "M13", "M14", "M16"):
        assert pr[m]["category"] == "structural" and pr[m]["hyp4_structural"]


def test_fixture_outputs_written(D, fixture_run, tmp_path):
    dec = D.decide("XP", [fixture_run])
    D.write_outputs("XP", [fixture_run], dec, str(tmp_path))
    D.write_report([dec], str(tmp_path))
    with open(tmp_path / "points_XP.csv") as f:
        rows = list(csv.DictReader(f))
    assert rows and set(rows[0]) == set(D.POINT_COLS)
    with open(tmp_path / "pairs_XP.csv") as f:
        prs = list(csv.DictReader(f))
    assert {r["mutant"] for r in prs} == set(O.MUTANTS) and set(prs[0]) == set(D.PAIR_COLS)
    js = json.loads((tmp_path / "results_XP.json").read_text())
    assert js["decision"]["decides"] is False
    assert "XP can never yield PASS" in json.dumps(js["decision"]) or js["decision"]["scopes"]["physical"]["decision"] != "PASS"
    assert (tmp_path / "coupling_scaling_v2_report.md").exists()


def test_determinism_on_fixture(D, fixture_cfg):
    assert D.determinism([fixture_cfg], D._solvers(), grid=[F(1), F(1, 2**20)])


# --------------------------------------------------------------------------------------------- first deviation
def test_first_deviation_logic(D):
    with mpmath.workdps(60):
        a, b = mpmath.matrix([[1]]), mpmath.matrix([[2]])
        ok = [("X", a), ("Y", a), ("Z", a)]
        mut = [("X", a), ("Y", b), ("Z", b)]
        tol = mpmath.mpf(10) ** -50
        assert D.first_deviation(mut, ok, ("Y", b), tol)["ok"]
        assert not D.first_deviation(mut, ok, ("Z", b), tol)["ok"]          # wrong slice
        assert not D.first_deviation(mut, ok, ("Y", a), tol)["ok"]          # wrong value
        assert D.first_deviation(ok, ok, (None, None), tol)["ok"]           # identity: nothing may deviate
        assert not D.first_deviation(mut, ok, (None, None), tol)["ok"]


def test_category_rules(D):
    c = lambda N, cls: dict(N=N, cls=cls)  # noqa: E731
    assert D.category("M11", c(2, "R")) == "off_path" and D.category("M1", c(2, "R")) == "off_path"
    assert D.category("M5", c(3, "R")) == "identical_by_algebra" and D.category("M8", c(2, "H")) == "identical_by_algebra"
    assert D.category("M12", c(4, "R")) == "equivalent" and D.category("M1", c(4, "D")) == "class_D_order"
    assert D.category("M9", c(3, "D")) == "perturbation-control" and D.category("M16", c(2, "R")) == "structural"
    assert D.category("M1", c(3, "H")) == "candidate"


def _res(cls="R", invalid=(), pairs=None, pts=None):
    return dict(config=dict(set="HO", idx=0, ny=2, cls=cls, N=3, n_open=1, E="1"), invalid=list(invalid),
                points=pts or [], pairs=pairs or {}, noise_trend={})


def _pt(e_glob=1e-16, ref_bound=1e-60, tau=1e-13, judge=True):
    blocks = {b: dict(e_glob=e_glob, ref_bound=ref_bound) for b in ("G_11", "G_N1", "G_1N", "G_NN")}
    A = dict(blocks=blocks, dS_max=1e-16, dS_ref_bound=1e-60, tau_S=1e-12, intermediates=[], U=0, U_bound=1, Q_value_judged=0,
             Q_bound=1, sym_applicable=False, sym_residual=0, sym_bound=1, tau_end=tau)
    return dict(s="1", aborted=False, A=A, judgeable=judge, tau_end=tau)


def test_decision_precedence(D):
    assert D.decide("HO", [_res(invalid=["x"])])["scopes"]["physical"]["decision"] == "INVALID"
    fw = D.decide("HO", [_res(pts=[_pt(e_glob=1e-9)])])["scopes"]["physical"]
    assert fw["decision"] == "FAIL-W"
    nr = D.decide("HO", [_res(pts=[_pt(e_glob=1e-9, ref_bound=None)])])["scopes"]["physical"]
    assert nr["decision"] == "DESIGN INSUFFICIENT" and "not reliably measured" in nr["design_insufficient"][0]
    an = D.decide("HO", [_res(pts=[_pt(e_glob=1e-12, tau=1e-13)])])["scopes"]["physical"]
    assert any("anomaly" in x for x in an["design_insufficient"])
    nj = D.decide("HO", [_res(pts=[_pt(tau=1e-10, judge=False)])])["scopes"]["physical"]
    assert any("non-judgeable" in x for x in nj["design_insufficient"])
    fdr = D.decide("HO", [_res(pairs={"M1": dict(category="ordering-meaningful", regression=dict(valid=True),
                                                  window_adequate=True, hyp2=True, hyp1=False)})])["scopes"]["physical"]
    assert fdr["fail_d"] and fdr["decision"] in ("FAIL-D",)
    assert D.decide("HO", [], determinism_ok=False)["scopes"]["algebraic"]["decision"] == "INVALID"
    assert D.decide("XP", [])["decides"] is False


def test_driver_imports_solvers_only_lazily(D):
    import ast
    tree = ast.parse((HERE / "coupling_scaling_v2.py").read_text())
    top = {(n.names[0].name if isinstance(n, ast.Import) else (n.module or "")).split(".")[0] for n in tree.body
           if isinstance(n, (ast.Import, ast.ImportFrom))}
    assert not (top & {"wz_rgf", "wz_mutants", "stationary", "coupled_blocks"}), top


# ------------------------------------------------------------------------------------------------------ analysis
def test_zero_test_thresholds():
    with mpmath.workdps(220):
        mf = mpmath.mpf
        assert AN.zero_test(mf("1e-3"), mf("1e-3") * (1 + mf("1e-25"))) == "nonzero"
        assert AN.zero_test(mf("1e-3"), mf("1.1e-3")) == "unresolved"
        assert AN.zero_test(mf("1e-90"), mf("1e-180")) == "zero"
        assert AN.zero_test(mf("1e-79"), mf("1e-180")) == "unresolved"
        assert AN.zero_test(mf("1e-59"), mf("1e-59")) == "nonzero" and AN.zero_test(mf("1e-61"), mf("1e-61")) == "unresolved"
    assert AN.identify_order(["zero", "zero", "nonzero"]) == ("q", 2)
    assert AN.identify_order(["zero", "unresolved", "nonzero"]) == ("unresolved", 1)
    assert AN.identify_order(["zero"] * 13) == ("all_zero", None)


def test_b_tail_exact_and_halving():
    B, s, x, z = AN.b_tail(3, 2, 1, F(1, 4), F(1))        # x = 1/4, zeta = 15/14 -> zeta x < 1
    assert s == 1 and x == F(1, 4) and z == F(15, 14)
    assert B == 3 * math.comb(14, 1) * F(1, 4) ** 13 / (1 - F(15, 14) * F(1, 4))
    B2, s2, x2, _ = AN.b_tail(1, 30, 1, F(1, 2), F(1))     # zeta = 43/14; needs halving
    assert s2 < 1 and F(43, 14) * x2 < 1
    assert AN.b_tail(5, 3, 2, 0, F(1))[0] == 0
    assert AN.s_star("+inf") == F(1, 2**60) and AN.s_star(F(1, 2**70)) == F(1, 2**90)


def test_series_backend_matches_direct_inverse():
    ctx = mpmath.MPContext()
    ctx.dps = 50
    A0 = ctx.matrix([[3, 1], [0.5, 4]])
    A1 = ctx.matrix([[0.2, -0.1], [0.3, 0.1]])
    S = AN.neumann_series(ctx, A0, A1, 32, K=40)
    s = ctx.mpf(1) / 16
    approx = sum((S[k] * s**k for k in range(1, 41)), S[0])
    exact = 32 * ctx.inverse(A0 + s * A1)
    assert max(abs(approx[i, j] - exact[i, j]) for i in range(2) for j in range(2)) < ctx.mpf(10) ** -40
    be = AN.SeriesBackend(ctx, K=3)
    X = be.mul([A0, A1, A0 * 0, A0 * 0], [A1, A0, A0 * 0, A0 * 0])
    assert X[1] == A0 * A0 + A1 * A1 and X[0] == A0 * A1


def test_ols_slopes_crossover_and_noise():
    svals = [F(1, 2**k) for k in range(40, -1, -1)]
    D = [2 * float(s) ** 2 for s in svals]
    o = AN.ols([math.log(float(s)) for s in svals], [math.log(d) for d in D])
    assert abs(o["slope"] - 2) < 1e-12 and abs(sum(o["weights"])) < 1e-12
    phi = AN.local_slopes(svals, D)
    assert all(abs(v - 2) < 1e-9 for v in phi.values()) and F(1) in phi
    sx, sec = AN.crossover(svals, D, 2.0, 2, phi)
    assert sx == F(1) and sec == F(1)
    pts = [(s, d, d, 1e6) for s, d in zip(svals, D)]
    reg = AN.regression(pts, 2)
    assert reg["valid"] and abs(reg["p_hat"] - 2) < 1e-9 and reg["hyp2"] and reg["interval"] == 0
    assert AN.window_adequacy(reg, list(phi.values()), 2)
    nt = AN.noise_trend(svals, [1e-16] * len(svals), [1e-17] * len(svals))
    assert nt["ok"] and nt["n"] == 29
    assert AN.regression([(s, d, d, 1.0) for s, d in zip(svals, D)], 2)["n_points"] == 0


def test_bootstrap_deterministic():
    a = AN.bootstrap_median([1.0, 2.0, 3.0, 4.0], key=0, n=200)
    b = AN.bootstrap_median([1.0, 2.0, 3.0, 4.0], key=0, n=200)
    assert a == b and a["median"] == 2.5


def test_magnitude_monomials_consistent():
    """monomial_factors has the same (coef, degree) structure as monomials, for every polynomial mutant and N."""
    class Sym:
        def block(self, X, a, b, ny): return ("blk", a, b)
        def neg(self, A): return A
        def mul(self, A, B): return ("mul",)
    for m in O.MUTANTS:
        for N in (2, 3, 4, 6):
            if not O.on_path(m, N) or m in ("M11", "M10"):
                continue
            facs = O.monomial_factors(m, N, 2, None, {j: None for j in range(1, N)}, Sym())
            assert sorted((c, len(f)) for c, f in facs) == sorted(O.monomials(m, N)), (m, N)
