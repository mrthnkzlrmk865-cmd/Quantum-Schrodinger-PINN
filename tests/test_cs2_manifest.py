"""Manifest generator (coupling_scaling_v2_manifest): §2.6 order, §2.3 draws, §2.4 energy/conditioning/certificate
plumbing, redraw cap and abort, solver guard, determinism and hash. Uses NON-frozen test roots only (never the XP/HO
roots), so no manifest configuration is generated here."""
import ast
import hashlib
import pathlib
import sys
from fractions import Fraction as F

import numpy as np
import pytest

HERE = pathlib.Path(__file__).resolve().parent.parent / "experiments" / "wz_rgf_diagnostic"
sys.path.insert(0, str(HERE))
import coupling_scaling_v2_manifest as M  # noqa: E402

TEST_ROOT = 987654321


def _fresh(code):
    """Run code in a fresh interpreter (the solver guard checks sys.modules, which other test files populate)."""
    import subprocess
    r = subprocess.run([sys.executable, "-c", "import sys; sys.path.insert(0, %r)\n" % str(HERE) + code],
                       cwd=str(HERE), capture_output=True, text=True, timeout=1200)
    assert r.returncode == 0, r.stderr[-3000:]
    return r.stdout


@pytest.fixture
def test_root(monkeypatch):
    monkeypatch.setitem(M.ROOTS, "TEST", TEST_ROOT)
    monkeypatch.setitem(M.DRAWS, "TEST", {"R": 1, "H": 1, "D": 1})
    return "TEST"


def test_index_order_and_counts():
    xp, ho = M.configurations("XP"), M.configurations("HO")
    assert len(xp) == 60 and len(ho) == 120
    assert [c[0] for c in ho] == list(range(120))
    assert xp[:3] == [(0, 2, "R", 2, 0), (1, 2, "R", 2, 1), (2, 2, "R", 3, 0)]
    assert xp[8] == (8, 2, "H", 2, 0) and xp[16] == (16, 2, "D", 2, 0) and xp[20] == (20, 3, "R", 2, 0)
    for s, n in (("XP", {"R": 24, "H": 24, "D": 12}), ("HO", {"R": 48, "H": 48, "D": 24})):
        cf = M.configurations(s)
        assert {k: sum(1 for c in cf if c[2] == k) for k in "RHD"} == n
    assert M.ROOTS["XP"] == 202611110001 and M.ROOTS["HO"] == 202611110002 and M.BOOTSTRAP_ROOT == 202611110003


def test_draw_order_matches_spec_generator():
    n_open, u, Jr, Ji = M.draw(TEST_ROOT, 7, 3, 3, 4)
    g = np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy=TEST_ROOT, spawn_key=(7, 3))))
    assert n_open == int(g.integers(1, 4)) and u == float(g.random())
    for j in range(4):
        assert Jr[j] == g.integers(-8, 9, size=(3, 3)).tolist() and Ji[j] == g.integers(-8, 9, size=(3, 3)).tolist()
    assert M.draw(TEST_ROOT, 7, 3, 3, 4) == (n_open, u, Jr, Ji)
    assert M.draw(TEST_ROOT, 7, 4, 3, 4) != (n_open, u, Jr, Ji)


@pytest.mark.parametrize("cls", "RHD")
def test_v1_blocks_structure(cls):
    _, _, Jr, Ji = M.draw(TEST_ROOT, 1, 0, 4, 2)
    V = M.v1_blocks(4, cls, Jr, Ji)
    for B in V:
        for a in range(4):
            for b in range(4):
                re, im = B[a][b]
                assert (re, -im) == B[b][a]                       # Hermitian
                assert (2 * re).denominator == 1 and -8 <= re <= 8 and (2 * im).denominator == 1
                if cls != "H":
                    assert im == 0
                if cls == "D" and a != b:
                    assert re == 0
        assert B[0][0][0] == F(Jr[0][0][0]) if B is V[0] else True


def test_energy_rule_exact_rhe_and_window():
    ny = 3
    Ep = [F(e) for e in M.eperp_doubles(ny)]
    for n_open in (1, 2, 3):
        for u in (0.0, 0.25, 0.5, 0.999999):
            E, lo, up, d, ok, why = M.energy(ny, n_open, u)
            assert E.denominator in (1, 2, 4, 8, 16, 32, 64)
            raw = 64 * (lo + F(u) * (up - lo))
            assert E * 64 == round(raw) and abs(E * 64 - raw) <= F(1, 2)
            assert up == (Ep[n_open] if n_open < ny else Ep[0] + 64 * (2 - F(1, 50)))
            assert d == [(E - e) / 64 for e in Ep]
            assert ok == all(abs(x) >= F(1, 50) and x <= 2 - F(1, 50) for x in d) or why == \
                "energy: open-channel count differs from n_open"
    assert M.energy(ny, 1, 0.0)[4] is False                      # E = E_perp,1 -> delta_1 = 0 fails


def test_rhe_half_to_even():
    assert round(F(5, 2)) == 2 and round(F(7, 2)) == 4 and round(F(-5, 2)) == -2


def test_conditioning_quantities_present_and_bounded(test_root):
    ny, cls, N = 2, "R", 3
    for attempt in range(50):
        n_open, u, Jr, Ji = M.draw(TEST_ROOT, 0, attempt, ny, N)
        E, *_, ok, _ = M.energy(ny, n_open, u)
        if ok:
            break
    V1 = M.v1_blocks(ny, cls, Jr, Ji)
    c = M.conditioning_at("L", E, ny, cls, V1, F(1, 2), N)
    assert set(c) == {"M", "M^(1)", "M^(2)", "SC_1", "SC_2", "Z", "Z-(Lambda-I)"}
    c2 = M.conditioning_at("R_W", E, ny, cls, V1, F(1, 2), N)
    assert set(c2) == set(c) and all(v >= 1 for v in c2.values())


def test_conditioning_rejects_singular_and_nonfinite():
    assert M._cond(np.zeros((2, 2))) == float("inf")
    assert M._cond(np.array([[1.0, 0], [0, 1e-20]])) > 1e3


def test_solver_guard_blocks_and_restores():
    out = _fresh("""
import coupling_scaling_v2_manifest as M
stationary = M._transverse_basis()
orig = stationary.scattering_matrix
with M.solver_guard():
    for f in (lambda: stationary.scattering_matrix(1.0, None), lambda: stationary.transverse_basis(1.0, 2)):
        try:
            f(); print('NOT BLOCKED')
        except M.SolverGuardViolation:
            print('blocked')
    assert M.eperp_doubles(2)
assert stationary.scattering_matrix is orig
print('restored')
""")
    assert out.split() == ["blocked", "blocked", "restored"]


def test_solver_guard_detects_forbidden_module(monkeypatch):
    monkeypatch.setitem(sys.modules, "wz_rgf", object())
    with pytest.raises(M.SolverGuardViolation):
        with M.solver_guard():
            pass


def test_manifest_module_imports_no_solver():
    tree = ast.parse((HERE / "coupling_scaling_v2_manifest.py").read_text())
    mods = {(n.names[0].name if isinstance(n, ast.Import) else (n.module or "")).split(".")[0]
            for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))}
    assert not (mods & set(M.FORBIDDEN_MODULES)), mods
    assert mods <= {"__future__", "hashlib", "json", "os", "platform", "sys", "fractions", "numpy",
                    "coupling_scaling_v2_certify", "coupling_scaling_v2_covering", "stationary", "mpmath",
                    "multiprocessing", "argparse"}, mods
    calls = {getattr(n.func, "id", getattr(n.func, "attr", None)) for n in ast.walk(tree) if isinstance(n, ast.Call)}
    assert not (calls & {"__import__", "import_module", "exec", "eval"})


def test_generate_configuration_deterministic_and_hash():
    code = """
import hashlib, json, sys
import coupling_scaling_v2_manifest as M
M.ROOTS['TEST'] = %d
txt = M.dumps(M._worker(('TEST', 0, 2, 'R', 2, 0, 50)))
assert 'wz_rgf' not in sys.modules and 'wz_mutants' not in sys.modules and 'coupled_blocks' not in sys.modules
print(hashlib.sha256(txt.encode()).hexdigest())
print(txt)
""" % TEST_ROOT
    a, b = _fresh(code), _fresh(code)
    assert a == b
    import json
    r = json.loads(a.split("\n", 1)[1])
    cert = r["certificate"]
    assert cert["accepted"] and set(cert["matrices"]) == {"M_L", "M_L^(1)"} and cert["r0"] is not None
    for rec in cert["matrices"].values():
        assert rec["covered"] and rec["points"][0]["s"] == "0/1" and rec["r0"] is not None
        ss = [F(p["s"]) for p in rec["points"]]
        assert ss == sorted(set(ss))
    assert all(isinstance(x["reason"], str) for x in r["rejection_log"])
    assert r["set"] == "TEST" and r["idx"] == 0 and len(r["Jr"]) == 2


def test_redraw_cap_aborts_without_partial_output(test_root, monkeypatch, tmp_path):
    monkeypatch.setattr(M, "energy", lambda ny, n, u: (F(0), F(0), F(0), [], False, "energy: forced (test)"))
    with pytest.raises(M.ManifestAbort):
        M.generate_configuration("TEST", 0, 2, "R", 2, 0, max_attempts=3)
    monkeypatch.setitem(M.ROOTS, "XP", TEST_ROOT)
    import contextlib
    monkeypatch.setattr(M, "solver_guard", contextlib.nullcontext)     # abort logic only; the guard is tested separately
    with pytest.raises(M.ManifestAbort):
        M.build_manifest("XP", indices=[0], max_attempts=2)
    assert not (tmp_path / "manifest_XP.json").exists()


def test_main_refuses_to_overwrite_frozen_manifest(tmp_path):
    (tmp_path / "manifest.sha256").write_text("x\n")
    with pytest.raises(SystemExit):
        M.main(["--out", str(tmp_path)])


def test_previous_files_listing_excludes_outputs_and_pycache():
    txt = M.write_previous_files()
    paths = [ln.split("  ", 1)[1] for ln in txt.splitlines()]
    assert paths == sorted(paths)
    assert not any("__pycache__" in p or p.startswith("results_coupling_scaling_v2") for p in paths)
    assert "wz_rgf.py" in paths and "wz_mutants.py" in paths


def test_previous_files_step_is_separate_and_refuses_overwrite(tmp_path):
    assert M.main(["--out", str(tmp_path), "--previous-files"]) == 0
    assert (tmp_path / "previous_files.sha256").exists() and not (tmp_path / "manifest.sha256").exists()
    with pytest.raises(SystemExit):
        M.main(["--out", str(tmp_path), "--previous-files"])
