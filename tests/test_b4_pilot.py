"""B4 limited diagnostic pilot — tests, phases T1–T3 (including T3.1) and T5.

Classification: DIAGNOSTIC / NOT B4. Nothing here is a B4 result, and no test marks any pilot gate as passed.

Protocol: docs/B4_PILOT_PROTOCOL.md as committed in the immutable baseline commit BASELINE below (decision A6).
Expected values are parsed from that committed protocol text, so every check is traceable to the protocol.
No threshold is defined or changed here.

T1  design-time arithmetic and protected-artifact integrity:
    §1.2 point table, §4.1 moved-node table, §2.4 frame-time equality, §3.3.1 dx-equality constants,
    protected source files / documents / B1–B3Y (+ validation, transport) artifacts against BASELINE.
T2  same realization: L1 disorder-parameter hashes (RGF path, TD production path, B3/B3X records) and
    L2 sampler fidelity on the production TD grid (bitwise). Shifted variants and L3 are NOT implemented here.
T3  b4_pilot_rgf: §1.4 observed-order classification, §1.6 admissibility and deterministic resolution selection,
    and (T3.1) the §1.4 point flag (status mismatch OR |ΔT1 + ΔR1| above the recorded unitarity residuals).
    Synthetic inputs only; passing these tests does NOT mean gate G-RGF passed.
T5  b4_pilot_io: §5.1–§5.5 output guard (allowlist/denylist), DIAGNOSTIC / NOT B4 labels, metadata (pinned protocol
    SHA-256 reference separate from the A6 baseline commit; HEAD; per-file code fingerprints), clean-tree abort and
    atomic publication via staging outside the repository. Publication tests use a throwaway git repository under
    tmp_path; the real results/ directory is never written. No G0/G-REAL is executed.

No solver propagation, RGF energy, or pilot driver is run. Git must be available: the baseline is read from the
repository object database; a missing git or commit is a test failure, never a skip.
"""
import hashlib
import inspect
import json
import os
import re
import subprocess
import sys

import numpy as np
import pytest

import validate_rgf_b2 as b2
import validate_rgf_b3 as b3
from config import Config, waveguide_config
from disorder import disorder_from_config
from potential import Potential
from reference_solver import SplitStepDST, build_grid

CLASSIFICATION = "DIAGNOSTIC / NOT B4"
BASELINE = "f3ab9accf04e8a4927fac3f276cd88910ac1bdd3"
PROTOCOL = "docs/B4_PILOT_PROTOCOL.md"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# Protected source files and documents (Phase 0 audit §5; B3Y protected list extended to the TD stack and the
# B4 documents). provenance.py stays protected (decision A1: no change to it).
PROTECTED_SOURCES = (
    "stationary.py", "smatrix_spec.py", "sparse_reference.py", "sparse_wavefunction.py", "sparse_exact_lattice.py",
    "reference_solver.py", "transport.py", "potential.py", "disorder.py", "config.py", "initial_condition.py",
    "observables.py", "analytic.py", "ode1d.py", "provenance.py", "experiment_disorder.py",
    "validate_solver.py", "validate_rgf_b1.py", "validate_rgf_b2.py", "validate_rgf_b3.py", "validate_rgf_b3x.py",
    "validate_rgf_b3y.py", "audit_rgf_phase.py",
)
PROTECTED_DOCS = (
    "docs/PHASE2_METHOD_SELECTION.md", "docs/RGF_B1_VALIDATION.md", "docs/RGF_B2_VALIDATION.md",
    "docs/RGF_B3_VALIDATION.md", "docs/RGF_B3X_VALIDATION.md", "docs/RGF_B3Y_VALIDATION.md",
    "docs/RGF_B3Y_TIER3_AUDIT.md", "docs/rgf_b3x_protocol.json", "docs/rgf_b3y_protocol.json",
    "docs/B4_IMPLEMENTATION_READINESS_AUDIT.md", PROTOCOL,
)
# Protected result trees: B1–B3Y artifacts (protocol §5.3) and the two other official directories of the §5.3 denylist.
PROTECTED_RESULT_DIRS = (
    "results/rgf_b1_validation", "results/rgf_b2_validation", "results/rgf_b3_validation",
    "results/rgf_b3x_validation", "results/rgf_b3y_validation", "results/validation", "results/transport",
)

# The four B4 realizations (B3 cases; protocol §1.1) and the committed production TD files of the same realizations.
CASES = ("A", "B", "C", "D")
TD_FILES = {"A": "results/transport/corr_length/zeta0.02.npz", "B": "results/transport/corr_length/zeta0.05.npz",
            "C": "results/transport/corr_length/zeta0.2.npz", "D": "results/transport/strength/W200.npz"}
B3X_NAMES = {"A": "W100_z0.02", "B": "W100_z0.05", "C": "W100_z0.2", "D": "W200_z0.05"}


# ----------------------------------------------------------------------
# baseline access (git object database; no network; failures are failures)
# ----------------------------------------------------------------------
def _git(*args, input_bytes=None):
    try:
        r = subprocess.run(["git", "-C", REPO, *args], input=input_bytes, capture_output=True, check=False)
    except FileNotFoundError as exc:                                 # pragma: no cover - environment failure
        pytest.fail(f"git is required to read the A6 baseline {BASELINE}: {exc}")
    if r.returncode != 0:
        pytest.fail(f"git {' '.join(args)} failed: {r.stderr.decode(errors='replace').strip()}")
    return r.stdout


def _baseline_ls(path):
    out = _git("ls-tree", "-r", "--name-only", BASELINE, "--", path).decode()
    return sorted(p for p in out.splitlines() if p)


def _baseline_contents(paths):
    """Bytes of each path as represented in BASELINE (git cat-file --batch). Missing paths -> None."""
    req = "".join(f"{BASELINE}:{p}\n" for p in paths).encode()
    raw = _git("cat-file", "--batch", input_bytes=req)
    out, pos = {}, 0
    for p in paths:
        nl = raw.index(b"\n", pos)
        header = raw[pos:nl].decode()
        pos = nl + 1
        if header.endswith(" missing"):
            out[p] = None
            continue
        size = int(header.split()[2])
        out[p] = raw[pos:pos + size]
        pos += size + 1                                              # content + trailing LF
    return out


def _sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def _protocol_text():
    data = _baseline_contents([PROTOCOL])[PROTOCOL]
    assert data is not None, f"{PROTOCOL} is not in baseline {BASELINE}"
    return data.decode("utf-8")


def _section(text, start, end):
    i = text.index(start)
    return text[i:text.index(end, i)]


_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def _num(s):
    """Protocol number -> canonical string: unicode minus -> '-', 'a×10⁻ᵇ' -> 'a e-0b' (Python %.Ne style)."""
    s = s.strip().replace("−", "-")
    m = re.fullmatch(r"([+-]?\d+(?:\.\d+)?)×10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻]+)", s)
    if m:
        return f"{m.group(1)}e{int(m.group(2).translate(_SUP)):+03d}"
    return s


def _rows(section, first_cell_regex):
    rows = []
    for line in section.splitlines():
        if not line.startswith("|"):
            continue
        cells = [c.strip().replace("**", "").strip() for c in line.strip().strip("|").split("|")]
        if re.fullmatch(first_cell_regex, cells[0]):
            rows.append(cells)
    return rows


# ----------------------------------------------------------------------
# T1.1 — §1.2 predeclared energy points (exact reproduction at printed precision)
# ----------------------------------------------------------------------
def _point_rows():
    sec = _section(_protocol_text(), "### 1.2 Predeclared energy points", "### 1.3")
    rows = _rows(sec, r"(L[123]|T(3|10|11)[−+]|U[123]|C[12])")
    assert len(rows) == 14, f"expected 14 point rows in protocol §1.2, found {len(rows)}"
    return rows


def test_section_1_2_points_reproduced_exactly():
    k = b2.kx_window(20.0, 0.2)                                     # coarse B4 grid (protocol §1.1/§1.2)
    Ep = b2.e_perp()
    sk = 1.0 / (np.sqrt(2.0) * 0.2)
    assert k.size == 1274
    for cells in _point_rows():
        pid, gi, kx_s, E_s, oc, thr, dE, z = cells[:8]
        if gi.startswith("—"):                                      # controls: exact B3 energies kx0 + j sigma_k
            j = {"C1": -1, "C2": 2}[pid]
            kx = 20.0 + j * sk
        else:
            kx = float(k[int(gi)])
        E = 0.5 * (kx ** 2 + np.pi ** 2)
        d = np.abs(E - Ep)
        n = int(np.argmin(d))
        n_open = int(np.sum(Ep < E))
        assert f"{kx:.6f}" == kx_s, (pid, kx, kx_s)
        assert f"{E:.5f}" == E_s, (pid, E, E_s)
        assert oc == f"{n_open} ({63 - n_open})", (pid, oc)
        assert thr == f"{n + 1} ({Ep[n]:.4f})", (pid, thr)
        dE_s, rel_s = re.fullmatch(r"(\S+) \((\S+)\)", dE).groups()
        assert f"{E - Ep[n]:+.3f}" == _num(dE_s), (pid, E - Ep[n], dE_s)
        assert f"{d[n] / E:.1e}" == _num(rel_s), (pid, d[n] / E, rel_s)
        assert f"{(kx - 20.0) / sk:+.2f}" == _num(z), (pid, z)
        a = b2.adjust_energy(E, Ep)
        assert a["adjusted"] is False, f"{pid} would be moved by adjust_energy (protocol: none of the 14 is moved)"


# ----------------------------------------------------------------------
# T1.2 — §4.1 threshold-moved nodes
# ----------------------------------------------------------------------
def _moved(dkx):
    k = b2.kx_window(20.0, 0.2, dkx)
    h = k[1] - k[0]
    Ep = b2.e_perp()
    out = []
    for i, kx in enumerate(k):
        E = 0.5 * (kx ** 2 + np.pi ** 2)
        a = b2.adjust_energy(E, Ep)
        if a["adjusted"]:
            km = float(np.sqrt(2.0 * a["E"] - np.pi ** 2))
            w = h * (0.5 if i in (0, k.size - 1) else 1.0) * float(b2.packet_weights(kx, 20.0, 0.2))
            out.append(dict(i=i, kx=float(kx), km=km, d_h=(km - kx) / h, E0=a["E_original"], E1=a["E"],
                            n=a["nearest_threshold_n"], side="above" if E > a["nearest_threshold_E"] else "below",
                            kept=a["n_open"] == int(np.sum(Ep < E)), w=w))
    return out


def test_section_4_1_coarse_moved_nodes_reproduced_exactly():
    sec = _section(_protocol_text(), "### 4.1 Design-time enumeration", "### 4.2")
    assert "11 nodes are moved" in sec
    rows = _rows(sec, r"\d+")
    rows = [r for r in rows if len(r) == 9]                         # the coarse table (9 columns)
    got = _moved(b2.DKX)
    assert len(got) == 11 and len(rows) == 11
    for r, g in zip(rows, got):
        i, kx, km, dh, ee, n, side, kept, w = r
        assert int(i) == g["i"]
        assert f"{g['kx']:.6f}" == kx and f"{g['km']:.6f}" == km, (i, g)
        assert f"{g['d_h']:+.3f}" == _num(dh), (i, g["d_h"], dh)
        e0, e1 = [x.strip() for x in ee.split("→")]
        assert f"{g['E0']:.4f}" == e0 and f"{g['E1']:.4f}" == e1, (i, ee)
        assert int(n) == g["n"] and side == g["side"] and kept == "yes" and g["kept"], (i, side, kept)
        assert f"{g['w']:.1e}" == _num(w), (i, g["w"], w)
    total = sum(g["w"] for g in got)
    assert f"{total:.1e}" == "8.1e-03" and "8.1×10⁻³ of the packet weight" in sec


def test_section_4_1_fine_moved_nodes_and_shared_energies():
    sec = _section(_protocol_text(), "### 4.1 Design-time enumeration", "### 4.2")
    assert "18 nodes are moved" in sec and "**2.44 grid steps**" in sec
    got = _moved(b2.DKX / 2)
    assert len(got) == 18
    assert all(g["kept"] for g in got), "a fine-grid moved node changed side"
    assert f"{max(abs(g['d_h']) for g in got):.2f}" == "2.44"
    by_E = {}
    for g in got:
        by_E.setdefault(round(g["E1"], 10), []).append(g["i"])
    shared = {tuple(v): e for e, v in by_E.items() if len(v) > 1}
    expected = {(int(a), int(b_)): e for a, b_, e in re.findall(r"\| (\d+), (\d+) \| E' = ([\d.]+) \|", sec)}
    assert len(expected) == 3 and set(shared) == set(expected), (shared, expected)
    for pair, e in expected.items():
        assert f"{shared[pair]:.5f}" == e, (pair, shared[pair], e)
    # every coarse node is a fine node, and the coarse grid is the even subset (halving grid nested)
    k, kf = b2.kx_window(20.0, 0.2), b2.kx_window(20.0, 0.2, b2.DKX / 2)
    assert kf.size == 2547 and np.allclose(kf[::2], k, rtol=0, atol=1e-12)


# ----------------------------------------------------------------------
# T1.3 — §2.4 frame-time equality (verified production formula; no solver run)
# ----------------------------------------------------------------------
def _frame_times(t_end, dt, n_frames):
    """Exact replica of SplitStepDST.run's step/frame bookkeeping (guarded by the source check below)."""
    n_steps = max(1, int(np.ceil(t_end / dt - 1e-9)))
    dt = t_end / n_steps
    record = sorted(set(np.unique(np.round(np.linspace(0, n_steps, n_frames)).astype(int)).tolist()))
    return np.array([s * dt for s in record]), record


def test_run_uses_the_replicated_frame_formula():
    src = inspect.getsource(SplitStepDST.run)
    for line in ("n_steps = max(1, int(np.ceil(t_end / dt - 1e-9)))", "dt = t_end / n_steps",
                 "record = set(np.unique(np.round(np.linspace(0, n_steps, n_frames)).astype(int)).tolist())",
                 'out["times"].append(step * dt)'):
        assert line in src, f"SplitStepDST.run no longer contains: {line}"


def test_section_2_4_frame_times_bitwise_equal():
    t_max, dtf = 6.0, 0.005
    n_frames = int(round(t_max / dtf)) + 1                          # 1201 frames, t_j = j*0.005 (protocol §2.1/§2.4)
    prod, rp = _frame_times(t_max, 2e-4, n_frames)                  # production and TD-R2 (dt unchanged)
    r1, r1s = _frame_times(t_max, 1e-4, n_frames)                   # TD-R1
    assert n_frames == 1201 and len(prod) == len(r1) == 1201
    assert rp == list(range(0, 30001, 25)) and r1s == list(range(0, 60001, 50))
    assert np.array_equal(prod, r1), "production and TD-R1 frame times are not bitwise equal"
    assert np.allclose(prod, np.arange(1201) * dtf, rtol=0, atol=1e-12)


# ----------------------------------------------------------------------
# T1.4 — §3.3.1 dx-equality constants (traceable to the protocol table)
# ----------------------------------------------------------------------
DX_PROD = 3.9 / 250


def _variant_rows():
    sec = _section(_protocol_text(), "#### 3.3.1 Variants and the shifted frame", "#### 3.3.2")
    rows = _rows(sec, r"V[0-4] \(.*\)")
    assert len(rows) == 5
    return sec, [(r[0], int(r[1]), r[2], r[6], int(r[7])) for r in rows]


def test_section_3_3_1_variant_table_constants():
    _, rows = _variant_rows()
    for name, m, s_txt, lx_txt, nx in rows:
        s = m * DX_PROD
        Lx = 3.9 + 2 * m * DX_PROD                                  # the construction of §3.3.1 (all variants)
        assert (m == 0 and s_txt == "0") or f"{s:.4f}" == s_txt, (name, s, s_txt)
        assert f"{Lx:.4f}".rstrip("0").rstrip(".") == lx_txt.rstrip("0").rstrip("."), (name, Lx, lx_txt)
        assert nx == 250 + 2 * m - 1, (name, nx)


@pytest.mark.parametrize("m,expected", [(39, True), (64, True), (115, True), (38, False)])
def test_section_3_3_1_dx_bitwise_equality(m, expected):
    """L_x is constructed as 3.9 + 2*m*dx (equivalently (nx+1)*dx); build_grid's dx = L_x/(nx+1) must equal the
    production dx bitwise for m = 39, 64, 115, and must not for m = 38 (protocol §3.3.1)."""
    sec, _ = _variant_rows()
    assert "For m = 38 it does not" in sec
    nx = 250 + 2 * m - 1
    Lx = 3.9 + 2 * m * DX_PROD
    cfg = waveguide_config().replace(Lx=Lx, nx=nx)
    _, _, dx, _ = build_grid(cfg)
    assert (dx == DX_PROD) is expected, (m, repr(dx), repr(DX_PROD))
    assert ((((nx + 1) * DX_PROD) / (nx + 1)) == DX_PROD) is expected     # equivalent construction (nx+1)*dx


def test_production_grid_dx_is_the_reference_dx():
    cfg = waveguide_config()
    _, _, dx, dy = build_grid(cfg)
    assert cfg.nx == 249 and cfg.ny == 63 and dx == DX_PROD and dy == 1.0 / 64


def test_printed_decimal_lx_is_not_a_safe_construction_for_m115():
    """Disclosure (not a protocol change): the table prints L_x = 7.4880 for V4, but the decimal literal 7.488
    gives L_x/480 != dx; only the formula construction 3.9 + 2*m*dx is bitwise-safe. Implementations must use it."""
    assert 7.488 / 480 != DX_PROD
    assert (3.9 + 2 * 115 * DX_PROD) / 480 == DX_PROD


# ----------------------------------------------------------------------
# T1.5 — protected artifacts against the A6 baseline (SHA-256 of file bytes)
# ----------------------------------------------------------------------
def _protected_paths():
    paths = list(PROTECTED_SOURCES) + list(PROTECTED_DOCS)
    baseline_tests = _baseline_ls("tests")                         # every test file tracked at the A6 baseline
    assert baseline_tests and "tests/test_b4_pilot.py" not in baseline_tests
    paths += baseline_tests
    for d in PROTECTED_RESULT_DIRS:
        tracked = _baseline_ls(d)
        assert tracked, f"no tracked files under {d} in baseline {BASELINE}; refusing to omit it silently"
        paths += tracked
    return paths


def test_baseline_commit_exists_and_is_an_ancestor_of_head():
    assert _git("cat-file", "-t", BASELINE).decode().strip() == "commit"
    r = subprocess.run(["git", "-C", REPO, "merge-base", "--is-ancestor", BASELINE, "HEAD"], capture_output=True)
    assert r.returncode == 0, f"baseline {BASELINE} is not an ancestor of HEAD"


def test_protected_files_unchanged_against_baseline():
    paths = _protected_paths()
    base = _baseline_contents(paths)
    missing_in_base = [p for p, v in base.items() if v is None]
    assert not missing_in_base, f"protected paths absent from baseline {BASELINE}: {missing_in_base}"
    problems = []
    for p in paths:
        full = os.path.join(REPO, p)
        if not os.path.isfile(full):
            problems.append(f"MISSING in working tree: {p}")
            continue
        with open(full, "rb") as f:
            now = _sha256(f.read())
        if now != _sha256(base[p]):
            problems.append(f"CHANGED: {p} (baseline {_sha256(base[p])}, now {now})")
    assert not problems, "\n".join(problems)
    n_results = sum(1 for p in paths if p.startswith("results/"))
    assert n_results == 109, f"protected result-file count changed: {n_results}"   # b1 8, b2 10, b3 13, b3x 8, b3y 11, validation 7, transport 52


def test_protected_result_dirs_contain_no_extra_files():
    for d in PROTECTED_RESULT_DIRS:
        tracked = set(_baseline_ls(d))
        present = set()
        for root, _, files in os.walk(os.path.join(REPO, d)):
            for f in files:
                present.add(os.path.relpath(os.path.join(root, f), REPO))
        extra = sorted(present - tracked)
        assert not extra, f"files not in baseline {BASELINE} under protected {d}: {extra}"


# ----------------------------------------------------------------------
# T2.1 — L1 disorder-parameter identity
# ----------------------------------------------------------------------
def disorder_params_hash(d):
    """Pilot L1 recipe: SHA-256 over the bytes of concat(qx, qy, phi) as contiguous little-endian IEEE-754 float64
    ('<f8'), prefixed 'sha256:'. On a little-endian platform this equals the B3/B3X record recipe
    (validate_rgf_b3.sha(np.concatenate([qx, qy, phi])), native float64)."""
    a = np.ascontiguousarray(np.concatenate([d.qx, d.qy, d.phi]), dtype="<f8")
    return "sha256:" + hashlib.sha256(a.tobytes()).hexdigest()


def _fields(d):
    return dict(W=d.strength, zeta=d.corr_length, seed=d.seed, M=d.n_modes, edge=d.edge,
                region=None if d.region is None else list(d.region))


def _td_production_cfg(case):
    data = np.load(os.path.join(REPO, TD_FILES[case]), allow_pickle=False)
    seeds = list(data["seeds"])
    assert 0 in seeds, f"{TD_FILES[case]} does not contain seed 0"
    return Config.from_dict(json.loads(str(data["config_json"])))


def test_L1_recipe_equals_b3_record_recipe():
    d = b3.case_potential("A")[1].disorder
    assert sys.byteorder == "little"
    assert disorder_params_hash(d) == b3.sha(np.concatenate([d.qx, d.qy, d.phi]))


@pytest.mark.parametrize("case", CASES)
def test_L1_parameter_identity_rgf_td_and_b3_records(case):
    rgf = b3.case_potential(case)[1].disorder                      # RGF path (B3 construction)
    td = disorder_from_config(_td_production_cfg(case), seed=0)    # TD production path (committed config)
    h_rgf, h_td = disorder_params_hash(rgf), disorder_params_hash(td)
    assert h_rgf == h_td, (case, h_rgf, h_td)
    assert _fields(rgf) == _fields(td), (case, _fields(rgf), _fields(td))
    p = b3.CASES[case]
    assert (rgf.strength, rgf.corr_length, rgf.seed) == (p["W"], p["zeta"], p["seed"])

    rec3 = json.load(open(os.path.join(REPO, "results/rgf_b3_validation/run_metadata.json")))["cases"][case]
    rec3x = json.load(open(os.path.join(REPO, "results/rgf_b3x_validation/run_metadata.json")))["configurations"][
        B3X_NAMES[case]]
    for name, rec in (("B3", rec3), ("B3X", rec3x)):
        assert rec["disorder_params_hash"] == h_rgf, (case, name, rec["disorder_params_hash"], h_rgf)
        assert (rec["W"], rec["zeta"], rec["seed"], rec["disorder_n_modes"], rec["disorder_edge"],
                rec["disorder_region"]) == (rgf.strength, rgf.corr_length, rgf.seed, rgf.n_modes, rgf.edge,
                                            list(rgf.region)), (case, name)


def test_L1_b3y_records_carry_no_disorder_hash():
    """Availability check, not an equality claim: B3Y metadata has no disorder_params_hash, so L1 compares against
    B3 and B3X records only."""
    text = open(os.path.join(REPO, "results/rgf_b3y_validation/run_metadata.json")).read()
    assert "disorder_params_hash" not in text


def test_L1_detects_a_different_realization():
    a = disorder_params_hash(b3.case_potential("A")[1].disorder)
    cfg = b3.case_config("A")
    assert disorder_params_hash(disorder_from_config(cfg, seed=1)) != a
    assert disorder_params_hash(disorder_from_config(cfg.replace(disorder_corr_length=0.021), seed=0)) != a


# ----------------------------------------------------------------------
# T2.2 — L2 sampler fidelity on the production TD grid (bitwise)
# ----------------------------------------------------------------------
def _first_difference(a, b):
    if a.shape != b.shape or a.dtype != b.dtype:
        return f"shape/dtype differ: {a.shape} {a.dtype} vs {b.shape} {b.dtype}"
    bad = np.argwhere(~((a == b) | (np.isnan(a) & np.isnan(b))))
    if bad.size == 0:
        return None
    idx = tuple(int(i) for i in bad[0])
    return f"{bad.shape[0]} differing entries; first at {idx}: {a[idx]!r} vs {b[idx]!r}"


@pytest.mark.parametrize("case", CASES)
def test_L2_production_td_grid_samples_bitwise(case):
    cfg = _td_production_cfg(case)
    assert (cfg.nx, cfg.ny, cfg.Lx, cfg.Ly) == (249, 63, 3.9, 1.0)
    solver = SplitStepDST(cfg, [Potential(cfg, disorder_from_config(cfg, seed=0))])
    # independent construction of the production coordinates: x_i = (i+1) dx, dx = Lx/(nx+1); same for y
    dx, dy = cfg.Lx / (cfg.nx + 1), cfg.Ly / (cfg.ny + 1)
    X, Y = np.meshgrid(np.arange(1, cfg.nx + 1) * dx, np.arange(1, cfg.ny + 1) * dy, indexing="ij")
    msg = _first_difference(solver.X, X) or _first_difference(solver.Y, Y)
    assert msg is None, f"case {case}: solver grid differs from the production coordinates: {msg}"
    direct = np.asarray(Potential(cfg, disorder_from_config(cfg, seed=0)).numpy(X, Y), dtype=np.float64)
    stored = solver.V[0]
    assert solver.V.shape == (1, cfg.nx, cfg.ny) and stored.dtype == np.float64
    msg = _first_difference(stored, direct)
    assert msg is None, f"case {case}: SplitStepDST.V is not bitwise equal to a direct Potential.numpy call: {msg}"
    assert np.all(np.isfinite(stored)) and float(np.max(np.abs(stored))) > 0.0


def test_L2_comparison_detects_a_one_ulp_change():
    a = np.linspace(-1.0, 1.0, 12).reshape(3, 4)
    b = a.copy()
    b[1, 2] = np.nextafter(b[1, 2], np.inf)
    msg = _first_difference(a, b)
    assert msg is not None and "(1, 2)" in msg


def test_module_is_labelled_diagnostic():
    assert CLASSIFICATION == "DIAGNOSTIC / NOT B4" and CLASSIFICATION in (__doc__ or "")


# ======================================================================
# T3 — RGF observed-order classification (§1.4) and deterministic resolution selection (§1.6)
# Synthetic, hand-constructed inputs only; no RGF evaluation. Passing these tests does NOT mean G-RGF passed.
# Dyadic inputs (Q(1/512) = 0 or 1, changes = powers of two) make every subtraction exact, so the expected Δ, p and
# Richardson values below are computed by hand from the literals, independently of the implementation.
# ======================================================================
import math  # noqa: E402

import b4_pilot_rgf as rgf  # noqa: E402


def _tri(d1, d2, q0=0.0):
    """(Q(1/512), Q(1/1024), Q(1/2048)) with exact consecutive changes d1, d2 (dyadic inputs)."""
    q1 = q0 + d1
    q2 = q1 + d2
    assert q1 - q0 == d1 and q2 - q1 == d2, "synthetic input is not exact"
    return (q0, q1, q2)


def _mirror(t):
    return tuple(1.0 - v for v in t)


# p = 2 cases (|Δ₁| = 4 |Δ₂|); hand values: e(1/2048) = |Δ₂|/3, e(1/1024) = 4|Δ₂|/3, e(1/512) = 16|Δ₂|/3
CASE_ALL = _tri(2.0 ** -15, 2.0 ** -17)        # |Δ₂| = 7.63e-6:  e = 4.07e-5 / 1.02e-5 / 2.54e-6 -> all admissible
CASE_1024 = _tri(2.0 ** -12, 2.0 ** -14)       # |Δ₂| = 6.10e-5:  e(512) = 3.26e-4 > 1e-4; e(1024) = 8.14e-5; e(2048) = 2.03e-5
CASE_2048 = _tri(2.0 ** -11, 2.0 ** -13)       # |Δ₂| = 1.22e-4:  e(1024) = 1.63e-4 > 1e-4; e(2048) = 4.07e-5, |Δ₂| <= 2e-4
CASE_NONE = _tri(2.0 ** -10, 2.0 ** -12)       # |Δ₂| = 2.44e-4 > 2e-4 (e(2048) = 8.14e-5 <= 1e-4): nothing admissible


def test_T3_protocol_values_and_label():
    assert rgf.CLASSIFICATION == "DIAGNOSTIC / NOT B4" and "DIAGNOSTIC / NOT B4" in rgf.__doc__
    assert rgf.ORDER_FLOOR == 1e-8 and rgf.ORDER_WINDOW == (1.8, 2.2)
    assert rgf.CHANGE_LIMIT == 2e-4 and rgf.ERROR_LIMIT == 1e-4
    assert rgf.DX_CANDIDATES == ("1/512", "1/1024", "1/2048") and rgf.RICHARDSON_K == {"1/512": 2, "1/1024": 1,
                                                                                         "1/2048": 0}
    sec = _section(_protocol_text(), "### 1.6 Deterministic resolution decision", "## 2. Part 2")
    assert "k = 2, 1, 0 for dx_c = 1/512, 1/1024, 1/2048" in sec and "**coarsest** candidate" in sec


def test_T3_both_changes_below_floor():
    r = rgf.classify_order(0.5, 0.5 + 3e-9, 0.5 + 3e-9 + 2e-9)
    assert r.status == "BELOW FLOOR" and r.p is None


def test_T3_floor_is_strict_exactly_one_below():
    # |Δ₁| = 1e-8 exactly (Q(1/512) = 0) is NOT below the floor; |Δ₂| = 1e-9 is -> exactly one below
    r = rgf.classify_order(0.0, 1e-8, 1e-8 + 1e-9)
    assert r.delta1 == 1e-8 and r.status == "NOT RESOLVABLE" and r.p is None
    r = rgf.classify_order(0.0, 1e-6, 1e-6 + 1e-9)
    assert r.status == "NOT RESOLVABLE" and r.p is None


def test_T3_opposite_signs_not_resolvable():
    r = rgf.classify_order(0.0, 2.0 ** -15, 2.0 ** -15 - 2.0 ** -18)
    assert r.delta1 > 0 > r.delta2 and r.status == "NOT RESOLVABLE" and r.p is None


def test_T3_p_below_window():
    r = rgf.classify_order(*_tri(2.0 ** -16, 2.0 ** -17))      # p = 1 exactly
    assert r.p == 1.0 and r.status == "OUTSIDE WINDOW"


def test_T3_p_exactly_1_8_is_resolved():
    q = (0.25, 0.26360235255150194, 0.26750860255150194)       # exact: Δ₁ = 0.01360235255150194, Δ₂ = 2^-8
    assert q[1] - q[0] == 0.01360235255150194 and q[2] - q[1] == 2.0 ** -8
    assert math.log2(0.01360235255150194 / 2.0 ** -8) == 1.8
    r = rgf.classify_order(*q)
    assert r.p == 1.8 and r.status == "RESOLVED"


def test_T3_p_inside_window():
    r = rgf.classify_order(*CASE_ALL)
    assert r.p == 2.0 and r.status == "RESOLVED"


def test_T3_p_exactly_2_2_is_resolved():
    q = (0.0, 6.846775740367383e-08, 8.336891859752149e-08)    # exact: Δ₂ = 2^-26 = 1.4901161193847656e-08
    assert q[2] - q[1] == 2.0 ** -26 and math.log2((q[1] - q[0]) / 2.0 ** -26) == 2.2
    r = rgf.classify_order(*q)
    assert r.p == 2.2 and r.status == "RESOLVED"


def test_T3_p_above_window():
    r = rgf.classify_order(*_tri(2.0 ** -14, 2.0 ** -17))      # p = 3 exactly
    assert r.p == 3.0 and r.status == "OUTSIDE WINDOW"


def test_T3_negative_same_sign_changes_classified_like_positive():
    r = rgf.classify_order(*_mirror(CASE_ALL))
    assert r.delta1 < 0 and r.delta2 < 0 and r.p == 2.0 and r.status == "RESOLVED"


def test_T3_T1_and_R1_classified_independently():
    t = rgf.classify_order(*CASE_ALL)
    r = rgf.classify_order(1.0, 1.0 - 2.0 ** -20, 1.0 - 2.0 ** -20 + 2.0 ** -21)   # opposite signs
    assert t.status == "RESOLVED" and r.status == "NOT RESOLVABLE"


def test_T3_non_numeric_or_non_finite_input_rejected():
    with pytest.raises(ValueError):
        rgf.classify_order(0.0, float("nan"), 1.0)
    with pytest.raises(TypeError):
        rgf.classify_order(0.0, None, 1.0)


def test_T3_richardson_estimates_by_hand():
    o = rgf.classify_order(*CASE_ALL)
    d2 = 2.0 ** -17
    assert rgf.richardson_error(o, "1/2048") == pytest.approx(d2 / 3, rel=1e-15)
    assert rgf.richardson_error(o, "1/1024") == pytest.approx(4 * d2 / 3, rel=1e-15)
    assert rgf.richardson_error(o, "1/512") == pytest.approx(16 * d2 / 3, rel=1e-15)
    # k = 2 identity of §1.6: e(1/512) = 2^p |Δ₁| / (2^p - 1)
    assert rgf.richardson_error(o, "1/512") == pytest.approx(4 * 2.0 ** -15 / 3, rel=1e-15)


@pytest.mark.parametrize("q", [(0.5, 0.5 + 3e-9, 0.5 + 5e-9), (0.0, 1e-6, 1e-6 + 1e-9), _tri(2.0 ** -16, 2.0 ** -17)])
def test_T3_no_richardson_estimate_for_non_resolved(q):
    o = rgf.classify_order(*q)
    assert o.status != "RESOLVED"
    with pytest.raises(ValueError):
        rgf.richardson_error(o, "1/1024")
    assert all(c.error_estimate is None for c in rgf.admissibility(o).values())


def _labels(q):
    a = rgf.admissibility(rgf.classify_order(*q))
    return tuple(a[dx].label for dx in ("1/512", "1/1024", "1/2048"))


def test_T3_admissibility_rule_1_512():
    assert _labels(CASE_ALL)[0] == "ADMISSIBLE"              # |Δ₁| = 3.05e-5, e = 4.07e-5
    assert _labels(CASE_1024)[0] == "INADMISSIBLE"           # e(1/512) = 3.26e-4 > 1e-4


def test_T3_admissibility_rule_1_1024():
    assert _labels(CASE_1024)[1] == "ADMISSIBLE"             # |Δ₂| = 6.10e-5, e = 8.14e-5
    assert _labels(CASE_2048)[1] == "INADMISSIBLE"           # e(1/1024) = 1.63e-4 > 1e-4


def test_T3_admissibility_rule_1_2048():
    assert _labels(CASE_2048)[2] == "ADMISSIBLE-R"           # |Δ₂| = 1.22e-4 <= 2e-4, e = 4.07e-5
    assert _labels(CASE_NONE)[2] == "INADMISSIBLE"           # |Δ₂| = 2.44e-4 > 2e-4 although e = 8.14e-5 <= 1e-4


def test_T3_below_floor_and_unclassifiable_labels():
    assert _labels((0.5, 0.5 + 3e-9, 0.5 + 5e-9)) == ("ADMISSIBLE-F",) * 3
    assert _labels((0.0, 1e-6, 1e-6 + 1e-9)) == ("INADMISSIBLE-U",) * 3
    assert _labels(_tri(2.0 ** -16, 2.0 ** -17)) == ("INADMISSIBLE-U",) * 3        # OUTSIDE WINDOW


def test_T3_admissible_F_is_never_an_error_bound():
    meaning = rgf.LABEL_MEANING["ADMISSIBLE-F"]
    assert "not an error bound" in meaning.lower()
    a = rgf.admissibility(rgf.classify_order(0.5, 0.5 + 3e-9, 0.5 + 5e-9))
    assert all(c.label == "ADMISSIBLE-F" and c.error_estimate is None for c in a.values())


@pytest.mark.parametrize("q", [CASE_ALL, CASE_1024, CASE_2048, CASE_NONE, _mirror(CASE_2048)])
def test_T3_admissible_R_only_at_1_2048(q):
    lab = _labels(q)
    assert "ADMISSIBLE-R" not in lab[:2] and lab[2] != "ADMISSIBLE"


def _vals(spec):
    """spec: {point: (T1 triple, R1 triple)} -> values mapping."""
    out = {}
    for x, (t, r) in spec.items():
        out[(x, "T1")] = t
        out[(x, "R1")] = r
    return out


def test_T3_selection_coarsest_admissible():
    spec = {"x1": (CASE_ALL, _mirror(CASE_ALL)), "x2": (CASE_1024, _mirror(CASE_1024))}
    s = rgf.select_resolution(_vals(spec), ["x1", "x2"], solver_agreement=True)
    assert s.outcome == "SELECTED" and s.dx_selected == "1/1024"   # 1/512 fails at x2; 1/1024 admissible everywhere
    assert set(s.labels_at_selected.values()) == {"ADMISSIBLE"} and len(s.labels_at_selected) == 4
    assert s.statements == ()


def test_T3_selection_finest_with_richardson_only_statement():
    spec = {"x1": (CASE_ALL, _mirror(CASE_ALL)), "x2": (CASE_2048, _mirror(CASE_2048))}
    s = rgf.select_resolution(_vals(spec), ["x1", "x2"], solver_agreement=True)
    assert s.outcome == "SELECTED" and s.dx_selected == "1/2048"
    assert set(s.labels_at_selected.values()) == {"ADMISSIBLE-R"}
    assert any("no next-finer change" in t for t in s.statements)


def test_T3_selection_preserves_mixed_labels_and_F_statement():
    below = (0.5, 0.5 + 3e-9, 0.5 + 5e-9)
    spec = {"x1": (CASE_ALL, _mirror(CASE_ALL)), "x2": (below, _mirror(below))}
    s = rgf.select_resolution(_vals(spec), ["x1", "x2"], solver_agreement=True)
    assert s.outcome == "SELECTED" and s.dx_selected == "1/512"
    assert s.labels_at_selected[("x1", "T1")] == "ADMISSIBLE" and s.labels_at_selected[("x2", "R1")] == "ADMISSIBLE-F"
    assert any("NOT an error bound" in t for t in s.statements)


def test_T3_selection_requires_both_quantities_admissible():
    # T1 admissible everywhere at 1/512, but R1 at x1 only at 1/1024 -> 1/1024
    spec = {"x1": (CASE_ALL, _mirror(CASE_1024)), "x2": (CASE_ALL, _mirror(CASE_ALL))}
    s = rgf.select_resolution(_vals(spec), ["x1", "x2"], solver_agreement=True)
    assert s.dx_selected == "1/1024"


def test_T3_not_evaluable_when_any_unclassifiable():
    nr = (1.0, 1.0 - 2.0 ** -20, 1.0 - 2.0 ** -20 + 2.0 ** -21)   # opposite signs
    spec = {"x1": (CASE_ALL, _mirror(CASE_ALL)), "x2": (CASE_ALL, nr)}
    s = rgf.select_resolution(_vals(spec), ["x1", "x2"], solver_agreement=True)
    assert s.outcome == "G-RGF NOT EVALUABLE" and s.dx_selected is None and s.labels_at_selected == {}
    assert s.status_mismatch_points == ("x2",)


def test_T3_fail_when_all_classified_but_limits_exceeded():
    spec = {"x1": (CASE_ALL, _mirror(CASE_ALL)), "x2": (CASE_NONE, _mirror(CASE_NONE))}
    s = rgf.select_resolution(_vals(spec), ["x1", "x2"], solver_agreement=True)
    assert s.outcome == "G-RGF FAIL" and s.dx_selected is None


def test_T3_not_evaluable_takes_precedence_over_fail():
    spec = {"x1": (CASE_NONE, _mirror(CASE_NONE)), "x2": (_tri(2.0 ** -14, 2.0 ** -17), _mirror(CASE_ALL))}
    s = rgf.select_resolution(_vals(spec), ["x1", "x2"], solver_agreement=True)
    assert s.outcome == "G-RGF NOT EVALUABLE"


def test_T3_step0_solver_agreement_is_mandatory():
    v = _vals({"x1": (CASE_ALL, _mirror(CASE_ALL))})
    assert rgf.select_resolution(v, ["x1"], solver_agreement=False).outcome == "STOP: SOLVER AGREEMENT NOT ESTABLISHED"
    with pytest.raises(TypeError):
        rgf.select_resolution(v, ["x1"], solver_agreement=None)
    with pytest.raises(TypeError):
        rgf.select_resolution(v, ["x1"])                                        # keyword is required


def test_T3_missing_quantity_point_or_extra_input_rejected():
    v = _vals({"x1": (CASE_ALL, _mirror(CASE_ALL)), "x2": (CASE_ALL, _mirror(CASE_ALL))})
    no_r = {k: q for k, q in v.items() if k != ("x2", "R1")}
    with pytest.raises(ValueError, match="missing"):
        rgf.select_resolution(no_r, ["x1", "x2"], solver_agreement=True)
    with pytest.raises(ValueError, match="missing"):
        rgf.select_resolution(v, ["x1", "x2", "x3"], solver_agreement=True)    # required point without data
    with pytest.raises(ValueError, match="unexpected"):
        rgf.select_resolution(v, ["x1"], solver_agreement=True)                # data for a non-required point
    with pytest.raises(ValueError):
        rgf.select_resolution(v, [], solver_agreement=True)
    with pytest.raises(ValueError):
        rgf.select_resolution(v, ["x1", "x1", "x2"], solver_agreement=True)
    bad = dict(v)
    bad[("x1", "T1")] = (0.0, 1.0)
    with pytest.raises(ValueError):
        rgf.select_resolution(bad, ["x1", "x2"], solver_agreement=True)


def test_T3_selection_is_deterministic():
    spec = {"x1": (CASE_ALL, _mirror(CASE_ALL)), "x2": (CASE_1024, _mirror(CASE_1024)),
            "x3": ((0.5, 0.5 + 3e-9, 0.5 + 5e-9), _mirror((0.5, 0.5 + 3e-9, 0.5 + 5e-9)))}
    a = rgf.select_resolution(_vals(spec), ["x1", "x2", "x3"], solver_agreement=True)
    b = rgf.select_resolution(dict(reversed(list(_vals(spec).items()))), ["x1", "x2", "x3"], solver_agreement=True)
    assert a == b and a.dx_selected == "1/1024" and a.outcome == "SELECTED"


def test_T3_module_is_pure():
    import ast
    tree = ast.parse(open(os.path.join(REPO, "b4_pilot_rgf.py"), encoding="utf-8").read())
    names = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, (ast.Import, ast.ImportFrom))
             for a in n.names} | {n.module.split(".")[0] for n in ast.walk(tree)
                                  if isinstance(n, ast.ImportFrom) and n.module}
    assert names <= {"__future__", "math", "dataclasses", "typing", "annotations", "dataclass", "field", "Mapping",
                     "Sequence"}, names
    calls = {n.func.id for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)}
    assert "open" not in calls


# ----------------------------------------------------------------------
# T3.1 — §1.4 point flag: status mismatch OR |ΔT1 + ΔR1| > recorded unitarity residuals
# Bound (derived in b4_pilot_rgf): |Δ₁T1 + Δ₁R1| <= U(1/512) + U(1/1024), |Δ₂T1 + Δ₂R1| <= U(1/1024) + U(1/2048),
# because T1 + R1 - 1 = (S†S - I)_11 and |(S†S - I)_11| <= U. Dyadic inputs; expected sums and bounds by hand.
# ----------------------------------------------------------------------
_T_OK = CASE_ALL                                                           # (0, 2^-15, 2^-15 + 2^-17), RESOLVED p = 2
_R_OK = _mirror(CASE_ALL)                                                  # ΔR = -ΔT exactly -> sums 0
_R_SHIFTED = (1.0, 1.0 - 2.0 ** -15 + 2.0 ** -25, 1.0 - 2.0 ** -15 + 2.0 ** -25 - 2.0 ** -17 + 2.0 ** -27)
#   ΔR = (-(2^-15 - 2^-25), -(2^-17 - 2^-27)): ratio 4 -> RESOLVED p = 2; sums = (2^-25, 2^-27)
_R_MISMATCH = (1.0, 1.0 - 2.0 ** -15, 1.0 - 2.0 ** -15 + 2.0 ** -17)     # opposite signs -> NOT RESOLVABLE; sums (0, 2^-16)


def test_T3_1_inputs_are_exact():
    assert _R_SHIFTED[1] - _R_SHIFTED[0] == -(2.0 ** -15 - 2.0 ** -25)
    assert _R_SHIFTED[2] - _R_SHIFTED[1] == -(2.0 ** -17 - 2.0 ** -27)
    assert _R_MISMATCH[2] - _R_MISMATCH[1] == 2.0 ** -17


def test_T3_1_neither_condition():
    f = rgf.point_flags({("x", "T1"): _T_OK, ("x", "R1"): _R_OK}, {"x": (2.0 ** -40,) * 3}, ["x"])["x"]
    assert f.unitarity.sum_delta1 == 0.0 and f.unitarity.sum_delta2 == 0.0
    assert not f.status_mismatch and not f.unitarity_inconsistent and not f.flagged and f.reasons == ()
    assert f.classification == "DIAGNOSTIC / NOT B4"


def test_T3_1_residual_flag_alone():
    f = rgf.point_flags({("x", "T1"): _T_OK, ("x", "R1"): _R_SHIFTED}, {"x": (2.0 ** -40,) * 3}, ["x"])["x"]
    assert f.t1_status == f.r1_status == "RESOLVED" and not f.status_mismatch
    u = f.unitarity
    assert (u.sum_delta1, u.sum_delta2) == (2.0 ** -25, 2.0 ** -27)
    assert (u.bound1, u.bound2) == (2.0 ** -39, 2.0 ** -39)
    assert u.exceeds1 and u.exceeds2 and f.unitarity_inconsistent and f.flagged
    assert len(f.reasons) == 2 and all("U(1/" in s for s in f.reasons)


def test_T3_1_status_mismatch_alone():
    # sums (0, 2^-16); U = (0, 2^-17, 2^-17): bound2 = 2^-16 = |sum2| (equality is consistent), bound1 = 2^-17 >= 0
    f = rgf.point_flags({("x", "T1"): _T_OK, ("x", "R1"): _R_MISMATCH}, {"x": (0.0, 2.0 ** -17, 2.0 ** -17)},
                        ["x"])["x"]
    assert (f.t1_status, f.r1_status) == ("RESOLVED", "NOT RESOLVABLE") and f.status_mismatch
    assert f.unitarity.sum_delta2 == 2.0 ** -16 and f.unitarity.bound2 == 2.0 ** -16
    assert not f.unitarity_inconsistent and f.flagged and len(f.reasons) == 1 and "status mismatch" in f.reasons[0]


def test_T3_1_both_conditions_keep_both_reasons():
    f = rgf.point_flags({("x", "T1"): _T_OK, ("x", "R1"): _R_MISMATCH}, {"x": (0.0, 0.0, 0.0)}, ["x"])["x"]
    assert f.status_mismatch and f.unitarity_inconsistent and f.flagged
    assert f.unitarity.exceeds2 and not f.unitarity.exceeds1                # sum1 = 0 = bound1
    assert len(f.reasons) == 2


def test_T3_1_equality_at_the_derived_boundary_is_not_flagged():
    t, r = (0.5, 0.5 + 2.0 ** -20, 0.5 + 2.0 ** -20), (0.5, 0.5, 0.5)     # sums (2^-20, 0)
    eq = rgf.unitarity_consistency(t, r, (2.0 ** -21, 2.0 ** -21, 0.0))   # bound1 = 2^-20 exactly
    assert eq.sum_delta1 == 2.0 ** -20 == eq.bound1 and not eq.exceeds1 and not eq.flagged
    over = rgf.unitarity_consistency(t, r, (2.0 ** -22, 2.0 ** -21, 0.0))  # bound1 = 0.75 * 2^-20
    assert over.bound1 == 0.75 * 2.0 ** -20 and over.exceeds1 and over.flagged
    # report-only per-dx defect |1 - T1 - R1| - U, by hand: 1/512: 0 - 2^-21; 1/1024: 2^-20 - 2^-21; 1/2048: 2^-20 - 0.
    # Positive entries mean these synthetic data violate |T1 + R1 - 1| <= U (they are deliberately inconsistent).
    assert eq.per_dx_defect == (-(2.0 ** -21), 2.0 ** -21, 2.0 ** -20)


def test_T3_1_missing_or_invalid_residuals_rejected():
    v = {("x", "T1"): _T_OK, ("x", "R1"): _R_OK}
    with pytest.raises(ValueError, match="missing residuals"):
        rgf.point_flags(v, {}, ["x"])
    with pytest.raises(ValueError, match="unexpected residuals"):
        rgf.point_flags(v, {"x": (0.0,) * 3, "y": (0.0,) * 3}, ["x"])
    with pytest.raises(ValueError, match="missing values"):
        rgf.point_flags({("x", "T1"): _T_OK}, {"x": (0.0,) * 3}, ["x"])
    for bad in ((0.0, -1e-16, 0.0), (0.0, float("nan"), 0.0), (0.0, 0.0), None):
        with pytest.raises((ValueError, TypeError)):
            rgf.point_flags(v, {"x": bad}, ["x"])
    with pytest.raises(TypeError):
        rgf.unitarity_consistency(_T_OK, _R_OK, (0.0, "1e-16", 0.0))


def test_T3_1_status_mismatch_agrees_with_select_resolution():
    v = {("x1", "T1"): _T_OK, ("x1", "R1"): _R_OK, ("x2", "T1"): _T_OK, ("x2", "R1"): _R_MISMATCH}
    flags = rgf.point_flags(v, {"x1": (0.0,) * 3, "x2": (0.0,) * 3}, ["x1", "x2"])
    sel = rgf.select_resolution(v, ["x1", "x2"], solver_agreement=True)
    assert tuple(x for x in ("x1", "x2") if flags[x].status_mismatch) == sel.status_mismatch_points == ("x2",)


# ======================================================================
# T5 — output guard, labels, metadata, clean-tree abort, atomic publication (b4_pilot_io)
# No physical calculation; G0 is NOT executed. G0Outcome objects below are synthetic test records.
# ======================================================================
import shutil  # noqa: E402

import b4_pilot_io as pio  # noqa: E402

_PINNED = "sha256:3cd620d2affc6ba5e59e6bd26bdd407afdc6eed7269ee91dd5448a30a0a035ab"   # independent literal (plan R1)


def _git_tmp(repo, *args):
    r = subprocess.run(["git", "-C", str(repo), "-c", "user.email=pilot@example.invalid", "-c", "user.name=pilot",
                        "-c", "commit.gpgsign=false", *args], capture_output=True)
    assert r.returncode == 0, r.stderr.decode()
    return r.stdout.decode().strip()


@pytest.fixture
def tmp_repo(tmp_path):
    """Throwaway git repository containing a byte copy of the committed protocol and an official results dir."""
    repo = tmp_path / "repo"
    (repo / "docs").mkdir(parents=True)
    (repo / "docs" / "B4_PILOT_PROTOCOL.md").write_bytes(_baseline_contents([PROTOCOL])[PROTOCOL])
    (repo / "results" / "rgf_b3_validation").mkdir(parents=True)
    (repo / "results" / "rgf_b3_validation" / "summary.json").write_text("{}")
    _git_tmp(repo, "init", "-q")
    _git_tmp(repo, "add", "-A")
    _git_tmp(repo, "commit", "-q", "-m", "init")
    return repo


def _scratch(tmp_path):
    d = tmp_path / "scratch"
    d.mkdir(exist_ok=True)
    return str(d)


def test_T5_pinned_protocol_reference_is_independent_and_correct():
    assert pio.PROTOCOL_SHA256_REFERENCE == _PINNED
    assert pio.A6_BASELINE_COMMIT == BASELINE and pio.A6_BASELINE_COMMIT != pio.PROTOCOL_SHA256_REFERENCE
    # consistency of the pinned literal with the committed file (computed here from the baseline blob, not by pio)
    assert _sha256(_baseline_contents([PROTOCOL])[PROTOCOL]) == _PINNED
    src = open(os.path.join(REPO, "b4_pilot_io.py"), encoding="utf-8").read()
    assert f'PROTOCOL_SHA256_REFERENCE = "{_PINNED}"' in src        # a literal, not computed at import


def test_T5_check_protocol_on_the_real_repository_is_read_only():
    assert pio.check_protocol(REPO) == _PINNED


def test_T5_check_protocol_detects_a_changed_byte(tmp_repo):
    assert pio.check_protocol(str(tmp_repo)) == _PINNED
    p = tmp_repo / "docs" / "B4_PILOT_PROTOCOL.md"
    p.write_bytes(p.read_bytes() + b" ")
    with pytest.raises(pio.PilotGuardError, match="protocol SHA-256 mismatch"):
        pio.check_protocol(str(tmp_repo))
    with pytest.raises(pio.PilotGuardError):
        pio.check_protocol(str(tmp_repo), reference="sha256:" + "0" * 64)


def test_T5_clean_tree_abort_on_untracked_and_modified(tmp_repo):
    head = pio.require_clean_tree(str(tmp_repo))
    assert len(head) == 40
    (tmp_repo / "new_file.py").write_text("x = 1\n")
    with pytest.raises(pio.PilotGuardError, match="not clean"):
        pio.require_clean_tree(str(tmp_repo))
    (tmp_repo / "new_file.py").unlink()
    (tmp_repo / "results" / "rgf_b3_validation" / "summary.json").write_text('{"x": 1}')
    with pytest.raises(pio.PilotGuardError, match="not clean"):
        pio.require_clean_tree(str(tmp_repo))


def test_T5_output_guard_allowlist_and_denylist(tmp_repo, tmp_path):
    r = str(tmp_repo)
    assert pio.check_output_path(os.path.join(r, "results/b4_pilot/rgf_resolution/x.csv"), r)
    for bad in ("results/rgf_b3_validation/x.csv", "results/rgf_b3y_validation", "results/transport/a.npz",
                "results/validation/v.csv", "results/stationary_validation", "results/rgf_b4_validation/x",
                "results/other/x.csv", "docs/x.md", "x.csv", "results/b4_pilot_evil/x.csv"):
        with pytest.raises(pio.PilotGuardError):
            pio.check_output_path(os.path.join(r, bad), r)
    with pytest.raises(pio.PilotGuardError):
        pio.check_output_path(str(tmp_path / "outside.csv"), r)
    # a symlink inside results/b4_pilot that resolves into an official directory is refused
    (tmp_repo / "results" / "b4_pilot").mkdir()
    os.symlink(tmp_repo / "results" / "rgf_b3_validation", tmp_repo / "results" / "b4_pilot" / "link")
    with pytest.raises(pio.PilotGuardError, match="official result directory"):
        pio.check_output_path(os.path.join(r, "results/b4_pilot/link/x.csv"), r)


def test_T5_staging_must_be_outside_the_repository(tmp_repo, tmp_path):
    with pytest.raises(pio.PilotGuardError, match="outside the repository"):
        pio.Staging(str(tmp_repo), scratch_root=str(tmp_repo / "results"))
    with pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path)) as st:
        assert not st.path.startswith(os.path.realpath(str(tmp_repo)) + os.sep)
        with pytest.raises(pio.PilotGuardError, match="escapes"):
            st.write_csv("../diagnostic_x.csv", ["a"], [])


def test_T5_labelled_writers(tmp_repo, tmp_path):
    with pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path)) as st:
        c = st.write_csv("rgf_resolution/diagnostic_rows.csv", ["a", "b"], [{"a": 1, "b": 2}])
        with open(c, encoding="utf-8") as f:
            assert f.readline().rstrip("\n") == "# DIAGNOSTIC / NOT B4 — docs/B4_PILOT_PROTOCOL.md"
            assert f.readline().strip() == "a,b"
        j = st.write_json("diagnostic_gates.json", {"x": 1})
        assert json.load(open(j))["classification"] == "DIAGNOSTIC / NOT B4"
        n = st.write_npz("diagnostic_S.npz", S=np.eye(2, dtype=complex))
        with np.load(n) as z:
            assert str(z["classification"]) == "DIAGNOSTIC / NOT B4" and z["S"].shape == (2, 2)
        for bad in ("rows.csv", "diagnostic_rows.txt"):
            with pytest.raises(pio.PilotGuardError):
                st.write_csv(bad, ["a"], [])
        with pytest.raises(pio.PilotGuardError):
            st.write_json("diagnostic_bad.json", {"classification": "B4 RESULT"})
        with pytest.raises(pio.PilotGuardError):
            st.write_npz("diagnostic_bad.npz", classification=np.array("B4"))


def _meta(repo, head):
    return pio.build_metadata(str(repo), head=head, protected_before={"a": "sha256:x"},
                              code_files=("docs/B4_PILOT_PROTOCOL.md",), command="test")


def test_T5_metadata_fields(tmp_repo):
    head = pio.require_clean_tree(str(tmp_repo))
    m = _meta(tmp_repo, head)
    assert m["classification"] == "DIAGNOSTIC / NOT B4"
    assert m["protocol_sha256"] == _PINNED and m["a6_baseline_commit"] == BASELINE
    assert m["git_commit"] == head and m["clean_tree"] is True
    fp = m["provenance"]["code_fingerprint"]
    assert fp["code_files"] == {"docs/B4_PILOT_PROTOCOL.md": _PINNED} and fp["code_fingerprint"].startswith("sha256:")
    assert {"python", "numpy", "OPENBLAS_NUM_THREADS", "cpu_count"} <= set(m["environment"])
    with pytest.raises(pio.PilotGuardError):
        pio.build_metadata(str(tmp_repo), head=head, protected_before={}, code_files=("missing.py",))
    with pytest.raises(pio.PilotGuardError):
        pio.build_metadata(str(tmp_repo), head=head, protected_before={}, code_files=(),
                           extra={"protocol_sha256": "x"})


def test_T5_real_code_fingerprint_covers_the_pilot_modules():
    fp = pio.code_fingerprint(REPO)
    assert {"b4_pilot_rgf.py", "b4_pilot_io.py", "stationary.py", "reference_solver.py"} <= set(fp["code_files"])
    with open(os.path.join(REPO, "stationary.py"), "rb") as f:
        assert fp["code_files"]["stationary.py"] == _sha256(f.read())


def test_T5_failed_g0_publishes_nothing(tmp_repo, tmp_path):
    head = pio.require_clean_tree(str(tmp_repo))
    st = pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path))
    st.write_csv("diagnostic_controls.csv", ["x"], [{"x": 1}])
    with pytest.raises(pio.PilotGuardError, match="G0 has not passed"):
        pio.publish(st, pio.G0Outcome(False, _PINNED, head), _meta(tmp_repo, head))
    assert not (tmp_repo / "results" / "b4_pilot").exists()
    assert not os.path.exists(st.path) and st.discarded
    assert [n for n in os.listdir(tmp_repo / "results") if n.startswith(".b4_pilot")] == []
    assert pio.dirty_entries(str(tmp_repo)) == []                  # repository untouched


def test_T5_exception_before_publication_leaves_no_output(tmp_repo, tmp_path):
    with pytest.raises(RuntimeError, match="boom"):
        with pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path)) as st:
            st.write_csv("diagnostic_controls.csv", ["x"], [{"x": 1}])
            path = st.path
            raise RuntimeError("boom")
    assert not os.path.exists(path) and not (tmp_repo / "results" / "b4_pilot").exists()


@pytest.mark.parametrize("bad", ["wrong_hash", "head_mismatch", "unlabelled_file", "not_a_record"])
def test_T5_publication_refusals_leave_no_output(tmp_repo, tmp_path, bad):
    head = pio.require_clean_tree(str(tmp_repo))
    st = pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path))
    st.write_csv("diagnostic_controls.csv", ["x"], [{"x": 1}])
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    if bad == "wrong_hash":
        outcome = pio.G0Outcome(True, "sha256:" + "0" * 64, head)
    elif bad == "head_mismatch":
        meta = dict(meta, git_commit="0" * 40)
    elif bad == "unlabelled_file":
        with open(os.path.join(st.path, "stray.txt"), "w") as f:
            f.write("x")
    else:
        outcome = {"passed": True}
    with pytest.raises(pio.PilotGuardError):
        pio.publish(st, outcome, meta)
    assert not (tmp_repo / "results" / "b4_pilot").exists() and not os.path.exists(st.path)


def test_T5_successful_publication_is_complete_and_labelled(tmp_repo, tmp_path):
    head = pio.require_clean_tree(str(tmp_repo))
    with pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path)) as st:
        st.write_readme("synthetic test publication")
        st.write_csv("rgf_resolution/diagnostic_controls.csv", ["x"], [{"x": 1}])
        st.write_npz("rgf_resolution/diagnostic_S.npz", S=np.eye(2))
        target = pio.publish(st, pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head))
    assert target == os.path.join(os.path.realpath(str(tmp_repo)), "results", "b4_pilot")
    files = sorted(os.path.relpath(os.path.join(d, f), target) for d, _, fs in os.walk(target) for f in fs)
    assert files == ["README.md", "diagnostic_run_metadata.json", "rgf_resolution/diagnostic_S.npz",
                     "rgf_resolution/diagnostic_controls.csv"]
    meta = json.load(open(os.path.join(target, "diagnostic_run_metadata.json")))
    assert meta["protocol_sha256"] == _PINNED and meta["git_commit"] == head
    assert [n for n in os.listdir(tmp_repo / "results") if n.startswith(".b4_pilot")] == []
    assert not os.path.exists(st.path)
    # official directory untouched; a second publication is refused (no overwrite or merge)
    assert (tmp_repo / "results" / "rgf_b3_validation" / "summary.json").read_text() == "{}"
    st2 = pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path))
    with pytest.raises(pio.PilotGuardError, match="already exists"):
        pio.publish(st2, pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head))


def test_T5_leftover_partial_directory_blocks_publication(tmp_repo, tmp_path):
    head = pio.require_clean_tree(str(tmp_repo))
    (tmp_repo / "results" / ".b4_pilot.partial-deadbeef").mkdir()
    st = pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path))
    with pytest.raises(pio.PilotGuardError, match="leftover partial"):
        pio.publish(st, pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head))
    assert not (tmp_repo / "results" / "b4_pilot").exists()


def test_T5_snapshots_detect_changes(tmp_repo):
    paths = ["docs/B4_PILOT_PROTOCOL.md", "results/rgf_b3_validation/summary.json", "missing.txt"]
    before = pio.protected_snapshot(str(tmp_repo), paths)
    assert before["missing.txt"] is None and before["docs/B4_PILOT_PROTOCOL.md"] == _PINNED
    assert pio.compare_snapshots(before, pio.protected_snapshot(str(tmp_repo), paths)) == []
    (tmp_repo / "results" / "rgf_b3_validation" / "summary.json").write_text("[]")
    diff = pio.compare_snapshots(before, pio.protected_snapshot(str(tmp_repo), paths))
    assert len(diff) == 1 and "summary.json" in diff[0]


def test_T5_baseline_hashes_use_file_bytes_at_the_a6_commit():
    h = pio.baseline_hashes(REPO, [PROTOCOL, "stationary.py", "no/such/file"])
    assert h[PROTOCOL] == _PINNED and h["no/such/file"] is None
    assert h["stationary.py"] == _sha256(_baseline_contents(["stationary.py"])["stationary.py"])


def test_T5_real_results_b4_pilot_does_not_exist():
    assert not os.path.exists(os.path.join(REPO, "results", "b4_pilot"))
    assert [n for n in os.listdir(os.path.join(REPO, "results")) if n.startswith(".b4_pilot")] == []
