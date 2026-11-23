"""B4 limited diagnostic pilot — tests, phases T1–T3 (including T3.1), T5 and T7.

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
T7  b4_pilot_rgf_records: §1.4 record extraction from ONE method's raw output (full complex S; T1, R1, U from both
    products, Q, channel counts, method extras, per-method code fingerprint), the complete protocol record set and
    same-problem checks, raw M1-M2 / M1-M3 / M2-M3 agreement rows and the §1.3 M4 requirement (M1-M2 or M1-M3 only;
    an M2-M3-only exceedance is reported as a protocol gap), raw M4 adjudication inputs (no classification), and
    lossless serialization (complex128; clongdouble as exact hi/lo pairs) through T5 Staging. Synthetic S only.

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
    env = m["provenance"]["environment_info"]
    assert {"python", "platform", "machine", "processor", "cpu_count", "numpy", "scipy",
            "OPENBLAS_NUM_THREADS"} <= set(env)
    assert "environment" not in m and m["dirty_entries"] == []
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
    meta0 = _meta(tmp_repo, head)
    with pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path)) as st:
        st.write_readme("synthetic test publication")
        st.write_csv("rgf_resolution/diagnostic_controls.csv", ["x"], [{"x": 1}])
        st.write_npz("rgf_resolution/diagnostic_S.npz", S=np.eye(2))
        target = pio.publish(st, pio.G0Outcome(True, _PINNED, head), meta0)
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
    # (the published, uncommitted output now makes the tree dirty, so fresh metadata is refused too)
    with pytest.raises(pio.PilotGuardError, match="not clean"):
        _meta(tmp_repo, head)
    st2 = pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path))
    with pytest.raises(pio.PilotGuardError, match="already exists"):
        pio.publish(st2, pio.G0Outcome(True, _PINNED, head), meta0)
    assert not os.path.exists(st2.path)


def _no_partial(repo):
    return [n for n in os.listdir(repo / "results") if n.startswith(".b4_pilot")] == []


def _staged(repo, tmp_path):
    st = pio.Staging(str(repo), scratch_root=_scratch(tmp_path))
    st.write_csv("diagnostic_controls.csv", ["x"], [{"x": 1}])
    return st


def _assert_nothing_published(repo, st):
    assert not os.path.lexists(repo / "results" / "b4_pilot") and _no_partial(repo)
    assert not os.path.exists(st.path) and st.discarded


def test_T5_head_change_after_g0_is_refused(tmp_repo, tmp_path):
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    (tmp_repo / "later.txt").write_text("x")
    _git_tmp(tmp_repo, "add", "later.txt")
    _git_tmp(tmp_repo, "commit", "-q", "-m", "after G0")                # clean tree, but HEAD moved
    st = _staged(tmp_repo, tmp_path)
    with pytest.raises(pio.PilotGuardError, match="HEAD changed since G0"):
        pio.publish(st, outcome, meta)
    _assert_nothing_published(tmp_repo, st)


@pytest.mark.parametrize("kind", ["untracked", "modified"])
def test_T5_dirty_tree_after_g0_is_refused(tmp_repo, tmp_path, kind):
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    if kind == "untracked":
        (tmp_repo / "stray.py").write_text("x = 1\n")
    else:
        (tmp_repo / "results" / "rgf_b3_validation" / "summary.json").write_text('{"changed": 1}')
    st = _staged(tmp_repo, tmp_path)
    with pytest.raises(pio.PilotGuardError, match="not clean"):
        pio.publish(st, outcome, meta)
    _assert_nothing_published(tmp_repo, st)


def test_T5_protocol_change_after_g0_is_refused(tmp_repo, tmp_path):
    p = tmp_repo / "docs" / "B4_PILOT_PROTOCOL.md"
    p.write_bytes(p.read_bytes() + b"\n")
    _git_tmp(tmp_repo, "commit", "-q", "-am", "protocol changed")
    head = pio.require_clean_tree(str(tmp_repo))                       # outcome claims the new HEAD
    outcome = pio.G0Outcome(True, _PINNED, head)
    meta = _meta(tmp_repo, head)
    st = _staged(tmp_repo, tmp_path)
    with pytest.raises(pio.PilotGuardError, match="protocol SHA-256 mismatch"):
        pio.publish(st, outcome, meta)
    _assert_nothing_published(tmp_repo, st)


def test_T5_final_revalidation_runs_last_before_the_rename(tmp_repo, tmp_path, monkeypatch):
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    st = _staged(tmp_repo, tmp_path)
    calls = []
    orig_json, orig_verify = st.write_json, pio._verify_labels
    orig_reval, orig_rename = pio.revalidate_for_publication, pio._rename_noreplace
    orig_resolve = pio._resolve_renameat2
    monkeypatch.setattr(pio, "_resolve_renameat2", lambda: (calls.append("resolve"), orig_resolve())[1])
    monkeypatch.setattr(st, "write_json", lambda *a, **k: (calls.append("metadata"), orig_json(*a, **k))[1])
    monkeypatch.setattr(pio, "_verify_labels", lambda *a: (calls.append("labels"), orig_verify(*a))[1])
    monkeypatch.setattr(pio, "revalidate_for_publication", lambda *a: (calls.append("revalidate"), orig_reval(*a))[1])
    monkeypatch.setattr(pio, "_rename_noreplace", lambda *a: (calls.append("rename"), orig_rename(*a))[1])
    target = pio.publish(st, outcome, meta)
    assert calls == ["metadata", "labels", "resolve", "revalidate", "rename"]
    assert os.path.isdir(target) and _no_partial(tmp_repo)


def test_T5_tree_dirtied_during_label_check_is_caught(tmp_repo, tmp_path, monkeypatch):
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    st = _staged(tmp_repo, tmp_path)
    orig_verify = pio._verify_labels

    def verify_then_dirty(root):
        orig_verify(root)
        (tmp_repo / "late_change.txt").write_text("x")

    monkeypatch.setattr(pio, "_verify_labels", verify_then_dirty)
    with pytest.raises(pio.PilotGuardError, match="not clean"):
        pio.publish(st, outcome, meta)
    _assert_nothing_published(tmp_repo, st)


def test_T5_concurrent_target_is_never_overwritten(tmp_repo, tmp_path, monkeypatch):
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    st = _staged(tmp_repo, tmp_path)
    orig_rename = pio._rename_noreplace
    target = tmp_repo / "results" / "b4_pilot"

    def race(fn, src, dst):
        os.mkdir(dst)                                               # a concurrent process creates an empty target
        orig_rename(fn, src, dst)

    monkeypatch.setattr(pio, "_rename_noreplace", race)
    with pytest.raises(pio.PilotGuardError, match="never overwritten"):
        pio.publish(st, outcome, meta)
    assert target.is_dir() and os.listdir(target) == []            # untouched and still empty
    assert not os.path.exists(st.path) and _no_partial(tmp_repo)


@pytest.mark.parametrize("mode", ["missing", "ENOSYS", "EINVAL", "EXDEV", "ENOTSUP", "EOPNOTSUPP", "EACCES",
                                  "UNKNOWN"])
def test_T5_no_replace_primitive_failures_fail_closed(tmp_repo, tmp_path, monkeypatch, mode):
    import ctypes
    import errno as _errno
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    st = _staged(tmp_repo, tmp_path)
    if mode == "missing":
        monkeypatch.setattr(pio, "_libc_renameat2", lambda: None)
    else:
        code = 9999 if mode == "UNKNOWN" else getattr(_errno, mode)   # 9999: an errno with no special handling

        def fake(*args):
            ctypes.set_errno(code)
            return -1

        monkeypatch.setattr(pio, "_libc_renameat2", lambda: fake)
    with pytest.raises(pio.PilotGuardError, match="fails? closed"):
        pio.publish(st, outcome, meta)
    _assert_nothing_published(tmp_repo, st)


def test_T5_cross_filesystem_staging_is_refused(tmp_repo, tmp_path, monkeypatch):
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    st = _staged(tmp_repo, tmp_path)
    real = pio._device
    monkeypatch.setattr(pio, "_device", lambda p: real(p) + (1 if os.path.realpath(p) == st.path else 0))
    with pytest.raises(pio.PilotGuardError, match="cross-filesystem"):
        pio.publish(st, outcome, meta)
    _assert_nothing_published(tmp_repo, st)
    assert not os.path.exists(os.path.join(st.path, "diagnostic_run_metadata.json"))


def test_T5_device_helper_and_real_no_replace_rename(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir()
    (a / "f").write_text("1")
    assert pio._device(str(a)) == os.stat(str(a)).st_dev
    fn = pio._resolve_renameat2()                                     # available on this platform (glibc 2.39)
    pio._rename_noreplace(fn, os.fsencode(str(a)), os.fsencode(str(b)))
    assert not a.exists() and (b / "f").read_text() == "1"
    a.mkdir()
    (a / "g").write_text("2")
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(pio.PilotGuardError, match="never overwritten"):
        pio._rename_noreplace(fn, os.fsencode(str(a)), os.fsencode(str(empty)))
    assert os.listdir(empty) == [] and (a / "g").exists()


def test_T5_same_filesystem_publication_adds_only_the_target(tmp_repo, tmp_path):
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    before = set(os.listdir(tmp_repo / "results"))
    st = _staged(tmp_repo, tmp_path)
    assert pio._device(st.path) == pio._device(str(tmp_repo / "results"))
    pio.publish(st, outcome, meta)
    assert set(os.listdir(tmp_repo / "results")) - before == {"b4_pilot"} and _no_partial(tmp_repo)
    assert not hasattr(pio, "PARTIAL_PREFIX")


def test_T5_build_metadata_refuses_dirty_tree_and_head_mismatch(tmp_repo):
    head = pio.require_clean_tree(str(tmp_repo))
    with pytest.raises(pio.PilotGuardError, match="HEAD is"):
        pio.build_metadata(str(tmp_repo), head="0" * 40, protected_before={}, code_files=())
    (tmp_repo / "stray.txt").write_text("x")
    with pytest.raises(pio.PilotGuardError, match="not clean"):
        pio.build_metadata(str(tmp_repo), head=head, protected_before={}, code_files=())


@pytest.mark.parametrize("bad", ["clean_tree_false", "clean_tree_missing", "env_missing", "provenance_missing",
                                 "env_not_dict"])
def test_T5_metadata_record_checks_refuse_publication(tmp_repo, tmp_path, bad):
    head = pio.require_clean_tree(str(tmp_repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(tmp_repo, head)
    if bad == "clean_tree_false":
        meta["clean_tree"] = False
    elif bad == "clean_tree_missing":
        del meta["clean_tree"]
    elif bad == "env_missing":
        del meta["provenance"]["environment_info"]
    elif bad == "provenance_missing":
        del meta["provenance"]
    else:
        meta["provenance"]["environment_info"] = "python 3"
    st = _staged(tmp_repo, tmp_path)
    with pytest.raises(pio.PilotGuardError, match="clean tree|environment_info"):
        pio.publish(st, outcome, meta)
    _assert_nothing_published(tmp_repo, st)


def test_T5_filesystem_mismatch_is_detected_before_any_staging_write(tmp_repo, tmp_path, monkeypatch):
    scratch = _scratch(tmp_path)
    real = pio._device
    monkeypatch.setattr(pio, "_device", lambda p: real(p) + (1 if os.path.realpath(p) == os.path.realpath(scratch)
                                                              else 0))
    with pytest.raises(pio.PilotGuardError, match="different filesystems"):
        pio.Staging(str(tmp_repo), scratch_root=scratch)
    assert os.listdir(scratch) == []                                    # not even the staging directory was created
    assert not os.path.lexists(tmp_repo / "results" / "b4_pilot") and _no_partial(tmp_repo)


def test_T5_staging_refuses_missing_results_directory(tmp_repo, tmp_path):
    import shutil
    shutil.rmtree(tmp_repo / "results")
    scratch = _scratch(tmp_path)
    with pytest.raises(pio.PilotGuardError, match="results/ does not exist"):
        pio.Staging(str(tmp_repo), scratch_root=scratch)
    assert os.listdir(scratch) == [] and not (tmp_repo / "results").exists()


class _PostRevalidationTracer:
    """Records what ACTUALLY runs in publish() after the final revalidation returns, up to the renameat2 foreign call.

    Two independent instruments are armed the moment revalidate_for_publication returns:
      * a CPython audit hook (PEP 578): sees subprocess/posix_spawn/exec, open, os.listdir/scandir/mkdir/rename/remove,
        ctypes.dlopen/dlsym, ... even when triggered from C code;
      * sys.setprofile: sees every Python-level call and every builtin/C call (os.stat, os.lstat, os.fsencode, ...).
    The window ends when the resolved renameat2 function is invoked. Calling a ctypes function pointer raises no audit
    event, so the REAL lookup (find_library, dlopen, dlsym) is left untouched and only its result is wrapped in a
    marker proxy that forwards to the real foreign function; entering the proxy closes the window.
    The full audit sequence from the start of publish() is kept as well, to check where the lookup happens.
    The profiler active before the tracer installs its own (sys.getprofile(), this thread) is saved and reinstated by
    close(); close() is idempotent and leaves profiling untouched if the tracer never installed its profiler."""

    _hook_installed = False
    current = None

    def __init__(self):
        self.armed = False
        self.window = []          # everything observed after revalidation, before the renameat2 call
        self.sequence = []        # audit events + markers for the whole publish() call
        self.rename_called = False
        self._profile_installed = False
        self._previous_profile = None
        if not _PostRevalidationTracer._hook_installed:   # audit hooks cannot be removed: install one dispatcher once
            sys.addaudithook(_PostRevalidationTracer._dispatch)
            _PostRevalidationTracer._hook_installed = True

    @staticmethod
    def _dispatch(event, args):
        t = _PostRevalidationTracer.current
        if t is not None:
            t._audit(event, args)

    def _audit(self, event, args):
        if event == "ctypes.dlsym":
            self.sequence.append(("dlsym", args[1] if len(args) > 1 else None))
        elif event in ("subprocess.Popen", "os.posix_spawn", "os.exec", "ctypes.dlopen", "open"):
            self.sequence.append((event, None))
        if not self.armed or event == "sys.setprofile":
            return
        if event == "object.__getattr__" and _is_tracer_introspection(args):
            return                                      # the profiler's own frame.f_code read, not publish()'s
        self.window.append(("audit", event))

    def _profile(self, frame, what, arg):
        if not self.armed:
            return
        if what == "call":
            if frame.f_code is _MARKER_CODE:                        # the renameat2 invocation: close the window
                self.armed = False
                self.rename_called = True
                self.sequence.append(("RENAME_CALL", None))
                return
            self.window.append(("call", frame.f_code.co_name))
        elif what == "c_call":
            self.window.append(("c_call", getattr(arg, "__name__", repr(arg))))

    def install(self, monkeypatch):
        orig = pio.revalidate_for_publication
        real_lookup = pio._libc_renameat2

        def traced_lookup():
            fn = real_lookup()                      # the real find_library / dlopen / dlsym
            return None if fn is None else _RenameMarker(fn)

        monkeypatch.setattr(pio, "_libc_renameat2", traced_lookup)

        def traced_revalidate(*a):
            self.sequence.append(("REVALIDATE_START", None))
            orig(*a)
            self.sequence.append(("REVALIDATE_END", None))
            if not self._profile_installed:
                self._previous_profile = sys.getprofile()   # preserved and reinstated by close()
                self._profile_installed = True
            sys.setprofile(self._profile)
            self.armed = True                       # nothing after this line but returning into publish()

        monkeypatch.setattr(pio, "revalidate_for_publication", traced_revalidate)
        _PostRevalidationTracer.current = self
        return self

    def close(self):
        self.armed = False
        if self._profile_installed:                     # restore exactly once; later calls are no-ops
            self._profile_installed = False
            previous, self._previous_profile = self._previous_profile, None
            sys.setprofile(previous)
        _PostRevalidationTracer.current = None


class _RenameMarker:
    """Forwards to the real resolved renameat2 foreign function; only marks the moment it is invoked."""

    def __init__(self, fn):
        self.fn = fn

    def __call__(self, *args):
        return self.fn(*args)


_MARKER_CODE = _RenameMarker.__call__.__code__


def _is_tracer_introspection(args):
    """CPython audits reads of frame.f_code; the profiler performs exactly those reads. Nothing else is ignored."""
    import types
    return len(args) == 2 and args[1] == "f_code" and isinstance(args[0], types.FrameType)


@pytest.fixture
def tracer(monkeypatch):
    t = _PostRevalidationTracer().install(monkeypatch)
    yield t
    t.close()


def _traced_publication(repo, tmp_path):
    head = pio.require_clean_tree(str(repo))
    outcome, meta = pio.G0Outcome(True, _PINNED, head), _meta(repo, head)
    return _staged(repo, tmp_path), outcome, meta


def test_T5_nothing_runs_between_final_revalidation_and_the_rename(tmp_repo, tmp_path, tracer):
    st, outcome, meta = _traced_publication(tmp_repo, tmp_path)
    target = pio.publish(st, outcome, meta)
    tracer.close()
    assert tracer.rename_called and os.path.isdir(target)
    # the only thing observed after the revalidation is entering _rename_noreplace, whose first operation is the call
    assert tracer.window == [("call", "_rename_noreplace")], tracer.window


def test_T5_renameat2_is_resolved_once_before_the_final_revalidation(tmp_repo, tmp_path, tracer):
    st, outcome, meta = _traced_publication(tmp_repo, tmp_path)
    pio.publish(st, outcome, meta)
    tracer.close()
    names = [n for n, _ in tracer.sequence]
    lookups = [i for i, (n, a) in enumerate(tracer.sequence) if n == "dlsym" and a == "renameat2"]
    assert len(lookups) == 1, tracer.sequence
    assert lookups[0] < names.index("REVALIDATE_START") < names.index("REVALIDATE_END") < names.index("RENAME_CALL")
    after = tracer.sequence[names.index("REVALIDATE_END") + 1:]
    assert after == [("RENAME_CALL", None)], after


@pytest.mark.parametrize("intervening", ["lookup", "stat", "open", "subprocess"])
def test_T5_post_revalidation_detector_catches_intervening_operations(tmp_repo, tmp_path, monkeypatch, tracer,
                                                                      intervening):
    """Negative controls: the detector is not a tautology. Re-inserting each kind of operation between the final
    revalidation and the rename (including the pre-correction library lookup) is caught."""
    orig = pio._rename_noreplace

    def with_intervening(fn, src, dst):
        if intervening == "lookup":
            fn = pio._libc_renameat2()              # the pre-correction behaviour (find_library -> ldconfig)
        elif intervening == "stat":
            os.path.lexists(dst)
        elif intervening == "open":
            with open(os.devnull, "rb"):
                pass
        else:
            subprocess.run(["true"], check=True)
        return orig(fn, src, dst)

    monkeypatch.setattr(pio, "_rename_noreplace", with_intervening)
    st, outcome, meta = _traced_publication(tmp_repo, tmp_path)
    pio.publish(st, outcome, meta)
    tracer.close()
    seen = set(tracer.window)
    assert tracer.window != [("call", "_rename_noreplace")]
    expected = {"lookup": ("audit", "ctypes.dlopen"), "stat": ("c_call", "lstat"), "open": ("audit", "open"),
                "subprocess": ("audit", "subprocess.Popen")}[intervening]
    assert expected in seen, tracer.window


def test_T5_tracer_restores_the_previously_active_profiler(tmp_repo, tmp_path, tracer):
    original_profiler = sys.getprofile()                            # restored in the finally block
    calls = []

    def previous_profiler(frame, event, arg):
        calls.append(event)

    st, outcome, meta = _traced_publication(tmp_repo, tmp_path)
    sys.setprofile(previous_profiler)
    try:
        pio.publish(st, outcome, meta)
        assert sys.getprofile() is not previous_profiler            # replaced by the tracer during the window
        tracer.close()
        assert sys.getprofile() is previous_profiler                # reinstated, not disabled
        tracer.close()                                              # idempotent: a second close changes nothing
        assert sys.getprofile() is previous_profiler
        assert tracer.rename_called and tracer.window == [("call", "_rename_noreplace")], tracer.window
        assert calls                                                # the previous profiler really was active
    finally:
        sys.setprofile(original_profiler)


def test_T5_tracer_close_without_profiling_leaves_the_profiler_untouched(tracer):
    original_profiler = sys.getprofile()                            # restored in the finally block

    def previous_profiler(frame, event, arg):
        pass

    sys.setprofile(previous_profiler)
    try:
        tracer.close()                                              # the tracer never installed its profiler
        assert sys.getprofile() is previous_profiler
    finally:
        sys.setprofile(original_profiler)


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


# ======================================================================================================================
# T7 — b4_pilot_rgf_records: RGF record extraction, completeness, raw agreement / M4 inputs, lossless serialization.
# Synthetic hand-built S matrices only: no solver call, no physical data, no pilot gate.
# ======================================================================================================================
import dataclasses  # noqa: E402
import hashlib as _hashlib  # noqa: E402

import b4_pilot_rgf_records as rr  # noqa: E402
import sparse_exact_lattice as _xp  # noqa: E402

_CLD = _xp.CLD
_A = 0j
_B = 1 + 0j
_ALL_CODE_FILES = sorted({f for fs in rr.METHOD_CODE_FILES.values() for f in fs} | set(rr.EXTRACTION_CODE_FILES))
_CODE_FILES = {f: "sha256:" + _hashlib.sha256(f.encode()).hexdigest() for f in _ALL_CODE_FILES}
_META = {"provenance": {"code_fingerprint": {"code_files": dict(_CODE_FILES)}}}
_EVAL = {"evaluable": True}
_NOT_EVAL = {"evaluable": False}
_NSL = {512: 719, 1024: 1437, 2048: 2873}


def _S0():
    """Exactly unitary, symmetric 2x2 S (n_open = 1) with exactly representable moduli: T1 = 1, R1 = 0."""
    return np.array([[_A, _B], [_B, _A]], dtype=np.complex128)


def _psha(config, inv, tv):
    return "sha256:" + _hashlib.sha256(f"V|{config}|{inv}|{tv}".encode()).hexdigest()


def _rec(config="A", point="L1", tv=(63, 127), inv=512, method="M1", S=None, E=None, x_start=-0.5,
         platform=_EVAL, code_files=_CODE_FILES, **m4kw):
    S = _S0() if S is None else S
    common = dict(config=config, point=point, inv_dx=inv, n_y=tv[0], n_q=tv[1],
                  E=(float(10 + rr.POINT_ORDER.index(point)) if point in rr.POINT_ORDER else 1.0) if E is None else E, kx=1.25, V_lead_max=0.0,
                  potential_sha256=_psha(config, inv, tv), n_slices=_NSL.get(inv, 1), x_start=x_start, wall_time_s=0.5,
                  code_files=code_files)
    if method == "M4":
        out = dict(S=np.asarray(S).astype(_CLD), n_open=S.shape[0] // 2, error_indicator=1e-13,
                   refinement_gate={"passed": True}, green_refinement_gate={"passed": True},
                   fl_vs_state_max_abs=0.0, eta_final=1e-20)
        out.update(m4kw)
        return rr.extract_m4_record(solver_out=out, platform_report=platform, **common)
    out = dict(S=S, n_open=S.shape[0] // 2)
    if method == "M1":
        out.update(min_rcond=0.125, argmin_rcond=7)
    elif method == "M3":
        out.update(residual_after_refinement=2.0 ** -60, residual_before_refinement=2.0 ** -50)
    return rr.extract_record(method=method, solver_out=out, **common)


def _full_set(config="A", overrides=None):
    """The complete 145-key non-M4 record set for one configuration; overrides maps key -> S."""
    overrides = overrides or {}
    out = []
    for k in rr.expected_record_keys([config]):
        c, p, ny, nq, inv, m = k
        out.append(_rec(config=c, point=p, tv=(ny, nq), inv=inv, method=m, S=overrides.get(k)))
    return out


def _bumped(delta, i=0, j=1):
    S = _S0()
    S[i, j] += delta
    return S


# ---- extraction ----------------------------------------------------------------------------------------------
def test_T7_extraction_exact_T1_R1_U_on_a_dyadic_unitary_S():
    r = _rec()
    assert r.T1 == 1.0 and r.R1 == 0.0 and r.one_minus_T1_R1 == 0.0
    assert r.U == 0.0 and r.U_StS == 0.0 and r.U_SSt == 0.0 and r.Q == 0.0
    assert (r.n_open, r.n_closed) == (1, 62) and r.classification == "DIAGNOSTIC / NOT B4"
    H = 0.5 * np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]], dtype=np.complex128)
    h = _rec(S=H)                                                      # n_open = 2, real dyadic, unitary, symmetric
    assert (h.T1, h.R1, h.one_minus_T1_R1, h.U, h.Q) == (0.5, 0.5, 0.0, 0.0, 0.0) and h.n_closed == 61


def test_T7_observables_are_the_existing_definitions():
    rng = np.random.default_rng(7)
    S = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / 4
    r = _rec(S=S)
    assert r.T1 == b3.T1(S, 2) == float(np.sum(np.abs(S[2:, 0]) ** 2))          # Σ_m |(t_RL)_m1|²
    assert r.R1 == float(np.sum(np.abs(b3.sm.split_S(S, 2)["r_L"][:, 0]) ** 2))  # B3Y expression
    assert r.R1 == float(np.sum(np.abs(S[:2, 0]) ** 2))
    u, q = b3.unitarity(S), b3.reciprocity(S, 2)
    assert (r.U, r.U_StS, r.U_SSt) == (u["U"], u["U_StS"], u["U_SSt"]) and r.U == max(r.U_StS, r.U_SSt)
    assert all(getattr(r, k) == q[k] for k in ("Q", "Q_full", "Q_rL", "Q_rR", "Q_t", "Q_prob",
                                                "same_dir_RL_report_only", "same_dir_LR_report_only"))


def test_T7_T1_R1_use_only_incident_left_channel_1():
    rng = np.random.default_rng(3)
    S = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / 4
    S2 = S.copy()
    S2[:, 1:] *= 3                                                     # other incident channels changed
    a, b = _rec(S=S), _rec(S=S2)
    assert (a.T1, a.R1) == (b.T1, b.R1) and a.U != b.U


def test_T7_caller_supplied_observables_are_refused():
    r = _rec()
    for f, v in (("T1", 0.25), ("R1", 0.75), ("U", 0.0 + 1e-30), ("Q", 1.0), ("one_minus_T1_R1", 1e-16)):
        if getattr(r, f) == v:
            continue
        with pytest.raises(rr.PilotRecordError, match="recomputed from S"):
            dataclasses.replace(r, **{f: v})
    with pytest.raises(rr.PilotRecordError, match="n_closed"):
        dataclasses.replace(r, n_closed=61)


# ---- full complex S ------------------------------------------------------------------------------------------
def test_T7_full_complex_S_is_preserved_bitwise_and_read_only():
    S = _S0()
    r = _rec(S=S)
    assert r.S.dtype == np.complex128 and r.S.tobytes() == S.tobytes()
    S[0, 0] = 0.25
    assert r.S[0, 0] == _A                                             # a copy, not a view
    with pytest.raises(ValueError):
        r.S[0, 0] = 1
    with pytest.raises(ValueError):
        r.S.setflags(write=True)                                       # immutable backing buffer


def test_T7_a_phase_change_alters_the_stored_S_not_the_probabilities():
    a, b = _rec(), _rec(S=_S0() * 1j)                                  # exact global phase i
    assert (a.T1, a.R1, a.U) == (b.T1, b.R1, b.U) and a.s_sha256 != b.s_sha256


# ---- record validation and placement --------------------------------------------------------------------------
@pytest.mark.parametrize("kw", [dict(point="B3-E1", method="M2"), dict(point="B3-E2", method="M3"),
                                dict(point="B3-E1", method="M4"), dict(tv=(63, 254), point="U1", inv=1024,
                                                                         method="M2"),
                                dict(tv=(127, 254), point="U1", inv=1024, method="M4"),
                                dict(tv=(63, 254), point="L1", inv=1024), dict(tv=(63, 254), point="U1", inv=512),
                                dict(tv=(255, 254), point="U1", inv=1024), dict(point="X9"), dict(config="E"),
                                dict(inv=4096), dict(method="M5")])
def test_T7_placement_rules_refuse(kw):
    with pytest.raises(rr.PilotRecordError):
        _rec(**kw)


@pytest.mark.parametrize("bad", ["dtype", "shape", "nonfinite", "n_open_too_large", "S_not_array"])
def test_T7_invalid_S_is_refused(bad):
    r = _rec()
    S = {"dtype": _S0().astype(np.complex64), "shape": np.zeros((2, 3), complex),
         "nonfinite": np.array([[np.nan, _B], [_B, _A]]), "S_not_array": [[_A, _B], [_B, _A]]}.get(bad)
    if bad == "n_open_too_large":
        with pytest.raises(rr.PilotRecordError, match="n_open"):
            _rec(S=np.eye(126, dtype=np.complex128))                   # n_open = 63 = n_y
        return
    with pytest.raises(rr.PilotRecordError):
        dataclasses.replace(r, S=S)


@pytest.mark.parametrize("method,extras", [("M1", dict(min_rcond=0.1, argmin_rcond=1, error_indicator=0.0)),
                                           ("M1", dict(min_rcond=0.1)), ("M2", dict(min_rcond=0.1)),
                                           ("M3", dict(residual_after_refinement=0.0)),
                                           ("M1", dict(min_rcond=0.1, argmin_rcond=1.5))])
def test_T7_method_extras_must_match_exactly(method, extras):
    r = _rec(method=method)
    with pytest.raises(rr.PilotRecordError):
        dataclasses.replace(r, extras=extras)


def test_T7_extract_record_refuses_M4_and_incomplete_solver_output():
    with pytest.raises(rr.PilotRecordError, match="extract_m4_record"):
        rr.extract_record(config="A", point="L1", inv_dx=512, method="M4", n_y=63, n_q=127, E=1.0, kx=1.0,
                          solver_out={}, V_lead_max=0.0, potential_sha256=_psha("A", 512, (63, 127)), n_slices=1,
                          x_start=0.0, wall_time_s=0.0, code_files=_CODE_FILES)
    with pytest.raises(rr.PilotRecordError, match="lacks a required field"):
        rr.extract_record(config="A", point="L1", inv_dx=512, method="M1", n_y=63, n_q=127, E=1.0, kx=1.0,
                          solver_out={"S": _S0(), "n_open": 1}, V_lead_max=0.0,
                          potential_sha256=_psha("A", 512, (63, 127)), n_slices=1, x_start=0.0, wall_time_s=0.0,
                          code_files=_CODE_FILES)


# ---- completeness, separation and three-grid marshalling -------------------------------------------------------
def test_T7_expected_record_set_sizes_and_order():
    k1, k4 = rr.expected_record_keys(["A"]), rr.expected_record_keys(["D", "B", "A", "C"])
    assert len(k1) == 145 and len(set(k1)) == 145 and len(k4) == 580
    assert k4[:145] == k1 and k1[0] == ("A", "L1", 63, 127, 512, "M1")
    assert sum(1 for k in k1 if k[5] == "M1") == 14 * 3 + 9 + 10 and not any(k[5] == "M4" for k in k1)
    assert sum(1 for k in k1 if (k[2], k[3]) != (63, 127)) == 10


def test_T7_complete_set_validates_and_feeds_T3_unchanged():
    recs = _full_set()
    assert len(rr.validate_record_set(recs, ["A"])) == 145
    g = rr.three_grid_values(recs, config="A", method="M1")
    assert g["points"] == rr.POINT_ORDER and set(g["values"]) == {(p, q) for p in rr.POINT_ORDER for q in ("T1", "R1")}
    sel = rgf.select_resolution(g["values"], g["points"], solver_agreement=True)
    flags = rgf.point_flags(g["values"], g["residuals"], g["points"])
    assert sel.outcome == rgf.SELECTED and sel.dx_selected == "1/512" and not any(f.flagged for f in flags.values())
    g2 = rr.three_grid_values(recs, config="A", method="M2")
    assert g2["points"] == rr.POINTS_14


def test_T7_three_grid_values_never_mixes_methods_and_refuses_M4():
    key = ("A", "L1", 63, 127, 1024, "M2")
    recs = _full_set(overrides={key: _S0() * np.array([[1, 0.5], [0.5, 1]])})
    g1, g2 = (rr.three_grid_values(recs, config="A", method=m) for m in ("M1", "M2"))
    assert g1["values"][("L1", "T1")] == (1.0, 1.0, 1.0) and g2["values"][("L1", "T1")][1] == 0.25
    with pytest.raises(rr.PilotRecordError, match="M4"):
        rr.three_grid_values(recs, config="A", method="M4")


@pytest.mark.parametrize("method", ["M1", "M2", "M3"])
def test_T7_one_missing_method_record_at_one_point_is_refused(method):
    recs = _full_set()
    drop = ("A", "T10+", 63, 127, 1024, method)
    recs = [r for r in recs if r.key != drop]
    with pytest.raises(rr.PilotRecordError, match=r"missing \[\('A', 'T10\+', 63, 127, 1024, '" + method):
        rr.agreement_rows(recs, configs=["A"])
    with pytest.raises(rr.PilotRecordError, match="missing"):
        rr.three_grid_values(recs, config="A", method="M1")


@pytest.mark.parametrize("drop", [("A", "B3-E2", 63, 127, 2048, "M1"), ("A", "U3", 127, 254, 1024, "M1"),
                                  ("A", "T11−", 63, 254, 1024, "M1")])
def test_T7_missing_B3_or_transverse_record_is_refused(drop):
    recs = [r for r in _full_set() if r.key != drop]
    with pytest.raises(rr.PilotRecordError, match="missing"):
        rr.validate_record_set(recs, ["A"])


def test_T7_duplicate_foreign_or_empty_record_sets_are_refused():
    recs = _full_set()
    with pytest.raises(rr.PilotRecordError, match="duplicate"):
        rr.validate_record_set(recs + [recs[0]], ["A"])
    with pytest.raises(rr.PilotRecordError, match="unlisted configs"):
        rr.validate_record_set(recs + [_rec(config="B")], ["A"])
    with pytest.raises(rr.PilotRecordError, match="missing"):
        rr.validate_record_set([], ["A"])
    with pytest.raises(rr.PilotRecordError, match="configs is empty"):
        rr.validate_record_set(recs, [])


@pytest.mark.parametrize("field,value", [("E", 10.000000000000002), ("x_start", -0.5000000000000001),
                                         ("potential_sha256", "sha256:" + "f" * 64), ("n_slices", 720)])
def test_T7_check_same_problem_detects_mismatches(field, value):
    a = _rec(method="M1")
    b = dataclasses.replace(_rec(method="M2"), **{field: value})
    with pytest.raises(rr.PilotRecordError, match=field):
        rr.check_same_problem([a, b])


def test_T7_check_same_problem_detects_n_open_and_shape_and_transverse_mismatch():
    a = _rec(point="U1", inv=1024)
    with pytest.raises(rr.PilotRecordError, match="n_open"):
        rr.check_same_problem([a, _rec(point="U1", inv=1024, method="M2", S=np.eye(4, dtype=np.complex128))])
    with pytest.raises(rr.PilotRecordError, match="n_y|n_q|potential"):
        rr.check_same_problem([a, _rec(point="U1", inv=1024, tv=(63, 254))])
    with pytest.raises(rr.PilotRecordError, match="at least two"):
        rr.check_same_problem([a])


def test_T7_agreement_rows_reject_an_inconsistent_problem_inside_a_complete_set():
    recs = _full_set()
    i = next(n for n, r in enumerate(recs) if r.key == ("A", "U2", 63, 127, 2048, "M3"))
    recs[i] = dataclasses.replace(recs[i], E=recs[i].E + 1.0)
    with pytest.raises(rr.PilotRecordError, match="differ in E"):
        rr.agreement_rows(recs, configs=["A"])


# ---- raw agreement rows and the M4 trigger rule ----------------------------------------------------------------
def test_T7_agreement_rows_complete_and_within_limit():
    rows = rr.agreement_rows(_full_set(), configs=["A"])
    assert len(rows) == 14 * 3 * 3 and all(r["dS"] == 0.0 and not r["exceeds"] for r in rows)
    s = rr.agreement_summary(rows, configs=["A"])
    assert s.all_within_limit and s.exceedances == () and s.protocol_gap_keys == () and s.n_rows == 126
    assert rr.m4_required_keys(rows, configs=["A"]) == tuple(("A", p, 2048) for p in rr.M4_PROTOCOL_POINTS)
    names = {f.name for f in dataclasses.fields(rr.AgreementSummary)}
    assert not any(t in n for n in names for t in ("step", "solver_agreement", "passed", "gate"))


@pytest.mark.parametrize("method,trig_pairs", [("M2", {"M1-M2"}), ("M3", {"M1-M3"})])
def test_T7_M1_M2_or_M1_M3_exceedance_triggers_M4(method, trig_pairs):
    key = ("A", "L2", 63, 127, 1024, method)
    rows = rr.agreement_rows(_full_set(overrides={key: _bumped(2.0 ** -32)}), configs=["A"])
    hit = [r for r in rows if r["exceeds"]]
    assert {r["pair"] for r in hit if r["m4_trigger"]} == trig_pairs
    assert {r["pair"] for r in hit} == trig_pairs | {"M2-M3"} and all(r["point"] == "L2" for r in hit)
    req = rr.m4_required_keys(rows, configs=["A"])
    assert ("A", "L2", 1024) in req and len(req) == 6
    s = rr.agreement_summary(rows, configs=["A"])
    assert not s.all_within_limit and s.m4_triggered_keys == (("A", "L2", 1024),) and s.protocol_gap_keys == ()


def test_T7_M2_M3_only_exceedance_is_a_protocol_gap_not_an_M4_trigger():
    d = 2.0 ** -34                                                     # 5.8e-11 each way: M1-M2, M1-M3 <= 1e-10
    ov = {("A", "U1", 63, 127, 512, "M2"): _bumped(d), ("A", "U1", 63, 127, 512, "M3"): _bumped(-d)}
    rows = rr.agreement_rows(_full_set(overrides=ov), configs=["A"])
    hit = [r for r in rows if r["exceeds"]]
    assert [(r["pair"], r["m4_trigger"]) for r in hit] == [("M2-M3", False)]
    s = rr.agreement_summary(rows, configs=["A"])
    assert not s.all_within_limit and s.m4_triggered_keys == () and s.protocol_gap_keys == (("A", "U1", 512),)
    assert ("A", "U1", 512) not in rr.m4_required_keys(rows, configs=["A"])


def test_T7_agreement_limit_is_inclusive_at_exactly_1e_10():
    base = np.array([[_A, 0], [0, _A]], dtype=np.complex128)
    at, over = base.copy(), base.copy()
    at[0, 1], over[0, 1] = 1e-10, np.nextafter(1e-10, 1.0)
    for S2, expect in ((at, False), (over, True)):
        ov = {k: base for k in rr.expected_record_keys(["A"])}
        ov[("A", "C1", 63, 127, 512, "M2")] = S2
        rows = rr.agreement_rows(_full_set(overrides=ov), configs=["A"])
        r = next(x for x in rows if (x["point"], x["inv_dx"], x["pair"]) == ("C1", 512, "M1-M2"))
        assert r["exceeds"] is expect and r["m4_trigger"] is expect


def test_T7_incomplete_or_tampered_rows_are_refused():
    rows = rr.agreement_rows(_full_set(), configs=["A"])
    for bad in ([], rows[1:], rows + [rows[0]]):
        with pytest.raises(rr.PilotRecordError):
            rr.agreement_summary(bad, configs=["A"])
    with pytest.raises(rr.PilotRecordError, match="configs is empty"):
        rr.agreement_summary(rows, configs=[])
    with pytest.raises(rr.PilotRecordError, match="missing"):
        rr.agreement_summary(rows, configs=["A", "B"])                 # B rows absent: never all_within_limit
    t = [dict(r) for r in rows]
    t[5]["exceeds"] = True
    with pytest.raises(rr.PilotRecordError, match="exceeds"):
        rr.agreement_summary(t, configs=["A"])
    t = [dict(r) for r in rows]
    t[0]["dS"], t[0]["exceeds"], t[0]["m4_trigger"] = 1.0, True, False
    with pytest.raises(rr.PilotRecordError, match="m4_trigger"):
        rr.m4_required_keys(t, configs=["A"])


# ---- M4 adjudication inputs -------------------------------------------------------------------------------------
_M4_KEYS = tuple(("A", p, 2048) for p in rr.M4_PROTOCOL_POINTS)


def _m4_set(**kw):
    return [_rec(point=p, inv=2048, method="M4", **kw) for p in rr.M4_PROTOCOL_POINTS]


def test_T7_M4_inputs_exact_eps_and_statuses():
    S4 = _S0().astype(_CLD)
    S4[0, 0] += _CLD(2.0) ** -40
    m4 = [_rec(point=p, inv=2048, method="M4") for p in rr.M4_PROTOCOL_POINTS[1:]]
    m4.insert(0, rr.extract_m4_record(config="A", point="L1", inv_dx=2048, n_y=63, n_q=127, E=10.0, kx=1.25,
                                      solver_out=dict(S=S4, n_open=1, error_indicator=1e-13,
                                                      refinement_gate={"passed": True},
                                                      green_refinement_gate={"passed": True},
                                                      fl_vs_state_max_abs=0.0, eta_final=1e-20),
                                      V_lead_max=0.0, potential_sha256=_psha("A", 2048, (63, 127)), n_slices=2873,
                                      x_start=-0.5, wall_time_s=0.5, code_files=_CODE_FILES, platform_report=_EVAL))
    out = rr.m4_adjudication_inputs(_full_set() + m4, _EVAL, required_keys=_M4_KEYS)
    assert [o["point"] for o in out] == list(rr.M4_PROTOCOL_POINTS)
    assert out[0]["eps_M1"] == out[0]["eps_M2"] == out[0]["eps_M3"] == 2.0 ** -40
    assert all(o["m4_status"] == rr.M4_EVALUATED and o["eps_M1"] == 0.0 for o in out[1:])


def test_T7_M4_gate_failure_and_platform_not_evaluable_are_never_a_pass():
    bad = [_rec(point=p, inv=2048, method="M4", error_indicator=1e-9) for p in rr.M4_PROTOCOL_POINTS]
    assert not bad[0].extras["X7_passed"] and not bad[0].extras["gates_passed"]
    out = rr.m4_adjudication_inputs(_full_set() + bad, _EVAL, required_keys=_M4_KEYS)
    assert {o["m4_status"] for o in out} == {rr.M4_NE_GATE}
    out = rr.m4_adjudication_inputs(_full_set(), _NOT_EVAL, required_keys=_M4_KEYS)
    assert {o["m4_status"] for o in out} == {rr.M4_NE_PLATFORM} and all(o["eps_M1"] is None for o in out)
    with pytest.raises(rr.PilotRecordError, match="NOT EVALUABLE"):
        _rec(point="L1", inv=2048, method="M4", platform=_NOT_EVAL)
    with pytest.raises(rr.PilotRecordError, match="not evaluable"):
        rr.m4_adjudication_inputs(_full_set() + _m4_set(), _NOT_EVAL, required_keys=_M4_KEYS)


def test_T7_M4_set_must_equal_the_required_keys():
    full = _full_set()
    with pytest.raises(rr.PilotRecordError, match="missing"):
        rr.m4_adjudication_inputs(full + _m4_set()[1:], _EVAL, required_keys=_M4_KEYS)
    with pytest.raises(rr.PilotRecordError, match="unexpected"):
        rr.m4_adjudication_inputs(full + _m4_set() + [_rec(point="L2", inv=512, method="M4")], _EVAL,
                                  required_keys=_M4_KEYS)
    for bad in ((), _M4_KEYS + _M4_KEYS[:1], (("A", "B3-E1", 2048),), (("A", "L1", 4096),)):
        with pytest.raises(rr.PilotRecordError):
            rr.m4_adjudication_inputs(full + _m4_set(), _EVAL, required_keys=bad)


def test_T7_M4_must_describe_the_same_problem_and_have_grid_partners():
    m4 = _m4_set()
    m4[2] = dataclasses.replace(m4[2], E=m4[2].E + 1.0)
    with pytest.raises(rr.PilotRecordError, match="differ in E"):
        rr.m4_adjudication_inputs(_full_set() + m4, _EVAL, required_keys=_M4_KEYS)
    two = [_rec(point="T3+", inv=2048, method="M4", S=np.eye(4, dtype=np.complex128))]
    with pytest.raises(rr.PilotRecordError, match="n_open"):
        rr.m4_adjudication_inputs(_full_set() + two, _EVAL, required_keys=(("A", "T3+", 2048),))
    partial = [r for r in _full_set() if r.key != ("A", "U3", 63, 127, 2048, "M2")]
    with pytest.raises(rr.PilotRecordError, match="all required"):
        rr.m4_adjudication_inputs(partial + _m4_set(), _EVAL, required_keys=_M4_KEYS)


# ---- M4 exact storage -------------------------------------------------------------------------------------------
def test_T7_clongdouble_hi_lo_split_is_exact_and_refuses_lossy_values():
    S = _S0().astype(_CLD)
    S[0, 1] += _CLD(2.0) ** -60                                        # below float64 resolution of 1/2
    hi, lo = rr.split_hi_lo(S)
    assert hi.dtype == lo.dtype == np.complex128 and np.array_equal(hi.astype(_CLD) + lo.astype(_CLD), S)
    assert lo[0, 1] != 0
    U = _S0().astype(_CLD)
    U[1, 1] = _CLD(2.0) ** -1074 * _CLD(1.5)                           # lo would underflow in float64
    with pytest.raises(rr.PilotRecordError, match="lossy"):
        rr.split_hi_lo(U)


# ---- fixed-width strings ---------------------------------------------------------------------------------------
def test_T7_every_string_field_holds_its_longest_allowed_value():
    allowed = {"record_config": rr.CONFIGS, "record_point": rr.POINT_ORDER, "record_method": rr.METHODS,
               "record_storage": tuple(rr.STORAGE.values()), "record_s_sha256": ("sha256:" + "a" * 64,),
               "record_method_fingerprint": ("sha256:" + "b" * 64,), "schema_version": (rr.SCHEMA_VERSION,)}
    assert set(allowed) == set(rr.STRING_WIDTHS)
    assert len("clongdouble_hi_lo") == 17 <= rr.STRING_WIDTHS["record_storage"]
    for name, vals in allowed.items():
        arr = rr.fixed_str_array(name, list(vals))
        assert arr.tolist() == list(vals) and arr.dtype == np.dtype(f"<U{rr.STRING_WIDTHS[name]}")
        assert max(len(v) for v in vals) <= rr.STRING_WIDTHS[name]
        with pytest.raises(rr.PilotRecordError, match="exceeds the declared width"):
            rr.fixed_str_array(name, ["x" * (rr.STRING_WIDTHS[name] + 1)])


# ---- npz round trip through the T5 staging interface --------------------------------------------------------
def _write_npz(tmp_repo, tmp_path, recs):
    st = pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path))
    p = st.write_npz("rgf_resolution/diagnostic_rgf_S.npz", **rr.s_archive_arrays(recs))
    return st, p


def test_T7_S_archive_round_trip_is_bitwise(tmp_repo, tmp_path):
    S4 = _S0().astype(_CLD)
    S4[0, 1] += _CLD(2.0) ** -60
    m4 = rr.extract_m4_record(config="A", point="L1", inv_dx=2048, n_y=63, n_q=127, E=10.0, kx=1.25,
                              solver_out=dict(S=S4, n_open=1, error_indicator=1e-13, refinement_gate={"passed": True},
                                              green_refinement_gate={"passed": True}, fl_vs_state_max_abs=0.0,
                                              eta_final=1e-20),
                              V_lead_max=0.0, potential_sha256=_psha("A", 2048, (63, 127)), n_slices=2873,
                              x_start=-0.5, wall_time_s=0.5, code_files=_CODE_FILES, platform_report=_EVAL)
    rng = np.random.default_rng(11)
    S2 = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / 4
    recs = [_rec(method="M3", S=S2), m4, _rec(point="B3-E3", inv=2048), _rec(point="U1", tv=(127, 254), inv=1024)]
    st, p = _write_npz(tmp_repo, tmp_path, recs)
    try:
        ent = rr.load_s_archive(p, _META)
        by = {(e["point"], e["method"]): e for e in ent}
        assert by[("L1", "M3")]["S"].tobytes() == S2.tobytes()
        assert by[("L1", "M4")]["S"].dtype == _CLD and np.array_equal(by[("L1", "M4")]["S"], S4)
        assert [rr._entry_key(e) for e in ent] == [r.key for r in rr.sort_records(recs)]
        assert all(e["s_sha256"] == r.s_sha256 for e, r in zip(ent, rr.sort_records(recs)))
    finally:
        st.discard()


def _tamper(src, dst, fn):
    with np.load(src, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}
    fn(d)
    np.savez(dst, **d)
    return dst


@pytest.mark.parametrize("mutation", ["flip_S_bit", "drop_array", "extra_array", "schema", "classification",
                                      "wrong_int_dtype", "unsorted", "storage", "fingerprint"])
def test_T7_malformed_S_archives_are_refused(tmp_repo, tmp_path, mutation):
    recs = [_rec(), _rec(method="M2"), _rec(point="L2")]
    st, p = _write_npz(tmp_repo, tmp_path, recs)

    def fn(d):
        if mutation == "flip_S_bit":
            d["r0001_S"] = d["r0001_S"].copy()
            d["r0001_S"][0, 0] += 2.0 ** -52
        elif mutation == "drop_array":
            del d["r0002_S"]
        elif mutation == "extra_array":
            d["r0099_S"] = d["r0000_S"]
        elif mutation == "schema":
            d["schema_version"] = np.array("t7-v0", dtype="<U8")
        elif mutation == "classification":
            d["classification"] = np.array("B4")
        elif mutation == "wrong_int_dtype":
            d["record_inv_dx"] = d["record_inv_dx"].astype(np.int32)
        elif mutation == "unsorted":
            for k in [k for k in d if k.startswith("record_")]:
                d[k] = d[k][[1, 0, 2]]
            d["r0000_S"], d["r0001_S"] = d["r0001_S"], d["r0000_S"]
        elif mutation == "storage":
            d["record_storage"] = rr.fixed_str_array("record_storage", ["clongdouble_hi_lo"] * 3)
        else:
            d["record_method_fingerprint"] = rr.fixed_str_array("record_method_fingerprint",
                                                                ["sha256:" + "0" * 64] * 3)
    try:
        bad = _tamper(p, str(tmp_path / "diagnostic_bad.npz"), fn)
        with pytest.raises(rr.PilotRecordError):
            rr.load_s_archive(bad, _META)
    finally:
        st.discard()


# ---- CSV round trip ------------------------------------------------------------------------------------------
def _write_csv(tmp_repo, tmp_path, recs):
    st = pio.Staging(str(tmp_repo), scratch_root=_scratch(tmp_path))
    cols, rows = rr.points_csv_rows(recs)
    return st, st.write_csv("rgf_resolution/diagnostic_rgf_points.csv", cols, rows)


def _mixed_records():
    rng = np.random.default_rng(5)
    S2 = (rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4))) / 4
    return [_rec(method="M1", S=S2, x_start=-0.0), _rec(method="M2", S=S2, x_start=-0.0),
            _rec(method="M3", S=S2, x_start=-0.0), _rec(point="L1", inv=2048, method="M4")]


def test_T7_points_csv_round_trip_is_lossless(tmp_repo, tmp_path):
    recs = _mixed_records()
    st, p = _write_csv(tmp_repo, tmp_path, recs)
    try:
        ent = rr.parse_points_csv(p, _META)
        for e, r in zip(ent, rr.sort_records(recs)):
            for col in rr.CSV_COLUMNS:
                if col in rr.ALL_EXTRAS:
                    assert (col in e) == (col in r.extras) and (col not in e or e[col] == r.extras[col])
                    continue
                v = getattr(r, col)
                assert e[col] == v and (not isinstance(v, float) or np.float64(e[col]).tobytes()
                                        == np.float64(v).tobytes()), col
        assert math.copysign(1.0, ent[0]["x_start"]) == -1.0           # -0.0 preserved
        lines = open(p, encoding="utf-8").read().splitlines()
        assert lines[0] == pio.CSV_HEADER and lines[2].count(rr.NOT_APPLICABLE) == len(rr.ALL_EXTRAS) - 2
    finally:
        st.discard()


@pytest.mark.parametrize("mutation", ["label", "header", "noncanonical_float", "nan", "not_applicable",
                                      "bool", "order", "fingerprint", "short_row"])
def test_T7_malformed_points_csv_is_refused(tmp_repo, tmp_path, mutation):
    st, p = _write_csv(tmp_repo, tmp_path, _mixed_records())
    lines = open(p, encoding="utf-8").read().splitlines()
    cols = lines[1].split(",")
    row = lines[2].split(",")
    if mutation == "label":
        lines[0] = "# something else"
    elif mutation == "header":
        lines[1] = lines[1].replace("T1", "T_1", 1)
    elif mutation in ("noncanonical_float", "nan"):
        row[cols.index("T1")] = "nan" if mutation == "nan" else row[cols.index("T1")] + "0"
    elif mutation == "not_applicable":
        row[cols.index("eta_final")] = "0.0"
    elif mutation == "bool":
        lines[5] = lines[5].replace("True", "true", 1)
    elif mutation == "order":
        lines[2], lines[3] = lines[3], lines[2]
    elif mutation == "fingerprint":
        row[cols.index("method_fingerprint")] = "sha256:" + "1" * 64
    else:
        row = row[:-1]
    if mutation not in ("label", "header", "bool", "order"):
        lines[2] = ",".join(row)
    bad = tmp_path / "diagnostic_bad.csv"
    bad.write_text("\n".join(lines) + "\n", encoding="utf-8")
    try:
        with pytest.raises(rr.PilotRecordError):
            rr.parse_points_csv(str(bad), _META)
    finally:
        st.discard()


def test_T7_archive_and_csv_are_cross_validated(tmp_repo, tmp_path):
    recs = _mixed_records()
    st1, pn = _write_npz(tmp_repo, tmp_path, recs)
    st2, pc = _write_csv(tmp_repo, tmp_path, recs)
    try:
        a, c = rr.load_s_archive(pn, _META), rr.parse_points_csv(pc, _META)
        rr.verify_archive_csv_consistency(a, c)
        c2 = [dict(e) for e in c]
        c2[1]["T1"] = c2[1]["T1"] + 2.0 ** -40
        with pytest.raises(rr.PilotRecordError, match="recomputed"):
            rr.verify_archive_csv_consistency(a, c2)
        with pytest.raises(rr.PilotRecordError, match="keys"):
            rr.verify_archive_csv_consistency(a, c[:-1])
        c3 = [dict(e) for e in c]
        c3[0]["s_sha256"] = "sha256:" + "2" * 64
        with pytest.raises(rr.PilotRecordError, match="s_sha256"):
            rr.verify_archive_csv_consistency(a, c3)
    finally:
        st1.discard()
        st2.discard()


# ---- method fingerprints ------------------------------------------------------------------------------------------
def test_T7_method_fingerprints_are_derived_per_method():
    fps = {m: rr.method_fingerprint(m, _CODE_FILES) for m in rr.METHODS}
    assert len(set(fps.values())) == 4 and fps == {m: rr.method_fingerprint(m, dict(_CODE_FILES)) for m in rr.METHODS}
    changed = dict(_CODE_FILES, **{"stationary.py": "sha256:" + "3" * 64})
    after = {m: rr.method_fingerprint(m, changed) for m in rr.METHODS}
    assert [m for m in rr.METHODS if after[m] != fps[m]] == ["M1"]
    shared = dict(_CODE_FILES, **{"smatrix_spec.py": "sha256:" + "4" * 64})
    assert [m for m in rr.METHODS if rr.method_fingerprint(m, shared) != fps[m]] == ["M1", "M2", "M3", "M4"]
    assert "smatrix_spec.py" in rr.EXTRACTION_CODE_FILES and all("smatrix_spec.py" in rr.method_files(m)
                                                                  for m in rr.METHODS)
    assert _rec(method="M2").method_fingerprint == fps["M2"]


def test_T7_fingerprint_mismatches_are_refused(tmp_repo, tmp_path):
    r2 = _rec(method="M2")
    wrong = dataclasses.replace(r2, method_fingerprint=_rec(method="M1").method_fingerprint)
    with pytest.raises(rr.PilotRecordError, match="fingerprint mismatch"):
        rr.verify_method_fingerprints([wrong], _META)
    missing = {"provenance": {"code_fingerprint": {"code_files": {k: v for k, v in _CODE_FILES.items()
                                                                 if k != "sparse_wavefunction.py"}}}}
    with pytest.raises(rr.PilotRecordError, match="sparse_wavefunction.py"):
        rr.verify_method_fingerprints([r2], missing)
    with pytest.raises(rr.PilotRecordError, match="lacks"):
        rr.verify_method_fingerprints([r2], {})
    st, p = _write_npz(tmp_repo, tmp_path, [_rec(), r2])
    try:
        altered = {"provenance": {"code_fingerprint": {"code_files": dict(
            _CODE_FILES, **{"stationary.py": "sha256:" + "5" * 64})}}}
        with pytest.raises(rr.PilotRecordError, match="M1"):
            rr.load_s_archive(p, altered)
    finally:
        st.discard()


def test_T7_real_code_fingerprint_covers_every_method_file():
    assert set(_ALL_CODE_FILES) <= set(pio.CODE_FILES) and "b4_pilot_rgf_records.py" in pio.CODE_FILES
    fp = pio.code_fingerprint(REPO)
    assert {m: rr.method_fingerprint(m, fp["code_files"]) for m in rr.METHODS}


# ---- determinism ---------------------------------------------------------------------------------------------
def test_T7_serialization_is_deterministic_under_input_order(tmp_repo, tmp_path):
    recs = _mixed_records() + [_rec(point="C2", inv=1024), _rec(config="B", point="T3+")]
    rev = list(reversed(recs))
    assert rr.archive_content_digest(rr.s_archive_arrays(recs)) == rr.archive_content_digest(rr.s_archive_arrays(rev))
    assert rr.points_csv_rows(recs) == rr.points_csv_rows(rev)
    st1, p1 = _write_csv(tmp_repo, tmp_path, recs)
    st2, p2 = _write_csv(tmp_repo, tmp_path, rev)
    try:
        assert open(p1, "rb").read() == open(p2, "rb").read()
    finally:
        st1.discard()
        st2.discard()


# ---- module purity (read-only file access allowed; no writes, solvers, TD or energy averages) --------------
def test_T7_module_is_read_only_and_solver_free():
    import ast
    src = open(os.path.join(REPO, "b4_pilot_rgf_records.py"), encoding="utf-8").read()
    tree = ast.parse(src)
    imported = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    imported |= {n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module}
    assert imported <= {"__future__", "csv", "hashlib", "math", "re", "types", "collections", "dataclasses", "numpy",
                        "smatrix_spec", "sparse_exact_lattice", "validate_rgf_b3", "b4_pilot_io"}, imported
    from_names = {a.name for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) for a in n.names}
    assert from_names <= {"annotations", "Mapping", "dataclass", "field", "CLASSIFICATION", "CSV_HEADER"}, from_names
    forbidden_attr = {"save", "savez", "savez_compressed", "savetxt", "tofile", "remove", "unlink", "rename",
                      "replace", "mkdir", "makedirs", "rmdir", "rmtree", "write_text", "write_bytes", "system",
                      "popen", "scattering_matrix", "sparse_scattering_matrix", "wavefunction_scattering_matrix",
                      "exact_lattice_S", "exact_lattice_S_from_blocks", "build_device", "make_rgf_device",
                      "run_case", "energy_average", "write_csv", "write_npz", "write_json", "publish"}
    for n in ast.walk(tree):
        if isinstance(n, ast.Call):
            f = n.func
            name = f.attr if isinstance(f, ast.Attribute) else getattr(f, "id", None)
            assert name not in forbidden_attr, name
            if name == "open":
                mode = n.args[1] if len(n.args) > 1 else next((k.value for k in n.keywords if k.arg == "mode"), None)
                assert isinstance(mode, ast.Constant) and not set(mode.value) & set("wax+"), ast.dump(n)
            if isinstance(f, ast.Attribute) and f.attr == "load":
                kw = {k.arg: k.value for k in n.keywords}
                assert isinstance(kw.get("allow_pickle"), ast.Constant) and kw["allow_pickle"].value is False
    for s in [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)]:
        assert not re.search(r"results/(rgf_b|stationary_validation|transport|validation|b4_pilot)", s), s
    assert "energy_average" not in src and "split_step" not in src and "td_" not in src.lower().replace("std_", "")


# ---- T7 corrections: immutability, strict gate Booleans, row and CSV validation, strict M4 input types -------
def test_T7_extras_are_immutable_and_isolated_from_the_caller():
    failed = _rec(point="L1", inv=2048, method="M4", error_indicator=1e-9)          # X7 fails -> gates_passed False
    assert failed.extras["gates_passed"] is False
    with pytest.raises(TypeError):
        failed.extras["gates_passed"] = True
    with pytest.raises((TypeError, AttributeError)):
        failed.extras.clear()
    with pytest.raises(dataclasses.FrozenInstanceError):
        failed.extras = {}
    src = dict(min_rcond=0.125, argmin_rcond=7)
    r = dataclasses.replace(_rec(), extras=src)
    src["min_rcond"] = 99.0                                                       # caller's dict changed later
    assert r.extras["min_rcond"] == 0.125
    r2 = dataclasses.replace(r, wall_time_s=1.0)                                  # replace still works
    assert r2.wall_time_s == 1.0 and dict(r2.extras) == dict(r.extras) and type(r2.extras) is type(r.extras)
    m4 = [failed] + [_rec(point=p, inv=2048, method="M4", error_indicator=1e-9) for p in rr.M4_PROTOCOL_POINTS[1:]]
    out = rr.m4_adjudication_inputs(_full_set() + m4, _EVAL, required_keys=_M4_KEYS)
    assert {o["m4_status"] for o in out} == {rr.M4_NE_GATE}


@pytest.mark.parametrize("method", ["M1", "M4"])
def test_T7_S_cannot_be_made_writeable(method):
    S = _S0()
    r = _rec(method=method, inv=2048, S=S)
    before = r.s_sha256
    with pytest.raises(ValueError):
        r.S.setflags(write=True)
    with pytest.raises(ValueError):
        r.S[0, 0] = 5
    expect = S.astype(_CLD) if method == "M4" else S
    assert r.S.dtype == expect.dtype and r.S.tobytes() == expect.tobytes() and r.s_sha256 == before


@pytest.mark.parametrize("gate", ["refinement_gate", "green_refinement_gate"])
@pytest.mark.parametrize("value", ["False", "True", 0, 1, None, 1.0, {"passed": "False"}, {"passed": 1},
                                   {"passed": None}, "GATE_OBJECT_NONE", True, {}, [True], ("passed", True)],
                         ids=lambda v: repr(v))
def test_T7_malformed_M4_gate_values_never_produce_a_record(gate, value):
    obj = value if not isinstance(value, (str, int, float)) or isinstance(value, bool) else {"passed": value}
    if value is None:
        obj = {"passed": None}                                      # a None 'passed' value
    if value == "GATE_OBJECT_NONE":
        obj = None                                                  # the gate object itself is None
    with pytest.raises(rr.PilotRecordError):
        _rec(point="L1", inv=2048, method="M4", **{gate: obj})


def test_T7_numpy_bool_gate_values_are_accepted():
    r = _rec(point="L1", inv=2048, method="M4", refinement_gate={"passed": np.bool_(True)},
             green_refinement_gate={"passed": np.bool_(True)})
    assert r.extras["X4_state_passed"] is True and r.extras["gates_passed"] is True


@pytest.mark.parametrize("field,value", [("limit", 1e-9), ("limit", 1e-10 * (1 + 2.0 ** -52)), ("limit", True),
                                         ("limit", "1e-10"), ("limit", None), ("classification", "B4"),
                                         ("classification", None), ("limit", "DROP"),
                                         ("classification", "DROP")])
def test_T7_agreement_rows_with_altered_limit_or_classification_are_refused(field, value):
    rows = [dict(r) for r in rr.agreement_rows(_full_set(), configs=["A"])]
    if value == "DROP":
        del rows[7][field]
    else:
        rows[7][field] = value
    with pytest.raises(rr.PilotRecordError):
        rr.agreement_summary(rows, configs=["A"])
    with pytest.raises(rr.PilotRecordError):
        rr.m4_required_keys(rows, configs=["A"])


@pytest.mark.parametrize("col,value,msg", [("n_open", "0", "n_open"), ("n_open", "63", "n_open"),
                                           ("n_closed", "61", "n_closed"), ("n_slices", "0", "n_slices"),
                                           ("wall_time_s", "-0.5", "wall_time_s")])
def test_T7_points_csv_structural_violations_are_refused(tmp_repo, tmp_path, col, value, msg):
    st, p = _write_csv(tmp_repo, tmp_path, [_rec()])
    try:
        lines = open(p, encoding="utf-8").read().splitlines()
        cols = lines[1].split(",")
        row = lines[2].split(",")
        row[cols.index(col)] = value
        lines[2] = ",".join(row)
        bad = tmp_path / "diagnostic_bad.csv"
        bad.write_text("\n".join(lines) + "\n", encoding="utf-8")
        assert rr.parse_points_csv(p, _META)[0][col] != (float(value) if col == "wall_time_s" else int(value))
        with pytest.raises(rr.PilotRecordError, match=msg):
            rr.parse_points_csv(str(bad), _META)
    finally:
        st.discard()


@pytest.mark.parametrize("platform", [_EVAL, _NOT_EVAL])
@pytest.mark.parametrize("intruder", [{"method": "M4"}, None, "M4"], ids=["dict", "None", "str"])
def test_T7_m4_adjudication_inputs_refuse_non_record_objects(platform, intruder):
    recs = _full_set() + (_m4_set() if platform["evaluable"] else []) + [intruder]
    with pytest.raises(rr.PilotRecordError, match="RGFRecord objects only"):
        rr.m4_adjudication_inputs(recs, platform, required_keys=_M4_KEYS)
