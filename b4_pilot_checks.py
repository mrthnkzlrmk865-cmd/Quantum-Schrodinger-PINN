"""B4 limited diagnostic pilot — gate G0 and gate G-REAL check logic (phase T6).

Classification: DIAGNOSTIC / NOT B4.

Pure check functions and fail-closed orchestration for docs/B4_PILOT_PROTOCOL.md §6 (G0, G-REAL), §1.2, §1.4, §1.6,
§3.3.1, §3.3.2 and §5.3–§5.5. Every physical dependency (M1 control solves, device construction, TD solver
construction, potential sampling) is INJECTED by the runner modules b4_pilot_run_g0 / b4_pilot_run_greal; this module
imports no solver and performs no file write.

Pass predicates (never weakened):
* G0 is PASS only if every G0 check is PASS (any failure stops the pilot; no partial output).
* G-REAL is PASS only if every blocking L1, L2 and L3 check is PASS. A FAIL or NOT EVALUABLE is "not passed".
  The shifted-grid dx-equality result (§3.3.1) is NON-BLOCKING evidence: recorded with the actual dx values and
  listed in nonblocking_failures; it is never dropped. NOT APPLICABLE marks L3 pairs without common physical points
  (RGF <-> TD, different transverse grids); it is not a claim and is never counted.
* NOT EVALUABLE never counts as PASS.
"""
from __future__ import annotations

import hashlib
import math
import os
import re
import types
from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

import b4_pilot_io as pio
import b4_pilot_rgf_records as rr

CLASSIFICATION = pio.CLASSIFICATION

PASS, FAIL, NOT_EVALUABLE, NOT_APPLICABLE = "PASS", "FAIL", "NOT EVALUABLE", "NOT APPLICABLE"
STATUSES = (PASS, FAIL, NOT_EVALUABLE, NOT_APPLICABLE)
LEVELS = ("G0", "L1", "L2", "L3", "DX")
GATES = ("G0", "G-REAL")


class PilotCheckError(RuntimeError):
    """Refusal of the gate orchestration (fail closed)."""


class NotEvaluable(RuntimeError):
    """Raised by an injected dependency when an item cannot be evaluated (e.g. the lead-tail guard)."""


# ---- protected artifacts (protocol §5.3; the same lists as tests/test_b4_pilot.py T1.5) -------------------------
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
    "docs/B4_IMPLEMENTATION_READINESS_AUDIT.md", pio.PROTOCOL_PATH,
)
PROTECTED_RESULT_DIRS = (
    "results/rgf_b1_validation", "results/rgf_b2_validation", "results/rgf_b3_validation",
    "results/rgf_b3x_validation", "results/rgf_b3y_validation", "results/validation", "results/transport",
)
PROTECTED_RESULT_COUNT = 109


@dataclass(frozen=True)
class ProtectedSpec:
    sources: tuple = PROTECTED_SOURCES
    docs: tuple = PROTECTED_DOCS
    result_dirs: tuple = PROTECTED_RESULT_DIRS
    tests_dir: str | None = "tests"            # every test file tracked at the baseline is protected
    expected_result_count: int | None = PROTECTED_RESULT_COUNT


# ---- control references (§1.4, §1.6; decision: B3X dx ladder primary) --------------------------------------------
CONFIGS = rr.CONFIGS
B3X_NAME = {"A": "W100_z0.02", "B": "W100_z0.05", "C": "W100_z0.2", "D": "W200_z0.05"}
CONTROL_TOL = 1e-12                                            # §0.1 control reproduction, <= (inclusive)
LADDER_POINTS = ("B3-E1", "C1", "B3-E2", "B3-E3", "C2")         # the five B3 energies in ascending E (j = -2..2)
CONTROL_INV = {"C1": (512, 1024), "C2": (512, 1024), "B3-E1": (512, 1024, 2048), "B3-E2": (512, 1024, 2048),
               "B3-E3": (512, 1024, 2048)}
CONTROL_KEYS = tuple((c, p, inv) for c in CONFIGS for p in LADDER_POINTS for inv in CONTROL_INV[p])   # 52
B3X_LADDER = "results/rgf_b3x_validation/dx_ladder_diagnostic.csv"
B3_COMPARISON = "results/rgf_b3_validation/rgf_sparse_comparison.csv"
B3Y_COMPARISON = "results/rgf_b3y_validation/per_energy_comparison.csv"
B3_METADATA = "results/rgf_b3_validation/run_metadata.json"
B3X_METADATA = "results/rgf_b3x_validation/run_metadata.json"

# ---- §3.3.1 variants (cross-checked against the protocol table at run time) ---------------------------------------
@dataclass(frozen=True)
class VariantSpec:
    name: str
    m: int
    cap_width_extends: bool        # width = 0.6 + s (V1, V2, V4) or the production 0.6 (V0, V3)
    cap_strength: float
    cap_power: int = 2


VARIANTS = (VariantSpec("V0", 0, False, 300.0), VariantSpec("V1", 39, True, 300.0),
            VariantSpec("V2", 39, True, 150.0), VariantSpec("V3", 64, False, 300.0),
            VariantSpec("V4", 115, True, 300.0))
PRODUCTION = dict(Lx=3.9, nx=249, ny=63, Ly=1.0, cap_width=0.6, cap_strength=300.0, cap_power=2)

# ---- RGF devices (decision 3: 5 per configuration) ---------------------------------------------------------------
RGF_DEVICES = tuple((inv, 63, 127) for inv in rr.INV_DX) + ((1024, 63, 254), (1024, 127, 254))
RGF_NESTED_PAIRS = ((512, 1024), (512, 2048), (1024, 2048))


# ======================================================================================================================
# immutable result types
# ======================================================================================================================
def _freeze(v):
    if isinstance(v, Mapping):
        return types.MappingProxyType({str(k): _freeze(x) for k, x in v.items()})
    if isinstance(v, (list, tuple)):
        return tuple(_freeze(x) for x in v)
    if isinstance(v, np.ndarray):
        return tuple(_freeze(x) for x in v.tolist())
    if isinstance(v, (np.bool_,)):
        return bool(v)
    if isinstance(v, np.integer):
        return int(v)
    if isinstance(v, np.floating):
        return float(v)
    if v is None or isinstance(v, (bool, int, float, str)):
        return v
    raise PilotCheckError(f"evidence value of type {type(v).__name__} is not allowed")


def _thaw(v):
    if isinstance(v, Mapping):
        return {k: _thaw(x) for k, x in v.items()}
    if isinstance(v, tuple):
        return [_thaw(x) for x in v]
    return v


@dataclass(frozen=True)
class CheckResult:
    check_id: str
    level: str
    status: str
    blocking: bool = True
    evidence: Mapping = field(default_factory=dict)
    message: str = ""
    classification: str = CLASSIFICATION

    def __post_init__(self):
        if self.level not in LEVELS:
            raise PilotCheckError(f"unknown level {self.level!r}")
        if self.status not in STATUSES:
            raise PilotCheckError(f"unknown status {self.status!r}")
        if not isinstance(self.blocking, bool):
            raise PilotCheckError("blocking must be a bool")
        if self.level == "DX" and self.blocking:
            raise PilotCheckError("the shifted-grid dx-equality result is non-blocking evidence (§3.3.1)")
        if self.level != "DX" and self.status != NOT_APPLICABLE and not self.blocking:
            raise PilotCheckError(f"{self.check_id}: only DX evidence and NOT APPLICABLE entries are non-blocking")
        if self.status == NOT_APPLICABLE and (self.level != "L3" or self.blocking):
            raise PilotCheckError("NOT APPLICABLE is allowed only for non-blocking L3 entries")
        if self.classification != CLASSIFICATION:
            raise PilotCheckError("classification must be DIAGNOSTIC / NOT B4")
        object.__setattr__(self, "evidence", _freeze(dict(self.evidence)))

    def to_dict(self):
        return dict(check_id=self.check_id, level=self.level, status=self.status, blocking=self.blocking,
                    evidence=_thaw(self.evidence), message=self.message, classification=self.classification)


def _derived_status(results, not_run):
    blocking = [r for r in results if r.blocking]
    if results and not not_run and all(r.status == PASS for r in blocking):
        return PASS
    if any(r.status == FAIL for r in blocking):
        return FAIL
    return NOT_EVALUABLE


@dataclass(frozen=True)
class GateReport:
    """Gate outcome. status is PASS only if every blocking result is PASS and nothing was left unrun; the given
    status must equal the derived one (an inconsistent report is refused)."""
    gate: str
    status: str
    results: tuple
    not_run: tuple
    head: str
    protocol_sha256: str
    stopped_at: str | None = None
    classification: str = CLASSIFICATION

    def __post_init__(self):
        if self.gate not in GATES:
            raise PilotCheckError(f"unknown gate {self.gate!r}")
        res = tuple(self.results)
        if not all(isinstance(r, CheckResult) for r in res):
            raise PilotCheckError("results must be CheckResult objects")
        object.__setattr__(self, "results", res)
        object.__setattr__(self, "not_run", tuple(self.not_run))
        if self.status not in (PASS, FAIL, NOT_EVALUABLE):
            raise PilotCheckError(f"invalid gate status {self.status!r}")
        derived = _derived_status(res, self.not_run)
        if self.status != derived:
            raise PilotCheckError(f"inconsistent {self.gate} report: status {self.status!r}, derived {derived!r}")
        if self.classification != CLASSIFICATION:
            raise PilotCheckError("classification must be DIAGNOSTIC / NOT B4")

    @property
    def passed(self) -> bool:
        return self.status == PASS

    @property
    def nonblocking_failures(self) -> tuple:
        return tuple(r.check_id for r in self.results if r.level == "DX" and r.status != PASS)

    def result(self, check_id):
        for r in self.results:
            if r.check_id == check_id:
                return r
        raise KeyError(check_id)

    def to_dict(self):
        return dict(gate=self.gate, status=self.status, head=self.head, protocol_sha256=self.protocol_sha256,
                    stopped_at=self.stopped_at, not_run=list(self.not_run),
                    nonblocking_failures=list(self.nonblocking_failures),
                    results=[r.to_dict() for r in self.results], classification=self.classification)


class _Sequence:
    """Fail-closed sequential execution: the first blocking non-PASS stops the gate; the rest is listed not_run."""

    def __init__(self, gate, planned):
        self.gate, self.planned, self.results, self.stopped = gate, list(planned), [], None

    def add(self, res: CheckResult) -> bool:
        self.results.append(res)
        if res.blocking and res.status != PASS:
            self.stopped = res.check_id
            return False
        return True

    def guard(self, check_id, level, fn, *args, **kw) -> CheckResult:
        try:
            return fn(*args, **kw)
        except NotEvaluable as exc:
            return CheckResult(check_id, level, NOT_EVALUABLE, message=str(exc))
        except Exception as exc:                                      # any error is a failure, never a pass
            return CheckResult(check_id, level, FAIL, message=f"{type(exc).__name__}: {exc}")

    def report(self, head, protocol_sha256):
        done = {r.check_id for r in self.results}
        not_run = tuple(c for c in self.planned if c not in done) if self.stopped else ()
        return GateReport(self.gate, _derived_status(tuple(self.results), not_run), tuple(self.results), not_run,
                          head, protocol_sha256, self.stopped)


# ======================================================================================================================
# protocol-table parsing (ports of the T1 test logic; the protocol hash is verified first by G0-1)
# ======================================================================================================================
_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")


def _section(text, start, end):
    i = text.index(start)
    return text[i:text.index(end, i)]


def _num(s):
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


def _expect(cond, msg, problems):
    if not cond:
        problems.append(msg)


def g0_point_table(text, grid) -> CheckResult:
    """G0-4: protocol §1.2 reproduced exactly. grid provides kx_window, e_perp, adjust_energy (validate_rgf_b2)."""
    sec = _section(text, "### 1.2 Predeclared energy points", "### 1.3")
    rows = _rows(sec, r"(L[123]|T(3|10|11)[−+]|U[123]|C[12])")
    problems = []
    _expect(len(rows) == 14, f"expected 14 point rows, found {len(rows)}", problems)
    k = grid.kx_window(20.0, 0.2)
    Ep = grid.e_perp()
    sk = 1.0 / (np.sqrt(2.0) * 0.2)
    _expect(k.size == 1274, f"coarse grid has {k.size} nodes, expected 1274", problems)
    for cells in rows:
        pid, gi, kx_s, E_s, oc, thr, dE, z = cells[:8]
        if gi.startswith("—"):
            kx = 20.0 + {"C1": -1, "C2": 2}[pid] * sk
        else:
            kx = float(k[int(gi)])
        E = 0.5 * (kx ** 2 + np.pi ** 2)
        d = np.abs(E - Ep)
        n = int(np.argmin(d))
        n_open = int(np.sum(Ep < E))
        _expect(f"{kx:.6f}" == kx_s, f"{pid}: kx {kx:.6f} != {kx_s}", problems)
        _expect(f"{E:.5f}" == E_s, f"{pid}: E {E:.5f} != {E_s}", problems)
        _expect(oc == f"{n_open} ({63 - n_open})", f"{pid}: open/closed {oc}", problems)
        _expect(thr == f"{n + 1} ({Ep[n]:.4f})", f"{pid}: threshold {thr}", problems)
        m = re.fullmatch(r"(\S+) \((\S+)\)", dE)
        if m is None:
            problems.append(f"{pid}: malformed E - E_n cell {dE!r}")
        else:
            _expect(f"{E - Ep[n]:+.3f}" == _num(m.group(1)), f"{pid}: E - E_n", problems)
            _expect(f"{d[n] / E:.1e}" == _num(m.group(2)), f"{pid}: relative distance", problems)
        _expect(f"{(kx - 20.0) / sk:+.2f}" == _num(z), f"{pid}: (kx - k0)/sigma_k", problems)
        _expect(grid.adjust_energy(E, Ep)["adjusted"] is False, f"{pid}: moved by adjust_energy", problems)
    return CheckResult("G0-4", "G0", FAIL if problems else PASS, evidence=dict(n_points=len(rows)),
                       message="; ".join(problems))


def _moved(grid, dkx):
    k = grid.kx_window(20.0, 0.2, dkx)
    h = k[1] - k[0]
    Ep = grid.e_perp()
    out = []
    for i, kx in enumerate(k):
        E = 0.5 * (kx ** 2 + np.pi ** 2)
        a = grid.adjust_energy(E, Ep)
        if a["adjusted"]:
            km = float(np.sqrt(2.0 * a["E"] - np.pi ** 2))
            w = h * (0.5 if i in (0, k.size - 1) else 1.0) * float(grid.packet_weights(kx, 20.0, 0.2))
            out.append(dict(i=i, kx=float(kx), km=km, d_h=(km - kx) / h, E0=a["E_original"], E1=a["E"],
                            n=a["nearest_threshold_n"], side="above" if E > a["nearest_threshold_E"] else "below",
                            kept=a["n_open"] == int(np.sum(Ep < E)), w=w))
    return out


def g0_moved_nodes(text, grid) -> CheckResult:
    """G0-5: protocol §4.1 reproduced exactly (coarse and fine tables, shared energies, sides, nesting)."""
    sec = _section(text, "### 4.1 Design-time enumeration", "### 4.2")
    problems = []
    _expect("11 nodes are moved" in sec and "18 nodes are moved" in sec and "**2.44 grid steps**" in sec,
            "§4.1 summary statements missing", problems)
    rows = [r for r in _rows(sec, r"\d+") if len(r) == 9]
    got = _moved(grid, grid.DKX)
    _expect(len(got) == 11 and len(rows) == 11, f"coarse moved nodes {len(got)}, table rows {len(rows)}", problems)
    for r, g in zip(rows, got):
        i, kx, km, dh, ee, n, side, kept, w = r
        _expect(int(i) == g["i"], f"node {i}", problems)
        _expect(f"{g['kx']:.6f}" == kx and f"{g['km']:.6f}" == km, f"node {i}: kx", problems)
        _expect(f"{g['d_h']:+.3f}" == _num(dh), f"node {i}: dkx/h", problems)
        parts = [x.strip() for x in ee.split("→")]
        _expect(len(parts) == 2 and f"{g['E0']:.4f}" == parts[0] and f"{g['E1']:.4f}" == parts[-1],
                f"node {i}: E", problems)
        _expect(int(n) == g["n"] and side == g["side"] and kept == "yes" and g["kept"], f"node {i}: side", problems)
        _expect(f"{g['w']:.1e}" == _num(w), f"node {i}: weight", problems)
    total = sum(g["w"] for g in got)
    _expect(f"{total:.1e}" == "8.1e-03" and "8.1×10⁻³ of the packet weight" in sec, "total weight", problems)
    fine = _moved(grid, grid.DKX / 2)
    _expect(len(fine) == 18 and all(g["kept"] for g in fine), "fine moved nodes", problems)
    _expect(bool(fine) and f"{max(abs(g['d_h']) for g in fine):.2f}" == "2.44", "fine max displacement", problems)
    by_E = {}
    for g in fine:
        by_E.setdefault(round(g["E1"], 10), []).append(g["i"])
    shared = {tuple(v): e for e, v in by_E.items() if len(v) > 1}
    expected = {(int(a), int(b)): e for a, b, e in re.findall(r"\| (\d+), (\d+) \| E' = ([\d.]+) \|", sec)}
    _expect(len(expected) == 3 and set(shared) == set(expected), "shared fine-grid energies", problems)
    for pair, e in expected.items():
        _expect(pair in shared and f"{shared[pair]:.5f}" == e, f"shared energy {pair}", problems)
    k, kf = grid.kx_window(20.0, 0.2), grid.kx_window(20.0, 0.2, grid.DKX / 2)
    _expect(kf.size == 2547 and np.allclose(kf[::2], k, rtol=0, atol=1e-12), "fine grid nesting", problems)
    return CheckResult("G0-5", "G0", FAIL if problems else PASS,
                       evidence=dict(n_coarse=len(got), n_fine=len(fine)), message="; ".join(problems))


def parse_variant_table(text) -> tuple:
    """§3.3.1 table rows -> (name, m, width text, strength, power, nx); cross-checked against VARIANTS."""
    sec = _section(text, "#### 3.3.1 Variants and the shifted frame", "#### 3.3.2")
    rows = _rows(sec, r"V[0-4] \(.*\)")
    if len(rows) != 5:
        raise PilotCheckError(f"§3.3.1 table has {len(rows)} variant rows, expected 5")
    out = []
    for r in rows:
        name = r[0].split()[0]
        out.append(dict(name=name, m=int(r[1]), width=r[3], strength=float(r[4]), power=int(r[5]), nx=int(r[7])))
    dx = PRODUCTION["Lx"] / (PRODUCTION["nx"] + 1)
    for spec, row in zip(VARIANTS, out):
        width = float(re.findall(r"\d+\.\d+|\d+", row["width"])[-1])        # "0.6 + s = 1.2084" or "1.2084"
        expect = PRODUCTION["cap_width"] + spec.m * dx if spec.cap_width_extends else PRODUCTION["cap_width"]
        if (spec.name, spec.m, f"{expect:.4f}", spec.cap_strength, spec.cap_power,
                PRODUCTION["nx"] + 2 * spec.m) != (row["name"], row["m"], f"{width:.4f}", row["strength"],
                                                   row["power"], row["nx"]):
            raise PilotCheckError(f"§3.3.1 row {row} disagrees with the variant specification {spec}")
    return tuple(out)


# ======================================================================================================================
# G0 helpers
# ======================================================================================================================
def g0_protocol(repo) -> CheckResult:
    sha = pio.check_protocol(repo)
    return CheckResult("G0-1", "G0", PASS, evidence=dict(protocol_sha256=sha))


def g0_clean_tree(repo, check_id="G0-2", expect_head=None) -> CheckResult:
    head = pio.require_clean_tree(repo)
    if expect_head is not None and head != expect_head:
        return CheckResult(check_id, "G0", FAIL, evidence=dict(head=head, expected=expect_head),
                           message="HEAD changed during G0")
    return CheckResult(check_id, "G0", PASS, evidence=dict(head=head))


def _git_lines(repo, *args):
    return [p for p in pio._git(repo, *args).decode().splitlines() if p]


def protected_paths(repo, spec: ProtectedSpec, baseline) -> tuple:
    paths = list(spec.sources) + list(spec.docs)
    if spec.tests_dir:
        tests = _git_lines(repo, "ls-tree", "-r", "--name-only", baseline, "--", spec.tests_dir)
        if not tests:
            raise PilotCheckError(f"no tracked test files under {spec.tests_dir} at baseline {baseline}")
        paths += tests
    n_results = 0
    for d in spec.result_dirs:
        tracked = _git_lines(repo, "ls-tree", "-r", "--name-only", baseline, "--", d)
        if not tracked:
            raise PilotCheckError(f"no tracked files under {d} at baseline {baseline}; refusing to omit it")
        paths += tracked
        n_results += len(tracked)
    if spec.expected_result_count is not None and n_results != spec.expected_result_count:
        raise PilotCheckError(f"protected result-file count {n_results} != {spec.expected_result_count}")
    if len(set(paths)) != len(paths):
        raise PilotCheckError("protected path list contains duplicates")
    return tuple(paths)


def g0_baseline(repo, baseline, spec) -> tuple:
    """G0-3a: the baseline commit exists, is an ancestor of HEAD, and the protected list resolves at it."""
    if pio._git(repo, "cat-file", "-t", baseline).decode().strip() != "commit":
        raise PilotCheckError(f"baseline {baseline} is not a commit")
    pio._git(repo, "merge-base", "--is-ancestor", baseline, "HEAD")      # non-zero exit raises PilotGuardError
    paths = protected_paths(repo, spec, baseline)
    base = pio.baseline_hashes(repo, paths, commit=baseline)
    missing = sorted(p for p, v in base.items() if v is None)
    if missing:
        raise PilotCheckError(f"protected paths absent from baseline {baseline}: {missing}")
    res = CheckResult("G0-3a", "G0", PASS, evidence=dict(baseline=baseline, n_paths=len(paths)))
    return res, paths, base


def _extra_files(repo, spec, baseline):
    extra = []
    for d in spec.result_dirs:
        tracked = set(_git_lines(repo, "ls-tree", "-r", "--name-only", baseline, "--", d))
        for root, _, files in os.walk(os.path.join(repo, d)):
            for f in files:
                rel = os.path.relpath(os.path.join(root, f), repo)
                if rel not in tracked:
                    extra.append(rel)
    return sorted(extra)


def snapshot_check(repo, check_id, paths, base, spec, baseline, pre=None) -> tuple:
    """Snapshot of the protected files; must equal the baseline (and `pre`, if given); no extra result files."""
    snap = pio.protected_snapshot(repo, paths)
    problems = pio.compare_snapshots(base, snap)
    extra = _extra_files(repo, spec, baseline)
    if extra:
        problems.append(f"extra files in protected result directories: {extra}")
    if pre is not None:
        problems += [f"changed since S_pre: {p}" for p in pio.compare_snapshots(pre, snap)]
    res = CheckResult(check_id, "G0", FAIL if problems else PASS,
                      evidence=dict(n_paths=len(paths), snapshot=snap), message="; ".join(problems))
    return res, snap


def _f(s, what):
    try:
        v = float(s)
    except (TypeError, ValueError) as exc:
        raise PilotCheckError(f"{what}: not a number {s!r}") from exc
    if not math.isfinite(v):
        raise PilotCheckError(f"{what}: non-finite {s!r}")
    return v


def _same_bits(a, b):
    return np.float64(a).tobytes() == np.float64(b).tobytes()


def control_references(read_csv) -> dict:
    """{(config, point, inv_dx): (E, T1)} from the B3X dx ladder (primary). Every overlapping B3 (dx = 1/512) and
    B3Y value must equal it bitwise in E and T1, otherwise the reference is inconsistent (raise)."""
    ladder = read_csv(B3X_LADDER)
    refs, by_E = {}, {}
    for c in CONFIGS:
        rows = [r for r in ladder if r.get("config") == B3X_NAME[c]]
        if len(rows) != 5 or any(r.get("label") != "diagnostic" for r in rows):
            raise PilotCheckError(f"B3X ladder: expected 5 diagnostic rows for {B3X_NAME[c]}, got {len(rows)}")
        rows = sorted(rows, key=lambda r: _f(r["E"], "B3X E"))
        for p, r in zip(LADDER_POINTS, rows):
            E = _f(r["E"], "B3X E")
            by_E[(c, E)] = (p, r)
            for inv in CONTROL_INV[p]:
                refs[(c, p, inv)] = (E, _f(r[f"T1_dx_1_{inv}"], f"B3X T1 {c} {p} 1/{inv}"))
    b3 = read_csv(B3_COMPARISON)
    for c in CONFIGS:
        rows = [r for r in b3 if r.get("case") == c]
        if len(rows) != 5:
            raise PilotCheckError(f"B3: expected 5 rows for case {c}, got {len(rows)}")
        for r in rows:
            E = _f(r["E"], "B3 E")
            if (c, E) not in by_E:
                raise PilotCheckError(f"B3 energy {r['E']} of case {c} is not bitwise a B3X ladder energy")
            p, lr = by_E[(c, E)]
            if not _same_bits(_f(r["T1_RGF"], "B3 T1"), _f(lr["T1_dx_1_512"], "B3X T1")):
                raise PilotCheckError(f"B3 and B3X T1 differ for {c} {p} at dx = 1/512")
    name_to_c = {v: k for k, v in B3X_NAME.items()}
    n_b3y = 0
    for r in read_csv(B3Y_COMPARISON):
        c = name_to_c.get(r.get("config"))
        if c is None:
            continue
        E = _f(r["E"], "B3Y E")
        if (c, E) not in by_E:
            raise PilotCheckError(f"B3Y energy {r['E']} of {r['config']} is not bitwise a B3X ladder energy")
        p, lr = by_E[(c, E)]
        inv = int(r["inv_dx"])
        if not _same_bits(_f(r["T1_M1"], "B3Y T1"), _f(lr[f"T1_dx_1_{inv}"], "B3X T1")):
            raise PilotCheckError(f"B3Y and B3X T1 differ for {c} {p} at dx = 1/{inv}")
        n_b3y += 1
    if set(refs) != set(CONTROL_KEYS):
        raise PilotCheckError("control reference keys incomplete")
    return refs, n_b3y


def g0_references(read_csv, expected_energies) -> tuple:
    """G0-REF: reference consistency plus bitwise equality of every reference E with the B3 construction."""
    refs, n_b3y = control_references(read_csv)
    problems = []
    for c in CONFIGS:
        exp = expected_energies(c)
        for p in LADDER_POINTS:
            E = refs[(c, p, 512)][0]
            if p not in exp or not _same_bits(exp[p][1], E):
                problems.append(f"{c} {p}: reference E {E!r} != construction {exp.get(p)!r}")
    res = CheckResult("G0-REF", "G0", FAIL if problems else PASS,
                      evidence=dict(n_keys=len(refs), n_b3y_overlap=n_b3y, primary=B3X_LADDER),
                      message="; ".join(problems))
    return res, refs


def within_tolerance(value, ref, tol=CONTROL_TOL) -> bool:
    return abs(value - ref) <= tol


def g0_controls(refs, records, code_files) -> CheckResult:
    """G0-6/7: every control key has exactly one M1 record (T7 RGFRecord, base transverse) with E bitwise equal to
    the reference, the M1 method fingerprint, and |T1 - T1_committed| <= 1e-12."""
    problems, rows = [], []
    got = {}
    for r in records:
        if not isinstance(r, rr.RGFRecord):
            raise PilotCheckError(f"control record of type {type(r).__name__}")
        k = (r.config, r.point, r.inv_dx)
        if k in got:
            raise PilotCheckError(f"duplicate control record {k}")
        got[k] = r
    missing = sorted(set(CONTROL_KEYS) - set(got))
    extra = sorted(set(got) - set(CONTROL_KEYS))
    if missing or extra:
        raise PilotCheckError(f"control records incomplete: missing {missing}; unexpected {extra}")
    fp = rr.method_fingerprint("M1", code_files)
    for k in CONTROL_KEYS:
        r, (E, T1) = got[k], refs[k]
        if r.method != "M1" or (r.n_y, r.n_q) != rr.BASE_TRANSVERSE:
            problems.append(f"{k}: not a base-transverse M1 record")
            continue
        if r.method_fingerprint != fp:
            problems.append(f"{k}: M1 fingerprint mismatch")
        if not _same_bits(r.E, E):
            problems.append(f"{k}: E {r.E!r} != reference {E!r}")
        ok = within_tolerance(r.T1, T1)
        if not ok:
            problems.append(f"{k}: |T1 - committed| = {abs(r.T1 - T1):.3e} > {CONTROL_TOL}")
        rows.append(dict(config=k[0], point=k[1], inv_dx=k[2], T1=r.T1, T1_committed=T1,
                         abs_diff=abs(r.T1 - T1), within=ok))
    return CheckResult("G0-6/7", "G0", FAIL if problems else PASS,
                       evidence=dict(n=len(rows), tolerance=CONTROL_TOL, rows=rows), message="; ".join(problems))


# ======================================================================================================================
# G0 orchestration
# ======================================================================================================================
@dataclass(frozen=True)
class G0Deps:
    protocol_text: object          # repo -> protocol text (working tree; hash verified by G0-1 first)
    grid: object                   # validate_rgf_b2-like: kx_window, e_perp, adjust_energy, packet_weights, DKX
    read_csv: object               # relative path -> list of dict rows (committed artifact)
    expected_energies: object      # config -> {point: (kx, E)} from the B3 construction
    code_files: object             # repo -> provenance.code_fingerprint.code_files mapping
    m1_control: object             # (config, point, inv_dx, E, kx, code_files) -> RGFRecord  (PHYSICAL: M1 solve)


G0_PLAN = ("G0-1", "G0-2", "G0-4", "G0-5", "G0-REF", "G0-3a", "G0-3b", "G0-CTRL", "G0-3c", "G0-6/7", "G0-8")


@dataclass(frozen=True)
class G0Run:
    report: GateReport
    control_records: tuple = ()


def run_g0(repo, deps: G0Deps, *, baseline=pio.A6_BASELINE_COMMIT, spec: ProtectedSpec = ProtectedSpec()) -> G0Run:
    """Gate G0, fail closed, in the order of G0_PLAN. S_pre is taken immediately before the first control solve
    and S_post immediately after the last; nothing else runs in between."""
    seq = _Sequence("G0", G0_PLAN)
    head = sha = ""
    r = seq.guard("G0-1", "G0", g0_protocol, repo)
    if not seq.add(r):
        return G0Run(seq.report(head, sha))
    sha = r.evidence["protocol_sha256"]
    r = seq.guard("G0-2", "G0", g0_clean_tree, repo)
    if not seq.add(r):
        return G0Run(seq.report(head, sha))
    head = r.evidence["head"]
    text = deps.protocol_text(repo)
    if not seq.add(seq.guard("G0-4", "G0", g0_point_table, text, deps.grid)):
        return G0Run(seq.report(head, sha))
    if not seq.add(seq.guard("G0-5", "G0", g0_moved_nodes, text, deps.grid)):
        return G0Run(seq.report(head, sha))
    holder = {}

    def _refs():
        res, holder["refs"] = g0_references(deps.read_csv, deps.expected_energies)
        return res
    if not seq.add(seq.guard("G0-REF", "G0", _refs)):
        return G0Run(seq.report(head, sha))

    def _base():
        res, holder["paths"], holder["base"] = g0_baseline(repo, baseline, spec)
        return res
    if not seq.add(seq.guard("G0-3a", "G0", _base)):
        return G0Run(seq.report(head, sha))
    code_files = deps.code_files(repo)
    energies = {c: deps.expected_energies(c) for c in CONFIGS}

    def _pre():
        res, holder["pre"] = snapshot_check(repo, "G0-3b", holder["paths"], holder["base"], spec, baseline)
        return res
    if not seq.add(seq.guard("G0-3b", "G0", _pre)):
        return G0Run(seq.report(head, sha))
    records = []

    def _controls():                                   # S_pre ... controls ... S_post: nothing else in between
        for (c, p, inv) in CONTROL_KEYS:
            kx, E = energies[c][p]
            records.append(deps.m1_control(c, p, inv, E, kx, code_files))
        return CheckResult("G0-CTRL", "G0", PASS, evidence=dict(n_solves=len(records)))
    ctrl = seq.guard("G0-CTRL", "G0", _controls)

    def _post():
        res, _ = snapshot_check(repo, "G0-3c", holder["paths"], holder["base"], spec, baseline, pre=holder["pre"])
        return res
    post = seq.guard("G0-3c", "G0", _post)                # always taken once the control phase was entered
    if not seq.add(ctrl):
        seq.results.append(post)
        return G0Run(seq.report(head, sha))
    if not seq.add(post):
        return G0Run(seq.report(head, sha))
    if not seq.add(seq.guard("G0-6/7", "G0", g0_controls, holder["refs"], records, code_files)):
        return G0Run(seq.report(head, sha))
    seq.add(seq.guard("G0-8", "G0", g0_clean_tree, repo, "G0-8", head))
    return G0Run(seq.report(head, sha), tuple(records))


# ======================================================================================================================
# G-REAL primitives
# ======================================================================================================================
def disorder_params_hash(d) -> str:
    """sha256 over contiguous little-endian float64 concat(qx, qy, phi) (equals the B3/B3X record recipe)."""
    a = np.ascontiguousarray(np.concatenate([d.qx, d.qy, d.phi]), dtype="<f8")
    return "sha256:" + hashlib.sha256(a.tobytes()).hexdigest()


def disorder_fields(d) -> dict:
    return dict(W=float(d.strength), zeta=float(d.corr_length), region=None if d.region is None else
                [float(x) for x in d.region], edge=float(d.edge), M=int(d.n_modes), seed=int(d.seed))


def record_fields(rec) -> dict:
    return dict(W=float(rec["W"]), zeta=float(rec["zeta"]), region=[float(x) for x in rec["disorder_region"]],
                edge=float(rec["disorder_edge"]), M=int(rec["disorder_n_modes"]), seed=int(rec["seed"]))


def l1_identity(config, sources, records) -> CheckResult:
    """L1 for one configuration. sources: {name: disorder object}; records: {name: committed record dict}."""
    if not sources:
        raise PilotCheckError("no disorder sources")
    ids = {n: (disorder_params_hash(d), disorder_fields(d)) for n, d in sources.items()}
    ids.update({n: (r["disorder_params_hash"], record_fields(r)) for n, r in records.items()})
    ref_name = next(iter(sources))
    ref = ids[ref_name]
    bad = sorted(n for n, v in ids.items() if v != ref)
    return CheckResult(f"L1/{config}", "L1", FAIL if bad else PASS,
                       evidence=dict(reference=ref_name, hash=ref[0], fields=ref[1], sources=sorted(ids),
                                     mismatching=bad, b3y="no disorder hash recorded (availability only)"),
                       message=f"identity differs for {bad}" if bad else "")


def _first_difference(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or a.dtype != b.dtype:
        return f"shape/dtype differ: {a.shape} {a.dtype} vs {b.shape} {b.dtype}"
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return "non-finite values present"
    bad = np.argwhere(a != b)
    if bad.size == 0:
        return None
    idx = tuple(int(i) for i in bad[0])
    return f"{bad.shape[0]} differing entries; first at {idx}: {a[idx]!r} vs {b[idx]!r}"


def array_sha(a) -> str:
    return "sha256:" + hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


@dataclass(frozen=True)
class GridSample:
    """What one grid stores and what a direct Potential.numpy call gives at independently built coordinates."""
    name: str
    V_stored: np.ndarray
    V_direct: np.ndarray
    x_stored: np.ndarray
    x_direct: np.ndarray
    y_stored: np.ndarray
    y_direct: np.ndarray
    psi0: np.ndarray | None = None
    info: Mapping = field(default_factory=dict)


def l2_fidelity(check_id, s: GridSample) -> CheckResult:
    diffs = {k: _first_difference(a, b) for k, (a, b) in dict(
        V=(s.V_stored, s.V_direct), x=(s.x_stored, s.x_direct), y=(s.y_stored, s.y_direct)).items()}
    bad = {k: v for k, v in diffs.items() if v is not None}
    ev = dict(grid=s.name, potential_sha256=array_sha(s.V_stored), shape=list(np.shape(s.V_stored)),
              bitwise_equal=not bad, **{k: v for k, v in dict(s.info).items()})
    return CheckResult(check_id, "L2", FAIL if bad else PASS, evidence=ev,
                       message="; ".join(f"{k}: {v}" for k, v in bad.items()))


@dataclass(frozen=True)
class VariantGeometry:
    spec: VariantSpec
    cfg: object
    x_phys: np.ndarray
    dx_prod: float
    dx_variant: float
    evidence: Mapping


def variant_geometry(cfg_td, spec: VariantSpec) -> VariantGeometry:
    """§3.3.1 construction: x_phys,j = (j + 1 - m) dx from integer indices (never x_j - s); L_x = 3.9 + 2 m dx;
    nx = 249 + 2 m; CAP width 0.6 + s for the widened variants."""
    for k in ("Lx", "nx", "ny", "Ly", "cap_width", "cap_strength", "cap_power"):
        if getattr(cfg_td, k) != PRODUCTION[k]:
            raise PilotCheckError(f"TD production config has {k} = {getattr(cfg_td, k)!r}, expected {PRODUCTION[k]!r}")
    m = spec.m
    dx = cfg_td.Lx / (cfg_td.nx + 1)
    s = m * dx
    Lx = cfg_td.Lx + 2 * m * dx
    nx = cfg_td.nx + 2 * m
    width = cfg_td.cap_width + s if spec.cap_width_extends else cfg_td.cap_width
    cfg_s = cfg_td.replace(Lx=Lx, nx=nx, cap_width=width, cap_strength=spec.cap_strength, cap_power=spec.cap_power)
    x_phys = (np.arange(nx, dtype=np.int64) + 1 - m) * dx
    dxv = cfg_s.Lx / (cfg_s.nx + 1)
    ev = dict(variant=spec.name, m=m, s=s, L_x=Lx, nx=nx, cap_width=width, cap_strength=spec.cap_strength,
              dx_production=dx, dx_variant=dxv, dx_equal=bool(_same_bits(dxv, dx)))
    return VariantGeometry(spec, cfg_s, x_phys, dx, dxv, types.MappingProxyType(ev))


def dx_evidence(geo: VariantGeometry, solver_dx=None) -> CheckResult:
    ev = dict(geo.evidence)
    if solver_dx is not None:
        ev["dx_solver"] = float(solver_dx)
        ev["dx_solver_equal"] = bool(_same_bits(solver_dx, geo.dx_prod))
    ok = ev["dx_equal"] and ev.get("dx_solver_equal", True)
    return CheckResult(f"DX/{geo.spec.name}", "DX", PASS if ok else FAIL, blocking=False, evidence=ev,
                       message="" if ok else "L_x/(nx + 1) differs from the production dx (recorded; §3.3.1)")


def l3_variant(check_id, prod: GridSample, var: GridSample, m: int) -> CheckResult:
    """x_variant[i + m] == x_production[i], V_variant[i + m, :] == V_production[i, :] and psi0 likewise, bitwise,
    for every production index i; y grids equal. Points outside the production domain: max |V| recorded (L2 only)."""
    n = np.shape(prod.V_stored)[0]
    nv = np.shape(var.V_stored)[0]
    problems = []
    if m < 0 or m + n > nv:
        raise PilotCheckError(f"offset m = {m} with {n} production points does not fit a variant of {nv} points")
    sl = slice(m, m + n)
    for name, a, b in (("x", var.x_stored[sl], prod.x_stored), ("V", var.V_stored[sl, :], prod.V_stored),
                       ("psi0", None if var.psi0 is None else var.psi0[sl, :], prod.psi0),
                       ("y", var.y_stored, prod.y_stored)):
        if a is None or b is None:
            problems.append(f"{name}: missing")
            continue
        d = _first_difference(a, b)
        if d is not None:
            problems.append(f"{name}: {d}")
    outside = np.concatenate([np.ravel(var.V_stored[:m, :]), np.ravel(var.V_stored[m + n:, :])])
    ev = dict(variant=var.name, m=m, variant_index_range=[m, m + n - 1], production_index_range=[0, n - 1],
              max_abs_V_outside=float(np.max(np.abs(outside))) if outside.size else 0.0)
    return CheckResult(check_id, "L3", FAIL if problems else PASS, evidence=ev, message="; ".join(problems))


def l3_full_equal(check_id, a: GridSample, b: GridSample) -> CheckResult:
    problems = []
    for name in ("x_stored", "V_stored", "psi0", "y_stored"):
        x, y = getattr(a, name), getattr(b, name)
        d = "missing" if x is None or y is None else _first_difference(x, y)
        if d is not None:
            problems.append(f"{name}: {d}")
    return CheckResult(check_id, "L3", FAIL if problems else PASS, evidence=dict(a=a.name, b=b.name),
                       message="; ".join(problems))


def l3_rgf_nested(config, coarse_inv, fine_inv, coarse: GridSample, fine: GridSample) -> CheckResult:
    """Coincident RGF slices across nested dx (base transverse grid): coordinates and samples bitwise."""
    if fine_inv % coarse_inv:
        raise PilotCheckError("resolutions are not nested")
    step = fine_inv // coarse_inv
    n = np.shape(coarse.V_stored)[0]
    problems = []
    xf = np.asarray(fine.x_stored)[::step][:n]
    Vf = np.asarray(fine.V_stored)[::step][:n]
    for name, a, b in (("x", xf, coarse.x_stored), ("V", Vf, coarse.V_stored), ("y", fine.y_stored,
                                                                                 coarse.y_stored)):
        d = _first_difference(a, b)
        if d is not None:
            problems.append(f"{name}: {d}")
    return CheckResult(f"L3/{config}/RGF/{coarse_inv}-{fine_inv}", "L3", FAIL if problems else PASS,
                       evidence=dict(step=step, n_coincident=int(n)), message="; ".join(problems))


def not_applicable(check_id, reason) -> CheckResult:
    return CheckResult(check_id, "L3", NOT_APPLICABLE, blocking=False, evidence=dict(reason=reason))


# ======================================================================================================================
# G-REAL orchestration
# ======================================================================================================================
@dataclass(frozen=True)
class GRealDeps:
    protocol_text: object          # repo -> protocol text
    case_disorder: object          # config -> disorder object used by the RGF devices (case_potential path)
    case_config: object            # config -> B3 configuration
    td_config: object              # config -> committed TD production configuration (seed 0 present)
    disorder_from_config: object   # (cfg, seed) -> disorder object
    record_identities: object      # config -> {name: committed record dict with disorder_params_hash and fields}
    rgf_device: object             # (config, inv_dx, n_y, n_q) -> GridSample  (may raise NotEvaluable)
    td_production: object          # (config, cfg_td) -> (GridSample, solver_dx)        constructor only
    td_variant: object             # (config, cfg_td, VariantGeometry) -> (GridSample, solver_dx)  constructor only


def _require_passed(report, gate):
    if not isinstance(report, GateReport) or report.gate != gate:
        raise PilotCheckError(f"a {gate} GateReport is required")
    if report.status != PASS:
        raise PilotCheckError(f"{gate} has not passed (status {report.status}); nothing later may run")


def run_greal(repo, deps: GRealDeps, *, g0: GateReport) -> GateReport:
    """Gate G-REAL after a PASS G0 (refused otherwise): L1 -> L2 -> L3, fail closed; DX evidence always recorded."""
    _require_passed(g0, "G0")
    head = pio.head_commit(repo)
    if head != g0.head:
        raise PilotCheckError(f"HEAD changed since G0: {g0.head} -> {head}")
    parse_variant_table(deps.protocol_text(repo))        # §3.3.1 table must agree with VARIANTS (else refuse)
    planned = ([f"L1/{c}" for c in CONFIGS]
               + [f"L2/{c}/RGF/{i}/{y}x{q}" for c in CONFIGS for (i, y, q) in RGF_DEVICES]
               + [f"L2/{c}/TD/production" for c in CONFIGS]
               + [f"L2/{c}/TD/{v.name}" for c in CONFIGS for v in VARIANTS]
               + [f"L3/{c}/RGF/{a}-{b}" for c in CONFIGS for (a, b) in RGF_NESTED_PAIRS]
               + [f"L3/{c}/TD/{v.name}" for c in CONFIGS for v in VARIANTS]
               + [f"L3/{c}/TD/V1=V2" for c in CONFIGS])
    seq = _Sequence("G-REAL", planned)

    def done():
        return seq.report(head, g0.protocol_sha256)

    cfg_td, geos = {}, {}
    for c in CONFIGS:                                                    # ---- L1
        def _l1(c=c):
            cfg_td[c] = deps.td_config(c)
            geos[c] = {v.name: variant_geometry(cfg_td[c], v) for v in VARIANTS}
            src = {"rgf_case_potential": deps.case_disorder(c),
                   "rgf_fresh": deps.disorder_from_config(deps.case_config(c), 0),
                   "td_production": deps.disorder_from_config(cfg_td[c], 0)}
            src.update({f"td_{n}": deps.disorder_from_config(g.cfg, 0) for n, g in geos[c].items()})
            return l1_identity(c, src, deps.record_identities(c))
        if not seq.add(seq.guard(f"L1/{c}", "L1", _l1)):
            return done()
    rgf, prod, var, dxres = {}, {}, {}, []
    for c in CONFIGS:                                                    # ---- L2
        for (inv, ny, nq) in RGF_DEVICES:
            cid = f"L2/{c}/RGF/{inv}/{ny}x{nq}"

            def _dev(c=c, inv=inv, ny=ny, nq=nq, cid=cid):
                s = deps.rgf_device(c, inv, ny, nq)
                rgf[(c, inv, ny, nq)] = s
                return l2_fidelity(cid, s)
            if not seq.add(seq.guard(cid, "L2", _dev)):
                return done()

        def _prod(c=c):
            s, _ = deps.td_production(c, cfg_td[c])
            prod[c] = s
            return l2_fidelity(f"L2/{c}/TD/production", s)
        if not seq.add(seq.guard(f"L2/{c}/TD/production", "L2", _prod)):
            return done()
        for v in VARIANTS:
            def _var(c=c, v=v):
                s, solver_dx = deps.td_variant(c, cfg_td[c], geos[c][v.name])
                var[(c, v.name)] = s
                if v.m:
                    dxres.append(dx_evidence(geos[c][v.name], solver_dx))
                return l2_fidelity(f"L2/{c}/TD/{v.name}", s)
            if not seq.add(seq.guard(f"L2/{c}/TD/{v.name}", "L2", _var)):
                return done()
    for c in CONFIGS:                                                    # ---- L3
        for (a, b) in RGF_NESTED_PAIRS:
            cid = f"L3/{c}/RGF/{a}-{b}"
            if not seq.add(seq.guard(cid, "L3", l3_rgf_nested, c, a, b, rgf[(c, a, 63, 127)], rgf[(c, b, 63, 127)])):
                return done()
        seq.results.append(not_applicable(f"L3/{c}/RGF/Nq254", "N_q = 254 devices share no transverse sampling "
                                          "points with N_q = 127; no L3 claim (L2 only)"))
        seq.results.append(not_applicable(f"L3/{c}/RGF-TD", "RGF slices and the TD grid share no physical points; "
                                          "RGF <-> TD is judged on L1 + L2 only"))
        for v in VARIANTS:
            cid = f"L3/{c}/TD/{v.name}"
            if not seq.add(seq.guard(cid, "L3", l3_variant, cid, prod[c], var[(c, v.name)], v.m)):
                return done()
        if not seq.add(seq.guard(f"L3/{c}/TD/V1=V2", "L3", l3_full_equal, f"L3/{c}/TD/V1=V2", var[(c, "V1")],
                                 var[(c, "V2")])):
            return done()
    for r in dxres:                                                      # ---- DX (non-blocking evidence)
        seq.results.append(r)
    return done()


# ======================================================================================================================
# clearance and T5 integration
# ======================================================================================================================
@dataclass(frozen=True)
class Clearance:
    head: str
    protocol_sha256: str
    classification: str = CLASSIFICATION


def require_clearance(g0, greal) -> Clearance:
    """The only entry for later pilot stages: both gates must be PASS (NOT EVALUABLE is never a pass)."""
    _require_passed(g0, "G0")
    _require_passed(greal, "G-REAL")
    if g0.head != greal.head or g0.protocol_sha256 != greal.protocol_sha256:
        raise PilotCheckError("G0 and G-REAL reports refer to different HEADs or protocol hashes")
    return Clearance(g0.head, g0.protocol_sha256)


def to_g0_outcome(g0, greal, repo) -> pio.G0Outcome:
    """T5 G0Outcome from two validated PASS reports with the pinned protocol hash and the live HEAD. T5 publish
    then applies all of its own checks (nothing is bypassed)."""
    cl = require_clearance(g0, greal)
    if cl.protocol_sha256 != pio.PROTOCOL_SHA256_REFERENCE:
        raise PilotCheckError("G0 report does not carry the pinned protocol reference")
    live = pio.head_commit(repo)
    if live != cl.head:
        raise PilotCheckError(f"HEAD changed since the gates: {cl.head} -> {live}")
    checks = {"G0": g0.status, "G-REAL": greal.status, "G-REAL_nonblocking_failures": list(greal.nonblocking_failures)}
    return pio.G0Outcome(passed=True, protocol_sha256=cl.protocol_sha256, head=cl.head, checks=checks)


def snapshots(g0) -> tuple:
    """(S_pre, S_post) protected-file snapshots of a PASS G0 report (for build_metadata protected_before/after)."""
    _require_passed(g0, "G0")
    return _thaw(g0.result("G0-3b").evidence["snapshot"]), _thaw(g0.result("G0-3c").evidence["snapshot"])


def metadata_extra(g0, greal) -> dict:
    require_clearance(g0, greal)
    l1 = {r.check_id.split("/")[1]: dict(disorder_params_hash=r.evidence["hash"], fields=_thaw(r.evidence["fields"]))
          for r in greal.results if r.level == "L1"}
    return {"gates": {"G0": g0.to_dict(), "G-REAL": greal.to_dict()},
            "disorder_identity": l1,
            "potential_samples": [dict(check=r.check_id, **_thaw(r.evidence)) for r in greal.results
                                  if r.level == "L2"],
            "common_coordinate_checks": [dict(check=r.check_id, status=r.status, **_thaw(r.evidence))
                                         for r in greal.results if r.level == "L3"],
            "shifted_grid_dx": [dict(check=r.check_id, status=r.status, **_thaw(r.evidence))
                                for r in greal.results if r.level == "DX"],
            "nonblocking_failures": list(greal.nonblocking_failures)}


def run_gates(repo, staging, g0_deps, greal_deps, **g0_kw):
    """G0 then G-REAL. Any failure discards the staging directory (outside the repository); nothing is published
    here. Returns (G0Run, G-REAL report or None)."""
    try:
        g0 = run_g0(repo, g0_deps, **g0_kw)
        if not g0.report.passed:
            staging.discard()
            return g0, None
        greal = run_greal(repo, greal_deps, g0=g0.report)
        if not greal.passed:
            staging.discard()
        return g0, greal
    except Exception:
        staging.discard()
        raise
