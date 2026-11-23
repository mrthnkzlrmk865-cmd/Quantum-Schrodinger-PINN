"""B4 limited diagnostic pilot — RGF record extraction and serialization helpers (phase T7).

Classification: DIAGNOSTIC / NOT B4.

Pure helpers for docs/B4_PILOT_PROTOCOL.md §1.3–§1.6 and §5:

* extraction of one protocol §1.4 record (full complex S, T1, R1, 1 - T1 - R1, U from both products, Q and its
  components, channel counts, method-specific extras) from the raw output of ONE method; every observable is computed
  from that record's own S with the existing definitions (validate_rgf_b3.T1 / unitarity / reciprocity /
  s_difference, smatrix_spec.split_S; R1 as in validate_rgf_b3y). Caller-supplied observables are never accepted;
* the complete protocol record set and strict completeness / same-problem validation;
* raw M1-M2, M1-M3, M2-M3 complex-S agreement rows, the M4 requirement set (protocol §1.3: an M1-M2 or M1-M3
  difference > 1e-10, plus the predeclared M4 points) and raw M4 adjudication inputs. No exceedance classification,
  precision consistency, decomposition, order classification, Richardson analysis or gate assembly is done here
  (that is the separate summarizer's job; protocol §1.4 "computed ... by a separate summarizer");
* lossless, deterministic serialization: arrays for b4_pilot_io.Staging.write_npz, rows for Staging.write_csv, and
  read-only, verifying loaders.

This module never writes a file and never calls a solver. M4 is never a fourth grid: it is excluded from the
three-grid marshalling used by the observed-order functions.
"""
from __future__ import annotations

import csv
import hashlib
import math
import re
from dataclasses import dataclass, field

import numpy as np

import smatrix_spec as sm
import sparse_exact_lattice as xp
import validate_rgf_b3 as b3
from b4_pilot_io import CLASSIFICATION, CSV_HEADER

SCHEMA_VERSION = "t7-v1"


class PilotRecordError(ValueError):
    """Malformed, incomplete or inconsistent pilot RGF data (fail closed)."""


# ---- protocol identifiers and values (§1.2, §1.3, §1.6, §0.1); quoted, never tuned -----------------------------
METHODS = ("M1", "M2", "M3", "M4")
GRID_METHODS = ("M1", "M2", "M3")
INV_DX = (512, 1024, 2048)
CONFIGS = ("A", "B", "C", "D")
POINTS_14 = ("L1", "L2", "L3", "T3−", "T3+", "T10−", "T10+", "T11−", "T11+", "U1", "U2", "U3",
             "C1", "C2")
B3_POINTS = ("B3-E1", "B3-E2", "B3-E3")                 # E = 88.51, 204.93, 281.90 (§1.6), ascending
POINT_ORDER = POINTS_14 + B3_POINTS
M4_PROTOCOL_POINTS = ("L1", "T3+", "T10+", "T11+", "U3")
M4_PROTOCOL_INV_DX = 2048
NYNQ_POINTS = ("U1", "U2", "U3", "T11−", "T11+")
NYNQ_INV_DX = 1024
BASE_TRANSVERSE = (63, 127)
TRANSVERSE_VARIANTS = ((63, 127), (63, 254), (127, 254))   # base; N_q 127 -> 254; N_y 63 -> 127 at N_q = 254
S_AGREEMENT = 1e-10
PAIRS = (("M1", "M2"), ("M1", "M3"), ("M2", "M3"))
M4_TRIGGER_PAIRS = (("M1", "M2"), ("M1", "M3"))            # §1.3: "an M1-M2/M3 difference exceeds 1e-10"

S_DTYPE = {"M1": np.dtype(np.complex128), "M2": np.dtype(np.complex128), "M3": np.dtype(np.complex128),
           "M4": np.dtype(xp.CLD)}

EXTRAS = {                                                    # exact key sets and types per method (§1.4)
    "M1": {"min_rcond": float, "argmin_rcond": int},
    "M2": {},
    "M3": {"residual_after_refinement": float, "residual_before_refinement": float},
    "M4": {"eta_final": float, "error_indicator": float, "X4_state_passed": bool, "X4_green_passed": bool,
           "X5_passed": bool, "X7_passed": bool, "gates_passed": bool},
}
ALL_EXTRAS = tuple(k for m in METHODS for k in EXTRAS[m])

METHOD_CODE_FILES = {                                         # repository imports of each method (verified)
    "M1": ("stationary.py", "smatrix_spec.py"),
    "M2": ("sparse_reference.py", "smatrix_spec.py"),
    "M3": ("sparse_wavefunction.py",),
    "M4": ("sparse_exact_lattice.py",),
}
EXTRACTION_CODE_FILES = ("validate_rgf_b3.py", "b4_pilot_rgf_records.py")

OBSERVABLES = ("T1", "R1", "one_minus_T1_R1", "U", "U_StS", "U_SSt", "Q", "Q_full", "Q_rL", "Q_rR", "Q_t", "Q_prob",
               "same_dir_RL_report_only", "same_dir_LR_report_only")

_SHA_RE = re.compile(r"^sha256:[0-9a-f]{64}$")
_INT_RE = re.compile(r"^-?[0-9]+$")
NOT_APPLICABLE = "NOT APPLICABLE"


# ----------------------------------------------------------------------------------------------------------------
# small validators
# ----------------------------------------------------------------------------------------------------------------
def _real(x, name):
    if isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, float, np.floating, np.integer)):
        raise PilotRecordError(f"{name} must be a real number, got {type(x).__name__}")
    v = float(x)
    if not math.isfinite(v):
        raise PilotRecordError(f"{name} must be finite, got {v!r}")
    return v


def _int(x, name):
    if isinstance(x, (bool, np.bool_)) or not isinstance(x, (int, np.integer)):
        raise PilotRecordError(f"{name} must be an integer, got {type(x).__name__}")
    return int(x)


def _sha(x, name):
    if not isinstance(x, str) or not _SHA_RE.match(x):
        raise PilotRecordError(f"{name} must be 'sha256:' + 64 lowercase hex digits, got {x!r}")
    return x


def _same_bits(a: float, b: float) -> bool:
    return np.float64(a).tobytes() == np.float64(b).tobytes()


def _check_configs(configs):
    cs = list(configs)
    if not cs:
        raise PilotRecordError("configs is empty")
    if len(set(cs)) != len(cs):
        raise PilotRecordError(f"configs contains duplicates: {cs}")
    bad = [c for c in cs if c not in CONFIGS]
    if bad:
        raise PilotRecordError(f"unknown configs {bad}; expected a subset of {CONFIGS}")
    return tuple(c for c in CONFIGS if c in cs)


# ----------------------------------------------------------------------------------------------------------------
# observables (existing definitions only)
# ----------------------------------------------------------------------------------------------------------------
def _observables(S, n_open: int) -> dict:
    """T1 = Σ_m |(t_RL)_m1|², R1 = Σ_m |(r_L)_m1|² (column 1 = incidence in left channel 1), U from both products,
    Q and its block-wise components — all computed from S with the existing B3/B3Y definitions."""
    u = b3.unitarity(S)
    q = b3.reciprocity(S, n_open)
    T1 = float(b3.T1(S, n_open))
    R1 = float(np.sum(np.abs(sm.split_S(S, n_open)["r_L"][:, 0]) ** 2))       # validate_rgf_b3y.run_methods
    out = dict(T1=T1, R1=R1, one_minus_T1_R1=1.0 - T1 - R1, U=float(u["U"]), U_StS=float(u["U_StS"]),
               U_SSt=float(u["U_SSt"]))
    out.update({k: float(q[k]) for k in ("Q", "Q_full", "Q_rL", "Q_rR", "Q_t", "Q_prob", "same_dir_RL_report_only",
                                         "same_dir_LR_report_only")})
    return out


def split_hi_lo(S) -> tuple:
    """Exact two-double representation of a clongdouble S: hi = S rounded to complex128, lo = (S - hi) in complex128.
    Refuses if hi + lo does not reproduce S exactly (e.g. lo underflow)."""
    S = np.asarray(S)
    if S.dtype != S_DTYPE["M4"]:
        raise PilotRecordError(f"hi/lo split needs dtype {S_DTYPE['M4']}, got {S.dtype}")
    hi = S.astype(np.complex128)
    lo = (S - hi.astype(xp.CLD)).astype(np.complex128)
    if not (np.all(np.isfinite(hi)) and np.all(np.isfinite(lo))):
        raise PilotRecordError("hi/lo split produced non-finite values")
    if not np.array_equal(hi.astype(xp.CLD) + lo.astype(xp.CLD), S):
        raise PilotRecordError("clongdouble S is not exactly representable as hi + lo complex128 (refused: lossy)")
    return hi, lo


def s_digest(S) -> str:
    """Content hash of the stored S: dtype tag, shape and C-contiguous bytes (complex128), or of hi and lo (M4)."""
    S = np.asarray(S)
    h = hashlib.sha256()
    if S.dtype == np.complex128:
        h.update(b"complex128\0" + repr(S.shape).encode() + b"\0" + np.ascontiguousarray(S).tobytes())
    elif S.dtype == S_DTYPE["M4"]:
        hi, lo = split_hi_lo(S)
        h.update(b"clongdouble_hi_lo\0" + repr(S.shape).encode() + b"\0" + np.ascontiguousarray(hi).tobytes()
                 + b"\0" + np.ascontiguousarray(lo).tobytes())
    else:
        raise PilotRecordError(f"unsupported S dtype {S.dtype}")
    return "sha256:" + h.hexdigest()


def _digest_hi_lo(hi, lo) -> str:
    h = hashlib.sha256()
    h.update(b"clongdouble_hi_lo\0" + repr(hi.shape).encode() + b"\0" + np.ascontiguousarray(hi).tobytes()
             + b"\0" + np.ascontiguousarray(lo).tobytes())
    return "sha256:" + h.hexdigest()


# ----------------------------------------------------------------------------------------------------------------
# method fingerprint (protocol §1.4 "method version (code fingerprint)")
# ----------------------------------------------------------------------------------------------------------------
def method_files(method: str) -> tuple:
    if method not in METHODS:
        raise PilotRecordError(f"unknown method {method!r}")
    return tuple(sorted(set(METHOD_CODE_FILES[method]) | set(EXTRACTION_CODE_FILES)))


def method_fingerprint(method: str, code_files) -> str:
    """sha256 over (file name, per-file sha256) of the method's files, taken from the T5 metadata map
    provenance.code_fingerprint.code_files. A missing or malformed entry is refused."""
    if not isinstance(code_files, dict):
        raise PilotRecordError("code_files must be the provenance.code_fingerprint.code_files mapping")
    h = hashlib.sha256()
    for name in method_files(method):
        if name not in code_files:
            raise PilotRecordError(f"code fingerprint for {name!r} (method {method}) missing from metadata")
        h.update(name.encode() + b"\0" + _sha(code_files[name], f"code_files[{name!r}]").encode() + b"\0")
    return "sha256:" + h.hexdigest()


def _metadata_code_files(metadata) -> dict:
    try:
        cf = metadata["provenance"]["code_fingerprint"]["code_files"]
    except (KeyError, TypeError) as exc:
        raise PilotRecordError("metadata lacks provenance.code_fingerprint.code_files") from exc
    if not isinstance(cf, dict):
        raise PilotRecordError("metadata provenance.code_fingerprint.code_files is not a mapping")
    return cf


def verify_method_fingerprints(entries, metadata) -> None:
    """Every entry (an RGFRecord or a loaded/parsed mapping with 'method' and 'method_fingerprint') must carry exactly
    the fingerprint recomputed from the run metadata."""
    cf = _metadata_code_files(metadata)
    expected = {m: method_fingerprint(m, cf) for m in METHODS}
    for e in entries:
        m = e.method if isinstance(e, RGFRecord) else e["method"]
        fp = e.method_fingerprint if isinstance(e, RGFRecord) else e["method_fingerprint"]
        if m not in expected:
            raise PilotRecordError(f"unknown method {m!r}")
        if fp != expected[m]:
            raise PilotRecordError(f"method fingerprint mismatch for {m} at {_entry_key(e)}: record {fp}, "
                                   f"metadata {expected[m]}")


# ----------------------------------------------------------------------------------------------------------------
# the record
# ----------------------------------------------------------------------------------------------------------------
def _check_placement(config, point, n_y, n_q, inv_dx, method):
    if config not in CONFIGS:
        raise PilotRecordError(f"unknown config {config!r}")
    if point not in POINT_ORDER:
        raise PilotRecordError(f"unknown point {point!r}")
    if method not in METHODS:
        raise PilotRecordError(f"unknown method {method!r}")
    if inv_dx not in INV_DX:
        raise PilotRecordError(f"unknown inv_dx {inv_dx!r}; expected one of {INV_DX}")
    tv = (n_y, n_q)
    if tv not in TRANSVERSE_VARIANTS:
        raise PilotRecordError(f"unknown transverse variant {tv}; expected one of {TRANSVERSE_VARIANTS}")
    if tv != BASE_TRANSVERSE:
        if method != "M1" or point not in NYNQ_POINTS or inv_dx != NYNQ_INV_DX:
            raise PilotRecordError(f"transverse variant {tv} is allowed only for M1 at {NYNQ_POINTS}, "
                                   f"inv_dx {NYNQ_INV_DX}; got {method} at {point}, inv_dx {inv_dx}")
    if point in B3_POINTS and method != "M1":
        raise PilotRecordError(f"B3 energy points are re-run with M1 only (§1.6); got {method} at {point}")


@dataclass(frozen=True, eq=False)
class RGFRecord:
    """One protocol §1.4 row: one configuration x point x transverse variant x dx x method. Validated on
    construction; the observables must equal those recomputed from S (bitwise)."""
    config: str
    point: str
    n_y: int
    n_q: int
    inv_dx: int
    method: str
    E: float
    kx: float
    S: np.ndarray
    n_open: int
    n_closed: int
    T1: float
    R1: float
    one_minus_T1_R1: float
    U: float
    U_StS: float
    U_SSt: float
    Q: float
    Q_full: float
    Q_rL: float
    Q_rR: float
    Q_t: float
    Q_prob: float
    same_dir_RL_report_only: float
    same_dir_LR_report_only: float
    extras: dict
    V_lead_max: float
    potential_sha256: str
    n_slices: int
    x_start: float
    wall_time_s: float
    method_fingerprint: str
    classification: str = CLASSIFICATION
    s_sha256: str = field(default="", init=False)

    def __post_init__(self):
        for k in ("config", "point", "method"):
            if not isinstance(getattr(self, k), str):
                raise PilotRecordError(f"{k} must be a string")
        n_y, n_q, inv = _int(self.n_y, "n_y"), _int(self.n_q, "n_q"), _int(self.inv_dx, "inv_dx")
        _check_placement(self.config, self.point, n_y, n_q, inv, self.method)
        if self.classification != CLASSIFICATION:
            raise PilotRecordError("classification must be DIAGNOSTIC / NOT B4")
        for k in ("E", "kx", "V_lead_max", "x_start", "wall_time_s"):
            _real(getattr(self, k), k)
        if self.wall_time_s < 0:
            raise PilotRecordError("wall_time_s must be >= 0")
        if _int(self.n_slices, "n_slices") <= 0:
            raise PilotRecordError("n_slices must be > 0")
        _sha(self.potential_sha256, "potential_sha256")
        _sha(self.method_fingerprint, "method_fingerprint")
        n_open, n_closed = _int(self.n_open, "n_open"), _int(self.n_closed, "n_closed")
        if not 1 <= n_open <= n_y - 1:
            raise PilotRecordError(f"n_open = {n_open} must satisfy 1 <= n_open <= n_y - 1 = {n_y - 1}")
        if n_closed != n_y - n_open:
            raise PilotRecordError(f"n_closed = {n_closed} != n_y - n_open = {n_y - n_open}")
        S = self.S
        if not isinstance(S, np.ndarray):
            raise PilotRecordError("S must be a numpy array")
        if S.dtype != S_DTYPE[self.method]:
            raise PilotRecordError(f"S dtype for {self.method} must be {S_DTYPE[self.method]}, got {S.dtype}")
        if S.ndim != 2 or S.shape != (2 * n_open, 2 * n_open):
            raise PilotRecordError(f"S shape {S.shape} != (2 n_open, 2 n_open) = {(2 * n_open, 2 * n_open)}")
        if not np.all(np.isfinite(S)):
            raise PilotRecordError("S contains non-finite entries")
        Sc = np.array(S, copy=True, order="C")
        Sc.setflags(write=False)
        object.__setattr__(self, "S", Sc)
        ext = self.extras
        if not isinstance(ext, dict) or set(ext) != set(EXTRAS[self.method]):
            got = sorted(ext) if isinstance(ext, dict) else ext
            raise PilotRecordError(f"extras for {self.method} must be exactly {sorted(EXTRAS[self.method])}, got {got}")
        clean = {}
        for k, typ in EXTRAS[self.method].items():
            v = ext[k]
            if typ is bool:
                if not isinstance(v, (bool, np.bool_)):
                    raise PilotRecordError(f"extras[{k!r}] must be a bool")
                clean[k] = bool(v)
            elif typ is int:
                clean[k] = _int(v, f"extras[{k!r}]")
            else:
                clean[k] = _real(v, f"extras[{k!r}]")
        if self.method == "M4" and clean["gates_passed"] != all(
                clean[k] for k in ("X4_state_passed", "X4_green_passed", "X5_passed", "X7_passed")):
            raise PilotRecordError("extras['gates_passed'] must equal X4_state and X4_green and X5 and X7")
        object.__setattr__(self, "extras", clean)
        obs = _observables(Sc, n_open)
        for k in OBSERVABLES:
            v = _real(getattr(self, k), k)
            if not _same_bits(v, obs[k]):
                raise PilotRecordError(f"{k} = {v!r} differs from the value recomputed from S ({obs[k]!r}); "
                                       "observables are never accepted from the caller")
        object.__setattr__(self, "s_sha256", s_digest(Sc))

    @property
    def key(self):
        return record_key(self)


def record_key(rec) -> tuple:
    return (rec.config, rec.point, rec.n_y, rec.n_q, rec.inv_dx, rec.method)


def _entry_key(e):
    if isinstance(e, RGFRecord):
        return record_key(e)
    return (e["config"], e["point"], e["n_y"], e["n_q"], e["inv_dx"], e["method"])


def _order(key):
    c, p, ny, nq, inv, m = key
    return (CONFIGS.index(c), POINT_ORDER.index(p), TRANSVERSE_VARIANTS.index((ny, nq)), INV_DX.index(inv),
            METHODS.index(m))


def sort_records(records) -> list:
    recs = list(records)
    for r in recs:
        if not isinstance(r, RGFRecord):
            raise PilotRecordError(f"expected RGFRecord, got {type(r).__name__}")
    keys = [record_key(r) for r in recs]
    if len(set(keys)) != len(keys):
        dup = sorted({k for k in keys if keys.count(k) > 1}, key=_order)
        raise PilotRecordError(f"duplicate record keys: {dup}")
    return sorted(recs, key=lambda r: _order(record_key(r)))


def _build(method, *, config, point, inv_dx, n_y, n_q, E, kx, S, n_open, extras, V_lead_max, potential_sha256,
           n_slices, x_start, wall_time_s, code_files):
    S = np.asarray(S)
    n_open = _int(n_open, "n_open")
    if S.ndim != 2 or S.shape != (2 * n_open, 2 * n_open):
        raise PilotRecordError(f"solver S shape {S.shape} inconsistent with solver n_open = {n_open}")
    obs = _observables(S, n_open)
    return RGFRecord(config=config, point=point, n_y=n_y, n_q=n_q, inv_dx=inv_dx, method=method, E=E, kx=kx, S=S,
                     n_open=n_open, n_closed=_int(n_y, "n_y") - n_open, extras=extras, V_lead_max=V_lead_max,
                     potential_sha256=potential_sha256, n_slices=n_slices, x_start=x_start, wall_time_s=wall_time_s,
                     method_fingerprint=method_fingerprint(method, code_files), **obs)


def extract_record(*, config, point, inv_dx, method, n_y, n_q, E, kx, solver_out, V_lead_max, potential_sha256,
                   n_slices, x_start, wall_time_s, code_files) -> RGFRecord:
    """Record for M1, M2 or M3 from the raw solver dict (stationary.scattering_matrix,
    sparse_reference.sparse_scattering_matrix or sparse_wavefunction.wavefunction_scattering_matrix)."""
    if method not in GRID_METHODS:
        raise PilotRecordError(f"extract_record handles {GRID_METHODS}; use extract_m4_record for M4")
    try:
        S, n_open = solver_out["S"], solver_out["n_open"]
        if method == "M1":
            extras = dict(min_rcond=solver_out["min_rcond"], argmin_rcond=solver_out["argmin_rcond"])
        elif method == "M2":
            extras = {}
        else:
            extras = dict(residual_after_refinement=solver_out["residual_after_refinement"],
                          residual_before_refinement=solver_out["residual_before_refinement"])
    except (KeyError, TypeError) as exc:
        raise PilotRecordError(f"solver output for {method} lacks a required field: {exc}") from exc
    return _build(method, config=config, point=point, inv_dx=inv_dx, n_y=n_y, n_q=n_q, E=E, kx=kx, S=S,
                  n_open=n_open, extras=extras, V_lead_max=V_lead_max, potential_sha256=potential_sha256,
                  n_slices=n_slices, x_start=x_start, wall_time_s=wall_time_s, code_files=code_files)


def _platform_evaluable(platform_report) -> bool:
    if not isinstance(platform_report, dict) or not isinstance(platform_report.get("evaluable"), bool):
        raise PilotRecordError("platform_report must be sparse_exact_lattice.platform_precision() output")
    return platform_report["evaluable"]


def extract_m4_record(*, config, point, inv_dx, n_y, n_q, E, kx, solver_out, V_lead_max, potential_sha256,
                      n_slices, x_start, wall_time_s, code_files, platform_report) -> RGFRecord:
    """Record for M4 (sparse_exact_lattice.exact_lattice_S). Refused unless the platform gate is evaluable."""
    if not _platform_evaluable(platform_report):
        raise PilotRecordError("M4 is NOT EVALUABLE on this platform (extended-precision gate failed); no record")
    try:
        S, n_open = solver_out["S"], solver_out["n_open"]
        e_ind = _real(solver_out["error_indicator"], "error_indicator")
        x4s = bool(solver_out["refinement_gate"]["passed"])
        x4g = bool(solver_out["green_refinement_gate"]["passed"])
        x5 = bool(_real(solver_out["fl_vs_state_max_abs"], "fl_vs_state_max_abs") <= xp.FL_WF_TOL)
        x7 = bool(e_ind <= xp.ERROR_INDICATOR_MAX)
        extras = dict(eta_final=solver_out["eta_final"], error_indicator=e_ind, X4_state_passed=x4s,
                      X4_green_passed=x4g, X5_passed=x5, X7_passed=x7, gates_passed=bool(x4s and x4g and x5 and x7))
    except (KeyError, TypeError) as exc:
        raise PilotRecordError(f"solver output for M4 lacks a required field: {exc}") from exc
    return _build("M4", config=config, point=point, inv_dx=inv_dx, n_y=n_y, n_q=n_q, E=E, kx=kx, S=S, n_open=n_open,
                  extras=extras, V_lead_max=V_lead_max, potential_sha256=potential_sha256, n_slices=n_slices,
                  x_start=x_start, wall_time_s=wall_time_s, code_files=code_files)


# ----------------------------------------------------------------------------------------------------------------
# completeness and same-problem validation
# ----------------------------------------------------------------------------------------------------------------
def expected_record_keys(configs) -> tuple:
    """The complete non-M4 protocol record set (145 keys per configuration), in deterministic order."""
    keys = []
    for c in _check_configs(configs):
        for p in POINTS_14:
            for inv in INV_DX:
                for m in GRID_METHODS:
                    keys.append((c, p, *BASE_TRANSVERSE, inv, m))
        for p in B3_POINTS:
            for inv in INV_DX:
                keys.append((c, p, *BASE_TRANSVERSE, inv, "M1"))
        for tv in TRANSVERSE_VARIANTS[1:]:
            for p in NYNQ_POINTS:
                keys.append((c, p, *tv, NYNQ_INV_DX, "M1"))
    return tuple(sorted(keys, key=_order))


def validate_record_set(records, configs) -> dict:
    """Refuse unless the non-M4 records equal the complete expected set for `configs` (no missing, extra, duplicate
    or foreign-config record). Returns {key: record} for the non-M4 records."""
    cs = _check_configs(configs)
    recs = [r for r in records if not (isinstance(r, RGFRecord) and r.method == "M4")]
    sort_records(recs)                                          # type and duplicate check
    foreign = sorted({r.config for r in recs if r.config not in cs})
    if foreign:
        raise PilotRecordError(f"records for unlisted configs {foreign}")
    got = {record_key(r): r for r in recs}
    exp = set(expected_record_keys(cs))
    missing = sorted(exp - set(got), key=_order)
    extra = sorted(set(got) - exp, key=_order)
    if missing or extra:
        raise PilotRecordError(f"incomplete or inconsistent record set: missing {missing}; unexpected {extra}")
    return got


_SAME_FIELDS = ("config", "point", "inv_dx", "n_y", "n_q", "n_open", "potential_sha256", "n_slices")


def check_same_problem(recs) -> None:
    """The given records must describe the same lattice problem: equal config, point, E (bitwise), inv_dx, n_y, n_q,
    n_open, S shape, potential sha256, n_slices and x_start (bitwise)."""
    recs = list(recs)
    if len(recs) < 2:
        raise PilotRecordError("check_same_problem needs at least two records")
    r0 = recs[0]
    for r in recs[1:]:
        for f in _SAME_FIELDS:
            if getattr(r, f) != getattr(r0, f):
                raise PilotRecordError(f"records differ in {f}: {getattr(r0, f)!r} ({r0.method}) vs "
                                       f"{getattr(r, f)!r} ({r.method}) at {record_key(r0)[:5]}")
        for f in ("E", "x_start"):
            if not _same_bits(getattr(r, f), getattr(r0, f)):
                raise PilotRecordError(f"records differ in {f} (bitwise): {getattr(r0, f)!r} vs {getattr(r, f)!r}")
        if r.S.shape != r0.S.shape:
            raise PilotRecordError(f"records differ in S shape: {r0.S.shape} vs {r.S.shape}")


def three_grid_values(records, *, config, method="M1") -> dict:
    """Marshal ONE method's T1, R1 and U at dx = 1/512, 1/1024, 1/2048 (base transverse) into the exact input form of
    b4_pilot_rgf.select_resolution (values) and point_flags (residuals). No classification is done here.
    M1 covers the §1.6 point set (14 points + 3 B3 energies); M2/M3 the 14 points. M4 is refused (never a grid)."""
    if method not in GRID_METHODS:
        raise PilotRecordError(f"three-grid values are defined for {GRID_METHODS} only, not {method!r}")
    recs = [r for r in records if isinstance(r, RGFRecord) and r.config == config]
    got = validate_record_set(recs, [config])
    points = POINT_ORDER if method == "M1" else POINTS_14
    values, residuals = {}, {}
    for p in points:
        rs = [got[(config, p, *BASE_TRANSVERSE, inv, method)] for inv in INV_DX]
        values[(p, "T1")] = tuple(r.T1 for r in rs)
        values[(p, "R1")] = tuple(r.R1 for r in rs)
        residuals[p] = tuple(r.U for r in rs)
    return dict(config=config, method=method, points=tuple(points), values=values, residuals=residuals)


# ----------------------------------------------------------------------------------------------------------------
# raw agreement rows, M4 requirement and adjudication inputs (no classification)
# ----------------------------------------------------------------------------------------------------------------
def _pair_name(pair):
    return f"{pair[0]}-{pair[1]}"


def _expected_rows(configs):
    return [(c, p, inv, _pair_name(pr)) for c in configs for p in POINTS_14 for inv in INV_DX for pr in PAIRS]


def agreement_rows(records, *, configs) -> list:
    """max|S_a - S_b| of the complex S for every pair (M1-M2, M1-M3, M2-M3) at every (config, point, dx) of the 14
    points at base transverse. The record set must be complete for every listed config."""
    cs = _check_configs(configs)
    got = validate_record_set(records, cs)
    rows = []
    for c in cs:
        for p in POINTS_14:
            for inv in INV_DX:
                rs = {m: got[(c, p, *BASE_TRANSVERSE, inv, m)] for m in GRID_METHODS}
                check_same_problem(rs.values())
                for a, b in PAIRS:
                    dS = float(np.max(b3.s_difference(rs[a].S, rs[b].S)))
                    exceeds = bool(dS > S_AGREEMENT)
                    rows.append(dict(config=c, point=p, inv_dx=inv, pair=_pair_name((a, b)), dS=dS, limit=S_AGREEMENT,
                                     exceeds=exceeds, m4_trigger=bool(exceeds and (a, b) in M4_TRIGGER_PAIRS),
                                     classification=CLASSIFICATION))
    return rows


def _check_rows(rows, configs) -> dict:
    cs = _check_configs(configs)
    rows = list(rows)
    if not rows:
        raise PilotRecordError("agreement rows are empty")
    keyed = {}
    for r in rows:
        try:
            k = (r["config"], r["point"], r["inv_dx"], r["pair"])
            dS, exceeds, trig = r["dS"], r["exceeds"], r["m4_trigger"]
        except (KeyError, TypeError) as exc:
            raise PilotRecordError(f"malformed agreement row: {r!r}") from exc
        if k in keyed:
            raise PilotRecordError(f"duplicate agreement row {k}")
        dS = _real(dS, "dS")
        if exceeds is not bool(dS > S_AGREEMENT):
            raise PilotRecordError(f"row {k}: 'exceeds' inconsistent with dS")
        pr = tuple(k[3].split("-"))
        if trig is not bool(exceeds and pr in M4_TRIGGER_PAIRS):
            raise PilotRecordError(f"row {k}: 'm4_trigger' inconsistent with the §1.3 rule")
        keyed[k] = r
    exp = set(_expected_rows(cs))
    missing = sorted(exp - set(keyed), key=repr)
    extra = sorted(set(keyed) - exp, key=repr)
    if missing or extra:
        raise PilotRecordError(f"incomplete agreement rows: missing {missing}; unexpected {extra}")
    return keyed


def _key3_order(k):
    return (CONFIGS.index(k[0]), POINT_ORDER.index(k[1]), INV_DX.index(k[2]))


def m4_required_keys(rows, *, configs) -> tuple:
    """(config, point, inv_dx) keys where M4 is required: the predeclared M4 points at dx = 1/2048 for every config,
    plus every key with an M1-M2 or M1-M3 exceedance (protocol §1.3). Rows must be complete."""
    cs = _check_configs(configs)
    keyed = _check_rows(rows, cs)
    req = {(c, p, M4_PROTOCOL_INV_DX) for c in cs for p in M4_PROTOCOL_POINTS}
    req |= {(k[0], k[1], k[2]) for k, r in keyed.items() if r["m4_trigger"]}
    return tuple(sorted(req, key=_key3_order))


@dataclass(frozen=True)
class AgreementSummary:
    """Facts only. Contains NO step-0 decision: step 0 may be True without adjudication only if all_within_limit;
    any exceedance needs the summarizer's adjudication, and protocol-gap keys have no protocol-mandated M4."""
    configs: tuple
    n_rows: int
    all_within_limit: bool
    exceedances: tuple            # (config, point, inv_dx, pair, dS)
    m4_triggered_keys: tuple      # (config, point, inv_dx)
    protocol_gap_keys: tuple      # M2-M3 exceeds while M1-M2 and M1-M3 do not
    classification: str = CLASSIFICATION


def agreement_summary(rows, *, configs) -> AgreementSummary:
    cs = _check_configs(configs)
    keyed = _check_rows(rows, cs)
    ordered = sorted(keyed, key=lambda k: (_key3_order(k[:3]), [_pair_name(p) for p in PAIRS].index(k[3])))
    exc = tuple((k[0], k[1], k[2], k[3], keyed[k]["dS"]) for k in ordered if keyed[k]["exceeds"])
    trig = tuple(sorted({k[:3] for k in keyed if keyed[k]["m4_trigger"]}, key=_key3_order))
    gap = []
    for c in cs:
        for p in POINTS_14:
            for inv in INV_DX:
                e = {pr: keyed[(c, p, inv, pr)]["exceeds"] for pr in ("M1-M2", "M1-M3", "M2-M3")}
                if e["M2-M3"] and not e["M1-M2"] and not e["M1-M3"]:
                    gap.append((c, p, inv))
    return AgreementSummary(configs=cs, n_rows=len(keyed), all_within_limit=not exc, exceedances=exc,
                            m4_triggered_keys=trig, protocol_gap_keys=tuple(gap))


M4_EVALUATED = "EVALUATED"
M4_NE_PLATFORM = "NOT EVALUABLE: platform"
M4_NE_GATE = "NOT EVALUABLE: M4 gate failed"


def m4_adjudication_inputs(records, platform_report, *, required_keys) -> list:
    """Raw M4 adjudication inputs per required key: eps_m = max|S_m - S_M4| (m = M1, M2, M3, compared in clongdouble)
    and the M4 status. Not evaluable platform: M4 records must be absent and every key is NOT EVALUABLE (never a
    pass). Evaluable platform: the M4 records must equal required_keys exactly. No classification is done here."""
    req = list(required_keys)
    if not req:
        raise PilotRecordError("required_keys is empty")
    if len(set(req)) != len(req):
        raise PilotRecordError("required_keys contains duplicates")
    for k in req:
        if (not isinstance(k, tuple) or len(k) != 3 or k[0] not in CONFIGS or k[1] not in POINTS_14
                or k[2] not in INV_DX):
            raise PilotRecordError(f"invalid required key {k!r}")
    req = sorted(req, key=_key3_order)
    recs = [r for r in records if isinstance(r, RGFRecord)]
    sort_records(recs)
    m4 = {(r.config, r.point, r.inv_dx): r for r in recs if r.method == "M4"}
    if not _platform_evaluable(platform_report):
        if m4:
            raise PilotRecordError("M4 records present although the platform gate is not evaluable")
        return [dict(config=k[0], point=k[1], inv_dx=k[2], m4_status=M4_NE_PLATFORM, eps_M1=None, eps_M2=None,
                     eps_M3=None, classification=CLASSIFICATION) for k in req]
    missing = sorted(set(req) - set(m4), key=_key3_order)
    extra = sorted(set(m4) - set(req), key=_key3_order)
    if missing or extra:
        raise PilotRecordError(f"M4 record set inconsistent with the required keys: missing {missing}; "
                               f"unexpected {extra}")
    grid = {record_key(r): r for r in recs if r.method != "M4"}
    out = []
    for k in req:
        try:
            rs = {m: grid[(k[0], k[1], *BASE_TRANSVERSE, k[2], m)] for m in GRID_METHODS}
        except KeyError as exc:
            raise PilotRecordError(f"M1, M2 and M3 records are all required at M4 key {k}: missing {exc}") from exc
        s4 = m4[k]
        check_same_problem([rs["M1"], rs["M2"], rs["M3"], s4])
        eps = {f"eps_{m}": float(np.max(np.abs(rs[m].S.astype(xp.CLD) - s4.S))) for m in GRID_METHODS}
        status = M4_EVALUATED if s4.extras["gates_passed"] else M4_NE_GATE
        out.append(dict(config=k[0], point=k[1], inv_dx=k[2], m4_status=status, **eps, classification=CLASSIFICATION))
    return out


# ----------------------------------------------------------------------------------------------------------------
# serialization: npz arrays (for Staging.write_npz) and CSV rows (for Staging.write_csv); read-only loaders
# ----------------------------------------------------------------------------------------------------------------
STRING_WIDTHS = {"record_config": 1, "record_point": 8, "record_method": 2, "record_storage": 32,
                 "record_s_sha256": 71, "record_method_fingerprint": 71, "schema_version": 8}
STORAGE = {"M1": "complex128", "M2": "complex128", "M3": "complex128", "M4": "clongdouble_hi_lo"}
INDEX_INT = ("record_inv_dx", "record_n_y", "record_n_q", "record_n_open")
INDEX_FLOAT = ("record_E", "record_kx")


def fixed_str_array(name: str, values) -> np.ndarray:
    """Fixed-width unicode array with the declared width; refuses any value that would be truncated."""
    if name not in STRING_WIDTHS:
        raise PilotRecordError(f"no declared width for {name!r}")
    w = STRING_WIDTHS[name]
    vals = list(values)
    for v in vals:
        if not isinstance(v, str):
            raise PilotRecordError(f"{name}: non-string value {v!r}")
        if len(v) > w:
            raise PilotRecordError(f"{name}: value {v!r} ({len(v)} characters) exceeds the declared width {w}")
    arr = np.array(vals, dtype=f"<U{w}")
    if arr.tolist() != vals:
        raise PilotRecordError(f"{name}: values were truncated")
    return arr


def s_archive_arrays(records) -> dict:
    """Arrays for b4_pilot_io.Staging.write_npz (which adds the classification array). Deterministic record order."""
    recs = sort_records(records)
    if not recs:
        raise PilotRecordError("no records to serialize")
    arrays = {
        "schema_version": fixed_str_array("schema_version", [SCHEMA_VERSION])[0:1].reshape(()),
        "record_config": fixed_str_array("record_config", [r.config for r in recs]),
        "record_point": fixed_str_array("record_point", [r.point for r in recs]),
        "record_method": fixed_str_array("record_method", [r.method for r in recs]),
        "record_storage": fixed_str_array("record_storage", [STORAGE[r.method] for r in recs]),
        "record_s_sha256": fixed_str_array("record_s_sha256", [r.s_sha256 for r in recs]),
        "record_method_fingerprint": fixed_str_array("record_method_fingerprint", [r.method_fingerprint for r in recs]),
        "record_inv_dx": np.array([r.inv_dx for r in recs], dtype=np.int64),
        "record_n_y": np.array([r.n_y for r in recs], dtype=np.int64),
        "record_n_q": np.array([r.n_q for r in recs], dtype=np.int64),
        "record_n_open": np.array([r.n_open for r in recs], dtype=np.int64),
        "record_E": np.array([r.E for r in recs], dtype=np.float64),
        "record_kx": np.array([r.kx for r in recs], dtype=np.float64),
    }
    for i, r in enumerate(recs):
        if r.method == "M4":
            hi, lo = split_hi_lo(r.S)
            arrays[f"r{i:04d}_S_hi"], arrays[f"r{i:04d}_S_lo"] = hi, lo
        else:
            arrays[f"r{i:04d}_S"] = np.array(r.S, dtype=np.complex128, copy=True)
    for k, a in arrays.items():
        if a.dtype == object:
            raise PilotRecordError(f"object array {k!r} refused")
    return arrays


def archive_content_digest(arrays) -> str:
    """sha256 over sorted (name, dtype, shape, bytes). The npz FILE bytes are not reproducible (zip timestamps);
    this digest is the declared reproducibility invariant of the S archive."""
    h = hashlib.sha256()
    for k in sorted(arrays):
        a = np.ascontiguousarray(np.asarray(arrays[k]))
        if a.dtype == object:
            raise PilotRecordError(f"object array {k!r} refused")
        h.update(k.encode() + b"\0" + a.dtype.str.encode() + b"\0" + repr(a.shape).encode() + b"\0" + a.tobytes()
                 + b"\0")
    return "sha256:" + h.hexdigest()


def load_s_archive(path, metadata) -> list:
    """Read-only load of diagnostic_rgf_S.npz (allow_pickle=False). Verifies the classification, schema, index
    arrays, every S hash, shapes and dtypes, and every method fingerprint against the run metadata. Returns one dict
    per record; 'S' is complex128 for M1-M3, clongdouble for M4 only if the platform gate is evaluable (else None,
    with 'S_hi'/'S_lo' always present)."""
    with np.load(path, allow_pickle=False) as z:
        data = {k: z[k] for k in z.files}
    for k in ("classification", "schema_version", *STRING_WIDTHS, *INDEX_INT, *INDEX_FLOAT):
        if k not in data:
            raise PilotRecordError(f"S archive lacks {k!r}")
    if str(data["classification"]) != CLASSIFICATION:
        raise PilotRecordError("S archive classification is not DIAGNOSTIC / NOT B4")
    if data["schema_version"].shape != () or str(data["schema_version"]) != SCHEMA_VERSION:
        raise PilotRecordError(f"S archive schema_version {data['schema_version']!r} != {SCHEMA_VERSION!r}")
    for k, w in STRING_WIDTHS.items():
        if data[k].dtype != np.dtype(f"<U{w}"):
            raise PilotRecordError(f"{k}: dtype {data[k].dtype} != <U{w}")
    for k in INDEX_INT:
        if data[k].dtype != np.int64:
            raise PilotRecordError(f"{k}: dtype {data[k].dtype} != int64")
    for k in INDEX_FLOAT:
        if data[k].dtype != np.float64:
            raise PilotRecordError(f"{k}: dtype {data[k].dtype} != float64")
    idx = [k for k in STRING_WIDTHS if k != "schema_version"] + list(INDEX_INT) + list(INDEX_FLOAT)
    n = data["record_config"].shape[0]
    if n == 0 or any(data[k].shape != (n,) for k in idx):
        raise PilotRecordError("S archive index arrays are empty or of unequal length")
    expected_names = {"classification", "schema_version", *idx}
    evaluable = bool(xp.platform_precision()["evaluable"])
    entries = []
    for i in range(n):
        e = dict(config=str(data["record_config"][i]), point=str(data["record_point"][i]),
                 method=str(data["record_method"][i]), n_y=int(data["record_n_y"][i]), n_q=int(data["record_n_q"][i]),
                 inv_dx=int(data["record_inv_dx"][i]), n_open=int(data["record_n_open"][i]),
                 E=float(data["record_E"][i]), kx=float(data["record_kx"][i]),
                 storage=str(data["record_storage"][i]), s_sha256=str(data["record_s_sha256"][i]),
                 method_fingerprint=str(data["record_method_fingerprint"][i]))
        _check_placement(e["config"], e["point"], e["n_y"], e["n_q"], e["inv_dx"], e["method"])
        if e["storage"] != STORAGE[e["method"]]:
            raise PilotRecordError(f"record {i}: storage {e['storage']!r} invalid for {e['method']}")
        side = 2 * e["n_open"]
        if not 1 <= e["n_open"] <= e["n_y"] - 1:
            raise PilotRecordError(f"record {i}: invalid n_open {e['n_open']}")
        if e["method"] == "M4":
            names = (f"r{i:04d}_S_hi", f"r{i:04d}_S_lo")
            for nm in names:
                if nm not in data or data[nm].dtype != np.complex128 or data[nm].shape != (side, side):
                    raise PilotRecordError(f"record {i}: {nm} missing or of wrong dtype/shape")
            hi, lo = data[names[0]], data[names[1]]
            if _digest_hi_lo(hi, lo) != e["s_sha256"]:
                raise PilotRecordError(f"record {i}: S hash mismatch")
            e["S_hi"], e["S_lo"] = hi, lo
            e["S"] = (hi.astype(xp.CLD) + lo.astype(xp.CLD)) if evaluable else None
        else:
            names = (f"r{i:04d}_S",)
            S = data.get(names[0])
            if S is None or S.dtype != np.complex128 or S.shape != (side, side):
                raise PilotRecordError(f"record {i}: {names[0]} missing or of wrong dtype/shape")
            if s_digest(S) != e["s_sha256"]:
                raise PilotRecordError(f"record {i}: S hash mismatch")
            e["S"] = S
        if not all(np.all(np.isfinite(data[nm])) for nm in names):
            raise PilotRecordError(f"record {i}: non-finite S")
        expected_names.update(names)
        entries.append(e)
    extra = sorted(set(data) - expected_names)
    if extra:
        raise PilotRecordError(f"S archive contains unexpected arrays: {extra}")
    keys = [_entry_key(e) for e in entries]
    if len(set(keys)) != len(keys):
        raise PilotRecordError("S archive contains duplicate record keys")
    if keys != sorted(keys, key=_order):
        raise PilotRecordError("S archive records are not in the deterministic order")
    verify_method_fingerprints(entries, metadata)
    return entries


# ---- CSV -------------------------------------------------------------------------------------------------------
_STR_COLS = ("config", "point", "method", "potential_sha256", "s_sha256", "method_fingerprint")
_INT_COLS = ("n_y", "n_q", "inv_dx", "n_open", "n_closed", "argmin_rcond", "n_slices")
_BOOL_COLS = ("X4_state_passed", "X4_green_passed", "X5_passed", "X7_passed", "gates_passed")
CSV_COLUMNS = (("config", "point", "n_y", "n_q", "inv_dx", "method", "E", "kx", "n_open", "n_closed")
               + OBSERVABLES + ALL_EXTRAS
               + ("V_lead_max", "potential_sha256", "n_slices", "x_start", "wall_time_s", "s_sha256",
                  "method_fingerprint"))
TIMING_COLUMNS = ("wall_time_s",)                             # the only columns excluded from the §5.5 byte check


def _fmt(col, v):
    if col in _STR_COLS:
        return v
    if col in _BOOL_COLS:
        return "True" if v else "False"
    if col in _INT_COLS:
        return str(int(v))
    x = float(v)
    if not math.isfinite(x):
        raise PilotRecordError(f"{col}: non-finite value refused")
    return repr(x)


def points_csv_rows(records) -> tuple:
    """(columns, rows) for b4_pilot_io.Staging.write_csv: one row per record, deterministic order, every value
    already formatted (floats as repr, lossless). Another method's extras are 'NOT APPLICABLE'."""
    recs = sort_records(records)
    if not recs:
        raise PilotRecordError("no records to serialize")
    rows = []
    for r in recs:
        row = {}
        for col in CSV_COLUMNS:
            if col in ALL_EXTRAS:
                row[col] = _fmt(col, r.extras[col]) if col in r.extras else NOT_APPLICABLE
            else:
                row[col] = _fmt(col, getattr(r, col))
        rows.append(row)
    return CSV_COLUMNS, rows


def _parse(col, s):
    if col in _STR_COLS:
        return s
    if col in _BOOL_COLS:
        if s not in ("True", "False"):
            raise PilotRecordError(f"{col}: invalid boolean {s!r}")
        return s == "True"
    if col in _INT_COLS:
        if not _INT_RE.match(s):
            raise PilotRecordError(f"{col}: invalid integer {s!r}")
        return int(s)
    try:
        x = float(s)
    except ValueError as exc:
        raise PilotRecordError(f"{col}: invalid float {s!r}") from exc
    if not math.isfinite(x) or repr(x) != s:
        raise PilotRecordError(f"{col}: {s!r} is not a canonical finite float")
    return x


def parse_points_csv(path, metadata) -> list:
    """Read-only parse of diagnostic_rgf_points.csv: label line, exact header, canonical values, NOT APPLICABLE
    pattern per method, placement, deterministic order and method fingerprints."""
    with open(path, "r", encoding="utf-8", newline="") as f:
        first = f.readline().rstrip("\n").rstrip("\r")
        if first != CSV_HEADER:
            raise PilotRecordError("points CSV lacks the DIAGNOSTIC / NOT B4 label line")
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None or tuple(header) != CSV_COLUMNS:
            raise PilotRecordError("points CSV header does not match the schema")
        out = []
        for ln, row in enumerate(reader, start=3):
            if len(row) != len(CSV_COLUMNS):
                raise PilotRecordError(f"line {ln}: {len(row)} fields, expected {len(CSV_COLUMNS)}")
            raw = dict(zip(CSV_COLUMNS, row))
            m = raw["method"]
            if m not in METHODS:
                raise PilotRecordError(f"line {ln}: unknown method {m!r}")
            e = {}
            for col in CSV_COLUMNS:
                if col in ALL_EXTRAS and col not in EXTRAS[m]:
                    if raw[col] != NOT_APPLICABLE:
                        raise PilotRecordError(f"line {ln}: {col} must be {NOT_APPLICABLE!r} for {m}")
                    continue
                e[col] = _parse(col, raw[col])
            _check_placement(e["config"], e["point"], e["n_y"], e["n_q"], e["inv_dx"], m)
            for col in ("potential_sha256", "s_sha256", "method_fingerprint"):
                _sha(e[col], col)
            out.append(e)
    if not out:
        raise PilotRecordError("points CSV has no rows")
    keys = [_entry_key(e) for e in out]
    if len(set(keys)) != len(keys):
        raise PilotRecordError("points CSV contains duplicate record keys")
    if keys != sorted(keys, key=_order):
        raise PilotRecordError("points CSV rows are not in the deterministic order")
    verify_method_fingerprints(out, metadata)
    return out


def verify_archive_csv_consistency(archive_entries, csv_entries) -> None:
    """The S archive and the points CSV must describe the same records: same ordered keys, s_sha256, n_open, E, kx and
    method fingerprint; and every observable in the CSV must equal the value recomputed from the archived S
    (bitwise; M4 only where the clongdouble S could be reconstructed)."""
    ka = [_entry_key(e) for e in archive_entries]
    kc = [_entry_key(e) for e in csv_entries]
    if ka != kc:
        raise PilotRecordError("S archive and points CSV have different record keys or order")
    for a, c in zip(archive_entries, csv_entries):
        for f in ("s_sha256", "n_open", "method_fingerprint"):
            if a[f] != c[f]:
                raise PilotRecordError(f"{_entry_key(a)}: {f} differs between archive and CSV")
        for f in ("E", "kx"):
            if not _same_bits(a[f], c[f]):
                raise PilotRecordError(f"{_entry_key(a)}: {f} differs between archive and CSV")
        if a["S"] is not None:
            obs = _observables(a["S"], a["n_open"])
            for k in OBSERVABLES:
                if not _same_bits(obs[k], c[k]):
                    raise PilotRecordError(f"{_entry_key(a)}: CSV {k} = {c[k]!r} != recomputed {obs[k]!r}")
