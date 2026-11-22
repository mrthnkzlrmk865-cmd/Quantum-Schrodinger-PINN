"""B4 limited diagnostic pilot — output guard, labelling, metadata and atomic publication (phase T5).

Classification: DIAGNOSTIC / NOT B4.

Implements docs/B4_PILOT_PROTOCOL.md §5.1–§5.5 infrastructure and the T4-plan rules R1/R2 and §3.1 (G0 atomicity):

* every pilot file is written first into a STAGING directory OUTSIDE the repository;
* results/b4_pilot/ is created only by `publish`, and only with a passed G0 outcome; it is never created or populated
  before that, and a failed G0 (or any exception) leaves no pilot output in the repository;
* the publication target must be exactly <repo>/results/b4_pilot, which is checked against the §5.3 allowlist/denylist;
* every CSV/JSON/NPZ file carries the label DIAGNOSTIC / NOT B4 and every file name starts with "diagnostic_";
* G0 helpers: the protocol SHA-256 is compared with a PINNED literal reference (never derived from the file being
  checked), the A6 baseline commit is a separate identifier, and the run aborts on modified or untracked files.

This module performs no physical calculation and never executes G0 or G-REAL; it only provides the guarded machinery.
A G0Outcome is a record supplied by the (future) G0 runner; constructing one here does not mean G0 has passed.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import uuid
from dataclasses import dataclass, field

import numpy as np

CLASSIFICATION = "DIAGNOSTIC / NOT B4"
CSV_HEADER = "# DIAGNOSTIC / NOT B4 — docs/B4_PILOT_PROTOCOL.md"
FILE_PREFIX = "diagnostic_"

PROTOCOL_PATH = "docs/B4_PILOT_PROTOCOL.md"
# G0 protocol reference (plan R1): SHA-256 of the bytes of the committed protocol file (last changed in commit f3ab9ac).
# PINNED LITERAL — never recomputed from the file being validated.
PROTOCOL_SHA256_REFERENCE = "sha256:3cd620d2affc6ba5e59e6bd26bdd407afdc6eed7269ee91dd5448a30a0a035ab"
# A6 baseline commit for the protected-artifact comparison (decision A6). A different identifier: not a file hash.
A6_BASELINE_COMMIT = "f3ab9accf04e8a4927fac3f276cd88910ac1bdd3"

PILOT_DIR = "results/b4_pilot"                                   # §5.1 publication target (the only allowed one)
FORBIDDEN_RESULT_DIRS = (                                        # §5.3 denylist (explicit + the rgf_b*_validation family)
    "results/stationary_validation", "results/rgf_b4_validation", "results/transport", "results/validation",
    "results/rgf_b1_validation", "results/rgf_b2_validation", "results/rgf_b3_validation",
    "results/rgf_b3x_validation", "results/rgf_b3y_validation",
)
PARTIAL_PREFIX = ".b4_pilot.partial-"                            # hidden same-filesystem staging inside results/

# §5.4 modules whose per-file SHA-256 is recorded (decision A1: pilot-local, under provenance.code_fingerprint)
CODE_FILES = (
    "stationary.py", "smatrix_spec.py", "sparse_reference.py", "sparse_wavefunction.py", "sparse_exact_lattice.py",
    "reference_solver.py", "potential.py", "disorder.py", "config.py", "initial_condition.py",
    "validate_rgf_b2.py", "validate_rgf_b3.py", "b4_pilot_rgf.py", "b4_pilot_io.py",
)


class PilotGuardError(RuntimeError):
    """A guard refused an action (path, label, tree state, protocol hash or publication)."""


# ----------------------------------------------------------------------------------------------------------------
# hashing and repository state
# ----------------------------------------------------------------------------------------------------------------
def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def sha256_file(path: str) -> str:
    with open(path, "rb") as f:
        return sha256_bytes(f.read())


def _git(repo: str, *args: str, input_bytes: bytes | None = None) -> bytes:
    try:
        r = subprocess.run(["git", "-C", repo, *args], input=input_bytes, capture_output=True, check=False)
    except FileNotFoundError as exc:
        raise PilotGuardError(f"git is required: {exc}") from exc
    if r.returncode != 0:
        raise PilotGuardError(f"git {' '.join(args)} failed: {r.stderr.decode(errors='replace').strip()}")
    return r.stdout


def head_commit(repo: str) -> str:
    return _git(repo, "rev-parse", "HEAD").decode().strip()


def dirty_entries(repo: str) -> list:
    """Every modified, staged, deleted or untracked (non-ignored) path, from `git status --porcelain`."""
    out = _git(repo, "status", "--porcelain", "--untracked-files=all").decode()
    return [line for line in out.splitlines() if line.strip()]


def require_clean_tree(repo: str) -> str:
    """Plan R2: abort on any modified or untracked repository file. Returns HEAD on success."""
    dirty = dirty_entries(repo)
    if dirty:
        raise PilotGuardError("repository is not clean (commit all pilot code and tests first): " + "; ".join(dirty))
    return head_commit(repo)


def check_protocol(repo: str, reference: str = PROTOCOL_SHA256_REFERENCE) -> str:
    """Plan R1: the working-tree protocol AND the protocol at HEAD must both hash to the pinned reference."""
    p = os.path.join(repo, PROTOCOL_PATH)
    if not os.path.isfile(p):
        raise PilotGuardError(f"protocol file missing: {PROTOCOL_PATH}")
    wt = sha256_file(p)
    head = sha256_bytes(_git(repo, "show", f"HEAD:{PROTOCOL_PATH}"))
    if wt != reference or head != reference:
        raise PilotGuardError(f"protocol SHA-256 mismatch: working tree {wt}, HEAD {head}, reference {reference}")
    return reference


def code_fingerprint(repo: str, files=CODE_FILES) -> dict:
    """Decision A1: per-file SHA-256 (pilot-local; provenance.py is not modified). Missing files raise."""
    per_file, h = {}, hashlib.sha256()
    for name in files:
        p = os.path.join(repo, name)
        if not os.path.isfile(p):
            raise PilotGuardError(f"code file missing for fingerprint: {name}")
        with open(p, "rb") as f:
            data = f.read()
        per_file[name] = sha256_bytes(data)
        h.update(name.encode() + b"\0" + data + b"\0")
    return {"code_fingerprint": "sha256:" + h.hexdigest(), "code_files": per_file}


def protected_snapshot(repo: str, paths) -> dict:
    """SHA-256 of each protected path in the working tree (missing -> None)."""
    out = {}
    for p in paths:
        full = os.path.join(repo, p)
        out[p] = sha256_file(full) if os.path.isfile(full) else None
    return out


def baseline_hashes(repo: str, paths, commit: str = A6_BASELINE_COMMIT) -> dict:
    """SHA-256 of each path's bytes as represented in the A6 baseline commit (git cat-file --batch; missing -> None)."""
    paths = list(paths)
    raw = _git(repo, "cat-file", "--batch", input_bytes="".join(f"{commit}:{p}\n" for p in paths).encode())
    out, pos = {}, 0
    for p in paths:
        nl = raw.index(b"\n", pos)
        header = raw[pos:nl].decode()
        pos = nl + 1
        if header.endswith(" missing"):
            out[p] = None
            continue
        size = int(header.split()[2])
        out[p] = sha256_bytes(raw[pos:pos + size])
        pos += size + 1
    return out


def compare_snapshots(expected: dict, actual: dict) -> list:
    """Human-readable differences (missing on either side, or changed). Empty list = identical."""
    problems = []
    for p in sorted(set(expected) | set(actual)):
        e, a = expected.get(p, "<absent>"), actual.get(p, "<absent>")
        if e != a:
            problems.append(f"{p}: expected {e}, found {a}")
    return problems


# ----------------------------------------------------------------------------------------------------------------
# path guards (§5.3)
# ----------------------------------------------------------------------------------------------------------------
def _real(p: str) -> str:
    return os.path.realpath(os.path.abspath(p))


def _within(path: str, root: str) -> bool:
    return path == root or path.startswith(root + os.sep)


def pilot_target(repo: str) -> str:
    return os.path.join(_real(repo), PILOT_DIR)


def check_output_path(path: str, repo: str) -> str:
    """Allowlist: the resolved path must lie in <repo>/results/b4_pilot. Denylist: never in a §5.3 official directory
    (also catches symlinks resolving into them)."""
    rp, repo_r = _real(path), _real(repo)
    for d in FORBIDDEN_RESULT_DIRS:
        if _within(rp, os.path.join(repo_r, d)):
            raise PilotGuardError(f"refusing to write into the official result directory {d}")
    if not _within(rp, os.path.join(repo_r, PILOT_DIR)):
        raise PilotGuardError(f"refusing to write outside {PILOT_DIR}: {path}")
    return rp


def check_staging_path(path: str, repo: str) -> str:
    """Pre-G0 staging and scratch must be OUTSIDE the repository (plan §3.1)."""
    rp = _real(path)
    if _within(rp, _real(repo)):
        raise PilotGuardError(f"staging/scratch must be outside the repository: {path}")
    return rp


# ----------------------------------------------------------------------------------------------------------------
# labelled writers (staging only)
# ----------------------------------------------------------------------------------------------------------------
def _check_name(name: str, ext: str) -> None:
    base = os.path.basename(name)
    if not base.startswith(FILE_PREFIX) or not base.endswith(ext):
        raise PilotGuardError(f"pilot file names must start with {FILE_PREFIX!r} and end with {ext!r}: {name}")


@dataclass
class G0Outcome:
    """Record produced by the (future) G0 runner. passed=True is required for publication; protocol_sha256 must be the
    pinned reference. Constructing this object does not execute or pass G0."""
    passed: bool
    protocol_sha256: str
    head: str
    checks: dict = field(default_factory=dict)


class Staging:
    """Staging directory OUTSIDE the repository. All pilot writes go here; nothing reaches the repository until
    `publish` with a passed G0 outcome. Use as a context manager: an exception discards the staging directory."""

    def __init__(self, repo: str, scratch_root: str | None = None):
        self.repo = _real(repo)
        root = check_staging_path(scratch_root or tempfile.gettempdir(), self.repo)
        self.path = check_staging_path(tempfile.mkdtemp(prefix="b4_pilot_staging_", dir=root), self.repo)
        self.published = False
        self.discarded = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is not None or not self.published:
            self.discard()
        return False

    def _target(self, rel: str) -> str:
        if self.discarded or self.published:
            raise PilotGuardError("staging is closed")
        full = _real(os.path.join(self.path, rel))
        if not _within(full, self.path):
            raise PilotGuardError(f"path escapes the staging directory: {rel}")
        os.makedirs(os.path.dirname(full), exist_ok=True)
        return full

    def write_csv(self, rel: str, columns, rows) -> str:
        _check_name(rel, ".csv")
        p = self._target(rel)
        with open(p, "w", newline="", encoding="utf-8") as f:
            f.write(CSV_HEADER + "\n")
            w = csv.DictWriter(f, fieldnames=list(columns))
            w.writeheader()
            for r in rows:
                w.writerow(r)
        return p

    def write_json(self, rel: str, obj: dict) -> str:
        _check_name(rel, ".json")
        if not isinstance(obj, dict):
            raise PilotGuardError("JSON payload must be a dict")
        if obj.get("classification", CLASSIFICATION) != CLASSIFICATION:
            raise PilotGuardError("JSON classification must be DIAGNOSTIC / NOT B4")
        p = self._target(rel)
        with open(p, "w", encoding="utf-8") as f:
            json.dump(dict(obj, classification=CLASSIFICATION), f, indent=2, sort_keys=True, default=str)
        return p

    def write_npz(self, rel: str, **arrays) -> str:
        _check_name(rel, ".npz")
        if "classification" in arrays:
            raise PilotGuardError("the classification array is added by the writer")
        p = self._target(rel)
        np.savez_compressed(p, classification=np.array(CLASSIFICATION), **arrays)
        return p

    def write_readme(self, text: str) -> str:
        p = self._target("README.md")
        with open(p, "w", encoding="utf-8") as f:
            f.write(f"# {CLASSIFICATION}\n\nProtocol: {PROTOCOL_PATH} ({PROTOCOL_SHA256_REFERENCE})\n\n{text}\n")
        return p

    def discard(self) -> None:
        if not self.discarded and not self.published:
            shutil.rmtree(self.path, ignore_errors=True)
            self.discarded = True


def _verify_labels(root: str) -> None:
    for d, _, files in os.walk(root):
        for f in files:
            p = os.path.join(d, f)
            if f == "README.md":
                with open(p, encoding="utf-8") as fh:
                    ok = CLASSIFICATION in fh.readline()
            elif f.endswith(".csv") and f.startswith(FILE_PREFIX):
                with open(p, encoding="utf-8") as fh:
                    ok = fh.readline().rstrip("\n") == CSV_HEADER
            elif f.endswith(".json") and f.startswith(FILE_PREFIX):
                with open(p, encoding="utf-8") as fh:
                    ok = json.load(fh).get("classification") == CLASSIFICATION
            elif f.endswith(".npz") and f.startswith(FILE_PREFIX):
                with np.load(p, allow_pickle=False) as z:
                    ok = "classification" in z.files and str(z["classification"]) == CLASSIFICATION
            else:
                ok = False
            if not ok:
                raise PilotGuardError(f"unlabelled or non-pilot file in staging: {os.path.relpath(p, root)}")


def publish(staging: Staging, outcome: G0Outcome, metadata: dict) -> str:
    """Atomically publish the staging contents as <repo>/results/b4_pilot. Requires a passed G0 outcome with the pinned
    protocol reference, consistent metadata, an absent target and labelled files only. On any refusal the staging
    directory is discarded and the repository is left without pilot output."""
    try:
        if not isinstance(outcome, G0Outcome) or outcome.passed is not True:
            raise PilotGuardError("G0 has not passed: nothing is published")
        if outcome.protocol_sha256 != PROTOCOL_SHA256_REFERENCE:
            raise PilotGuardError("G0 outcome does not carry the pinned protocol reference")
        if metadata.get("protocol_sha256") != PROTOCOL_SHA256_REFERENCE or metadata.get("git_commit") != outcome.head:
            raise PilotGuardError("metadata inconsistent with the G0 outcome")
        target = pilot_target(staging.repo)
        check_output_path(target, staging.repo)
        results = os.path.dirname(target)
        if os.path.lexists(target):
            raise PilotGuardError(f"{PILOT_DIR} already exists; refusing to overwrite or merge")
        leftovers = [n for n in os.listdir(results) if n.startswith(PARTIAL_PREFIX)] if os.path.isdir(results) else []
        if leftovers:
            raise PilotGuardError(f"leftover partial publication directories in results/: {leftovers}")
        staging.write_json("diagnostic_run_metadata.json", metadata)
        _verify_labels(staging.path)
        os.makedirs(results, exist_ok=True)
        partial = os.path.join(results, PARTIAL_PREFIX + uuid.uuid4().hex)
        try:
            shutil.copytree(staging.path, partial)
            os.rename(partial, target)                 # atomic on one filesystem: target appears complete or not at all
        finally:
            if os.path.lexists(partial):
                shutil.rmtree(partial, ignore_errors=True)
        staging.published = True
        shutil.rmtree(staging.path, ignore_errors=True)
        return target
    except Exception:
        staging.discard()
        raise


# ----------------------------------------------------------------------------------------------------------------
# metadata (§5.4)
# ----------------------------------------------------------------------------------------------------------------
def build_metadata(repo: str, *, head: str, protected_before: dict, protected_after: dict | None = None,
                   code_files=CODE_FILES, command: str | None = None, extra: dict | None = None) -> dict:
    """§5.4 metadata. head must come from require_clean_tree (plan R2). Environment fields are read-only."""
    env = dict(python=platform.python_version(), platform=platform.platform(), machine=platform.machine(),
               cpu_count=os.cpu_count(), numpy=np.__version__,
               OPENBLAS_NUM_THREADS=os.environ.get("OPENBLAS_NUM_THREADS"))
    try:
        import scipy
        env["scipy"] = scipy.__version__
    except Exception:  # noqa: BLE001
        env["scipy"] = None
    meta = {
        "classification": CLASSIFICATION,
        "protocol_path": PROTOCOL_PATH,
        "protocol_sha256": PROTOCOL_SHA256_REFERENCE,
        "a6_baseline_commit": A6_BASELINE_COMMIT,
        "git_commit": head,
        "clean_tree": True,
        "provenance": {"code_fingerprint": code_fingerprint(repo, code_files)},
        "environment": env,
        "protected_before": protected_before,
        "protected_after": protected_after,
        "command": command if command is not None else " ".join(sys.argv),
    }
    if extra:
        overlap = set(extra) & set(meta)
        if overlap:
            raise PilotGuardError(f"extra metadata may not override core fields: {sorted(overlap)}")
        meta.update(extra)
    return meta
