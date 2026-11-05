"""
Provenance and cache validation for experiment outputs.

Every cached result (an .npz file or a PINN run directory) carries a
metadata record. The record contains:

  spec              the scientifically relevant inputs: the relevant Config
                    fields plus run parameters (end time, number of saved
                    frames, explicit seed list, number of realizations, ...)
  fingerprint       "sha256:<hex>" of the canonical JSON of `spec`
  code_fingerprint  "sha256:<hex>" of the source files that determine the
                    numerical result, plus per-file hashes
  git_commit / git_dirty, created_utc, library versions, platform,
                    CPU count, wall time (not part of the fingerprint)

Cache policy (check_cache):
  * no file                       -> "miss" (compute)
  * fingerprint matches           -> "hit" (reuse). If the code fingerprint
                                     differs, a warning is printed; with
                                     strict_code=True this is an error.
  * fingerprint differs, or the   -> CacheMismatchError, unless
    file has no metadata (legacy)    recompute=True ("recompute": the file
                                     will be overwritten by the caller).
A result produced with a different configuration is therefore never
accepted silently, and existing results are never overwritten unless the
caller asks for it explicitly.

Fields that do not change the numbers (device, log interval, batch size
for propagating realizations together, FFT worker count) are recorded but
excluded from the fingerprint. The PyTorch intra-op thread count DOES change
PINN training results (reduction order; docs/REPRODUCIBILITY_AUDIT.md) and is
part of the PINN spec.

Environment: an environment fingerprint (Python, NumPy, SciPy, [PyTorch],
CPU architecture) is stored with every result. A cache produced in a
different numerical environment is reused with a warning by default, and
refused with strict_env. Scientific runs should use both strict_code and
strict_env (command-line flag --strict).
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import math
import os
import platform
import subprocess
import warnings

REPO = os.path.dirname(os.path.abspath(__file__))
SCHEMA_VERSION = 2   # 2: environment fingerprint, PINN thread count in spec

# Source files whose content determines the numbers of each experiment kind.
_CORE = ["config.py", "disorder.py", "potential.py", "initial_condition.py", "reference_solver.py",
         "observables.py"]
CODE_FILES = {
    "transport": _CORE + ["transport.py", "experiment_disorder.py"],
    "localization": _CORE + ["experiment_localization.py"],
    "pinn": _CORE + ["embeddings.py", "model.py", "pinn_core.py", "train.py", "evaluate.py",
                     "experiment_pinn.py"],
    "validation": _CORE + ["analytic.py", "transport.py", "validate_solver.py"],
}

# Config fields that never influence a reference-solver result.
_PINN_ONLY_FIELDS = {
    "hidden_layers", "hidden_units", "use_fourier", "fourier_num_freq", "fourier_max_freq",
    "n_interior", "n_norm_side", "n_norm_times", "sampling", "rad_pool_factor", "epochs_adam",
    "lr_adam", "lr_decay_every", "lr_decay_gamma", "use_lbfgs", "epochs_lbfgs", "loss_weights",
    "pde_warmup_epochs", "grad_clip_norm", "resample_every", "log_every", "seed", "device",
    "rad_resample_every", "lbfgs_points_factor", "causal_eps", "causal_bins",
}
_NON_NUMERICAL_FIELDS = {"device", "log_every"}


class CacheMismatchError(RuntimeError):
    pass


# ----------------------------------------------------------------------
def _normalize(obj):
    """Make an object JSON-canonical: tuples -> lists, numpy scalars ->
    Python, floats exact (repr round-trips), non-finite floats as strings."""
    try:
        import numpy as np
        if isinstance(obj, np.generic):
            obj = obj.item()
        elif isinstance(obj, np.ndarray):
            obj = obj.tolist()
    except ImportError:  # pragma: no cover
        pass
    if isinstance(obj, dict):
        return {str(k): _normalize(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, range)):
        return [_normalize(v) for v in obj]
    if isinstance(obj, bool) or obj is None or isinstance(obj, str):
        return obj
    if isinstance(obj, (int, float)):
        # numbers are compared by value: 100 and 100.0 are the same input
        # (exact for |int| < 2**53, which covers seeds and grid sizes)
        obj = float(obj)
        return obj if math.isfinite(obj) else repr(obj)
    raise TypeError(f"cannot fingerprint object of type {type(obj)}")


def canonical_json(obj) -> str:
    return json.dumps(_normalize(obj), sort_keys=True, separators=(",", ":"))


def fingerprint(spec: dict) -> str:
    return "sha256:" + hashlib.sha256(canonical_json(spec).encode()).hexdigest()


def solver_config_dict(cfg) -> dict:
    """Config fields relevant to reference-solver (non-PINN) experiments."""
    return {k: v for k, v in cfg.to_dict().items() if k not in _PINN_ONLY_FIELDS}


def pinn_config_dict(cfg) -> dict:
    return {k: v for k, v in cfg.to_dict().items() if k not in _NON_NUMERICAL_FIELDS}


def code_fingerprint(kind: str) -> dict:
    files = {}
    h = hashlib.sha256()
    for name in CODE_FILES[kind]:
        with open(os.path.join(REPO, name), "rb") as f:
            data = f.read()
        files[name] = "sha256:" + hashlib.sha256(data).hexdigest()
        h.update(name.encode() + b"\0" + data + b"\0")
    return {"code_fingerprint": "sha256:" + h.hexdigest(), "code_files": files}


def _git(*args):
    try:
        return subprocess.run(["git", *args], cwd=REPO, capture_output=True, text=True,
                              timeout=10).stdout.strip()
    except Exception:  # noqa: BLE001 - git is optional
        return None


# Libraries whose versions can change the floating-point results of each
# experiment kind (the reference solver does not use PyTorch).
ENV_PACKAGES = {
    "transport": ("python", "numpy", "scipy", "machine"),
    "localization": ("python", "numpy", "scipy", "machine"),
    "validation": ("python", "numpy", "scipy", "machine"),
    "pinn": ("python", "numpy", "scipy", "torch", "machine"),
}


def environment_info() -> dict:
    info = dict(python=platform.python_version(), platform=platform.platform(),
                machine=platform.machine(), cpu_count=os.cpu_count())
    for mod in ("numpy", "scipy", "torch"):
        try:
            info[mod] = __import__(mod).__version__
        except Exception:  # noqa: BLE001
            info[mod] = None
    commit = _git("rev-parse", "HEAD")
    status = _git("status", "--porcelain", "--untracked-files=no")
    info["git_commit"] = commit or None
    info["git_dirty"] = None if status is None else bool(status)
    return info


def _machine_of(env: dict):
    """CPU architecture of a recorded environment. Records written before
    schema 2 have no 'machine' entry; it is then read from the platform
    string (e.g. 'Linux-...-x86_64-with-glibc2.39')."""
    if env.get("machine"):
        return env["machine"]
    plat = env.get("platform") or ""
    for arch in ("x86_64", "aarch64", "arm64", "AMD64", "ppc64le", "s390x"):
        if arch in plat:
            return arch
    return None


def environment_fingerprint(kind: str, env: dict | None = None) -> dict:
    """Fingerprint of the numerically relevant environment for `kind`:
    the versions of the libraries that do the arithmetic, the Python
    version and the CPU architecture. Not included (recorded only): OS
    kernel, CPU model, BLAS/FFT build details, which cannot be read
    portably."""
    env = environment_info() if env is None else env
    fields = {k: (_machine_of(env) if k == "machine" else env.get(k)) for k in ENV_PACKAGES[kind]}
    return {"env_fingerprint": fingerprint(fields), "env_fields": fields}


def make_metadata(kind: str, spec: dict, **extra) -> dict:
    meta = dict(schema_version=SCHEMA_VERSION, kind=kind, spec=_normalize(spec),
                fingerprint=fingerprint(spec),
                created_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"))
    meta.update(code_fingerprint(kind))
    meta["environment"] = environment_info()
    meta.update(environment_fingerprint(kind, meta["environment"]))
    meta.update(_normalize(extra))
    return meta


# ----------------------------------------------------------------------
def read_npz_metadata(path):
    import numpy as np
    with np.load(path, allow_pickle=False) as d:
        if "metadata_json" not in d.files:
            return None
        return json.loads(str(d["metadata_json"]))


def read_dir_metadata(run_dir):
    p = os.path.join(run_dir, "run_metadata.json")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def write_dir_metadata(run_dir, meta):
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, "run_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2, sort_keys=True)


def check_cache(stored_meta, exists: bool, spec: dict, kind: str, where: str,
                recompute=False, strict_code=False, strict_env=False) -> str:
    """Return "miss", "hit" or "recompute"; raise CacheMismatchError for
    an incompatible cache unless recompute=True.

    Configuration mismatch -> always an error. Source-code or numerical-
    environment mismatch -> warning by default, error with strict_code /
    strict_env (use both, i.e. --strict, for results that enter a report)."""
    if not exists:
        return "miss"
    want = fingerprint(spec)
    if stored_meta is None:
        if recompute:
            return "recompute"
        raise CacheMismatchError(
            f"{where}: existing result has no provenance metadata (produced before cache "
            f"validation existed). It cannot be verified against the requested configuration. "
            f"Use a fresh output directory, or pass --recompute to overwrite it.")
    have = stored_meta.get("fingerprint")
    if have != want:
        if recompute:
            return "recompute"
        diff = _spec_diff(stored_meta.get("spec", {}), _normalize(spec))
        raise CacheMismatchError(
            f"{where}: cached result was produced with a different configuration "
            f"(stored {have}, requested {want}). Differing entries: {diff}. "
            f"Use a fresh output directory, or pass --recompute to overwrite it.")
    code_now = code_fingerprint(kind)["code_fingerprint"]
    if stored_meta.get("code_fingerprint") != code_now:
        msg = (f"{where}: configuration matches but the source code that produced the cached "
               f"result has changed ({stored_meta.get('code_fingerprint')} -> {code_now}).")
        if strict_code:
            if recompute:
                return "recompute"
            raise CacheMismatchError(msg + " Refusing to reuse it (strict_code).")
        warnings.warn(msg + " Reusing it; pass --strict-code to forbid this.", stacklevel=2)
    stored_env = stored_meta.get("env_fields")
    if stored_env is None:                     # schema-1 record: derive from recorded versions
        stored_env = environment_fingerprint(kind, stored_meta.get("environment", {}))["env_fields"]
    env_now = environment_fingerprint(kind)["env_fields"]
    if _normalize(stored_env) != _normalize(env_now):
        diff = _spec_diff(_normalize(stored_env), _normalize(env_now))
        msg = (f"{where}: configuration matches but the numerical environment differs: {diff}.")
        if strict_env:
            if recompute:
                return "recompute"
            raise CacheMismatchError(msg + " Refusing to reuse it (strict_env).")
        warnings.warn(msg + " Reusing it; pass --strict-env to forbid this.", stacklevel=2)
    return "hit"


def _spec_diff(a, b, prefix=""):
    out = []
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            out += _spec_diff(a.get(k, "<missing>"), b.get(k, "<missing>"), f"{prefix}{k}.")
    elif a != b:
        sa, sb = json.dumps(a), json.dumps(b)
        if len(sa) > 60:
            sa = sa[:57] + "..."
        if len(sb) > 60:
            sb = sb[:57] + "..."
        out.append(f"{prefix[:-1]}: {sa} -> {sb}")
    return out[:10]
