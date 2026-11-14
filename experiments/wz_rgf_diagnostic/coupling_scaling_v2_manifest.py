r"""Manifest generator for the W/Z coupling-scaling experiment v2 (PREDICTIONS_coupling_scaling_v2.md, Revision 3.7).

Scope: §2.3 draws, the §2.4 acceptance filter (energy rule, grid conditioning, continuous certificate on [0, 1] with the
common Taylor radius), the §2.4 item 5 redraw cap, §2.6 index order, and the §8.1 manifest files. Filter and certificate
steps ONLY (§8.2 phase 2): this module never runs a W/Z solver, the production RGF solver, or a main-experiment mpmath
reference.

Solver protection:
  * imports: numpy, the production certificate/covering modules and `stationary.transverse_basis` (E_perp source values,
    §2.1). No import of wz_rgf, wz_mutants, wz_variants, coupled_blocks, the driver, the oracle or the transcription
    (AST-checked in tests);
  * at generation time every other callable of the `stationary` module is replaced by a raising stub (`solver_guard`),
    and the absence of the solver modules from sys.modules is asserted before and after generation.

Implementer choices not fixed by the frozen spec (flagged for R0, see R0_REVIEW_PACKET.md Part I):
  * C-1 grid conditioning (§2.4 item 2) is evaluated in binary64 (numpy SVD, `numpy.linalg.cond`, 2-norm) on the exact
    L(s) / R_W(s) matrices rounded once to complex128; inverses (prefix Ĝ, Ĝ_NN) by `numpy.linalg.inv`. A non-finite
    value counts as > 1e3;
  * C-2 E_upper for n_open = N_y is E_perp,1 (double, imported exactly) + 64·(2 − δ_min) in exact rationals;
  * C-3 filter steps are evaluated in the order energy → grid conditioning → certificate, stopping at the first failure
    (the filter is a conjunction, so acceptance does not depend on the order; only the logged reason does);
  * C-4 Λ in the M10 matrix Z − (Λ − I) is the lattice's own lead root (enclosure midpoint), as assembled in L / R_W.
"""
from __future__ import annotations

import hashlib
import json
import os
import platform
import sys
from fractions import Fraction

import numpy as np

import coupling_scaling_v2_certify as C
import coupling_scaling_v2_covering as COV

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT_DIR = os.path.join(HERE, "results_coupling_scaling_v2")

ROOTS = {"XP": 202611110001, "HO": 202611110002}
BOOTSTRAP_ROOT = 202611110003
NYS = (2, 3, 4)
CLASSES = ("R", "H", "D")
NS = (2, 3, 4, 6)
DRAWS = {"XP": {"R": 2, "H": 2, "D": 1}, "HO": {"R": 4, "H": 4, "D": 2}}
V0 = 8
DX = 1.0 / 8
DX2 = Fraction(1, 64)
DELTA_MIN = Fraction(1, 50)
COND_MAX = 1e3
MAX_ATTEMPTS = 10000
GRID = COV.GRID                                  # the 81 §3 points, increasing
FORBIDDEN_MODULES = ("wz_rgf", "wz_mutants", "wz_variants", "coupled_blocks", "coupling_scaling_v2",
                     "wz_mutants_oracle", "wz_mutants_mp")
SCHEMA = "coupling_scaling_v2_manifest/1"


class ManifestAbort(RuntimeError):
    """§2.4 item 5: redraw cap reached; no partial manifest is written."""


class SolverGuardViolation(RuntimeError):
    pass


# --------------------------------------------------------------------------------------------- E_perp (read-only)
def _transverse_basis():
    if ROOT_DIR not in sys.path:
        sys.path.insert(0, ROOT_DIR)
    import stationary                                       # noqa: E402  (E_perp source values only)
    return stationary


def eperp_doubles(ny, _cache={}):
    if ny not in _cache:
        st = _transverse_basis()
        _cache[ny] = tuple(float(e) for e in st.transverse_basis(1.0, ny).E_perp)
    return _cache[ny]


class solver_guard:
    """Context manager: every callable attribute of `stationary` except transverse_basis-derived data is replaced by a
    raising stub, and the solver modules must not be imported. Restores the module on exit."""

    def __enter__(self):
        for ny in NYS:
            eperp_doubles(ny)                               # E_perp fixed before the guard
        self._check_modules()
        st = _transverse_basis()
        self._st, self._saved = st, {}
        for name in dir(st):
            obj = getattr(st, name)
            if callable(obj) and getattr(obj, "__module__", None) == st.__name__:
                self._saved[name] = obj
                setattr(st, name, self._stub(name))
        return self

    @staticmethod
    def _stub(name):
        def stub(*a, **k):
            raise SolverGuardViolation(f"stationary.{name} called during manifest generation")
        return stub

    @staticmethod
    def _check_modules():
        bad = [m for m in FORBIDDEN_MODULES if m in sys.modules]
        if bad:
            raise SolverGuardViolation(f"solver/driver modules imported during manifest generation: {bad}")

    def __exit__(self, *exc):
        try:
            self._check_modules()
        finally:
            for name, obj in self._saved.items():
                setattr(self._st, name, obj)
        return False


# ------------------------------------------------------------------------------------------------ §2.6 index order
def configurations(set_name):
    """[(idx, N_y, class, N, d)] in §2.6 order: N_y → class → N → draw d."""
    out, idx = [], 0
    for ny in NYS:
        for cls in CLASSES:
            for N in NS:
                for d in range(DRAWS[set_name][cls]):
                    out.append((idx, ny, cls, N, d))
                    idx += 1
    return out


# --------------------------------------------------------------------------------------------------- §2.3 draws
def draw(root, idx, attempt, ny, N):
    gen = np.random.Generator(np.random.PCG64(np.random.SeedSequence(entropy=root, spawn_key=(idx, attempt))))
    n_open = int(gen.integers(1, ny + 1))
    u = float(gen.random())
    Jr, Ji = [], []
    for _ in range(N):
        Jr.append(gen.integers(-8, 9, size=(ny, ny)).astype(int).tolist())
        Ji.append(gen.integers(-8, 9, size=(ny, ny)).astype(int).tolist())
    return n_open, u, Jr, Ji


def v1_blocks(ny, cls, Jr, Ji):
    """V_j(1) = v0 (B_j + B_j^H)/2, B_j = (Jr_j + i Ji_j [class H])/8; class D: off-diagonal zeroed. Exact pairs."""
    out = []
    for jr, ji in zip(Jr, Ji):
        blk = []
        for a in range(ny):
            row = []
            for b in range(ny):
                re = Fraction(V0, 8) * Fraction(jr[a][b] + jr[b][a], 2)
                im = Fraction(V0, 8) * Fraction(ji[a][b] - ji[b][a], 2) if cls == "H" else Fraction(0)
                if cls == "D" and a != b:
                    re = im = Fraction(0)
                row.append((re, im))
            blk.append(tuple(row))
        out.append(tuple(blk))
    return out


# --------------------------------------------------------------------------------------------- §2.4 item 1 energy
def energy(ny, n_open, u):
    """Returns (E, E_lo, E_upper, deltas, ok, reason). Exact rationals from the double E_lo, E_upper and u."""
    Ep = [Fraction(e) for e in eperp_doubles(ny)]
    E_lo = Ep[n_open - 1]
    E_up = Ep[n_open] if n_open < ny else Ep[0] + 64 * (2 - DELTA_MIN)                # choice C-2
    E = Fraction(round(64 * (E_lo + Fraction(u) * (E_up - E_lo))), 64)             # round(Fraction): half-to-even
    deltas = [(E - e) * DX2 for e in Ep]
    bad = [n for n, d in enumerate(deltas) if not (abs(d) >= DELTA_MIN and d <= 2 - DELTA_MIN)]
    reason = None if not bad else f"energy: channel(s) {bad} violate |delta| >= 1/50 and delta <= 2 - 1/50"
    if not bad and sum(1 for d in deltas if d > 0) != n_open:
        reason = "energy: open-channel count differs from n_open"
    return E, E_lo, E_up, deltas, reason is None, reason


# ------------------------------------------------------------------------------------ §2.4 item 2 grid conditioning
def _to_np(A):
    return np.array([[complex(float(z[0]), float(z[1])) for z in row] for row in A], dtype=complex)


def _cond(X):
    with np.errstate(all="ignore"):
        try:
            c = float(np.linalg.cond(X, 2))
        except np.linalg.LinAlgError:
            return float("inf")
    return c if np.isfinite(c) else float("inf")


def _inv(X):
    return np.linalg.inv(X)


def _lattice_matrices(lat, E, Ep_q, V1, s, Dq=None, dq=None, N=None):
    """Exact full and prefix matrices of L or R_W at s, rounded once to complex128 (choice C-1)."""
    mats = {}
    for key in ["full"] + list(range(1, N)):
        kind, j = ("full", None) if key == "full" else ("prefix", key)
        if lat == "L":
            A = C.assemble_L(E, Ep_q, V1, s, kind, j)
        else:
            A = C.assemble_R_W(Dq, dq, kind, j)
        mats[key] = _to_np(A)
    return mats


def _rw_inputs(E, ny, cls, V1, s, N):
    """Production W/Z inputs in binary64 (the formation used by the driver; §4.1), imported exactly."""
    Epd = np.asarray(eperp_doubles(ny), dtype=float)
    V1d = np.array([[[complex(float(z[0]), float(z[1])) for z in r] for r in B] for B in V1])
    if cls in "RD":
        V1d = V1d.real.copy()
    Vs = float(s) * V1d
    Delta = (float(E) * np.eye(ny) - np.diag(Epd) - Vs) * (DX * DX)
    dlead = (float(E) - Epd) * (DX * DX)
    Dq = [C.import_matrix(Delta[j]) for j in range(N)]
    dq = [C.import_binary64(float(x)) for x in dlead]
    return Dq, dq


def conditioning_at(lat, E, ny, cls, V1, s, N):
    """§2.4 item 2 quantities at one grid point on one lattice: dict name -> cond_2 (binary64)."""
    Ep_q = [C.import_binary64(e) for e in eperp_doubles(ny)]
    Dq = dq = None
    if lat == "R_W":
        Dq, dq = _rw_inputs(E, ny, cls, V1, s, N)
        gam = [1 - d for d in dq]
    else:
        gam = [1 - (E - e) * DX2 for e in Ep_q]
    mats = _lattice_matrices(lat, E, Ep_q, V1, s, Dq, dq, N)
    t = float(C.T_HOP)
    out = {"M": _cond(mats["full"])}
    Gp = {}
    for j in range(1, N):
        out[f"M^({j})"] = _cond(mats[j])
        Gp[j] = _inv(mats[j] / t)
    for j in range(1, N):
        if j == 1:
            SC = mats[1] / t
        else:
            Ajj = mats["full"][(j - 1) * ny:j * ny, (j - 1) * ny:j * ny] / t
            SC = Ajj - Gp[j - 1][(j - 2) * ny:(j - 1) * ny, (j - 2) * ny:(j - 1) * ny]
        out[f"SC_{j}"] = _cond(SC)
    G = _inv(mats["full"] / t)
    GNN = G[(N - 1) * ny:, (N - 1) * ny:]
    Z = _inv(GNN)
    out["Z"] = _cond(Z)
    lam = [C.lead_enclosure(g)["lam_mid"] for g in gam]
    Lm1 = np.diag([complex(float(l[0]) - 1.0, float(l[1])) for l in lam])
    out["Z-(Lambda-I)"] = _cond(Z - Lm1)
    return out


def grid_conditioning(E, ny, cls, V1, N):
    """Returns (ok, reason, max_cond). Stops at the first violation."""
    worst = 0.0
    for s in GRID:
        for lat in ("L", "R_W"):
            try:
                c = conditioning_at(lat, E, ny, cls, V1, s, N)
            except np.linalg.LinAlgError as exc:
                return False, f"conditioning: singular matrix on {lat} at s={s} ({exc})", None
            for name, v in c.items():
                worst = max(worst, v)
                if not v <= COND_MAX:
                    return False, f"conditioning: cond2({name}) = {v:.6g} > 1e3 on {lat} at s={s}", None
    return True, None, worst


# ------------------------------------------------------------------------------------ §2.4 item 3 certificate
def certificate(E, ny, V1, N):
    Ep_q = [C.import_binary64(e) for e in eperp_doubles(ny)]
    mats = []
    for key in ["full"] + list(range(1, N)):
        kind, j = ("full", None) if key == "full" else ("prefix", key)
        name = "M_L" if key == "full" else f"M_L^({key})"
        mats.append((name, C.assemble_L(E, Ep_q, V1, Fraction(0), kind, j), C.assemble_A1(V1, kind, j),
                     C.eps_A_L(E, Ep_q, N, kind)))
    return COV.certify_configuration(mats)


# ------------------------------------------------------------------------------------------- serialisation helpers
def fs(x):
    if x is COV.INF:
        return "+inf"
    if x is None:
        return None
    x = Fraction(x)
    return f"{x.numerator}/{x.denominator}"


def _cert_record(cert):
    recs = {}
    for name, r in cert["records"].items():
        recs[name] = dict(
            covered=r["covered"], b=fs(r["b"]), b_zero=r["b_zero"], n_points=r["n_points"],
            budget_exhausted=r["budget_exhausted"], cond_bound=fs(r["cond_bound"]), min_r=fs(r["min_r"]),
            points=[dict(s=fs(p["s"]), depth=p["depth"], dps=p["dps"], r=fs(p["r"])) for p in r["certified"]],   # a = 1/(2 r b)
            uncertified=[dict(s=fs(p["s"]), depth=p["depth"]) for p in r["uncertified"]],
            r0=fs(COV.taylor_radius(r)))
    return dict(accepted=cert["accepted"], r0=fs(cert["r0"]), beta=fs(cert["beta"]), matrices=recs)


def _failed_cert_reason(cert):
    parts = []
    for name, r in cert["records"].items():
        if not r["covered"]:
            parts.append(f"{name} not covered (points={r['n_points']}, budget_exhausted={r['budget_exhausted']})")
        z = r["zero"]
        if z is None or z["status"] != "certified":
            parts.append(f"{name}: s = 0 uncertified (R3.20)")
    return "certificate: " + "; ".join(parts)


# ---------------------------------------------------------------------------------------------- one configuration
def generate_configuration(set_name, idx, ny, cls, N, d, max_attempts=MAX_ATTEMPTS):
    root = ROOTS[set_name]
    log = []
    for attempt in range(max_attempts):
        n_open, u, Jr, Ji = draw(root, idx, attempt, ny, N)
        V1 = v1_blocks(ny, cls, Jr, Ji)
        E, E_lo, E_up, deltas, ok, reason = energy(ny, n_open, u)
        if ok:
            ok, reason, worst = grid_conditioning(E, ny, cls, V1, N)
        if ok:
            cert = certificate(E, ny, V1, N)
            ok = cert["accepted"]
            if not ok:
                reason = _failed_cert_reason(cert)
        if not ok:
            log.append(dict(attempt=attempt, reason=reason))
            continue
        return dict(
            set=set_name, idx=idx, draw=d, attempt=attempt, N_y=ny, cls=cls, N=N, n_open=n_open, u=float.hex(u),
            E=fs(E), E_lo=fs(E_lo), E_upper=fs(E_up), deltas=[fs(x) for x in deltas],
            E_perp=[float.hex(e) for e in eperp_doubles(ny)],
            Jr=Jr, Ji=Ji, V1=[[[[fs(z[0]), fs(z[1])] for z in row] for row in B] for B in V1],
            max_grid_cond_binary64=float.hex(worst), rejection_log=log, certificate=_cert_record(cert))
    raise ManifestAbort(f"{set_name} idx {idx}: redraw cap of {max_attempts} attempts reached; manifest generation "
                        "aborts (DESIGN INSUFFICIENT, pre-run)")


def _worker(args):
    with solver_guard():
        return generate_configuration(*args)


# ---------------------------------------------------------------------------------------------------- versions
def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def versions():
    import mpmath
    files = ["coupling_scaling_v2_manifest.py", "coupling_scaling_v2_certify.py", "coupling_scaling_v2_covering.py",
             "PREDICTIONS_coupling_scaling_v2.md"]
    return dict(python=platform.python_version(), numpy=np.__version__, mpmath=mpmath.__version__,
                platform=platform.platform(), machine=platform.machine(),
                OPENBLAS_NUM_THREADS=os.environ.get("OPENBLAS_NUM_THREADS"),
                file_sha256={f: _sha256(os.path.join(HERE, f)) for f in files},
                stationary_sha256=_sha256(os.path.join(ROOT_DIR, "stationary.py")))


def dumps(obj):
    return json.dumps(obj, sort_keys=True, indent=1, separators=(",", ": ")) + "\n"


def build_manifest(set_name, workers=1, indices=None, max_attempts=MAX_ATTEMPTS):
    cfgs = configurations(set_name)
    if indices is not None:
        cfgs = [c for c in cfgs if c[0] in set(indices)]
    jobs = [(set_name, idx, ny, cls, N, d, max_attempts) for idx, ny, cls, N, d in cfgs]
    if workers > 1:
        import multiprocessing as mpc
        with mpc.get_context("spawn").Pool(workers) as pool:
            recs = pool.map(_worker, jobs, chunksize=1)     # map preserves input (manifest) order
    else:
        recs = [_worker(j) for j in jobs]
    assert [r["idx"] for r in recs] == [j[1] for j in jobs]
    return dict(schema=SCHEMA, set=set_name, root_entropy=ROOTS[set_name], spec="PREDICTIONS_coupling_scaling_v2.md "
                "Revision 3.7", index_order="N_y (2,3,4) -> class (R,H,D) -> N (2,3,4,6) -> draw d",
                n_configurations=len(recs), configurations=recs, versions=versions())


def write_previous_files(out_dir=OUT_DIR):
    """§8.3: sha256 of every existing file under experiments/wz_rgf_diagnostic/ (excluding __pycache__ and the output
    directory itself), sorted by path."""
    lines = []
    for dp, dns, fns in os.walk(HERE):
        dns[:] = sorted(d for d in dns if d != "__pycache__" and os.path.join(dp, d) != out_dir)
        for fn in sorted(fns):
            p = os.path.join(dp, fn)
            lines.append(f"{_sha256(p)}  {os.path.relpath(p, HERE)}")
    lines.sort(key=lambda s: s.split("  ", 1)[1])
    return "\n".join(lines) + "\n"


def main(argv=None):
    """Default: generate manifest_XP.json, manifest_HO.json and manifest.sha256 (filter and certificate only).
    --previous-files: write previous_files.sha256 (§8.3). That is a separate step, run immediately before the first solver
    run (after R0). It is not part of the manifest freeze, so later review documents do not invalidate it."""
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--out", default=OUT_DIR)
    ap.add_argument("--previous-files", action="store_true")
    a = ap.parse_args(argv)
    if a.previous_files:
        p = os.path.join(a.out, "previous_files.sha256")
        if os.path.exists(p):
            raise SystemExit("previous_files.sha256 already exists; refusing to overwrite")
        txt = write_previous_files(a.out)
        os.makedirs(a.out, exist_ok=True)
        with open(p, "w") as f:
            f.write(txt)
        print(hashlib.sha256(txt.encode()).hexdigest(), "previous_files.sha256")
        return 0
    if os.path.exists(os.path.join(a.out, "manifest.sha256")):
        raise SystemExit("manifest.sha256 already exists: the manifest is frozen; refusing to regenerate")
    for m in FORBIDDEN_MODULES:
        if m in sys.modules:
            raise SolverGuardViolation(m)
    texts = {}
    for set_name in ("XP", "HO"):
        texts[set_name] = dumps(build_manifest(set_name, a.workers))       # ManifestAbort propagates: nothing written
    os.makedirs(a.out, exist_ok=True)
    names = {}
    for set_name, txt in texts.items():
        fn = f"manifest_{set_name}.json"
        with open(os.path.join(a.out, fn), "w") as f:
            f.write(txt)
        names[fn] = hashlib.sha256(txt.encode()).hexdigest()
    with open(os.path.join(a.out, "manifest.sha256"), "w") as f:
        f.write("".join(f"{h}  {fn}\n" for fn, h in sorted(names.items())))
    for fn, h in sorted(names.items()):
        print(h, fn)
    return 0


if __name__ == "__main__":
    sys.exit(main())
