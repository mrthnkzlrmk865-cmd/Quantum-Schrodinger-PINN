"""
High-precision phase audit of the B1 clean-waveguide RGF results, for all
six B1 configurations (docs/RGF_B1_VALIDATION.md §6).

What is compared
----------------
For a clean guide with exact leads, the diagonal transmission amplitudes
under convention P are exactly t_nn = lambda_n^(N-1), with
lambda_n = exp(i theta_n), theta_n = k_n dx = arccos(c_n),
c_n = 1 - (E - E_perp_n) dx^2  (spec §2.7.2).

* RGF value: the diagonal of t_RL (incidence from the left) and of t_LR
  (incidence from the right) from stationary.scattering_matrix, double
  precision. This is unchanged production code.
* High-precision reference: the EXACT value of lambda_n^(N-1) for the lattice
  problem with the parameters the RGF was given. E (the spec-grid double),
  E_perp_n (the same double array used by the RGF), dx and N are taken
  exactly as given (doubles are exact binary rationals), and then
      delta = (E - E_perp_n) dx^2,   theta = arccos(1 - delta),
      phi_ref = (N - 1) theta
  is evaluated in mpmath at `dps` significant digits. No double-precision
  intermediate of the RGF or of smatrix_spec.lead_modes (lambda, k, angles)
  enters the reference. The reference is NOT the continuum phase.
* Cross-checks of the reference itself: the same theta from the
  cancellation-free form theta = 2 arcsin(sqrt(delta/2)), and the full
  computation at dps in {30, 50, 100}. The reported reference is the
  dps = 100 arccos form.

Phase-difference convention
---------------------------
dphi = arg( t_RGF * exp(-i phi_ref) ) in (-pi, pi]: the principal value of
phi_RGF - phi_ref. It is independent of the 2 pi branch of each phase, and
valid while the true difference is below pi in magnitude (here ~1e-10).
Reported phases phi_rgf = arg(t_RGF) and phi_ref (reduced mod 2 pi into
(-pi, pi]) are both principal values. The unwrapped (N-1) theta, which can
exceed 10^3 rad, is also stored. Magnitude errors ||t| - 1| are reported
separately.

B1 criterion
------------
The predeclared B1 phase criterion (docs/PHASE2_METHOD_SELECTION.md §3,
threshold 1e-10) is defined on |t_nn - lambda_double^(N-1)|, i.e. against
the double-precision reference. That quantity (b1_phase_quantity) is
recomputed here, and the predeclared pass/fail is decided on it. The same
1e-10 value applied to |dphi| against the high-precision reference is
reported as a DIAGNOSTIC, not as a redefinition of B1.

Outputs (results/rgf_b1_validation/):
  phase_audit_hp.csv                  one row per (config, energy, open channel, direction)
  phase_audit_hp_summary.json         per-configuration maxima, RMS, argmax, pass/fail
  phase_audit_hp_metadata.json        commit, dirty flag, versions, precision, command, time
Existing B1 artifacts (b1_grid.csv, phase_exact.csv, ...) are not modified.

Usage:
  python audit_rgf_phase.py                      # all six configurations (~12 min, 1 BLAS thread)
  python audit_rgf_phase.py --configs A F --out DIR
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")    # performance; before NumPy (stationary.py docstring)

import argparse
import csv
import datetime as _dt
import json
import sys
import time

import mpmath as mp
import numpy as np

from provenance import code_fingerprint, environment_info
from stationary import scattering_matrix
from validate_rgf_b1 import LIMITS, device, spec_energies

CONFIGS = {   # id: (N_y, dx, L_d) -- the six B1 configurations of validate_rgf_b1.py
    "A": (63, 1 / 512, 1.0),
    "B": (127, 1 / 512, 1.0),
    "C": (63, 1 / 1024, 1.0),
    "D": (127, 1 / 1024, 1.0),
    "E": (63, 1 / 512, 3.0),
    "F": (63, 1 / 1024, 3.0),
}
REF_DPS = 100
CHECK_DPS = (30, 50)
PHASE_THRESHOLD = LIMITS["phase"]          # 1e-10, unchanged


def wrap_to_pi(x, dps=None):
    """Reduce a real (float or mpf) into (-pi, pi].
    For mpf input the reduction is done at `dps` digits (required: mpmath
    arithmetic otherwise uses the ambient precision, which would make the
    result depend on hidden global state)."""
    if isinstance(x, mp.mpf):
        if dps is None:
            raise ValueError("wrap_to_pi: pass dps for mpmath input")
        with mp.workdps(dps):
            return _wrap_mp(x)
    y = (x + np.pi) % (2 * np.pi) - np.pi
    return np.where(y == -np.pi, np.pi, y) if np.ndim(y) else (np.pi if y == -np.pi else y)


def _wrap_mp(x):
    """(-pi, pi] reduction of an mpf at the current mpmath precision."""
    two_pi = 2 * mp.pi
    y = x - two_pi * mp.floor((x + mp.pi) / two_pi)              # [-pi, pi)
    return y if y != -mp.pi else +mp.pi


def reference_phase_wrapped(E, Eperp, dx, N, dps):
    """Exact lattice phase (N-1) theta reduced into (-pi, pi], computed and
    reduced entirely at `dps` digits, returned as float."""
    return float(wrap_to_pi(reference_phase(E, Eperp, dx, N, dps), dps=dps))


def wrapped_difference(z, phi_ref):
    """arg(z * exp(-i phi_ref)) in (-pi, pi]: principal value of arg(z) - phi_ref,
    independent of 2 pi branches. z complex (nonzero), phi_ref real."""
    return float(np.angle(z * np.exp(-1j * float(phi_ref))))


def reference_phase(E, Eperp, dx, N, dps, form="acos"):
    """Unwrapped exact lattice phase (N-1) theta with theta = k dx, at `dps` digits."""
    with mp.workdps(dps):
        delta = (mp.mpf(float(E)) - mp.mpf(float(Eperp))) * mp.mpf(float(dx)) ** 2
        if not (0 < delta < 2):
            raise ValueError("channel not open in the lattice dispersion")
        if form == "acos":
            theta = mp.acos(1 - delta)
        elif form == "asin":
            theta = 2 * mp.asin(mp.sqrt(delta / 2))               # 1 - cos(theta) = 2 sin^2(theta/2)
        else:
            raise ValueError(form)
        return (N - 1) * theta


def audit_config(cid):
    ny, dx, Ld = CONFIGS[cid]
    dev = device(ny, dx, Ld)
    N, Ep = dev.n_slices, dev.basis.E_perp
    rows = []
    for E in spec_energies(Ep):
        r = scattering_matrix(E, dev)
        m = r["modes"]
        open_idx = np.where(m["open"])[0]
        for k, n in enumerate(open_idx):
            phi_u = {d: reference_phase(E, Ep[n], dx, N, d) for d in (REF_DPS, *CHECK_DPS)}
            phi_asin = reference_phase(E, Ep[n], dx, N, REF_DPS, form="asin")
            with mp.workdps(REF_DPS):
                phi_ref = wrap_to_pi(phi_u[REF_DPS], dps=REF_DPS)
                ch_dps = {d: float(abs(phi_u[d] - phi_u[REF_DPS])) for d in CHECK_DPS}
                ch_form = float(abs(phi_asin - phi_u[REF_DPS]))
                ref_complex = complex(mp.cos(phi_u[REF_DPS]), mp.sin(phi_u[REF_DPS]))
            lam_d = m["lam"][n]
            dbl_ref = lam_d ** (N - 1)
            for direction, t in (("L->R", r["t_RL"][k, k]), ("R->L", r["t_LR"][k, k])):
                dphi = wrapped_difference(t, phi_ref)
                rows.append(dict(
                    config=cid, N_y=ny, dx=dx, L_d=Ld, N=N, E=float(E), channel=int(n) + 1, direction=direction,
                    sin_kdx=float(np.sin(m["k"][n] * dx)),
                    phase_rgf=float(np.angle(t)), phase_ref_hp=float(phi_ref),
                    phase_ref_hp_unwrapped=float(phi_u[REF_DPS]),
                    dphase_wrapped=dphi, abs_phase_err=abs(dphi),
                    mag_err=abs(abs(t) - 1.0),
                    err_complex_vs_hp=abs(t - ref_complex),
                    dphase_double_ref_vs_hp=wrapped_difference(dbl_ref, phi_ref),
                    b1_phase_quantity=abs(t - dbl_ref),
                    ref_change_dps30_vs_dps100=ch_dps[30], ref_change_dps50_vs_dps100=ch_dps[50],
                    ref_acos_vs_asin_dps100=ch_form,
                    ref_dps=REF_DPS))
    return rows


def summarize(rows, cids):
    out = {}
    for cid in cids:
        rr = [r for r in rows if r["config"] == cid]
        err = np.array([r["abs_phase_err"] for r in rr])
        b1q = np.array([r["b1_phase_quantity"] for r in rr])
        dref = np.array([abs(r["dphase_double_ref_vs_hp"]) for r in rr])
        i, j = int(np.argmax(err)), int(np.argmax(b1q))
        out[cid] = dict(
            N_y=rr[0]["N_y"], dx=rr[0]["dx"], L_d=rr[0]["L_d"], N=rr[0]["N"],
            n_energies=len({r["E"] for r in rr}), n_rows=len(rr),
            hp_max_abs_phase_err=float(err[i]), hp_rms_phase_err=float(np.sqrt(np.mean(err**2))),
            hp_argmax=dict(E=rr[i]["E"], channel=rr[i]["channel"], direction=rr[i]["direction"],
                           sin_kdx=rr[i]["sin_kdx"]),
            hp_diagnostic_pass_at_threshold=bool(err[i] <= PHASE_THRESHOLD),
            b1_phase_quantity_max=float(b1q[j]),
            b1_argmax=dict(E=rr[j]["E"], channel=rr[j]["channel"], direction=rr[j]["direction"]),
            b1_phase_criterion_pass=bool(b1q[j] <= PHASE_THRESHOLD),
            double_ref_vs_hp_max_abs_phase=float(dref.max()),
            max_mag_err=float(max(r["mag_err"] for r in rr)),
            max_ref_change_dps30_vs_dps100=float(max(r["ref_change_dps30_vs_dps100"] for r in rr)),
            max_ref_change_dps50_vs_dps100=float(max(r["ref_change_dps50_vs_dps100"] for r in rr)),
            max_ref_acos_vs_asin=float(max(r["ref_acos_vs_asin_dps100"] for r in rr)))
    return out


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--configs", nargs="*", default=list(CONFIGS), choices=list(CONFIGS))
    ap.add_argument("--out", default="results/rgf_b1_validation")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    rows = []
    for cid in a.configs:
        t1 = time.time()
        rr = audit_config(cid)
        rows += rr
        print(f"  config {cid} {CONFIGS[cid]}: {len(rr)} rows in {time.time() - t1:.1f}s", flush=True)
    summary = dict(
        phase_threshold=PHASE_THRESHOLD,
        b1_criterion_definition="max |t_nn - lambda_double^(N-1)| over energies, open channels, both directions",
        hp_reference="exact lattice phase (N-1) arccos(1 - (E - E_perp) dx^2), mpmath",
        phase_difference_convention="arg(t_RGF exp(-i phi_ref)) in (-pi, pi]",
        ref_dps=REF_DPS, check_dps=list(CHECK_DPS),
        configs=summarize(rows, a.configs),
        wall_time_s=time.time() - t0)
    summary["all_configs_b1_phase_pass"] = all(v["b1_phase_criterion_pass"] for v in summary["configs"].values())
    write_csv(os.path.join(a.out, "phase_audit_hp.csv"), rows)
    with open(os.path.join(a.out, "phase_audit_hp_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    env = environment_info()
    meta = dict(
        generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        command=" ".join([os.path.basename(sys.executable)] + sys.argv),
        source_commit=env.get("git_commit"), working_tree_dirty=env.get("git_dirty"),
        python=env.get("python"), numpy=env.get("numpy"), scipy=env.get("scipy"),
        mpmath=mp.__version__, platform=env.get("platform"), machine=env.get("machine"),
        arbitrary_precision=dict(library="mpmath", reference_dps=REF_DPS, check_dps=list(CHECK_DPS)),
        openblas_num_threads=os.environ.get("OPENBLAS_NUM_THREADS"),
        configs={c: dict(zip(("N_y", "dx", "L_d"), CONFIGS[c])) for c in a.configs},
        wall_time_s=summary["wall_time_s"],
        **code_fingerprint("rgf"))
    with open(os.path.join(a.out, "phase_audit_hp_metadata.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps({c: {k: v[k] for k in ("N", "n_rows", "hp_max_abs_phase_err", "b1_phase_quantity_max",
                                            "b1_phase_criterion_pass", "double_ref_vs_hp_max_abs_phase",
                                            "max_ref_change_dps50_vs_dps100")}
                      for c, v in summary["configs"].items()}, indent=1))
