r"""
EXPLORATORY follow-up of the B3 dx-refinement failure (not part of the
predeclared B3 protocol; docs/RGF_B3_VALIDATION.md §8.2).

The predeclared criterion |T_1(dx = 1/1024) - T_1(dx = 1/512)| <= 2e-4 failed
for case D (W = 200, zeta = 0.05) at 3 of its 5 energies. This script only
diagnoses that failure; it changes no threshold and replaces no B3 result.
For the failing case it computes, at the same five energies and with N_y = 63,
N_q = 127, buffer 0.2012 (all fixed):
  * RGF T_1 at dx = 1/256, 1/512, 1/1024, 1/2048  (nested lattices),
  * the sparse reference T_1 at dx = 1/1024  (is the dx dependence seen by
    both methods, i.e. a property of the discretization and not of the RGF?),
  * the observed orders and a Richardson estimate IF the order is in [1.8, 2.2].
Writes results/rgf_b3_validation/exploratory_dx_case_D.csv and
exploratory_dx_case_D_summary.json.
Usage: python explore_b3_dx.py [--case D]
"""

from __future__ import annotations

import os

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import argparse
import csv
import datetime as _dt
import json
import sys
import time

import numpy as np

import sparse_reference as sr
import validate_rgf_b3 as b3
from provenance import environment_info
from stationary import scattering_matrix

DXS = (1 / 256, 1 / 512, 1 / 1024, 1 / 2048)
SPARSE_DX = 1 / 1024


def order(a, b, c):
    """log2(|a - b| / |b - c|) for T(h), T(h/2), T(h/4)."""
    return float(np.log2(abs(a - b) / abs(b - c)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="D")
    ap.add_argument("--out", default="results/rgf_b3_validation")
    a = ap.parse_args()
    case = a.case
    cfg, pot = b3.case_potential(case)
    E_list = b3.case_energies(cfg)
    t0 = time.time()
    T = {i: {} for i in range(len(E_list))}
    sp_T = {}
    for dx in DXS:
        var = dict(b3.BASE, dx=dx)
        dev, geo, info, _ = b3.make_rgf_device(case, cfg, pot, var)
        assert info["bitwise_equal"]
        for i, e in enumerate(E_list):
            r = scattering_matrix(e["E"], dev)
            T[i][dx] = b3.T1(r["S"], r["n_open"])
            if dx == SPARSE_DX:
                s = sr.sparse_scattering_matrix(e["E"], pot.numpy, geo["x_start"], dx, dev.n_slices, cfg.Ly,
                                                var["Ny"], var["Nq"])
                sp_T[i] = dict(T1=b3.T1(s["S"], s["n_open"]), max_abs_dS=float(np.max(b3.s_difference(r["S"], s["S"]))))
        print(f"  dx = 1/{round(1 / dx)}: N = {dev.n_slices} done ({time.time() - t0:.0f}s)", flush=True)
    rows = []
    for i, e in enumerate(E_list):
        t = [T[i][dx] for dx in DXS]
        p1, p2 = order(t[0], t[1], t[2]), order(t[1], t[2], t[3])
        d512 = t[2] - t[1]
        rich_ok = bool(1.8 <= p1 <= 2.2 and 1.8 <= p2 <= 2.2)
        T_ext = t[3] + (t[3] - t[2]) / 3.0 if rich_ok else float("nan")
        rows.append(dict(
            case=case, kx=e["kx"], E=e["E"], n_open=e["n_open_exact"], **{f"T1_dx_1_{round(1 / dx)}": T[i][dx] for dx in DXS},
            delta_512_to_1024=d512, delta_1024_to_2048=t[3] - t[2], order_256_512_1024=p1, order_512_1024_2048=p2,
            richardson_applied=rich_ok, T1_richardson_from_1024_2048=T_ext,
            est_error_dx512_vs_extrapolated=(t[1] - T_ext) if rich_ok else float("nan"),
            est_error_dx1024_vs_extrapolated=(t[2] - T_ext) if rich_ok else float("nan"),
            sparse_T1_dx_1_1024=sp_T[i]["T1"], sparse_minus_rgf_T1_dx_1_1024=sp_T[i]["T1"] - T[i][SPARSE_DX],
            sparse_vs_rgf_max_abs_dS_dx_1_1024=sp_T[i]["max_abs_dS"], label="exploratory"))
    with open(os.path.join(a.out, f"exploratory_dx_case_{case}.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    ords = [r["order_256_512_1024"] for r in rows] + [r["order_512_1024_2048"] for r in rows]
    summary = dict(
        label="exploratory (diagnosis of the failed B3 dx criterion); not part of the predeclared protocol",
        case=case, params=b3.CASES[case], dxs=list(DXS), N_y=63, N_q=127, buffer=b3.buffer_length(0.2),
        max_abs_delta_512_to_1024=float(max(abs(r["delta_512_to_1024"]) for r in rows)),
        max_abs_delta_1024_to_2048=float(max(abs(r["delta_1024_to_2048"]) for r in rows)),
        orders=dict(min=float(min(ords)), max=float(max(ords)), median=float(np.median(ords))),
        n_richardson_applied=int(sum(r["richardson_applied"] for r in rows)),
        max_abs_est_error_dx512=float(np.nanmax([abs(r["est_error_dx512_vs_extrapolated"]) for r in rows])) if any(r["richardson_applied"] for r in rows) else None,
        max_abs_est_error_dx1024=float(np.nanmax([abs(r["est_error_dx1024_vs_extrapolated"]) for r in rows])) if any(r["richardson_applied"] for r in rows) else None,
        max_abs_sparse_minus_rgf_T1_dx1024=float(max(abs(r["sparse_minus_rgf_T1_dx_1_1024"]) for r in rows)),
        max_sparse_vs_rgf_S_dx1024=float(max(r["sparse_vs_rgf_max_abs_dS_dx_1_1024"] for r in rows)),
        dT1_limit_of_the_failed_criterion=b3.LIMITS["dT1"], wall_time_s=time.time() - t0)
    env = environment_info()
    summary.update(generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
                   command=" ".join([os.path.basename(sys.executable)] + sys.argv),
                   source_commit=env.get("git_commit"), working_tree_dirty=env.get("git_dirty"),
                   python=env.get("python"), numpy=np.__version__)
    with open(os.path.join(a.out, f"exploratory_dx_case_{case}_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    print(json.dumps(summary, indent=1))
