"""
Spreading of a wave packet launched inside a closed box filled with
correlated disorder: a FEASIBILITY test for localization diagnostics.

Observables (observables.py), per realization and time:
  MSD(t) = <|r - r0|^2>   mean-square displacement from the launch point
  PR(t)  = (int rho)^2 / int rho^2   participation ratio (effective area)

In a closed box both saturate at late times even without localization,
because the box is finite. The finite-size test compares two box sizes
L = 1 and L = 2 (same disorder statistics, same packet):
  * extended (delocalized) states: the late-time PR / L^2 is about the same
    as in the CLEAN box of the same size, for both sizes
  * localized states with localization length xi << L: the late-time PR
    is set by xi, not by L, so PR / L^2 falls when L doubles
    (PR_sat(L=2) / PR_sat(L=1) -> 1 instead of 4).
The test only works if xi is clearly smaller than the box. Section
"Localization feasibility" of the README discusses the estimates.

Usage:  python experiment_localization.py [--n-real 16]
Output: results/localization/<point>.npz and summary.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time

import numpy as np

from config import closed_box_config
from disorder import disorder_from_config
from observables import mean_square_displacement, participation_ratio, uniform_box_msd
from potential import Potential
from reference_solver import SplitStepDST
from stats_utils import summarize

OUT = "results/localization"
STRENGTHS = [0, 100, 200, 300]
SIZES = [1.0, 2.0]


def run_point(L, W, n_real, t_end, batch=8):
    cfg = closed_box_config(L, disorder_strength=float(W), T=t_end)
    seeds = list(range(1 if W == 0 else n_real))
    msd, pr = [], []
    t0 = time.time()
    for s in range(0, len(seeds), batch):
        pots = [Potential(cfg, disorder_from_config(cfg, seed=sd)) for sd in seeds[s:s + batch]]
        solver = SplitStepDST(cfg, pots)

        def obs(psi, slv):
            return dict(msd=mean_square_displacement(psi, slv.X, slv.Y, (cfg.x0, cfg.y0)),
                        pr=participation_ratio(psi, slv.dx, slv.dy))

        res = solver.run(n_frames=101, store_frames=False, observables=obs)
        msd.append(res["msd"])
        pr.append(res["pr"])
    return cfg, res["times"], np.concatenate(msd, axis=1), np.concatenate(pr, axis=1), time.time() - t0


def main(n_real):
    os.makedirs(OUT, exist_ok=True)
    rows = []
    for L in SIZES:
        t_end = 0.25 * L + 0.25   # several box traversals at v ~ 20
        for W in STRENGTHS:
            name = f"L{L}_W{W}"
            path = os.path.join(OUT, name + ".npz")
            if os.path.exists(path):
                d = np.load(path)
                cfg = closed_box_config(L, disorder_strength=float(W), T=t_end)
                times, msd, pr = d["times"], d["msd"], d["pr"]
            else:
                cfg, times, msd, pr, wt = run_point(L, W, n_real, t_end)
                np.savez_compressed(path, times=times, msd=msd, pr=pr, config_json=json.dumps(cfg.to_dict()))
                print(f"  {name}: done ({wt:.0f}s)", flush=True)
            late = times >= 0.5 * times[-1]          # time average over the second half
            pr_late = pr[late].mean(axis=0)
            msd_late = msd[late].mean(axis=0)
            row = dict(L=L, W=W, zeta=cfg.disorder_corr_length, kx0=cfg.kx0, n=pr.shape[1],
                       area=L * L, msd_uniform=uniform_box_msd(L, L, (cfg.x0, cfg.y0)))
            for key, v in (("PR_late", pr_late), ("PR_late_over_area", pr_late / L**2),
                           ("MSD_late", msd_late), ("MSD_late_over_uniform", msd_late / row["msd_uniform"])):
                s = summarize(v)
                for sk in ("mean", "std", "se", "ci_boot_lo", "ci_boot_hi"):
                    row[f"{key}_{sk}"] = s[sk]
            rows.append(row)
    with open(os.path.join(OUT, "summary.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-real", type=int, default=16)
    a = ap.parse_args()
    for r in main(a.n_real):
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if "std" not in k})
