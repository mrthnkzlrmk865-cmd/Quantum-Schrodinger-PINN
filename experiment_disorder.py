"""
Transport of a quantum wave packet through a correlated disordered slab in
a 2D hard-wall waveguide (open boundaries via CAPs). Reference solver only.

Research question
-----------------
How do the disorder strength W and the spatial correlation length zeta of
a static, Gaussian-correlated random potential control the transmission
and reflection of a wave packet, and in which parameter range is the
transport ballistic, diffusive, or strongly suppressed? Are the observed
trends consistent with the Born-approximation transport mean free path
l_tr(W, zeta, k) (disorder.born_rates), and where does that perturbative
picture fail?

Sweeps (python experiment_disorder.py --sweep NAME [--n-real N]):
  strength     W in {0,25,50,75,100,150,200,300}, zeta = 0.05, Ld = 1
  corr_length  zeta in {0.02,...,0.2}, W = 100, Ld = 1
  momentum     kx0 in {12,16,20,25,30}, W = 100, zeta = 0.05, Ld = 1
  length       Ld in {0.25,0.5,1,2,3}, W in {100, 200}, zeta = 0.05
  barrier      barrier (width 0.1, centre of slab) of height 150 or 250
               with W in {0, 50, 100, 200}
  convergence  W = 100, zeta = 0.05, many realizations (running mean / SE)
  all          everything above

Each configuration point uses seeds 0..N-1 (independent realizations; the
same seeds are reused across points, which makes neighbouring points
positively correlated and differences between them less noisy). Outputs:
  results/transport/<sweep>/<point>.npz   raw T(t), R(t), P_in(t) per seed + config json
  results/transport/<sweep>/summary.csv   statistics per point
Caching: an existing .npz is reused only if its configuration fingerprint
(provenance.py) matches the requested configuration exactly (config,
t_end, n_frames, seeds, number of realizations). Otherwise the script stops
with CacheMismatchError instead of silently using or overwriting it; pass
--recompute to overwrite, or --out DIR to write to a fresh directory.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time

import numpy as np

from config import waveguide_config
from disorder import born_rates
from provenance import (check_cache, make_metadata, read_npz_metadata, solver_config_dict,
                        write_dir_metadata)
from stats_utils import summarize
from transport import run_transport

OUT = "results/transport"
STRENGTHS = [0, 25, 50, 75, 100, 150, 200, 300]
ZETAS = [0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2]
MOMENTA = [12, 16, 20, 25, 30]
LENGTHS = [0.25, 0.5, 1.0, 2.0, 3.0]


def sweep_points(name, n_real):
    pts = []
    if name == "strength":
        for W in STRENGTHS:
            pts.append((f"W{W}", dict(disorder_strength=float(W)), 1 if W == 0 else n_real, 1.0))
    elif name == "corr_length":
        for z in ZETAS:
            pts.append((f"zeta{z}", dict(disorder_strength=100.0, disorder_corr_length=z), n_real, 1.0))
    elif name == "momentum":
        for k in MOMENTA:
            t_end = max(1.0, 1.0 * 20.0 / k)
            pts.append((f"k{k}", dict(disorder_strength=100.0, kx0=float(k)), n_real, t_end))
            pts.append((f"k{k}_clean", dict(disorder_strength=0.0, kx0=float(k)), 1, t_end))
    elif name == "length":
        for W in (100, 200):
            for Ld in LENGTHS:
                pts.append((f"W{W}_L{Ld}", dict(disorder_strength=float(W), slab_length=Ld), n_real, 1.0 + 0.5 * Ld))
    elif name == "barrier":
        for Vb in (150.0, 250.0):
            for W in (0, 50, 100, 200):
                pts.append((f"Vb{int(Vb)}_W{W}", dict(disorder_strength=float(W), barrier_height=Vb),
                            1 if W == 0 else n_real, 1.5))
    elif name == "convergence":
        pts.append(("W100_zeta0.05", dict(disorder_strength=100.0), n_real, 1.0))
    else:
        raise ValueError(name)
    return pts


N_FRAMES = 41


def transport_spec(cfg, n_real, t_end, n_frames=N_FRAMES):
    """Everything that determines the numbers in one sweep point (see
    provenance.py). Seeds are listed explicitly: realization i uses
    disorder seed i."""
    return dict(experiment="transport", config=solver_config_dict(cfg),
                run=dict(t_end=t_end, n_frames=n_frames, n_real=n_real, seeds=list(range(n_real))))


def run_point(sweep, name, overrides, n_real, t_end, out=None, recompute=False, strict_code=False,
              strict_env=False):
    d = os.path.join(out or OUT, sweep)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, f"{name}.npz")
    cfg = waveguide_config(**dict(overrides))
    spec = transport_spec(cfg, n_real, t_end)
    exists = os.path.exists(path)
    status = check_cache(read_npz_metadata(path) if exists else None, exists, spec, "transport",
                         path, recompute=recompute, strict_code=strict_code, strict_env=strict_env)
    if status == "hit":
        return cfg, dict(np.load(path, allow_pickle=False))
    r = run_transport(cfg, range(n_real), t_end=t_end, n_frames=N_FRAMES)
    meta = make_metadata("transport", spec, sweep=sweep, point=name, wall_time_s=r["wall_time_s"])
    np.savez_compressed(path, times=r["times"], T=r["T"], R=r["R"], P_in=r["P_in"], seeds=r["seeds"],
                        wall_time_s=r["wall_time_s"], config_json=json.dumps(cfg.to_dict()),
                        metadata_json=json.dumps(meta), fingerprint=meta["fingerprint"])
    print(f"  {sweep}/{name}: n={n_real} <T>={r['T'][-1].mean():.4f} <R>={r['R'][-1].mean():.4f} "
          f"<P_in>={r['P_in'][-1].mean():.2e}  ({r['wall_time_s']:.0f}s)", flush=True)
    return cfg, dict(np.load(path, allow_pickle=False))


def summarize_point(sweep, name, cfg, data):
    T, R, P = data["T"][-1], data["R"][-1], data["P_in"][-1]
    k = float(np.hypot(cfg.kx0, cfg.k_transverse))
    b = born_rates(cfg.disorder_strength, cfg.disorder_corr_length, k) if cfg.disorder_strength > 0 else None
    row = dict(sweep=sweep, point=name, W=cfg.disorder_strength, zeta=cfg.disorder_corr_length,
               kx0=cfg.kx0, k=k, E0=cfg.mean_energy, Ld=cfg.disorder_region[1] - cfg.disorder_region[0],
               barrier=cfg.barrier_height, t_end=float(data["times"][-1]),
               l_s_born=float(b["l_s"]) if b else np.inf, l_tr_born=float(b["l_tr"]) if b else np.inf,
               k_zeta=k * cfg.disorder_corr_length,
               fingerprint=str(data["fingerprint"]) if "fingerprint" in data else "")
    for key, v in (("T", T), ("R", R), ("P_in", P), ("lnT", np.log(np.maximum(T, 1e-300)))):
        s = summarize(v)
        for sk in ("n", "mean", "std", "se", "ci_t_lo", "ci_t_hi", "ci_boot_lo", "ci_boot_hi", "median"):
            row[f"{key}_{sk}"] = s[sk]
    return row


def run_sweep(sweep, n_real, out=None, recompute=False, strict_code=False, strict_env=False):
    out = out or OUT
    rows, points = [], []
    t0 = time.time()
    for name, over, n, t_end in sweep_points(sweep, n_real):
        cfg, data = run_point(sweep, name, over, n, t_end, out, recompute, strict_code, strict_env)
        rows.append(summarize_point(sweep, name, cfg, data))
        points.append(dict(point=name, fingerprint=rows[-1]["fingerprint"], n_real=n, t_end=t_end))
    with open(os.path.join(out, sweep, "summary.csv"), "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)
    write_dir_metadata(os.path.join(out, sweep), make_metadata(
        "transport", dict(experiment="transport_sweep", sweep=sweep, n_real_default=n_real, points=points),
        wall_time_s=time.time() - t0))
    return rows


DEFAULT_N = dict(strength=48, corr_length=48, momentum=32, length=24, barrier=24, convergence=192)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", default="all")
    ap.add_argument("--n-real", type=int, default=None)
    ap.add_argument("--out", default=OUT, help="output root directory (default results/transport)")
    ap.add_argument("--recompute", action="store_true",
                    help="overwrite cached results whose configuration does not match")
    ap.add_argument("--strict-code", action="store_true",
                    help="also refuse cached results produced by different source code")
    ap.add_argument("--strict-env", action="store_true",
                    help="also refuse cached results produced with different numerical library versions")
    ap.add_argument("--strict", action="store_true", help="--strict-code and --strict-env")
    a = ap.parse_args()
    names = list(DEFAULT_N) if a.sweep == "all" else [a.sweep]
    for nm in names:
        print(f"== sweep {nm}", flush=True)
        run_sweep(nm, a.n_real or DEFAULT_N[nm], a.out, a.recompute, a.strict_code or a.strict,
                  a.strict_env or a.strict)
