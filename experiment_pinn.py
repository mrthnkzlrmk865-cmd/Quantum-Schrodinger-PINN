"""
PINN performance study on the PINN problem (config.pinn_config): unit box,
Gaussian packet kx0 = 15 (E0 = 112.5), T = 0.02, Dirichlet walls.

Every run is compared with the DST reference solution of the identical
configuration (same disorder realization) by evaluate.evaluate_run.

Runs (name -> overrides). "s<k>" = network initialization / sampling seed,
"d<k>" = disorder realization seed.
  clean_*                  V = 0 (an analytic reference also exists)
  W28/W56/W112             disorder RMS W = E0/4, E0/2, E0 at zeta = 0.08
  W56_zeta0.04 / 0.16      correlation length
  W56_d1, W56_d2           other disorder realizations
  *_nofourier              plain MLP on normalized inputs
  *_rad                    residual-adaptive sampling
  clean_T0.01 / T0.04      half / twice the time window
  clean_long               3x training budget
All runs: 2000 Adam steps (fresh 4000-point batch per step) + 200 L-BFGS
steps on 16000 fixed points unless stated. Diagnostic runs that led to this
protocol are in results/pinn/diag_* (see README).

Usage:  python experiment_pinn.py [--only NAME ...] [--quick] [--out DIR] [--recompute]
Output: <out>/<run>/ (model.pt, config.json, run_metadata.json, ...) and
<out>/summary.csv. A trained model is reused only if run_metadata.json has
the same configuration fingerprint (provenance.py); otherwise the script
stops instead of silently reusing or overwriting it. --quick runs go to
<run>_quick directories.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import time

import torch

from config import pinn_config
from evaluate import evaluate_run
from potential import Potential
from reference_solver import SplitStepDST
from provenance import check_cache, make_metadata, pinn_config_dict, read_dir_metadata, write_dir_metadata
from train import save_run, train

OUT = "results/pinn"
E0 = 112.5

RUNS = {
    # network-seed variability, clean potential
    "clean_fourier_uniform_s0": dict(),
    "clean_fourier_uniform_s1": dict(seed=1),
    "clean_fourier_uniform_s2": dict(seed=2),
    # disorder strength (zeta = 0.08)
    "W28_fourier_uniform_s0": dict(disorder_strength=E0 / 4),
    "W56_fourier_uniform_s0": dict(disorder_strength=E0 / 2),
    "W112_fourier_uniform_s0": dict(disorder_strength=E0),
    "W56_fourier_uniform_s1": dict(disorder_strength=E0 / 2, seed=1),
    "W56_fourier_uniform_s2": dict(disorder_strength=E0 / 2, seed=2),
    # correlation length and realization
    "W56_zeta0.04_fourier_uniform_s0": dict(disorder_strength=E0 / 2, disorder_corr_length=0.04),
    "W56_zeta0.16_fourier_uniform_s0": dict(disorder_strength=E0 / 2, disorder_corr_length=0.16),
    "W56_d1_fourier_uniform_s0": dict(disorder_strength=E0 / 2, disorder_seed=1),
    # architecture / sampling
    "clean_nofourier_uniform_s0": dict(use_fourier=False),
    "W56_nofourier_uniform_s0": dict(disorder_strength=E0 / 2, use_fourier=False),
    "clean_fourier_rad_s0": dict(sampling="rad"),
    "W56_fourier_rad_s0": dict(disorder_strength=E0 / 2, sampling="rad"),
    # problem difficulty: time window
    "clean_T0.01_fourier_uniform_s0": dict(T=0.01),
    "clean_T0.04_fourier_uniform_s0": dict(T=0.04),
    # training budget (3x Adam, 2.5x L-BFGS)
    "clean_long_fourier_uniform_s0": dict(epochs_adam=6000, epochs_lbfgs=500, lr_decay_every=2000),
}


def pinn_spec(cfg):
    """Everything that determines a PINN run's numbers: the full PINN
    configuration (problem, potential, disorder seed, architecture,
    training schedule, network seed). The reference solution used for
    evaluation is defined by the same config."""
    return dict(experiment="pinn", config=pinn_config_dict(cfg))


def run_one(name, overrides, quick=False, out=None, recompute=False, strict_code=False):
    over = dict(overrides)
    if quick:
        over.update(epochs_adam=300, epochs_lbfgs=20)
        name = name + "_quick"          # never share a directory with a full run
    d = os.path.join(out or OUT, name)
    cfg = pinn_config(**over)
    spec = pinn_spec(cfg)
    exists = os.path.exists(os.path.join(d, "model.pt"))
    status = check_cache(read_dir_metadata(d), exists, spec, "pinn", d,
                         recompute=recompute, strict_code=strict_code)
    if status != "hit":
        print(f"== training {name}", flush=True)
        model, hist, info = train(cfg, verbose=False)
        # cost of the reference solution for the same problem
        t0 = time.time()
        SplitStepDST(cfg, Potential(cfg)).run(store_frames=False)
        info["reference_time_s"] = time.time() - t0
        save_run(d, cfg, model, hist, info)
        write_dir_metadata(d, make_metadata("pinn", spec, run=name, torch_threads=torch.get_num_threads(),
                                            device=cfg.device, train_time_s=info["train_time_s"]))
    summ = evaluate_run(d)
    with open(os.path.join(d, "train_info.json")) as f:
        info = json.load(f)
    row = dict(run=name, fingerprint=read_dir_metadata(d)["fingerprint"], W=cfg.disorder_strength, zeta=cfg.disorder_corr_length, disorder_seed=cfg.disorder_seed,
               seed=cfg.seed, T=cfg.T, fourier=cfg.use_fourier, sampling=cfg.sampling, **summ,
               final_loss=info.get("final_loss"), train_time_s=info["train_time_s"],
               reference_time_s=info.get("reference_time_s"))
    print(f"   {name}: relL2 final={summ['rel_l2_final']:.3e} mean={summ['rel_l2_mean']:.3e} "
          f"rho={summ['rel_l2_rho_final']:.3e} norm_dev={summ['norm_max_dev']:.2e} "
          f"train={info['train_time_s']:.0f}s", flush=True)
    return row


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--recompute", action="store_true")
    ap.add_argument("--strict-code", action="store_true")
    a = ap.parse_args()
    if a.threads:
        torch.set_num_threads(a.threads)
    names = a.only or list(RUNS)
    rows = [run_one(n, RUNS[n], a.quick, a.out, a.recompute, a.strict_code) for n in names]
    # merge with existing summary (rows carry their configuration fingerprint)
    path = os.path.join(a.out, "summary.csv")
    old = {}
    if os.path.exists(path):
        with open(path) as f:
            old = {r["run"]: r for r in csv.DictReader(f)}
    for r in rows:
        old[r["run"]] = r
    allrows = [old[k] for k in list(RUNS) + [k + "_quick" for k in RUNS] if k in old]
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(allrows)
