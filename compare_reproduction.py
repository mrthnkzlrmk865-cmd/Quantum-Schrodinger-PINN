"""
Compare a reproduction run (fresh output directories) with the committed
baseline results, array by array.

    python compare_reproduction.py [--repro results/reproduction] [--base results]

For every quantity it reports whether the reproduction is bit-for-bit
identical and the maximum absolute difference, plus the summary-statistic
differences. It also checks that the scientific trends (monotonic decrease
of <T> with W, minimum of <T> vs zeta) are the same in both. Output:
<repro>/comparison.json and a printed table.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os

import numpy as np


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def compare_arrays(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return dict(identical=False, shape_base=list(a.shape), shape_repro=list(b.shape), max_abs_diff=None)
    same = bool(np.array_equal(a, b, equal_nan=True))
    d = float(np.nanmax(np.abs(a.astype(float) - b.astype(float)))) if a.size else 0.0
    return dict(identical=same, max_abs_diff=d)


def compare_npz_dir(base_dir, rep_dir, keys):
    out = {}
    for rp in sorted(glob.glob(os.path.join(rep_dir, "*.npz"))):
        name = os.path.basename(rp)
        bp = os.path.join(base_dir, name)
        if not os.path.exists(bp):
            out[name] = dict(missing_in_baseline=True)
            continue
        with np.load(bp, allow_pickle=False) as B, np.load(rp, allow_pickle=False) as Rr:
            out[name] = {k: compare_arrays(B[k], Rr[k]) for k in keys if k in B.files and k in Rr.files}
    return out


def compare_csv(base, rep, key_cols):
    if not (os.path.exists(base) and os.path.exists(rep)):
        return dict(skipped=f"missing {'baseline' if not os.path.exists(base) else 'reproduction'} file")
    with open(base) as f:
        b = list(csv.DictReader(f))
    with open(rep) as f:
        r = list(csv.DictReader(f))
    res = dict(rows_base=len(b), rows_repro=len(r), max_abs_diff={}, identical_numeric=True)
    for rb, rr in zip(b, r):
        for k, vb in rb.items():
            if k in key_cols or k not in rr:
                continue
            x, y = _num(vb), _num(rr[k])
            if x is None or y is None:
                continue
            if not (x == y or (np.isnan(x) and np.isnan(y))):
                res["identical_numeric"] = False
                d = abs(x - y)
                res["max_abs_diff"][k] = max(res["max_abs_diff"].get(k, 0.0), d)
    return res


def trend_checks(base_root, rep_root):
    out = {}
    for label, root in (("base", base_root), ("repro", rep_root)):
        p = os.path.join(root, "transport", "strength", "summary.csv")
        if os.path.exists(p):
            with open(p) as f:
                rows = sorted(csv.DictReader(f), key=lambda r: float(r["W"]))
            T = [float(r["T_mean"]) for r in rows]
            out[f"{label}_strength_T_monotonic_decreasing"] = bool(np.all(np.diff(T) < 0))
        p = os.path.join(root, "transport", "corr_length", "summary.csv")
        if os.path.exists(p):
            with open(p) as f:
                rows = sorted(csv.DictReader(f), key=lambda r: float(r["zeta"]))
            T = [float(r["T_mean"]) for r in rows]
            out[f"{label}_corr_length_argmin_zeta"] = float(rows[int(np.argmin(T))]["zeta"])
    return out


def main(rep, base):
    res = {}
    # validation
    vb, vr = os.path.join(base, "validation"), os.path.join(rep, "validation")
    if os.path.isdir(vr):
        res["validation"] = {os.path.basename(p): compare_csv(os.path.join(vb, os.path.basename(p)), p, ("case", "quantity"))
                             for p in sorted(glob.glob(os.path.join(vr, "v*.csv")))}
    # transport
    for sw in ("strength", "corr_length", "momentum", "length", "barrier", "convergence"):
        d = os.path.join(rep, "transport", sw)
        if os.path.isdir(d):
            res[f"transport/{sw}/npz"] = compare_npz_dir(os.path.join(base, "transport", sw), d,
                                                         ("times", "T", "R", "P_in", "seeds"))
            res[f"transport/{sw}/summary"] = compare_csv(os.path.join(base, "transport", sw, "summary.csv"),
                                                         os.path.join(d, "summary.csv"), ("sweep", "point", "fingerprint"))
    # localization
    d = os.path.join(rep, "localization")
    if os.path.isdir(d):
        res["localization/npz"] = compare_npz_dir(os.path.join(base, "localization"), d, ("times", "msd", "pr"))
        res["localization/summary"] = compare_csv(os.path.join(base, "localization", "summary.csv"),
                                                  os.path.join(d, "summary.csv"), ("fingerprint",))
    # PINN
    for rd in sorted(glob.glob(os.path.join(rep, "pinn*", "*", "eval_summary.json"))):
        run = os.path.basename(os.path.dirname(rd))
        label = os.path.basename(os.path.dirname(os.path.dirname(rd)))
        bd = os.path.join(base, "pinn", run)
        with open(rd) as f:
            er = json.load(f)
        with open(os.path.join(bd, "eval_summary.json")) as f:
            eb = json.load(f)
        entry = dict(eval_metrics={k: dict(base=eb[k], repro=er[k], abs_diff=abs(eb[k] - er[k])) for k in eb})
        import torch
        sb = torch.load(os.path.join(bd, "model.pt"), map_location="cpu")
        sr = torch.load(os.path.join(os.path.dirname(rd), "model.pt"), map_location="cpu")
        entry["weights"] = dict(identical=all(torch.equal(sb[k], sr[k]) for k in sb),
                                max_abs_diff=max(float((sb[k] - sr[k]).abs().max()) for k in sb))
        res[f"{label}/{run}"] = entry
    res["trends"] = trend_checks(base, rep)
    with open(os.path.join(rep, "comparison.json"), "w") as f:
        json.dump(res, f, indent=2)
    return res


def _print(res):
    for k, v in res.items():
        if k.endswith("/npz"):
            allid = all(c.get("identical", False) for e in v.values() for c in e.values() if isinstance(c, dict))
            mx = max((c.get("max_abs_diff") or 0.0) for e in v.values() for c in e.values() if isinstance(c, dict))
            print(f"{k:40s} files={len(v):3d}  bit-identical={allid}  max|diff|={mx:.3e}")
        elif k == "validation":
            for f, c in v.items():
                if "skipped" in c:
                    continue
                print(f"validation/{f:28s} identical={c['identical_numeric']}  max|diff|={c['max_abs_diff']}")
        elif k.endswith("/summary") and "skipped" in v:
            print(f"{k:40s} skipped ({v['skipped']})")
        elif k.endswith("/summary"):
            print(f"{k:40s} identical={v['identical_numeric']}  max|diff|={v['max_abs_diff']}")
        elif k.startswith("pinn"):
            print(f"{k:40s} weights identical={v['weights']['identical']} "
                  f"(max|dw|={v['weights']['max_abs_diff']:.3e}); rel_l2_final "
                  f"{v['eval_metrics']['rel_l2_final']['base']:.6f} -> {v['eval_metrics']['rel_l2_final']['repro']:.6f}")
        else:
            print(k, v)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--repro", default="results/reproduction")
    ap.add_argument("--base", default="results")
    a = ap.parse_args()
    _print(main(a.repro, a.base))
