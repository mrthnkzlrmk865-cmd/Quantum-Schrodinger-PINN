"""
Seed-wise paired analysis of the correlation-length sweep.

In experiment_disorder.py, realization i of every zeta point uses disorder
seed i. For a given seed, make_disorder() draws the same normalized wave
vectors and phases at every zeta and rescales q by 1/zeta. Realizations
of different zeta points with the same index are therefore paired, and
differences between zeta points must be analysed as paired differences
d_i = T_i(zeta_a) - T_i(zeta_b), resampling seed indices (not the two
groups independently).

Contrasts (fixed before this script was run, but AFTER the unpaired results
had been seen, so this is a post-hoc re-analysis, not a confirmatory test;
see docs/PHASE2_METHOD_SELECTION.md §4):
  C1  T(0.05) - T(0.02)
  C2  T(0.05) - T(0.20)
  C3  curvature at the minimum: K = [T(0.03) + T(0.07)]/2 - T(0.05)
      (K > 0 means zeta = 0.05 lies below the mean of its neighbours),
      with its components T(0.03) - T(0.05) and T(0.07) - T(0.05)
  S1  bootstrap distribution of the grid argmin of <T>(zeta)
  S2  vertex of a parabola in ln(zeta) through the three points around the
      minimum (0.03, 0.05, 0.07), reported as k*zeta*

For each contrast: mean, SD of the paired differences, paired-bootstrap
percentile CI (95 % and Bonferroni 98.33 % for 3 primary contrasts), paired
Student-t CI, the unpaired (Welch) CI for comparison, the per-seed
correlation r, effect sizes d_z = mean(d)/sd(d) and
d_av = mean(d)/[(sd_a + sd_b)/2], a two-sided bootstrap p-value, and the
paired t and Wilcoxon signed-rank p-values, and a finite-resampling-corrected,
null-centered bootstrap tail estimate (b+1)/(B+1) (bootstrap_p_value; not an
exact finite-sample p-value), all three Holm-adjusted over the primary
contrasts C1-C3.

Input: the fingerprinted reproduction of the sweep (bit-identical to the
committed baseline). Before analysing, the script checks that all points
share the same seeds and that their configuration specs differ ONLY in
disorder_corr_length.

Usage:  python paired_analysis.py [--data results/reproduction/transport/corr_length]
Output: results/analysis/paired_corr_length.json and .csv,
        figures/paired_corr_length.png
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os

import numpy as np
from scipy import stats

from provenance import read_npz_metadata

N_BOOT = 20000
BOOT_SEED = 20261105
PRIMARY = ("C1", "C2", "C3")


# ----------------------------------------------------------------------
# Generic paired statistics (unit-tested in tests/test_paired.py)
# ----------------------------------------------------------------------
def paired_boot_indices(n, n_boot=N_BOOT, seed=BOOT_SEED):
    """One shared set of resampled seed indices, so that all contrasts and
    the argmin distribution use the same bootstrap replicates."""
    return np.random.default_rng(seed).integers(0, n, size=(n_boot, n))


def paired_contrast(values: dict, weights: dict, idx, levels=(0.95, 1 - 0.05 / 3)):
    """Contrast sum_j w_j T_j(per seed), analysed seed-wise.
    values: {label: array (n,)} aligned by seed; weights: {label: w}."""
    d = sum(w * np.asarray(values[k], float) for k, w in weights.items())
    n = d.size
    mean, sd = d.mean(), d.std(ddof=1)
    boot = d[idx].mean(axis=1)
    out = dict(n=n, mean=float(mean), sd_diff=float(sd), se=float(sd / np.sqrt(n)))
    for lv in levels:
        a = (1 - lv) / 2
        lo, hi = np.quantile(boot, [a, 1 - a])
        tq = stats.t.ppf(1 - a, n - 1)
        out[f"boot_ci_{lv:.4f}"] = [float(lo), float(hi)]
        out[f"t_ci_{lv:.4f}"] = [float(mean - tq * sd / np.sqrt(n)), float(mean + tq * sd / np.sqrt(n))]
    out.update(bootstrap_p_value(boot, mean))
    out["p_t"] = float(stats.ttest_1samp(d, 0.0).pvalue)
    out["p_wilcoxon"] = float(stats.wilcoxon(d).pvalue) if np.any(d != 0) else 1.0
    out["d_z"] = float(mean / sd) if sd > 0 else float("inf")
    bz = boot / (d[idx].std(axis=1, ddof=1) + 1e-300)
    out["d_z_boot_ci_0.95"] = [float(x) for x in np.quantile(bz, [0.025, 0.975])]
    if len(weights) == 2:
        (ka, wa), (kb, wb) = weights.items()
        a_, b_ = np.asarray(values[ka], float), np.asarray(values[kb], float)
        out["r_pairs"] = float(np.corrcoef(a_, b_)[0, 1])
        out["d_av"] = float(mean / (0.5 * (a_.std(ddof=1) + b_.std(ddof=1))))
        se_u = np.sqrt(a_.var(ddof=1) / n + b_.var(ddof=1) / n)
        df = se_u**4 / ((a_.var(ddof=1) / n) ** 2 / (n - 1) + (b_.var(ddof=1) / n) ** 2 / (n - 1))
        tq = stats.t.ppf(0.975, df)
        out["unpaired_welch_ci_0.95"] = [float(mean - tq * se_u), float(mean + tq * se_u)]
    return out


def bootstrap_p_value(boot_means, observed):
    """Finite-resampling-corrected, null-centered bootstrap tail estimate
    (two-sided), reported under the key p_boot. The (+1) correction follows
    the Monte Carlo / permutation literature (North, Curtis & Sham 2002;
    Phipson & Smyth 2010; see also Davison & Hinkley 1997). It prevents a
    zero estimate from finitely many resamples. It does NOT make this an
    exact finite-sample p-value: the null distribution is itself
    approximated by the centered bootstrap, so the bootstrap approximation
    error remains. It is also not an unbiased estimate of the true p.

    Null-centering procedure: the B bootstrap means boot_b are centered at
    the observed mean, boot_b - observed. This approximates the sampling
    distribution of (mean - true mean), and hence the null distribution of
    the mean when the true mean is 0. A replicate is "at least as extreme" as
    the observation if |boot_b - observed| >= |observed|. With b such
    replicates,
        p = (b + 1) / (B + 1).
    The smallest value this can report is 1/(B + 1) (about 5.0e-5 for
    B = 20000), so it never reports 0. Smaller tail probabilities cannot
    be resolved with B replicates."""
    boot_means = np.asarray(boot_means, float)
    B = boot_means.size
    b = int(np.sum(np.abs(boot_means - observed) >= abs(observed)))
    return dict(p_boot=(b + 1) / (B + 1), p_boot_count=b, p_boot_B=B, p_boot_resolution=1.0 / (B + 1))


def holm(pvals: dict):
    keys = sorted(pvals, key=pvals.get)
    m = len(keys)
    adj, running = {}, 0.0
    for i, k in enumerate(keys):
        running = max(running, min(1.0, (m - i) * pvals[k]))
        adj[k] = running
    return adj


def parabola_vertex(x, y):
    """Vertex of the parabola through three points; nan if not convex."""
    c = np.polyfit(x, y, 2)
    return -c[1] / (2 * c[0]) if c[0] > 0 else np.nan


# ----------------------------------------------------------------------
def load_sweep(data_dir):
    files = sorted(glob.glob(os.path.join(data_dir, "zeta*.npz")), key=lambda f: float(f.split("zeta")[-1][:-4]))
    T, seeds, specs = {}, None, {}
    for f in files:
        z = float(f.split("zeta")[-1][:-4])
        with np.load(f, allow_pickle=False) as d:
            T[z] = d["T"][-1].copy()
            s = d["seeds"].copy()
        if seeds is None:
            seeds = s
        elif not np.array_equal(seeds, s):
            raise ValueError(f"{f}: seeds differ from the other points; cannot pair")
        meta = read_npz_metadata(f)
        if meta is None:
            raise ValueError(f"{f}: no provenance metadata")
        specs[z] = meta["spec"]
    # the specs must differ only in the correlation length
    ref = next(iter(specs.values()))
    for z, sp in specs.items():
        diff = {k for k in set(ref["config"]) | set(sp["config"]) if ref["config"].get(k) != sp["config"].get(k)}
        if diff - {"disorder_corr_length"} or ref["run"] != sp["run"]:
            raise ValueError(f"zeta={z}: configuration differs in {sorted(diff)} or run parameters")
        if sp["config"]["disorder_corr_length"] != z:
            raise ValueError(f"zeta={z}: file name and stored correlation length disagree")
    k = float(np.hypot(ref["config"]["kx0"], np.pi / ref["config"]["Ly"]))
    return T, seeds, k


def analyse(data_dir):
    T, seeds, k = load_sweep(data_dir)
    zs = sorted(T)
    n = len(seeds)
    idx = paired_boot_indices(n)
    res = dict(data=data_dir, n_seeds=n, k=k, n_boot=N_BOOT, boot_seed=BOOT_SEED, zetas=zs,
               T_mean={str(z): float(T[z].mean()) for z in zs})
    c = {}
    c["C1"] = paired_contrast(T, {0.05: 1.0, 0.02: -1.0}, idx)
    c["C2"] = paired_contrast(T, {0.05: 1.0, 0.2: -1.0}, idx)
    c["C3"] = paired_contrast(T, {0.03: 0.5, 0.07: 0.5, 0.05: -1.0}, idx)
    c["C3a_T0.03_minus_T0.05"] = paired_contrast(T, {0.03: 1.0, 0.05: -1.0}, idx)
    c["C3b_T0.07_minus_T0.05"] = paired_contrast(T, {0.07: 1.0, 0.05: -1.0}, idx)
    for name, p in holm({k_: c[k_]["p_t"] for k_ in PRIMARY}).items():
        c[name]["p_t_holm"] = p
    for name, p in holm({k_: c[k_]["p_wilcoxon"] for k_ in PRIMARY}).items():
        c[name]["p_wilcoxon_holm"] = p
    for name, p in holm({k_: c[k_]["p_boot"] for k_ in PRIMARY}).items():
        c[name]["p_boot_holm"] = p
    res["contrasts"] = c
    # S1: argmin of the bootstrap mean curve
    M = np.stack([T[z][idx].mean(axis=1) for z in zs], axis=1)   # (n_boot, n_zeta)
    am = np.asarray(zs)[np.argmin(M, axis=1)]
    res["argmin_distribution"] = {str(z): float(np.mean(am == z)) for z in zs}
    # S2: parabola vertex in ln(zeta) through 0.03, 0.05, 0.07
    lx = np.log([0.03, 0.05, 0.07])
    verts = np.array([parabola_vertex(lx, M[b, [zs.index(0.03), zs.index(0.05), zs.index(0.07)]])
                      for b in range(M.shape[0])])
    ok = np.isfinite(verts)
    v0 = parabola_vertex(lx, [T[0.03].mean(), T[0.05].mean(), T[0.07].mean()])
    res["vertex_k_zeta"] = dict(estimate=float(k * np.exp(v0)),
                                boot_ci_95=[float(x) for x in k * np.exp(np.quantile(verts[ok], [0.025, 0.975]))],
                                fraction_convex=float(ok.mean()))
    return res, T


def write_outputs(res, T, out_json, out_csv, fig_path):
    os.makedirs(os.path.dirname(out_json), exist_ok=True)
    with open(out_json, "w") as f:
        json.dump(res, f, indent=2)
    rows = []
    for name, c in res["contrasts"].items():
        rows.append(dict(contrast=name, n=c["n"], mean=c["mean"], sd_diff=c["sd_diff"],
                         boot_lo95=c["boot_ci_0.9500"][0], boot_hi95=c["boot_ci_0.9500"][1],
                         boot_lo9833=c["boot_ci_0.9833"][0], boot_hi9833=c["boot_ci_0.9833"][1],
                         t_lo95=c["t_ci_0.9500"][0], t_hi95=c["t_ci_0.9500"][1],
                         welch_lo95=c.get("unpaired_welch_ci_0.95", [np.nan])[0],
                         welch_hi95=c.get("unpaired_welch_ci_0.95", [np.nan, np.nan])[1],
                         r_pairs=c.get("r_pairs", np.nan), d_z=c["d_z"], d_av=c.get("d_av", np.nan),
                         p_boot=c["p_boot"], p_boot_count=c["p_boot_count"], p_boot_B=c["p_boot_B"],
                         p_boot_holm=c.get("p_boot_holm", np.nan), p_t=c["p_t"], p_wilcoxon=c["p_wilcoxon"],
                         p_t_holm=c.get("p_t_holm", np.nan), p_wilcoxon_holm=c.get("p_wilcoxon_holm", np.nan)))
    with open(out_csv, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wr.writeheader()
        wr.writerows(rows)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    blue, orange, grey = "#2a78d6", "#eb6834", "#7a7974"
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    names = ["C1", "C2", "C3", "C3a_T0.03_minus_T0.05", "C3b_T0.07_minus_T0.05"]
    labels = ["T(.05)−T(.02)", "T(.05)−T(.20)", "curvature K", "T(.03)−T(.05)", "T(.07)−T(.05)"]
    y = np.arange(len(names))
    for i, nm in enumerate(names):
        c = res["contrasts"][nm]
        lo, hi = c["boot_ci_0.9500"]
        ax[0].errorbar(c["mean"], i, xerr=[[c["mean"] - lo], [hi - c["mean"]]], color=blue, marker="o",
                       capsize=3, label="paired bootstrap 95 %" if i == 0 else None)
        if "unpaired_welch_ci_0.95" in c:
            wl, wh = c["unpaired_welch_ci_0.95"]
            ax[0].errorbar(c["mean"], i + 0.25, xerr=[[c["mean"] - wl], [wh - c["mean"]]], color=orange,
                           marker="s", capsize=3, label="unpaired Welch 95 %" if i == 0 else None)
    ax[0].axvline(0, color=grey, lw=1, ls="--")
    ax[0].set_yticks(y, labels)
    ax[0].invert_yaxis()
    ax[0].set(xlabel="difference in transmission T", title="(a) seed-wise paired contrasts (n = %d)" % res["n_seeds"])
    ax[0].legend(fontsize=8)
    zs = res["zetas"]
    ax[1].bar([f"{z:g}" for z in zs], [res["argmin_distribution"][str(z)] for z in zs], color=blue)
    ax[1].set(xlabel="correlation length ζ", ylabel="bootstrap probability",
              title="(b) grid argmin of ⟨T⟩(ζ) on the 7-point ζ grid")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="results/reproduction/transport/corr_length")
    ap.add_argument("--out", default="results/analysis")
    a = ap.parse_args()
    res, T = analyse(a.data)
    write_outputs(res, T, os.path.join(a.out, "paired_corr_length.json"),
                  os.path.join(a.out, "paired_corr_length.csv"), "figures/paired_corr_length.png")
    for name, c in res["contrasts"].items():
        print(f"{name:24s} mean={c['mean']:+.4f}  paired95={np.round(c['boot_ci_0.9500'], 4).tolist()}  "
              f"Bonf98.3={np.round(c['boot_ci_0.9833'], 4).tolist()}  d_z={c['d_z']:+.2f}  "
              f"r={c.get('r_pairs', float('nan')):+.2f}  p_boot={c['p_boot']:.2e} (b={c['p_boot_count']}/{c['p_boot_B']})  "
              f"p_boot_holm={c.get('p_boot_holm', float('nan')):.2e}  p_t_holm={c.get('p_t_holm', float('nan')):.2g}")
    print("argmin distribution:", res["argmin_distribution"])
    print("vertex k*zeta:", res["vertex_k_zeta"])
