"""
Turkish-labelled, print-size versions of the six main-report figures (TÜBİTAK report R1).

Plotting only. The data, data selection, error bars and plotted quantities are exactly those of
make_figures.py (functions fig_potentials, fig_corr_length, fig_born_ratio, fig_strength, fig_length,
fig_pinn). This script reads the same stored result files; it runs no solver and no training, and it does
not modify any existing file. Changes relative to the originals: Turkish labels and titles, Turkish decimal
comma on linear axes, larger fonts relative to the canvas, and canvas sizes chosen to match the width at
which each figure is printed in the report (so 7-8 pt text in the file is 7-8 pt on paper).

The disorder field of the first figure is evaluated deterministically from (seed, W, zeta) with
disorder.make_disorder, exactly as fig_potentials does; Born mean free paths are evaluated with
disorder.born_rates, exactly as fig_corr_length / fig_strength do. Neither is a simulation.

Usage (repository root):  python visualizations/tubitak_r1/make_figures_tr.py
"""
import csv
import json
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FuncFormatter, FixedLocator, NullFormatter  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from config import waveguide_config  # noqa: E402
from disorder import born_rates, make_disorder  # noqa: E402

RES = os.path.join(ROOT, "results")
OUT = os.path.join(ROOT, "visualizations", "tubitak_r1")
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MK = ["o", "s", "^", "D", "v", "P", "X", "*"]
CM = 1 / 2.54

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 6.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#e4e3df", "grid.linewidth": 0.5, "lines.linewidth": 1.4, "lines.markersize": 3.8,
    "savefig.dpi": 300, "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
    "legend.framealpha": 0.9, "errorbar.capsize": 2.5,
})


def tr_num(v, nd=None):
    s = f"{v:.{nd}f}" if nd is not None else f"{v:g}"
    return s.replace(".", ",").replace("-", "−")


TRFMT = FuncFormatter(lambda v, _: tr_num(round(v, 6)))


def comma_axes(*axes, x=True, y=True):
    for a in axes:
        if x and a.get_xscale() == "linear":
            a.xaxis.set_major_formatter(TRFMT)
        if y and a.get_yscale() == "linear":
            a.yaxis.set_major_formatter(TRFMT)


def read_csv(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k, v in r.items():
            try:
                r[k] = float(v)
            except (TypeError, ValueError):
                pass
    return rows


def errbar(ax, x, rows, key, color, marker, label):
    m = np.array([r[f"{key}_mean"] for r in rows])
    lo = np.array([r[f"{key}_ci_boot_lo"] for r in rows])
    hi = np.array([r[f"{key}_ci_boot_hi"] for r in rows])
    yerr = np.vstack([np.nan_to_num(m - lo), np.nan_to_num(hi - m)])
    ax.errorbar(x, m, yerr=yerr, color=color, marker=marker, label=label)
    return m


def theory_T(Ld, l_tr, factor=np.pi / 2):
    l = factor * l_tr
    return l / (Ld + l)


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    print("wrote", name)


# 1 ------------------------------------------------------------------ disorder examples
def fig_potentials():
    cfg = waveguide_config(disorder_strength=100.0)
    x = np.linspace(0, cfg.Lx, 500)
    y = np.linspace(0, cfg.Ly, 130)
    X, Y = np.meshgrid(x, y, indexing="ij")
    fig, ax = plt.subplots(3, 1, figsize=(11 * CM, 8.6 * CM), sharex=True)
    for i, z in enumerate((0.02, 0.05, 0.15)):
        d = make_disorder(100.0, z, seed=0, region=cfg.disorder_region, edge=cfg.disorder_edge)
        im = ax[i].imshow(d.numpy(X, Y).T, origin="lower", extent=[0, cfg.Lx, 0, cfg.Ly], cmap="RdBu_r",
                          vmin=-300, vmax=300, aspect="equal")
        ax[i].set(ylabel="y")
        ax[i].set_title(f"ζ = {tr_num(z)} (aynı tohum, W = 100)", pad=2)
        ax[i].grid(False)
        for xc in (cfg.cap_width, cfg.Lx - cfg.cap_width):
            ax[i].axvline(xc, color="#52514e", ls=":", lw=0.9)
        comma_axes(ax[i])
    ax[-1].set_xlabel("x")
    cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cb.set_label("düzensizlik potansiyeli V$_{dis}$(x, y)")
    cb.ax.yaxis.set_major_formatter(TRFMT)
    save(fig, "disorder_examples_tr.png")


# 2 ------------------------------------------------------------------ correlation length
def fig_corr_length():
    rows = sorted(read_csv(os.path.join(RES, "transport", "corr_length", "summary.csv")), key=lambda r: r["zeta"])
    kz = np.array([r["k_zeta"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(16 * CM, 6.4 * CM))
    errbar(ax[0], kz, rows, "T", C[0], "o", "⟨T⟩ (simülasyon)")
    errbar(ax[0], kz, rows, "R", C[1], "s", "⟨R⟩ (simülasyon)")
    zt = np.linspace(min(r["zeta"] for r in rows), max(r["zeta"] for r in rows), 200)
    k = rows[0]["k"]
    l_tr = born_rates(rows[0]["W"], zt, k)["l_tr"]
    ax[0].plot(zt * k, theory_T(rows[0]["Ld"], l_tr), color="#7a7974", ls="--", lw=1.1,
               label="T için Born + difüzyon tahmini")
    ax[0].set(xscale="log", xlabel="kζ (dalga sayısı × korelasyon uzunluğu)", ylabel="t$_{end}$ anında olasılık")
    ax[0].set_title(f"(a) W = {rows[0]['W']:g} (W/E$_0$ = {tr_num(rows[0]['W'] / rows[0]['E0'], 2)}), "
                    f"L$_d$ = {rows[0]['Ld']:g}")
    ax[0].legend(loc="lower center", bbox_to_anchor=(0.42, 0.0), fontsize=6.3)
    ax[1].loglog(zt * k, born_rates(rows[0]["W"], zt, k)["l_tr"], color=C[2], label="ℓ$_{tr}$ (taşınım)")
    ax[1].loglog(zt * k, born_rates(rows[0]["W"], zt, k)["l_s"], color=C[6], ls="-.", label="ℓ$_s$ (saçılma)")
    ax[1].axhline(rows[0]["Ld"], color="#7a7974", lw=0.9, ls=":", label="L$_d$ (bölge uzunluğu)")
    ax[1].set(xlabel="kζ", ylabel="ortalama serbest yol (kılavuz genişliği birimi)")
    ax[1].set_title("(b) Born yaklaşımında ortalama serbest yollar (kuram)")
    ax[1].legend(loc="lower left")
    ticks = [0.4, 0.6, 1, 2, 3, 4]
    for a in ax:
        a.xaxis.set_major_locator(FixedLocator(ticks))
        a.xaxis.set_major_formatter(FuncFormatter(lambda v, _: tr_num(v)))
        a.xaxis.set_minor_formatter(NullFormatter())
    ax[1].yaxis.set_major_formatter(FuncFormatter(lambda v, _: tr_num(v)))
    comma_axes(ax[0], x=False)
    fig.tight_layout(w_pad=1.2)
    save(fig, "transport_vs_correlation_length_tr.png")


# 3 ------------------------------------------------------------------ Born ratio
def fig_born_ratio():
    with open(os.path.join(RES, "transport", "analysis.json")) as f:
        a = json.load(f)
    fig, ax = plt.subplots(figsize=(11 * CM, 7.0 * CM))
    labels = dict(strength="W taraması (ζ = 0,05; k$_{x0}$ = 20)", corr_length="ζ taraması (W = 100; k$_{x0}$ = 20)",
                  momentum="k taraması (W = 100; ζ = 0,05)")
    for i, sw in enumerate(("strength", "corr_length", "momentum")):
        rr = a[sw]
        x = np.array([r["k_l_tr"] for r in rr])
        y = np.array([r["ratio_avg"] for r in rr])
        e = np.array([r["T_se"] / r["T_est_avg"] * 1.96 for r in rr])
        ax.errorbar(x, y, yerr=e, color=C[i], marker=MK[i], ls="none", label=labels[sw])
    ax.axhline(1.0, color="#7a7974", lw=0.9, ls="--")
    ax.set(xscale="log", xlabel="kℓ$_{tr}$ (Born yaklaşımı)", ylabel="⟨T⟩ simülasyon / tahmin")
    ax.set_title("Born + difüzyon tahmininden sapma")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: tr_num(v)))
    comma_axes(ax, x=False)
    ax.legend(loc="upper right")
    fig.tight_layout()
    save(fig, "born_ratio_vs_kl_tr.png")


# 4 ------------------------------------------------------------------ disorder strength
def fig_strength():
    rows = sorted(read_csv(os.path.join(RES, "transport", "strength", "summary.csv")), key=lambda r: r["W"])
    W = np.array([r["W"] for r in rows])
    E0 = rows[0]["E0"]
    fig = plt.figure(figsize=(16 * CM, 11.5 * CM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])
    a0 = fig.add_subplot(gs[0, :])
    a1 = fig.add_subplot(gs[1, 0])
    a2 = fig.add_subplot(gs[1, 1])
    errbar(a0, W / E0, rows, "T", C[0], "o", "iletim ⟨T⟩")
    errbar(a0, W / E0, rows, "R", C[1], "s", "yansıma ⟨R⟩")
    errbar(a0, W / E0, rows, "P_in", C[2], "^", "bölgede kalan ⟨P$_{in}$⟩")
    Wt = np.linspace(W[W > 0].min(), W.max(), 100)
    l_tr = born_rates(Wt, rows[0]["zeta"], rows[0]["k"])["l_tr"]
    a0.plot(Wt / E0, theory_T(rows[0]["Ld"], l_tr), color="#7a7974", ls="--", lw=1.1,
            label="Born + difüzyon tahmini ℓ/(L+ℓ)")
    a0.set(xlabel="düzensizlik gücü W / E$_0$", ylabel="t$_{end}$ anında olasılık")
    a0.set_title(f"(a) L$_d$ = {rows[0]['Ld']:g}, ζ = {tr_num(rows[0]['zeta'])} (kζ = {tr_num(rows[0]['k_zeta'], 2)})")
    a0.legend(loc="center right")
    d = os.path.join(RES, "transport", "strength")
    for i, Wv in enumerate([50, 100, 200, 300]):
        f = os.path.join(d, f"W{Wv}.npz")
        T = np.load(f)["T"][-1]
        a1.hist(T, bins=np.linspace(0, 1, 41), color=C[i], alpha=0.55,
                label=f"W/E$_0$ = {tr_num(Wv / E0, 2)} (n = {T.size})")
    a1.set(xlabel="tek realizasyonun iletimi T", ylabel="realizasyon sayısı")
    a1.set_title("(b) realizasyonlar üzerindeki dağılım")
    a1.set_ylim(0, 1.55 * a1.get_ylim()[1])
    a1.legend(loc="upper center", ncol=2, fontsize=6.2, columnspacing=0.8, handlelength=1.2)
    for i, Wv in enumerate([0, 50, 100, 200, 300]):
        z = np.load(os.path.join(d, f"W{Wv}.npz"))
        a2.semilogy(z["times"], np.maximum(z["P_in"].mean(axis=1), 1e-8), color=C[i], marker=MK[i],
                    markevery=5, label=f"W/E$_0$ = {tr_num(Wv / E0, 2)}")
    a2.set(xlabel="zaman t", ylabel="bölgedeki olasılık ⟨P$_{in}$(t)⟩")
    a2.set_title("(c) düzensiz bölgede kalma")
    a2.legend(loc="upper right", ncol=2, fontsize=6.0, columnspacing=0.6, handlelength=1.4, borderaxespad=0.2)
    comma_axes(a0, a1, a2)
    fig.tight_layout(h_pad=1.0, w_pad=1.2)
    save(fig, "transport_vs_strength_tr.png")


# 5 ------------------------------------------------------------------ slab length
def fig_length():
    rows = read_csv(os.path.join(RES, "transport", "length", "summary.csv"))
    with open(os.path.join(RES, "transport", "analysis.json")) as f:
        fits = json.load(f)["length"]
    fig, ax = plt.subplots(1, 2, figsize=(16 * CM, 7.2 * CM), gridspec_kw=dict(width_ratios=[1.35, 1]))
    for i, W in enumerate(sorted({r["W"] for r in rows})):
        rr = sorted([r for r in rows if r["W"] == W], key=lambda r: r["Ld"])
        L = np.array([r["Ld"] for r in rr])
        errbar(ax[0], L, rr, "T", C[i], MK[i], f"⟨T⟩, W = {W:g}")
        ft = fits[str(W)]
        Lt = np.linspace(0.15, 3.2, 100)
        ax[0].plot(Lt, ft["ohmic_l"] / (Lt + ft["ohmic_l"]), color=C[i], ls="--", lw=0.9,
                   label=f"difüzif fit ℓ = {tr_num(ft['ohmic_l'], 2)} (χ²$_{{red}}$ = {tr_num(ft['ohmic_chi2_red'], 1)})")
        ax[0].plot(Lt, ft["exp_A"] * np.exp(-Lt / ft["exp_xi"]), color=C[i], ls=":", lw=1.3,
                   label=f"üstel fit ξ = {tr_num(ft['exp_xi'], 2)} (χ²$_{{red}}$ = {tr_num(ft['exp_chi2_red'], 1)})")
        ax[1].plot(L, ft["lnT_var"], color=C[i], marker=MK[i], label=f"var(ln T), W = {W:g}")
    ax[0].set(yscale="log", xlabel="bölge uzunluğu L$_d$", ylabel="t$_{end}$ anında ⟨T⟩ (log ölçek)")
    ax[0].set_title("(a) ortalama iletimin uzunlukla değişimi")
    ax[0].yaxis.set_major_formatter(FuncFormatter(lambda v, _: tr_num(v)))
    ax[0].legend(fontsize=6.2, loc="lower left")
    ax[1].set(xlabel="bölge uzunluğu L$_d$", ylabel="realizasyonlar üzerinde ln T varyansı")
    ax[1].set_title("(b) örnekten örneğe dalgalanmalar")
    ax[1].legend(loc="upper left")
    comma_axes(*ax)
    fig.tight_layout(w_pad=1.2)
    save(fig, "transport_vs_length_tr.png")


# 6 ------------------------------------------------------------------ PINN errors
def fig_pinn():
    rows = read_csv(os.path.join(RES, "pinn", "summary.csv"))
    fig, ax = plt.subplots(figsize=(13 * CM, 9.5 * CM))
    names = [r["run"].replace("_fourier", "").replace("_uniform", "") for r in rows]
    y = np.arange(len(rows))
    ax.barh(y - 0.2, [r["rel_l2_final"] for r in rows], height=0.38, color=C[0], label="ψ'nin bağıl L2 hatası (t = T)")
    ax.barh(y + 0.2, [r["rel_l2_rho_final"] for r in rows], height=0.38, color=C[1], label="|ψ|²'nin bağıl L2 hatası (t = T)")
    ax.set_yticks(y, names)
    ax.tick_params(axis="y", labelsize=6.8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator([0.2, 0.3, 0.5, 1.0]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: tr_num(v)))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set(xlabel="DST referansına göre bağıl L2 hatası (log ölçek)")
    ax.set_title("Yapılandırmaya göre PINN doğruluğu", pad=22)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False)
    fig.tight_layout()
    save(fig, "pinn_errors_tr.png")


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for fn in (fig_potentials, fig_corr_length, fig_born_ratio, fig_strength, fig_length, fig_pinn):
        fn()
