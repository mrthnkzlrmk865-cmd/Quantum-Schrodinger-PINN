"""
Regenerate every figure in figures/ from the raw data in results/.
Captions are written to figures/CAPTIONS.md. Figures whose data are
missing are skipped with a message.

Usage:  python make_figures.py
"""

from __future__ import annotations

import csv
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from config import Config, waveguide_config
from disorder import born_rates, make_disorder

FIG = "figures"
RES = "results"
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MK = ["o", "s", "^", "D", "v", "P", "X", "*"]
CAPTIONS = {}

plt.rcParams.update({
    "font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10, "legend.fontsize": 8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
    "grid.color": "#e4e3df", "grid.linewidth": 0.6, "lines.linewidth": 1.8,
    "lines.markersize": 5, "savefig.dpi": 200, "figure.dpi": 100,
})


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


def save(fig, name, caption):
    os.makedirs(FIG, exist_ok=True)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, name), bbox_inches="tight")
    plt.close(fig)
    CAPTIONS[name] = caption
    print("wrote", name)


def errbar(ax, x, rows, key, color, marker, label):
    m = np.array([r[f"{key}_mean"] for r in rows])
    lo = np.array([r[f"{key}_ci_boot_lo"] for r in rows])
    hi = np.array([r[f"{key}_ci_boot_hi"] for r in rows])
    yerr = np.vstack([np.nan_to_num(m - lo), np.nan_to_num(hi - m)])
    ax.errorbar(x, m, yerr=yerr, color=color, marker=marker, capsize=3, label=label)
    return m


# ----------------------------------------------------------------------
def fig_validation():
    d = os.path.join(RES, "validation")
    if not os.path.exists(os.path.join(d, "v3_time_convergence.csv")):
        print("skip validation")
        return
    fig, ax = plt.subplots(1, 4, figsize=(15, 3.4))
    rows = read_csv(os.path.join(d, "v2_free_gaussian.csv"))
    for i, (n, dt) in enumerate(sorted({(r["n"], r["dt"]) for r in rows})):
        rr = [r for r in rows if r["n"] == n and r["dt"] == dt and r["t"] > 0]
        ax[0].semilogy([r["t"] for r in rr], [r["rel_l2"] for r in rr], color=C[i], marker=MK[i],
                       label=f"N={int(n)}, Δt={dt:g}")
    ax[0].set(xlabel="time t", ylabel="relative L2 error vs analytic", title="(a) free Gaussian packet")
    ax[0].legend()

    rows = read_csv(os.path.join(d, "v3_time_convergence.csv"))
    for i, (case, q) in enumerate([("pinn_W112", "rel_l2_psi(T)"), ("waveguide_W200", "|T-T_ref|"),
                                   ("waveguide_W200", "|R-R_ref|")]):
        rr = [r for r in rows if r["case"] == case and r["quantity"] == q]
        x = np.array([r["dt"] for r in rr])
        ax[1].loglog(x, [r["error"] for r in rr], color=C[i], marker=MK[i], label=f"{case}: {q}")
    xs = np.array([1e-5, 8e-4])
    ax[1].loglog(xs, 3e-5 * (xs / 1e-4) ** 2 * 0.3, color="#7a7974", ls="--", lw=1, label="slope 2")
    ax[1].set(xlabel="time step Δt", ylabel="error vs finest Δt", title="(b) time convergence (strong disorder)")
    ax[1].legend()

    rows = read_csv(os.path.join(d, "v4_space_convergence.csv"))
    i = 0
    for case in sorted({r["case"] for r in rows if r["case"].startswith("waveguide")}):
        for sd in (0, 1):
            rr = sorted([r for r in rows if r["case"] == case and r["seed"] == sd], key=lambda r: r["dx"])
            ref = rr[0]["T"]
            ax[2].loglog([r["dx"] for r in rr[1:]], [abs(r["T"] - ref) + 1e-16 for r in rr[1:]],
                         color=C[i % 8], marker=MK[i % 8], label=f"{case.replace('waveguide_', '')}, seed {sd}: |ΔT|")
            i += 1
    rr = sorted([r for r in rows if r["case"] == "pinn_W112_psi"], key=lambda r: r["dx"])
    ax[2].loglog([r["dx"] for r in rr], [r["T"] for r in rr], color=C[i % 8], marker=MK[i % 8],
                 label="PINN problem W=112: rel. L2 ψ(T)")
    ax[2].set(xlabel="grid spacing Δx", ylabel="difference to finest grid", title="(c) spatial convergence")
    ax[2].legend(fontsize=7)

    rows = read_csv(os.path.join(d, "v5_cap_reflection.csv"))
    ax[3].semilogy([r["kx0"] for r in rows], [max(r["spurious_return"], 1e-12) for r in rows],
                   color=C[0], marker="o", label="probability returned to the left")
    ax[3].semilogy([r["kx0"] for r in rows], [max(r["remaining"], 1e-12) for r in rows],
                   color=C[1], marker="s", label="probability left in domain")
    ax[3].set(xlabel="incident wave number k_x0", ylabel="probability", title="(d) CAP quality, clean waveguide")
    ax[3].legend()
    save(fig, "validation_solver.png",
         "Validation of the DST split-step reference solver. (a) Relative L2 error against the exact "
         "free Gaussian packet (2x2 box, walls far from the packet). (b) Self-convergence in the time step "
         "with strong disorder; errors relative to the finest step; dashed line: slope 2 (Strang splitting "
         "is second order). PINN problem: wavefunction at t = T, W = E0 = 112.5; waveguide: T and R at t = 1, "
         "W = 200. (c) Spatial convergence: the same continuum disorder realization sampled on grids with "
         "dx = 1/32, 1/64 (production) and 1/128 (reference). (d) Spurious return from the complex "
         "absorbing potential for a free packet in the waveguide, vs incident wave number. Data: results/validation/.")


# ----------------------------------------------------------------------
def theory_T(Ld, l_tr, factor=np.pi / 2):
    """Order-of-magnitude diffusive estimate T ~ l/(L + l), l = (pi/2) l_tr."""
    l = factor * l_tr
    return l / (Ld + l)


def fig_strength():
    p = os.path.join(RES, "transport", "strength", "summary.csv")
    if not os.path.exists(p):
        print("skip strength")
        return
    rows = sorted(read_csv(p), key=lambda r: r["W"])
    W = np.array([r["W"] for r in rows])
    E0 = rows[0]["E0"]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    errbar(ax[0], W / E0, rows, "T", C[0], "o", "transmission ⟨T⟩")
    errbar(ax[0], W / E0, rows, "R", C[1], "s", "reflection ⟨R⟩")
    errbar(ax[0], W / E0, rows, "P_in", C[2], "^", "still in slab ⟨P_in⟩")
    Wt = np.linspace(W[W > 0].min(), W.max(), 100)
    l_tr = born_rates(Wt, rows[0]["zeta"], rows[0]["k"])["l_tr"]
    ax[0].plot(Wt / E0, theory_T(rows[0]["Ld"], l_tr), color="#7a7974", ls="--", lw=1.2,
               label="Born + diffusion estimate ℓ/(L+ℓ)")
    ax[0].set(xlabel="disorder strength W / E₀", ylabel="probability at t_end",
              title=f"(a) L_d={rows[0]['Ld']:g}, ζ={rows[0]['zeta']:g} (kζ={rows[0]['k_zeta']:.2f})")
    ax[0].legend(loc="center right")

    # distributions
    d = os.path.join(RES, "transport", "strength")
    for i, Wv in enumerate([50, 100, 200, 300]):
        f = os.path.join(d, f"W{Wv}.npz")
        if os.path.exists(f):
            T = np.load(f)["T"][-1]
            ax[1].hist(T, bins=np.linspace(0, 1, 41), color=C[i], alpha=0.55, label=f"W/E₀={Wv/E0:.2f} (n={T.size})")
    ax[1].set(xlabel="transmission T of one realization", ylabel="number of realizations",
              title="(b) distribution over disorder realizations")
    ax[1].legend()

    for i, Wv in enumerate([0, 50, 100, 200, 300]):
        f = os.path.join(d, f"W{Wv}.npz")
        if os.path.exists(f):
            z = np.load(f)
            ax[2].semilogy(z["times"], np.maximum(z["P_in"].mean(axis=1), 1e-8), color=C[i], marker=MK[i],
                           markevery=5, label=f"W/E₀={Wv/E0:.2f}")
    ax[2].set(xlabel="time t", ylabel="⟨P_in(t)⟩ (probability inside slab)", title="(c) dwell in the disordered slab")
    ax[2].legend()
    save(fig, "transport_vs_strength.png",
         "Wave-packet transport through a disordered slab (length L_d = 1, correlation length ζ = 0.05) in a "
         "hard-wall waveguide of width 1; incident transverse mode 1, kx0 = 20 (E₀ ≈ 205). (a) Ensemble means of "
         "transmission T, reflection R and probability remaining in the slab at t_end = 1, with 95% bootstrap "
         "confidence intervals over independent disorder realizations (n given in summary.csv). Dashed: "
         "order-of-magnitude estimate ℓ/(L_d+ℓ) with ℓ = (π/2)ℓ_tr from the 2D Born approximation (valid only "
         "for kℓ ≫ 1; not a fit). (b) Distribution of T over realizations. (c) Mean probability inside the slab "
         "vs time. Data: results/transport/strength/.")


def fig_corr_length():
    p = os.path.join(RES, "transport", "corr_length", "summary.csv")
    if not os.path.exists(p):
        print("skip corr_length")
        return
    rows = sorted(read_csv(p), key=lambda r: r["zeta"])
    kz = np.array([r["k_zeta"] for r in rows])
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
    errbar(ax[0], kz, rows, "T", C[0], "o", "⟨T⟩ (simulation)")
    errbar(ax[0], kz, rows, "R", C[1], "s", "⟨R⟩ (simulation)")
    zt = np.linspace(min(r["zeta"] for r in rows), max(r["zeta"] for r in rows), 200)
    l_tr = born_rates(rows[0]["W"], zt, rows[0]["k"])["l_tr"]
    ax[0].plot(zt * rows[0]["k"], theory_T(rows[0]["Ld"], l_tr), color="#7a7974", ls="--", lw=1.2,
               label="Born + diffusion estimate for T")
    ax[0].set(xscale="log", xlabel="k ζ (wave number × correlation length)", ylabel="probability at t_end",
              title=f"(a) W={rows[0]['W']:g} (W/E₀={rows[0]['W']/rows[0]['E0']:.2f}), L_d={rows[0]['Ld']:g}")
    ax[0].legend()
    ax[1].loglog(zt * rows[0]["k"], born_rates(rows[0]["W"], zt, rows[0]["k"])["l_tr"], color=C[2],
                 label="transport mean free path ℓ_tr")
    ax[1].loglog(zt * rows[0]["k"], born_rates(rows[0]["W"], zt, rows[0]["k"])["l_s"], color=C[6], ls="-.",
                 label="scattering mean free path ℓ_s")
    ax[1].axhline(rows[0]["Ld"], color="#7a7974", lw=1, ls=":", label="slab length L_d")
    ax[1].set(xlabel="k ζ", ylabel="length (units of waveguide width)", title="(b) Born-approximation mean free paths")
    ax[1].legend()
    save(fig, "transport_vs_correlation_length.png",
         "Effect of the disorder correlation length ζ at fixed RMS strength W = 100 (L_d = 1, kx0 = 20). "
         "(a) Ensemble-mean transmission and reflection with 95% bootstrap CIs vs kζ. Dashed: Born + "
         "diffusion estimate (not a fit). (b) Born-approximation mean free paths for Gaussian-correlated "
         "disorder in unbounded 2D: ℓ_tr is non-monotonic in ζ (white-noise limit ∝ ζ⁻², smooth limit ∝ ζ). "
         "Data: results/transport/corr_length/.")


def fig_length():
    p = os.path.join(RES, "transport", "length", "summary.csv")
    pa = os.path.join(RES, "transport", "analysis.json")
    if not (os.path.exists(p) and os.path.exists(pa)):
        print("skip length")
        return
    rows = read_csv(p)
    with open(pa) as f:
        fits = json.load(f)["length"]
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 3.8))
    for i, W in enumerate(sorted({r["W"] for r in rows})):
        rr = sorted([r for r in rows if r["W"] == W], key=lambda r: r["Ld"])
        L = np.array([r["Ld"] for r in rr])
        errbar(ax[0], L, rr, "T", C[i], MK[i], f"⟨T⟩, W={W:g}")
        ft = fits[str(W)]
        Lt = np.linspace(0.15, 3.2, 100)
        ax[0].plot(Lt, ft["ohmic_l"] / (Lt + ft["ohmic_l"]), color=C[i], ls="--", lw=1,
                   label=f"Ohmic fit ℓ={ft['ohmic_l']:.2f} (χ²ᵣ={ft['ohmic_chi2_red']:.1f})")
        ax[0].plot(Lt, ft["exp_A"] * np.exp(-Lt / ft["exp_xi"]), color=C[i], ls=":", lw=1.4,
                   label=f"exp. fit ξ={ft['exp_xi']:.2f} (χ²ᵣ={ft['exp_chi2_red']:.1f})")
        ax[1].plot(L, ft["lnT_var"], color=C[i], marker=MK[i], label=f"var(ln T), W={W:g}")
    ax[0].set(yscale="log", xlabel="slab length L_d", ylabel="⟨T⟩ at t_end (log scale)",
              title="(a) mean transmission vs length")
    ax[0].legend(fontsize=7)
    ax[1].set(xlabel="slab length L_d", ylabel="variance of ln T over realizations",
              title="(b) sample-to-sample fluctuations")
    ax[1].legend()
    save(fig, "transport_vs_length.png",
         "Transmission vs disordered-slab length (ζ = 0.05, kx0 = 20, 24 realizations per point, t_end = 1 + 0.5 L_d). "
         "(a) ⟨T⟩ with 95% bootstrap CIs (log scale) and weighted least-squares fits of the diffusive form ℓ/(L+ℓ) "
         "(dashed) and of an exponential A e^(−L/ξ) (dotted); χ²ᵣ = reduced chi-square. At W = 100 the diffusive form "
         "fits better; at W = 200 the exponential fits better, with ξ of the same order as the quasi-1D estimate Nℓ. "
         "(b) var(ln T): it stays ≪ |⟨ln T⟩|, unlike the broad distributions expected deep in the localized regime "
         "(the packet averages over a band of energies). Data: results/transport/length/, analysis.json.")


def fig_momentum():
    p = os.path.join(RES, "transport", "momentum", "summary.csv")
    if not os.path.exists(p):
        print("skip momentum")
        return
    rows = read_csv(p)
    dis = sorted([r for r in rows if r["W"] > 0], key=lambda r: r["kx0"])
    cl = sorted([r for r in rows if r["W"] == 0], key=lambda r: r["kx0"])
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    k = np.array([r["kx0"] for r in dis])
    errbar(ax, k, dis, "T", C[0], "o", "⟨T⟩ disordered (W=100)")
    ax.plot([r["kx0"] for r in cl], [r["T_mean"] for r in cl], color=C[1], marker="s", label="T clean waveguide")
    ax.plot(k, [theory_T(r["Ld"], r["l_tr_born"]) for r in dis], color="#7a7974", ls="--", lw=1.2,
            label="Born + diffusion estimate")
    ax.set(xlabel="incident wave number k_x0", ylabel="transmission at t_end", title="Transmission vs incident momentum")
    ax.legend()
    save(fig, "transport_vs_momentum.png",
         "Transmission through the disordered slab (W = 100, ζ = 0.05, L_d = 1) vs incident wave number, with "
         "95% bootstrap CIs; clean-waveguide values for comparison. Faster packets are scattered less (ℓ_tr grows "
         "with k). Data: results/transport/momentum/.")


def fig_barrier():
    p = os.path.join(RES, "transport", "barrier", "summary.csv")
    if not os.path.exists(p):
        print("skip barrier")
        return
    rows = read_csv(p)
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    for i, Vb in enumerate(sorted({r["barrier"] for r in rows})):
        rr = sorted([r for r in rows if r["barrier"] == Vb], key=lambda r: r["W"])
        errbar(ax, np.array([r["W"] for r in rr]), rr, "T", C[i], MK[i], f"barrier V_b={Vb:g}")
    ax.set(xlabel="disorder strength W", ylabel="⟨T⟩ at t_end", title="Barrier inside the disordered slab")
    ax.legend()
    save(fig, "transport_barrier.png",
         "Transmission through a smooth rectangular barrier (width 0.1, centred in the slab) of height below "
         "(150) or above (250) the longitudinal energy kx0²/2 = 200, with added disorder (ζ = 0.05, L_d = 1). "
         "95% bootstrap CIs. Data: results/transport/barrier/.")


def fig_convergence():
    p = os.path.join(RES, "transport", "convergence", "W100_zeta0.05.npz")
    if not os.path.exists(p):
        print("skip convergence")
        return
    from stats_utils import running_mean
    z = np.load(p)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.4))
    for i, key in enumerate(("T", "R")):
        n, m, se = running_mean(z[key][-1])
        ax[0].plot(n, m, color=C[i], label=f"running mean ⟨{key}⟩")
        ax[0].fill_between(n, m - 1.96 * se, m + 1.96 * se, color=C[i], alpha=0.2)
        ax[1].loglog(n[1:], se[1:], color=C[i], label=f"SE of ⟨{key}⟩")
    nn = np.arange(2, len(n) + 1)
    ax[1].loglog(nn, se[-1] * np.sqrt(len(n) / nn), color="#7a7974", ls="--", lw=1, label="∝ 1/√n")
    ax[0].set(xlabel="number of realizations n", ylabel="probability", title="(a) convergence of ensemble mean (±1.96 SE)")
    ax[1].set(xlabel="number of realizations n", ylabel="standard error", title="(b) standard error")
    ax[0].legend()
    ax[1].legend()
    save(fig, "realization_convergence.png",
         "Convergence of ensemble averages with the number of independent disorder realizations "
         "(W = 100, ζ = 0.05, L_d = 1). (a) Running mean with ±1.96 SE band. (b) Standard error vs n "
         "follows 1/√n; the curves for T and R coincide because T + R ≈ 1 (P_in < 10⁻³). "
         "Data: results/transport/convergence/.")


def fig_potentials():
    cfg = waveguide_config(disorder_strength=100.0)
    x = np.linspace(0, cfg.Lx, 500)
    y = np.linspace(0, cfg.Ly, 130)
    X, Y = np.meshgrid(x, y, indexing="ij")
    fig, ax = plt.subplots(3, 1, figsize=(8, 5.6), sharex=True)
    for i, z in enumerate((0.02, 0.05, 0.15)):
        d = make_disorder(100.0, z, seed=0, region=cfg.disorder_region, edge=cfg.disorder_edge)
        im = ax[i].imshow(d.numpy(X, Y).T, origin="lower", extent=[0, cfg.Lx, 0, cfg.Ly], cmap="RdBu_r",
                          vmin=-300, vmax=300, aspect="equal")
        ax[i].set(ylabel="y", title=f"ζ = {z} (same seed, W = 100)")
        ax[i].grid(False)
        for xc in (cfg.cap_width, cfg.Lx - cfg.cap_width):
            ax[i].axvline(xc, color="#52514e", ls=":", lw=1)
    ax[-1].set_xlabel("x")
    fig.colorbar(im, ax=ax, label="V_dis(x, y)", shrink=0.8)
    fig.savefig(os.path.join(FIG, "disorder_examples.png"), bbox_inches="tight")
    plt.close(fig)
    CAPTIONS["disorder_examples.png"] = (
        "Disorder realizations (seed 0, RMS W = 100) for three correlation lengths in the waveguide geometry. "
        "Disorder is confined to the slab 1.9 < x < 2.9 by a smooth window; dotted lines mark the inner edges of "
        "the absorbing layers. The same seed gives the same random draw rescaled in ζ.")
    print("wrote disorder_examples.png")


def fig_snapshots():
    from disorder import disorder_from_config
    from potential import Potential
    from reference_solver import SplitStepDST
    cfg = waveguide_config(disorder_strength=150.0)
    s = SplitStepDST(cfg, Potential(cfg, disorder_from_config(cfg, seed=0)))
    res = s.run(t_end=0.24, n_frames=4)
    fig, ax = plt.subplots(len(res["times"]), 1, figsize=(8, 6.5), sharex=True)
    for i, (t, f) in enumerate(zip(res["times"], res["frames"][:, 0])):
        ax[i].imshow(np.abs(f).T ** 2, origin="lower", extent=[0, cfg.Lx, 0, cfg.Ly], cmap="magma", aspect="equal")
        a, b = cfg.disorder_region
        for xc in (a, b):
            ax[i].axvline(xc, color="#8fd1ff", lw=1)
        ax[i].set(ylabel="y", title=f"|ψ(x,y,t)|², t = {t:.2f}")
        ax[i].grid(False)
    ax[-1].set_xlabel("x")
    save(fig, "waveguide_snapshots.png",
         "Probability density of one realization (W = 150, ζ = 0.05, seed 0) at four times. Light lines: slab "
         "edges. The packet is partly reflected and partly transmitted, and scattered into higher transverse "
         "modes inside the slab; outgoing waves are absorbed at both ends. Each panel has its own colour "
         "scale (zero = black). Reference solver, dx = 1/64.")


def fig_localization():
    p = os.path.join(RES, "localization", "summary.csv")
    if not os.path.exists(p):
        print("skip localization")
        return
    rows = read_csv(p)
    fig, ax = plt.subplots(1, 3, figsize=(14, 3.6))
    for j, L in enumerate(sorted({r["L"] for r in rows})):
        for i, W in enumerate(sorted({r["W"] for r in rows})):
            f = os.path.join(RES, "localization", f"L{L}_W{int(W)}.npz")
            if not os.path.exists(f):
                continue
            z = np.load(f)
            ax[j].plot(z["times"], z["pr"].mean(axis=1) / L**2, color=C[i], marker=MK[i], markevery=10,
                       label=f"W={W:g}")
        ax[j].set(xlabel="time t", ylabel="⟨PR(t)⟩ / box area", title=f"({'ab'[j]}) box L = {L:g}")
        ax[j].legend()
    for j, L in enumerate(sorted({r["L"] for r in rows})):
        rr = sorted([r for r in rows if r["L"] == L], key=lambda r: r["W"])
        errbar(ax[2], np.array([r["W"] for r in rr]), rr, "PR_late_over_area", C[j], MK[j], f"L = {L:g}")
    ax[2].set(xlabel="disorder strength W", ylabel="late-time PR / area", title="(c) finite-size test")
    ax[2].legend()
    save(fig, "localization_closed_box.png",
         "Participation ratio PR = (∫ρ)²/∫ρ² (effective occupied area) of a packet launched in a closed box "
         "with disorder (ζ = 0.05, kx0 = 20) for box sizes L = 1 and 2. (a,b) ensemble-mean PR/area vs time. "
         "(c) late-time (second half) average with 95% bootstrap CIs. Localization with ξ ≪ L would make "
         "PR/area drop by ~4x when L doubles; equal values indicate states extended over the box. W = 0 is a single "
         "deterministic run (no ensemble average), hence its large recurrences. Data: results/localization/.")


def fig_pinn():
    p = os.path.join(RES, "pinn", "summary.csv")
    if not os.path.exists(p):
        print("skip pinn")
        return
    rows = read_csv(p)
    fig, ax = plt.subplots(figsize=(8, 4.2))
    names = [r["run"].replace("_fourier", "").replace("_uniform", "") for r in rows]
    y = np.arange(len(rows))
    ax.barh(y - 0.2, [r["rel_l2_final"] for r in rows], height=0.38, color=C[0], label="rel. L2 of ψ at t = T")
    ax.barh(y + 0.2, [r["rel_l2_rho_final"] for r in rows], height=0.38, color=C[1], label="rel. L2 of |ψ|² at t = T")
    ax.set_yticks(y, names)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set(xlabel="relative L2 error vs DST reference (log scale)", title="PINN accuracy by configuration")
    ax.legend()
    save(fig, "pinn_errors.png",
         "Final-time (t = T) relative L2 errors of every PINN run against the DST reference solution of the same "
         "configuration and disorder realization. W in units of the potential (E₀ = 112.5); s = network seed, "
         "d = disorder seed; nofourier = plain MLP; rad = residual-adaptive sampling; T0.04 = doubled time window. "
         "Data: results/pinn/summary.csv.")

    fig, ax = plt.subplots(figsize=(6, 3.8))
    sel = ["clean_fourier_uniform_s0", "W28_fourier_uniform_s0", "W56_fourier_uniform_s0",
           "W112_fourier_uniform_s0", "clean_T0.04_fourier_uniform_s0"]
    for i, nm in enumerate(sel):
        f = os.path.join(RES, "pinn", nm, "eval_timeseries.csv")
        if os.path.exists(f):
            rr = read_csv(f)
            ax.semilogy([r["t"] for r in rr[1:]], [r["rel_l2"] for r in rr[1:]], color=C[i], marker=MK[i],
                        markevery=3, label=nm.replace("_fourier_uniform_s0", ""))
    ax.set(xlabel="time t", ylabel="relative L2 error of ψ", title="PINN error growth in time")
    ax.legend()
    save(fig, "pinn_error_vs_time.png",
         "Relative L2 error of the PINN wavefunction vs time for clean and disordered potentials "
         "(ζ = 0.08) and for a doubled time window. Data: results/pinn/<run>/eval_timeseries.csv.")

    f = os.path.join(RES, "pinn", "W56_fourier_uniform_s0", "eval_fields.npz")
    if os.path.exists(f):
        z = np.load(f)
        X, Y = z["X"], z["Y"]
        ext = [0, 1, 0, 1]
        fig, ax = plt.subplots(1, 4, figsize=(15, 3.4))
        im = ax[0].imshow(z["V"].T, origin="lower", extent=ext, cmap="RdBu_r")
        fig.colorbar(im, ax=ax[0], shrink=0.8)
        ax[0].set(title="potential V(x,y), W=56", xlabel="x", ylabel="y")
        ax[0].grid(False)
        ax[3].grid(False)
        rp, rr_ = np.abs(z["psi_pinn"][-1]) ** 2, np.abs(z["psi_ref"][-1]) ** 2
        vmax = rr_.max()
        for a, dat, ttl in ((ax[1], rr_, "reference |ψ|², t=T"), (ax[2], rp, "PINN |ψ|², t=T")):
            im = a.imshow(dat.T, origin="lower", extent=ext, cmap="magma", vmin=0, vmax=vmax)
            fig.colorbar(im, ax=a, shrink=0.8)
            a.set(title=ttl, xlabel="x")
            a.grid(False)
        im = ax[3].imshow((rp - rr_).T, origin="lower", extent=ext, cmap="RdBu_r",
                          vmin=-np.abs(rp - rr_).max(), vmax=np.abs(rp - rr_).max())
        fig.colorbar(im, ax=ax[3], shrink=0.8)
        ax[3].set(title="PINN − reference", xlabel="x")
        save(fig, "pinn_vs_reference_W56.png",
             "PINN vs DST reference for the disordered PINN problem (W = 56 = E₀/2, ζ = 0.08, disorder seed 0) "
             "at the final time T = 0.02. Both solvers use the identical potential realization.")


def fig_born_ratio():
    p = os.path.join(RES, "transport", "analysis.json")
    if not os.path.exists(p):
        print("skip born ratio")
        return
    with open(p) as f:
        a = json.load(f)
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    labels = dict(strength="vary W (ζ=0.05, k=20)", corr_length="vary ζ (W=100, k=20)",
                  momentum="vary k (W=100, ζ=0.05)")
    for i, sw in enumerate(("strength", "corr_length", "momentum")):
        if sw not in a:
            continue
        rr = a[sw]
        x = np.array([r["k_l_tr"] for r in rr])
        y = np.array([r["ratio_avg"] for r in rr])
        e = np.array([r["T_se"] / r["T_est_avg"] * 1.96 for r in rr])
        ax.errorbar(x, y, yerr=e, color=C[i], marker=MK[i], ls="none", capsize=3, label=labels[sw])
    ax.axhline(1.0, color="#7a7974", lw=1, ls="--")
    ax.set(xscale="log", xlabel="k ℓ_tr (Born approximation)",
           ylabel="⟨T⟩ simulated / estimate", title="Departure from the perturbative estimate")
    ax.legend()
    save(fig, "born_ratio_vs_kl.png",
         "Ratio of the simulated ensemble-mean transmission to the parameter-free Born + diffusion estimate "
         "⟨ℓ(k)/(L_d+ℓ(k))⟩ (averaged over the packet's momentum distribution, ℓ = (π/2)ℓ_tr), for all three "
         "sweeps, plotted against kℓ_tr. Error bars: 95% CI of the simulated mean (estimate has no error bar). "
         "The ratio approaches ~1.05 for weak scattering and grows to ~1.3 as kℓ_tr → 1, where the Born "
         "approximation is not expected to hold. Smooth disorder (kζ ≳ 2, points with ratio < 1) deviates in "
         "the opposite direction. The overall O(1) prefactor of the estimate is convention-dependent; the trend "
         "is the robust result. Data: results/transport/analysis.json.")


def write_captions():
    with open(os.path.join(FIG, "CAPTIONS.md"), "w") as f:
        f.write("# Figure captions\n\nGenerated by `python make_figures.py`.\n\n")
        for k in sorted(CAPTIONS):
            f.write(f"## {k}\n\n![{k}]({k})\n\n{CAPTIONS[k]}\n\n")


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    for fn in (fig_validation, fig_potentials, fig_snapshots, fig_strength, fig_corr_length, fig_length,
               fig_momentum, fig_barrier, fig_convergence, fig_born_ratio, fig_localization, fig_pinn):
        try:
            fn()
        except Exception as e:  # keep going, report
            print(f"FAILED {fn.__name__}: {e!r}")
    write_captions()
