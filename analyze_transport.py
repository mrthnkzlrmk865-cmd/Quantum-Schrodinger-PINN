"""
Quantitative analysis of the transport sweeps (run after experiment_disorder.py).

1. Strength / correlation-length / momentum sweeps: compare <T> with the
   parameter-free Born + diffusion estimate T_est = l/(L + l),
   l = (pi/2) l_tr(W, zeta, k). Reports T_sim - T_est and the ratio.
2. Length sweep: for each W, weighted least-squares fits of <T>(L) to
     Ohmic / diffusive:  T = l / (L + l)          (1 parameter)
     exponential:        T = A exp(-L / xi)        (2 parameters)
   with weights 1/SE^2, and the reduced chi^2 of each. Also reports
   var(ln T) vs L (it grows with L in the localized regime because T
   becomes broadly distributed).
   With only five lengths and L_max / l_tr ~ 6-20, this analysis can only
   say which form fits better in the measured range. It cannot prove
   localization.

Output: results/transport/analysis.json (and printed tables).
"""

from __future__ import annotations

import csv
import json
import os

import numpy as np
from scipy.optimize import curve_fit

RES = "results/transport"


def read(sweep):
    with open(os.path.join(RES, sweep, "summary.csv")) as f:
        rows = list(csv.DictReader(f))
    return [{k: (float(v) if k not in ("sweep", "point") else v) for k, v in r.items()} for r in rows]


def est(r):
    if not np.isfinite(r["l_tr_born"]):
        return 1.0
    l = 0.5 * np.pi * r["l_tr_born"]
    return l / (r["Ld"] + l)


def est_energy_averaged(r, sigma_x=0.2, Ly=1.0, n=401):
    """Same estimate averaged over the packet's longitudinal momentum
    distribution |phi(kx)|^2 ~ exp(-(kx - kx0)^2 sigma_x^2) (transverse mode 1
    fixed). Removes the confound that T(k) is convex and the packet is not
    monochromatic."""
    from disorder import born_rates
    if r["W"] == 0:
        return 1.0
    sk = 1.0 / (np.sqrt(2) * sigma_x)
    kx = np.linspace(max(r["kx0"] - 5 * sk, 0.5), r["kx0"] + 5 * sk, n)
    w = np.exp(-((kx - r["kx0"]) ** 2) / (2 * sk**2))
    k = np.hypot(kx, np.pi / Ly)
    l = 0.5 * np.pi * born_rates(r["W"], r["zeta"], k)["l_tr"]
    return float(np.sum(w * l / (r["Ld"] + l)) / np.sum(w))


def compare_born(sweep):
    out = []
    for r in read(sweep):
        if r["W"] == 0:
            continue
        e = est(r)
        out.append(dict(point=r["point"], W=r["W"], zeta=r["zeta"], kx0=r["kx0"], k_l_tr=r["k"] * r["l_tr_born"],
                        T_sim=r["T_mean"], T_se=r["T_se"], T_est=e, diff=r["T_mean"] - e, ratio=r["T_mean"] / e,
                        T_est_avg=est_energy_averaged(r), ratio_avg=r["T_mean"] / est_energy_averaged(r)))
    return out


def fit_length():
    rows = read("length")
    res = {}
    for W in sorted({r["W"] for r in rows}):
        rr = sorted([r for r in rows if r["W"] == W], key=lambda r: r["Ld"])
        L = np.array([r["Ld"] for r in rr])
        T = np.array([r["T_mean"] for r in rr])
        se = np.array([r["T_se"] for r in rr])
        ohm = lambda L, l: l / (L + l)
        ex = lambda L, A, xi: A * np.exp(-L / xi)
        p1, c1 = curve_fit(ohm, L, T, p0=[0.5], sigma=se, absolute_sigma=True)
        p2, c2 = curve_fit(ex, L, T, p0=[1.0, 1.0], sigma=se, absolute_sigma=True)
        chi1 = float(np.sum(((T - ohm(L, *p1)) / se) ** 2) / (len(L) - 1))
        chi2 = float(np.sum(((T - ex(L, *p2)) / se) ** 2) / (len(L) - 2))
        res[str(W)] = dict(
            L=L.tolist(), T_mean=T.tolist(), T_se=se.tolist(),
            lnT_mean=[r["lnT_mean"] for r in rr], lnT_var=[r["lnT_std"] ** 2 for r in rr],
            P_in=[r["P_in_mean"] for r in rr],
            ohmic_l=float(p1[0]), ohmic_l_err=float(np.sqrt(c1[0, 0])), ohmic_chi2_red=chi1,
            exp_A=float(p2[0]), exp_xi=float(p2[1]), exp_xi_err=float(np.sqrt(c2[1, 1])), exp_chi2_red=chi2,
            born_l=float(0.5 * np.pi * rr[0]["l_tr_born"]), born_l_tr=float(rr[0]["l_tr_born"]),
        )
    return res


if __name__ == "__main__":
    out = {}
    for sw in ("strength", "corr_length", "momentum"):
        if os.path.exists(os.path.join(RES, sw, "summary.csv")):
            out[sw] = compare_born(sw)
            print(f"\n== {sw}: simulated <T> vs Born+diffusion estimate")
            for r in out[sw]:
                print(f"  {r['point']:>12s}  k*l_tr={r['k_l_tr']:7.1f}  T_sim={r['T_sim']:.3f}±{r['T_se']:.3f}  "
                      f"T_est={r['T_est']:.3f}  ratio={r['ratio']:.2f}  T_est_avg={r['T_est_avg']:.3f}  ratio_avg={r['ratio_avg']:.2f}")
    if os.path.exists(os.path.join(RES, "length", "summary.csv")):
        out["length"] = fit_length()
        print("\n== length sweep fits")
        for W, r in out["length"].items():
            print(f"  W={W}: Ohmic l={r['ohmic_l']:.3f}±{r['ohmic_l_err']:.3f} (Born l={r['born_l']:.3f}) "
                  f"chi2_red={r['ohmic_chi2_red']:.2f} | exp xi={r['exp_xi']:.3f}±{r['exp_xi_err']:.3f} "
                  f"chi2_red={r['exp_chi2_red']:.2f}")
            print(f"      var(lnT) = {np.round(r['lnT_var'], 3).tolist()}  P_in = {np.format_float_scientific(max(r['P_in']), 2)} max")
    with open(os.path.join(RES, "analysis.json"), "w") as f:
        json.dump(out, f, indent=2)
