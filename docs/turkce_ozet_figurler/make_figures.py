"""Türkçe teknik özet için grafikler. YALNIZCA commit edilmiş sonuç dosyalarını okur; hiçbir solver, referans hesabı
veya driver çalıştırmaz. Çalıştırma: python3 docs/turkce_ozet_figurler/make_figures.py (depo kökünden)."""
import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
OUT = os.path.dirname(os.path.abspath(__file__))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.grid": True, "grid.alpha": 0.3,
                     "savefig.dpi": 200, "figure.dpi": 100})
C_A, C_B, C_T = "#1f4e79", "#c0504d", "#555555"

# --- Şekil 1: B1 faz hatası (6 konfigürasyon) ---
s = json.load(open(os.path.join(ROOT, "results/rgf_b1_validation/phase_audit_hp_summary.json")))["configs"]
names = ["A", "B", "C", "D", "E", "F"]
hp = [s[k]["hp_max_abs_phase_err"] for k in names]
b1 = [s[k]["b1_phase_quantity_max"] for k in names]
lab = [f"{k}\nN_y={s[k]['N_y']}, N={s[k]['N']}" for k in names]
fig, ax = plt.subplots(figsize=(6.6, 2.9))
x = range(len(names))
ax.bar([i - 0.2 for i in x], b1, 0.4, color=C_B, label="(i) B1 faz büyüklüğü (double vs double referans)")
ax.bar([i + 0.2 for i in x], hp, 0.4, color=C_A, label="(ii) RGF vs 100 basamaklı kafes fazı")
ax.axhline(1e-10, color="k", ls="--", lw=1, label="önceden belirlenmiş kriter 1e-10")
ax.set_yscale("log"); ax.set_ylim(1e-11, 4e-10)
ax.set_xticks(list(x)); ax.set_xticklabels(lab, fontsize=7)
ax.set_ylabel("maks. |faz hatası| (rad)")
ax.legend(fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "sekil1_b1_faz.png")); plt.close(fig)

# --- Şekil 2: 63 kanal, açık kanal başına faz ve genlik hatası ---
a = json.load(open(os.path.join(ROOT, "experiments/wz_rgf_diagnostic/results_amplitude_roundoff/amplitude_roundoff_results.json")))
wz = a["runs"]["63x63:solve"]["open_channels"]
pr = a["production_original"]["channels"]
ch = sorted(wz, key=int)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.6, 2.7))
xs = [int(c) for c in ch]
a1.semilogy([i - 0.15 for i in xs], [abs(pr[c]["eps_phi_vs_L"]) for c in ch], "s", color=C_B, label="orijinal RGF (L'ye göre)")
a1.semilogy([i + 0.15 for i in xs], [max(abs(wz[c]["eps_phi_tRL"]), 1e-18) for c in ch], "o", color=C_A, label="W/Z (kendi kafesi R_W)")
a1.axhline(1e-10, color="k", ls="--", lw=1)
a1.set_title("|ε_φ| (t_RL fazı)", fontsize=9); a1.set_xlabel("açık kanal n"); a1.set_xticks(xs); a1.legend(fontsize=6.5)
a2.semilogy([i - 0.15 for i in xs], [abs(pr[c]["eps_A_vs_L"]) for c in ch], "s", color=C_B, label="orijinal RGF")
a2.semilogy([i + 0.15 for i in xs], [abs(wz[c]["eps_A_tRL"]) for c in ch], "o", color=C_A, label="W/Z")
a2.set_title("|ε_A| (|t_RL| genliği)", fontsize=9); a2.set_xlabel("açık kanal n"); a2.set_xticks(xs); a2.legend(fontsize=6.5)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "sekil2_63kanal.png")); plt.close(fig)

# --- Şekil 3: coupled-block, W/Z (solver A) uç blok hatası vs τ_end ---
rows = list(csv.DictReader(open(os.path.join(ROOT, "experiments/wz_rgf_diagnostic/results_coupled_blocks/coupled_blocks_comparisons.csv"))))
A = [r for r in rows if r["row"] == "solver A"]
fig, ax = plt.subplots(figsize=(6.6, 2.6))
ax.semilogy(range(len(A)), [float(r["max_e_glob"]) for r in A], "o", ms=3, color=C_A, label="W/Z uç blok e_glob (maks. 4 blok)")
ax.semilogy(range(len(A)), [float(r["tau_end"]) for r in A], "_", ms=6, color=C_T, label="sezgisel sınır τ_end")
ax.axhline(1e-10, color="k", ls="--", lw=1, label="1e-10 doğruluk kriteri")
ax.set_xlabel(f"execution sırası (n = {len(A)})"); ax.set_ylabel("Ĝ birimlerinde hata"); ax.set_ylim(1e-17, 1e-7)
ax.legend(fontsize=7, loc="upper center", ncol=3)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "sekil3_coupled.png")); plt.close(fig)
print("ok", len(A))
