"""
TÜBİTAK 2204-A presentation figures: white-background 3D Plotly surface views of
EXISTING project data. Visualization only. This script does not run any solver,
does not train or evaluate a PINN, and does not modify any data file.

Inputs (read-only):
  results/pinn/W56_fourier_uniform_s0/eval_fields.npz   X, Y, times, V, psi_ref, psi_pinn
  results/pinn/W56_fourier_uniform_s0/eval_summary.json rel_l2_final, rel_l2_rho_final
  results/pinn/W56_fourier_uniform_s0/config.json       run parameters
  wavefunction_animation.html (repository root)          embedded PINN/DST animation frames

Outputs: visualizations/tubitak/*.html (standalone, Plotly.js inlined).
PNG/SVG are exported from these HTML files by export_static.js (headless Chromium).

Usage (from the repository root):
  python visualizations/tubitak/make_tubitak_3d.py
  node visualizations/tubitak/export_static.js
"""

import base64
import json
import os

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
RUN = os.path.join(ROOT, "results", "pinn", "W56_fourier_uniform_s0")
OUT = os.path.join(ROOT, "visualizations", "tubitak")
ANIM_SRC = os.path.join(ROOT, "wavefunction_animation.html")

FONT = "Times New Roman, Liberation Serif, Times, serif"
F_TITLE, F_AXIS, F_TICK = 26, 22, 18
F_TICK3D = 16
INK = "#1a1a1a"
GRID = "#d0d0d0"
PANE = "#f4f4f4"

with open(os.path.join(RUN, "config.json")) as f:
    CFG = json.load(f)
with open(os.path.join(RUN, "eval_summary.json")) as f:
    SUMMARY = json.load(f)

REL_PSI = SUMMARY["rel_l2_final"]        # relative L2 error of psi at t = T (stored value)
REL_RHO = SUMMARY["rel_l2_rho_final"]    # relative L2 error of |psi|^2 at t = T (stored value)
W = CFG["disorder_strength"]
ZETA = CFG["disorder_corr_length"]
SEED = CFG["disorder_seed"]
T_END = CFG["T"]


def tr(v, nd=3):
    """Turkish decimal comma for text labels."""
    return f"{v:.{nd}f}".replace(".", ",")


def nice_ticks(zmax, n=4):
    """Positive tick values for a [0, zmax] axis (the 0 tick is omitted: it collides with the
    horizontal-axis labels at the back corner of the 3D box)."""
    raw = zmax / n
    mag = 10 ** np.floor(np.log10(raw))
    step = min((m * mag for m in (1, 2, 2.5, 5, 10) if m * mag >= raw))
    return [float(v) for v in np.arange(step, zmax * 1.0001, step)]


def axis(title, rng=None, nticks=6, tickvals=None):
    d = dict(
        title=dict(text=title, font=dict(family=FONT, size=F_AXIS, color=INK)),
        tickfont=dict(family=FONT, size=F_TICK3D, color=INK),
        backgroundcolor=PANE, showbackground=True, gridcolor=GRID, zerolinecolor="#9a9a9a",
        linecolor="#555555", showline=True, nticks=nticks,
    )
    if rng is not None:
        d["range"] = rng
    if tickvals is not None:
        d["tickmode"] = "array"
        d["tickvals"] = tickvals
    return d


def scene(ztitle, zrange=None, eye=(-1.25, -1.35, 0.85), positive=True):
    return dict(
        xaxis=axis("x", [0, 1]),
        yaxis=axis("y", [0, 1]),
        zaxis=axis(ztitle, zrange, nticks=5,
                   tickvals=nice_ticks(zrange[1]) if (positive and zrange and zrange[0] == 0) else None),
        aspectmode="manual",
        aspectratio=dict(x=1, y=1, z=0.6),
        camera=dict(eye=dict(x=eye[0], y=eye[1], z=eye[2]), center=dict(x=0, y=0, z=-0.08)),
    )


def colorbar(title, x=None, length=0.72):
    cb = dict(
        title=dict(text=title, side="right", font=dict(family=FONT, size=F_AXIS, color=INK)),
        tickfont=dict(family=FONT, size=F_TICK, color=INK),
        thickness=18, len=length, outlinewidth=0.6, outlinecolor="#555555",
    )
    if x is not None:
        cb["x"] = x
    return cb


def base_layout(fig, title, subtitle, width, height):
    fig.update_layout(
        title=dict(
            text=f"{title}<br><span style='font-size:{F_TICK}px;color:#444'>{subtitle}</span>",
            font=dict(family=FONT, size=F_TITLE, color=INK), x=0.5, xanchor="center", y=0.97,
        ),
        font=dict(family=FONT, color=INK),
        paper_bgcolor="white", plot_bgcolor="white", separators=",.",
        width=width, height=height,
        margin=dict(l=10, r=10, t=95, b=10),
    )


LIGHT = dict(ambient=0.55, diffuse=0.75, specular=0.15, roughness=0.6, fresnel=0.1)


def load_fields():
    z = np.load(os.path.join(RUN, "eval_fields.npz"))
    X, Y, V = z["X"], z["Y"], z["V"]
    times_saved = z["times"][::5]           # evaluate.py saves psi_*[::5]  -> 0, 0.005, ..., 0.02
    assert np.isclose(times_saved[-1], T_END)
    rho_ref = np.abs(z["psi_ref"][-1]) ** 2
    rho_pinn = np.abs(z["psi_pinn"][-1]) ** 2
    # consistency with the stored summary (pure re-reading, no new computation of physics)
    chk_psi = np.linalg.norm(z["psi_pinn"][-1] - z["psi_ref"][-1]) / np.linalg.norm(z["psi_ref"][-1])
    chk_rho = np.linalg.norm(rho_pinn - rho_ref) / np.linalg.norm(rho_ref)
    assert abs(chk_psi - REL_PSI) < 1e-12 and abs(chk_rho - REL_RHO) < 1e-12
    return X, Y, V, rho_ref, rho_pinn


SUB_COMMON = (f"2B model, 3B yüzey gösterimi · tek düzensizlik realizasyonu "
              f"(W = {tr(W, 2)}, ζ = {tr(ZETA, 2)}, tohum {SEED})")


def write(fig, name):
    path = os.path.join(OUT, name)
    fig.write_html(path, include_plotlyjs=True, full_html=True,
                   config=dict(displaylogo=False, responsive=False,
                               toImageButtonOptions=dict(format="png", scale=3)))
    print("wrote", os.path.relpath(path, ROOT))


def fig1(X, Y, V):
    vabs = float(np.abs(V).max())
    fig = go.Figure(go.Surface(
        x=X, y=Y, z=V, colorscale="RdBu_r", cmid=0.0, cmin=-vabs, cmax=vabs,
        colorbar=colorbar("V(x, y)"), lighting=LIGHT, lightposition=dict(x=-1e4, y=-1e4, z=2e4),
        contours=dict(z=dict(show=False)), hovertemplate="x=%{x:.3f}<br>y=%{y:.3f}<br>V=%{z:.1f}<extra></extra>",
    ))
    fig.update_layout(scene=scene("V(x, y)", [-vabs, vabs], positive=False))
    base_layout(fig, "İki Boyutlu Düzensizlik Potansiyelinin 3B Yüzey Gösterimi", SUB_COMMON, 900, 760)
    write(fig, "fig1_duzensizlik_potansiyeli_3b.html")


def fig2(X, Y, rho_ref, rho_pinn):
    zmax = float(max(rho_ref.max(), rho_pinn.max())) * 1.05
    fig = make_subplots(
        rows=1, cols=2, specs=[[{"type": "scene"}, {"type": "scene"}]],
        subplot_titles=("DST referans çözümü: |ψ<sub>DST</sub>|²",
                        "PINN tahmini: |ψ<sub>PINN</sub>|²"),
        horizontal_spacing=0.02,
    )
    common = dict(colorscale="Viridis", cmin=0.0, cmax=zmax, lighting=LIGHT,
                  lightposition=dict(x=-1e4, y=-1e4, z=2e4))
    fig.add_trace(go.Surface(x=X, y=Y, z=rho_ref, showscale=False, name="DST referans",
                             hovertemplate="DST<br>x=%{x:.3f}<br>y=%{y:.3f}<br>|ψ|²=%{z:.3f}<extra></extra>",
                             **common), 1, 1)
    fig.add_trace(go.Surface(x=X, y=Y, z=rho_pinn, name="PINN",
                             colorbar=colorbar("|ψ|²", x=1.0, length=0.62),
                             hovertemplate="PINN<br>x=%{x:.3f}<br>y=%{y:.3f}<br>|ψ|²=%{z:.3f}<extra></extra>",
                             **common), 1, 2)
    sc = scene("|ψ|²", [0, zmax])
    fig.update_layout(scene=dict(sc, domain=dict(x=[0.0, 0.47], y=[0.1, 0.84])),
                      scene2=dict(sc, domain=dict(x=[0.49, 0.96], y=[0.1, 0.84])))
    for a in fig.layout.annotations:
        a.font = dict(family=FONT, size=F_AXIS, color=INK)
        a.y = 0.86
    base_layout(fig, "PINN ve DST Referans Dalga Fonksiyonu Karşılaştırması",
                f"{SUB_COMMON} · t = T = {tr(T_END, 2)} · aynı ızgara, z aralığı, renk ölçeği ve kamera",
                1500, 760)
    fig.add_annotation(
        text=(f"Tepe değerleri: DST {tr(rho_ref.max(), 2)}, PINN {tr(rho_pinn.max(), 2)} · "
              f"bağıl L2 hatası: ψ için {tr(REL_PSI)}, |ψ|² için {tr(REL_RHO)}"),
        x=0.5, y=0.02, xref="paper", yref="paper", showarrow=False,
        font=dict(family=FONT, size=F_TICK, color=INK))
    write(fig, "fig2_pinn_dst_karsilastirma_3b.html")


def fig3(X, Y, rho_ref, rho_pinn):
    absdiff = np.abs(rho_pinn - rho_ref)
    zmax = float(absdiff.max()) * 1.05
    fig = go.Figure(go.Surface(
        x=X, y=Y, z=absdiff, colorscale="YlOrRd", cmin=0.0, cmax=zmax,
        colorbar=colorbar("| |ψ<sub>PINN</sub>|² − |ψ<sub>DST</sub>|² |"),
        lighting=LIGHT, lightposition=dict(x=-1e4, y=-1e4, z=2e4),
        hovertemplate="x=%{x:.3f}<br>y=%{y:.3f}<br>mutlak fark=%{z:.3f}<extra></extra>",
    ))
    fig.update_layout(scene=scene("mutlak fark |Δρ|", [0, zmax]))
    base_layout(fig, "PINN Tahmininin DST Referansından Mutlak Farkı",
                f"{SUB_COMMON} · t = T = {tr(T_END, 2)} · ρ = |ψ|²", 900, 760)
    fig.add_annotation(
        text=(f"Bağıl L2 hatası (t = T): ψ için {tr(REL_PSI)}, ρ = |ψ|² için {tr(REL_RHO)} · "
              f"en büyük mutlak fark {tr(absdiff.max(), 2)}"),
        x=0.5, y=0.02, xref="paper", yref="paper", showarrow=False,
        font=dict(family=FONT, size=F_TICK, color=INK))
    write(fig, "fig3_pinn_dst_mutlak_fark_3b.html")
    return absdiff


# ------------------------------------------------------------------ animation
def _arr(o):
    a = np.frombuffer(base64.b64decode(o["bdata"]), dtype=o["dtype"])
    sh = o.get("shape")
    return a.reshape([int(v) for v in str(sh).split(",")]) if sh else a


def load_animation_frames():
    """Read the frames embedded in the existing root animation (no recomputation)."""
    s = open(ANIM_SRC, encoding="utf-8").read()
    dec = json.JSONDecoder()
    data, _ = dec.raw_decode(s, s.index("[", s.index("Plotly.newPlot(")))
    frames, _ = dec.raw_decode(s, s.index("[", s.index("addFrames")))
    assert [d.get("name") for d in data] == ["PINN", "Reference (DST split-step)"]
    X, Y = _arr(data[0]["x"]), _arr(data[0]["y"])
    pinn = [_arr(fr["data"][0]["z"]) for fr in frames]
    ref = [_arr(fr["data"][1]["z"]) for fr in frames]
    n = len(frames)
    times = np.linspace(0.0, T_END, n)
    # the original frame labels are f"t={t:.3f}"; confirm the time axis matches them exactly
    assert [fr["name"] for fr in frames] == [f"t={t:.3f}" for t in times]
    return X, Y, times, pinn, ref


def animation():
    X, Y, times, pinn, ref = load_animation_frames()
    zmax = max(max(p.max() for p in pinn), max(r.max() for r in ref)) * 1.05
    common = dict(colorscale="Viridis", cmin=0.0, cmax=zmax, lighting=LIGHT,
                  lightposition=dict(x=-1e4, y=-1e4, z=2e4))

    def traces(i):
        return [go.Surface(x=X, y=Y, z=ref[i], showscale=False, scene="scene", name="DST referans",
                           **common),
                go.Surface(x=X, y=Y, z=pinn[i], scene="scene2", name="PINN",
                           colorbar=colorbar("|ψ|²", x=1.0, length=0.6), **common)]

    def frame_traces(i):
        # only z changes between frames; colour/axis ranges stay as chosen with the scale buttons
        return [go.Surface(z=ref[i]), go.Surface(z=pinn[i])]

    # close-up scale: the largest value reached by either solution for t >= T/2 (no data are clipped
    # in that half of the animation; earlier, sharper frames exceed it and are cut by the z range)
    half = len(times) // 2
    z_close = max(max(p.max() for p in pinn[half:]), max(r.max() for r in ref[half:])) * 1.05

    def time_note(t):
        return dict(text=f"<b>t = {t:.4f}</b>".replace(".", ","), x=0.5, y=0.83, xref="paper", yref="paper",
                    showarrow=False, font=dict(family=FONT, size=F_AXIS + 2, color=INK),
                    bgcolor="white", bordercolor="#555555", borderwidth=1, borderpad=4)

    info = dict(
        text=(f"W = {tr(W, 2)}, ζ = {tr(ZETA, 2)}, tohum {SEED} · k<sub>x0</sub> = {tr(CFG['kx0'], 0)} · "
              f"T = {tr(T_END, 2)} · tek realizasyon · 2B model, 3B yüzey gösterimi<br>"
              f"DST = referans çözüm. PINN farkları gizlenmemiştir; son anda bağıl L2 hatası "
              f"ψ için {tr(REL_PSI)}, |ψ|² için {tr(REL_RHO)} (results/pinn/W56_fourier_uniform_s0)"),
        x=0.5, y=0.005, xref="paper", yref="paper", showarrow=False, align="center",
        font=dict(family=FONT, size=F_TICK - 2, color="#333333"))
    panel = [dict(text="<b>DST referans çözümü</b> |ψ<sub>DST</sub>|²", x=0.24, y=0.92, xref="paper",
                  yref="paper", showarrow=False, font=dict(family=FONT, size=F_AXIS, color=INK)),
             dict(text="<b>PINN tahmini</b> |ψ<sub>PINN</sub>|²", x=0.74, y=0.92, xref="paper",
                  yref="paper", showarrow=False, font=dict(family=FONT, size=F_AXIS, color=INK))]

    fig = go.Figure(data=traces(0))
    names = [f"k{i:02d}" for i in range(len(times))]     # unique names (original labels collide at 3 decimals)
    fig.frames = [go.Frame(data=frame_traces(i), traces=[0, 1], name=names[i],
                           layout=dict(annotations=panel + [info, time_note(t)]))
                  for i, t in enumerate(times)]
    sc = scene("|ψ|²", [0, zmax])
    sc_l = dict(sc, domain=dict(x=[0.0, 0.48], y=[0.24, 0.86]))
    sc_r = dict(sc, domain=dict(x=[0.50, 0.98], y=[0.24, 0.86]))
    btn = dict(family=FONT, size=F_TICK, color=INK)
    fig.update_layout(
        scene=sc_l, scene2=sc_r, annotations=panel + [info, time_note(times[0])],
        updatemenus=[dict(
            type="buttons", showactive=False, x=0.02, y=0.10, xanchor="left", yanchor="bottom",
            bgcolor="white", bordercolor="#555555", font=btn, direction="left",
            buttons=[dict(label="▶ Oynat", method="animate",
                          args=[None, dict(frame=dict(duration=120, redraw=True), fromcurrent=True,
                                           transition=dict(duration=0))]),
                     dict(label="⏸ Duraklat", method="animate",
                          args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")])],
        ), dict(
            type="buttons", showactive=True, active=0, x=0.98, y=0.93, xanchor="right", yanchor="bottom",
            bgcolor="white", bordercolor="#555555", font=dict(family=FONT, size=F_TICK - 2, color=INK),
            direction="left",
            buttons=[dict(label=f"z ölçeği: tam (0–{zmax:.0f})", method="update",
                          args=[dict(cmax=[zmax, zmax]),
                                {"scene.zaxis.range": [0, zmax], "scene2.zaxis.range": [0, zmax],
                                 "scene.zaxis.tickvals": nice_ticks(zmax), "scene2.zaxis.tickvals": nice_ticks(zmax)}]),
                     dict(label=f"z ölçeği: yakın (0–{z_close:.1f})".replace(".", ","), method="update",
                          args=[dict(cmax=[z_close, z_close]),
                                {"scene.zaxis.range": [0, z_close], "scene2.zaxis.range": [0, z_close],
                                 "scene.zaxis.tickvals": nice_ticks(z_close),
                                 "scene2.zaxis.tickvals": nice_ticks(z_close)}])],
        )],
        sliders=[dict(
            active=0, x=0.18, len=0.8, y=0.10, yanchor="bottom", pad=dict(t=0, b=0),
            currentvalue=dict(prefix="zaman t = ", font=btn), font=dict(family=FONT, size=F_TICK - 4, color=INK),
            bgcolor="#eeeeee", bordercolor="#555555", tickcolor="#555555",
            steps=[dict(method="animate", label=f"{t:.4f}".replace(".", ","),
                        args=[[names[i]], dict(mode="immediate", frame=dict(duration=0, redraw=True),
                                               transition=dict(duration=0))])
                   for i, t in enumerate(times)],
        )],
    )
    base_layout(fig, "Dalga Paketinin Zamanla Değişimi: DST Referans ve PINN", SUB_COMMON, 1500, 820)
    fig.update_layout(margin=dict(l=10, r=10, t=95, b=10))
    path = os.path.join(OUT, "wavefunction_animation_white.html")
    fig.write_html(path, include_plotlyjs=True, full_html=True, auto_play=False,
                   config=dict(displaylogo=False, responsive=False))
    print("wrote", os.path.relpath(path, ROOT), "frames:", len(times))
    return X.shape, len(times), float(zmax), float(z_close)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    X, Y, V, rho_ref, rho_pinn = load_fields()
    fig1(X, Y, V)
    fig2(X, Y, rho_ref, rho_pinn)
    d = fig3(X, Y, rho_ref, rho_pinn)
    info = animation()
    meta = dict(
        grid=list(X.shape), x_range=[float(X.min()), float(X.max())],
        t_figures=float(T_END), V_min=float(V.min()), V_max=float(V.max()),
        rho_ref_max=float(rho_ref.max()), rho_pinn_max=float(rho_pinn.max()),
        absdiff_max=float(d.max()), rel_l2_psi_final=REL_PSI, rel_l2_rho_final=REL_RHO,
        animation_grid=list(info[0]), animation_frames=info[1], animation_zmax=info[2],
        animation_z_close=info[3],
    )
    with open(os.path.join(OUT, "figure_values.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(json.dumps(meta, indent=2))
