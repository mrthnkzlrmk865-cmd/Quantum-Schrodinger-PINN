"""
Interactive 3D Plotly animation of |psi(x, y, t)|^2 from a trained PINN run,
optionally overlaid with the DST reference solution of the same
configuration.

Usage:
    python visualize.py results/pinn/W56_fourier_uniform_s0
    -> writes <run>/wavefunction_animation.html and
       <run>/wavefunction_comparison.html
(plotly is only needed for this script.)
"""

import os
import sys

import numpy as np
import plotly.graph_objects as go
import torch

from evaluate import pinn_on_grid
from potential import Potential
from reference_solver import SplitStepDST
from train import load_run


def compute_probability_density(model, cfg, n_grid=90, n_frames=41):
    xs = np.linspace(0, cfg.Lx, n_grid)
    ys = np.linspace(0, cfg.Ly, n_grid)
    ts = np.linspace(0, cfg.T, n_frames)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    psi = pinn_on_grid(model, X, Y, ts)
    return X, Y, ts, [np.abs(p) ** 2 for p in psi]


# ----------------------------------------------------------------------------------
# 3. DARK / NEON THEMED 3D PLOTLY ANIMATION
# ----------------------------------------------------------------------------------
NEON_COLORSCALE = [
    [0.0, "#05010d"],   # near-black - base
    [0.2, "#170b3b"],   # dark purple
    [0.45, "#3a1c7a"],  # purple
    [0.65, "#6e2fb5"],  # vivid purple
    [0.8, "#8f3ed6"],   # pink-purple
    [0.9, "#39a6ff"],   # neon blue
    [1.0, "#7ef9ff"],   # bright cyan (peaks)
]


def build_animation(X, Y, ts, prob_frames, z_max=None, out_html="wavefunction_animation.html"):
    if z_max is None:
        z_max = max(p.max() for p in prob_frames) * 1.05

    # First frame
    fig = go.Figure(
        data=[go.Surface(
            x=X, y=Y, z=prob_frames[0],
            colorscale=NEON_COLORSCALE,
            cmin=0, cmax=z_max,
            showscale=False,
            lighting=dict(ambient=0.35, diffuse=0.9, specular=1.0, roughness=0.4, fresnel=0.2),
            lightposition=dict(x=100, y=200, z=300),
            contours_z=dict(show=True, usecolormap=True, project_z=True, highlightwidth=1),
        )]
    )

    # Remaining frames (animation frames)
    frames = []
    for i, (t_val, prob) in enumerate(zip(ts, prob_frames)):
        frames.append(go.Frame(
            data=[go.Surface(x=X, y=Y, z=prob, colorscale=NEON_COLORSCALE,
                              cmin=0, cmax=z_max, showscale=False,
                              lighting=dict(ambient=0.35, diffuse=0.9, specular=1.0,
                                            roughness=0.4, fresnel=0.2),
                              lightposition=dict(x=100, y=200, z=300),
                              contours_z=dict(show=True, usecolormap=True,
                                               project_z=True, highlightwidth=1))],
            name=f"t={t_val:.3f}"
        ))
    fig.frames = frames

    # --- Dark theme / 3Blue1Brown feel ---
    fig.update_layout(
        title=dict(
            text="|Ψ(x, y, t)|² — Quantum Wavefunction Probability Density",
            font=dict(color="#e6e6fa", size=20, family="Arial"),
            x=0.5,
        ),
        paper_bgcolor="#03000a",
        plot_bgcolor="#03000a",
        scene=dict(
            xaxis=dict(title="x", backgroundcolor="#03000a", gridcolor="#2a1a4a",
                       showbackground=True, zerolinecolor="#4a2a7a", color="#c9c9e8"),
            yaxis=dict(title="y", backgroundcolor="#03000a", gridcolor="#2a1a4a",
                       showbackground=True, zerolinecolor="#4a2a7a", color="#c9c9e8"),
            zaxis=dict(title="|Ψ|²", backgroundcolor="#03000a", gridcolor="#2a1a4a",
                       showbackground=True, zerolinecolor="#4a2a7a", color="#c9c9e8",
                       range=[0, z_max]),
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.9)),
            aspectratio=dict(x=1, y=1, z=0.6),
        ),
        margin=dict(l=0, r=0, t=60, b=0),
        updatemenus=[dict(
            type="buttons",
            showactive=False,
            y=0.02, x=0.02, xanchor="left", yanchor="bottom",
            bgcolor="#1a0f33",
            bordercolor="#7a4fd6",
            font=dict(color="#e6e6fa"),
            buttons=[
                dict(label="▶ Play", method="animate",
                     args=[None, dict(frame=dict(duration=60, redraw=True),
                                       fromcurrent=True, transition=dict(duration=0))]),
                dict(label="⏸ Pause", method="animate",
                     args=[[None], dict(frame=dict(duration=0, redraw=False),
                                         mode="immediate")]),
            ],
        )],
        sliders=[dict(
            active=0,
            currentvalue=dict(prefix="t = ", font=dict(color="#e6e6fa")),
            pad=dict(t=40),
            bgcolor="#1a0f33",
            bordercolor="#7a4fd6",
            font=dict(color="#e6e6fa"),
            steps=[dict(method="animate",
                        args=[[f"t={t_val:.3f}"],
                              dict(mode="immediate",
                                   frame=dict(duration=0, redraw=True))],
                        label=f"{t_val:.2f}")
                   for t_val in ts],
        )],
    )

    fig.write_html(out_html, auto_play=False, include_plotlyjs="cdn")
    print(f"Animation saved -> {out_html}")
    return fig


def build_comparison_animation(
    X, Y, ts, prob_pinn, prob_ref, out_html="wavefunction_comparison.html"
):
    """Side-by-side neon surface (PINN) + translucent wireframe (reference)."""
    z_max = max(max(p.max() for p in prob_pinn), max(p.max() for p in prob_ref)) * 1.05

    fig = go.Figure(
        data=[
            go.Surface(
                x=X, y=Y, z=prob_pinn[0], colorscale=NEON_COLORSCALE,
                cmin=0, cmax=z_max, showscale=False, opacity=0.95, name="PINN",
            ),
            go.Surface(
                x=X, y=Y, z=prob_ref[0], colorscale=[[0, "#888888"], [1, "#dddddd"]],
                cmin=0, cmax=z_max, showscale=False, opacity=0.35, name="Reference (DST split-step)",
            ),
        ]
    )
    frames = []
    for t_val, p_pinn, p_ref in zip(ts, prob_pinn, prob_ref):
        frames.append(
            go.Frame(
                data=[
                    go.Surface(x=X, y=Y, z=p_pinn, colorscale=NEON_COLORSCALE,
                               cmin=0, cmax=z_max, showscale=False, opacity=0.95),
                    go.Surface(x=X, y=Y, z=p_ref, colorscale=[[0, "#888888"], [1, "#dddddd"]],
                               cmin=0, cmax=z_max, showscale=False, opacity=0.35),
                ],
                name=f"t={t_val:.3f}",
            )
        )
    fig.frames = frames
    fig.update_layout(
        title=dict(
            text="PINN (colour) vs validated DST reference (grey) - disordered potential W=56, zeta=0.08; PINN final rel. L2 error 0.59",
            font=dict(color="#e6e6fa", size=16, family="Arial"), x=0.5,
        ),
        paper_bgcolor="#03000a", plot_bgcolor="#03000a",
        scene=dict(
            xaxis=dict(title="x", backgroundcolor="#03000a", gridcolor="#2a1a4a",
                       showbackground=True, color="#c9c9e8"),
            yaxis=dict(title="y", backgroundcolor="#03000a", gridcolor="#2a1a4a",
                       showbackground=True, color="#c9c9e8"),
            zaxis=dict(title="|Ψ|²", backgroundcolor="#03000a", gridcolor="#2a1a4a",
                       showbackground=True, color="#c9c9e8", range=[0, z_max]),
            camera=dict(eye=dict(x=1.5, y=1.5, z=0.9)),
            aspectratio=dict(x=1, y=1, z=0.6),
        ),
        margin=dict(l=0, r=0, t=60, b=0),
        updatemenus=[dict(
            type="buttons", showactive=False, y=0.02, x=0.02,
            xanchor="left", yanchor="bottom", bgcolor="#1a0f33",
            bordercolor="#7a4fd6", font=dict(color="#e6e6fa"),
            buttons=[
                dict(label="▶ Play", method="animate",
                     args=[None, dict(frame=dict(duration=60, redraw=True),
                                       fromcurrent=True, transition=dict(duration=0))]),
                dict(label="⏸ Pause", method="animate",
                     args=[[None], dict(frame=dict(duration=0, redraw=False), mode="immediate")]),
            ],
        )],
        sliders=[dict(
            active=0, currentvalue=dict(prefix="t = ", font=dict(color="#e6e6fa")),
            pad=dict(t=40), bgcolor="#1a0f33", bordercolor="#7a4fd6",
            font=dict(color="#e6e6fa"),
            steps=[dict(method="animate", args=[[f"t={t_val:.3f}"],
                        dict(mode="immediate", frame=dict(duration=0, redraw=True))],
                        label=f"{t_val:.2f}") for t_val in ts],
        )],
    )
    fig.write_html(
    out_html,
    auto_play=False,
    include_plotlyjs="cdn",   # Plotly.js’i CDN’den yükler, dosya boyutu küçülür
    full_html=True
)


if __name__ == "__main__":
    run = sys.argv[1] if len(sys.argv) > 1 else "results/pinn/clean_fourier_uniform_s0"
    cfg, model = load_run(run)
    X, Y, ts, prob = compute_probability_density(model, cfg)
    build_animation(X, Y, ts, prob, out_html=os.path.join(run, "wavefunction_animation.html"))
    solver = SplitStepDST(cfg.replace(nx=89, ny=89), Potential(cfg))
    ref = solver.run(n_frames=41)
    prob_ref = [np.abs(f[0]) ** 2 for f in ref["frames"]]
    prob_pinn = [np.abs(p) ** 2 for p in pinn_on_grid(model, solver.X, solver.Y, ref["times"])]
    build_comparison_animation(solver.X, solver.Y, ref["times"], prob_pinn, prob_ref,
                               out_html=os.path.join(run, "wavefunction_comparison.html"))
    print("done")
