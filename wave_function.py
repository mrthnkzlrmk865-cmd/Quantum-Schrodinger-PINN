"""
Renders the trained Schrödinger PINN's |Psi(x,y,t)|^2 (probability density)
output as a 3Blue1Brown-inspired, dark/neon-themed, interactive 3D surface
animation.

Usage:
    python visualize_wavefunction.py
    -> open the generated "wavefunction_animation.html" file in your browser.
"""

import torch
import numpy as np
import plotly.graph_objects as go

from pinn_schrodinger import SchrodingerPINN, cfg


# ----------------------------------------------------------------------------------
# 1. LOAD THE MODEL
# ----------------------------------------------------------------------------------
def load_model(path="schrodinger_pinn.pt"):
    model = SchrodingerPINN().to(cfg.device)
    model.load_state_dict(torch.load(path, map_location=cfg.device))
    model.eval()
    return model


# ----------------------------------------------------------------------------------
# 2. COMPUTE |Psi|^2 ON A GRID (for every time step)
# ----------------------------------------------------------------------------------
def compute_probability_density(model, n_grid=90, n_frames=60):
    xs = np.linspace(0, cfg.L, n_grid)
    ys = np.linspace(0, cfg.L, n_grid)
    ts = np.linspace(0, cfg.T, n_frames)

    X, Y = np.meshgrid(xs, ys)
    x_flat = torch.tensor(X.reshape(-1, 1), dtype=torch.float32, device=cfg.device)
    y_flat = torch.tensor(Y.reshape(-1, 1), dtype=torch.float32, device=cfg.device)

    prob_frames = []
    with torch.no_grad():
        for t_val in ts:
            t_flat = torch.full_like(x_flat, float(t_val))
            u, v = model(x_flat, y_flat, t_flat)
            prob = (u ** 2 + v ** 2).cpu().numpy().reshape(n_grid, n_grid)
            prob_frames.append(prob)

    return X, Y, ts, prob_frames


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

    fig.write_html(out_html, auto_play=False)
    print(f"Animation saved -> {out_html}")
    return fig


if __name__ == "__main__":
    model = load_model("schrodinger_pinn.pt")
    X, Y, ts, prob_frames = compute_probability_density(model, n_grid=90, n_frames=60)
    build_animation(X, Y, ts, prob_frames)
