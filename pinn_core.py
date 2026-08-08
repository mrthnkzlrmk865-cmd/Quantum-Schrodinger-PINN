"""
Physics components of the PINN: initial condition, PDE residual (via
autograd), collocation-point sampling, and the total training loss.

Key addition vs. the original draft: a NORMALIZATION loss term.

Why it matters
---------------
Without it, u = v = 0 everywhere is a near-optimal solution: the PDE
residual is exactly zero, the boundary loss is exactly zero, and only the
initial-condition loss pushes back against it. If w_ic isn't large enough
relative to the rest, gradient descent can collapse toward this trivial
"empty universe" solution, especially early in training. The normalization
term directly penalizes any global collapse (or blow-up) of total
probability by enforcing integral(|Psi|^2) dx dy = 1 at several random time
slices throughout training, independent of the initial condition.

This does not replace the post-hoc conservation check in evaluate.py --
that script verifies the *trained* model actually conserves probability
over time. This loss term is what makes it plausible that it will.
"""

import torch

from config import cfg
from potential import potential_torch
from initial_condition import initial_condition_torch as initial_condition





def pde_residual(model, x, y, t):
    x.requires_grad_(True)
    y.requires_grad_(True)
    t.requires_grad_(True)

    u, v = model(x, y, t)

    u_t = torch.autograd.grad(u, t, grad_outputs=torch.ones_like(u), create_graph=True)[0]
    v_t = torch.autograd.grad(v, t, grad_outputs=torch.ones_like(v), create_graph=True)[0]

    u_x = torch.autograd.grad(u, x, grad_outputs=torch.ones_like(u), create_graph=True)[0]
    u_y = torch.autograd.grad(u, y, grad_outputs=torch.ones_like(u), create_graph=True)[0]
    v_x = torch.autograd.grad(v, x, grad_outputs=torch.ones_like(v), create_graph=True)[0]
    v_y = torch.autograd.grad(v, y, grad_outputs=torch.ones_like(v), create_graph=True)[0]

    u_xx = torch.autograd.grad(u_x, x, grad_outputs=torch.ones_like(u_x), create_graph=True)[0]
    u_yy = torch.autograd.grad(u_y, y, grad_outputs=torch.ones_like(u_y), create_graph=True)[0]
    v_xx = torch.autograd.grad(v_x, x, grad_outputs=torch.ones_like(v_x), create_graph=True)[0]
    v_yy = torch.autograd.grad(v_y, y, grad_outputs=torch.ones_like(v_y), create_graph=True)[0]

    lap_u = u_xx + u_yy
    lap_v = v_xx + v_yy

    V = potential_torch(x, y)
    f_u = u_t + 0.5 * lap_v - V * v
    f_v = v_t - 0.5 * lap_u + V * u
    return f_u, f_v


def sample_points():
    dev = cfg.device

    x_int = torch.rand(cfg.n_interior, 1, device=dev) * cfg.L
    y_int = torch.rand(cfg.n_interior, 1, device=dev) * cfg.L
    t_int = torch.rand(cfg.n_interior, 1, device=dev) * cfg.T

    # --- Initial-condition points: importance-sampled around the packet ---
    # A purely uniform sample over [0,L]x[0,L] mostly misses the localized
    # Gaussian packet (sigma=0.06 in a unit box): only ~5% of uniform draws
    # land where the envelope is non-negligible. That starves the IC loss of
    # gradient signal about the packet's actual shape, and the network can
    # reach a deceptively low IC loss by predicting ~0 everywhere -- which
    # is exactly the "collapse to the trivial zero solution" failure mode.
    # Fix: draw most IC points from a Gaussian centered on the packet
    # (concentrated where it matters) and the rest uniformly (to still
    # enforce u=v=0 in the background, away from the packet).
    n_focus = int(cfg.n_ic * 0.7)
    n_bg = cfg.n_ic - n_focus

    x_focus = torch.randn(n_focus, 1, device=dev) * (2 * cfg.sigma) + cfg.x0
    y_focus = torch.randn(n_focus, 1, device=dev) * (2 * cfg.sigma) + cfg.y0
    x_focus = x_focus.clamp(0.0, cfg.L)
    y_focus = y_focus.clamp(0.0, cfg.L)

    x_bg = torch.rand(n_bg, 1, device=dev) * cfg.L
    y_bg = torch.rand(n_bg, 1, device=dev) * cfg.L

    x_ic = torch.cat([x_focus, x_bg], dim=0)
    y_ic = torch.cat([y_focus, y_bg], dim=0)
    t_ic = torch.zeros(cfg.n_ic, 1, device=dev)

    n_edge = cfg.n_bc // 4
    t_bc = torch.rand(cfg.n_bc, 1, device=dev) * cfg.T
    edges_x = torch.cat(
        [
            torch.zeros(n_edge, 1, device=dev),
            torch.full((n_edge, 1), cfg.L, device=dev),
            torch.rand(n_edge, 1, device=dev) * cfg.L,
            torch.rand(n_edge, 1, device=dev) * cfg.L,
        ]
    )
    edges_y = torch.cat(
        [
            torch.rand(n_edge, 1, device=dev) * cfg.L,
            torch.rand(n_edge, 1, device=dev) * cfg.L,
            torch.zeros(n_edge, 1, device=dev),
            torch.full((n_edge, 1), cfg.L, device=dev),
        ]
    )

    # Normalization-loss points: a FIXED spatial grid (not random MC
    # samples). A naive uniform Monte-Carlo estimate of integral(|Psi|^2)
    # with only a few hundred random points badly under-estimates a
    # narrow, localized function like our packet (sigma=0.06) purely from
    # sampling noise -- most random draws land where the packet has near-
    # zero amplitude. A fixed grid resolves the packet deterministically
    # regardless of the RNG. (This was misdiagnosed as a training/collapse
    # problem before the real bug -- an unnormalized IC amplitude, see
    # initial_condition.py -- was found; a coarse MC estimate of the
    # integral made it hard to tell the two apart. Both are fixed now.)
    n_side = max(4, int(cfg.n_norm**0.5))
    grid_1d = torch.linspace(0, cfg.L, n_side, device=dev)
    gx, gy = torch.meshgrid(grid_1d, grid_1d, indexing="ij")
    x_norm = gx.reshape(-1, 1)
    y_norm = gy.reshape(-1, 1)
    t_norm = torch.rand(cfg.n_norm_times, device=dev) * cfg.T

    return (
        (x_int, y_int, t_int),
        (x_ic, y_ic, t_ic),
        (edges_x, edges_y, t_bc),
        (x_norm, y_norm, t_norm),
    )


def normalization_loss(model, norm_points):
    """Monte-Carlo estimate of integral(|Psi|^2) dx dy at several times, penalized against 1."""
    x_n, y_n, t_slices = norm_points
    losses = []
    for t_val in t_slices:
        t_n = torch.full_like(x_n, float(t_val))
        u, v = model(x_n, y_n, t_n)
        prob_integral = torch.mean(u**2 + v**2) * (cfg.L**2)
        losses.append((prob_integral - 1.0) ** 2)
    return torch.stack(losses).mean()


def total_loss(model, interior, ic, bc, norm_points, weights=None):
    if weights is None:
        weights = cfg.loss_weights

    x_i, y_i, t_i = interior
    x_ic, y_ic, t_ic = ic
    x_bc, y_bc, t_bc = bc

    f_u, f_v = pde_residual(model, x_i, y_i, t_i)
    loss_pde = torch.mean(f_u**2) + torch.mean(f_v**2)

    u_pred_ic, v_pred_ic = model(x_ic, y_ic, t_ic)
    u0, v0 = initial_condition(x_ic, y_ic)
    loss_ic = torch.mean((u_pred_ic - u0) ** 2) + torch.mean((v_pred_ic - v0) ** 2)

    u_pred_bc, v_pred_bc = model(x_bc, y_bc, t_bc)
    loss_bc = torch.mean(u_pred_bc**2) + torch.mean(v_pred_bc**2)

    loss_norm = normalization_loss(model, norm_points)

    loss = (
        weights["pde"] * loss_pde
        + weights["ic"] * loss_ic
        + weights["bc"] * loss_bc
        + weights["norm"] * loss_norm
    )
    parts = dict(
        pde=loss_pde.item(), ic=loss_ic.item(), bc=loss_bc.item(), norm=loss_norm.item()
    )
    return loss, parts
