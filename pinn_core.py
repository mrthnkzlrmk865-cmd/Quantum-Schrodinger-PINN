"""
Physics part of the PINN: PDE residual (autograd), collocation sampling
(uniform or residual-adaptive), normalization loss and total loss.

Residual of i psi_t = -1/2 lap psi + V psi with psi = u + i v:

    real part of the equation:  f_v = v_t - 1/2 lap u + V u = 0
    imag part of the equation:  f_u = u_t + 1/2 lap v - V v = 0

The loss uses the scaled residual  (T / A) * f,  where T is the time window
and A = max|psi0|. This is the residual of the equation written in the
normalized time tau = t/T and normalized amplitude, so its natural size is
O(E0 * T), independent of units. (The original code used the raw residual,
whose mean square was ~1e7 at initialization.)

Normalization loss: (int |psi|^2 dx dy - 1)^2 averaged over a few random
times, with the integral evaluated by the midpoint rule on a fixed grid.
The hard-constrained IC already has unit norm at t = 0, so this term only
penalizes norm drift during the evolution. It is a soft constraint:
whether the trained model conserves probability is checked separately by
evaluate.py.

Sampling strategies for interior collocation points:
  "uniform": i.i.d. uniform in [0,Lx]x[0,Ly]x[0,T], redrawn every
             resample_every epochs.
  "rad":     residual-based adaptive distribution (Wu et al., CMAME 403,
             115671 (2023)): draw a pool of rad_pool_factor * n uniform
             candidates and pick n of them with probability
             proportional to |f|^2 / mean|f|^2 + 1.
"""

from __future__ import annotations

import torch

from potential import Potential


def _grad(out, inp):
    return torch.autograd.grad(out, inp, grad_outputs=torch.ones_like(out), create_graph=True)[0]


def pde_residual(model, potential: Potential, x, y, t, V=None, scaled=True):
    """Returns (f_u, f_v). V may be passed precomputed (same points)."""
    x = x.detach().requires_grad_(True)
    y = y.detach().requires_grad_(True)
    t = t.detach().requires_grad_(True)
    u, v = model(x, y, t)
    u_t, v_t = _grad(u, t), _grad(v, t)
    u_x, u_y = _grad(u, x), _grad(u, y)
    v_x, v_y = _grad(v, x), _grad(v, y)
    lap_u = _grad(u_x, x) + _grad(u_y, y)
    lap_v = _grad(v_x, x) + _grad(v_y, y)
    if V is None:
        V = potential.torch(x.detach(), y.detach())
    f_u = u_t + 0.5 * lap_v - V * v
    f_v = v_t - 0.5 * lap_u + V * u
    if scaled:
        s = model.cfg.T / model.amp
        f_u, f_v = s * f_u, s * f_v
    return f_u, f_v


def _uniform(cfg, n, gen, device):
    x = torch.rand(n, 1, generator=gen, device=device, dtype=torch.float64).float() * cfg.Lx
    y = torch.rand(n, 1, generator=gen, device=device, dtype=torch.float64).float() * cfg.Ly
    t = torch.rand(n, 1, generator=gen, device=device, dtype=torch.float64).float() * cfg.T
    return x, y, t


def sample_interior(cfg, model, potential, gen):
    dev = cfg.device
    if cfg.sampling == "uniform":
        x, y, t = _uniform(cfg, cfg.n_interior, gen, dev)
    elif cfg.sampling == "rad":
        pool = cfg.n_interior * cfg.rad_pool_factor
        xp, yp, tp = _uniform(cfg, pool, gen, dev)
        res2 = []
        for s in range(0, pool, 8192):
            fu, fv = pde_residual(model, potential, xp[s:s + 8192], yp[s:s + 8192], tp[s:s + 8192])
            res2.append((fu**2 + fv**2).detach())
        r = torch.cat(res2).squeeze(1)
        p = r / r.mean() + 1.0
        idx = torch.multinomial(p / p.sum(), cfg.n_interior, replacement=False, generator=gen)
        x, y, t = xp[idx], yp[idx], tp[idx]
    else:
        raise ValueError(cfg.sampling)
    V = potential.torch(x, y).detach()
    return x, y, t, V


def norm_quadrature_points(cfg, device):
    n = cfg.n_norm_side
    xs = (torch.arange(n, device=device, dtype=torch.float32) + 0.5) * cfg.Lx / n
    ys = (torch.arange(n, device=device, dtype=torch.float32) + 0.5) * cfg.Ly / n
    gx, gy = torch.meshgrid(xs, ys, indexing="ij")
    w = cfg.Lx * cfg.Ly / n**2
    return gx.reshape(-1, 1), gy.reshape(-1, 1), w


def normalization_loss(model, xq, yq, w, times):
    losses = []
    for t_val in times:
        u, v = model(xq, yq, torch.full_like(xq, float(t_val)))
        losses.append((torch.sum(u**2 + v**2) * w - 1.0) ** 2)
    return torch.stack(losses).mean()


def total_loss(model, potential, interior, norm_pts, weights):
    x, y, t, V = interior
    f_u, f_v = pde_residual(model, potential, x, y, t, V=V)
    loss_pde = torch.mean(f_u**2) + torch.mean(f_v**2)
    xq, yq, w, times = norm_pts
    loss_norm = normalization_loss(model, xq, yq, w, times)
    loss = weights["pde"] * loss_pde + weights["norm"] * loss_norm
    return loss, dict(pde=loss_pde.item(), norm=loss_norm.item())
