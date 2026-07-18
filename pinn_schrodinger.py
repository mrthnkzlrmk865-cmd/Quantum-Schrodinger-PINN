"""
Physics-Informed Neural Network (PINN) for the 2D Time-Dependent Schrödinger Equation (TDSE)
==============================================================================================

Scenario: A particle's wavefunction Psi(x, y, t) = u(x,y,t) + i*v(x,y,t) evolves
inside a 2D infinite potential well (or tunnels through a potential barrier).

Equation (hbar = 1, m = 1, natural units):
    i * dPsi/dt = -0.5 * (d2Psi/dx2 + d2Psi/dy2) + V(x,y) * Psi

After splitting into real/imaginary parts, the PDE residuals (which must go to 0) are:
    f_u = du/dt + 0.5 * (d2v/dx2 + d2v/dy2) - V * v   -> should be 0
    f_v = dv/dt - 0.5 * (d2u/dx2 + d2u/dy2) + V * u   -> should be 0

Author: Draft prepared with Claude (Scope & speed optimized)
"""

import torch
import torch.nn as nn
import numpy as np


# ----------------------------------------------------------------------------------
# 1. CONFIGURATION
# ----------------------------------------------------------------------------------
class Config:
    # Spatial and temporal domain (2D infinite well: [0, L] x [0, L])
    L = 1.0          # box width
    T = 1.0          # total simulation time

    # Initial Gaussian wave packet parameters
    x0, y0 = 0.3, 0.5      # initial center (away from the barrier, left side)
    sigma = 0.06           # packet width
    kx0, ky0 = 40.0, 0.0   # initial momentum (in x direction -> tunneling scenario)

    # Potential: SCENARIO = "well" (flat floor, only walls) or
    #            SCENARIO = "tunneling" (rectangular barrier in the middle)
    SCENARIO = "tunneling"
    barrier_center = 0.55
    barrier_width = 0.06
    barrier_height = 600.0   # moderate height so tunneling is visibly noticeable

    # Network architecture (CPU dostu hafif model)
    hidden_layers = 4
    hidden_units = 64

    # Training (Hızlı test ve donmaları önlemek için optimize edildi)
    n_interior = 2000   # PDE kayıp noktası (20000 -> 2000)
    n_ic = 1000          # Başlangıç durumu noktası (4000 -> 1000)
    n_bc = 1000          # Sınır noktası (4000 -> 1000)
    epochs_adam = 3000   # Toplam adım sayısı (15000 -> 3000)
    lr = 1e-3

    device = "cuda" if torch.cuda.is_available() else "cpu"


# ÇOK ÖNEMLİ: cfg nesnesini sınıflardan önce tanımlıyoruz ki NameError hatası çözülsün.
cfg = Config()
torch.manual_seed(42)
np.random.seed(42)


# ----------------------------------------------------------------------------------
# 2. POTENTIAL FUNCTION V(x, y)
# ----------------------------------------------------------------------------------
def potential(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
    """
    Since the infinite well itself forbids the region outside the box, we
    take V=0 inside the box; the walls are enforced separately via a
    Dirichlet boundary condition (Psi=0). For the tunneling scenario we add
    an extra rectangular barrier.
    """
    V = torch.zeros_like(x)
    if cfg.SCENARIO == "tunneling":
        mask = (x > cfg.barrier_center - cfg.barrier_width / 2) & \
               (x < cfg.barrier_center + cfg.barrier_width / 2)
        V = torch.where(mask, torch.full_like(x, cfg.barrier_height), V)
    return V


# ----------------------------------------------------------------------------------
# 3. INITIAL WAVE PACKET (t = 0)
# ----------------------------------------------------------------------------------
def initial_condition(x: torch.Tensor, y: torch.Tensor):
    """Normalized 2D Gaussian wave packet, moving in the kx0 direction."""
    envelope = torch.exp(-((x - cfg.x0) ** 2 + (y - cfg.y0) ** 2) / (2 * cfg.sigma ** 2))
    phase = cfg.kx0 * x + cfg.ky0 * y
    u0 = envelope * torch.cos(phase)
    v0 = envelope * torch.sin(phase)
    return u0, v0


# ----------------------------------------------------------------------------------
# 4. PINN NETWORK ARCHITECTURE
# ----------------------------------------------------------------------------------
class SchrodingerPINN(nn.Module):
    """
    Input:  (x, y, t)  -> 3 dimensions
    Output: (u, v)      -> real and imaginary parts of Psi
    """

    def __init__(self, in_dim=3, out_dim=2, hidden=cfg.hidden_units, layers=cfg.hidden_layers):
        super().__init__()
        net = [nn.Linear(in_dim, hidden), nn.Tanh()]
        for _ in range(layers - 1):
            net += [nn.Linear(hidden, hidden), nn.Tanh()]
        net += [nn.Linear(hidden, out_dim)]
        self.net = nn.Sequential(*net)

        # Xavier init (works well with tanh)
        for m in self.net:
            if isinstance(m, nn.Linear):
                nn.init.xavier_normal_(m.weight)
                nn.init.zeros_(m.bias)

    def forward(self, x, y, t):
        inp = torch.cat([x, y, t], dim=1)
        out = self.net(inp)
        u, v = out[:, 0:1], out[:, 1:2]
        return u, v


# ----------------------------------------------------------------------------------
# 5. PHYSICS LOSS (PDE Residual) — derivatives via autograd
# ----------------------------------------------------------------------------------
def pde_residual(model: SchrodingerPINN, x, y, t):
    x.requires_grad_(True)
    y.requires_grad_(True)
    t.requires_grad_(True)

    u, v = model(x, y, t)

    # First derivatives
    u_t = torch.autograd.grad(u, t, grad_outputs=torch.ones_like(u), create_graph=True)[0]
    v_t = torch.autograd.grad(v, t, grad_outputs=torch.ones_like(v), create_graph=True)[0]

    u_x = torch.autograd.grad(u, x, grad_outputs=torch.ones_like(u), create_graph=True)[0]
    u_y = torch.autograd.grad(u, y, grad_outputs=torch.ones_like(u), create_graph=True)[0]
    v_x = torch.autograd.grad(v, x, grad_outputs=torch.ones_like(v), create_graph=True)[0]
    v_y = torch.autograd.grad(v, y, grad_outputs=torch.ones_like(v), create_graph=True)[0]

    # Second derivatives (for the Laplacian)
    u_xx = torch.autograd.grad(u_x, x, grad_outputs=torch.ones_like(u_x), create_graph=True)[0]
    u_yy = torch.autograd.grad(u_y, y, grad_outputs=torch.ones_like(u_y), create_graph=True)[0]
    v_xx = torch.autograd.grad(v_x, x, grad_outputs=torch.ones_like(v_x), create_graph=True)[0]
    v_yy = torch.autograd.grad(v_y, y, grad_outputs=torch.ones_like(v_y), create_graph=True)[0]

    lap_u = u_xx + u_yy
    lap_v = v_xx + v_yy

    V = potential(x, y)

    f_u = u_t + 0.5 * lap_v - V * v
    f_v = v_t - 0.5 * lap_u + V * u

    return f_u, f_v


# ----------------------------------------------------------------------------------
# 6. SAMPLING COLLOCATION POINTS
# ----------------------------------------------------------------------------------
def sample_points():
    dev = cfg.device

    # Interior points (where the PDE residual is enforced)
    x_int = torch.rand(cfg.n_interior, 1, device=dev) * cfg.L
    y_int = torch.rand(cfg.n_interior, 1, device=dev) * cfg.L
    t_int = torch.rand(cfg.n_interior, 1, device=dev) * cfg.T

    # Initial condition points (t=0)
    x_ic = torch.rand(cfg.n_ic, 1, device=dev) * cfg.L
    y_ic = torch.rand(cfg.n_ic, 1, device=dev) * cfg.L
    t_ic = torch.zeros(cfg.n_ic, 1, device=dev)

    # Boundary points (4 walls, random times)
    n_edge = cfg.n_bc // 4
    t_bc = torch.rand(cfg.n_bc, 1, device=dev) * cfg.T

    edges_x = torch.cat([
        torch.zeros(n_edge, 1, device=dev),                # left wall x=0
        torch.full((n_edge, 1), cfg.L, device=dev),        # right wall x=L
        torch.rand(n_edge, 1, device=dev) * cfg.L,         # bottom wall (random x)
        torch.rand(n_edge, 1, device=dev) * cfg.L,         # top wall (random x)
    ])
    edges_y = torch.cat([
        torch.rand(n_edge, 1, device=dev) * cfg.L,
        torch.rand(n_edge, 1, device=dev) * cfg.L,
        torch.zeros(n_edge, 1, device=dev),
        torch.full((n_edge, 1), cfg.L, device=dev),
    ])

    return (x_int, y_int, t_int), (x_ic, y_ic, t_ic), (edges_x, edges_y, t_bc)


# ----------------------------------------------------------------------------------
# 7. TOTAL LOSS FUNCTION
# ----------------------------------------------------------------------------------
def total_loss(model, interior, ic, bc, weights=(1.0, 10.0, 10.0)):
    w_pde, w_ic, w_bc = weights
    x_i, y_i, t_i = interior
    x_ic, y_ic, t_ic = ic
    x_bc, y_bc, t_bc = bc

    # --- PDE loss ---
    f_u, f_v = pde_residual(model, x_i, y_i, t_i)
    loss_pde = torch.mean(f_u ** 2) + torch.mean(f_v ** 2)

    # --- Initial condition loss ---
    u_pred_ic, v_pred_ic = model(x_ic, y_ic, t_ic)
    u0, v0 = initial_condition(x_ic, y_ic)
    loss_ic = torch.mean((u_pred_ic - u0) ** 2) + torch.mean((v_pred_ic - v0) ** 2)

    # --- Boundary loss (Dirichlet, Psi=0) ---
    u_pred_bc, v_pred_bc = model(x_bc, y_bc, t_bc)
    loss_bc = torch.mean(u_pred_bc ** 2) + torch.mean(v_pred_bc ** 2)

    loss = w_pde * loss_pde + w_ic * loss_ic + w_bc * loss_bc
    return loss, dict(pde=loss_pde.item(), ic=loss_ic.item(), bc=loss_bc.item())


# ----------------------------------------------------------------------------------
# 8. TRAINING LOOP
# ----------------------------------------------------------------------------------
def train():
    model = SchrodingerPINN().to(cfg.device)
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.lr)
    # scheduler adım sayısını da yeni epoch limitine göre ayarladık
    scheduler = torch.optim.lr_scheduler.StepLR(optimizer, step_size=1000, gamma=0.5)

    interior, ic, bc = sample_points()

    for epoch in range(1, cfg.epochs_adam + 1):
        optimizer.zero_grad()
        loss, parts = total_loss(model, interior, ic, bc)
        loss.backward()
        optimizer.step()
        scheduler.step()

        # Noktaları belirli aralıklarla yenile
        if epoch % 1000 == 0:
            interior, ic, bc = sample_points()

        # İlerleme durumunu terminale bas
        if epoch % 100 == 0 or epoch == 1:
            print(f"[{epoch:6d}/{cfg.epochs_adam}] "
                  f"total={loss.item():.6f}  pde={parts['pde']:.6f}  "
                  f"ic={parts['ic']:.6f}  bc={parts['bc']:.6f}")

    torch.save(model.state_dict(), "schrodinger_pinn.pt")
    print("Model saved -> schrodinger_pinn.pt")
    return model


if __name__ == "__main__":
    print(f"Device: {cfg.device} | Scenario: {cfg.SCENARIO}")
    train()
