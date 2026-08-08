"""
Central configuration for the Schrodinger PINN project.

All scripts (training, reference solver, evaluation, visualization) import
this single source of truth so that the PINN and the classical reference
solver are always compared on an identical physical setup.
"""

import torch


class Config:
    # ------------------------------------------------------------------
    # Domain (2D infinite square well: [0, L] x [0, L], time in [0, T])
    # ------------------------------------------------------------------
    L: float = 1.0
    T: float = 1.0

    # ------------------------------------------------------------------
    # Initial Gaussian wave packet
    # ------------------------------------------------------------------
    x0, y0 = 0.3, 0.5
    sigma = 0.06
    kx0, ky0 = 40.0, 0.0

    # ------------------------------------------------------------------
    # Potential. "tunneling" = single rectangular barrier.
    # "disordered" = rectangular barrier + spatially random roughness,
    # used in the disordered-transport experiment (see experiment_disorder.py)
    # ------------------------------------------------------------------
    SCENARIO: str = "tunneling"
    barrier_center = 0.55
    barrier_width = 0.06
    barrier_height = 600.0

    # Disorder experiment parameters (only used when SCENARIO == "disordered")
    disorder_strength = 150.0
    disorder_correlation_length = 0.03
    disorder_seed = 0

    # ------------------------------------------------------------------
    # Network architecture
    # ------------------------------------------------------------------
    hidden_layers = 4
    hidden_units = 64

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------
    n_interior = 4000
    n_ic = 1500
    n_bc = 1500
    n_norm = 1600           # -> 40x40 fixed grid for the normalization-loss quadrature
    n_norm_times = 4        # number of random time-slices used to estimate norm(t)

    epochs_adam = 4000
    lr_adam = 1e-3

    use_lbfgs = True
    epochs_lbfgs = 400      # L-BFGS steps (each step does several internal line-search evals)

    # NOTE: the values above are the recommended settings for a full run
    # (best on GPU, ~15-20 min). For the demo run baked into this repo's
    # figures/ and REPORT.md, a CPU-only sandbox was used with reduced
    # epochs (see README "Reproducing the results" section for exact
    # numbers) so it completes in a few minutes. Convergence quality scales
    # with epoch count -- for coursework/competition-grade results, run
    # with the defaults above on a GPU.

    # ic weight raised from 10 -> 50: with importance-sampled IC points
    # (see pinn_core.sample_points), the packet's peak carries much more
    # gradient signal than before, and a higher weight lets the network
    # prioritize matching it precisely instead of settling for a low
    # "mostly background" MSE. See README "Debugging note" for the story.
    # ic/bc are now satisfied by construction (see model.py hard-constraint
    # parametrization) so their loss weight is 0 -- they are still computed
    # and logged during training purely as a sanity check (they should sit
    # at/near machine precision for ic, and near the packet's negligible
    # tail amplitude for bc).
    loss_weights = dict(pde=1.0, ic=0.0, bc=0.0, norm=10.0)

    # PDE-loss warm-up: with Fourier features, autograd's second derivatives
    # through the embedding start out very large for a randomly initialized
    # network (empirically, mean(f_u^2+f_v^2) ~ 1e7 at init -- see README
    # "Debugging note"), which can dominate the total loss before IC/BC/norm
    # get a chance to shape the solution. Ramping the PDE weight in linearly
    # over the first `pde_warmup_epochs` gives the network a chance to first
    # match the initial condition and boundary before being pushed hard to
    # also satisfy the PDE everywhere.
    pde_warmup_epochs = 300
    grad_clip_norm = 5.0

    resample_every = 1000   # refresh collocation points every N Adam epochs

    seed = 42
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # ------------------------------------------------------------------
    # Evaluation / reference solver grid
    # ------------------------------------------------------------------
    eval_n_grid = 128       # spatial grid points per axis for reference solver
    eval_n_steps = 4000     # time steps for split-step propagation
    eval_n_frames = 60      # number of frames sampled for comparison/visualization


cfg = Config()

torch.manual_seed(cfg.seed)
