"""B4 limited diagnostic pilot — real wiring of gate G-REAL (phase T6).

Classification: DIAGNOSTIC / NOT B4.

Supplies b4_pilot_checks.GRealDeps. It constructs and samples ONLY what G-REAL needs (protocol §3.3.1, §3.3.2):
disorder objects, RGF devices (stationary.build_device: construction and projection, no solve), the TD production
grid and the V0-V4 variant grids (reference_solver.SplitStepDST CONSTRUCTOR only, to read the stored potential), the
initial state psi0, and direct Potential.numpy samples on independently built coordinates. It never calls
SplitStepDST.run, scattering_matrix or any propagation. NOT executed during T6; `main` requires --execute.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

import b4_pilot_checks as chk
import b4_pilot_io as pio
import validate_rgf_b3 as b3
from config import Config
from disorder import disorder_from_config
from initial_condition import initial_condition_numpy
from potential import Potential
from reference_solver import SplitStepDST
from stationary import build_device

TD_FILES = {"A": "results/transport/corr_length/zeta0.02.npz", "B": "results/transport/corr_length/zeta0.05.npz",
            "C": "results/transport/corr_length/zeta0.2.npz", "D": "results/transport/strength/W200.npz"}


def read_protocol(repo):
    with open(os.path.join(repo, pio.PROTOCOL_PATH), "r", encoding="utf-8") as f:
        return f.read()


def make_td_config(repo):
    def td_config(config):
        with np.load(os.path.join(repo, TD_FILES[config]), allow_pickle=False) as data:
            if 0 not in list(data["seeds"]):
                raise chk.PilotCheckError(f"{TD_FILES[config]} does not contain seed 0")
            return Config.from_dict(json.loads(str(data["config_json"])))
    return td_config


def make_record_identities(repo):
    def record_identities(config):
        with open(os.path.join(repo, chk.B3_METADATA), encoding="utf-8") as f:
            r3 = json.load(f)["cases"][config]
        with open(os.path.join(repo, chk.B3X_METADATA), encoding="utf-8") as f:
            r3x = json.load(f)["configurations"][chk.B3X_NAME[config]]
        return {"B3_record": r3, "B3X_record": r3x}                  # B3Y records carry no disorder hash
    return record_identities


def rgf_device(config, inv_dx, n_y, n_q):
    cfg, pot = b3.case_potential(config)
    dx = 1.0 / inv_dx
    geo = b3.geometry(cfg, dx, b3.BASE["buffer"])
    rec = b3.Recorder(pot.numpy)
    try:
        dev = build_device(rec, geo["x_start"], geo["x_end"], dx, cfg.Ly, n_y, n_q)   # construction only
    except ValueError as exc:                                                       # lead-tail guard
        raise chk.NotEvaluable(f"RGF device {config} 1/{inv_dx} {n_y}x{n_q}: {exc}") from exc
    xs, y, direct = b3.direct_reference_samples(config, geo["x_start"], dx, dev.n_slices, cfg.Ly, n_q)
    return chk.GridSample(name=f"RGF {config} 1/{inv_dx} N_y={n_y} N_q={n_q}", V_stored=rec.calls[0], V_direct=direct,
                          x_stored=dev.xs, x_direct=xs, y_stored=dev.basis.y, y_direct=y,
                          info=dict(V_lead_max=dev.V_lead_max, n_slices=dev.n_slices, x_start=geo["x_start"]))


def td_production(config, cfg):
    solver = SplitStepDST(cfg, [Potential(cfg, disorder_from_config(cfg, seed=0))])   # constructor only
    dx, dy = cfg.Lx / (cfg.nx + 1), cfg.Ly / (cfg.ny + 1)
    x, y = np.arange(1, cfg.nx + 1) * dx, np.arange(1, cfg.ny + 1) * dy
    X, Y = np.meshgrid(x, y, indexing="ij")
    direct = np.asarray(Potential(cfg, disorder_from_config(cfg, seed=0)).numpy(X, Y), dtype=np.float64)
    psi0 = initial_condition_numpy(cfg, solver.X, solver.Y)
    s = chk.GridSample(name=f"TD production {config}", V_stored=solver.V[0], V_direct=direct,
                       x_stored=solver.X[:, 0], x_direct=x, y_stored=solver.Y[0, :], y_direct=y, psi0=psi0,
                       info=dict(nx=cfg.nx, ny=cfg.ny, dx=solver.dx))
    return s, solver.dx


def td_variant(config, cfg, geo):
    m, dx = geo.spec.m, geo.dx_prod
    y = np.arange(1, cfg.ny + 1) * (cfg.Ly / (cfg.ny + 1))                    # production y-grid
    Xp, Y = np.meshgrid(geo.x_phys, y, indexing="ij")
    V = np.asarray(Potential(cfg, disorder_from_config(cfg, seed=0)).numpy(Xp, Y), dtype=np.float64)
    psi0 = initial_condition_numpy(cfg, Xp, Y)                                   # production cfg (x0 = 1.2)
    solver = SplitStepDST(geo.cfg, potentials=V)                                 # constructor only
    xd = np.array([float(j + 1 - m) for j in range(geo.cfg.nx)]) * dx            # independently rebuilt x_phys
    Xd, Yd = np.meshgrid(xd, y, indexing="ij")
    direct = np.asarray(Potential(cfg, disorder_from_config(cfg, seed=0)).numpy(Xd, Yd), dtype=np.float64)
    s = chk.GridSample(name=f"TD {geo.spec.name} {config}", V_stored=solver.V[0], V_direct=direct,
                       x_stored=geo.x_phys, x_direct=xd, y_stored=solver.Y[0, :], y_direct=y, psi0=psi0,
                       info=dict(m=m, nx=geo.cfg.nx, solver_x_first=float(solver.X[0, 0]),
                                 solver_x_last=float(solver.X[-1, 0])))
    return s, solver.dx


def real_deps(repo) -> chk.GRealDeps:
    return chk.GRealDeps(protocol_text=read_protocol, case_disorder=lambda c: b3.case_potential(c)[1].disorder,
                         case_config=b3.case_config, td_config=make_td_config(repo),
                         disorder_from_config=lambda cfg, seed: disorder_from_config(cfg, seed=seed),
                         record_identities=make_record_identities(repo), rgf_device=rgf_device,
                         td_production=td_production, td_variant=td_variant)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--execute" not in argv:
        print("b4_pilot_run_greal: refusing to run without --execute (gate G-REAL samples real potentials)",
              file=sys.stderr)
        return 2
    print("b4_pilot_run_greal: run G-REAL through b4_pilot_checks.run_gates after G0", file=sys.stderr)
    return 2


if __name__ == "__main__":                                    # pragma: no cover
    sys.exit(main())
