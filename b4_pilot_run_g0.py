"""B4 limited diagnostic pilot — real wiring of gate G0 (phase T6).

Classification: DIAGNOSTIC / NOT B4.

Supplies b4_pilot_checks.G0Deps with the repository's real dependencies. The ONLY solver call is the M1 control
solve (stationary.scattering_matrix) inside `m1_control`, used for the G0 control reproduction (§1.4, §1.6). This
module is NOT executed during T6; `main` refuses to run without the explicit --execute flag. It writes nothing:
evidence stays in memory until the single T5 publication.
"""
from __future__ import annotations

import csv
import os
import sys
import time

import b4_pilot_checks as chk
import b4_pilot_io as pio
import b4_pilot_rgf_records as rr
import validate_rgf_b2 as b2
import validate_rgf_b3 as b3
from stationary import build_device, scattering_matrix

J_TO_POINT = {-2: "B3-E1", -1: "C1", 0: "B3-E2", 1: "B3-E3", 2: "C2"}


def read_protocol(repo):
    with open(os.path.join(repo, pio.PROTOCOL_PATH), "r", encoding="utf-8") as f:
        return f.read()


def make_csv_reader(repo):
    def read_csv(rel):
        with open(os.path.join(repo, rel), "r", encoding="utf-8", newline="") as f:
            return list(csv.DictReader(f))
    return read_csv


def expected_energies(config):
    """{point: (kx, E)} from the existing B3 construction (validate_rgf_b3.case_energies)."""
    return {J_TO_POINT[e["j"]]: (e["kx"], e["E"]) for e in b3.case_energies(b3.case_config(config))}


def m1_control(config, point, inv_dx, E, kx, code_files):
    """One M1 control solve on the B3 device (buffer 0.2 -> 103/512, N_y 63, N_q 127) -> T7 RGFRecord."""
    cfg, pot = b3.case_potential(config)
    dx = 1.0 / inv_dx
    geo = b3.geometry(cfg, dx, b3.BASE["buffer"])
    rec = b3.Recorder(pot.numpy)
    dev = build_device(rec, geo["x_start"], geo["x_end"], dx, cfg.Ly, 63, 127)
    t0 = time.time()
    out = scattering_matrix(E, dev)                       # the only solver call in T6 (G0 controls)
    wall = time.time() - t0
    return rr.extract_record(config=config, point=point, inv_dx=inv_dx, method="M1", n_y=63, n_q=127, E=E, kx=kx,
                             solver_out=out, V_lead_max=dev.V_lead_max, potential_sha256=b3.sha(rec.calls[0]),
                             n_slices=dev.n_slices, x_start=geo["x_start"], wall_time_s=wall, code_files=code_files)


def real_deps(repo) -> chk.G0Deps:
    return chk.G0Deps(protocol_text=read_protocol, grid=b2, read_csv=make_csv_reader(repo),
                      expected_energies=expected_energies,
                      code_files=lambda r: pio.code_fingerprint(r)["code_files"], m1_control=m1_control)


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--execute" not in argv:
        print("b4_pilot_run_g0: refusing to run without --execute (gate G0 performs 52 physical M1 solves)",
              file=sys.stderr)
        return 2
    repo = os.path.dirname(os.path.abspath(__file__))
    run = chk.run_g0(repo, real_deps(repo))
    print(run.report.status)
    return 0 if run.report.passed else 1


if __name__ == "__main__":                                    # pragma: no cover
    sys.exit(main())
