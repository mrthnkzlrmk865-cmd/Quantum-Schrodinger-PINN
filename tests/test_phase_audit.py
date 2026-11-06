"""High-precision phase audit of the B1 RGF results (audit_rgf_phase.py).

Checks the audit machinery (phase wrapping, reference convergence) and the
committed artifacts (results/rgf_b1_validation/phase_audit_hp*): all six
B1 configurations present, the energy grid identical to the B1 grid in
b1_grid.csv, expected row counts, configuration F evaluated, and summary
values consistent with the CSV. One row of configuration F is recomputed
live. Thresholds are not changed here.
"""
import csv
import json
import os

import mpmath as mp
import numpy as np
import pytest

import audit_rgf_phase as ap

RES = "results/rgf_b1_validation"
CSV = os.path.join(RES, "phase_audit_hp.csv")
SUMMARY = os.path.join(RES, "phase_audit_hp_summary.json")
META = os.path.join(RES, "phase_audit_hp_metadata.json")
HAVE = os.path.exists(CSV) and os.path.exists(SUMMARY)


# ----------------------------------------------------------------------
# machinery
# ----------------------------------------------------------------------
def test_six_configurations_match_validate_rgf_b1():
    assert ap.CONFIGS == {"A": (63, 1 / 512, 1.0), "B": (127, 1 / 512, 1.0), "C": (63, 1 / 1024, 1.0),
                          "D": (127, 1 / 1024, 1.0), "E": (63, 1 / 512, 3.0), "F": (63, 1 / 1024, 3.0)}
    assert ap.PHASE_THRESHOLD == 1e-10          # B1 threshold, unchanged


def test_wrap_to_pi():
    for x in (0.0, 1.0, np.pi - 1e-12, -np.pi + 1e-12, 7.0, -7.0, 1000.3):
        y = ap.wrap_to_pi(x)
        assert -np.pi < y <= np.pi and np.isclose(np.exp(1j * y), np.exp(1j * x), atol=1e-12)
    assert ap.wrap_to_pi(-np.pi) == np.pi
    with mp.workdps(50):
        x = mp.mpf(1000) + mp.mpf("0.3")
    y = ap.wrap_to_pi(x, dps=50)
    with mp.workdps(50):
        assert -mp.pi < y <= mp.pi and abs(mp.sin(y) - mp.sin(mp.mpf("1000.3"))) < mp.mpf(10) ** -45
    with pytest.raises(ValueError):
        ap.wrap_to_pi(x)                                   # precision must be explicit for mpf input


def test_wrapped_reference_independent_of_ambient_precision():
    args = (125.69231487269683, 0.5 * (5 * np.pi) ** 2, 1 / 1024, 3483, 100)
    a = ap.reference_phase_wrapped(*args)
    with mp.workdps(15):
        b = ap.reference_phase_wrapped(*args)
    with mp.workdps(300):
        c = ap.reference_phase_wrapped(*args)
    assert a == b == c


def test_wrapped_difference_is_branch_independent_and_correct_near_pi():
    for k in (-3, 0, 5, 400):
        phi = 0.7 + 2 * np.pi * k
        assert abs(ap.wrapped_difference(np.exp(1j * 0.7), phi)) < 1e-12
    # across the +/- pi cut: arg(z) just below +pi, reference just above -pi
    z = np.exp(1j * (np.pi - 1e-9))
    d = ap.wrapped_difference(z, -np.pi + 1e-9)
    assert abs(d - (-2e-9)) < 1e-15                       # naive subtraction would give ~2 pi
    assert abs(np.angle(z) - (-np.pi + 1e-9)) > 6.0
    # sign convention: positive when the RGF phase leads the reference
    assert ap.wrapped_difference(np.exp(1j * 0.3), 0.3 - 1e-10) == pytest.approx(1e-10, abs=1e-16)


@pytest.mark.parametrize("E,n,dx,N", [(125.69231487269683, 5, 1 / 1024, 3483),   # F worst-case channel
                                      (9.579321918704377, 1, 1 / 512, 718),
                                      (240.0, 6, 1 / 1024, 1435)])
def test_reference_converges_with_precision_and_formula(E, n, dx, N):
    Ep = 0.5 * (n * np.pi) ** 2
    ref = ap.reference_phase(E, Ep, dx, N, 150)
    for dps, tol in ((30, 1e-20), (50, 1e-40), (100, 1e-90)):
        assert abs(ap.reference_phase(E, Ep, dx, N, dps) - ref) < tol
    assert abs(ap.reference_phase(E, Ep, dx, N, 100, form="asin") - ref) < 1e-90


def test_reference_rejects_closed_channel():
    with pytest.raises(ValueError):
        ap.reference_phase(10.0, 0.5 * (2 * np.pi) ** 2, 1 / 512, 718, 50)


# ----------------------------------------------------------------------
# committed artifacts
# ----------------------------------------------------------------------
def _rows():
    with open(CSV) as f:
        return list(csv.DictReader(f))


def _b1_grid():
    with open(os.path.join(RES, "b1_grid.csv")) as f:
        return list(csv.DictReader(f))


@pytest.mark.skipif(not HAVE, reason="audit artifacts not generated")
def test_all_six_configurations_present_and_F_evaluated():
    rows = _rows()
    assert sorted({r["config"] for r in rows}) == list("ABCDEF")
    F = [r for r in rows if r["config"] == "F"]
    assert F and {int(r["N"]) for r in F} == {3483}
    with open(SUMMARY) as f:
        s = json.load(f)
    assert sorted(s["configs"]) == list("ABCDEF") and s["configs"]["F"]["N"] == 3483


@pytest.mark.skipif(not HAVE, reason="audit artifacts not generated")
def test_energy_grid_identical_to_b1_grid_and_row_counts():
    rows, grid = _rows(), _b1_grid()
    for cid, (ny, dx, Ld) in ap.CONFIGS.items():
        g = [r for r in grid if int(r["N_y"]) == ny and float(r["dx"]) == dx and float(r["L_d"]) == Ld]
        a = [r for r in rows if r["config"] == cid]
        assert len(g) == 50
        assert [float(r["E"]) for r in g] == sorted({float(r["E"]) for r in a})       # exact equality
        # one row per open channel and direction
        assert len(a) == 2 * sum(int(r["n_open"]) for r in g)
        assert {int(r["N"]) for r in a} == {int(g[0]["N"])}


@pytest.mark.skipif(not HAVE, reason="audit artifacts not generated")
def test_b1_quantity_matches_b1_grid():
    """The recomputed predeclared B1 phase quantity equals the per-energy
    'phase' column of b1_grid.csv (max over channels and directions)."""
    rows, grid = _rows(), _b1_grid()
    for cid, (ny, dx, Ld) in ap.CONFIGS.items():
        for g in (r for r in grid if int(r["N_y"]) == ny and float(r["dx"]) == dx and float(r["L_d"]) == Ld):
            E = float(g["E"])
            m = max(float(r["b1_phase_quantity"]) for r in rows if r["config"] == cid and float(r["E"]) == E)
            # b1_grid.csv evaluated lambda^(N-1) as an array power, the audit as a
            # scalar power: agreement to rounding is required, not bit equality
            assert abs(m - float(g["phase"])) <= 1e-12 * float(g["phase"]), (cid, E)


@pytest.mark.skipif(not HAVE, reason="audit artifacts not generated")
def test_summary_agrees_with_csv():
    rows = _rows()
    with open(SUMMARY) as f:
        s = json.load(f)
    assert s["phase_threshold"] == 1e-10
    for cid, v in s["configs"].items():
        rr = [r for r in rows if r["config"] == cid]
        err = np.array([float(r["abs_phase_err"]) for r in rr])
        b1q = np.array([float(r["b1_phase_quantity"]) for r in rr])
        assert v["n_rows"] == len(rr) and v["n_energies"] == 50
        assert v["hp_max_abs_phase_err"] == err.max()
        assert np.isclose(v["hp_rms_phase_err"], np.sqrt(np.mean(err**2)), rtol=1e-12)
        i = int(np.argmax(err))
        assert v["hp_argmax"]["E"] == float(rr[i]["E"]) and v["hp_argmax"]["channel"] == int(rr[i]["channel"])
        assert v["b1_phase_quantity_max"] == b1q.max()
        assert v["b1_phase_criterion_pass"] == bool(b1q.max() <= 1e-10)
        assert v["max_ref_change_dps50_vs_dps100"] == max(float(r["ref_change_dps50_vs_dps100"]) for r in rr)
        for r in rr:   # each row: abs error is |wrapped difference|
            assert float(r["abs_phase_err"]) == abs(float(r["dphase_wrapped"]))
    assert s["all_configs_b1_phase_pass"] == all(v["b1_phase_criterion_pass"] for v in s["configs"].values())


@pytest.mark.skipif(not HAVE, reason="audit artifacts not generated")
def test_metadata_complete():
    with open(META) as f:
        m = json.load(f)
    for k in ("generated_utc", "command", "source_commit", "working_tree_dirty", "python", "numpy", "scipy",
              "mpmath", "arbitrary_precision", "code_fingerprint"):
        assert k in m
    assert m["arbitrary_precision"]["library"] == "mpmath" and m["arbitrary_precision"]["reference_dps"] == 100


@pytest.mark.skipif(not HAVE, reason="audit artifacts not generated")
def test_configuration_F_worst_row_reproduces_live():
    """Recompute the worst row of configuration F from scratch and compare
    with the committed CSV (same code, same machine class: bitwise for the
    reference, and to rounding for the RGF value)."""
    rows = [r for r in _rows() if r["config"] == "F"]
    w = max(rows, key=lambda r: float(r["abs_phase_err"]))
    from stationary import scattering_matrix
    from validate_rgf_b1 import device
    dev = device(63, 1 / 1024, 3.0)
    E, n = float(w["E"]), int(w["channel"])
    res = scattering_matrix(E, dev)
    k = list(np.where(res["modes"]["open"])[0]).index(n - 1)
    t = res["t_RL"][k, k] if w["direction"] == "L->R" else res["t_LR"][k, k]
    phi = ap.reference_phase_wrapped(E, dev.basis.E_perp[n - 1], dev.dx, dev.n_slices, ap.REF_DPS)
    assert phi == float(w["phase_ref_hp"])                 # bitwise: same code path and precision
    assert abs(ap.wrapped_difference(t, phi) - float(w["dphase_wrapped"])) <= 1e-15
