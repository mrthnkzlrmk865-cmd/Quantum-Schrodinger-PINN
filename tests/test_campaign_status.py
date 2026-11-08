"""
Fail-closed campaign status: gate registry, driver orchestration, independent verifier (docs/PHYSICAL_CAMPAIGN_PLAN.md section 8).

All campaign data here are SYNTHETIC tables (no simulation is run). The only real data used is the read-only regression against the
existing B3Y result directory, whose files are hashed before and after to prove they are unchanged.
"""
import ast
import csv
import hashlib
import json
import math
import os
import shutil

import pytest

import campaign_status as cs
import physical_campaign_driver as drv
import verify_campaign_outputs as ver

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
B3Y_DIR = os.path.join(HERE, "results", "rgf_b3y_validation")
B3Y_PROTOCOL = os.path.join(HERE, "docs", "rgf_b3y_protocol.json")
N_TEST = 100                                     # the protocol minimum N_min (any admissible N works)


# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def rows_of(path):
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        return list(r), r.fieldnames


def write_rows(path, rows, cols=None):
    cols = cols or list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def edit_csv(path, fn):
    rows, cols = rows_of(path)
    rows = fn(rows)
    cols = list(dict.fromkeys(list(cols) + [k for r in rows for k in r]))
    write_rows(path, rows, cols)


def tree_hashes(d):
    out = {}
    for root, _, files in os.walk(d):
        for fn in sorted(files):
            p = os.path.join(root, fn)
            with open(p, "rb") as f:
                out[os.path.relpath(p, d)] = hashlib.sha256(f.read()).hexdigest()
    return out


# ----------------------------------------------------------------------
# synthetic physical campaign: passing stage implementations
# ----------------------------------------------------------------------
DRAFT = cs.load_draft()


def passing_stages(N=N_TEST, mutate=None):
    D = DRAFT
    cells = [(w, z) for w in D["grid"]["W"] for z in D["grid"]["zeta"]]
    corners = [tuple(c) for c in D["grid"]["corner_cells"]]
    E = range(D["energies"]["n_energies"])
    seeds = list(range(1, N + 1))
    cid = cs.cell_id

    def control_rows(cid_, cells_, n_k):
        out = []
        cset = set(cells_)
        for o in D["controls"][cid_]["observables"]:
            for (w, z) in sorted(cset):
                out.append(dict(control_id=cid_, observable=o, unit_type="cell", unit_id=cid(w, z), n_k=n_k, d_mean=1e-6, se_d=1e-6, h=0.01))
            for c in D["contrasts"]:
                if all(tuple(x) in cset for x in c["cells"]):
                    out.append(dict(control_id=cid_, observable=o, unit_type="contrast", unit_id=c["id"], n_k=n_k, d_mean=-1e-6, se_d=1e-6, h=0.01))
        return out

    st = {
        "C6_clean_guide": lambda ctx: {"clean_guide": [dict(energy_index=e, max_abs_Tnn_minus_1=1e-14) for e in E]},
        "C5_b3y_reproduction": lambda ctx: {"b3y_reproduction": [dict(zeta=z, energy_index=e, abs_dT1=1e-11) for z in D["controls"]["C5"]["zeta"] for e in E]},
        "main_ensemble": lambda ctx: {"realizations": [dict(W=w, zeta=z, seed=s, energy_index=e, status="OK", unitarity=1e-13, reciprocity=1e-13, n_open=6,
                                                            n_open_analytic=6, lead_tail_max=1e-12, v_hash="sha256:x", tau=0.5, T1=0.4)
                                                       for (w, z) in cells for s in seeds for e in E]},
        "C0_reproducibility": lambda ctx: {"reproducibility": [dict(W=cells[i % 16][0], zeta=cells[i % 16][1], seed=1 + i, bitwise_equal=True)
                                                               for i in range(math.ceil(0.01 * 16 * N))]},
        "C2_solver_cross_check": lambda ctx: {"solver_cross_check": [dict(W=w, zeta=z, seed=s, energy_index=e, max_abs_dS=1e-12)
                                                                     for (w, z) in cells for s in seeds[:10] for e in E], "adjudication": []},
        "C1_dx_1024": lambda ctx: {"controls": control_rows("C1", cells, 20)},
        "C1b_dx_2048": lambda ctx: {"controls": control_rows("C1b", corners, 5)},
        "C3_basis": lambda ctx: {"controls": control_rows("C3", [(25.0, 0.025), (200.0, 0.025)], 10)},
        "C4_buffer": lambda ctx: {"controls": control_rows("C4", corners, 10)},
        "statistics": lambda ctx: {"tests": [dict(observable=o, contrast_id=c["id"], p_hat=0.3, p_holm_adj_from_phi=1.0, decision="NOT REJECTED")
                                             for o in D["statistics"]["observables"] for c in D["contrasts"]]},
        "C7_precision": lambda ctx: {"precision": [dict(observable="tau_bar", unit_type="cell", unit_id=cid(w, z), half_width=0.015) for (w, z) in cells] +
                                     [dict(observable="tau_bar", unit_type="contrast", unit_id=c["id"], half_width=0.018) for c in D["contrasts"]]},
    }
    if mutate:
        mutate(st)
    return st


def run_phys(tmp_path, mutate=None, N=N_TEST, name="run"):
    out = str(tmp_path / name)
    s = drv.run_campaign(DRAFT, out, N, passing_stages(N, mutate), log=lambda m: None)
    return out, s


def wrap(st, stage, table, fn):
    """Replace the rows that one stage returns for one table by fn(rows)."""
    inner = st[stage]

    def new(ctx):
        r = inner(ctx)
        r[table] = fn(list(r[table]))
        return r
    st[stage] = new


def reaggregate_phys(out, N=N_TEST):
    return cs.aggregate(cs.physical_gate_spec(DRAFT, N), out)["overall_status"]


# ----------------------------------------------------------------------
# synthetic B3Y-shaped outputs (passing), derived from the real B3Y protocol keys
# ----------------------------------------------------------------------
def make_b3y_like(d):
    P = json.load(open(B3Y_PROTOCOL))
    S = P["point_sets"]
    cfg = {p["id"]: p["config"] for p in P["points"]}
    i1, i2, i3 = 512, 1024, 2048
    os.makedirs(d, exist_ok=True)
    comp = []
    for p, i in [(p, i) for p in S["tier1"] for i in (i1, i2)] + [(p, i3) for p in S["tier2_dx2048"]]:
        r = dict(point=p, inv_dx=i, status="OK", V_lead_max=1e-7, S1_dS=1e-12, S2_dS=1e-12, S3_dS=1e-12, S4_S1_max=1e-12, S4_S2_max=1e-12,
                 R1_backward_error_M2=1e-16, R1_backward_error_M3=1e-16)
        for m in ("M1", "M2", "M3"):
            r[f"U_{m}"] = 1e-13
            r[f"Q_{m}"] = 1e-13
        for c in ("P4_lead_guard_ok", "N1_passed", "P1_rgf_device_bitwise", "P1_M2_samples_bitwise", "P1_M3_samples_bitwise", "P1_Eperp_bitwise_across_methods"):
            r[c] = True
        comp.append(r)
    write_rows(os.path.join(d, "per_energy_comparison.csv"), comp)
    t3 = [dict(point=p, inv_dx=i, X4_state_passed=True, X4_green_passed=True, X5_passed=True, X7_passed=True, gates_passed=True, X7_error_indicator=1e-13,
               X5_fl_vs_state_max_abs=1e-17, eps_M1=1e-12, eps_M2=1e-12, eps_M3=1e-12)
          for p, i in [(p, i1) for p in S["tier3_xp_512_validation"]] + [(p, i) for p in S["tier3_xp_1024_2048"] for i in (i2, i3)]]
    write_rows(os.path.join(d, "tier3_exact_lattice.csv"), t3)
    dec = [dict(point=p, inv_dx=i, method=m, total=1e-12, bulk_rounding=1e-12, lead_data=1e-13, solver=1e-13, identity_residual_three_term=0.0,
                e_ind_combined=1e-14, gate_fl_passed=True, gate_hybrid_passed=True, mirror_source_bitwise=True, mirror_extraction_bitwise=True,
                M2_recomputed_bitwise=True, M3_recomputed_bitwise=True)
           for p in S["tier3_xp_1024_2048"] for i in (i2, i3) for m in ("M2", "M3")]
    write_rows(os.path.join(d, "error_decomposition.csv"), dec)
    allpts = list(dict.fromkeys(S["tier3_xp_512_validation"] + S["tier3_xp_1024_2048"]))
    xg = [dict(gate="G0", case=p, inv_dx="", passed=True) for p in S["tier1"]]
    xg += [dict(gate="X1", case=p, inv_dx=i, passed=True) for p in allpts for i in (i1, i2, i3)]
    xg += [dict(gate=g, case=p, inv_dx=i, passed=True) for g in ("X2", "X3") for p in S["tier3_xp_1024_2048"] for i in (i1, i2, i3)]
    xg += [dict(gate="X6", case=p, inv_dx=i, passed=True) for p in S["tier3_xp_1024_2048"] for i in (i1, i2)]
    xg += [dict(gate="X0", case="platform", inv_dx="", passed=True), dict(gate="X8", case="all", inv_dx="", passed=True)]
    write_rows(os.path.join(d, "xp_gates.csv"), xg)
    write_rows(os.path.join(d, "buffer_control.csv"), [dict(point=p, K1_passed_dT1=True, S_within_geometry_passed=True, abs_dT1_M1=1e-11, abs_dT1_M3=1e-11)
                                                       for p in S["buffer_control"]])
    write_rows(os.path.join(d, "solver_spread.csv"), [dict(point=p, inv_dx=i, passed=True, T1_spread=1e-12) for p in S["tier2_dx2048"] for i in (i1, i2, i3)])
    write_rows(os.path.join(d, "refinement_convergence.csv"), [dict(point=p, method=m, observed_order=2.0) for p in S["tier2_dx2048"] for m in ("M1", "M2", "M3", "M4")])
    cf = lambda ids: sorted({cfg[i] for i in ids})
    pi = [dict(config=c, inv_dx=i, buffer_requested=b, P1_rgf_device_vs_direct_bitwise=True, P4_lead_guard_ok=True, V_lead_max=1e-7)
          for c, i, b in [(c, i1, 0.2) for c in cf(S["tier1"])] + [(c, i2, 0.2) for c in cf(S["tier1"])] + [(c, i3, 0.2) for c in cf(S["tier2_dx2048"])] +
          [(c, i2, 0.4) for c in cf(S["buffer_control"])]]
    pi += [dict(config=c, inv_dx=f"{i}_vs_{i1}", buffer_requested="", P2_potential_bitwise=True, coincident_x_bitwise=True) for c in cf(S["tier2_dx2048"]) for i in (i2, i3)]
    write_rows(os.path.join(d, "potential_identity.csv"), pi)
    summ = dict(overall_status="PASS", stopped_at=None,
                criteria=dict(S1_M1_vs_M2=dict(status="PASS"), XP_gates=dict(X2=dict(status="PASS"), X3=dict(status="PASS")),
                              O1_convergence_orders=dict(status="CLASSIFICATION ONLY")))
    json.dump(summ, open(os.path.join(d, "summary.json"), "w"))
    return d


def b3y_status(d):
    return cs.aggregate_b3y(d, B3Y_PROTOCOL)


def gate_state(res, gid):
    return [g["state"] for g in res["registry"] if g["gate_id"] == gid][0]


# ======================================================================
# baselines: synthetic outputs that satisfy every gate PASS, and the verifier agrees
# ======================================================================
def test_physical_baseline_passes_and_verifier_agrees(tmp_path):
    out, s = run_phys(tmp_path)
    assert s["overall_status"] == "PASS" and s["campaign_status"] == "COMPLETED" and s["stopped_at"] is None
    assert s["n_mandatory_gates"] == len(cs.physical_gate_spec(DRAFT, N_TEST)) == s["n_gates_pass"]
    assert all(g["state"] == "PASS" for g in s["registry"].values())
    v = ver.verify_physical(out, N_TEST)
    assert v["recomputed"] == "PASS" and v["agree"] and v["summary_consistent"], v["failures"][:3]
    assert s["row_counts"]["realizations"] == 16 * N_TEST * 5


def test_b3y_like_baseline_passes_and_verifier_agrees(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))
    r = b3y_status(d)
    assert r["overall_status"] == "PASS", r["problems"]
    v = ver.verify_b3y(d, B3Y_PROTOCOL)
    assert v["recomputed"] == "PASS" and v["agree"], v["failures"][:3]


# ======================================================================
# mutation 1: a failing nested gate
# ======================================================================
def test_m1_failing_nested_gate_in_summary_prevents_pass(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))
    s = json.load(open(os.path.join(d, "summary.json")))
    s["criteria"]["XP_gates"]["X3"]["status"] = "FAIL"                      # hidden one level down, overall still claims PASS
    json.dump(s, open(os.path.join(d, "summary.json"), "w"))
    r = b3y_status(d)
    assert r["overall_status"] != "PASS" and gate_state(r, "B3Y.summary.nested_statuses") == "FAIL"
    v = ver.verify_b3y(d, B3Y_PROTOCOL)
    assert not v["summary_consistent"] and not v["agree"]


def test_m1_failing_nested_row_level_gate_prevents_pass_and_verifier_flags_summary(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))
    edit_csv(os.path.join(d, "xp_gates.csv"), lambda rows: [dict(r, passed=False) if (r["gate"] == "X3" and r["inv_dx"] == "2048") else r for r in rows][:])
    r = b3y_status(d)
    assert r["overall_status"] == "FAIL" and gate_state(r, "B3Y.xp_gate.X3") == "FAIL"
    v = ver.verify_b3y(d, B3Y_PROTOCOL)
    assert v["recomputed"] == "NOT PASS" and not v["agree"]                  # raw data disagree with the summary's PASS


def test_m1_deeply_nested_and_list_statuses_are_flattened():
    flat = cs.flatten_statuses({"a": {"b": [{"c": {"status": "FAIL"}}, {"passed": False}]}, "status": "PASS"})
    assert ("/a/b[0]/c/status", "FAIL") in flat and ("/a/b[1]/passed", False) in flat and ("/status", "PASS") in flat


# ======================================================================
# mutation 2: a failed extended-precision or solver-cross-check row
# ======================================================================
def test_m2_failed_extended_precision_row_prevents_pass(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))
    edit_csv(os.path.join(d, "tier3_exact_lattice.csv"),
             lambda rows: [dict(r, X7_passed=False, gates_passed=False, X7_error_indicator=1e-9) if r["inv_dx"] == "2048" else r for r in rows])
    r = b3y_status(d)
    assert r["overall_status"] == "FAIL" and gate_state(r, "B3Y.tier3.X4_X5_X7_eps") == "FAIL"
    assert ver.verify_b3y(d, B3Y_PROTOCOL)["recomputed"] == "NOT PASS"


def test_m2_failed_solver_cross_check_row_prevents_pass(tmp_path):
    out, s = run_phys(tmp_path, mutate=lambda st: wrap(st, "C2_solver_cross_check", "solver_cross_check",
                                                       lambda rows: [dict(rows[0], max_abs_dS=2e-10)] + rows[1:]))
    assert s["overall_status"] == "FAIL" and s["registry"]["PHYS.C2.solver_cross_check"]["state"] == "FAIL"
    rows, _ = rows_of(os.path.join(out, "solver_cross_check.csv"))
    assert len(rows) == 16 * 10 * 5 and any(float(r["max_abs_dS"]) == 2e-10 for r in rows)          # the failing row is kept
    v = ver.verify_physical(out, N_TEST)
    assert v["recomputed"] == "NOT PASS" and v["agree"]


# ======================================================================
# mutation 3: a NOT EVALUABLE decomposition row
# ======================================================================
def test_m3_not_evaluable_decomposition_row_prevents_pass(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))

    def ne(rows):
        r0 = rows[0]
        return [dict(point=r0["point"], inv_dx=r0["inv_dx"], method=r0["method"], status="NOT EVALUABLE")] + rows[1:]
    edit_csv(os.path.join(d, "error_decomposition.csv"), ne)
    r = b3y_status(d)
    assert gate_state(r, "B3Y.decomposition.evaluable") == "NOT EVALUABLE" and r["overall_status"] == "NOT EVALUABLE"
    assert ver.verify_b3y(d, B3Y_PROTOCOL)["recomputed"] == "NOT PASS"


def test_m3_extra_not_evaluable_row_beyond_protocol_keys_is_not_ignored(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))
    edit_csv(os.path.join(d, "error_decomposition.csv"), lambda rows: rows + [dict(point="W200_z0.2@E88.51", inv_dx=512, method="M2", status="NOT EVALUABLE")])
    assert b3y_status(d)["overall_status"] == "NOT EVALUABLE"


def test_m3_not_evaluable_adjudication_row_in_physical_campaign(tmp_path):
    def m(st):
        inner = st["C2_solver_cross_check"]
        st["C2_solver_cross_check"] = lambda ctx: dict(inner(ctx), adjudication=[dict(W=25.0, zeta=0.025, seed=1, energy_index=0, status="NOT EVALUABLE")])
    out, s = run_phys(tmp_path, mutate=m)
    assert s["registry"]["PHYS.C2.adjudication_evaluable"]["state"] == "NOT EVALUABLE" and s["overall_status"] == "NOT EVALUABLE"
    assert ver.verify_physical(out, N_TEST)["recomputed"] == "NOT PASS"


# ======================================================================
# mutation 4: a missing expected row
# ======================================================================
def test_m4_missing_realization_row_prevents_pass(tmp_path):
    out, s = run_phys(tmp_path, mutate=lambda st: wrap(st, "main_ensemble", "realizations", lambda rows: rows[:-1]))
    assert s["overall_status"] == "FAIL" and "missing" in s["registry"]["PHYS.C0.realizations"]["reason"]
    v = ver.verify_physical(out, N_TEST)
    assert v["recomputed"] == "NOT PASS" and v["agree"]


def test_m4_missing_b3y_comparison_row_prevents_pass(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))
    edit_csv(os.path.join(d, "per_energy_comparison.csv"), lambda rows: [r for r in rows if not (r["point"] == "W200_z0.2@E88.51" and r["inv_dx"] == "2048")])
    r = b3y_status(d)
    assert r["overall_status"] == "FAIL" and "missing" in [g for g in r["registry"] if g["gate_id"].startswith("B3Y.comparison")][0]["reason"]
    assert ver.verify_b3y(d, B3Y_PROTOCOL)["recomputed"] == "NOT PASS"


def test_m4_missing_table_and_duplicate_row(tmp_path):
    out, _ = run_phys(tmp_path)
    os.remove(os.path.join(out, "clean_guide.csv"))
    res = cs.aggregate(cs.physical_gate_spec(DRAFT, N_TEST), out)
    assert res["overall_status"] != "PASS" and [g for g in res["registry"] if g["gate_id"] == "PHYS.C6.clean_guide"][0]["state"] == "NOT RUN"
    out2, s2 = run_phys(tmp_path, mutate=lambda st: wrap(st, "C5_b3y_reproduction", "b3y_reproduction", lambda rows: rows + rows[:1]), name="dup")
    assert s2["overall_status"] == "FAIL" and "duplicated" in s2["registry"]["PHYS.C5.b3y_reproduction"]["reason"]


def test_m4_wrong_count_of_random_reproducibility_subset(tmp_path):
    out, s = run_phys(tmp_path, mutate=lambda st: wrap(st, "C0_reproducibility", "reproducibility", lambda rows: rows[:-1]))
    assert s["overall_status"] == "FAIL" and s["registry"]["PHYS.C0.reproducibility"]["state"] == "FAIL"


# ======================================================================
# mutation 5: a missing or unknown status
# ======================================================================
def test_m5_unknown_row_status_prevents_pass(tmp_path):
    out, s = run_phys(tmp_path, mutate=lambda st: wrap(st, "main_ensemble", "realizations", lambda rows: [dict(rows[0], status="MAYBE")] + rows[1:]))
    assert s["overall_status"] == "FAIL" and "MAYBE" in s["registry"]["PHYS.C0.realizations"]["reason"] + str(
        [g for g in cs.aggregate(cs.physical_gate_spec(DRAFT, N_TEST), out)["registry"] if g["gate_id"] == "PHYS.C0.realizations"][0]["first_failures"])
    assert ver.verify_physical(out, N_TEST)["recomputed"] == "NOT PASS"


def test_m5_missing_row_status_and_missing_stage_status_prevent_pass(tmp_path):
    out, s = run_phys(tmp_path, mutate=lambda st: wrap(st, "main_ensemble", "realizations", lambda rows: [dict(rows[0], status="")] + rows[1:]))
    assert s["overall_status"] == "FAIL"
    out2, _ = run_phys(tmp_path, name="stage_blank")
    edit_csv(os.path.join(out2, "stages.csv"), lambda rows: [dict(r, status="") if r["stage"] == "statistics" else r for r in rows])
    assert reaggregate_phys(out2) == "FAIL"
    assert ver.verify_physical(out2, N_TEST)["recomputed"] == "NOT PASS"


def test_m5_unknown_or_missing_registry_state_and_unknown_decision():
    good = [dict(gate_id="a", state="PASS"), dict(gate_id="b", state="PASS")]
    assert cs.combine(good, ["a", "b"])["overall_status"] == "PASS"
    assert cs.combine([dict(gate_id="a", state="PASS"), dict(gate_id="b", state="PASSED")], ["a", "b"])["overall_status"] == "FAIL"
    assert cs.combine([dict(gate_id="a", state="PASS"), dict(gate_id="b", state=None)], ["a", "b"])["overall_status"] == "FAIL"
    assert cs.combine([dict(gate_id="a", state="PASS")], ["a", "b"])["overall_status"] == "FAIL"                 # missing entry
    assert cs.combine(good + [dict(gate_id="a", state="PASS")], ["a", "b"])["overall_status"] == "FAIL"           # duplicated entry
    assert cs.combine([], [])["overall_status"] == "FAIL"                                                       # nothing declared
    assert cs.combine([dict(gate_id="a", state="PASS"), dict(gate_id="b", state="NOT RUN")], ["a", "b"])["overall_status"] == "INCOMPLETE"


def test_m5_unknown_statistical_decision_prevents_pass(tmp_path):
    out, s = run_phys(tmp_path, mutate=lambda st: wrap(st, "statistics", "tests", lambda rows: [dict(rows[0], decision="SIGNIFICANT")] + rows[1:]))
    assert s["overall_status"] == "FAIL" and s["registry"]["PHYS.statistics.tests_complete"]["state"] == "FAIL"


def test_m5_unknown_nested_status_string_in_summary(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))
    s = json.load(open(os.path.join(d, "summary.json")))
    s["criteria"]["S1_M1_vs_M2"]["status"] = "pass"                          # wrong case: not the vocabulary
    json.dump(s, open(os.path.join(d, "summary.json"), "w"))
    assert b3y_status(d)["overall_status"] == "FAIL"


# ======================================================================
# mutation 6: a stage that stopped early
# ======================================================================
def test_m6_exception_in_a_stage_stops_the_run_and_prevents_pass(tmp_path):
    def m(st):
        def boom(ctx):
            raise RuntimeError("solver crashed")
        st["C1_dx_1024"] = boom
    out, s = run_phys(tmp_path, mutate=m)
    assert s["overall_status"] != "PASS" and s["stopped_at"] == "C1_dx_1024" and s["campaign_status"].startswith("STOPPED")
    stages, _ = rows_of(os.path.join(out, "stages.csv"))
    st = {r["stage"]: r["status"] for r in stages}
    assert st["C1_dx_1024"] == "FAILED" and all(st[k] == "NOT RUN" for k in ("C1b_dx_2048", "C3_basis", "C4_buffer", "statistics", "C7_precision"))
    assert len(rows_of(os.path.join(out, "realizations.csv"))[0]) == 16 * N_TEST * 5                    # completed stages' rows preserved
    assert s["registry"]["PHYS.stages.completed"]["state"] == "FAIL"
    v = ver.verify_physical(out, N_TEST)
    assert v["recomputed"] == "NOT PASS" and v["agree"]


def test_m6_missing_stage_implementation_is_not_run_and_never_pass(tmp_path):
    out, s = run_phys(tmp_path, mutate=lambda st: st.pop("statistics"))
    assert s["overall_status"] != "PASS" and s["stopped_at"] == "statistics"
    assert s["registry"]["PHYS.statistics.tests_complete"]["state"] in ("NOT RUN", "FAIL")


def test_m6_undeclared_table_returned_by_a_stage_stops_the_run(tmp_path):
    def m(st):
        inner = st["C6_clean_guide"]
        st["C6_clean_guide"] = lambda ctx: dict(inner(ctx), surprise=[{"x": 1}])
    out, s = run_phys(tmp_path, mutate=m)
    assert s["stopped_at"] == "C6_clean_guide" and s["overall_status"] != "PASS"


def test_m6_upstream_summary_stopped_early(tmp_path):
    d = make_b3y_like(str(tmp_path / "b3y"))
    s = json.load(open(os.path.join(d, "summary.json")))
    s["stopped_at"] = dict(stage="X8/V1", reason="test")
    json.dump(s, open(os.path.join(d, "summary.json"), "w"))
    assert b3y_status(d)["overall_status"] == "FAIL"
    assert not ver.verify_b3y(d, B3Y_PROTOCOL)["agree"]


# ======================================================================
# mutation 7: an under-powered control
# ======================================================================
def test_m7_underpowered_control_prevents_pass_even_if_tost_interval_is_inside_the_margin(tmp_path):
    from scipy.stats import t as tdist
    h = 0.01
    eps = 0.25 * h
    se = 0.6 * eps / float(tdist.ppf(0.95, 19))                              # t*SE = 0.6 eps > eps/2 (under-powered) but |d| + t*SE < eps
    def m(st):
        wrap(st, "C1_dx_1024", "controls", lambda rows: [dict(rows[0], d_mean=0.0, se_d=se, h=h)] + rows[1:])
    out, s = run_phys(tmp_path, mutate=m)
    assert s["overall_status"] == "FAIL" and s["registry"]["PHYS.C1.equivalence"]["state"] == "FAIL"
    reg = {g["gate_id"]: g for g in cs.aggregate(cs.physical_gate_spec(DRAFT, N_TEST), out)["registry"]}
    assert "UNDER-POWERED" in reg["PHYS.C1.equivalence"]["first_failures"]
    v = ver.verify_physical(out, N_TEST)
    assert v["recomputed"] == "NOT PASS" and any("under-powered" in f for f in v["failures"])


def test_m7_tost_rule_boundaries():
    ev = cs.tost_evaluator(0.25, 0.05, 0.5, 100)
    base = dict(n_k=20, h=0.01)
    assert ev(dict(base, d_mean=0.0, se_d=1e-6)) is None                                          # accepted, powered
    assert "TOST not accepted" in ev(dict(base, d_mean=0.0024, se_d=1e-4))                        # interval crosses +eps (powered)
    assert "UNDER-POWERED" in ev(dict(base, d_mean=0.0, se_d=0.002))
    assert ev(dict(base, d_mean=0.0, se_d=1e-6, n_k=101)) is not None                             # beyond the extension cap
    assert ev(dict(base, d_mean="nan", se_d=1e-6)) is not None and ev(dict(base, d_mean=0.0, se_d=1e-6, h=0)) is not None


# ======================================================================
# protocol derivation, row preservation, independence, scope
# ======================================================================
def test_registry_is_protocol_derived_not_run_derived(tmp_path):
    specs = cs.physical_gate_spec(DRAFT, N_TEST)
    ids = [s.gate_id for s in specs]
    assert len(ids) == len(set(ids)) and {"PHYS.stages.completed", "PHYS.C0.realizations", "PHYS.C2.solver_cross_check", "PHYS.C1.equivalence",
                                          "PHYS.C1b.equivalence", "PHYS.C3.equivalence", "PHYS.C4.equivalence", "PHYS.C5.b3y_reproduction",
                                          "PHYS.C6.clean_guide", "PHYS.C7.precision", "PHYS.statistics.tests_complete"} <= set(ids)
    exp = {s.gate_id: (len(s.expected_keys) if s.expected_keys is not None else s.expected_count) for s in specs}
    assert exp["PHYS.C0.realizations"] == 16 * N_TEST * 5 and exp["PHYS.C2.solver_cross_check"] == 16 * 10 * 5
    assert exp["PHYS.C1.equivalence"] == 2 * (16 + 34) and exp["PHYS.C7.precision"] == 16 + 34 and exp["PHYS.statistics.tests_complete"] == 68
    assert exp["PHYS.C1b.equivalence"] == 2 * 4 and exp["PHYS.C3.equivalence"] == 2 * 2 and exp["PHYS.C0.reproducibility"] == math.ceil(0.01 * 16 * N_TEST)
    for bad in (None, 99, 1001, 150.5):
        with pytest.raises(ValueError):
            cs.physical_gate_spec(DRAFT, bad)
    out = str(tmp_path / "empty")
    os.makedirs(out)
    res = cs.aggregate(specs, out)                                           # nothing ran: every gate present in the registry, none PASS
    assert res["overall_status"] != "PASS" and len(res["registry"]) == len(specs) and all(g["state"] == "NOT RUN" for g in res["registry"])


def test_invalid_n_seeds_in_driver_is_fail_closed(tmp_path):
    s = drv.run_campaign(DRAFT, str(tmp_path / "badN"), 50, passing_stages(50), log=lambda m: None)
    assert s["overall_status"] == "FAIL" and "could not be built" in s["problems"][0]


def test_failed_rows_are_written_not_dropped(tmp_path):
    out, s = run_phys(tmp_path, mutate=lambda st: wrap(st, "main_ensemble", "realizations",
                                                       lambda rows: [dict(rows[0], status="NOT EVALUABLE", unitarity="", tau="")] + [dict(rows[1], unitarity=1e-3)] + rows[2:]))
    rows, _ = rows_of(os.path.join(out, "realizations.csv"))
    assert len(rows) == 16 * N_TEST * 5 and rows[0]["status"] == "NOT EVALUABLE" and float(rows[1]["unitarity"]) == 1e-3
    assert s["overall_status"] == "FAIL"


def test_driver_refuses_historical_result_directories():
    for d in ("rgf_b3y_validation", "rgf_b3x_validation", "rgf_b1_validation"):
        with pytest.raises(SystemExit):
            drv.check_output_dir(os.path.join(HERE, "results", d))


def _imports(path):
    tree = ast.parse(open(path).read())
    names = set()
    for n in ast.walk(tree):
        if isinstance(n, ast.Import):
            names |= {a.name.split(".")[0] for a in n.names}
        elif isinstance(n, ast.ImportFrom) and n.module:
            names.add(n.module.split(".")[0])
    return names


def test_independence_of_modules():
    v = _imports(os.path.join(HERE, "verify_campaign_outputs.py"))
    assert "campaign_status" not in v and "physical_campaign_driver" not in v and "validate_rgf_b3y" not in v
    for f in ("campaign_status.py", "physical_campaign_driver.py"):
        im = _imports(os.path.join(HERE, f))
        assert "validate_rgf_b3y" not in im and "stationary" not in im and "sparse_exact_lattice" not in im     # no physics, no old validator


def test_draft_spec_matches_the_approved_plan_text():
    plan = open(os.path.join(HERE, "docs", "PHYSICAL_CAMPAIGN_PLAN.md")).read()
    D = DRAFT
    assert "NOT frozen" in D["status"] and D["source_plan_commit"] == "c4d1dc1" and D["seeds"]["N"] is None
    assert D["grid"]["W"] == [25.0, 50.0, 100.0, 200.0] and D["grid"]["zeta"] == [0.025, 0.05, 0.1, 0.2]
    for frag in ("W | 25, 50, 100, 200", "ζ | 0.025, 0.05, 0.10, 0.20", "W·ζ = 5", "N_max = 1000", "δ = 0.02", "first 20 seeds", "first 5 seeds",
                 "first 10 seeds", "410/512 vs 205/512", "**68 in total**", "unitarity ≤ 10⁻¹⁰", "reciprocity ≤ 10⁻⁸", "lead-tail ≤ 10⁻⁶",
                 "\\|ΔT1\\| ≤ 10⁻⁶", "T_nn = 1 within 10⁻¹²", "1 % random subset", "up to 100"):
        assert frag in plan, frag
    L = D["limits"]
    assert (L["unitarity"], L["reciprocity"], L["lead_tail_guard"], L["delta"], L["equivalence_fraction_of_h"], L["tost_alpha"],
            L["underpowered_fraction_of_eps"], L["control_seed_extension_max"]) == (1e-10, 1e-8, 1e-6, 0.02, 0.25, 0.05, 0.5, 100)
    assert len(D["contrasts"]) == 34 and D["statistics"]["n_tests"] == 68
    fam = [c["family"] for c in D["contrasts"]]
    assert (fam.count("C-W"), fam.count("C-zeta"), fam.count("C-int"), fam.count("C-iso")) == (12, 12, 9, 1)
    iso = [c for c in D["contrasts"] if c["family"] == "C-iso"][0]
    assert all(abs(w * z - 5.0) < 1e-12 for w, z in iso["cells"])


# ======================================================================
# B3Y regression (read-only): the new aggregator and the verifier reproduce PASS on the existing results without modifying them
# ======================================================================
def test_b3y_regression_reproduces_pass_without_modifying_any_file():
    before = tree_hashes(B3Y_DIR)
    proto_before = hashlib.sha256(open(B3Y_PROTOCOL, "rb").read()).hexdigest()
    r = cs.aggregate_b3y(B3Y_DIR, B3Y_PROTOCOL)
    v = ver.verify_b3y(B3Y_DIR, B3Y_PROTOCOL)
    assert r["overall_status"] == "PASS", r["problems"]
    assert all(g["state"] == "PASS" for g in r["registry"]) and len(r["registry"]) == 16
    assert v["recomputed"] == "PASS" and v["summary_overall"] == "PASS" and v["agree"] and v["summary_consistent"]
    assert tree_hashes(B3Y_DIR) == before and hashlib.sha256(open(B3Y_PROTOCOL, "rb").read()).hexdigest() == proto_before
    assert "PASS" == json.load(open(os.path.join(B3Y_DIR, "summary.json")))["overall_status"]


def test_b3y_regression_copy_with_a_single_nested_failure_is_not_pass(tmp_path):
    """The three B3Y-driver gaps, demonstrated on a COPY of the real results (the originals are untouched)."""
    d = str(tmp_path / "copy")
    shutil.copytree(B3Y_DIR, d)
    s = json.load(open(os.path.join(d, "summary.json")))
    s["criteria"]["XP_gates"]["X6"]["status"] = "FAIL"
    json.dump(s, open(os.path.join(d, "summary.json"), "w"))
    assert cs.aggregate_b3y(d, B3Y_PROTOCOL)["overall_status"] == "FAIL"                              # gap (i)
    d2 = str(tmp_path / "copy2")
    shutil.copytree(B3Y_DIR, d2)
    edit_csv(os.path.join(d2, "tier3_exact_lattice.csv"), lambda rows: [dict(r, gates_passed=False) if r["inv_dx"] == "2048" else r for r in rows])
    assert cs.aggregate_b3y(d2, B3Y_PROTOCOL)["overall_status"] == "FAIL"                             # gap (ii)
    d3 = str(tmp_path / "copy3")
    shutil.copytree(B3Y_DIR, d3)
    edit_csv(os.path.join(d3, "error_decomposition.csv"), lambda rows: [dict(rows[0], status="NOT EVALUABLE", bulk_rounding="")] + rows[1:])
    assert cs.aggregate_b3y(d3, B3Y_PROTOCOL)["overall_status"] == "NOT EVALUABLE"                    # gap (iii)
