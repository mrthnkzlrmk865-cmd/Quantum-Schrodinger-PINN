r"""
Independent verifier of campaign outputs (docs/PHYSICAL_CAMPAIGN_PLAN.md section 8, item 4).

Recomputes a binary verdict (PASS / NOT PASS) from the raw CSV files and the protocol, and checks it against summary.json.
Deliberately independent of campaign_status.py: it does not import it, uses its own CSV loading, its own derivation of the expected rows
from the protocol, its own row predicates (plain functions, not the declarative checks of the aggregator), its own nested-summary walk,
and its own TOST computation. Agreement of the two is the cross-check; disagreement is a FAIL.

  python verify_campaign_outputs.py --kind b3y --out results/rgf_b3y_validation             (read-only)
  python verify_campaign_outputs.py --kind physical --out DIR --n-seeds N
Exit status 0 iff the recomputed verdict agrees with summary.json and summary.json is internally consistent.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def _rows(out, name):
    p = os.path.join(out, name)
    if not os.path.isfile(p):
        return None
    with open(p, newline="") as f:
        return [dict(r) for r in csv.DictReader(f)]


def _f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    return v


def _k(*vals):
    out = []
    for v in vals:
        try:
            out.append(float(v))
        except (TypeError, ValueError):
            out.append(str(v).strip())
    return tuple(out)


def _le(x, lim):
    v = _f(x)
    return math.isfinite(v) and v <= lim


def _true(x):
    return str(x).strip() == "True"


def _complete(rows, keyfun, expected, failures, label, allow_dup=False):
    """Every expected key exactly once (or at least once if allow_dup)."""
    if rows is None:
        failures.append(f"{label}: table missing")
        return []
    got = {}
    for r in rows:
        got.setdefault(keyfun(r), []).append(r)
    for k in expected:
        n = len(got.get(k, []))
        if n == 0 or (n > 1 and not allow_dup):
            failures.append(f"{label}: key {k} present {n} times")
            return rows
    return rows


def _walk(o, path, acc):
    if isinstance(o, dict):
        for k, v in o.items():
            if k in ("status", "state", "passed") and not isinstance(v, (dict, list)):
                acc.append((f"{path}/{k}", v))
            else:
                _walk(v, f"{path}/{k}", acc)
    elif isinstance(o, list):
        for i, v in enumerate(o):
            _walk(v, f"{path}[{i}]", acc)
    return acc


# ----------------------------------------------------------------------
def verify_b3y(out, protocol_path=None):
    with open(protocol_path or os.path.join(HERE, "docs", "rgf_b3y_protocol.json")) as f:
        P = json.load(f)
    L, X, S = P["limits"], P["xp_gates"], P["point_sets"]
    cfg_of = {p["id"]: p["config"] for p in P["points"]}
    i1, i2, i3 = (int(round(1 / d)) for d in P["discretization"]["dx_levels"])
    fails = []
    # comparison
    rows = _complete(_rows(out, "per_energy_comparison.csv"), lambda r: _k(r["point"], r["inv_dx"]),
                     [_k(p, i) for p in S["tier1"] for i in (i1, i2)] + [_k(p, i3) for p in S["tier2_dx2048"]], fails, "comparison")
    for r in rows:
        ok = (r.get("status") == "OK" and all(_le(r.get(c), L["S_abs"]) for c in ("S1_dS", "S2_dS", "S3_dS"))
              and all(_le(r.get(c), L["TR_offdiag_abs"]) for c in ("S4_S1_max", "S4_S2_max"))
              and all(_le(r.get(f"U_{m}"), L["unitarity"]) and _le(r.get(f"Q_{m}"), L["reciprocity"]) for m in ("M1", "M2", "M3"))
              and all(_le(r.get(f"R1_backward_error_{m}"), L["sparse_backward_error"]) for m in ("M2", "M3"))
              and _le(r.get("V_lead_max"), L["lead_tail_guard"])
              and all(_true(r.get(c)) for c in ("P4_lead_guard_ok", "N1_passed", "P1_rgf_device_bitwise", "P1_M2_samples_bitwise",
                                                "P1_M3_samples_bitwise", "P1_Eperp_bitwise_across_methods")))
        if not ok:
            fails.append(f"comparison row {r.get('point')} {r.get('inv_dx')} fails")
    # tier 3
    rows = _complete(_rows(out, "tier3_exact_lattice.csv"), lambda r: _k(r["point"], r["inv_dx"]),
                     [_k(p, i1) for p in S["tier3_xp_512_validation"]] + [_k(p, i) for p in S["tier3_xp_1024_2048"] for i in (i2, i3)], fails, "tier3")
    for r in rows:
        if not (all(_true(r.get(c)) for c in ("X4_state_passed", "X4_green_passed", "X5_passed", "X7_passed", "gates_passed"))
                and _le(r.get("X7_error_indicator"), X["error_indicator_max"]) and _le(r.get("X5_fl_vs_state_max_abs"), X["fisher_lee_vs_wavefunction_forms_abs"])
                and all(_le(r.get(f"eps_{m}"), L["S_abs"]) for m in ("M1", "M2", "M3"))):
            fails.append(f"tier3 row {r.get('point')} {r.get('inv_dx')} fails")
    # decomposition: protocol-required rows present; every row (required or extra) evaluable
    rows = _complete(_rows(out, "error_decomposition.csv"), lambda r: _k(r["point"], r["inv_dx"], r["method"]),
                     [_k(p, i, m) for p in S["tier3_xp_1024_2048"] for i in (i2, i3) for m in ("M2", "M3")], fails, "decomposition")
    for r in rows:
        st = (r.get("status") or "").strip()
        terms = [_f(r.get(c)) for c in ("total", "bulk_rounding", "lead_data", "solver", "identity_residual_three_term", "e_ind_combined")]
        if st not in ("", "OK") or not all(math.isfinite(v) for v in terms) or not all(
                _true(r.get(c)) for c in ("gate_fl_passed", "gate_hybrid_passed", "mirror_source_bitwise", "mirror_extraction_bitwise",
                                          "M2_recomputed_bitwise", "M3_recomputed_bitwise")):
            fails.append(f"decomposition row {r.get('point')} {r.get('inv_dx')} {r.get('method')} not evaluable/failing ({st or 'no status'})")
    # xp gates
    xr = _rows(out, "xp_gates.csv")
    if xr is None:
        fails.append("xp_gates: table missing")
        xr = []
    allpts = list(dict.fromkeys(S["tier3_xp_512_validation"] + S["tier3_xp_1024_2048"]))
    exp = {"G0": [_k(p, "") for p in S["tier1"]], "X1": [_k(p, i) for p in allpts for i in (i1, i2, i3)],
           "X2": [_k(p, i) for p in S["tier3_xp_1024_2048"] for i in (i1, i2, i3)], "X3": [_k(p, i) for p in S["tier3_xp_1024_2048"] for i in (i1, i2, i3)],
           "X6": [_k(p, i) for p in S["tier3_xp_1024_2048"] for i in (i1, i2)]}
    for g, keys in exp.items():
        _complete([r for r in xr if r.get("gate") == g], lambda r: _k(r["case"], r["inv_dx"]), keys, fails, f"xp {g}")
    for g in ("X0", "X8"):
        if sum(r.get("gate") == g for r in xr) != 1:
            fails.append(f"xp {g}: not exactly one row")
    for r in xr:
        if not _true(r.get("passed")):
            fails.append(f"xp {r.get('gate')} {r.get('case')} {r.get('inv_dx')} not passed")
    # K1, H1, O1
    for r in _complete(_rows(out, "buffer_control.csv"), lambda r: _k(r["point"]), [_k(p) for p in S["buffer_control"]], fails, "K1"):
        if not (_true(r.get("K1_passed_dT1")) and _true(r.get("S_within_geometry_passed")) and _le(r.get("abs_dT1_M1"), L["dT1_buffer"])
                and _le(r.get("abs_dT1_M3"), L["dT1_buffer"])):
            fails.append(f"K1 {r.get('point')} fails")
    for r in _complete(_rows(out, "solver_spread.csv"), lambda r: _k(r["point"], r["inv_dx"]),
                       [_k(p, i) for p in S["tier2_dx2048"] for i in (i1, i2, i3)], fails, "H1"):
        if not (_true(r.get("passed")) and _le(r.get("T1_spread"), L["T1_method_spread_for_solver_independence"])):
            fails.append(f"H1 {r.get('point')} {r.get('inv_dx')} fails")
    for r in _complete(_rows(out, "refinement_convergence.csv"), lambda r: _k(r["point"], r["method"]),
                       [_k(p, m) for p in S["tier2_dx2048"] for m in ("M1", "M2", "M3", "M4")], fails, "O1"):
        if not math.isfinite(_f(r.get("observed_order"))):
            fails.append(f"O1 {r.get('point')} {r.get('method')} order missing")
    # P1/P4 devices and P2
    pi = _rows(out, "potential_identity.csv")
    cfgs = lambda ids: sorted({cfg_of[i] for i in ids})
    dev = [_k(c, i1, 0.2) for c in cfgs(S["tier1"])] + [_k(c, i2, 0.2) for c in cfgs(S["tier1"])] + \
        [_k(c, i3, 0.2) for c in cfgs(S["tier2_dx2048"])] + [_k(c, i2, 0.4) for c in cfgs(S["buffer_control"])]
    devrows = [r for r in (pi or []) if "_vs_" not in str(r.get("inv_dx"))]
    _complete(devrows if pi is not None else None, lambda r: _k(r["config"], r["inv_dx"], r["buffer_requested"]), dev, fails, "P1", allow_dup=True)
    for r in devrows:
        if not (_true(r.get("P1_rgf_device_vs_direct_bitwise")) and _true(r.get("P4_lead_guard_ok")) and _le(r.get("V_lead_max"), L["lead_tail_guard"])):
            fails.append(f"P1/P4 {r.get('config')} {r.get('inv_dx')} fails")
    p2 = [r for r in (pi or []) if "_vs_" in str(r.get("inv_dx"))]
    _complete(p2 if pi is not None else None, lambda r: _k(r["config"], r["inv_dx"]),
              [_k(c, f"{i}_vs_{i1}") for c in cfgs(S["tier2_dx2048"]) for i in (i2, i3)], fails, "P2")
    for r in p2:
        if not (_true(r.get("P2_potential_bitwise")) and _true(r.get("coincident_x_bitwise"))):
            fails.append(f"P2 {r.get('config')} {r.get('inv_dx')} fails")
    return _finish(out, fails, nested_subtree="criteria", informational=("CLASSIFICATION ONLY",))


# ----------------------------------------------------------------------
def verify_physical(out, n_seeds, draft_path=None):
    from scipy.stats import t as tdist
    with open(draft_path or os.path.join(HERE, "docs", "physical_campaign_gates_draft.json")) as f:
        D = json.load(f)
    fails = []
    G, L, C, sd = D["grid"], D["limits"], D["controls"], D["seeds"]
    if not (isinstance(n_seeds, int) and sd["N_min"] <= n_seeds <= sd["N_max"]):
        return dict(recomputed="NOT PASS", summary_overall=None, agree=False, failures=[f"N = {n_seeds!r} not admissible"], summary_consistent=False)
    allcells = [(w, z) for w in G["W"] for z in G["zeta"]]
    corners = [tuple(c) for c in G["corner_cells"]]
    sel = lambda spec: allcells if spec == "all" else corners if spec == "corners" else [tuple(c) for c in spec]
    seeds = list(range(sd["confirmatory_first"], sd["confirmatory_first"] + n_seeds))
    E = list(range(D["energies"]["n_energies"]))
    st = _rows(out, "stages.csv")
    _complete(st, lambda r: _k(r["stage"]), [_k(s) for s in D["stages"]], fails, "stages")
    for r in st or []:
        if r.get("status") != "COMPLETED":
            fails.append(f"stage {r.get('stage')}: {r.get('status')!r}")
    if st is not None and len(st) != len(D["stages"]):
        fails.append("stages: unexpected row count")
    for r in _complete(_rows(out, "realizations.csv"), lambda r: _k(r["W"], r["zeta"], r["seed"], r["energy_index"]),
                       [_k(w, z, s, e) for (w, z) in allcells for s in seeds for e in E], fails, "realizations"):
        if not (r.get("status") == "OK" and _le(r.get("unitarity"), L["unitarity"]) and _le(r.get("reciprocity"), L["reciprocity"])
                and _le(r.get("lead_tail_max"), L["lead_tail_guard"]) and _f(r.get("n_open")) == _f(r.get("n_open_analytic"))
                and str(r.get("v_hash", "")).strip() and math.isfinite(_f(r.get("tau"))) and math.isfinite(_f(r.get("T1")))):
            fails.append(f"realization {r.get('W')} {r.get('zeta')} {r.get('seed')} {r.get('energy_index')} fails ({r.get('status')})")
    rp = _rows(out, "reproducibility.csv")
    if rp is None or len(rp) != math.ceil(L["reproducibility_subset_fraction"] * len(allcells) * n_seeds) or not all(_true(r.get("bitwise_equal")) for r in rp):
        fails.append("reproducibility subset incomplete or not bitwise")
    c2 = C["C2"]
    for r in _complete(_rows(out, "solver_cross_check.csv"), lambda r: _k(r["W"], r["zeta"], r["seed"], r["energy_index"]),
                       [_k(w, z, s, e) for (w, z) in sel(c2["cells"]) for s in seeds[:c2["n_seeds"]] for e in E], fails, "C2"):
        if not _le(r.get("max_abs_dS"), L["solver_cross_check_S_abs"]):
            fails.append(f"C2 {r.get('W')} {r.get('zeta')} {r.get('seed')} {r.get('energy_index')} exceeds")
    adj = _rows(out, "adjudication.csv")
    if adj is None:
        fails.append("adjudication table missing")
    elif any(r.get("status") != "ADJUDICATED" for r in adj):
        fails.append("adjudication row not ADJUDICATED")
    for r in _complete(_rows(out, "b3y_reproduction.csv"), lambda r: _k(r["zeta"], r["energy_index"]),
                       [_k(z, e) for z in C["C5"]["zeta"] for e in E], fails, "C5"):
        if not _le(r.get("abs_dT1"), L["b3y_reproduction_abs_dT1"]):
            fails.append("C5 exceeds")
    for r in _complete(_rows(out, "clean_guide.csv"), lambda r: _k(r["energy_index"]), [_k(e) for e in E], fails, "C6"):
        if not _le(r.get("max_abs_Tnn_minus_1"), L["clean_guide_Tnn_abs"]):
            fails.append("C6 exceeds")
    cid_of = lambda w, z: f"W{float(w)!r}_z{float(z)!r}"
    ctr = _rows(out, "controls.csv")
    for cid in ("C1", "C1b", "C3", "C4"):
        cells = set(sel(C[cid]["cells"]))
        keys = []
        for o in C[cid]["observables"]:
            keys += [_k(cid, o, "cell", cid_of(w, z)) for (w, z) in sorted(cells)]
            keys += [_k(cid, o, "contrast", c["id"]) for c in D["contrasts"] if all(tuple(x) in cells for x in c["cells"])]
        mine = [r for r in (ctr or []) if r.get("control_id") == cid]
        _complete(mine if ctr is not None else None, lambda r: _k(r["control_id"], r["observable"], r["unit_type"], r["unit_id"]), keys, fails, cid)
        for r in mine:
            n, d, se, h = (_f(r.get(c)) for c in ("n_k", "d_mean", "se_d", "h"))
            if not (all(math.isfinite(v) for v in (n, d, se, h)) and n == int(n) and 2 <= n <= L["control_seed_extension_max"] and h > 0 and se >= 0):
                fails.append(f"{cid} {r.get('unit_id')}: invalid statistics")
                continue
            eps = L["equivalence_fraction_of_h"] * h
            half = tdist.ppf(1 - L["tost_alpha"], int(n) - 1) * se
            if half > L["underpowered_fraction_of_eps"] * eps:
                fails.append(f"{cid} {r.get('observable')} {r.get('unit_id')}: under-powered")
            elif not (abs(d) + half < eps):
                fails.append(f"{cid} {r.get('observable')} {r.get('unit_id')}: equivalence not shown")
    o7 = C["C7"]["observable"]
    pr = _rows(out, "precision.csv")
    mine = [r for r in (pr or []) if r.get("observable") == o7]
    _complete(mine if pr is not None else None, lambda r: _k(r["observable"], r["unit_type"], r["unit_id"]),
              [_k(o7, "cell", cid_of(w, z)) for (w, z) in allcells] + [_k(o7, "contrast", c["id"]) for c in D["contrasts"]], fails, "C7")
    for r in mine:
        if not _le(r.get("half_width"), L["delta"]):
            fails.append(f"C7 {r.get('unit_id')}: half-width above delta")
    dec = D["statistics"]["decisions"]
    tr = _rows(out, "tests.csv")
    _complete(tr, lambda r: _k(r["observable"], r["contrast_id"]), [_k(o, c["id"]) for o in D["statistics"]["observables"] for c in D["contrasts"]], fails, "tests")
    for r in tr or []:
        p = _f(r.get("p_holm_adj_from_phi"))
        d = r.get("decision")
        if d not in dec or not (0 <= p <= 1) or (d == dec[0] and p > L["holm_alpha"]) or (d == dec[1] and p <= L["holm_alpha"]):
            fails.append(f"test {r.get('observable')} {r.get('contrast_id')}: invalid decision/p")
    res = _finish(out, fails, nested_subtree=None, informational=())
    try:
        with open(os.path.join(out, "summary.json")) as f:
            s = json.load(f)
        if s.get("n_seeds") != n_seeds:
            res["summary_consistent"] = False
            res["agree"] = False
            res["failures"].append("summary n_seeds differs from the protocol addendum value")
    except (OSError, ValueError):
        pass
    return res


def _finish(out, fails, nested_subtree, informational):
    try:
        with open(os.path.join(out, "summary.json")) as f:
            s = json.load(f)
    except (OSError, ValueError) as e:
        return dict(recomputed="NOT PASS", summary_overall=None, agree=False, summary_consistent=False, failures=fails + [f"summary unreadable: {e}"])
    consistent = True
    if s.get("overall_status") == "PASS":
        if s.get("stopped_at") is not None:
            consistent = False
            fails.append("summary says PASS but stopped_at is not null")
        if nested_subtree is not None:
            for p, v in _walk(s.get(nested_subtree, {}), nested_subtree, []):
                if (p.endswith("/passed") and v is not True) or (not p.endswith("/passed") and v != "PASS" and v not in informational):
                    consistent = False
                    fails.append(f"summary says PASS but {p} = {v!r}")
        reg = s.get("registry")
        if isinstance(reg, dict) and any(isinstance(g, dict) and g.get("state") != "PASS" for g in reg.values()):
            consistent = False
            fails.append("summary says PASS but a registry gate is not PASS")
    recomputed = "PASS" if not fails else "NOT PASS"
    agree = consistent and ((recomputed == "PASS") == (s.get("overall_status") == "PASS"))
    return dict(recomputed=recomputed, summary_overall=s.get("overall_status"), agree=agree, summary_consistent=consistent, failures=fails)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--kind", choices=("b3y", "physical"), required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-seeds", type=int, default=None)
    a = ap.parse_args()
    r = verify_b3y(a.out) if a.kind == "b3y" else verify_physical(a.out, a.n_seeds)
    print(json.dumps(dict(recomputed=r["recomputed"], summary_overall=r["summary_overall"], agree=r["agree"],
                          summary_consistent=r["summary_consistent"], n_failures=len(r["failures"]), first_failures=r["failures"][:5]), indent=1))
    sys.exit(0 if r["agree"] else 1)
