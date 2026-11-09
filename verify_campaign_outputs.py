r"""
Independent verifier of campaign outputs (docs/PHYSICAL_CAMPAIGN_PLAN.md section 8, item 4).

Recomputes the gate outcomes from the raw CSV files and the protocol, and checks them against summary.json:
  physical: a four-state outcome per mandatory gate (PASS / FAIL / NOT RUN / NOT EVALUABLE), the exact overall status, the summary's
            gate registry (same ids, same states), the mandatory-gate and PASS counts; an unknown overall status is never accepted.
  b3y     : the binary verdict (PASS / NOT PASS) against the B3Y summary (whose format predates the registry), its exact "PASS" string,
            and the internal consistency of its nested statuses.
Deliberately independent of campaign_status.py and physical_campaign_driver.py: neither is imported. CSV loading, the derivation of the
exact expected key sets from the protocol, row predicates, TOST, Clopper-Pearson / Holm recomputation and the nested-summary walk are
implemented here separately. Key sets are exact: missing, duplicated and unexpected keys all fail; rows are never filtered away.

  python verify_campaign_outputs.py --kind b3y --out results/rgf_b3y_validation             (read-only)
  python verify_campaign_outputs.py --kind physical --out DIR --n-seeds N
Exit status 0 iff the verifier agrees with summary.json.
"""

from __future__ import annotations

import argparse
import csv
import decimal
import hashlib
import json
import math
import os
import sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
SELECTION_ALGORITHMS_V: dict = {}          # independent registry; empty while the protocol leaves the selection rule open
RNG_RULES_V: dict = {}                     # independent registry of bootstrap stream-rule validators; empty: no rule approved
C2_CRITERIA_V: dict = {}                   # independent registry of C2 resolution criteria; empty: none approved
V_RTOL, V_ATOL = 1e-9, 0.0                 # |a-b| <= V_ATOL + V_RTOL*max(|a|,|b|), finite only (docs/CAMPAIGN_STATUS_CONTRACT.md section 2)


def _same(a, b):
    """Independent implementation of the reported-vs-recomputed float agreement rule."""
    x, y = _f(a), _f(b)
    if not (math.isfinite(x) and math.isfinite(y)):
        return False
    return abs(x - y) <= V_ATOL + V_RTOL * max(abs(x), abs(y))
PHYS_OVERALL = ("PASS", "FAIL", "NOT EVALUABLE", "INCOMPLETE")


def _rows(out, name):
    p = os.path.join(out, name)
    if not os.path.isfile(p):
        return None
    with open(p, newline="") as f:
        return [dict(r) for r in csv.DictReader(f)]


def _f(x):
    if x is None or isinstance(x, bool):
        return float("nan")
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _k(*vals):
    out = []
    for v in vals:
        if v is None:
            out.append(None)
            continue
        try:
            out.append(float(v) + 0.0)
        except (TypeError, ValueError):
            out.append(str(v).strip())
    return tuple(out)


def _le(x, lim):
    v = _f(x)
    return math.isfinite(v) and v <= lim


def _true(x):
    return isinstance(x, str) and x.strip() == "True"


def _filled(x):
    return isinstance(x, str) and x.strip() != ""


def _broken(r):
    return None in r or any(v is None for v in r.values())


def _keyset(observed, required, optional=(), multiplicity=None):
    """Exact key-set comparison: returns a list of problems (missing, too many, unexpected)."""
    c = Counter(observed)
    probs = []
    if multiplicity is not None:
        for k, n in multiplicity.items():
            if c.get(k, 0) != n:
                probs.append(f"key {k}: {c.get(k, 0)} rows, protocol {n}")
        probs += [f"unexpected key {k}" for k in c if k not in multiplicity]
        return probs
    req = list(required)
    allowed = set(req) | set(optional)
    probs += [f"missing key {k}" for k in req if c.get(k, 0) == 0]
    probs += [f"key {k} {n} times" for k, n in c.items() if n > 1]
    probs += [f"unexpected key {k}" for k in c if k not in allowed]
    return probs


def _gate(rows, keyfun, required, row_ok, optional=(), multiplicity=None, status_col=None, ok=(), ne=(), open_decision=None, exact_distinct=None):
    """(state, failures) for one gate. exact_distinct: protocol-derived number of distinct keys (all within `optional`)."""
    if rows is None:
        return "NOT RUN", ["table missing"]
    fails, n_ne = [], 0
    try:
        keys = [keyfun(r) for r in rows]
        fails += _keyset(keys, required, optional, multiplicity)
        if exact_distinct is not None and len(set(keys)) != exact_distinct:
            fails.append(f"{len(set(keys))} distinct keys, protocol {exact_distinct}")
    except (KeyError, TypeError) as e:
        fails.append(f"key extraction failed: {e!r}")
    for r in rows:
        if _broken(r):
            fails.append("malformed CSV row")
            continue
        if status_col is not None:
            st = r.get(status_col, "") if status_col in r else ""
            st = st.strip()
            if st in ne:
                n_ne += 1
                continue
            if st not in ok:
                fails.append(f"status {st!r}")
                continue
        why = row_ok(r)
        if why:
            fails.append(why)
    if fails:
        return "FAIL", fails
    if n_ne:
        return "NOT EVALUABLE", [f"{n_ne} NOT EVALUABLE rows"]
    if open_decision:
        return "NOT EVALUABLE", [open_decision]
    return "PASS", []


def _overall(states):
    if states and all(s == "PASS" for s in states):
        return "PASS"
    if not states or "FAIL" in states or any(s not in ("PASS", "NOT RUN", "NOT EVALUABLE") for s in states):
        return "FAIL"
    if "NOT EVALUABLE" in states:
        return "NOT EVALUABLE"
    return "INCOMPLETE"


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


def _summary(out):
    try:
        with open(os.path.join(out, "summary.json")) as f:
            return json.load(f), None
    except (OSError, ValueError) as e:
        return None, f"summary unreadable: {e}"


# ======================================================================
# B3Y (read-only)
# ======================================================================
def verify_b3y(out, protocol_path=None):
    with open(protocol_path or os.path.join(HERE, "docs", "rgf_b3y_protocol.json")) as f:
        P = json.load(f)
    L, X, S = P["limits"], P["xp_gates"], P["point_sets"]
    cfg_of = {p["id"]: p["config"] for p in P["points"]}
    i1, i2, i3 = (int(round(1 / d)) for d in P["discretization"]["dx_levels"])
    fails = []

    def run(label, rows, keyfun, required, row_ok, **kw):
        st, why = _gate(rows, keyfun, required, row_ok, **kw)
        if st != "PASS":
            fails.append(f"{label}: {st}: {why[:2]}")

    run("comparison", _rows(out, "per_energy_comparison.csv"), lambda r: _k(r["point"], r["inv_dx"]),
        [_k(p, i) for p in S["tier1"] for i in (i1, i2)] + [_k(p, i3) for p in S["tier2_dx2048"]],
        lambda r: None if (all(_le(r.get(c), L["S_abs"]) for c in ("S1_dS", "S2_dS", "S3_dS"))
                           and all(_le(r.get(c), L["TR_offdiag_abs"]) for c in ("S4_S1_max", "S4_S2_max"))
                           and all(_le(r.get(f"U_{m}"), L["unitarity"]) and _le(r.get(f"Q_{m}"), L["reciprocity"]) for m in ("M1", "M2", "M3"))
                           and all(_le(r.get(f"R1_backward_error_{m}"), L["sparse_backward_error"]) for m in ("M2", "M3"))
                           and _le(r.get("V_lead_max"), L["lead_tail_guard"])
                           and all(_true(r.get(c)) for c in ("P4_lead_guard_ok", "N1_passed", "P1_rgf_device_bitwise", "P1_M2_samples_bitwise",
                                                             "P1_M3_samples_bitwise", "P1_Eperp_bitwise_across_methods"))) else "comparison row fails",
        status_col="status", ok=("OK",), ne=("NOT EVALUABLE",))
    run("tier3", _rows(out, "tier3_exact_lattice.csv"), lambda r: _k(r["point"], r["inv_dx"]),
        [_k(p, i1) for p in S["tier3_xp_512_validation"]] + [_k(p, i) for p in S["tier3_xp_1024_2048"] for i in (i2, i3)],
        lambda r: None if (all(_true(r.get(c)) for c in ("X4_state_passed", "X4_green_passed", "X5_passed", "X7_passed", "gates_passed"))
                           and _le(r.get("X7_error_indicator"), X["error_indicator_max"])
                           and _le(r.get("X5_fl_vs_state_max_abs"), X["fisher_lee_vs_wavefunction_forms_abs"])
                           and all(_le(r.get(f"eps_{m}"), L["S_abs"]) for m in ("M1", "M2", "M3"))) else "tier3 row fails")
    t3rows = _rows(out, "tier3_exact_lattice.csv") or []
    t3ok = Counter(str(r.get("point")).strip() for r in t3rows if _f(r.get("inv_dx")) == float(i1) and _true(r.get("gates_passed")))
    run("decomposition", _rows(out, "error_decomposition.csv"), lambda r: _k(r["point"], r["inv_dx"], r["method"]),
        [_k(p, i, m) for p in S["tier3_xp_1024_2048"] for i in (i2, i3) for m in ("M2", "M3")],
        lambda r: None if (all(math.isfinite(_f(r.get(c))) for c in ("total", "bulk_rounding", "lead_data", "solver", "identity_residual_three_term", "e_ind_combined"))
                           and all(_true(r.get(c)) for c in ("gate_fl_passed", "gate_hybrid_passed", "mirror_source_bitwise", "mirror_extraction_bitwise",
                                                             "M2_recomputed_bitwise", "M3_recomputed_bitwise"))) else "decomposition row not evaluable/failing",
        optional=[_k(p, i1, m) for p in S["tier3_xp_1024_2048"] if p in S["tier3_xp_512_validation"] and t3ok.get(p, 0) == 1 for m in ("M2", "M3")],
        status_col="status", ok=("", "OK"), ne=("NOT EVALUABLE",))
    xr = _rows(out, "xp_gates.csv")
    allpts = list(dict.fromkeys(S["tier3_xp_512_validation"] + S["tier3_xp_1024_2048"]))
    exp = {"G0": [_k(p, "") for p in S["tier1"]], "X1": [_k(p, i) for p in allpts for i in (i1, i2, i3)],
           "X2": [_k(p, i) for p in S["tier3_xp_1024_2048"] for i in (i1, i2, i3)], "X3": [_k(p, i) for p in S["tier3_xp_1024_2048"] for i in (i1, i2, i3)],
           "X6": [_k(p, i) for p in S["tier3_xp_1024_2048"] for i in (i1, i2)]}
    passed = lambda r: None if _true(r.get("passed")) else "not passed"
    for g, keys in exp.items():
        run(f"xp {g}", None if xr is None else [r for r in xr if r.get("gate") == g], lambda r: _k(r["case"], r["inv_dx"]), keys, passed)
    for g in ("X0", "X8"):
        run(f"xp {g}", None if xr is None else [r for r in xr if r.get("gate") == g], lambda r: _k(r["gate"]), [_k(g)], passed)
    if xr is not None and any(r.get("gate") not in set(exp) | {"X0", "X8"} for r in xr):
        fails.append("xp_gates: rows of an undeclared gate")
    run("K1", _rows(out, "buffer_control.csv"), lambda r: _k(r["point"]), [_k(p) for p in S["buffer_control"]],
        lambda r: None if (_true(r.get("K1_passed_dT1")) and _true(r.get("S_within_geometry_passed")) and _le(r.get("abs_dT1_M1"), L["dT1_buffer"])
                           and _le(r.get("abs_dT1_M3"), L["dT1_buffer"])) else "K1 fails")
    run("H1", _rows(out, "solver_spread.csv"), lambda r: _k(r["point"], r["inv_dx"]), [_k(p, i) for p in S["tier2_dx2048"] for i in (i1, i2, i3)],
        lambda r: None if (_true(r.get("passed")) and _le(r.get("T1_spread"), L["T1_method_spread_for_solver_independence"])) else "H1 fails")
    run("O1", _rows(out, "refinement_convergence.csv"), lambda r: _k(r["point"], r["method"]),
        [_k(p, m) for p in S["tier2_dx2048"] for m in ("M1", "M2", "M3", "M4")],
        lambda r: None if math.isfinite(_f(r.get("observed_order"))) else "order missing")
    pi = _rows(out, "potential_identity.csv")
    mult = Counter()
    for p in S["tier1"]:
        mult[_k(cfg_of[p], i1, 0.2)] += 1
        mult[_k(cfg_of[p], i2, 0.2)] += 1
    for p in S["tier2_dx2048"]:
        mult[_k(cfg_of[p], i3, 0.2)] += 1
    for p in S["buffer_control"]:
        mult[_k(cfg_of[p], i2, 0.4)] += 1
    p2 = lambda r: "_vs_" in str(r.get("inv_dx"))
    run("P1/P4", None if pi is None else [r for r in pi if not p2(r)], lambda r: _k(r["config"], r["inv_dx"], r["buffer_requested"]), [],
        lambda r: None if (_true(r.get("P1_rgf_device_vs_direct_bitwise")) and _true(r.get("P4_lead_guard_ok")) and _le(r.get("V_lead_max"), L["lead_tail_guard"])) else "P1/P4 fails",
        multiplicity=dict(mult))
    run("P2", None if pi is None else [r for r in pi if p2(r)], lambda r: _k(r["config"], r["inv_dx"]),
        [_k(c, f"{i}_vs_{i1}") for c in sorted({cfg_of[p] for p in S["tier2_dx2048"]}) for i in (i2, i3)],
        lambda r: None if (_true(r.get("P2_potential_bitwise")) and _true(r.get("coincident_x_bitwise"))) else "P2 fails")
    recomputed = "PASS" if not fails else "NOT PASS"
    s, err = _summary(out)
    if s is None:
        return dict(recomputed=recomputed, summary_overall=None, agree=False, summary_consistent=False, failures=fails + [err])
    consistent, why = True, []
    so = s.get("overall_status")
    if not isinstance(so, str) or not so.strip():
        consistent = False
        why.append(f"malformed overall_status {so!r}")
    if so == "PASS":
        if s.get("stopped_at") is not None:
            consistent = False
            why.append("PASS with stopped_at set")
        for p, v in _walk(s.get("criteria", {}), "criteria", []):
            if (p.endswith("/passed") and v is not True) or (not p.endswith("/passed") and v not in ("PASS", "CLASSIFICATION ONLY")):
                consistent = False
                why.append(f"PASS but {p} = {v!r}")
    agree = consistent and ((recomputed == "PASS") == (so == "PASS"))
    return dict(recomputed=recomputed, summary_overall=so, agree=agree, summary_consistent=consistent, failures=fails + why)


# ======================================================================
# physical campaign
# ======================================================================
def _cp(K, B, a):
    from scipy.stats import beta
    lo = 0.0 if K == 0 else float(beta.ppf(a / 2, K, B - K + 1))
    hi = 1.0 if K == B else float(beta.ppf(1 - a / 2, K + 1, B - K))
    return lo, hi


def _holm(p: dict, alpha):
    """Rejection set: the largest prefix j of the ascending order such that every p_(i) <= alpha/(m - i + 1) for i <= j."""
    vals = sorted(p.items(), key=lambda kv: kv[1])
    m = len(vals)
    j = 0
    while j < m and vals[j][1] <= alpha / (m - j):
        j += 1
    return {k for k, _ in vals[:j]}


def _stats_ok(rows, D, order):
    st = D["statistics"]
    al, B1, B2, a1, a2 = st["holm_alpha"], st["B_stage1"], st["B_stage2"], st["mc_alpha_per_test_stage1"], st["mc_alpha_per_test_stage2"]
    REJ, NREJ, AMB = st["decisions"]
    close = _same
    T, out = {}, []
    for r in rows:
        key = (str(r.get("observable")).strip(), str(r.get("contrast_id")).strip())
        try:
            K1, K, B, stg = (int(_f(r[c])) for c in ("K_stage1", "K", "B", "stage"))
            if any(_f(r[c]) != v for c, v in (("K_stage1", K1), ("K", K), ("B", B), ("stage", stg))):
                raise ValueError
        except (KeyError, ValueError, OverflowError):
            out.append(f"{key}: counts not integers")
            continue
        g = {c: _f(r.get(c)) for c in ("theta_hat", "ci95_sampling_lo", "ci95_sampling_hi", "p_hat", "mc_se", "mc_ci_sim99_lo", "mc_ci_sim99_hi",
                                         "p_lo_final", "p_hi_final", "holm_rank", "holm_threshold", "p_holm_adj_from_phi")}
        if not all(math.isfinite(v) for v in g.values()) or g["ci95_sampling_lo"] > g["ci95_sampling_hi"]:
            out.append(f"{key}: missing/non-finite values or reversed sampling CI")
            continue
        if not ((stg == 1 and B == B1 and K == K1) or (stg == 2 and B == B2)) or not (0 <= K1 <= B1 and 0 <= K <= B):
            out.append(f"{key}: stage/B/K inconsistent")
            continue
        ph = (K + 1) / (B + 1)
        lo, hi = _cp(K, B, a1 if stg == 1 else a2)
        checks = [close(g["p_hat"], ph), close(g["mc_se"], math.sqrt(ph * (1 - ph) / B)), close(g["mc_ci_sim99_lo"], lo),
                  close(g["mc_ci_sim99_hi"], hi), close(g["p_lo_final"], lo), close(g["p_hi_final"], max(hi, ph))]
        if not all(checks):
            out.append(f"{key}: MC quantities inconsistent")
            continue
        if K == 0 and not (g["mc_ci_sim99_lo"] == 0.0 and g["p_lo_final"] == 0.0):
            out.append(f"{key}: K = 0 requires exact zero lower endpoint")
            continue
        if K == B and not (g["mc_ci_sim99_hi"] == 1.0 and g["p_hi_final"] == 1.0):
            out.append(f"{key}: K = B requires exact unit upper endpoint")
            continue
        if not (0.0 <= g["p_lo_final"] and g["p_lo_final"] <= g["p_hat"] and g["p_hat"] <= g["p_hi_final"] and g["p_hi_final"] <= 1.0):
            out.append(f"{key}: endpoint ordering violated")
            continue
        l1, h1 = _cp(K1, B1, a1)
        # decision quantities from the integer counts (never from the reported floats)
        T[key] = (stg, l1, max(h1, (K1 + 1) / (B1 + 1)), g, r.get("decision"), lo, max(hi, ph))
    if out:
        return out
    if set(T) != set(order):
        return ["test set differs from the 68-test family"]
    A1 = _holm({k: v[1] for k, v in T.items()}, al) - _holm({k: v[2] for k, v in T.items()}, al)
    out += [f"{k}: Stage-2 selection violated" for k, v in T.items() if (v[0] == 2) != (k in A1)]
    phi = {k: v[6] for k, v in T.items()}
    R = _holm(phi, al)
    A = _holm({k: v[5] for k, v in T.items()}, al) - R
    pos = {k: i for i, k in enumerate(order)}
    ranked = sorted(T, key=lambda k: (phi[k], pos[k]))
    m = len(ranked)
    adj_prev = 0.0
    for i, k in enumerate(ranked):
        g = T[k][3]
        adj = max(adj_prev, min(1.0, (m - i) * phi[k]))
        adj_prev = adj
        if int(g["holm_rank"]) != i + 1 or not close(g["holm_threshold"], al / (m - i)) or not close(g["p_holm_adj_from_phi"], adj):
            out.append(f"{k}: rank/threshold/adjusted p wrong")
        label = T[k][4]
        if k in R:
            if label != REJ:
                out.append(f"{k}: in R_claim but labelled {label!r}")
        elif k in A:
            if label != AMB:
                out.append(f"{k}: in A_final but labelled {label!r}")
        elif label != NREJ:
            out.append(f"{k}: outside R_claim and A_final but labelled {label!r}")
    return out


def _digits(x):
    return isinstance(x, str) and len(x.strip()) > 0 and all(ch in "0123456789" for ch in x.strip())


def _np_version():
    import numpy
    return numpy.__version__


def _subset_count(fraction, n_cells, n_seeds):
    """Independent exact rule: ceiling of fraction * cells * N in decimal arithmetic."""
    v = decimal.Decimal(str(fraction)) * n_cells * n_seeds
    return int(v.to_integral_value(rounding=decimal.ROUND_CEILING))


def _c0_select_v(G, seeds, n, seed):
    """Independent implementation of the approved C0 selection: units ordered by nested ascending loops over W, zeta, seed."""
    import numpy as np
    units = [(w, z, s) for w in sorted(G["W"]) for z in sorted(G["zeta"]) for s in sorted(seeds)]
    rng = np.random.Generator(np.random.PCG64(np.random.SeedSequence(seed)))
    chosen = sorted(rng.permutation(len(units))[:n].tolist())
    return [units[i] for i in chosen]


SELECTION_ALGORITHMS_V["SEEDSEQ-PCG64-PERMUTATION-C0"] = _c0_select_v


def _rule_a_v(meta, spec, n_seeds):
    """Independent validator of the approved bootstrap rule (SeedSequence root, children via spawn(2), PCG64)."""
    import numpy as np
    bad = []
    E = int(spec["root_entropy"])
    if meta["bitgen"] != "PCG64" or spec.get("bit_generator") != "PCG64":
        bad.append("bit generator not PCG64")
    if meta["numpy"] != spec.get("numpy_version"):
        bad.append("numpy version metadata not the pin")
    if meta["api"] != spec.get("resampling_api"):
        bad.append("resampling API metadata not the approved call")
    children = np.random.SeedSequence(E).spawn(2)
    for st, rec in ((1, meta["one"]), (2, meta["two"])):
        if rec is None:
            continue
        sid, seed, probe = rec
        want = (spec["stream_id_template"].replace("{root_entropy}", str(E)).replace("{spawn_index}", str(children[st - 1].spawn_key[0]))
                .replace("{numpy_version}", spec["numpy_version"]))
        if seed != E:
            bad.append(f"Stage-{st} seed is not the root entropy")
        if sid != want:
            bad.append(f"Stage-{st} stream identifier not of the approved form")
        if _np_version() == spec["numpy_version"]:
            g = np.random.Generator(np.random.PCG64(children[st - 1]))
            first = g.integers(low=0, high=n_seeds, size=n_seeds, dtype=np.int64, endpoint=False)
            if probe != hashlib.sha256(first.astype("<i8").tobytes()).hexdigest():
                bad.append(f"Stage-{st} probe not reproduced")
    return bad


RNG_RULES_V["SEEDSEQ-PCG64-SPAWN-A"] = _rule_a_v


def _digits(x):
    return isinstance(x, str) and len(x.strip()) > 0 and all(ch in "0123456789" for ch in x.strip())


def _rng_gate(rows, D, order, n_seeds):
    """Independent check of the bootstrap stream metadata (stage-specific streams, consistency, distinctness, approved rule)."""
    inv_cols = ("bootstrap_rng_rule", "bootstrap_bit_generator", "bootstrap_numpy_version", "bootstrap_resampling_api")
    cols = inv_cols + ("bootstrap_rng_stream_stage1", "bootstrap_seed_stage1", "bootstrap_rng_stream_stage2", "bootstrap_seed_stage2",
                       "bootstrap_probe_stage1", "bootstrap_probe_stage2")
    st, why = _gate(rows, lambda r: _k(r["observable"], r["contrast_id"]), [_k(*k) for k in order], lambda r: None)
    if st != "PASS":
        return st, why
    bad, one, two, invs = [], [], [], []
    for r in rows:
        if not all(c in r for c in cols):
            bad.append("RNG columns missing")
            break
        stage = _f(r.get("stage"))
        if not (all(_filled(r[c]) for c in inv_cols) and _filled(r["bootstrap_rng_stream_stage1"]) and _digits(r["bootstrap_seed_stage1"])
                and _filled(r["bootstrap_probe_stage1"])):
            bad.append("Stage-1 stream/seed/probe missing or malformed")
            continue
        invs.append(tuple(r[c].strip() for c in inv_cols))
        one.append((r["bootstrap_rng_stream_stage1"].strip(), int(r["bootstrap_seed_stage1"].strip()), r["bootstrap_probe_stage1"].strip()))
        any2 = any(_filled(r[c]) for c in ("bootstrap_rng_stream_stage2", "bootstrap_seed_stage2", "bootstrap_probe_stage2"))
        if stage == 1:
            if any2:
                bad.append("Stage-1 test with Stage-2 stream metadata")
        elif stage == 2:
            if not (_filled(r["bootstrap_rng_stream_stage2"]) and _digits(r["bootstrap_seed_stage2"]) and _filled(r["bootstrap_probe_stage2"])):
                bad.append("Stage-2 stream/seed/probe missing or malformed")
                continue
            two.append((r["bootstrap_rng_stream_stage2"].strip(), int(r["bootstrap_seed_stage2"].strip()), r["bootstrap_probe_stage2"].strip()))
        else:
            bad.append("invalid stage")
    if not bad:
        if len(set(invs)) != 1:
            bad.append("rule / bit generator / numpy / API metadata not unique")
        if len(set(one)) != 1:
            bad.append("Stage-1 stream not shared by all tests")
        if len(set(two)) > 1:
            bad.append("Stage-2 stream not shared by all Stage-2 tests")
        if one and two and one[0][0] == two[0][0]:
            bad.append("Stage-2 stream identifier identical to Stage 1")
    if bad:
        return "FAIL", bad
    spec = D.get("bootstrap_rng")
    if spec is None:
        return "NOT EVALUABLE", ["protocol decision required: bootstrap RNG rule"]
    if invs[0][0] != str(spec.get("rule")):
        return "FAIL", ["unapproved stream rule"]
    fn = RNG_RULES_V.get(spec.get("rule"))
    if fn is None:
        return "NOT EVALUABLE", ["no validator for the approved stream rule"]
    if spec.get("root_entropy") is None:
        return "NOT EVALUABLE", ["protocol decision required: root entropy not recorded"]
    meta = dict(bitgen=invs[0][1], numpy=invs[0][2], api=invs[0][3], one=one[0], two=two[0] if two else None)
    try:
        msgs = fn(meta, spec, n_seeds)
    except (KeyError, TypeError, ValueError) as e:                            # incomplete rule specification in the protocol: fail closed
        return "FAIL", [f"approved rule specification incomplete: {type(e).__name__}: {e}"]
    if msgs:
        return "FAIL", msgs
    if _np_version() != spec.get("numpy_version"):
        return "NOT EVALUABLE", ["installed numpy differs from the pin: streams not reproducible here"]
    return "PASS", []


def _template_problems(b):
    """Independent schema check of the rule-A identifier template (hand-written brace scan; no str.format machinery)."""
    t = b.get("stream_id_template")
    if not isinstance(t, str) or t.strip() == "":
        return ["identifier template missing"]
    names, out, i = [], [], 0
    while i < len(t):
        ch = t[i]
        if ch == "}":
            return ["stray closing brace in identifier template"]
        if ch == "{":
            j = t.find("}", i + 1)
            if j < 0:
                return ["unterminated brace in identifier template"]
            inner = t[i + 1:j]
            if inner == "" or "{" in inner or not all(c.isascii() and (c.isalnum() or c == "_") for c in inner):
                return [f"malformed placeholder {{{inner}}}"]
            names.append(inner)
            i = j + 1
            continue
        i += 1
    for f in ("root_entropy", "spawn_index", "numpy_version"):
        if names.count(f) != 1:
            out.append(f"placeholder {f} occurs {names.count(f)} times")
    out += [f"unknown placeholder {n}" for n in sorted(set(names)) if n not in ("root_entropy", "spawn_index", "numpy_version")]
    comps = t.split("|")
    expect = (b.get("rule"), "entropy={root_entropy}", "spawn_key=({spawn_index},)", "bitgen=" + str(b.get("bit_generator")), "numpy={numpy_version}")
    if not (isinstance(b.get("rule"), str) and isinstance(b.get("bit_generator"), str)) or len(comps) != 5 or tuple(comps) != expect:
        out.append("fixed identifier components malformed")
    return out


def _protocol_problems(D):
    """Inadmissible recorded protocol values (the aggregator refuses to build its registry for these: overall FAIL)."""
    bad = []
    b = D.get("bootstrap_rng") or {}
    c = D.get("reproducibility_selection") or {}
    if b and b.get("rule") == "SEEDSEQ-PCG64-SPAWN-A":
        bad += _template_problems(b)                              # independent of whether root_entropy has been recorded
    if b.get("root_entropy") is not None and not _digits(str(b["root_entropy"])):
        bad.append("bootstrap root_entropy malformed")
    if c.get("seed") is not None:
        if not _digits(str(c["seed"])):
            bad.append("C0 selection seed malformed")
        elif b.get("root_entropy") is not None and _digits(str(b["root_entropy"])) and int(c["seed"]) == int(b["root_entropy"]):
            bad.append("C0 selection seed not dedicated (equals the bootstrap root entropy)")
    return bad


PHYS_TABLES = ("stages", "realizations", "reproducibility", "solver_cross_check", "adjudication", "b3y_reproduction", "clean_guide", "controls",
               "precision", "tests")


def recompute_physical(out, n_seeds, D):
    """Independent four-state outcome per mandatory gate (gate ids are the published contract of the registry)."""
    from scipy.stats import t as tdist
    G, L, C, sd = D["grid"], D["limits"], D["controls"], D["seeds"]
    cells = [(w, z) for w in G["W"] for z in G["zeta"]]
    corners = [tuple(c) for c in G["corner_cells"]]
    pick = lambda spec: cells if spec == "all" else corners if spec == "corners" else [tuple(c) for c in spec]
    seeds = list(range(sd["confirmatory_first"], sd["confirmatory_first"] + n_seeds))
    E = list(range(D["energies"]["n_energies"]))
    tb = {t: _rows(out, f"{t}.csv") for t in PHYS_TABLES}
    R = {}
    R["PHYS.stages.completed"] = _gate(tb["stages"], lambda r: _k(r["stage"]), [_k(s) for s in D["stages"]], lambda r: None,
                                       status_col="status", ok=("COMPLETED",))
    R["PHYS.C0.realizations"] = _gate(
        tb["realizations"], lambda r: _k(r["W"], r["zeta"], r["seed"], r["energy_index"]), [_k(w, z, s, e) for (w, z) in cells for s in seeds for e in E],
        lambda r: None if (_le(r.get("unitarity"), L["unitarity"]) and _le(r.get("reciprocity"), L["reciprocity"]) and _le(r.get("lead_tail_max"), L["lead_tail_guard"])
                           and math.isfinite(_f(r.get("n_open"))) and _f(r.get("n_open")) == _f(r.get("n_open_analytic")) and _filled(r.get("v_hash"))
                           and math.isfinite(_f(r.get("tau"))) and math.isfinite(_f(r.get("T1")))) else "realization fails",
        status_col="status", ok=("OK",), ne=("NOT EVALUABLE",))
    n_sub = _subset_count(L["reproducibility_subset_fraction"], len(cells), n_seeds)
    ens = {_k(w, z, s) for (w, z) in cells for s in seeds}
    rp = tb["reproducibility"]
    sel = D.get("reproducibility_selection")
    rep_ok = lambda r: None if (_true(r.get("bitwise_equal")) and _filled(r.get("selection_algorithm")) and _filled(r.get("selection_seed"))
                                and _filled(r.get("selection_numpy_version"))) else "repro row fails"
    alg = None if sel is None else SELECTION_ALGORITHMS_V.get(sel.get("algorithm"))
    rkey = lambda r: _k(r["W"], r["zeta"], r["seed"])
    sel_ok = lambda r: None if (str(r.get("selection_algorithm")).strip() == sel["algorithm"] and _digits(r.get("selection_seed"))
                                and int(r["selection_seed"].strip()) == int(sel["seed"])
                                and str(r.get("selection_numpy_version")).strip() == sel["numpy_version"]) else "selection not the protocol's"
    if alg is None or sel.get("seed") is None:
        R["PHYS.C0.reproducibility"] = _gate(rp, rkey, [], rep_ok, optional=sorted(ens), exact_distinct=n_sub,
                                             open_decision="protocol decision required: C0 selection rule or seed")
    elif _np_version() != sel.get("numpy_version"):
        R["PHYS.C0.reproducibility"] = _gate(rp, rkey, [], lambda r: rep_ok(r) or sel_ok(r), optional=sorted(ens), exact_distinct=n_sub,
                                             open_decision="installed numpy differs from the pin: the exact subset cannot be reproduced here")
    else:
        req = [_k(*k) for k in alg(G, seeds, n_sub, int(sel["seed"]))]
        R["PHYS.C0.reproducibility"] = _gate(rp, rkey, req, lambda r: rep_ok(r) or sel_ok(r))
    c2 = C["C2"]
    R["PHYS.C2.solver_cross_check"] = _gate(tb["solver_cross_check"], lambda r: _k(r["W"], r["zeta"], r["seed"], r["energy_index"]),
                                            [_k(w, z, s, e) for (w, z) in pick(c2["cells"]) for s in seeds[:c2["n_seeds"]] for e in E],
                                            lambda r: None if _le(r.get("max_abs_dS"), L["solver_cross_check_S_abs"]) else "C2 exceedance")
    flagged = list(dict.fromkeys(_k(r.get("W"), r.get("zeta"), r.get("seed"), r.get("energy_index")) for r in (tb["solver_cross_check"] or [])
                                 if not _le(r.get("max_abs_dS"), L["solver_cross_check_S_abs"])))
    adj_numbers = lambda r: None if (math.isfinite(_f(r.get("eps_M1"))) and math.isfinite(_f(r.get("eps_M3")))) else "adjudication numbers missing"
    akey = lambda r: _k(r["W"], r["zeta"], r["seed"], r["energy_index"])
    R["PHYS.C2.adjudication_recorded"] = _gate(tb["adjudication"], akey, flagged, adj_numbers, status_col="status", ok=("ADJUDICATED",), ne=("NOT EVALUABLE",))
    crit = D.get("C2_adjudication_criteria")
    cfun = None if crit is None else C2_CRITERIA_V.get(crit.get("criterion"))
    open_res = None
    if flagged:
        open_res = ("protocol decision required: no C2 resolution criterion" if crit is None
                    else None if cfun is not None else "C2 resolution criterion not implemented")
    R["PHYS.C2.adjudication_resolved"] = _gate(tb["adjudication"], akey, flagged,
                                               (lambda r: adj_numbers(r) or cfun(r, crit)) if cfun is not None else adj_numbers,
                                               status_col="status", ok=("ADJUDICATED",), ne=("NOT EVALUABLE",), open_decision=open_res)
    R["PHYS.C5.b3y_reproduction"] = _gate(tb["b3y_reproduction"], lambda r: _k(r["zeta"], r["energy_index"]), [_k(z, e) for z in C["C5"]["zeta"] for e in E],
                                          lambda r: None if _le(r.get("abs_dT1"), L["b3y_reproduction_abs_dT1"]) else "C5 exceeds")
    R["PHYS.C6.clean_guide"] = _gate(tb["clean_guide"], lambda r: _k(r["energy_index"]), [_k(e) for e in E],
                                     lambda r: None if _le(r.get("max_abs_Tnn_minus_1"), L["clean_guide_Tnn_abs"]) else "C6 exceeds")
    cid = lambda w, z: f"W{float(w)!r}_z{float(z)!r}"

    def tost(r):
        n, d, se, h = (_f(r.get(c)) for c in ("n_k", "d_mean", "se_d", "h"))
        if not (all(math.isfinite(v) for v in (n, d, se, h)) and n == int(n) and 2 <= n <= L["control_seed_extension_max"] and h > 0 and se >= 0):
            return "invalid control statistics"
        eps = L["equivalence_fraction_of_h"] * h
        half = float(tdist.ppf(1 - L["tost_alpha"], int(n) - 1)) * se
        if half > L["underpowered_fraction_of_eps"] * eps:
            return "under-powered"
        return None if abs(d) + half < eps else "equivalence not shown"
    ctr = tb["controls"]
    for c_ in ("C1", "C1b", "C3", "C4"):
        cs_ = set(pick(C[c_]["cells"]))
        keys = []
        for o in C[c_]["observables"]:
            keys += [_k(c_, o, "cell", cid(w, z)) for (w, z) in sorted(cs_)]
            keys += [_k(c_, o, "contrast", c["id"]) for c in D["contrasts"] if all(tuple(x) in cs_ for x in c["cells"])]
        R[f"PHYS.{c_}.equivalence"] = _gate(None if ctr is None else [r for r in ctr if str(r.get("control_id")).strip() == c_],
                                            lambda r: _k(r["control_id"], r["observable"], r["unit_type"], r["unit_id"]), keys, tost)
    o7 = C["C7"]["observable"]
    R["PHYS.C7.precision"] = _gate(tb["precision"], lambda r: _k(r["observable"], r["unit_type"], r["unit_id"]),
                                   [_k(o7, "cell", cid(w, z)) for (w, z) in cells] + [_k(o7, "contrast", c["id"]) for c in D["contrasts"]],
                                   lambda r: None if _le(r.get("half_width"), L["delta"]) else "half-width above delta")
    order = [(o, c["id"]) for o in D["statistics"]["observables"] for c in D["contrasts"]]
    st, why = _gate(tb["tests"], lambda r: _k(r["observable"], r["contrast_id"]), [_k(*k) for k in order], lambda r: None,
                    status_col="decision", ok=tuple(D["statistics"]["decisions"]))
    if st == "PASS":
        bad = _stats_ok(tb["tests"], D, order)
        if bad:
            st, why = "FAIL", bad
    R["PHYS.statistics.tests_complete"] = (st, why)
    R["PHYS.statistics.bootstrap_rng"] = _rng_gate(tb["tests"], D, order, n_seeds)
    # per-table coverage: every row owned by some gate
    for t in PHYS_TABLES:
        rows = tb[t]
        if rows is None:
            R[f"TABLE.{t}.csv.coverage"] = ("NOT RUN", ["table missing"])
        elif t == "controls" and any(str(r.get("control_id")).strip() not in ("C1", "C1b", "C3", "C4") for r in rows):
            R[f"TABLE.{t}.csv.coverage"] = ("FAIL", ["control rows of an undeclared control"])
        else:
            R[f"TABLE.{t}.csv.coverage"] = ("PASS", [])
    return R


def verify_physical(out, n_seeds, draft_path=None):
    with open(draft_path or os.path.join(HERE, "docs", "physical_campaign_gates_draft.json")) as f:
        D = json.load(f)
    sd = D["seeds"]
    s, err = _summary(out)
    if not (isinstance(n_seeds, int) and not isinstance(n_seeds, bool) and sd["N_min"] <= n_seeds <= sd["N_max"]):
        return dict(recomputed="FAIL", summary_overall=(s or {}).get("overall_status"), agree=False, summary_consistent=False,
                    gate_states={}, failures=[f"N = {n_seeds!r} not admissible"])
    bad_protocol = _protocol_problems(D)
    if bad_protocol:
        so = (s or {}).get("overall_status")
        return dict(recomputed="FAIL", summary_overall=so, agree=(so == "FAIL"), summary_consistent=(so == "FAIL"), gate_states={},
                    failures=bad_protocol)
    R = recompute_physical(out, n_seeds, D)
    states = {g: v[0] for g, v in R.items()}
    recomputed = _overall(list(states.values()))
    fails = [f"{g}: {v[0]}: {v[1][:2]}" for g, v in R.items() if v[0] != "PASS"]
    if s is None:
        return dict(recomputed=recomputed, summary_overall=None, agree=False, summary_consistent=False, gate_states=states, failures=fails + [err])
    why = []
    so = s.get("overall_status")
    if so not in PHYS_OVERALL:
        why.append(f"overall_status {so!r} is not one of {PHYS_OVERALL}")
    elif so != recomputed:
        why.append(f"overall_status {so!r} != recomputed {recomputed!r}")
    reg = s.get("registry")
    if not isinstance(reg, dict):
        why.append("registry missing or malformed")
    else:
        if set(reg) != set(states):
            why.append(f"registry ids differ: missing {sorted(set(states) - set(reg))[:3]}, unexpected {sorted(set(reg) - set(states))[:3]}")
        for g, st in states.items():
            e = reg.get(g)
            if not isinstance(e, dict) or e.get("state") != st:
                why.append(f"registry {g}: {None if not isinstance(e, dict) else e.get('state')!r} != recomputed {st!r}")
        if so == "PASS" and any(not isinstance(e, dict) or e.get("state") != "PASS" for e in reg.values()):
            why.append("summary says PASS with a non-PASS registry entry")
    if s.get("n_mandatory_gates") != len(states):
        why.append(f"n_mandatory_gates {s.get('n_mandatory_gates')!r} != {len(states)}")
    if s.get("n_gates_pass") != sum(v == "PASS" for v in states.values()):
        why.append(f"n_gates_pass {s.get('n_gates_pass')!r} != {sum(v == 'PASS' for v in states.values())}")
    if s.get("n_seeds") != n_seeds:
        why.append("summary n_seeds differs from the protocol addendum value")
    if so == "PASS" and (s.get("stopped_at") is not None or s.get("campaign_status") != "COMPLETED"):
        why.append("summary says PASS but the run did not complete")
    return dict(recomputed=recomputed, summary_overall=so, agree=not why, summary_consistent=not why, gate_states=states, failures=fails + why)


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
