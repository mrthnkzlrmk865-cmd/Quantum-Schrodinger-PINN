r"""
Fail-closed orchestration for the physical ensemble campaign (docs/PHYSICAL_CAMPAIGN_PLAN.md section 8).

Scope of this module (prerequisite work item): stage ordering, table writing, early-termination handling and the fail-closed overall
status. It contains NO physics: the stage implementations (RGF solves, controls, statistics) are injected as callables and are not
part of this work item. With no implementation a stage is recorded as NOT RUN and the run stops there (fail closed).

  run_campaign(draft, out_dir, n_seeds, stage_impls)   -> summary dict; writes <table>.csv for every declared table (header-only if
                                                           empty), stages.csv, gate_registry.csv, summary.json
A stage implementation receives a context dict and returns {table_name: [row dicts]}. Every returned row is written: rows are never
filtered, failed and NOT EVALUABLE rows stay in the tables. An exception or an undeclared table stops the run (remaining stages NOT RUN).
The overall status comes only from campaign_status.aggregate over the protocol-derived gate registry.
"""

from __future__ import annotations

import csv
import datetime as _dt
import json
import os

import campaign_status as cs

HERE = os.path.dirname(os.path.abspath(__file__))
HISTORICAL_RESULT_DIRS = tuple(os.path.normpath(os.path.join(HERE, "results", d)) for d in (
    "rgf_b1_validation", "rgf_b2_validation", "rgf_b3_validation", "rgf_b3x_validation", "rgf_b3y_validation",
    "ode1d_validation", "validation", "reproduction", "transport", "localization", "analysis", "pinn"))


def check_output_dir(out):
    o = os.path.normpath(os.path.abspath(out))
    for d in HISTORICAL_RESULT_DIRS:
        if o == d or o.startswith(d + os.sep):
            raise SystemExit(f"refusing to write into the historical result directory {d}")
    return o


def _write_table(path, columns, rows):
    cols = list(columns) + [k for k in dict.fromkeys(k for r in rows for k in r) if k not in columns]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: (repr(float(v)) if isinstance(v, float) else v) for k, v in r.items()})


def run_campaign(draft: dict, out_dir: str, n_seeds: int, stage_impls: dict, log=print) -> dict:
    out = check_output_dir(out_dir)
    os.makedirs(out, exist_ok=True)
    schemas = draft["tables"]
    tables = {t: [] for t in schemas if t != "stages"}
    stage_rows, stopped = [], None
    ctx = dict(draft=draft, n_seeds=n_seeds, out_dir=out)
    for st in draft["stages"]:
        if stopped is not None:
            stage_rows.append(dict(stage=st, status="NOT RUN", reason=f"run stopped at {stopped}"))
            continue
        impl = stage_impls.get(st)
        if impl is None:
            stage_rows.append(dict(stage=st, status="NOT RUN", reason="no implementation supplied"))
            stopped = st
            continue
        try:
            res = impl(ctx) or {}
            bad = [t for t in res if t not in tables]
            if bad:
                raise ValueError(f"stage returned undeclared tables {bad}")
            for t, rows in res.items():
                tables[t].extend(dict(r) for r in rows)
            stage_rows.append(dict(stage=st, status="COMPLETED", reason=""))
        except Exception as e:                                                # noqa: BLE001  (fail closed; recorded, run stops)
            stage_rows.append(dict(stage=st, status="FAILED", reason=f"{type(e).__name__}: {e}"))
            stopped = st
        log(f"stage {st}: {stage_rows[-1]['status']}")
    for t, rows in tables.items():
        _write_table(os.path.join(out, f"{t}.csv"), schemas[t], rows)
    _write_table(os.path.join(out, "stages.csv"), schemas["stages"], stage_rows)
    try:
        specs = cs.physical_gate_spec(draft, n_seeds)
        agg = cs.aggregate(specs, out)
    except Exception as e:                                                    # noqa: BLE001
        agg = dict(overall_status=cs.FAIL, n_mandatory=None, n_pass=None, problems=[f"gate registry could not be built: {type(e).__name__}: {e}"], registry=[])
    cs.write_registry(os.path.join(out, "gate_registry.csv"), agg["registry"])
    summary = dict(
        overall_status=agg["overall_status"], campaign_status="COMPLETED" if stopped is None else f"STOPPED at {stopped}", stopped_at=stopped,
        n_seeds=n_seeds, n_mandatory_gates=agg["n_mandatory"], n_gates_pass=agg["n_pass"], problems=agg["problems"],
        registry={g["gate_id"]: dict(state=g["state"], reason=g["reason"]) for g in agg["registry"]},
        row_counts={t: len(r) for t, r in tables.items()}, stages=stage_rows,
        protocol_status=draft.get("status"), source_plan_commit=draft.get("source_plan_commit"),
        generated_utc=_dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"))
    with open(os.path.join(out, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2)
    return summary


if __name__ == "__main__":
    raise SystemExit("physical_campaign_driver: the stage implementations (physics) are not part of this work item; "
                     "the module is exercised only by tests/test_campaign_status.py with synthetic stages.")
