"""Report writer (JSON + Markdown). The report never states that W/Z is validated; it reports software-test outcomes."""
import json
import pathlib
from fractions import Fraction


def _default(o):
    if isinstance(o, Fraction):
        return str(o)
    if isinstance(o, (set, frozenset)):
        return sorted(str(x) for x in o)
    if isinstance(o, tuple):
        return [str(x) for x in o]
    return str(o)


def write(out_dir, payload):
    out = pathlib.Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "t6_report.json").write_text(json.dumps(payload, indent=1, default=_default), encoding="utf-8")
    L = ["# T6 harness report (software tests only)", "",
         f"Mode: **{payload['mode']}**. Overall: **{payload['overall']}**.", "",
         "Scientific statuses are unchanged and are not affected by this report: B1-F FAIL; coupled-block W/Z DESIGN "
         "INSUFFICIENT; W/Z NOT VALIDATED. The C4 case is a numerical cross-check (not validated).", "",
         "| stage | checks | failed | status |", "|---|---|---|---|"]
    for st in payload["stages"]:
        L.append(f"| {st['stage']} | {st.get('total', '-')} | {st.get('failed', '-')} | {st['status']} |")
    if payload.get("t1t5"):
        L += ["", "## T1-T5", "", "| test | checks | failed |", "|---|---|---|"]
        for st in payload["t1t5"]:
            L.append(f"| {st['stage']} | {st['total']} | {st['failed']} |")
    if payload.get("mutations"):
        L += ["", "## Mutation suite", "", "| id | detected | mechanism | description |", "|---|---|---|---|"]
        for r in payload["mutations"]:
            L.append(f"| {r['id']} | {'yes' if r['detected'] else '**NO**'} | {', '.join(r['mechanism'])} | {r['description']} |")
    if payload.get("c4"):
        c = payload["c4"]
        L += ["", f"## C4 — {c['label']}", "", f"Field outcome: **{c['field_outcome']}** (blocking: {c['blocking']}); "
              f"{c['n_cases']} case/Y combinations."]
    if payload.get("anchors"):
        a = payload["anchors"]
        L += ["", "## Double-precision anchors (diagnostic only)", "",
              f"{sum(r['within'] for r in a)}/{len(a)} within 10·cond₂·ε; max rel. error "
              f"{max(r['rel_err'] for r in a):.2e}."]
    if payload.get("notes"):
        L += ["", "## Notes", ""] + [f"- {n}" for n in payload["notes"]]
    (out / "t6_report.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    return out
