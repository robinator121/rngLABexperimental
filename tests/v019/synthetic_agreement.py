"""Agreement of the Stage D pipeline with the independent recount on a synthetic campaign.

Runs every stage of the pipeline once on the 207 synthetic runs of tests/v019/_pilot_fixture.py
(audit-only namespace, no sampler, no pilot or smoke draw), recounts the campaign with
research/v019-audit and compares every result file. Writes docs/v020/evidence/synthetic_agreement.json.

Two campaigns are compared:
  frozen          the pipeline as it is: the oracle follows the frozen rule;
  hypothetical_A  the same campaign with the switch of proposed Amendment 01, option A, forced in
                  memory. No amendment record is written; this shows what the comparison would give.

Needs the measurement layer: RNGLAB_DPILOT_MEASURE=<path to measure.py>, or the file in the tree.

    PYTHONPATH=src python tests/v019/synthetic_agreement.py
"""
from __future__ import annotations

import importlib.util
import json
import platform
import sys
import tempfile
from collections import Counter
from pathlib import Path

import numpy as np
import scipy

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("_pilot_fixture", HERE / "_pilot_fixture.py")
fx = importlib.util.module_from_spec(_spec)
sys.modules["_pilot_fixture"] = fx
_spec.loader.exec_module(fx)
au, pc, pg, pa, pf = fx.au, fx.pc, fx.pg, fx.pa, fx.pf
vf = au.load("v019_audit_verify", au.REPO / "research/v019-audit/verify.py")
OUTPUT = au.REPO / "docs/v020/evidence/synthetic_agreement.json"


def campaign(root: Path, rule: str | None) -> dict:
    ctx = fx.context(fx.TEST_PLAN, root)
    generated = pg.generate(ctx, engine_spec=fx.ENGINE_SPEC, workers=0)
    pa.calibrate(ctx)
    original = pf.oracle_rule
    try:
        if rule is not None:
            pf.oracle_rule = lambda: rule
        pf.fit(ctx)
    finally:
        pf.oracle_rule = original
    pa.issue(ctx)
    scores = pf.evaluate(ctx)
    pa.healthy_summary(ctx)
    pa.fault_summary(ctx)
    pa.implementation_checks(ctx, pg.replay_check(ctx, fx.SyntheticEngine(ctx), {"healthy": [0], "fault/DUP": [1]}))
    pa.criteria(ctx)
    plan = {c: (f, n) for c, f, n in ctx.plan.conditions}
    views = {c: [vf.load_view(ctx.root, c, i) for i in range(n)] for c, (_, n) in plan.items()}
    saved = vf.saved_models_from_primary(pc.read_json(ctx.result("models.json")))
    independent = vf.recount_campaign(views, ctx.plan.namespace, plan=plan, fitting=ctx.plan.fitting,
                                      evaluation=ctx.plan.evaluation, saved_models=saved)
    comparison = vf.compare_with_primary(ctx.root, independent)
    cells = [cell["oracle_all_evaluation_runs"] for family in scores["families"].values() for cell in family.values()]
    return {"oracle_rule": scores["oracle_rule"], "runs": generated["counts"], "recount_status": independent["status"],
            "fields_compared": comparison["compared"], "files_absent": comparison["files_absent"],
            "discrepancies": len(comparison["discrepancies"]),
            "discrepancies_by_field": dict(sorted(Counter(d["field"] for d in comparison["discrepancies"]).items())),
            "discrepancy_list": comparison["discrepancies"],
            "oracle_cells": len(cells), "oracle_cells_with_an_unavailable_run": sum(1 for c in cells if c["unavailable"]),
            "oracle_runs_unavailable": sum(c["unavailable"] for c in cells),
            "oracle_runs_lower_term_replaced": sum(c.get("lower_term_replaced", 0) for c in cells),
            "oracle_cells_scored": sum(1 for c in cells if c["brier"] is not None)}


def main() -> int:
    if au.measure_path() is None:
        print("the measurement layer is not available: set RNGLAB_DPILOT_MEASURE")
        return 2
    out = {"schema": "rnglab-v020-synthetic-agreement-v1",
           "what": "synthetic fixtures in the audit-only namespace; not Monte Carlo output and not a pilot result",
           "namespace": fx.NAMESPACE, "plan": {c: n for c, _, n in fx.TEST_PLAN.conditions},
           "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}}
    for name, rule in (("frozen", None), ("hypothetical_A", "A01")):
        with tempfile.TemporaryDirectory() as directory:
            out[name] = campaign(Path(directory), rule)
        row = out[name]
        print(f"{name}: rule {row['oracle_rule']}, {row['fields_compared']} fields compared, {row['discrepancies']} differ "
              f"{row['discrepancies_by_field']}; oracle unavailable in {row['oracle_runs_unavailable']} run-lengths, "
              f"{row['oracle_cells_scored']} of {row['oracle_cells']} oracle cells scored")
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote", OUTPUT.relative_to(au.REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
