"""Timing of reading stored runs: what the analysis stages and the independent verification cost per run.

Generates twelve full-length runs through the pipeline's run_one in the audit-only namespace (six
worker processes, the pilot's mix of conditions; generation is not timed here), then times, per run:

  pipeline_classify       pilot_generate.classify with the trace hash (resume, run manifest, completeness gate)
  pipeline_load           pilot_analysis.load_run (every analysis stage; verifies file and array hashes)
  audit_classify          research/v019-audit/integrity.classify_run with the trace hash
  audit_load_view         research/v019-audit/verify.load_view (the independent recount's reader)
  audit_verify_shallow    verify.verify_run without regeneration
  audit_verify_deep       verify.verify_run with regeneration (the 25 sampled runs)

The store is deleted at the end. The runs are timing material, not pilot, smoke or campaign results.
Needs the measurement layer: RNGLAB_DPILOT_MEASURE=<path to measure.py>, or the file in the tree.

    PYTHONPATH=src python tests/v019/pipeline_read_timing.py      writes docs/v020/evidence/pipeline_read_timing.json
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import statistics
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("pipeline_benchmark", HERE / "pipeline_benchmark.py")
pb = importlib.util.module_from_spec(_spec)
sys.modules["pipeline_benchmark"] = pb
_spec.loader.exec_module(pb)   # sets the thread variables before NumPy is imported
fx, au, pc, pg = pb.fx, pb.au, pb.pc, pb.pg
pa = fx.pa
vf = au.load("v019_audit_verify", au.REPO / "research/v019-audit/verify.py")
integrity = au.load("v019_audit_integrity", au.REPO / "research/v019-audit/integrity.py")

OUTPUT = au.REPO / "docs/v020/evidence/pipeline_read_timing.json"
ROOT = au.REPO / ".devtmp/pipeline_read_timing"
WORKERS, JOBS = 6, 12
# how often each operation is applied to a run in a pilot: stages that read every run of the sets they use
PILOT_USES = {"pipeline_classify": 3 * 5999,          # resume scan, run manifest, completeness gate
              "pipeline_load": 999 + (999 + 1800) + 5999 + 1200 + 2000 + 3000 + 5999,
              # calibrate; fit (calibration + fitting); issue; evaluate; healthy; faults; checks
              "audit_classify": 5999, "audit_load_view": 5999, "audit_verify_shallow": 5999 - 25, "audit_verify_deep": 25}


def timed(function, jobs) -> list[float]:
    out = []
    for job in jobs:
        started = time.perf_counter()
        function(*job)
        out.append(time.perf_counter() - started)
    return out


def main() -> int:
    if au.measure_path() is None:
        print("the measurement layer is not available: set RNGLAB_DPILOT_MEASURE")
        return 2
    shutil.rmtree(ROOT, ignore_errors=True)
    ROOT.mkdir(parents=True)
    context = fx.context(pb.PLAN, ROOT, need_measure=False)
    jobs = pb.jobs_for(JOBS)
    pc.atomic_write(ROOT / ".pending.json", pc.canonical_json([list(j) for j in jobs]), overwrite=True)
    processes = [subprocess.Popen([sys.executable, str(HERE / "pipeline_benchmark.py"), "--worker", f"{k}/{WORKERS}",
                                   "--root", str(ROOT)], env=dict(os.environ)) for k in range(WORKERS)]
    if any(p.wait() for p in processes):
        print("a worker failed")
        return 1
    ids = [(condition, replica) for condition, _, replica in jobs]
    namespace = context.plan.namespace

    def shallow(condition, replica):
        return vf.verify_run(*vf.load_stored(ROOT, condition, replica), namespace, deep=False)

    def deep(condition, replica):
        return vf.verify_run(*vf.load_stored(ROOT, condition, replica), namespace, deep=True)

    operations = {
        "pipeline_classify": lambda c, r: pg.classify(context, c, r, verify_hash=True),
        "pipeline_load": lambda c, r: pa.load_run(context, c, r),
        "audit_classify": lambda c, r: integrity.classify_run(ROOT, c, r, verify_hash=True),
        "audit_load_view": lambda c, r: vf.load_view(ROOT, c, r),
        "audit_verify_shallow": shallow,
        "audit_verify_deep": deep,
    }
    out = {"schema": "rnglab-v020-pipeline-read-timing-v1",
           "what": "timing of reading stored runs in the audit-only namespace; the runs are deleted and are not pilot, "
                   "smoke or campaign results",
           "namespace": namespace, "runs": len(ids), "single_process": True, "operations": {}}
    issues = 0
    for name, function in operations.items():
        targets = ids[:6] if name == "audit_verify_deep" else ids
        if name.startswith("audit_verify"):
            reports = []
            seconds = timed(lambda c, r: reports.append(function(c, r)), targets)
            issues += sum(len(report["issues"]) for report in reports)
        else:
            seconds = timed(function, targets)
        mean = statistics.fmean(seconds)
        out["operations"][name] = {"runs_timed": len(targets), "mean_seconds": mean, "max_seconds": max(seconds),
                                   "uses_in_a_pilot": PILOT_USES[name], "pilot_minutes_single_process": PILOT_USES[name] * mean / 60}
        print(f"{name}: {mean * 1000:.0f} ms per run, pilot {PILOT_USES[name] * mean / 60:.1f} min")
    out["independent_verification_issues"] = issues
    out["pilot_minutes_single_process"] = {
        "pipeline_analysis_stages": out["operations"]["pipeline_load"]["pilot_minutes_single_process"],
        "pipeline_scans": out["operations"]["pipeline_classify"]["pilot_minutes_single_process"],
        "independent_verification": sum(out["operations"][name]["pilot_minutes_single_process"]
                                        for name in ("audit_load_view", "audit_verify_shallow", "audit_verify_deep"))}
    out["note"] = "reading only: the arithmetic of the stages (fits, bootstrap, summaries) is not included"
    shutil.rmtree(ROOT, ignore_errors=True)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("independent verification issues on these runs:", issues)
    print("wrote", OUTPUT.relative_to(au.REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
