"""Stage D pilot: command line. protocol.json governs; nothing is generated without the author's authorization.

    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py <command> [--smoke]

Order of a pilot (freeze procedure, steps 5 to 7):

    status                 state of the output tree; changes nothing
    manifest               write pre_generation_manifest.json (needs --interrupted-run-policy)
                           -> the author reviews it, pushes it and writes authorization.json
    generate               run every planned run that has no record; resumable
    replay                 I-4: regenerate the 25 selected runs and compare
    calibrate              thresholds, application constants, sigma_eff
    fit                    forecasters and oracle parameters, from fitting runs only
    issue                  alarm, eligibility and issuance endpoint of every run
    evaluate               predictions and scores; opens the evaluation runs once (separate process)
    summaries              healthy.json, faults.json
    checks                 implementation_checks.json (I-1 to I-4)
    criteria               criteria.json

`--smoke` selects the smoke namespace and the directory smoke/. Only status, manifest, generate
and replay accept it: no threshold, model or criterion is computed from a smoke.

authorization.json, written by the author next to the manifest:

    {"authorized_by": "<name>", "authorized_on": "YYYY-MM-DD", "namespace": "<namespace>",
     "pre_generation_manifest_sha256": "<sha256 of pre_generation_manifest.json>"}
"""
from __future__ import annotations

import os

# protocol.json storage.parallelism: set before NumPy or SciPy is imported, in this process and its workers
for _name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
              "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_name] = "1"

import argparse  # noqa: E402
import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))
import pilot_common as pc  # noqa: E402

SMOKE_COMMANDS = ("status", "manifest", "generate", "_worker", "replay", "smoke-timing")


def context_for(smoke: bool) -> pc.Context:
    plan, root = (pc.SMOKE_PLAN, HERE / "smoke") if smoke else (pc.FROZEN_PLAN, HERE)
    return pc.production_context(plan, root)


def smoke_timing(context: pc.Context) -> dict:
    """smoke/smoke_timing.json: run plan, timings, failures and cap usage of the stored smoke runs."""
    import pilot_generate as pg
    rows = []
    for condition, family, replica in context.plan.runs():
        state = pg.classify(context, condition, replica)
        row = {"condition": condition, "replica": replica, "state": state["state"], "status": state.get("status")}
        if state["state"] in ("COMPLETE", "FAILED_RECORD"):
            record = pc.read_json(context.record_path(condition, replica))
            row.update({key: record.get(key) for key in pg.TIMING_FIELDS})
            row["windows_cut"] = sum(1 for c in record.get("conventional") or [] if c.get("was_cut"))
            row["checks"] = {name: (record.get(name) or {}).get("passed") for name in ("I-1", "I-2", "I-3", "I-3-step")}
            row["failures"] = len(record.get("failures") or [])
        rows.append(row)
    body = {"schema": "rnglab-v019-d-pilot-smoke-timing-v1", "namespace": context.plan.namespace,
            "excluded_from_all_results": True, "runs": rows}
    pc.atomic_write(context.root / "smoke_timing.json", pc.canonical_json(body))
    return body


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", choices=("status", "manifest", "generate", "_worker", "replay", "smoke-timing", "calibrate",
                                            "fit", "issue", "evaluate", "summaries", "checks", "criteria"))
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--interrupted-run-policy", choices=("regenerate_same_id", "record_failure"), default=None,
                        help="manifest: what a restart does with a run that was started and left no record (no default)")
    parser.add_argument("--tdk-deadline", default=None, help="manifest: TDK submission deadline, YYYY-MM-DD")
    parser.add_argument("--workers", type=int, default=pc.PROTOCOL_WORKERS, help="manifest: worker processes (protocol: 8)")
    parser.add_argument("--shard", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--engineering-complete-on", default=None, help="criteria: date the engineering was complete")
    parser.add_argument("--results-available-on", default=None, help="criteria: date the pilot results were available")
    args = parser.parse_args(argv)
    if args.smoke and args.command not in SMOKE_COMMANDS:
        parser.error(f"{args.command} is not available for the smoke")
    try:
        context = context_for(args.smoke)
        import pilot_analysis as pa
        import pilot_forecast as pf
        import pilot_generate as pg
        if args.command == "status":
            states = pg.scan(context, verify_hash=False)
            counts = {name: sum(1 for s in states if s["state"] == name) for name in pg.STATES}
            print(json.dumps({"namespace": context.plan.namespace, "root": str(context.root), "planned": len(states),
                              "counts": counts,
                              "pre_generation_manifest": (context.root / "pre_generation_manifest.json").is_file(),
                              "authorization": (context.root / "authorization.json").is_file()}, indent=2))
        elif args.command == "manifest":
            if args.interrupted_run_policy is None:
                parser.error("--interrupted-run-policy is required: it is the author's decision and has no default")
            manifest = pg.build_pre_generation_manifest(context, interrupted_run_policy=args.interrupted_run_policy,
                                                        tdk_deadline=args.tdk_deadline, workers=args.workers)
            digest = pc.atomic_write(context.root / "pre_generation_manifest.json", pc.canonical_json(manifest))
            print(f"wrote pre_generation_manifest.json, sha256 {digest}")
            print("No draw is permitted until the author has reviewed and pushed it and written authorization.json.")
            if manifest["tdk_deadline"] is None:
                print("TDK deadline not recorded: stop criterion 5 will be INSUFFICIENT_INFORMATION.")
            if manifest["workers_deviate_from_protocol"]:
                print(f"workers = {manifest['workers']} differs from the protocol's {pc.PROTOCOL_WORKERS}.")
        elif args.command == "generate":
            manifest = pg.require_generation_allowed(context)
            command = [sys.executable, str(HERE / "run_pilot.py"), "_worker"] + (["--smoke"] if args.smoke else [])
            print(json.dumps(pg.generate(context, workers=manifest["workers"], worker_command=command), indent=2))
        elif args.command == "_worker":
            pg.require_generation_allowed(context)
            pending = [tuple(p) for p in pc.read_json(context.root / ".pending.json")]
            pending = [p for p in pending if pg.classify(context, p[0], p[2], verify_hash=False)["state"] == "MISSING"]
            pg.work(context, "production", pg.shard(pending, args.shard))
        elif args.command == "replay":
            pg.require_complete(context)
            result = pg.replay_check(context, pg.PilotEngine(context))
            pc.atomic_write(context.result("replay.json"), pc.canonical_json(result))
            print(json.dumps({"passed": result["passed"], "runs": len(result["runs"])}))
        elif args.command == "smoke-timing":
            smoke_timing(context)
        else:
            pg.require_complete(context)
            if args.command == "calibrate":
                print(pa.calibrate(context)["status"])
            elif args.command == "fit":
                pf.fit(context)
            elif args.command == "issue":
                pa.issue(context)
            elif args.command == "evaluate":
                pf.evaluate(context)
            elif args.command == "summaries":
                pa.healthy_summary(context)
                pa.fault_summary(context)
            elif args.command == "checks":
                replay = pc.read_json(context.result("replay.json")) if context.result("replay.json").is_file() else None
                print(json.dumps({"passed": pa.implementation_checks(context, replay)["passed"]}))
            elif args.command == "criteria":
                out = pa.criteria(context, engineering_complete_on=args.engineering_complete_on,
                                  results_available_on=args.results_available_on)
                print(json.dumps(out["results"], indent=2))
    except pc.PilotError as error:
        print(f"STOPPED: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
