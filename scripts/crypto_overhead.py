"""Separate observational microbenchmark; no claim of production overhead."""
import argparse
import json
import statistics
import tempfile
import time
from pathlib import Path

from rnglab_crypto.sources import ReplaySource
from rnglab_crypto.monitoring import source_scores
from rnglab_crypto.primitives.aead import GuardedAEAD, NonceLedger, cipher
from rnglab_crypto.reporting import environment


def measure():
    src = ReplaySource(bytes(32))
    raw = [src.read_bytes(32) for _ in range(256)]
    rows = []
    # Fixed alternating order limits a simple always-first warmup confounder.
    for rep in range(12):
        order = ["aead", "guarded_aead", "source_monitor", "sqlite_counter"]
        if rep % 2:
            order.reverse()
        for name in order:
            with tempfile.TemporaryDirectory(prefix="blue-bench-") as tmp:
                a = cipher("AES-256-GCM", bytes(32))
                guarded = GuardedAEAD("AES-256-GCM", bytes(32))
                ledger = NonceLedger(Path(tmp)/"ledger.db", bytes(32))
                start = time.perf_counter_ns()
                if name == "source_monitor":
                    source_scores(raw)
                else:
                    for op in range(256):
                        if name == "sqlite_counter":
                            ledger.allocate()
                        elif name == "guarded_aead":
                            guarded.encrypt(op.to_bytes(12, "big"), bytes(16))
                        else:
                            a.encrypt(op.to_bytes(12, "big"), bytes(16), None)
                elapsed = time.perf_counter_ns()-start
                ledger.close()
                rows.append({"replicate": rep, "operation": name, "nanoseconds_per_operation": elapsed/256})
    summary = {name: {"median_ns_per_operation": statistics.median(r["nanoseconds_per_operation"] for r in rows if r["operation"] == name),
                      "range_ns_per_operation": [min(r["nanoseconds_per_operation"] for r in rows if r["operation"] == name),
                                                 max(r["nanoseconds_per_operation"] for r in rows if r["operation"] == name)]}
               for name in order}
    return {"environment": environment(), "rows": rows, "summary": summary,
            "qualification": "12 repeated timings on one host, 256 operations each; startup excluded, no independent-machine CI; monitor uses batch vectorization"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    with args.output.open("x") as handle:
        json.dump(measure(), handle, indent=2)
        handle.write("\n")
