import csv
import hashlib
import importlib.metadata
import json
import platform
from collections import defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import t
from rnglab.providers import make_provider, SyntheticSource, bits_to_uniform
from rnglab.degradation import FaultSpec
from rnglab.simulations import WORKLOADS
from rnglab.diagnostics import diagnose


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def validate(config):
    for field in ("samples", "window", "diagnostic_bits", "trials"):
        if not isinstance(config[field], int) or config[field] < (2 if field in {"trials", "diagnostic_bits"} else 1):
            raise ValueError(f"Invalid {field}")
    if not isinstance(config["seed"], int) or config["seed"] < 0:
        raise ValueError("seed must be a nonnegative integer")
    if config["checkpoints"] != sorted(set(config["checkpoints"])) or not config["checkpoints"]:
        raise ValueError("checkpoints must be nonempty, sorted and unique")
    if any(not isinstance(n, int) or n < 1 or n > config["samples"] for n in config["checkpoints"]):
        raise ValueError("checkpoints must be valid sample counts")
    if config["checkpoints"][-1] != config["samples"]:
        raise ValueError("last checkpoint must equal samples")
    if not config["workloads"] or any(w not in WORKLOADS for w in config["workloads"]):
        raise ValueError("unknown/empty workload list")
    scenarios = config["scenarios"]
    if not scenarios or len({c["name"] for c in scenarios}) != len(scenarios):
        raise ValueError("scenario names must be nonempty and unique")
    for scenario in scenarios:
        if scenario["provider"] not in {"synthetic", "prng", "csprng", "hybrid"}:
            raise ValueError("unknown provider")
        if scenario["provider"] == "hybrid" and (not isinstance(scenario.get("reseed_interval"), int) or scenario["reseed_interval"] < 1):
            raise ValueError("hybrid needs positive integer reseed_interval")
        for field in ("faults", "initial_faults"):
            for fault in scenario.get(field, []):
                FaultSpec(**fault)
        if scenario["provider"] in {"prng", "csprng"} and scenario.get("faults"):
            raise ValueError("Use initial_faults for CSPRNG; PRNG faults require a source adapter")
        if scenario["provider"] in {"synthetic", "prng"} and scenario.get("initial_faults"):
            raise ValueError("initial_faults only applies to CSPRNG/hybrid")


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[(row["scenario"], row["workload"], row["n"])].append(row)
    summary = []
    for (scenario, workload, n), data in groups.items():
        values = np.array([r["estimate"] for r in data])
        target = data[0]["target"]
        errors = values - target
        absolute = np.abs(errors)
        critical = t.ppf(0.975, len(data)-1)
        half = float(critical * values.std(ddof=1)/np.sqrt(len(data)))
        summary.append({"scenario": scenario, "workload": workload, "n": n, "trials": len(data),
                        "target": target, "mean_estimate": float(values.mean()),
                        "bias": float(errors.mean()), "variance": float(values.var(ddof=1)),
                        "sd": float(values.std(ddof=1)), "mean_ci95_low": float(values.mean()-half),
                        "mean_ci95_high": float(values.mean()+half),
                        "mean_abs_error": float(absolute.mean()),
                        "mae_ci95_half": float(critical*absolute.std(ddof=1)/np.sqrt(len(data))),
                        "rmse": float(np.sqrt(np.mean(errors**2)))})
    return summary


def run(config, output):
    validate(config)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "config.json").write_text(json.dumps(config, indent=2)+"\n")
    rows, diagnostics, windows, histograms, events, fingerprints = [], [], [], [], [], []
    for scenario in config["scenarios"]:
        print(f"Running {scenario['name']}", flush=True)
        for workload_name in config["workloads"]:
            work = WORKLOADS[workload_name]
            # Stable identity independent of the order of workloads or scenarios.
            wid = int.from_bytes(hashlib.sha256(workload_name.encode()).digest()[:4], "big")
            for trial in range(config["trials"]):
                seed = [config["seed"], trial, wid]
                provider = make_provider(scenario, seed)
                identity = {"scenario": scenario["name"], "workload": workload_name, "trial": trial}
                if scenario["provider"] in {"csprng", "hybrid"}:
                    fingerprints.append({**identity, "initial_key_sha256": hashlib.sha256(provider.key).hexdigest()})
                total, start = 0.0, 0
                checkpoints = iter(config["checkpoints"])
                checkpoint = next(checkpoints, None)
                while start < config["samples"]:
                    count = min(config["window"], config["samples"]-start)
                    u, bits = provider.sample(count, work.dimension)
                    contributions = work.contributions(u)
                    cumulative = np.cumsum(contributions, dtype=float)
                    end = start + count
                    while checkpoint is not None and checkpoint <= end:
                        estimate = float((total + cumulative[checkpoint-start-1])/checkpoint)
                        rows.append({**identity, "n": checkpoint, "estimate": estimate,
                                     "target": work.target, "absolute_error": abs(estimate-work.target)})
                        checkpoint = next(checkpoints, None)
                    total += float(cumulative[-1])
                    windows.append({**identity, "start_sample": start, "end_sample": end,
                                    "window_estimate": float(contributions.mean()), "target": work.target})
                    # Diagnose a bounded PREFIX of each output window. Histograms use the full window.
                    stats, hist = diagnose(bits[:config["diagnostic_bits"]], u.ravel())
                    diagnostics.append({**identity, "stream": "output", "start_sample": start, "end_sample": end, **stats})
                    for bin_index, value in enumerate(hist):
                        histograms.append({**identity, "start_sample": start, "end_sample": end,
                                           "bin": bin_index, "count": int(value)})
                    start = end
                events.extend({**identity, **event} for event in provider.events)
                # Separate audit streams: never advance the stream feeding the simulation.
                audits = []
                if scenario["provider"] == "synthetic":
                    audits.append(("source_audit", scenario.get("faults", [])))
                if scenario["provider"] in {"csprng", "hybrid"}:
                    audits.append(("initial_source_audit", scenario.get("initial_faults", [])))
                if scenario["provider"] == "hybrid":
                    audits.append(("reseed_source_audit", scenario.get("faults", [])))
                for stream, faults in audits:
                    raw = SyntheticSource([*seed, 9182], faults).read_bits(config["diagnostic_bits"])
                    length = len(raw)//32*32
                    if length:
                        stats, _ = diagnose(raw, bits_to_uniform(raw[:length], length//32, 1).ravel())
                        diagnostics.append({**identity, "stream": stream, "start_sample": None, "end_sample": None, **stats})
    summary = summarize(rows)
    key_groups = defaultdict(set)
    for row in fingerprints:
        key_groups[(row["scenario"], row["workload"])].add(row["initial_key_sha256"])
    for row in summary:
        keys = key_groups.get((row["scenario"], row["workload"]))
        row["distinct_initial_keys"] = len(keys) if keys is not None else None
        row["uncertainty_note"] = "pointwise trial-mean interval"
        if keys is not None and len(keys) == 1:
            row["uncertainty_note"] = "one initial key; no independent initial-key replication"
            row["mean_ci95_low"] = row["mean_ci95_high"] = row["mae_ci95_half"] = None
    for name, data in [("trials", rows), ("summary", summary), ("diagnostics", diagnostics),
                       ("windows", windows), ("histograms", histograms), ("reseed_events", events), ("key_fingerprints", fingerprints)]:
        write_csv(output / f"{name}.csv", data)
    source_root = Path(__file__).parents[1]
    source_hash = hashlib.sha256()
    for file in sorted(source_root.rglob("*.py")):
        source_hash.update(str(file.relative_to(source_root)).encode())
        source_hash.update(file.read_bytes())
    metadata = {"python": platform.python_version(), "platform": platform.platform(),
                "libraries": {name: importlib.metadata.version(name) for name in ["numpy", "scipy", "matplotlib", "cryptography"]},
                "source_sha256": source_hash.hexdigest(),
                "seed_scheme": "[master_seed, trial_index, int.from_bytes(sha256(workload)[:4], 'big')]",
                "config_pairing": "Same trial seed across configurations; trials use distinct seeds.",
                "entropy_model": "Deterministic PCG64 Bernoulli model; no real physical entropy or security claim.",
                "uncertainty": "Student-t 95% intervals across trials; pointwise, unadjusted; degenerate seeds invalidate effective replication.",
                "diagnostic_null": "IID Bernoulli(1/2) bits / IID uniform samples; not valid entropy or security certificates."}
    (output / "metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    return summary, diagnostics, windows, histograms
