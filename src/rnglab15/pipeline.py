"""Registered measured-data analysis with bounded records and feature memory."""
from __future__ import annotations
from collections import Counter
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import numpy as np
from .assessment import Health, Moments, Timing, screen
from .conditioning import treatments
from .ingest import waveforms, sensors, events, seconds_per_tick
from .integrity import canonical, digest, inventory, snapshot, write_json

SENSOR_COLUMNS = ["ax", "ay", "az", "gx", "gy", "gz", "pressure_raw", "temp_raw", "humidity_raw", "ir_pulses", "ir_age_us"]
MAX_CHANNELS = 16
SCREEN_LIMIT = 2_000_000


class SymbolSink:
    def __init__(self, path):
        self.path = path
        self.file = path.open("xb")
        self.health = Health()
        self.n = 0

    def add(self, bit):
        self.file.write(bytes([bit]))
        self.health.step(bit)
        self.n += 1

    def close(self):
        self.file.close()

    def report(self):
        # Explicit prefix bound; full health test state covers the complete stream.
        with self.path.open("rb") as f:
            sample = f.read(SCREEN_LIMIT)
        return {"health": self.health.report(), "screen": screen(sample) if sample else None,
                "screen_prefix_limit_symbols": SCREEN_LIMIT, "screen_truncated": self.n > SCREEN_LIMIT,
                "symbol_file_sha256": digest(self.path), "symbols_total": self.n}


def describe_waveform(path, spec, out):
    n_samples = spec["samples"]
    if not 4 <= n_samples <= 16384:
        raise ValueError("waveform size must be 4..16384")
    pre = spec.get("baseline_window", [0, min(100,n_samples//4)])
    post = spec.get("post_window", [160,220] if n_samples >= 220 else [n_samples//2,3*n_samples//4])
    tail = spec.get("tail_window", [max(0,n_samples-100),n_samples] if n_samples >= 220 else [3*n_samples//4,n_samples])
    for bounds in (pre,post,tail):
        if len(bounds)!=2 or not 0 <= bounds[0] < bounds[1] <= n_samples or bounds[1]-bounds[0]<2:
            raise ValueError("invalid explicit feature windows")
    scale = seconds_per_tick(spec)
    channels = {}
    sinks = {}
    results = {}
    start = time.monotonic()
    try:
        with (out/"events.csv").open("x", newline="") as ef:
            writer = csv.writer(ef)
            writer.writerow(["row", "channel", "timestamp_raw", "gap_raw", "baseline_median", "baseline_sd",
                             "baseline_diff_mad", "peak", "ceiling_samples", "floor_samples", "crossing_index",
                             "post_median_minus_baseline", "tail_median_minus_baseline"])
            for event in waveforms(path, samples=n_samples):
                key = str(event.channel)
                if key not in channels:
                    if len(channels) == MAX_CHANNELS:
                        raise ValueError("too many channels for bounded analysis")
                    channels[key] = {"timing": Timing(), "rows": 0, "clipped_events": 0, "ceiling_codes": 0,
                                     "floor_codes": 0, "crossings": Counter(), "baseline": Moments(),
                                     "sum": np.zeros(n_samples), "sumsq": np.zeros(n_samples),
                                     "adjacent_waveform_duplicates": 0, "previous_digest": None,
                                     "pre_x": 0.0, "pre_y": 0.0, "pre_xx": 0.0, "pre_yy": 0.0,
                                     "pre_xy": 0.0, "pre_pairs": 0, "blocks": [], "block": Moments()}
                    for label in ("baseline_sample0_bit0", "timestamp_bit0", "gap_bit0", "gap_divisor_bit0", "peak_bit0"):
                        sinks[(key, label)] = SymbolSink(out/f"ch{key}-{label}.symbols")
                c = channels[key]
                x = np.asarray(event.samples, dtype=np.float64)
                prex = x[pre[0]:pre[1]]
                baseline = float(np.median(prex))
                diff = np.diff(prex)
                mad = float(1.4826*np.median(np.abs(diff-np.median(diff)))/np.sqrt(2))
                clipped = int(np.count_nonzero(x == 65535))
                floor = int(np.count_nonzero(x == 0))
                crossing = np.flatnonzero(x > spec.get("analysis_threshold", 10000))
                cross = int(crossing[0]) if len(crossing) else -1
                gap = c["timing"].add(event.timestamp)
                peak = int(x.max())
                writer.writerow([event.row, key, event.timestamp, gap, baseline, float(prex.std(ddof=1)),
                                 mad, peak, clipped, floor, cross, float(np.median(x[post[0]:post[1]])-baseline),
                                 float(np.median(x[tail[0]:tail[1]])-baseline)])
                c["rows"] += 1; c["clipped_events"] += clipped > 0
                c["ceiling_codes"] += clipped; c["floor_codes"] += floor
                c["crossings"][cross] += 1; c["baseline"].add(baseline); c["block"].add(baseline)
                if c["block"].n == 512:
                    c["blocks"].append(c["block"].report())
                    c["block"] = Moments()
                    if len(c["blocks"]) > 4096:
                        raise ValueError("block-report budget exceeded; split acquisition")
                c["sum"] += x; c["sumsq"] += x*x
                a, b = prex[:-1]-baseline, prex[1:]-baseline
                c["pre_pairs"] += len(a)
                c["pre_x"] += float(a.sum()); c["pre_y"] += float(b.sum())
                c["pre_xx"] += float(a@a); c["pre_yy"] += float(b@b); c["pre_xy"] += float(a@b)
                wh = hashlib.sha256(x.astype("<u2").tobytes()).digest()
                c["adjacent_waveform_duplicates"] += wh == c["previous_digest"]
                c["previous_digest"] = wh
                for label, bit in (("baseline_sample0_bit0", event.samples[0]&1),
                                   ("timestamp_bit0", event.timestamp&1), ("peak_bit0", peak&1)):
                    sinks[(key, label)].add(bit)
                if gap is not None and gap > 0:
                    sinks[(key, "gap_bit0")].add(gap&1)
                    divisor = spec.get("exploratory_gap_divisor", 1)
                    if not isinstance(divisor, int) or divisor < 1 or gap % divisor:
                        raise ValueError("declared exploratory divisor does not divide a positive gap")
                    sinks[(key, "gap_divisor_bit0")].add((gap//divisor)&1)
        for key, c in channels.items():
            n = c["rows"]
            if c["block"].n:
                c["blocks"].append(c["block"].report())
            pn = c["pre_pairs"]
            den = (pn*c["pre_xx"]-c["pre_x"]**2)*(pn*c["pre_yy"]-c["pre_y"]**2)
            corr = (pn*c["pre_xy"]-c["pre_x"]*c["pre_y"])/np.sqrt(den) if den > 0 else None
            with (out/f"ch{key}-sample_profile.csv").open("x", newline="") as f:
                w = csv.writer(f); w.writerow(["sample_index", "mean", "sd"])
                mean = c["sum"]/n
                sd = np.sqrt(np.maximum(0, (c["sumsq"]-n*mean**2)/max(1,n-1)))
                w.writerows(zip(range(n_samples), mean.tolist(), sd.tolist()))
            results[key] = {"records": n, "waveform_samples": n_samples,
                            "clipped_events": c["clipped_events"], "clipped_event_fraction": c["clipped_events"]/n,
                            "ceiling_codes": c["ceiling_codes"], "floor_codes": c["floor_codes"],
                            "first_crossing_all_samples": dict(c["crossings"]),
                            "analysis_threshold": spec.get("analysis_threshold",10000),
                            "baseline_medians": c["baseline"].report(), "baseline_512_event_blocks": c["blocks"],
                            "pooled_baseline_lag1_correlation": float(corr) if corr is not None else None,
                            "baseline_correlation_method": "event-median centering, no detrending, configured pre-window; descriptive",
                            "adjacent_waveform_duplicates": c["adjacent_waveform_duplicates"],
                            "timing": c["timing"].report(scale), "feature_windows": {"baseline":pre,"post":post,"tail":tail}, "physical_identity": "unverified"}
    finally:
        for sink in sinks.values():
            sink.close()
    return {"kind": "waveform", "channels": results,
            "symbol_streams": {f"ch{k}-{label}": sink.report() for (k,label),sink in sinks.items()},
            "selection": "one first ADC sample per event; full timestamp LSB; positive gap LSB; divisor is exploratory only",
            "seconds": time.monotonic()-start}


def describe_sensors(path, spec, out):
    columns = spec.get("columns", SENSOR_COLUMNS)
    if len(columns) > 64 or len(columns) != len(set(columns)):
        raise ValueError("invalid bounded column list")
    sinks = {key: SymbolSink(out/f"{key}.symbols") for key in columns}
    stats = {key: Moments() for key in columns}
    missing = Counter()
    timing = Timing()
    previous = None
    repeated = 0
    count = 0
    try:
        with (out/"sensor_features.csv").open("x", newline="") as f:
            w = csv.writer(f); w.writerow(["row", "timestamp_raw", "gap_raw", *columns])
            for row, tick, values, original in sensors(path, columns):
                count += 1
                gap = timing.add(tick)
                record = (tick, *values.values())
                repeated += record == previous
                previous = record
                w.writerow([row,tick,gap,*(values[k] for k in columns)])
                for key, value in values.items():
                    if value is None:
                        missing[key] += 1
                    else:
                        stats[key].add(value)
                        sinks[key].add(value&1)
    finally:
        for sink in sinks.values():
            sink.close()
    return {"kind": "sensor_csv", "records": count, "timing": timing.report(seconds_per_tick(spec)),
            "adjacent_duplicate_records": repeated,
            "columns": {key: {"values": stats[key].report(), "missing": missing[key],
                              "symbols": sinks[key].report()} for key in columns},
            "provenance": spec.get("provenance"),
            "interpretation": "sensor register readings; missing values retained; LSB stream omits NA and may bridge resets"}


def describe_events(path, spec, out):
    timings, sinks = {}, {}
    try:
        for row, channel, tick in events(path, spec.get("timestamp_column","timestamp"), spec.get("channel_column","channel")):
            # channel used in filenames only after bounded integer canonicalization
            channel = str(int(channel))
            if channel not in timings:
                if len(timings) == MAX_CHANNELS:
                    raise ValueError("too many channels")
                timings[channel] = Timing()
                sinks[channel] = SymbolSink(out/f"ch{channel}-gap_bit0.symbols")
            gap = timings[channel].add(tick)
            if gap is not None and gap > 0:
                sinks[channel].add(gap&1)
    finally:
        for sink in sinks.values():
            sink.close()
    return {"kind": "event_csv", "channels": {k:t.report(seconds_per_tick(spec)) for k,t in timings.items()},
            "symbol_streams": {k:s.report() for k,s in sinks.items()}}


def code_manifest(repo):
    # Installed wheels have neither a src/ layout nor repository metadata.
    in_checkout = (repo/"src/rnglab15").is_dir() and (repo/"pyproject.toml").is_file()
    root = repo if in_checkout else Path(__file__).resolve().parent.parent
    package_root = root/"src" if in_checkout else root
    paths = sorted((package_root/"rnglab15").glob("*.py")) + sorted((package_root/"rnglab14").rglob("*.py"))
    if in_checkout:
        paths += [repo/"pyproject.toml",repo/"requirements-physical.lock"]
        paths += sorted((repo/"scripts").glob("*physical15*.py"))
        paths += sorted((repo/"tests").glob("test_physical15.py"))
    return {str(p.relative_to(root)):digest(p) for p in paths}


def git_parent(repo):
    if not (repo/"src/rnglab15").is_dir() or not (repo/"pyproject.toml").is_file():
        return None
    try:
        return subprocess.check_output(["git","rev-parse","HEAD"],cwd=repo,text=True,stderr=subprocess.DEVNULL).strip()
    except (FileNotFoundError,subprocess.CalledProcessError):
        return None


def run(config_path, out, vault):
    config_path, out, vault = Path(config_path), Path(out), Path(vault)
    config = json.loads(config_path.read_text())
    if config.get("schema_version") != 1:
        raise ValueError("unsupported run schema")
    out.mkdir(parents=True, exist_ok=False)
    repo = Path(__file__).resolve().parents[2]
    registry = {"schema_version":1, "status":"running", "config":config,
                "config_sha256":digest(config_path), "code":code_manifest(repo),
                "python":sys.version, "platform":platform.platform(),
                "dependencies":{k:importlib.metadata.version(k) for k in ("numpy","scipy","cryptography")},
                "official_nist_non_iid_executable":shutil.which("ea_non_iid"),
                "git_parent":git_parent(repo),
                "sources":[], "validated_entropy_credit_bits":0,
                "qualification":"All measured outputs are replayable research candidates; no secure TRNG claim."}
    write_json(out/"registration.json", registry)
    try:
        for index, spec in enumerate(config["sources"]):
            if spec.get("privacy_review") != "non_sensitive_physical_measurement":
                raise ValueError("explicit non-sensitive physical-data declaration required")
            source = Path(spec["path"]).resolve(strict=True)
            raw, sha = snapshot(source, vault)
            source_out = out/f"source-{index:02d}"
            source_out.mkdir()
            if spec["kind"] == "waveform":
                report = describe_waveform(raw,spec,source_out)
            elif spec["kind"] == "sensor_csv":
                report = describe_sensors(raw,spec,source_out)
            elif spec["kind"] == "event_csv":
                report = describe_events(raw,spec,source_out)
            elif spec["kind"] == "rng_bytes":
                # Exact hardware identity/entropy status belongs in metadata, never inferred from bytes.
                report = {"kind":"rng_bytes", "bytes":raw.stat().st_size, "provenance":spec.get("provenance"),
                          "hardware_entropy_status":"unverified unless documented externally"}
            else:
                raise ValueError("unknown explicit adapter kind")
            report["treatments"] = {}
            for symbol in sorted(source_out.glob("*.symbols")):
                report["treatments"][symbol.stem] = treatments(symbol,source_out/(symbol.stem+"-treatments"))
            if digest(raw) != sha or digest(source) != sha:
                raise ValueError("raw evidence changed during analysis")
            write_json(source_out/"report.json", report)
            registry["sources"].append({"original":str(source), "raw_sha256":sha,
                                        "vault_object":str(raw), "output":source_out.name})
        registry["status"] = "complete"
    except Exception as exc:
        registry["status"] = "failed"
        registry["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        registry["output_hashes"] = {str(p.relative_to(out)):digest(p) for p in sorted(out.rglob("*")) if p.is_file()}
        write_json(out/"registry.json",registry)
    return registry
