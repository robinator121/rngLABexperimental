"""Stage D per-run measurements and storage; frozen protocol.json governs.

No campaign driver or calibration fitting. measure_run returns (record, arrays),
without writing; generation is timed by its caller, measurement by this module.
Statistical fields use parallel statistic_statuses/statistic_reasons mappings.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import platform
import time
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from scipy.stats import t

from rnglab11.statistics import evaluate_bits
from rnglab12.uncertainty import summarize_series
from rnglab14.models.pimc_harmonic import exact_finite_bead_energy

REPO = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("d_pilot_measure_dlib", Path(__file__).with_name("dlib.py"))
dlib = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dlib)
h1lib = dlib.h1lib
_spec = importlib.util.spec_from_file_location("C_analyze", REPO / "research/v018-h1-detection/analyze.py")
C_analyze = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(C_analyze)
H, LENGTHS = dlib.H, (10240, 20480, 40960)
BATTERY_ENDS, APPLICATION_ENDS = tuple(range(512, H + 1, 256)), tuple(range(1024, H + 1, 256))
P_KEYS = h1lib.CONV_P_KEYS
E_M = exact_finite_bead_energy(16, 2.0)
ARRAY_DTYPES = {"CV": "<f8", "TD": "<f8", "marks": "<i8", "severity": "<f8"}
STAT_FIELDS = (
    "frequency_signed_z ones block_frequency_chi2 block_size blocks unused_bits "
    "runs_residual runs_signed_z V pi prerequisite_failed serial_psi2 serial_delta1 serial_delta2 "
    "pattern_counts_6 pattern_counts_7 pattern_counts_8 pattern_counts_9 "
    "approximate_entropy_phi8 approximate_entropy_phi9 approximate_entropy approximate_entropy_chi2 "
    "cumulative_sums_max_abs longest_run_chi2 longest_run_categories longest_run_block_size "
    "longest_run_blocks longest_run_unused_bits longest_run_table_key longest_run_edges "
    "lag1_r lag1_test_z lag1_feature duplicate_count duplicate_rate duplicate_lambda"
).split()


def _finite(value: object) -> bool:
    return value is not None and bool(np.isfinite(value))


def _json_safe(value: object) -> object:
    """storage.json: preserve counts, replace unavailable nonfinite numbers by null."""
    if isinstance(value, dict):
        return {k: _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_json_safe(v) for v in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def _exception(error: Exception) -> dict:
    return {"class": type(error).__name__, "message": str(error)}


def _log(out: dict, name: str, value: object, reason: str | None = None,
         status: str = "NOT_APPLICABLE") -> None:
    good = value is not None and bool(np.isfinite(value).all())
    out[name] = _json_safe(value) if good else None
    out["statistic_statuses"][name] = "AVAILABLE" if good and reason is None else status
    out["statistic_reasons"][name] = reason if reason else (None if good else "NONFINITE_OR_MISSING")


def statistic(p_values: Iterable[float | None]) -> tuple[float, int]:
    """operational_alarm.statistic: finite in-range p only, floor, empty S=0."""
    ps = [float(p) for p in p_values if _finite(p) and 0 <= p <= 1]
    return (-math.log10(max(min(ps), h1lib.P_FLOOR)) if ps else 0.0), len(ps)


def battery_checkpoint(emitted: bytes, marks: np.ndarray, end: int) -> dict:
    """operational_alarm.bit_cap; signal_schema.checkpoint_metadata/battery_statistics.

    Shared intermediates execute the same Stage C calls as conventional_window;
    the word call always sees the full window. Exceptions are retained per test.
    """
    out = {"end": int(end), "statistic_statuses": {}, "statistic_reasons": {}, "failures": []}
    for key in ("byte_start byte_end original_bytes original_bits bit_used_bytes bit_used_bits "
                "word_used_bytes word_count pair_count was_cut").split():
        out[key] = None
    for key in STAT_FIELDS:
        _log(out, key, None, "INCOMPLETE_CHECKPOINT", "ERROR")
    raw = dict.fromkeys(P_KEYS)
    annotations: dict[str, tuple[str, str | None]] = {}

    def calculate(keys: Iterable[str], call: Callable) -> object:
        try:
            return call()
        except Exception as error:
            out["failures"].append({"phase": "battery", "end": int(end), "tests": list(keys),
                                    "exception": _exception(error)})
            for key in keys:
                annotations[key] = ("ERROR", f"{type(error).__name__}: {error}")
            return None

    def finish() -> dict:
        values, statuses, reasons, originals = [], [], [], []
        for key in P_KEYS:
            p = raw[key]
            available = _finite(p) and 0 <= p <= 1
            status, reason = annotations.get(key, ("PASS" if available and p >= .01 else
                                                   "FAIL" if available else "ERROR", None))
            if not available and reason is None:
                status, reason = "ERROR", "NONFINITE_OR_OUT_OF_RANGE_P"
            if status == "ERROR" and not any(key in f.get("tests", []) for f in out["failures"]):
                out["failures"].append({"phase": "battery", "end": int(end), "tests": [key], "reason": reason})
            values.append(float(p) if available else None)
            statuses.append(status)
            reasons.append(reason)
            originals.append(str(p) if p is not None and not _finite(p) else p)
        out.update(p_values=values, p_statuses=statuses, p_reasons=reasons,
                   p_original_values=originals, complete=not out["failures"],
                   complete_prime="ERROR" not in statuses[:8])
        out["S"], out["available_p_count"] = statistic(values)
        out["S_prime"], out["available_p_count_prime"] = statistic(values[:8])
        return out

    marks = np.asarray(marks)
    if marks.ndim != 1 or end < 512 or len(marks) < end:
        annotations.update({k: ("ERROR", "INCOMPLETE_MARKS") for k in P_KEYS})
        return finish()
    lo, hi = h1lib.window_bytes(marks, end - 512, end)
    boundaries = np.r_[lo, marks[end - 512:end]]
    if lo < 0 or hi > len(emitted) or np.any(np.diff(boundaries) <= 0) or np.any(boundaries % 8):
        annotations.update({k: ("ERROR", "INCOMPLETE_OR_INVALID_BYTE_RANGE") for k in P_KEYS})
        return finish()
    data = emitted[lo:hi]
    bits = np.unpackbits(np.frombuffer(data[:125000], dtype=np.uint8), bitorder="big")
    n, N = len(bits), len(data) // 8
    P = N - 1
    out.update(byte_start=lo, byte_end=hi, original_bytes=hi - lo, original_bits=8 * (hi - lo),
               bit_used_bytes=n // 8, bit_used_bits=n, word_used_bytes=hi - lo,
               word_count=N, pair_count=P, was_cut=8 * (hi - lo) > 1000000)
    for key in STAT_FIELDS:
        _log(out, key, None, "COMPUTATION_UNAVAILABLE", "ERROR")
    M = max(128, n // 99 + 1)
    base = calculate(P_KEYS[:3], lambda: evaluate_bits(bits, alpha=.01, block_size=M)["tests"])
    if base is not None:
        for key in P_KEYS[:3]:
            raw[key] = base[key]["p_value"]
            annotations[key] = (base[key]["status"], base[key].get("reason"))
        _log(out, "block_frequency_chi2", base["block_frequency"].get("statistic"),
             base["block_frequency"].get("reason"))
    for key, value in (("block_size", M), ("blocks", n // M), ("unused_bits", n % M)):
        _log(out, key, value)
    O = int(bits.sum(dtype=np.int64))
    pi = O / n
    V = 1 + int(np.count_nonzero(bits[1:] != bits[:-1]))
    residual = V - 2 * n * pi * (1 - pi)
    for key, value in (("ones", O), ("frequency_signed_z", (2 * O - n) / math.sqrt(n)),
                       ("V", V), ("pi", pi), ("runs_residual", residual),
                       ("runs_signed_z", residual / (2 * math.sqrt(n) * pi * (1 - pi)) if 0 < pi < 1 else None),
                       ("prerequisite_failed", abs(pi - .5) >= 2 / math.sqrt(n))):
        _log(out, key, value)
    counts = calculate(P_KEYS[3:6], lambda: h1lib.pattern_counts(bits, 9))
    if counts is not None:
        serial = calculate(P_KEYS[3:5], lambda: h1lib.serial_test(bits, 8, counts))
        if serial is not None:
            raw["serial_1"], raw["serial_2"] = serial
        raw["approximate_entropy"] = calculate(["approximate_entropy"],
                                               lambda: h1lib.approximate_entropy_test(bits, 8, counts))
        psi = [float(2**k / n * np.sum(counts[k].astype(float)**2) - n) for k in (6, 7, 8)]
        for k in (6, 7, 8, 9):
            _log(out, f"pattern_counts_{k}", counts[k])
        phi = []
        for k in (8, 9):
            c = counts[k] / n
            c = c[c > 0]
            phi.append(float(np.sum(c * np.log(c))))
        apen = phi[0] - phi[1]
        for key, value in (("serial_psi2", psi), ("serial_delta1", psi[2] - psi[1]),
                           ("serial_delta2", psi[2] - 2 * psi[1] + psi[0]),
                           ("approximate_entropy_phi8", phi[0]), ("approximate_entropy_phi9", phi[1]),
                           ("approximate_entropy", apen), ("approximate_entropy_chi2", 2 * n * (math.log(2) - apen))):
            _log(out, key, value)
    raw["cumulative_sums"] = calculate(["cumulative_sums"], lambda: h1lib.cumulative_sums_forward(bits))
    _log(out, "cumulative_sums_max_abs", int(np.max(np.abs(np.cumsum(2 * bits.astype(np.int64) - 1)))))
    longest = calculate(["longest_run"], lambda: h1lib.longest_run_test(bits))
    M = h1lib.longest_run_block_size(n)
    K, edges, probs = h1lib.LONGEST_RUN_TABLES[M]
    J = n // M
    if longest is not None:
        raw["longest_run"], nu = longest
        chi2 = sum((nu[j] - J * probs[j])**2 / (J * probs[j]) for j in range(K + 1))
        for key, value in (("longest_run_chi2", chi2), ("longest_run_categories", nu),
                           ("longest_run_block_size", M), ("longest_run_blocks", J),
                           ("longest_run_unused_bits", n - J * M), ("longest_run_table_key", M),
                           ("longest_run_edges", edges)):
            _log(out, key, value)
    with np.errstate(invalid="ignore", divide="ignore"):
        words = calculate(P_KEYS[8:], lambda: h1lib.word_tests(data))
    if words is not None:
        raw.update({k: words[k] for k in P_KEYS[8:]})
        r, duplicates = words["lag1_r"], words["duplicate_count"]
        u = (np.frombuffer(data, dtype=">u8") >> np.uint64(11)).astype(float)
        if len(u) < 3 or np.ptp(u[:-1]) == 0 or np.ptp(u[1:]) == 0:
            # signal_schema lag1: zero-variance pairs make r and its p unavailable.
            r, raw["lag1_correlation"] = float("nan"), None
        if not _finite(r):
            annotations["lag1_correlation"] = ("NOT_APPLICABLE", "UNDEFINED_CORRELATION")
        for key, value in (("lag1_r", r), ("lag1_test_z", math.sqrt(N) * r),
                           ("lag1_feature", math.sqrt(P) * r), ("duplicate_count", duplicates),
                           ("duplicate_rate", duplicates / P), ("duplicate_lambda", P * 2.0**-64)):
            _log(out, key, value)
    return finish()


def application_checkpoint(CV: np.ndarray, TD: np.ndarray, marks: np.ndarray, end: int) -> dict:
    """signal_schema.application_statistics: signed TD-CV and complete local ledger."""
    out = {"end": int(end), "count_D": 1024, "finite_input_count": 0,
           "statistic_statuses": {}, "statistic_reasons": {}}
    start = end - 1024
    cv, td, marks = np.asarray(CV), np.asarray(TD), np.asarray(marks)
    present = start >= 0 and cv.ndim == td.ndim == 1
    count = max(0, min(end, len(cv), len(td)) - start) if present else 0
    with np.errstate(over="ignore", invalid="ignore"):
        diff = td[start:start + count] - cv[start:start + count] if present else np.array([])
        out["finite_input_count"] = int(np.isfinite(diff).sum())
        mean = float(diff.mean()) if count == 1024 and np.isfinite(diff).all() else None
    for key in ("count_D", "finite_input_count"):
        _log(out, key, out[key])
    _log(out, "mean_D", mean, None if _finite(mean) else "INCOMPLETE_OR_NONFINITE_INPUT", "ERROR")
    ledger_ok = start >= 0 and marks.ndim == 1 and len(marks) >= end
    if ledger_ok:
        lo, _ = h1lib.window_bytes(marks, start, end)
        increments = np.diff(np.r_[lo, marks[start:end]])
        ledger_ok = bool(np.isfinite(increments).all() and np.all((increments >= 128) &
                              (increments <= 256) & ((increments - 128) % 8 == 0)))
    uphill = int(((increments - 128) // 8).sum()) if ledger_ok else None
    for key, value in (("uphill_count", uphill), ("proposal_count", 16 * 1024 if ledger_ok else None),
                       ("uphill_fraction", uphill / (16 * 1024) if ledger_ok else None)):
        _log(out, key, value, None if ledger_ok else "INVALID_OR_MISSING_MARK_LEDGER", "ERROR")
    return out


def standardise_application(checkpoint: dict, constants: dict) -> dict:
    """signal_schema.application_statistics.constants_call: signed calibrated signals."""
    out = {**checkpoint, "statistic_statuses": dict(checkpoint.get("statistic_statuses", {})),
           "statistic_reasons": dict(checkpoint.get("statistic_reasons", {}))}
    for key, raw, sd, center in (("A1_over_sigma_D", "mean_D", "sigma_D", 0.0),
                                ("A2_standardised", "uphill_fraction", "sigma_f", constants.get("mu_f"))):
        value, scale = checkpoint.get(raw), constants.get(sd)
        good = _finite(value) and _finite(scale) and scale > 0 and _finite(center)
        _log(out, key, (value - center) / scale if good else None,
             None if good else "UNAVAILABLE_SIGNAL_OR_CALIBRATION")
    return out


def endpoint(cv: np.ndarray, n: int, E_M: float) -> dict:
    """endpoint.call/_convention: selected 16-batch result, no rescue on exception."""
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n <= 0 or n % 16:
        raise ValueError("endpoint length must be a positive multiple of 16")
    n = int(n)
    values = np.asarray(cv, dtype=np.float64)
    prefix = values[:n] if values.ndim else values.reshape(1)
    out = {"L": n, "input_count": len(prefix), "source_status": None, "valid": False, "y": 0,
           "covered": None, "invalid_reason": None, "exception": None, "estimate": None,
           "signed_error": None, "mcse": None, "half_width": None, "interval_95": None,
           "batch_means": None, "severity_path": {"array": "severity", "start": 0, "stop": n}}
    if not np.isfinite(prefix).all():
        out["invalid_reason"] = "NONFINITE_INPUT"
        return out
    if values.ndim != 1 or len(prefix) != n:
        out["invalid_reason"] = "INCOMPLETE_TRACE"
        return out
    with np.errstate(over="ignore", invalid="ignore"):
        means = prefix.reshape(16, n // 16).mean(axis=1)
        estimate = float(means.mean())
        out.update(batch_means=_json_safe(means) if np.isfinite(means).all() else None,
                   estimate=_json_safe(estimate), signed_error=_json_safe(estimate - E_M))
    try:
        result = summarize_series(prefix, reference=E_M, batch_sizes=(n // 16,))
    except Exception as error:
        out.update(invalid_reason="EXCEPTION", exception=_exception(error))
        return out
    try:
        batch = result["batch_means"][str(n // 16)]
        status = batch["status"]
        out["source_status"] = status if isinstance(status, str) else None
        estimate, mcse, bounds, covered = (batch[k] for k in ("prefix_mean", "mcse", "interval_95", "covered"))
        if not isinstance(status, str) or (covered is not None and type(covered) is not bool):
            raise ValueError("invalid status or coverage type")
        if any(v is not None and not _finite(v) for v in (estimate, mcse)):
            raise ValueError("nonfinite endpoint value")
        if bounds is not None and (len(bounds) != 2 or not all(_finite(v) for v in bounds)):
            raise ValueError("malformed interval")
        if status == "ESTIMATED" and (estimate is None or mcse is None or bounds is None or type(covered) is not bool):
            raise ValueError("incomplete estimated endpoint")
        if estimate is not None:
            out.update(estimate=float(estimate), signed_error=_json_safe(estimate - E_M))
        out.update(mcse=mcse, interval_95=bounds, covered=covered,
                   half_width=float(t.ppf(.975, df=15)) * mcse if mcse is not None else None)
        if any(out[k] is not None and not _finite(out[k]) for k in ("signed_error", "half_width")):
            raise ValueError("nonfinite endpoint diagnostic")
    except (KeyError, TypeError, ValueError, IndexError) as error:
        out.update(invalid_reason="MALFORMED_RESULT", exception=_exception(error))
        return _json_safe(out)
    if status != "ESTIMATED":
        out["invalid_reason"] = status if status in ("DEGENERATE_BATCH_VARIANCE", "INSUFFICIENT_BATCHES") else "OTHER_STATUS"
    elif not covered:
        out["invalid_reason"] = "UNCOVERED"
    else:
        out.update(valid=True, y=1)
    return out


def issuance(checkpoints: list[dict], threshold: float, *, prime: bool = False) -> dict:
    """operational_alarm.issuance/C_prime: chronological strict exceedance, no onset input.

    A missing/failed earlier checkpoint makes first issuance unknown. A known
    first alarm survives later missing checkpoints; a whole-horizon maximum does not.
    """
    if not _finite(threshold):
        raise ValueError("finite calibrated threshold required")
    key = "S_prime" if prime else "S"
    by_end = {c["end"]: c for c in checkpoints}
    if len(by_end) != len(checkpoints) or any(e not in BATTERY_ENDS for e in by_end):
        raise ValueError("duplicate or off-grid checkpoint")
    def good(c: dict | None) -> bool:
        return (c is not None and c.get("complete_prime" if prime else "complete", c.get("complete", True))
                and _finite(c.get(key)))
    complete = all(good(by_end.get(e)) for e in BATTERY_ENDS)
    observed_alarm = any(good(c) and c[key] > threshold for c in checkpoints)
    out = {"whole_horizon_max": max(c[key] for c in checkpoints) if complete else None,
           "complete_horizon": complete, "alarm_sweep": None, "S_at_alarm": None,
           "whole_horizon_alarm": True if observed_alarm else None, "status": "UNKNOWN"}
    for end in BATTERY_ENDS:
        c = by_end.get(end)
        if not good(c):
            return out
        if c[key] > threshold:
            out.update(alarm_sweep=end, S_at_alarm=c[key], whole_horizon_alarm=True, status="ALARM")
            return out
    out.update(whole_horizon_alarm=False, status="NO_ALARM")
    return out


def onset_classification(alarm_sweep: int | None, t0: int | None, *, known: bool = True) -> dict:
    """operational_alarm.onset_classification: classify only after onset-blind issuance."""
    label = ("unknown" if not known else "no_alarm" if alarm_sweep is None else
             None if t0 is None else "pre_onset" if alarm_sweep <= t0 else "post_onset")
    return {"onset_classification": label,
            "alarm_minus_t0": int(alarm_sweep - t0) if known and alarm_sweep is not None and t0 is not None else None}


def eligibility(alarm_sweep: int | None, L: int, *, known: bool = True) -> dict:
    """operational_alarm.eligibility: strict alarm<L; unknown is not no alarm."""
    eligible = alarm_sweep is not None and alarm_sweep < L
    return {"L": int(L), "eligible": bool(eligible) if known else None,
            "reason": "UNKNOWN_ISSUANCE" if not known else None if eligible else
            "no_alarm" if alarm_sweep is None else "alarm_at_or_after_L"}


def issuance_endpoint(cv: np.ndarray, alarm_sweep: int | None, E_M: float = E_M) -> dict:
    """endpoint.issuance: no alarm has null n, endpoint and already_invalid."""
    n = 16 * (alarm_sweep // 16) if alarm_sweep is not None else None
    result = endpoint(cv, n, E_M) if n is not None else None
    return {"n": n, "endpoint": result, "already_invalid": not result["valid"] if result is not None else None}


def short_interval_check(cv: np.ndarray, E_M: float = E_M) -> list[dict]:
    """endpoint.issuance.short_interval_check: per-run fixed-checkpoint prefixes."""
    return [endpoint(cv, end, E_M) for end in BATTERY_ENDS]


def _source_manifests(run: dict) -> dict:
    """storage.run_record: reconstructed source descriptions, never invented timestamps.

    dlib returns seeds, not live manifests. Reconstruct constructor metadata with
    fresh, unread sources; created_at is unavailable for the original instances.
    """
    sources = {role: dlib.PCG64Source(run["seeds"][role]) for role in ("local", "centroid", "start")}
    local = dlib.FaultyWords(sources["local"], run["family"], run["profile"], run["seeds"]["fault"])
    sources.update(local_adapter=local, centroid_clock=dlib.SweepClock(sources["centroid"], local))
    manifests = {role: source.manifest.to_dict() for role, source in sources.items()}
    for manifest in manifests.values():
        manifest.update(created_at=None, provenance="reconstructed_from_seeds_and_constructor")
    manifests.update({role: {"generator": "numpy.random.default_rng", "seed": run["seeds"][role]}
                      for role in ("fault", "profile")})
    return manifests


def measure_run(run: dict) -> tuple[dict, dict[str, np.ndarray]]:
    """storage.run_record/failed_runs: all uncalibrated per-run measurements, in memory."""
    started = time.perf_counter()
    arrays = {k: np.asarray(run[k]) for k in ARRAY_DTYPES}
    record = {k: v for k, v in run.items() if k not in (*ARRAY_DTYPES, "emitted", "raw_base")}
    record.update(split="fitting" if run["family"] != "healthy" and run["replica"] < 600 else
                  "evaluation" if run["family"] != "healthy" else
                  "calibration" if run["condition"] == "calibration" else "validation",
                  source_manifests=_source_manifests(run), failures=list(run.get("failures") or []))
    record["conventional"] = [battery_checkpoint(run["emitted"], arrays["marks"], end) for end in BATTERY_ENDS]
    record["application"] = [application_checkpoint(arrays["CV"], arrays["TD"], arrays["marks"], end)
                             for end in APPLICATION_ENDS]
    record["endpoints"] = {str(L): endpoint(arrays["CV"], L, E_M) for L in LENGTHS}
    flags = [record["endpoints"][str(L)]["valid"] for L in LENGTHS]
    record["validity_flags"] = flags
    record["validity_pattern"] = "".join(str(int(v)) for v in flags)
    record["validity_transitions"] = {f"{a}:{b}": [flags[i], flags[j]]
                                      for i, a in enumerate(LENGTHS) for j, b in enumerate(LENGTHS) if i < j}
    if run["condition"] == "healthy":
        record["short_intervals"] = short_interval_check(arrays["CV"])
    for c in record["conventional"]:
        record["failures"].extend(c["failures"])
    for c in record["application"]:
        if "ERROR" in c["statistic_statuses"].values():
            record["failures"].append({"phase": "application", "end": c["end"], "reasons": c["statistic_reasons"]})
    for e in [*record["endpoints"].values(), *record.get("short_intervals", [])]:
        if e["invalid_reason"] not in (None, "UNCOVERED", "DEGENERATE_BATCH_VARIANCE", "INSUFFICIENT_BATCHES"):
            record["failures"].append({"phase": "endpoint", "L": e["L"], "reason": e["invalid_reason"],
                                       "exception": e["exception"]})
    for key, a in arrays.items():
        if a.shape != (H,) or a.dtype != np.dtype(ARRAY_DTYPES[key]) or not np.isfinite(a).all():
            record["failures"].append({"phase": "trace", "array": key, "reason": "SHAPE_DTYPE_OR_NONFINITE"})
    check = dlib.check_i1(arrays["marks"], run["local_bytes"], H, False,
                          current_sweep=run.get("current_sweep"))
    record["I-1"] = check
    if not check["passed"]:
        record["failures"].append({"phase": "I-1", "reason": check["mismatches"]})
    record.update(status="FAILED" if record["failures"] else "SUCCESS",
                  failure=record["failures"][0] if record["failures"] else None,
                  completed_sweeps=len(arrays["marks"]), last_completed_sweep=len(arrays["marks"]) - 1,
                  completed_battery_checkpoints=sum(c["complete"] for c in record["conventional"]),
                  completed_application_checkpoints=sum("ERROR" not in c["statistic_statuses"].values()
                                                        for c in record["application"]))
    record = _json_safe(record)
    record["measurement_seconds"] = time.perf_counter() - started
    return record, arrays


def array_hashes(arrays: dict[str, np.ndarray]) -> dict:
    """storage.hashes: canonical little-endian C-order bytes, marks stay integers."""
    return {key: {"sha256": hashlib.sha256(np.ascontiguousarray(arrays[key], dtype=dtype).tobytes(order="C")).hexdigest(),
                  "dtype": dtype, "shape": list(arrays[key].shape)} for key, dtype in ARRAY_DTYPES.items()}


def _paths(condition: str, i: int) -> tuple[Path, Path]:
    if condition not in ("calibration", "healthy", "fault/SHIFT", "fault/LAG1", "fault/DUP") or i < 0:
        raise ValueError("unknown condition or negative replica")
    return Path("traces") / condition / f"{i:05d}.npz", Path("results/runs") / condition / f"{i:05d}.json"


def write_run(record: dict, arrays: dict[str, np.ndarray], root: Path) -> dict:
    """storage: write canonical files under explicit root, return manifest entry.

    JSON's own file hash lives in the returned entry (it cannot hash itself).
    Serialization errors propagate with phase context for the campaign's failure
    recorder; this function never retries or substitutes a run.
    """
    root = Path(root)
    trace, result = _paths(record["condition"], record["replica"])
    stored = {}
    for key, dtype in ARRAY_DTYPES.items():
        a = np.asarray(arrays[key])
        if a.ndim != 1 or (a.shape != (H,) and record.get("status") != "FAILED"):
            raise ValueError(f"storage: invalid {key} shape")
        if a.dtype != np.dtype(dtype):
            raise ValueError(f"storage: invalid {key} dtype")
        if not np.isfinite(a).all() and record.get("status") != "FAILED":
            raise ValueError(f"storage: nonfinite {key} in successful run")
        stored[key] = np.ascontiguousarray(a, dtype=dtype)
    hashes = array_hashes(stored)
    saved = _json_safe({**record, "trace_path": trace.as_posix(), "array_hashes": hashes,
                        "compression": {"method": "numpy.savez_compressed ZIP_DEFLATED defaults",
                                        "numpy": np.__version__, "python": platform.python_version()}})
    # Validate JSON before touching disk; trace hash is a string of fixed length.
    json.dumps(saved, allow_nan=False)
    try:
        (root / trace).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(root / trace, **stored)
        saved["trace_sha256"] = hashlib.sha256((root / trace).read_bytes()).hexdigest()
        encoded = (json.dumps(saved, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")
        (root / result).parent.mkdir(parents=True, exist_ok=True)
        (root / result).write_bytes(encoded)
    except Exception as error:
        raise RuntimeError(f"serialization: {type(error).__name__}: {error}") from error
    return {"namespace": record["namespace"], "condition": record["condition"], "replica": record["replica"],
            "files": {trace.as_posix(): saved["trace_sha256"], result.as_posix(): hashlib.sha256(encoded).hexdigest()},
            "array_hashes": hashes, "status": record.get("status")}


def load_run(root: Path, condition: str, i: int, *, manifest: dict | None = None) -> tuple[dict, dict]:
    """storage.trace_file/failed_runs.detection: pickle-free load and hash validation."""
    root = Path(root)
    trace, result = _paths(condition, i)
    if manifest is not None:
        for path in (trace, result):
            if hashlib.sha256((root / path).read_bytes()).hexdigest() != manifest["files"][path.as_posix()]:
                raise ValueError(f"file hash mismatch: {path}")
    record = json.loads((root / result).read_text(encoding="utf-8"))
    if record["condition"] != condition or record["replica"] != i or record["trace_path"] != trace.as_posix():
        raise ValueError("run identity/path mismatch")
    if hashlib.sha256((root / trace).read_bytes()).hexdigest() != record["trace_sha256"]:
        raise ValueError("trace file hash mismatch")
    with np.load(root / trace, allow_pickle=False) as archive:
        if set(archive.files) != set(ARRAY_DTYPES):
            raise ValueError("unexpected trace arrays")
        arrays = {key: archive[key] for key in ARRAY_DTYPES}
    for key, a in arrays.items():
        if a.dtype.str != ARRAY_DTYPES[key] or a.ndim != 1 or (a.shape != (H,) and record.get("status") != "FAILED"):
            raise ValueError(f"invalid trace dtype/shape: {key}")
        if not np.isfinite(a).all() and record.get("status") != "FAILED":
            raise ValueError(f"nonfinite successful trace: {key}")
    if array_hashes(arrays) != record["array_hashes"]:
        raise ValueError("array hash mismatch")
    return record, arrays
