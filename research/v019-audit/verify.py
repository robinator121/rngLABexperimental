"""Independent verification of stored v0.19 Stage D outputs; writes machine-readable discrepancies.

    PYTHONPATH=src python research/v019-audit/verify.py --root <output root> --namespace <namespace> \
        [--deep all|sample|none] [--checkpoint-stride 8] [--workers 4] [--out verification.json]

Reads only what is stored (results/runs/{condition}/{i:05d}.json and
traces/{condition}/{i:05d}.npz, the layout of protocol.json storage) and
recomputes with research/v019-audit/recount.py and replay.py. It never calls
the pilot's own measurement or analysis functions and never writes into the
verified tree unless --out points there.

Two levels:
- per run: seeds, byte ledger, profile redraw, severity, onset offset, an
  independent regeneration of the whole run (marks, traces, stream hashes),
  battery statistics at a subset of checkpoints, application signals, endpoints;
- campaign: alarm thresholds, gate V2-D, issuance, eligibility, stop criteria
  1 to 3, the forecaster features, a reference fit, scores and criterion 4.

Every disagreement is a record {scope, id, field, stored, independent}. The
field names of the stored run record are read in `view_from_record` and
`verify_run` only; a change of the record schema is absorbed there.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import platform
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

_HERE = Path(__file__).resolve().parent


def _load(name: str, path: Path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


rc = _load("v019_audit_recount", _HERE / "recount.py")
rp = _load("v019_audit_replay", _HERE / "replay.py")

PLAN = {"calibration": ("healthy", 999), "healthy": ("healthy", 2000), "fault/SHIFT": ("SHIFT", 1000),
        "fault/LAG1": ("LAG1", 1000), "fault/DUP": ("DUP", 1000)}
FITTING, EVALUATION = range(0, 600), range(600, 1000)
ARRAYS = {"CV": "<f8", "TD": "<f8", "marks": "<i8", "severity": "<f8"}
BATTERY_FIELDS = ("frequency_signed_z", "block_frequency_chi2", "runs_residual", "runs_signed_z", "serial_delta1",
                  "serial_delta2", "approximate_entropy_chi2", "cumulative_sums_max_abs", "longest_run_chi2",
                  "longest_run_categories", "lag1_r", "lag1_feature", "duplicate_count", "duplicate_rate",
                  "original_bits", "bit_used_bits", "word_count", "pair_count", "was_cut", "S", "S_prime")


@dataclass
class RunView:
    """The stored quantities the campaign recount needs from one run."""
    condition: str
    family: str
    replica: int
    t0: int | None
    t_cap: float | None
    s_max: float | None
    statistic: np.ndarray        # S at the 159 battery checkpoints
    statistic_prime: np.ndarray  # S over the eight SP 800-22 columns
    target: np.ndarray | None    # the family's target statistic per checkpoint, NaN when unavailable
    mean_d: np.ndarray           # application checkpoints: signed mean of TD - CV
    uphill_fraction: np.ndarray
    cv: np.ndarray
    complete: bool = True
    notes: list[str] = field(default_factory=list)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(array: np.ndarray, dtype: str) -> str:
    return hashlib.sha256(np.ascontiguousarray(array, dtype=dtype).tobytes(order="C")).hexdigest()


def stored_paths(root: Path, condition: str, replica: int) -> tuple[Path, Path]:
    return (Path(root) / "results/runs" / condition / f"{replica:05d}.json",
            Path(root) / "traces" / condition / f"{replica:05d}.npz")


def load_stored(root: Path, condition: str, replica: int) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    record_path, trace_path = stored_paths(root, condition, replica)
    record = json.loads(record_path.read_text(encoding="utf-8"))
    with np.load(trace_path, allow_pickle=False) as archive:
        arrays = {key: archive[key] for key in archive.files}
    return record, arrays


def _nan(value: Any) -> float:
    return float("nan") if value is None else float(value)


def view_from_record(record: Mapping[str, Any], arrays: Mapping[str, np.ndarray]) -> RunView:
    family = record["family"]
    profile = record.get("profile") or {}
    battery = record["conventional"]
    target = None if family == "healthy" else np.array([_nan(c.get(rc.TARGETS[family])) for c in battery])
    return RunView(condition=record["condition"], family=family, replica=int(record["replica"]),
                   t0=profile.get("t0"), t_cap=profile.get("T_cap"), s_max=profile.get("s_max"),
                   statistic=np.array([_nan(c.get("S")) for c in battery]),
                   statistic_prime=np.array([_nan(c.get("S_prime")) for c in battery]), target=target,
                   mean_d=np.array([_nan(c.get("mean_D")) for c in record["application"]]),
                   uphill_fraction=np.array([_nan(c.get("uphill_fraction")) for c in record["application"]]),
                   cv=np.asarray(arrays["CV"], dtype=np.float64),
                   complete=record.get("status", "SUCCESS") == "SUCCESS")


def view_from_failure_record(record: Mapping[str, Any]) -> RunView:
    """A run that has only a failure record: no statistics, no trace. Its issuance is unknown and
    every endpoint is invalid (protocol.json failed_runs.endpoints)."""
    empty = np.zeros(0, dtype=np.float64)
    profile = record.get("profile") or {}
    return RunView(condition=record["condition"], family=record["family"], replica=int(record["replica"]),
                   t0=profile.get("t0"), t_cap=profile.get("T_cap"), s_max=profile.get("s_max"), statistic=empty,
                   statistic_prime=empty, target=None, mean_d=empty, uphill_fraction=empty, cv=empty, complete=False,
                   notes=["failure record"])


def load_view(root: Path, condition: str, replica: int) -> RunView:
    record = json.loads(stored_paths(root, condition, replica)[0].read_text(encoding="utf-8"))
    if record.get("trace_path") is None:
        return view_from_failure_record(record)
    return view_from_record(*load_stored(root, condition, replica))


# ---------------------------------------------------------------------------
# Per-run verification
# ---------------------------------------------------------------------------
def verify_run(record: Mapping[str, Any], arrays: Mapping[str, np.ndarray], namespace: str, *,
               deep: bool = True, checkpoint_stride: int = 8) -> dict[str, Any]:
    """Compare one stored run with its independent recomputation."""
    condition, family, replica = record["condition"], record["family"], int(record["replica"])
    run_id = f"{condition}/{replica:05d}"
    issues: list[dict[str, Any]] = []
    compared = 0

    def check(name: str, stored: Any, independent: Any, *, rtol: float = 1e-12, atol: float = 1e-12) -> None:
        nonlocal compared
        compared += 1
        if not rc.compare(stored, independent, rtol=rtol, atol=atol):
            issues.append({"scope": "run", "id": run_id, "field": name,
                           "stored": _brief(stored), "independent": _brief(independent)})

    check("namespace", record.get("namespace"), namespace)
    seeds = rc.run_seeds(namespace, condition, replica)
    check("seeds", {k: record["seeds"].get(k) for k in seeds}, seeds)
    for key, dtype in ARRAYS.items():
        array = np.asarray(arrays.get(key, np.zeros(0)))
        check(f"{key}.shape", list(array.shape), [rc.H])
        check(f"{key}.dtype", array.dtype.str, dtype)
        stored_hash = (record.get("array_hashes") or {}).get(key, {}).get("sha256")
        if stored_hash is not None:
            check(f"{key}.sha256", stored_hash, array_sha256(array, dtype))
    aligned = all(np.asarray(arrays.get(key, np.zeros(0))).shape == (rc.H,) for key in ARRAYS)
    marks = np.asarray(arrays["marks"])
    ledger = rc.check_ledger(marks, rc.H)
    check("ledger", [], ledger)
    check("local_bytes", record.get("local_bytes"), int(marks[-1]) if len(marks) else None)
    profile = record.get("profile") or {}
    if family == "healthy":
        expected = {"t0": None, "T_cap": None, "r": None, "s_max": None}
    else:
        drawn = rc.draw_profile(seeds["profile"], rc.protocol_caps()[family])
        expected = {k: (int(v) if k == "t0" else float(v)) for k, v in drawn.items()}
    check("profile", {k: profile.get(k) for k in expected}, expected, rtol=0.0, atol=0.0)
    severity = rc.severity_path(expected["t0"], expected["T_cap"], expected["s_max"])
    check("severity.bitwise", array_sha256(arrays["severity"], "<f8"), array_sha256(severity, "<f8"))
    onset = None if family == "healthy" or ledger else int(marks[expected["t0"] - 1])
    check("B0", record.get("B0"), onset)
    regenerated = None
    if deep and not ledger and aligned:
        try:
            regenerated = rp.replay(namespace, condition, family, replica)
        except ValueError as error:  # a LAG1 domain error is a run failure, not a discrepancy of the record
            issues.append({"scope": "run", "id": run_id, "field": "replay", "stored": record.get("status"),
                           "independent": f"{type(error).__name__}: {error}"})
        if regenerated is not None:
            for key in ("marks", "CV", "TD"):
                check(f"replay.{key}.bitwise", array_sha256(arrays[key], ARRAYS[key]), array_sha256(regenerated[key], ARRAYS[key]))
            for key in ("emitted_sha256", "base_sha256", "prefix_sha256", "start_sha256", "local_bytes", "B0"):
                check(f"replay.{key}", record.get(key), regenerated[key])
    battery = record.get("conventional") or []
    check("battery.ends", [c.get("end") for c in battery], rc.BATTERY_ENDS.tolist())
    for index, checkpoint in enumerate(battery):
        s, count = rc.alarm_statistic(checkpoint.get("p_values") or [])
        check(f"battery[{checkpoint.get('end')}].S_from_stored_p", checkpoint.get("S"), s)
        s_prime, _ = rc.alarm_statistic((checkpoint.get("p_values") or [])[:8])
        check(f"battery[{checkpoint.get('end')}].S_prime_from_stored_p", checkpoint.get("S_prime"), s_prime)
        check(f"battery[{checkpoint.get('end')}].available_p_count", checkpoint.get("available_p_count"), count)
        if regenerated is None or index % checkpoint_stride:
            continue
        end = int(checkpoint["end"])
        lo, hi = rc.window_words(marks, end - rc.BATTERY_WINDOW, end)
        mine = rc.battery_checkpoint(regenerated["emitted_words"][lo:hi])
        for key, stored, independent in zip(rc.P_KEYS, checkpoint["p_values"], mine["p_values"]):
            if stored is not None and independent is not None and max(stored, independent) <= 1e-290:
                compared += 1  # both at or below the floor: the alarm statistic is identical
                continue
            check(f"battery[{end}].p.{key}", stored, independent, rtol=1e-9, atol=1e-300)
        for key in BATTERY_FIELDS:
            if key in checkpoint and key in mine:
                check(f"battery[{end}].{key}", checkpoint[key], mine[key], rtol=1e-9, atol=1e-9)
    application = record.get("application") or []
    check("application.ends", [c.get("end") for c in application], rc.APPLICATION_ENDS.tolist())
    if not ledger and aligned:  # a malformed trace is already reported above; its windows are not comparable
        difference = np.asarray(arrays["TD"], dtype=np.float64) - np.asarray(arrays["CV"], dtype=np.float64)
        uphill = (np.diff(marks, prepend=0) - 128) // 8
        for checkpoint in application:
            end = int(checkpoint["end"])
            window = slice(end - rc.APPLICATION_WINDOW, end)
            check(f"application[{end}].mean_D", checkpoint.get("mean_D"), math.fsum(difference[window].tolist()) / 1024)
            count = int(uphill[window].sum())
            check(f"application[{end}].uphill_count", checkpoint.get("uphill_count"), count)
            check(f"application[{end}].uphill_fraction", checkpoint.get("uphill_fraction"), count / (16 * 1024))
    reference = rc.reference_energy()
    near_boundary = []
    for length in rc.LENGTHS:
        stored = (record.get("endpoints") or {}).get(str(length)) or {}
        mine = rc.endpoint(arrays["CV"], length, reference)
        check(f"endpoint[{length}].valid", stored.get("valid"), mine["valid"])
        check(f"endpoint[{length}].reason", stored.get("invalid_reason"), mine["reason"])
        for key in ("estimate", "signed_error", "mcse", "half_width", "interval_95", "covered", "batch_means"):
            check(f"endpoint[{length}].{key}", stored.get(key), mine[key])
        if mine["margin"] is not None and abs(mine["margin"]) <= 1e-12:
            near_boundary.append(length)
    for stored in record.get("short_intervals") or []:
        mine = rc.endpoint(arrays["CV"], int(stored["L"]), reference)
        check(f"short_interval[{stored['L']}].valid", stored.get("valid"), mine["valid"])
        check(f"short_interval[{stored['L']}].half_width", stored.get("half_width"), mine["half_width"])
    return {"id": run_id, "compared": compared, "issues": issues, "deep": regenerated is not None,
            "near_boundary_endpoints": near_boundary, "stored_status": record.get("status")}


def _brief(value: Any) -> Any:
    """A JSON-safe, bounded rendering of a compared value."""
    if isinstance(value, np.ndarray):
        value = value.tolist()
    if isinstance(value, (list, tuple)) and len(value) > 20:
        return {"length": len(value), "head": [_brief(v) for v in value[:5]]}
    if isinstance(value, (list, tuple)):
        return [_brief(v) for v in value]
    if isinstance(value, Mapping):
        return {str(k): _brief(v) for k, v in value.items()}
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    return value


# ---------------------------------------------------------------------------
# Campaign recount
# ---------------------------------------------------------------------------
def whole_horizon_max(view: RunView, prime: bool = False) -> float | None:
    s = view.statistic_prime if prime else view.statistic
    return float(s.max()) if view.complete and len(s) == len(rc.BATTERY_ENDS) and np.isfinite(s).all() else None


def issue(view: RunView, threshold: float, prime: bool = False) -> dict[str, Any]:
    s = view.statistic_prime if prime else view.statistic
    try:
        alarm: Any = rc.first_alarm(s, threshold)
    except rc.RecountError:
        return {"alarm_sweep": "UNKNOWN", "class": "unknown", "alarm_minus_t0": None}
    return {"alarm_sweep": alarm, "class": rc.onset_class(alarm, view.t0),
            "alarm_minus_t0": None if alarm is None or view.t0 is None else alarm - view.t0}


def features(view: RunView, alarm: int, constants: Mapping[str, float]) -> dict[str, Any]:
    """P forecasters.features_in_order for one eligible run."""
    j = int(np.searchsorted(rc.BATTERY_ENDS, alarm))
    values = [None if not math.isfinite(v) else float(v) for v in view.target]
    slope, count = rc.trailing_slope(rc.BATTERY_ENDS.tolist(), values, j)
    k = rc.application_index(alarm)
    a1 = a2 = None
    if k is not None and math.isfinite(view.mean_d[k]) and math.isfinite(view.uphill_fraction[k]):
        a1 = float(view.mean_d[k]) / constants["sigma_D"]
        a2 = (float(view.uphill_fraction[k]) - constants["mu_f"]) / constants["sigma_f"]
    return {"alarm_sweep": float(alarm), "target_at_alarm": values[j], "target_slope": slope, "slope_count": count,
            "application_index": k, "A1_over_sigma_D": a1, "A2_standardised": a2}


_COLUMNS = {"F0+": ("alarm_sweep",), "F-B-lite": ("alarm_sweep", "target_at_alarm", "target_slope"),
            "F-AB-lite": ("alarm_sweep", "target_at_alarm", "target_slope", "A1_over_sigma_D", "A2_standardised")}


def fit_and_predict(fit_rows: list[dict[str, Any]], fit_y: list[int], eval_rows: list[dict[str, Any]],
                    saved: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """P forecasters: F0, F0+, F-B-lite, F-AB-lite with the fallback rules of Q05.

    Without `saved` parameters the logistic models are fitted by the Newton reference
    optimiser; the pilot's own L-BFGS-B fit may differ from it by its stopping tolerance.
    """
    out: dict[str, Any] = {"models": {}, "predictions": {m: [] for m in rc.MODELS}, "fallbacks": {m: {} for m in rc.MODELS},
                           "reasons": {m: [] for m in rc.MODELS}}
    if not fit_y:
        return {**out, "q0": None}
    q0 = sum(fit_y) / len(fit_y)
    out["q0"] = q0
    fitted: dict[str, Any] = {}
    for model, columns in _COLUMNS.items():
        usable = [i for i, row in enumerate(fit_rows) if all(row[c] is not None and math.isfinite(row[c]) for c in columns)]
        y = np.array([fit_y[i] for i in usable], dtype=np.float64)
        info: dict[str, Any] = {"rows": len(usable), "classes": [int((y == 0).sum()), int((y == 1).sum())], "fallback": None}
        if saved is not None:
            if model not in saved:
                raise rc.RecountError(f"saved models hold no entry for {model}")
            if "theta" in saved[model]:
                fitted[model] = saved[model]
                info["source"] = "saved"
            else:
                info.update(source="saved", fallback=saved[model]["fallback"])
        elif len(usable) == 0 or y.min() == y.max():
            info["fallback"] = "EMPTY_SUBSET" if len(usable) == 0 else "ONE_CLASS"
        else:
            x = np.array([[fit_rows[i][c] for c in columns] for i in usable], dtype=np.float64)
            mu, scale = rc.standardise(x)
            result = rc.fit_logistic_newton((x - mu) / scale, y)
            if result["converged"]:
                fitted[model] = {"theta": result["theta"].tolist(), "mu": mu.tolist(), "scale": scale.tolist()}
                info.update(source="newton_reference", objective=result["objective"], gradient_inf=result["gradient_inf"])
            else:
                info["fallback"] = "NOT_CONVERGED"
        out["models"][model] = info

    def predict(model: str, row: Mapping[str, Any]) -> tuple[float, str | None]:
        columns = _COLUMNS[model]
        if model == "F-AB-lite":
            if row["application_index"] is None:
                return predict("F-B-lite", row)[0], "NO_APPLICATION_CHECKPOINT"
            if row["A1_over_sigma_D"] is None or row["A2_standardised"] is None:
                return predict("F-B-lite", row)[0], "MISSING_APPLICATION_FEATURE"
        if any(row[c] is None or not math.isfinite(row[c]) for c in columns if not c.startswith("A")):
            return q0, "MISSING_FEATURE"
        if model not in fitted:
            return q0, out["models"][model]["fallback"]
        parameters = fitted[model]
        value = rc.logistic_predict(parameters["theta"], parameters["mu"], parameters["scale"],
                                    np.array([[row[c] for c in columns]], dtype=np.float64))
        return float(value[0]), None

    for row in eval_rows:
        out["predictions"]["F0"].append(q0)
        out["reasons"]["F0"].append(None)
        for model in _COLUMNS:
            q, reason = predict(model, row)
            out["predictions"][model].append(q)
            out["reasons"][model].append(reason)
            if reason:
                out["fallbacks"][model][reason] = out["fallbacks"][model].get(reason, 0) + 1
    return out


def recount_campaign(views: Mapping[str, list[RunView]], namespace: str, *, plan: Mapping[str, tuple[str, int]] = PLAN,
                     fitting: Iterable[int] = FITTING, evaluation: Iterable[int] = EVALUATION,
                     resamples: int = 10000, saved_models: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Independent campaign analysis from stored per-run quantities.

    `plan` and the split ranges default to the frozen values; tests pass smaller ones, and
    the output records whether the frozen plan was used. The eligibility minimums of the
    stop criteria (300 and 240) are never scaled: a small campaign reports
    INSUFFICIENT_INFORMATION, as the protocol prescribes.
    """
    fitting, evaluation = set(fitting), set(evaluation)
    reference = rc.reference_energy()
    out: dict[str, Any] = {"plan": {k: list(v) for k, v in plan.items()}, "frozen_plan": dict(plan) == PLAN,
                           "missing": {}, "reference_energy": reference}
    for condition, (_, count) in plan.items():
        present = {v.replica for v in views.get(condition, [])}
        out["missing"][condition] = sorted(set(range(count)) - present)
    calibration = sorted(views.get("calibration", []), key=lambda v: v.replica)
    required = plan["calibration"][1]
    order = 950 if required == 999 else math.ceil(0.95 * (required + 1))
    thresholds: dict[str, float | None] = {}
    for label, prime in (("full", False), ("prime", True)):
        maxima = [whole_horizon_max(v, prime) for v in calibration]
        ok = len(maxima) == required and all(m is not None for m in maxima)
        thresholds[label] = rc.order_statistic_threshold(maxima, order, required) if ok else None
    out["thresholds"] = {**thresholds, "order_statistic": order, "calibration_runs": len(calibration)}
    if thresholds["full"] is None:
        return {**out, "status": "INSUFFICIENT_INFORMATION", "reason": "calibration maxima incomplete"}
    pooled_d = np.concatenate([v.mean_d for v in calibration])
    pooled_f = np.concatenate([v.uphill_fraction for v in calibration])
    constants = {"sigma_D": float(np.std(pooled_d, ddof=1)), "mu_f": float(np.mean(pooled_f)),
                 "sigma_f": float(np.std(pooled_f, ddof=1))}
    out["application_constants"] = constants
    endpoints = {(v.condition, v.replica): {L: rc.endpoint(v.cv, L, reference) for L in rc.LENGTHS}
                 for runs in views.values() for v in runs}
    out["sigma_eff"] = {L: (rc.sigma_eff(L, [endpoints[(v.condition, v.replica)][L]["mcse"] for v in calibration])
                            if all(endpoints[(v.condition, v.replica)][L]["mcse"] is not None for v in calibration) else None)
                        for L in rc.LENGTHS}
    healthy = sorted(views.get("healthy", []), key=lambda v: v.replica)
    covered = {L: sum(bool(endpoints[(v.condition, v.replica)][L]["valid"]) for v in healthy) for L in rc.LENGTHS}
    out["criterion_1"] = rc.criterion_1(covered, len(healthy)) if healthy else None
    healthy_issue = [issue(v, thresholds["full"]) for v in healthy]
    unknown = sum(1 for i in healthy_issue if i["alarm_sweep"] == "UNKNOWN")
    alarms = sum(1 for i in healthy_issue if i["alarm_sweep"] not in (None, "UNKNOWN"))
    out["V2_D"] = {**rc.gate_v2d(alarms, len(healthy)), "unknown": unknown} if healthy else None
    if thresholds["prime"] is not None and healthy:
        prime_alarms = sum(1 for v in healthy if issue(v, thresholds["prime"], True)["alarm_sweep"] not in (None, "UNKNOWN"))
        out["V2_D_prime_descriptive"] = rc.gate_v2d(prime_alarms, len(healthy))
    short = {int(end): 0 for end in rc.BATTERY_ENDS[::16]}
    for v in healthy:
        for end in short:
            short[end] += bool(rc.endpoint(v.cv, end, reference)["valid"])
    out["short_interval_coverage"] = {end: {"covered": k, "runs": len(healthy), "cp95": rc.clopper_pearson(k, len(healthy))}
                                      for end, k in short.items()}
    families: dict[str, Any] = {}
    criterion_input: dict[str, Any] = {}
    bootstrap: dict[str, Any] = {}
    for condition, (family, _) in plan.items():
        if family == "healthy":
            continue
        runs = sorted(views.get(condition, []), key=lambda v: v.replica)
        issued = {v.replica: issue(v, thresholds["full"]) for v in runs}
        row: dict[str, Any] = {"runs": len(runs), "classes": {}, "per_length": {}}
        for i in issued.values():
            row["classes"][i["class"]] = row["classes"].get(i["class"], 0) + 1
        alarm_list, valid_list, invalid_list = [], [], []
        for v in runs:
            alarm = issued[v.replica]["alarm_sweep"]
            alarm_list.append(alarm)
            valid_list.append(bool(endpoints[(v.condition, v.replica)][20480]["valid"]))
            n = rc.issuance_length(alarm) if alarm not in (None, "UNKNOWN") else None
            invalid_list.append(n is not None and not rc.endpoint(v.cv, n, reference)["valid"])
        criterion_input[family] = {"alarm": alarm_list, "valid": valid_list, "already_invalid": invalid_list}
        row["per_run"] = [{"replica": v.replica, "alarm_sweep": issued[v.replica]["alarm_sweep"],
                           "class": issued[v.replica]["class"], "alarm_minus_t0": issued[v.replica]["alarm_minus_t0"],
                           "already_invalid": invalid, "valid": {L: bool(endpoints[(v.condition, v.replica)][L]["valid"])
                                                                 for L in rc.LENGTHS}}
                          for v, invalid in zip(runs, invalid_list)]
        x, e = [], []
        for v in runs:
            if v.replica in fitting:
                for L in rc.LENGTHS:
                    error = endpoints[(v.condition, v.replica)][L]["signed_error"]
                    if error is not None:
                        x.append(rc.s_bar_closed_form(L, v.t0, v.t_cap, v.s_max))
                        e.append(error)
        row["oracle_beta"] = rc.slope_through_origin(x, e) if x else None
        for L in rc.LENGTHS:
            usable = [v for v in runs if issued[v.replica]["alarm_sweep"] not in (None, "UNKNOWN")
                      and rc.eligible(issued[v.replica]["alarm_sweep"], L)]
            fit = [v for v in usable if v.replica in fitting]
            held = [v for v in usable if v.replica in evaluation]
            rows = {v.replica: features(v, issued[v.replica]["alarm_sweep"], constants) for v in usable}
            y = {v.replica: int(bool(endpoints[(v.condition, v.replica)][L]["valid"])) for v in usable}
            model = fit_and_predict([rows[v.replica] for v in fit], [y[v.replica] for v in fit],
                                    [rows[v.replica] for v in held],
                                    None if saved_models is None else saved_models.get(family, {}).get(str(L)))
            y_eval = [y[v.replica] for v in held]
            seed = rc.analysis_seed(namespace, condition, f"bootstrap-L{L:05d}")
            scored = (rc.paired_bootstrap(y_eval, model["predictions"], seed, resamples=resamples)
                      if held and model["q0"] is not None else {"n": len(held), "models": {}})
            calibration_table = ({m: rc.calibration_bins(q, y_eval) for m, q in model["predictions"].items()}
                                 if held and model["q0"] is not None else {})
            oracle = {"frozen_unavailable": None, "brier_stable": None}
            if row["oracle_beta"] is not None and out["sigma_eff"][L]:
                evaluated = [v for v in runs if v.replica in evaluation]
                delta = [row["oracle_beta"] * rc.s_bar_closed_form(L, v.t0, v.t_cap, v.s_max)
                         / (out["sigma_eff"][L] / math.sqrt(L)) for v in evaluated]
                frozen = [rc.validity_probability_frozen(d) for d in delta]
                stable = [rc.validity_probability_stable(d) for d in delta]
                truth = [int(bool(endpoints[(v.condition, v.replica)][L]["valid"])) for v in evaluated]
                oracle = {"runs": len(evaluated), "frozen_unavailable": sum(1 for q in frozen if not math.isfinite(q)),
                          "brier_frozen": rc.brier(frozen, truth), "brier_stable": rc.brier(stable, truth),
                          "max_abs_delta": max(map(abs, delta), default=None)}
            row["per_length"][L] = {"eligible": len(usable), "eligible_fitting": len(fit), "eligible_evaluation": len(held),
                                    "valid_among_eligible": sum(y.values()), "q0": model["q0"], "models": model["models"],
                                    "fallbacks": model["fallbacks"], "scores": scored, "calibration": calibration_table,
                                    "oracle": oracle, "eligible_evaluation_ids": [v.replica for v in held],
                                    "predictions": model["predictions"], "reasons": model["reasons"], "y": y_eval}
            if L == 20480:
                bootstrap[family] = scored
        families[family] = row
    out["families"] = families
    c23 = rc.criteria_2_3(criterion_input)
    out["criteria_2_3"] = c23
    c4 = rc.criterion_4(bootstrap)
    out["criterion_4"] = c4
    gate = out["V2_D"]
    interpretable = bool(gate and gate["pass"] and gate["unknown"] == 0)
    out["issuance_interpretable"] = interpretable
    out["verdicts"] = {"criterion_1": out["criterion_1"]["result"] if out["criterion_1"] else None,
                       "criterion_2": c23["criterion_2"] if interpretable else "UNINTERPRETABLE",
                       "criterion_3": c23["criterion_3"] if interpretable else "UNINTERPRETABLE",
                       "criterion_4": c4["result"] if interpretable else "UNINTERPRETABLE",
                       "criterion_5": "INSUFFICIENT_INFORMATION"}
    out["status"] = "RECOUNTED"
    return out


# ---------------------------------------------------------------------------
# Comparison with the pilot's own analysis outputs (protocol.json implementation_checks.recount)
# ---------------------------------------------------------------------------
def saved_models_from_primary(models: Mapping[str, Any]) -> dict[str, Any]:
    """Parameters of results/models.json in the form fit_and_predict takes. The recount does not refit."""
    out: dict[str, Any] = {}
    for family, cells in models["families"].items():
        out[family] = {}
        for length, cell in cells.items():
            out[family][length] = {name: ({"theta": list(m["w"]) + [m["b"]], "mu": m["mu"], "scale": m["scale"]}
                                          if m["status"] == "FITTED" else {"fallback": m["fallback_reason"]})
                                   for name, m in cell["models"].items()}
    return out


def _jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def compare_with_primary(root: Path, independent: Mapping[str, Any]) -> dict[str, Any]:
    """Every quantity of the pilot's result files that the independent recount also has.

    Flags, counts, labels, identifiers and null patterns must be equal; finite floats must agree to
    1e-12 relative and absolute. A file that does not exist yet is listed, not treated as agreement.
    """
    results = Path(root) / "results"
    issues: list[dict[str, Any]] = []
    compared = 0
    absent: list[str] = []

    def check(where: str, name: str, stored: Any, mine: Any) -> None:
        nonlocal compared
        compared += 1
        if not rc.compare(stored, mine):
            issues.append({"scope": "campaign", "id": where, "field": name, "stored": _brief(stored), "independent": _brief(mine)})

    def read(name: str) -> Any:
        path = results / name
        if not path.is_file():
            absent.append(name)
            return None
        return _jsonl(path) if name.endswith(".jsonl") else json.loads(path.read_text(encoding="utf-8"))

    calibration = read("calibration.json")
    if calibration is not None:
        check("calibration", "threshold_full", calibration["threshold_full"], independent["thresholds"]["full"])
        check("calibration", "threshold_prime", calibration["threshold_prime"], independent["thresholds"]["prime"])
        check("calibration", "application_constants", calibration["application_constants"],
              independent.get("application_constants"))
        for length in rc.LENGTHS:
            check("calibration", f"sigma_eff[{length}]", calibration["sigma_eff"][str(length)],
                  (independent.get("sigma_eff") or {}).get(length))
    if independent.get("status") != "RECOUNTED":
        return {"compared": compared, "discrepancies": issues, "files_absent": absent, "agrees": not issues and not absent,
                "note": "calibration is incomplete: nothing downstream exists to compare"}
    healthy = read("healthy.json")
    if healthy is not None:
        for length in rc.LENGTHS:
            mine = independent["criterion_1"]["per_length"][length]
            check("healthy", f"coverage[{length}].covered", healthy["coverage"][str(length)]["covered"], mine["covered"])
            check("healthy", f"coverage[{length}].cp95", healthy["coverage"][str(length)]["cp95"], mine["cp95"])
        gate = independent["V2_D"]
        check("healthy", "V2_D.alarms", healthy["V2_D"]["alarms"], gate["alarms"])
        check("healthy", "V2_D.unknown", healthy["V2_D"]["unknown"], gate["unknown"])
        check("healthy", "V2_D.cp95", healthy["V2_D"]["cp95"], gate["cp95"])
        check("healthy", "V2_D.result", healthy["V2_D"]["result"],
              "INSUFFICIENT_INFORMATION" if gate["unknown"] else "PASS" if gate["pass"] else "FAIL")
    issuance = read("issuance.jsonl")
    if issuance is not None:
        stored_rows = {(r["family"], r["replica"]): r for r in issuance if r["family"] != "healthy"}
        for family, row in independent["families"].items():
            for run in row["per_run"]:
                stored = stored_rows.get((family, run["replica"]))
                where = f"issuance {family}/{run['replica']:05d}"
                if stored is None:
                    issues.append({"scope": "campaign", "id": where, "field": "row", "stored": None, "independent": "present"})
                    continue
                unknown = run["alarm_sweep"] == "UNKNOWN"
                check(where, "status", stored["full"]["status"],
                      "UNKNOWN" if unknown else "NO_ALARM" if run["alarm_sweep"] is None else "ALARM")
                check(where, "alarm_sweep", stored["full"]["alarm_sweep"], None if unknown else run["alarm_sweep"])
                check(where, "class", stored["full"]["onset_classification"], run["class"])
                check(where, "alarm_minus_t0", stored["full"]["alarm_minus_t0"], run["alarm_minus_t0"])
                if not unknown and run["alarm_sweep"] is not None:
                    check(where, "already_invalid", stored["full"]["issuance_endpoint"]["already_invalid"], run["already_invalid"])
                for length in rc.LENGTHS:
                    check(where, f"valid[{length}]", stored["valid"][str(length)], run["valid"][length])
                    check(where, f"eligible[{length}]", stored["full"]["eligibility"][str(length)]["eligible"],
                          None if unknown else rc.eligible(run["alarm_sweep"], length))
    models = read("models.json")
    if models is not None:
        for family, row in independent["families"].items():
            check(f"models {family}", "oracle.beta", models["oracle"]["beta"][family], row["oracle_beta"])
            for length in rc.LENGTHS:
                cell, mine = models["families"][family][str(length)], row["per_length"][length]
                check(f"models {family} L={length}", "eligible_fitting", cell["eligible_fitting"], mine["eligible_fitting"])
                check(f"models {family} L={length}", "q0", cell["q0"], mine["q0"])
    predictions = read("predictions.jsonl")
    if predictions is not None:
        for family, row in independent["families"].items():
            for length in rc.LENGTHS:
                mine = row["per_length"][length]
                stored = sorted((p for p in predictions if p["family"] == family and p["L"] == length and p["eligible"]),
                                key=lambda p: p["replica"])
                where = f"predictions {family} L={length}"
                check(where, "eligible evaluation ids", [p["replica"] for p in stored], mine["eligible_evaluation_ids"])
                if [p["replica"] for p in stored] != mine["eligible_evaluation_ids"]:
                    continue
                check(where, "y", [p["y"] for p in stored], mine["y"])
                for model in rc.MODELS:
                    check(where, f"q[{model}]", [p["predictions"][model]["q"] for p in stored], mine["predictions"][model])
                    check(where, f"fallback[{model}]", [p["predictions"][model]["fallback"] for p in stored], mine["reasons"][model])
    scores = read("scores.json")
    if scores is not None:
        for family, row in independent["families"].items():
            for length in rc.LENGTHS:
                cell, mine = scores["families"][family][str(length)], row["per_length"][length]
                where = f"scores {family} L={length}"
                check(where, "eligible", cell["counts"]["eligible"], mine["eligible_evaluation"])
                # Under the frozen rule the oracle is unavailable where scipy.stats.nct returns NaN. Under
                # Amendment A01 it is available wherever its parameters are, and its values are compared
                # with the quadrature form (they agree to better than 1e-14 in the Brier score).
                amended = scores.get("oracle_rule") == "A01"
                check(where, "oracle.unavailable", cell["oracle_all_evaluation_runs"]["unavailable"],
                      0 if amended and mine["oracle"].get("runs") else mine["oracle"]["frozen_unavailable"])
                check(where, "oracle.brier", cell["oracle_all_evaluation_runs"]["brier"],
                      mine["oracle"].get("brier_stable" if amended else "brier_frozen"))
                for model in rc.MODELS:
                    stored_model, my_model = cell["models"][model], mine["scores"]["models"].get(model)
                    if my_model is None:
                        check(where, f"{model}.brier", stored_model["brier"], None)
                        continue
                    check(where, f"{model}.brier", stored_model["brier"], my_model["brier"])
                    check(where, f"{model}.skill", stored_model["skill"], my_model["skill"])
                    check(where, f"{model}.brier_interval_95", stored_model.get("brier_interval_95"), my_model["brier_ci95"])
                    check(where, f"{model}.skill_interval_95", stored_model.get("skill_interval_95"), my_model["skill_ci95"])
                    check(where, f"{model}.skill_interval_criterion_4", stored_model.get("skill_interval_criterion_4"),
                          my_model["skill_ci_criterion"])
                    check(where, f"{model}.calibration", [[b["count"], b["mean_predicted"], b["observed_valid_fraction"]]
                                                           for b in stored_model["calibration"]],
                          [[b["count"], b["mean_predicted"], b["observed_fraction"]] for b in mine["calibration"][model]])
    criteria = read("criteria.json")
    if criteria is not None:
        before = criteria["results_before_precedence"]
        check("criteria", "1", before["1"], independent["criterion_1"]["result"])
        check("criteria", "2", before["2"], independent["criteria_2_3"]["criterion_2"])
        check("criteria", "3", before["3"], independent["criteria_2_3"]["criterion_3"])
        check("criteria", "4", before["4"], independent["criterion_4"]["result"])
        for key in ("1", "2", "3", "4"):
            check("criteria", f"final {key}", criteria["results"][key], independent["verdicts"][f"criterion_{key}"])
        for family, row in independent["criteria_2_3"]["per_family"].items():
            stored = criteria["criteria_2_3"]["per_family"][family]
            for name in ("eligible", "fraction_valid", "fraction_already_invalid", "median_remaining"):
                check(f"criteria {family}", name, stored[name], row[name])
    return {"compared": compared, "discrepancies": issues, "files_absent": absent, "agrees": not issues and not absent,
            "tolerance": {"rtol": 1e-12, "atol": 1e-12}}


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def _verify_one(job: tuple[str, str, int, str, bool, int]) -> dict[str, Any]:
    root, condition, replica, namespace, deep, stride = job
    try:
        record, arrays = load_stored(Path(root), condition, replica)
        return verify_run(record, arrays, namespace, deep=deep, checkpoint_stride=stride)
    except FileNotFoundError as error:
        return {"id": f"{condition}/{replica:05d}", "compared": 0, "deep": False, "missing": str(error), "issues": []}
    except Exception as error:  # a file that cannot be read or parsed is itself a finding
        return {"id": f"{condition}/{replica:05d}", "compared": 0, "deep": False,
                "issues": [{"scope": "run", "id": f"{condition}/{replica:05d}", "field": "load",
                            "stored": None, "independent": f"{type(error).__name__}: {error}"}]}


def deep_sample(namespace: str, plan: Mapping[str, tuple[str, int]] = PLAN) -> dict[str, list[int]]:
    """The I-4 replay sample of the frozen protocol: five planned replicas per condition."""
    chosen = {}
    for condition, (_, count) in plan.items():
        rng = np.random.default_rng(rc.analysis_seed(namespace, condition, "I4-replay-sample"))
        chosen[condition] = np.sort(rng.choice(np.arange(count, dtype=np.int64), size=min(5, count), replace=False)).tolist()
    return chosen


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--namespace", required=True)
    parser.add_argument("--deep", choices=("all", "sample", "none"), default="sample")
    parser.add_argument("--checkpoint-stride", type=int, default=8)
    parser.add_argument("--replicas", type=int, default=None, help="verify only replicas below this number (smoke: 2)")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--no-campaign", action="store_true", help="skip the campaign recount (smoke namespaces)")
    parser.add_argument("--compare-primary", action="store_true",
                        help="recount with the saved model parameters and compare with the pilot's result files")
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args(argv)
    plan = PLAN if args.replicas is None else {c: (f, min(n, args.replicas)) for c, (f, n) in PLAN.items()}
    sample = deep_sample(args.namespace, plan)
    jobs = [(str(args.root), condition, replica, args.namespace,
             args.deep == "all" or (args.deep == "sample" and replica in sample[condition]), args.checkpoint_stride)
            for condition, (_, count) in plan.items() for replica in range(count)]
    if args.workers > 1:
        with ProcessPoolExecutor(args.workers) as pool:
            runs = list(pool.map(_verify_one, jobs, chunksize=4))
    else:
        runs = [_verify_one(job) for job in jobs]
    issues = [i for r in runs for i in r["issues"]]
    report: dict[str, Any] = {
        "schema": "rnglab-v019-audit-verification-v1", "root": str(args.root), "namespace": args.namespace,
        "protocol_sha256": file_sha256(rc.PROTOCOL_PATH),
        "code_sha256": {name: file_sha256(_HERE / name) for name in ("recount.py", "replay.py", "verify.py")},
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": __import__("scipy").__version__, "platform": platform.platform(),
                        "libc": " ".join(platform.libc_ver())},
        "runs_planned": len(jobs), "runs_missing": [r["id"] for r in runs if "missing" in r],
        "runs_deep": sum(1 for r in runs if r["deep"]), "deep_sample": sample, "comparisons": sum(r["compared"] for r in runs),
        "near_boundary_endpoints": {r["id"]: r["near_boundary_endpoints"] for r in runs if r.get("near_boundary_endpoints")},
        "discrepancies": issues, "discrepancy_count": len(issues)}
    if not args.no_campaign and not report["runs_missing"]:
        views: dict[str, list[RunView]] = {}
        for condition, (_, count) in plan.items():
            views[condition] = [load_view(args.root, condition, replica) for replica in range(count)]
        saved = None
        models_path = Path(args.root) / "results/models.json"
        if args.compare_primary and models_path.is_file():
            saved = saved_models_from_primary(json.loads(models_path.read_text(encoding="utf-8")))
        report["campaign"] = recount_campaign(views, args.namespace, plan=plan, saved_models=saved)
        if args.compare_primary:
            report["primary_comparison"] = compare_with_primary(args.root, report["campaign"])
            issues = issues + report["primary_comparison"]["discrepancies"]
            report["discrepancies"], report["discrepancy_count"] = issues, len(issues)
    else:
        report["campaign"] = {"status": "NOT_RUN", "reason": "disabled" if args.no_campaign else "runs missing"}
    report["agrees"] = (not issues and not report["runs_missing"]
                        and not (args.compare_primary and report.get("primary_comparison", {}).get("files_absent")))
    text = json.dumps(report, sort_keys=True, indent=2, allow_nan=False, default=_brief) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
    print(f"{len(jobs)} planned runs, {len(report['runs_missing'])} missing, {report['runs_deep']} regenerated, "
          f"{report['comparisons']} comparisons, {len(issues)} discrepancies")
    return 0 if report["agrees"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
