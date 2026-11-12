"""Stage D pilot: feasibility forecasters, oracle reference, predictions and scores.

protocol.json `forecasters`, `oracle` and `scores` govern; expressions quoted there are used as
written. Two stages, run as two processes:

    fit       opens calibration results and fitting runs only, writes results/models.json and
              results/model_manifest.json (the hash of the model file and the list of inputs);
    evaluate  checks that hash, opens the evaluation runs once, writes results/predictions.jsonl
              and results/scores.json.

No model is refitted, tuned or replaced after an evaluation run has been opened. The fit stage
refuses to open any run that is not a fitting run.

Oracle. The probability is computed by the frozen call. That call can return NaN
(docs/v020/F01_NUMERICAL_RESOLUTION.md); such a value is reported unavailable, as the protocol
says. If the author approves Amendment A01 (pilot_common.approved_amendments), a non-finite
lower term at positive noncentrality is taken as 0; nothing else changes.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import scipy.optimize
import scipy.special
import scipy.stats

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))
import pilot_analysis as pa  # noqa: E402
import pilot_common as pc  # noqa: E402

FEATURES = {"F0": [], "F0+": ["alarm_sweep"], "F-B-lite": ["alarm_sweep", "target_at_alarm", "target_slope"],
            "F-AB-lite": ["alarm_sweep", "target_at_alarm", "target_slope", "A1_over_sigma_D", "A2_standardised"]}
BATTERY_FEATURES = ("alarm_sweep", "target_at_alarm", "target_slope")
OPTIMISER = {"method": "L-BFGS-B", "options": {"gtol": 1e-4, "ftol": 64 * float(np.finfo(np.float64).eps), "maxiter": 100,
                                               "maxfun": 15000, "maxls": 50, "maxcor": 10}}
EDGES = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]


# ---------------------------------------------------------------------------
# Features (protocol.json forecasters.features_in_order, slope, application_alignment)
# ---------------------------------------------------------------------------
def features(context: pc.Context, record: Mapping[str, Any], alarm_sweep: int, constants: Mapping[str, float] | None) -> dict:
    """Inputs of the forecasters for one run, from checkpoints at or before the alarm only."""
    target = pc.TARGETS[record["family"]]
    battery = record["conventional"]
    ends = [int(c["end"]) for c in battery]
    j = ends.index(int(alarm_sweep))
    values = [pc.finite_or_none(c.get(target)) for c in battery[max(0, j - 3):j + 1]]  # the last min(4, j+1), alarm included
    count = len(values)
    if any(v is None for v in values):
        slope = None
    elif count == 1:
        slope = 0.0
    else:
        x = np.asarray(ends[max(0, j - 3):j + 1], dtype=np.float64)
        v = np.asarray(values, dtype=np.float64)
        slope = float(np.sum((x - np.mean(x)) * (v - np.mean(v))) / np.sum((x - np.mean(x))**2))
    out: dict[str, Any] = {"alarm_sweep": float(alarm_sweep), "target_at_alarm": values[-1], "target_slope": slope,
                           "slope_count": count, "A1_over_sigma_D": None, "A2_standardised": None,
                           "application_end": None, "application_reason": None}
    application = [c for c in record["application"] if int(c["end"]) <= int(alarm_sweep)]
    if not application:
        out["application_reason"] = "NO_APPLICATION_CHECKPOINT"
        return out
    last = max(application, key=lambda c: int(c["end"]))
    standardised = context.measure.standardise_application(last, dict(constants or {}))
    out.update(application_end=int(last["end"]), A1_over_sigma_D=pc.finite_or_none(standardised.get("A1_over_sigma_D")),
               A2_standardised=pc.finite_or_none(standardised.get("A2_standardised")))
    if out["A1_over_sigma_D"] is None or out["A2_standardised"] is None:
        out["application_reason"] = "MISSING_APPLICATION_FEATURE"
    return out


# ---------------------------------------------------------------------------
# Logistic model (protocol.json forecasters.logistic, standardisation)
# ---------------------------------------------------------------------------
def fit_logistic(z: np.ndarray, y: np.ndarray) -> dict:
    """The frozen objective, gradient, call and convergence test."""
    d = z.shape[1]

    def objective(theta: np.ndarray) -> float:
        w, b = theta[:d], theta[d]
        eta = z @ w + b
        return float(np.sum(np.logaddexp(0.0, eta) - y * eta, dtype=np.float64) + 0.5 * np.dot(w, w))

    def gradient(theta: np.ndarray) -> np.ndarray:
        w, b = theta[:d], theta[d]
        residual = scipy.special.expit(z @ w + b) - y
        return np.concatenate([z.T @ residual + w, [np.sum(residual, dtype=np.float64)]])

    out: dict[str, Any] = {"accepted": False, "w": None, "b": None, "exception": None}
    try:
        result = scipy.optimize.minimize(fun=objective, x0=np.zeros(d + 1, dtype=np.float64), jac=gradient,
                                         method="L-BFGS-B", bounds=None, options=dict(OPTIMISER["options"]))
    except Exception as error:  # an exception is nonconvergence
        return {**out, "exception": {"class": type(error).__name__, "message": str(error)}}
    message = result.message.decode() if isinstance(result.message, bytes) else str(result.message)
    value, grad = objective(result.x), gradient(result.x)
    norm = float(np.max(np.abs(grad)))
    finite = bool(np.isfinite(result.x).all() and math.isfinite(float(result.fun)) and np.isfinite(result.jac).all()
                  and math.isfinite(value) and np.isfinite(grad).all())
    gradient_test = finite and norm <= 1e-4
    reduction_test = "RELATIVE REDUCTION OF F" in message.upper().replace("_", " ").replace("REL REDUCTION", "RELATIVE REDUCTION")
    accepted = bool(result.success is True and result.status == 0 and result.nit <= 100 and finite
                    and (gradient_test or reduction_test))
    out.update(accepted=accepted, success=bool(result.success), status=int(result.status), message=message,
               nit=int(result.nit), nfev=int(result.nfev), njev=int(getattr(result, "njev", result.nfev)),
               objective=value if finite else None, gradient_infinity_norm=norm if finite else None,
               terminated_by=("gradient" if gradient_test else "relative_reduction" if reduction_test else "other"))
    if accepted:
        out.update(w=[float(v) for v in result.x[:d]], b=float(result.x[d]))
    return out


def fit_model(name: str, rows: list[dict], labels: list[int], ids: list[int]) -> dict:
    """One logistic model on the rows whose required features are all finite; otherwise a fallback to F0."""
    columns = FEATURES[name]
    usable = [i for i, row in enumerate(rows) if all(row[c] is not None and math.isfinite(row[c]) for c in columns)]
    y = np.asarray([labels[i] for i in usable], dtype=np.float64)
    model: dict[str, Any] = {"features": columns, "training_rows": len(usable),
                             "training_classes": [int(np.sum(y == 0)), int(np.sum(y == 1))],
                             "excluded_ids": [ids[i] for i in range(len(rows)) if i not in set(usable)],
                             "status": "FALLBACK_F0", "fallback_reason": None, "mu": None, "scale": None,
                             "constant_feature": None, "w": None, "b": None, "optimiser": None}
    if not usable:
        return {**model, "fallback_reason": "EMPTY_SUBSET"}
    if y.min() == y.max():
        return {**model, "fallback_reason": "ONE_CLASS"}
    x = np.asarray([[rows[i][c] for c in columns] for i in usable], dtype=np.float64)
    mu = np.array([np.mean(x[:, j], dtype=np.float64) for j in range(x.shape[1])])
    sd = np.array([np.std(x[:, j], ddof=0, dtype=np.float64) for j in range(x.shape[1])])
    scale = np.where(sd == 0, 1.0, sd)
    fitted = fit_logistic((x - mu) / scale, y)
    model.update(mu=mu.tolist(), scale=scale.tolist(), constant_feature=[bool(v == 0) for v in sd],
                 optimiser={k: v for k, v in fitted.items() if k not in ("w", "b")})
    if not fitted["accepted"]:
        return {**model, "fallback_reason": "NOT_CONVERGED"}
    return {**model, "status": "FITTED", "w": fitted["w"], "b": fitted["b"]}


def predict(cell: Mapping[str, Any], row: Mapping[str, Any]) -> dict[str, dict]:
    """Probabilities of all four forecasters for one eligible run, with the reason for any fallback
    (protocol.json forecasters.training_and_missing_features)."""
    q0 = cell["q0"]
    if q0 is None:
        return {name: {"q": None, "fallback": "NO_ELIGIBLE_FITTING_RUNS"} for name in pc.MODELS}

    def logistic(name: str) -> tuple[float, str | None]:
        model = cell["models"][name]
        if any(row[c] is None or not math.isfinite(row[c]) for c in BATTERY_FEATURES if c in model["features"]):
            return q0, "MISSING_FEATURE"
        if model["status"] != "FITTED":
            return q0, model["fallback_reason"]
        z = (np.asarray([row[c] for c in model["features"]], dtype=np.float64) - np.asarray(model["mu"])) / np.asarray(model["scale"])
        return float(scipy.special.expit(z @ np.asarray(model["w"], dtype=np.float64) + model["b"])), None

    out = {"F0": {"q": q0, "fallback": None}}
    for name in ("F0+", "F-B-lite"):
        q, reason = logistic(name)
        out[name] = {"q": q, "fallback": reason}
    if row["application_reason"] is not None:  # no application checkpoint, or an unavailable application signal
        out["F-AB-lite"] = {"q": out["F-B-lite"]["q"], "fallback": row["application_reason"]}
    else:
        q, reason = logistic("F-AB-lite")
        out["F-AB-lite"] = {"q": q, "fallback": reason}
    return out


# ---------------------------------------------------------------------------
# Oracle (protocol.json oracle)
# ---------------------------------------------------------------------------
def s_bar(severity: np.ndarray, length: int) -> float:
    return float(np.sum(severity[:length], dtype=np.float64) / length)


def oracle_probability(beta: float | None, sigma_eff: float | None, mean_severity: float, length: int,
                       rule: str = "frozen") -> dict:
    """q = P(valid at L) from the frozen noncentral-t call. rule 'A01' applies the approved amendment."""
    if beta is None or sigma_eff is None:
        return {"q": None, "available": False, "reason": "ORACLE_PARAMETER_UNAVAILABLE", "delta": None, "lower_term_replaced": False}
    se = sigma_eff / math.sqrt(length)
    delta = beta * mean_severity / se
    c = scipy.stats.t.ppf(0.975, df=15)
    with np.errstate(all="ignore"):
        upper = float(scipy.stats.nct.cdf(c, df=15, nc=delta))
        lower = float(scipy.stats.nct.cdf(-c, df=15, nc=delta))
    q, replaced = float(upper - lower), False
    if rule == "A01" and not math.isfinite(lower) and delta > 0 and math.isfinite(upper):
        q, replaced = upper, True
    if not math.isfinite(q) or q < 0 or q > 1:
        return {"q": None, "available": False, "reason": "NONFINITE_OR_OUT_OF_RANGE", "delta": pc.finite_or_none(delta),
                "lower_term_replaced": False}
    return {"q": q, "available": True, "reason": None, "delta": float(delta), "lower_term_replaced": replaced}


def oracle_rule() -> str:
    return "A01" if "A01" in pc.approved_amendments() else "frozen"


# ---------------------------------------------------------------------------
# Fit stage
# ---------------------------------------------------------------------------
def _issued(context: pc.Context, record: Mapping[str, Any], arrays: Any, threshold: float | None) -> dict:
    if arrays is None or threshold is None:
        return {"status": "UNKNOWN", "alarm_sweep": None}
    return context.measure.issuance(record["conventional"], threshold)


def fit(context: pc.Context) -> dict:
    """Fit every forecaster and the oracle parameters on fitting runs; write the model file and its manifest."""
    pa._require_science_plan(context)
    plan = context.plan
    calibration, calibration_sha = pa.read_calibration(context)
    inputs = [{"path": "results/calibration.json", "sha256": calibration_sha}]
    constants, tau = calibration["application_constants"], calibration["threshold_full"]
    out: dict[str, Any] = {"schema": "rnglab-v019-d-pilot-models-v1", "protocol_sha256": pc.protocol_sha256(),
                           "design_sha256": pc.sha256_file(HERE / "design.json"), "calibration_sha256": calibration_sha,
                           "feature_order": FEATURES, "targets": pc.TARGETS, "optimiser": OPTIMISER, "families": {},
                           "oracle": {"rule": oracle_rule(), "sigma_eff": calibration["sigma_eff"], "beta": {}, "omitted_pairs": []}}
    for condition, family, _ in plan.conditions:
        if family == "healthy":
            continue
        rows: dict[int, list] = {length: [] for length in pc.LENGTHS}
        unknown, xs, es = 0, [], []
        for replica in plan.fitting:
            if plan.split(condition, replica) != "fitting":
                raise pc.PilotError("evaluation barrier: the fit stage was asked to open a run that is not a fitting run")
            record, arrays = pa.load_run(context, condition, replica)
            for relative in pc.run_paths(condition, replica):
                path = context.root / relative
                if path.exists():
                    inputs.append({"path": relative.as_posix(), "sha256": pc.sha256_file(path)})
            for length in pc.LENGTHS:  # oracle slope: every fitting run and length, whatever the alarm
                error = pc.finite_or_none(((record.get("endpoints") or {}).get(str(length)) or {}).get("signed_error"))
                if arrays is None or error is None:
                    out["oracle"]["omitted_pairs"].append({"condition": condition, "replica": replica, "L": length,
                                                           "reason": "signed error unavailable"})
                    continue
                xs.append(s_bar(arrays["severity"], length))
                es.append(error)
            issued = _issued(context, record, arrays, tau)
            if issued["status"] == "UNKNOWN":
                unknown += 1
                continue
            alarm = issued["alarm_sweep"]
            for length in pc.LENGTHS:
                if alarm is not None and alarm < length:
                    label = int(bool(record["endpoints"][str(length)]["valid"]))
                    rows[length].append((replica, features(context, record, alarm, constants), label))
        x, e = np.asarray(xs, dtype=np.float64), np.asarray(es, dtype=np.float64)
        denominator = float(np.sum(x * x)) if len(x) else 0.0
        out["oracle"]["beta"][family] = float(np.sum(x * e) / denominator) if denominator > 0 else None
        out["families"][family] = {}
        for length in pc.LENGTHS:
            ids = [r for r, _, _ in rows[length]]
            feature_rows = [f for _, f, _ in rows[length]]
            labels = [y for _, _, y in rows[length]]
            cell: dict[str, Any] = {"eligible_fitting": len(ids), "fitting_ids": ids,
                                    "classes": [labels.count(0), labels.count(1)], "issuance_unknown_fitting": unknown,
                                    "q0": (sum(labels) / len(labels)) if labels else None, "models": {}}
            for name in ("F0+", "F-B-lite", "F-AB-lite"):
                cell["models"][name] = fit_model(name, feature_rows, labels, ids)
            out["families"][family][str(length)] = cell
    digest = pc.atomic_write(context.result("models.json"), pc.canonical_json(out))
    manifest = {"schema": "rnglab-v019-d-pilot-model-manifest-v1", "models_sha256": digest,
                "written_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "fitting_inputs": inputs,
                "evaluation_barrier": {"evaluation_runs_opened": 0,
                                       "statement": "the fit stage opened calibration results and fitting runs only; "
                                                    "evaluation runs are first opened by the evaluate stage, after this file"}}
    pc.atomic_write(context.result("model_manifest.json"), pc.canonical_json(manifest))
    return out


# ---------------------------------------------------------------------------
# Scores (protocol.json scores)
# ---------------------------------------------------------------------------
def brier(q: np.ndarray, y: np.ndarray) -> float:
    return float(np.mean((q - y)**2, dtype=np.float64))


def calibration_table(context: pc.Context, q: np.ndarray, y: np.ndarray) -> list[dict]:
    index = np.minimum(np.searchsorted(EDGES, q, side="right") - 1, 4)
    rows = []
    for b in range(5):
        inside = index == b
        count, hits = int(inside.sum()), int(y[inside].sum())
        rows.append({"bin": b, "low": EDGES[b], "high": EDGES[b + 1], "count": count,
                     "mean_predicted": float(np.mean(q[inside], dtype=np.float64)) if count else None,
                     "observed_valid_fraction": hits / count if count else None, "cp95": pa.cp(context, hits, count)})
    return rows


def _interval(values: np.ndarray, level: float) -> list[float]:
    return np.quantile(values, [(1 - level) / 2, 1 - (1 - level) / 2], method="linear").tolist()


def score_cell(context: pc.Context, family: str, length: int, y: np.ndarray, forecasts: Mapping[str, list]) -> dict:
    """Scores of the four forecasters on the eligible evaluation runs of one family and length, in replica order."""
    from rnglab.seeds import namespace_seed
    n = len(y)
    out: dict[str, Any] = {"n": n, "models": {}, "contrasts": {}, "bootstrap": None}
    missing = {name: sum(1 for q in forecasts[name] if q is None) for name in pc.MODELS}
    complete = n > 0 and not any(missing.values())
    for name in pc.MODELS:
        available = [(q, label) for q, label in zip(forecasts[name], y) if q is not None]
        row: dict[str, Any] = {"missing": missing[name], "brier": None, "skill": None, "calibration": None,
                               "available_case_brier_descriptive": None}
        if n and not missing[name]:
            q = np.asarray(forecasts[name], dtype=np.float64)
            row.update(brier=brier(q, y), calibration=calibration_table(context, q, y))
        elif available:
            row["available_case_brier_descriptive"] = brier(np.asarray([a for a, _ in available], dtype=np.float64),
                                                            np.asarray([b for _, b in available], dtype=np.float64))
        out["models"][name] = row
    if not complete:
        out["status"] = "INSUFFICIENT_INFORMATION"
        return out
    q = {name: np.asarray(forecasts[name], dtype=np.float64) for name in pc.MODELS}
    base = out["models"]["F0+"]["brier"]
    for name in pc.MODELS:
        out["models"][name]["skill"] = None if base == 0 else 1 - out["models"][name]["brier"] / base
    fb, fab = out["models"]["F-B-lite"]["brier"], out["models"]["F-AB-lite"]["brier"]
    out["contrasts"] = {"FAB_versus_FB": {"delta_BS": fb - fab, "relative_skill": None if fb == 0 else 1 - fab / fb}}
    seed = namespace_seed(context.plan.namespace, "fault/" + family, "analysis", f"bootstrap-L{length:05d}")
    rng = np.random.default_rng(seed)
    loss = {name: (q[name] - y)**2 for name in pc.MODELS}
    resamples = 10000
    scores = {name: np.empty(resamples) for name in pc.MODELS}
    for b in range(resamples):
        idx = rng.integers(0, n, size=n)
        for name in pc.MODELS:
            scores[name][b] = np.mean(loss[name][idx], dtype=np.float64)
    undefined = int(np.count_nonzero(scores["F0+"] == 0))
    undefined_fb = int(np.count_nonzero(scores["F-B-lite"] == 0))
    out["bootstrap"] = {"seed": seed, "resamples": resamples, "undefined_skill_resamples": undefined,
                        "undefined_relative_skill_resamples": undefined_fb}
    for name in pc.MODELS:
        row = out["models"][name]
        row["brier_interval_95"] = _interval(scores[name], pc.OTHER_LEVEL)
        if undefined == 0:
            values = 1 - scores[name] / scores["F0+"]
            row["skill_interval_95"] = _interval(values, pc.OTHER_LEVEL)
            row["skill_interval_criterion_4"] = _interval(values, pc.CRITERION_4_LEVEL)
        else:
            row["skill_interval_95"] = row["skill_interval_criterion_4"] = None
    contrast = out["contrasts"]["FAB_versus_FB"]
    contrast["delta_BS_interval_95"] = _interval(scores["F-B-lite"] - scores["F-AB-lite"], pc.OTHER_LEVEL)
    contrast["relative_skill_interval_95"] = (_interval(1 - scores["F-AB-lite"] / scores["F-B-lite"], pc.OTHER_LEVEL)
                                              if undefined_fb == 0 else None)
    out["status"] = "SCORED"
    return out


def score_oracle(context: pc.Context, y: np.ndarray, rows: list[dict]) -> dict:
    q = [r["q"] for r in rows]
    unavailable = sum(1 for v in q if v is None)
    out = {"n": len(rows), "unavailable": unavailable, "lower_term_replaced": sum(1 for r in rows if r["lower_term_replaced"]),
           "brier": None, "calibration": None, "available_case_brier_descriptive": None}
    if rows and not unavailable:
        values = np.asarray(q, dtype=np.float64)
        out.update(brier=brier(values, y), calibration=calibration_table(context, values, y))
    elif len(rows) > unavailable:
        keep = [i for i, v in enumerate(q) if v is not None]
        out["available_case_brier_descriptive"] = brier(np.asarray([q[i] for i in keep], dtype=np.float64), y[keep])
    return out


# ---------------------------------------------------------------------------
# Evaluate stage
# ---------------------------------------------------------------------------
def evaluate(context: pc.Context) -> dict:
    """Open the evaluation runs once; write predictions and scores from the frozen models."""
    pa._require_science_plan(context)
    plan = context.plan
    manifest = pc.read_json(context.result("model_manifest.json"))
    models_sha = pc.sha256_file(context.result("models.json"))
    if models_sha != manifest["models_sha256"]:
        raise pc.PilotError("results/models.json does not match its manifest: evaluation refused")
    if context.result("predictions.jsonl").exists() or context.result("scores.json").exists():
        raise pc.PilotError("evaluation outputs already exist: the evaluation runs are opened once")
    models = pc.read_json(context.result("models.json"))
    calibration, calibration_sha = pa.read_calibration(context)
    if calibration_sha != models["calibration_sha256"]:
        raise pc.PilotError("calibration.json changed after the models were fitted")
    constants, tau, rule = calibration["application_constants"], calibration["threshold_full"], models["oracle"]["rule"]
    predictions: list[dict] = []
    scores: dict[str, Any] = {"schema": "rnglab-v019-d-pilot-scores-v1", "models_sha256": models_sha, "families": {},
                              "oracle_rule": rule}
    for condition, family, _ in plan.conditions:
        if family == "healthy":
            continue
        collected: dict[int, dict] = {length: {"y": [], "forecasts": {name: [] for name in pc.MODELS}, "oracle_all": [],
                                               "y_all": [], "oracle_eligible": [],
                                               "counts": {"planned": 0, "unknown": 0, "no_alarm": 0, "late_alarm": 0,
                                                          "eligible": 0}} for length in pc.LENGTHS}
        for replica in plan.evaluation:
            record, arrays = pa.load_run(context, condition, replica)
            issued = _issued(context, record, arrays, tau)
            for length in pc.LENGTHS:
                bucket = collected[length]
                bucket["counts"]["planned"] += 1
                valid = int(bool(((record.get("endpoints") or {}).get(str(length)) or {}).get("valid", False)))
                oracle = (oracle_probability(models["oracle"]["beta"][family], models["oracle"]["sigma_eff"][str(length)],
                                             s_bar(arrays["severity"], length), length, rule) if arrays is not None else
                          {"q": None, "available": False, "reason": "RUN_HAS_NO_TRACE", "delta": None, "lower_term_replaced": False})
                bucket["oracle_all"].append(oracle)
                bucket["y_all"].append(valid)
                row = {"family": family, "replica": replica, "L": length, "y": valid, "alarm_sweep": issued["alarm_sweep"],
                       "issuance_status": issued["status"], "eligible": False, "features": None, "predictions": None,
                       "oracle": oracle, "models_sha256": models_sha}
                if issued["status"] == "UNKNOWN":
                    bucket["counts"]["unknown"] += 1
                elif issued["alarm_sweep"] is None:
                    bucket["counts"]["no_alarm"] += 1
                elif issued["alarm_sweep"] >= length:
                    bucket["counts"]["late_alarm"] += 1
                else:
                    bucket["counts"]["eligible"] += 1
                    feature_row = features(context, record, issued["alarm_sweep"], constants)
                    predicted = predict(models["families"][family][str(length)], feature_row)
                    row.update(eligible=True, features=feature_row, predictions=predicted)
                    bucket["y"].append(valid)
                    bucket["oracle_eligible"].append(oracle)
                    for name in pc.MODELS:
                        bucket["forecasts"][name].append(predicted[name]["q"])
                predictions.append(row)
        scores["families"][family] = {}
        for length in pc.LENGTHS:
            bucket = collected[length]
            y = np.asarray(bucket["y"], dtype=np.float64)
            cell = score_cell(context, family, length, y, bucket["forecasts"])
            counts = bucket["counts"]
            counts["evaluable"] = counts["planned"] - counts["unknown"]
            cell["counts"] = counts
            cell["fallbacks"] = {name: {} for name in pc.MODELS}
            for row in predictions:
                if row["family"] == family and row["L"] == length and row["eligible"]:
                    for name in pc.MODELS:
                        reason = row["predictions"][name]["fallback"]
                        if reason:
                            cell["fallbacks"][name][reason] = cell["fallbacks"][name].get(reason, 0) + 1
            cell["oracle_all_evaluation_runs"] = score_oracle(context, np.asarray(bucket["y_all"], dtype=np.float64), bucket["oracle_all"])
            cell["oracle_eligible_subset"] = score_oracle(context, y, bucket["oracle_eligible"])
            scores["families"][family][str(length)] = cell
    pc.atomic_write(context.result("predictions.jsonl"), pc.canonical_jsonl(predictions))
    pc.atomic_write(context.result("scores.json"), pc.canonical_json(scores))
    return scores
