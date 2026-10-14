"""Small prospective RUS baseline with explicit censoring and abstention.

The one-sided split-conformal construction targets *marginal* lower-bound
coverage for exchangeable independent seed groups from one fixed population.
It gives no conditional safety guarantee, including among non-abstained cases.
Development censored labels are used as conservative surrogate regression
outcomes; calibration uses upper bounds on unavailable residuals. This makes
calibration conservative without silently selecting only observed failures.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.stats import beta


@dataclass(frozen=True)
class Observation:
    seed_id: str
    split: str
    application: str
    fault: str
    configuration: str
    issued_at: float
    features: tuple[float, ...]
    features_through: float
    event_at: float | None = None
    censored_at: float | None = None

    def __post_init__(self):
        if self.split not in {"development", "calibration", "final_holdout"}:
            raise ValueError("Unknown seed split")
        if not self.seed_id or not self.configuration:
            raise ValueError("Seed and configuration identities are required")
        if not math.isfinite(self.issued_at) or self.issued_at < 0:
            raise ValueError("Issuance must be a finite nonnegative time")
        if not math.isfinite(self.features_through) or self.features_through > self.issued_at:
            raise ValueError("Feature leakage: information after issuance")
        if self.features_through < 0 or not self.features:
            raise ValueError("Features need a valid past observation horizon")
        if not all(math.isfinite(float(v)) for v in self.features):
            raise ValueError("Features must be finite")
        if self.event_at is not None and self.censored_at is not None:
            raise ValueError("Supply an event OR a right censor, not both")
        for t in (self.event_at, self.censored_at):
            if t is not None and (not math.isfinite(t) or t < self.issued_at):
                raise ValueError("Outcome time cannot precede issuance")

    @property
    def remaining_observed(self) -> float:
        end = self.event_at if self.event_at is not None else self.censored_at
        if end is None:
            raise ValueError("Training/scoring needs an event or right censor")
        return float(end - self.issued_at)


@dataclass(frozen=True)
class Prediction:
    status: str
    lower_bound: float | None
    reason: str | None
    model_id: str


@dataclass(frozen=True)
class RUSPredictor:
    application: str
    fault: str
    configuration: str
    alpha: float
    coefficients: tuple[float, ...]
    feature_min: tuple[float, ...]
    feature_max: tuple[float, ...]
    residual_quantile: float
    development_seeds: tuple[str, ...]
    calibration_seeds: tuple[str, ...]
    development_censored: int
    calibration_censored: int
    model_id: str
    frozen: bool = True
    method: str = "log1p-ols/censor-conservative-one-sided-split-conformal-v1"

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "RUSPredictor":
        clean = dict(data)
        for key in ("coefficients", "feature_min", "feature_max", "development_seeds", "calibration_seeds"):
            clean[key] = tuple(clean[key])
        result = cls(**clean)
        if not result.frozen or _model_identity(result.to_dict()) != result.model_id:
            raise ValueError("Unfrozen or modified predictor")
        return result

    def save(self, path: str | Path) -> None:
        """Create a frozen model artifact; existing artifacts are never overwritten."""
        with Path(path).open("x") as stream:
            json.dump(self.to_dict(), stream, indent=2, allow_nan=False)
            stream.write("\n")


def _model_identity(data: dict) -> str:
    payload = {key: value for key, value in data.items() if key != "model_id"}
    return sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _validate_group(rows: list[Observation], split: str, population: tuple[str, str, str]) -> set[str]:
    seeds = set()
    dimensions = set()
    for row in rows:
        if row.split != split:
            raise ValueError(f"Expected {split} data")
        if (row.application, row.fault, row.configuration) != population:
            raise ValueError("Mixed or unsupported populations")
        if row.seed_id in seeds:
            raise ValueError("One issuance per independent seed group is required")
        seeds.add(row.seed_id)
        dimensions.add(len(row.features))
        row.remaining_observed
        if row.event_at is not None and row.event_at <= row.issued_at:
            raise ValueError("Already failed at issuance")
    if len(dimensions) > 1:
        raise ValueError("Inconsistent feature dimensions")
    return seeds


def train_predictor(
    development: Iterable[Observation], calibration: Iterable[Observation], *,
    application: str, fault: str, configuration: str, alpha: float = .1,
    min_development: int = 16, min_calibration: int = 20,
) -> RUSPredictor:
    """Fit on development only; calibrate once, then return an immutable model.

    For an event residual S = f(x)-log(1+T), censor followup C <= T yields
    S <= f(x)-log(1+C). The finite-sample upper order statistic of these upper
    residual bounds is conservative relative to the complete-label statistic.
    No event-time independence assumption for censoring is needed for this
    inequality; exchangeability of independent groups remains essential.
    """
    if not 0 < alpha < 1:
        raise ValueError("alpha must be in (0,1)")
    dev, cal = list(development), list(calibration)
    if len(dev) < min_development or len(cal) < min_calibration:
        raise ValueError("Insufficient independent development/calibration groups")
    population = (application, fault, configuration)
    ds = _validate_group(dev, "development", population)
    cs = _validate_group(cal, "calibration", population)
    if ds & cs:
        raise ValueError("Seed leakage between development and calibration")
    width = len(dev[0].features)
    if len(cal[0].features) != width:
        raise ValueError("Inconsistent feature dimensions")
    k = math.ceil((len(cal) + 1) * (1 - alpha))
    if k > len(cal):
        raise ValueError("Calibration too small for a nontrivial finite conformal bound")
    x = np.asarray([r.features for r in dev], dtype=float)
    design = np.column_stack([np.ones(len(dev)), x])
    # Censored development targets are lower-bound surrogates, not imputed events.
    y = np.log1p([r.remaining_observed for r in dev])
    coeff = np.linalg.lstsq(design, y, rcond=None)[0]
    cx = np.column_stack([np.ones(len(cal)), np.asarray([r.features for r in cal])])
    scores = cx @ coeff - np.log1p([r.remaining_observed for r in cal])
    q = float(np.sort(scores)[k - 1])
    data = dict(application=application, fault=fault, configuration=configuration,
                alpha=float(alpha), coefficients=tuple(float(v) for v in coeff),
                feature_min=tuple(float(v) for v in x.min(axis=0)),
                feature_max=tuple(float(v) for v in x.max(axis=0)),
                residual_quantile=q, development_seeds=tuple(sorted(ds)),
                calibration_seeds=tuple(sorted(cs)),
                development_censored=sum(r.event_at is None for r in dev),
                calibration_censored=sum(r.event_at is None for r in cal),
                model_id="", frozen=True,
                method="log1p-ols/censor-conservative-one-sided-split-conformal-v1")
    data["model_id"] = _model_identity(data)
    return RUSPredictor(**data)


def predict(model: RUSPredictor, row: Observation) -> Prediction:
    """Use only past features; future event/censor labels never affect the bound."""
    def abstain(reason):
        return Prediction("ABSTAIN", None, reason, model.model_id)
    if not model.frozen or _model_identity(model.to_dict()) != model.model_id:
        return abstain("UNFIT_OR_MODIFIED_MODEL")
    if (row.application, row.fault, row.configuration) != (model.application, model.fault, model.configuration):
        return abstain("UNSUPPORTED_POPULATION")
    if len(row.features) != len(model.feature_min):
        return abstain("UNSUPPORTED_FEATURE_DIMENSION")
    if row.event_at == row.issued_at:
        return abstain("ALREADY_FAILED")
    if any(v < lo or v > hi for v, lo, hi in zip(row.features, model.feature_min, model.feature_max)):
        return abstain("OUT_OF_DISTRIBUTION")
    value = model.coefficients[0] + sum(c * x for c, x in zip(model.coefficients[1:], row.features))
    if value - model.residual_quantile > 700:
        return abstain("NUMERICAL_RANGE")
    lower = max(0., math.expm1(value - model.residual_quantile))
    return Prediction("LOWER_BOUND", lower, None, model.model_id)


def exact_binomial_interval(events: int, total: int, confidence: float = .95) -> list[float] | None:
    """Two-sided pointwise Clopper-Pearson interval; no trials means unknown."""
    if total < 0 or events < 0 or events > total or not 0 < confidence < 1:
        raise ValueError("Invalid binomial count/confidence")
    if not total:
        return None
    tail = (1 - confidence) / 2
    lo = 0. if not events else float(beta.ppf(tail, events, total - events + 1))
    hi = 1. if events == total else float(beta.ppf(1 - tail, events + 1, total - events))
    return [lo, hi]


def validate(model: RUSPredictor, final_holdout: Iterable[Observation]) -> dict:
    """Score the frozen predictor once per independent final-holdout group.

    A censor with followup >= lower bound establishes safety; shorter followup
    is UNKNOWN. Risk bounds report worst/best cases across ALL emitted bounds,
    rather than treating unresolved censors as failures or silently dropping them.
    """
    rows = list(final_holdout)
    seen = set()
    outcomes = []
    for row in rows:
        if row.split != "final_holdout":
            raise ValueError("Validation requires final_holdout labels")
        if row.seed_id in seen or row.seed_id in model.development_seeds or row.seed_id in model.calibration_seeds:
            raise ValueError("Repeated or leaked final-holdout seed group")
        seen.add(row.seed_id)
        prediction = predict(model, row)
        outcome = "ABSTAINED"
        if prediction.status == "LOWER_BOUND":
            if row.event_at is not None:
                outcome = "UNSAFE" if prediction.lower_bound > row.remaining_observed else "SAFE"
            elif row.censored_at is not None:
                outcome = "SAFE" if row.remaining_observed >= prediction.lower_bound else "UNKNOWN"
            else:
                outcome = "UNKNOWN"
        outcomes.append({"seed_id": row.seed_id, "prediction": asdict(prediction), "outcome": outcome,
                         "event_at": row.event_at, "censored_at": row.censored_at})
    counts = {name: sum(r["outcome"] == name for r in outcomes) for name in ("SAFE", "UNSAFE", "UNKNOWN", "ABSTAINED")}
    emitted = len(rows) - counts["ABSTAINED"]
    unsafe_interval = exact_binomial_interval(counts["UNSAFE"], emitted)
    worst_interval = exact_binomial_interval(counts["UNSAFE"] + counts["UNKNOWN"], emitted)
    # This is an identification envelope: censorship uncertainty plus sampling CI.
    risk_envelope = None if not emitted else [unsafe_interval[0], worst_interval[1]]
    return {"model_id": model.model_id, "model_frozen": model.frozen,
            "population": {"application": model.application, "fault": model.fault, "configuration": model.configuration},
            "independent_seed_groups": len(rows), "emitted_bounds": emitted,
            "counts": counts, "abstention_fraction": counts["ABSTAINED"] / len(rows) if rows else None,
            "resolved_fraction_of_emitted": (counts["SAFE"] + counts["UNSAFE"]) / emitted if emitted else None,
            "unsafe_fraction_identification_bounds": [counts["UNSAFE"] / emitted, (counts["UNSAFE"] + counts["UNKNOWN"]) / emitted] if emitted else None,
            "unsafe_risk_95pct_identification_and_sampling_envelope": risk_envelope,
            "positive_bound_count": sum((r["prediction"]["lower_bound"] or 0) > 0 for r in outcomes),
            "censored_count": sum(r.censored_at is not None for r in rows),
            "missing_outcome_count": sum(r.event_at is None and r.censored_at is None for r in rows),
            "interpretation": "Exploratory marginal evaluation in one fixed population. No conditional risk guarantee; abstention can select a different risk population. Unknown censors remain explicit.",
            "outcomes": outcomes}
