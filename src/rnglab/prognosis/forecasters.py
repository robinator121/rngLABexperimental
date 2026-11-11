"""Forecasters of future validity and of the failure time, issued at an alarm.

Information rule. A forecaster receives a `History`: the monitoring statistics
up to and including the alarm checkpoint, and nothing else. `history_at`
copies that prefix out of a trajectory, so a forecaster holds no reference to
later values, to the severity path, or to any outcome. Parameters learnt from
labelled training runs (the bias slope, the noise scale, the sensitivity of
the statistic, logistic coefficients) are fixed by `fit` before any evaluation
history is seen.

Forecasters, in increasing use of the declared degradation family:

    BaseRate            the training fraction of valid runs (no run-specific information)
    AlarmTimeLogistic   logistic regression on the alarm sweep
    LiteLogistic        alarm sweep, statistic at the alarm, its slope over the last four checkpoints
    Persistence         the severity seen at the alarm stays as it is (no growth, no recovery)
    BayesRamp           posterior over the declared ramp parameters given the whole statistic
                        history, mixed over the validity model; also a failure-time distribution
    Oracle              the validity model at the true path (evaluation only: an upper reference)

Selection at the alarm. A forecast is issued at the first threshold crossing,
which depends on the data. BayesRamp uses the joint likelihood of every
statistic observed up to the alarm, so the stopping rule adds no factor that
depends on the parameters and no selection correction is needed. Statistics of
overlapping windows are correlated (1/2 for adjacent windows of a linear
statistic); the joint Gaussian likelihood carries that covariance.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Mapping, Protocol, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy import linalg, special

from .families import LinearRamp, RampPrior, SeverityPath, StepHold
from .surrogate import SurrogateSpec, Trajectories, simulate
from .validity import ValidityModel, coverage_probability_fast

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class History:
    """What is known at issuance."""
    alarm_index: int
    alarm_sweep: int
    statistic: FloatArray  # statistics of checkpoints 0..alarm_index, a private copy


def history_at(statistic: FloatArray, ends: NDArray[np.int64], alarm_index: int) -> History:
    """The admissible prefix of one run's statistic series."""
    if not 0 <= alarm_index < len(ends):
        raise ValueError("alarm index outside the checkpoint grid")
    prefix = np.array(statistic[:alarm_index + 1], dtype=np.float64, copy=True)
    prefix.setflags(write=False)
    return History(int(alarm_index), int(ends[alarm_index]), prefix)


class ValidityForecaster(Protocol):
    name: str

    def predict(self, histories: Sequence[History]) -> dict[int, FloatArray]:
        """Probability that the run is valid at each declared length, one value per history."""
        ...


# ---------------------------------------------------------------------------
# Parameter estimation from labelled training and healthy calibration runs
# ---------------------------------------------------------------------------
def estimate_sigma_eff(healthy: Trajectories) -> float:
    """Long-run noise scale from healthy runs: root mean square of sqrt(L) * standard error over runs and lengths."""
    spec = healthy.spec
    squares = [length * (healthy.interval(length)[1] / spec.critical)**2 for length in spec.lengths]
    return float(np.sqrt(np.mean(squares)))


def estimate_beta(train: Trajectories) -> float:
    """Bias slope: least squares through the origin of the signed error on the mean severity, pooled over lengths."""
    spec = train.spec
    x, e = [], []
    for length in spec.lengths:
        x.append(train.path.cumulative(np.array([length], dtype=np.float64))[..., 0] / length)
        e.append(train.interval(length)[0])
    x, e = np.concatenate(x), np.concatenate(e)
    return float(np.dot(x, e) / np.dot(x, x))


def estimate_kappa(train: Trajectories) -> float:
    """Sensitivity of the monitoring statistic: least squares through the origin on the window-mean severity."""
    x = train.spec.window_mean_severity(train.path)
    return float(np.sum(x * train.statistic) / np.sum(x * x))


# ---------------------------------------------------------------------------
# Baselines
# ---------------------------------------------------------------------------
@dataclass
class BaseRate:
    name: str = "base-rate"
    rate: dict[int, float] = field(default_factory=dict)

    def fit(self, histories: Sequence[History], outcomes: Mapping[int, NDArray[np.bool_]]) -> "BaseRate":
        self.rate = {length: float(np.mean(y)) for length, y in outcomes.items()}
        return self

    def predict(self, histories: Sequence[History]) -> dict[int, FloatArray]:
        return {length: np.full(len(histories), rate) for length, rate in self.rate.items()}


def _fit_logistic(x: FloatArray, y: FloatArray, iterations: int = 100) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Standardise, then minimise summed log-loss + ||w||^2/2 with a free intercept (Newton)."""
    mu, sd = x.mean(axis=0), x.std(axis=0)
    scale = np.where(sd == 0, 1.0, sd)
    design = np.column_stack([(x - mu) / scale, np.ones(len(y))])
    penalty = np.diag(np.r_[np.ones(x.shape[1]), 0.0])
    theta = np.zeros(design.shape[1])
    for _ in range(iterations):
        p = special.expit(design @ theta)
        gradient = design.T @ (p - y) + penalty @ theta
        if np.max(np.abs(gradient)) < 1e-9:
            break
        hessian = design.T @ (design * (p * (1 - p))[:, None]) + penalty + 1e-10 * np.eye(len(theta))
        theta = theta - np.linalg.solve(hessian, gradient)
    return theta, mu, scale


@dataclass
class LogisticForecaster:
    """Logistic regression on fixed features of the history; one model per declared length.

    features="alarm" uses the alarm sweep only; "lite" adds the statistic at the alarm and its
    least-squares slope over the last min(4, j+1) checkpoints (0 for a single point), as in the
    Stage D pilot's F-B-lite. A one-class training set falls back to the base rate.
    """
    features: str = "lite"
    stride: int = 256
    name: str = ""
    models: dict[int, tuple] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.features not in ("alarm", "lite"):
            raise ValueError("features must be 'alarm' or 'lite'")
        self.name = self.name or ("alarm-time" if self.features == "alarm" else "lite-logistic")

    def _design(self, histories: Sequence[History]) -> FloatArray:
        rows = []
        for h in histories:
            if self.features == "alarm":
                rows.append([h.alarm_sweep])
                continue
            tail = h.statistic[-4:]
            if len(tail) == 1:
                slope = 0.0
            else:
                t = (np.arange(len(tail)) - (len(tail) - 1) / 2) * self.stride
                slope = float(np.dot(t, tail - tail.mean()) / np.dot(t, t))
            rows.append([h.alarm_sweep, h.statistic[-1], slope])
        return np.array(rows, dtype=np.float64).reshape(len(histories), -1)

    def fit(self, histories: Sequence[History], outcomes: Mapping[int, NDArray[np.bool_]]) -> "LogisticForecaster":
        x = self._design(histories)
        for length, y in outcomes.items():
            y = np.asarray(y, dtype=np.float64)
            self.models[length] = ("rate", float(y.mean())) if y.min() == y.max() else ("fit", *_fit_logistic(x, y))
        return self

    def predict(self, histories: Sequence[History]) -> dict[int, FloatArray]:
        x = self._design(histories)
        out = {}
        for length, model in self.models.items():
            if model[0] == "rate":
                out[length] = np.full(len(histories), model[1])
            else:
                _, theta, mu, scale = model
                out[length] = special.expit(((x - mu) / scale) @ theta[:-1] + theta[-1])
        return out


# ---------------------------------------------------------------------------
# Model-based forecasters
# ---------------------------------------------------------------------------
@dataclass
class Persistence:
    """The severity implied by the alarm-window statistic holds from the middle of that window on;
    before it the input is taken as healthy. No growth and no recovery are assumed."""
    spec: SurrogateSpec
    model: ValidityModel
    kappa: float
    name: str = "persistence"

    def predict(self, histories: Sequence[History]) -> dict[int, FloatArray]:
        level = np.clip(np.array([h.statistic[-1] for h in histories]) / self.kappa, 0.0, self.spec.s_max)
        onset = np.array([h.alarm_sweep - self.spec.window // 2 for h in histories], dtype=np.float64)
        path = StepHold(onset, level)
        return {length: self.model.probability(path, length) for length in self.spec.lengths}


@dataclass
class BayesRamp:
    """Posterior predictive forecaster under the declared ramp family.

    fit() fixes a grid of ramp hypotheses with the declared prior weights, the mean of
    the statistic at every checkpoint under each hypothesis, and the validity
    probability of each hypothesis at each declared length. predict() weighs the
    hypotheses by the joint Gaussian likelihood of the observed statistic history.
    """
    spec: SurrogateSpec
    prior: RampPrior
    onsets: int = 49
    rates: int = 48
    name: str = "bayes-ramp"
    drift_widening: bool = True

    def fit(self, beta: float, sigma_eff: float, kappa: float) -> "BayesRamp":
        self.grid, self.log_prior = self.prior.grid(self.onsets, self.rates)
        self.log_prior = np.log(self.log_prior)
        self.mean = kappa * self.spec.window_mean_severity(self.grid)       # (hypotheses, checkpoints)
        self.model = ValidityModel(beta, sigma_eff, self.spec.batches, self.spec.level)
        self.p_valid = {}
        for length in self.spec.lengths:
            delta, lam = self.model.noncentralities(self.grid, length)
            self.p_valid[length] = coverage_probability_fast(delta, lam if self.drift_widening else 0.0,
                                                             self.spec.batches, self.spec.level)
        self.failure_pmf = None
        self._cholesky: dict[int, FloatArray] = {}
        return self

    def _banded_cholesky(self, n: int) -> FloatArray:
        """Cholesky factor (upper banded form) of the n x n covariance with 1 on the diagonal and 1/2 beside it."""
        if n not in self._cholesky:
            band = np.vstack([np.r_[0.0, np.full(n - 1, 0.5)], np.ones(n)])
            self._cholesky[n] = linalg.cholesky_banded(band)
        return self._cholesky[n]

    def posterior(self, history: History) -> FloatArray:
        """Posterior weights of the grid hypotheses given the statistic history."""
        n = history.alarm_index + 1
        residual = history.statistic[None, :] - self.mean[:, :n]
        solved = linalg.cho_solve_banded((self._banded_cholesky(n), False), residual.T)
        log_weight = self.log_prior - 0.5 * np.einsum("ij,ji->i", residual, solved)
        weight = np.exp(log_weight - log_weight.max())
        return weight / weight.sum()

    def predict(self, histories: Sequence[History]) -> dict[int, FloatArray]:
        weights = np.stack([self.posterior(h) for h in histories]) if histories else np.zeros((0, len(self.log_prior)))
        return {length: weights @ p for length, p in self.p_valid.items()}

    def fit_failure_times(self, rng: np.random.Generator, simulations: int = 100, chunk: int = 20000) -> "BayesRamp":
        """Failure-time distribution of every hypothesis by forward simulation under the declared model.

        The parameters of the simulated runs are the fitted ones (self.model, kappa is not
        needed: the failure time depends on the observable only). Column k < checkpoints
        is 'failed at checkpoint k'; the last column is 'not failed inside the horizon'.
        """
        spec = replace(self.spec, beta=self.model.beta, sigma_eff=self.model.sigma_eff)
        hypotheses, bins = len(self.log_prior), len(spec.ends) + 1
        counts = np.zeros((hypotheses, bins))
        index = np.repeat(np.arange(hypotheses), simulations)
        t0, t_cap = np.asarray(self.grid.t0), np.asarray(self.grid.t_cap)
        for start in range(0, len(index), chunk):
            rows = index[start:start + chunk]
            runs = simulate(spec, LinearRamp(t0[rows], t_cap[rows], self.grid.s_max), len(rows), rng)
            time, observed = runs.failure
            column = np.where(observed, np.searchsorted(spec.ends, time), bins - 1)
            np.add.at(counts, (rows, column), 1.0)
        self.failure_pmf = counts / simulations
        return self

    def predict_failure(self, histories: Sequence[History]) -> FloatArray:
        """(histories, checkpoints + 1): probability of failure at each checkpoint, last column = censored."""
        if self.failure_pmf is None:
            raise RuntimeError("fit_failure_times must run before predict_failure")
        return np.stack([self.posterior(h) @ self.failure_pmf for h in histories])


@dataclass
class Oracle:
    """Validity model at the true severity path. Uses evaluation truth: a reference, not a forecaster."""
    model: ValidityModel
    lengths: tuple[int, ...]
    name: str = "oracle"

    def predict_paths(self, path: SeverityPath) -> dict[int, FloatArray]:
        return {length: self.model.probability(path, length) for length in self.lengths}


@dataclass
class EmpiricalFailure:
    """Failure-time baseline: the training distribution of the failure time among runs whose alarm
    fell in the same quantile bin of the alarm sweep. Uses no statistic beyond the alarm time."""
    bins: int = 5
    name: str = "empirical-by-alarm-time"

    def fit(self, histories: Sequence[History], time: NDArray[np.int64], observed: NDArray[np.bool_],
            ends: NDArray[np.int64]) -> "EmpiricalFailure":
        alarm = np.array([h.alarm_sweep for h in histories])
        self.edges = np.quantile(alarm, np.linspace(0, 1, self.bins + 1)[1:-1])
        group = np.searchsorted(self.edges, alarm, side="right")
        column = np.where(observed, np.searchsorted(ends, time), len(ends))
        self.pmf = np.zeros((self.bins, len(ends) + 1))
        np.add.at(self.pmf, (group, column), 1.0)
        totals = self.pmf.sum(axis=1, keepdims=True)
        overall = self.pmf.sum(axis=0) / self.pmf.sum()
        self.pmf = np.where(totals > 0, self.pmf / np.maximum(totals, 1), overall)
        return self

    def predict_failure(self, histories: Sequence[History]) -> FloatArray:
        alarm = np.array([h.alarm_sweep for h in histories])
        return self.pmf[np.searchsorted(self.edges, alarm, side="right")]
