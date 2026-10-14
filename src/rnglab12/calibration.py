"""Finite-horizon healthy calibration on independent whole trajectories.

The rank bound is marginal over exchangeable calibration and future healthy
trajectories. It is not a conditional guarantee for a realized calibration set,
nor a guarantee for shifted source distributions, unbounded monitoring, or
arbitrarily changed observables, samplers, time grids, and configurations.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
import math
from typing import Sequence

import numpy as np
from scipy.stats import beta

from .seeds import SPLITS, assert_disjoint_lineages


def _probability(value: float, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 < value < 1:
        raise ValueError(f'{name} must lie strictly between zero and one')
    return float(value)


def _family(value: int) -> int:
    if type(value) is not int or value < 1:
        raise ValueError('family_size must be a positive integer')
    return value


def _rank(n: int, alpha: float, family_size: int) -> int:
    requested = (n + 1) * (1 - Fraction(str(alpha)) / family_size)
    return -(-requested.numerator // requested.denominator)


@dataclass(frozen=True)
class Trajectory:
    errors: Sequence[float]
    times: Sequence[int | float]
    application: str
    observable: str
    config_id: str
    lineage_id: str
    split: str
    reference: float | None = None
    seed: int | None = None

    def __post_init__(self):
        for name in ('application', 'observable', 'config_id', 'lineage_id'):
            if not isinstance(getattr(self, name), str) or not getattr(self, name):
                raise ValueError(f'A nonempty {name} is required')
        if self.split not in SPLITS:
            raise ValueError('Unknown trajectory partition')
        error = np.asarray(self.errors, dtype=float)
        times = np.asarray(self.times, dtype=float)
        if error.ndim != 1 or times.ndim != 1 or len(error) != len(times) or len(error) < 2:
            raise ValueError('At least two errors at matching frozen observation times are required')
        if not np.isfinite(error).all() or (error < 0).any():
            raise ValueError('Errors must be finite nonnegative absolute errors')
        if not np.isfinite(times).all() or (np.diff(times) <= 0).any():
            raise ValueError('Finite strictly increasing observation times are required')
        if self.reference is not None and (isinstance(self.reference, bool) or not math.isfinite(self.reference)):
            raise ValueError('Reference must be finite or None')
        if self.seed is not None and (isinstance(self.seed, bool) or not isinstance(self.seed, int) or self.seed < 0):
            raise ValueError('Seed must be a nonnegative integer or None')
        # Copy mutable inputs so a fitted identity cannot silently change later.
        object.__setattr__(self, 'errors', tuple(float(x) for x in error))
        object.__setattr__(self, 'times', tuple(float(x) for x in times))

    def lineage(self) -> dict:
        return {'lineage_id': self.lineage_id, 'split': self.split, 'seed': self.seed}

    def score(self) -> float:
        if self.reference is None:
            raise ValueError('ABSTAIN_UNKNOWN_REFERENCE: failure calibration requires a validated target')
        return max(min(a, b) for a, b in zip(self.errors[:-1], self.errors[1:]))


def _same_identity(first: Trajectory, other: Trajectory) -> None:
    fields = ('application', 'observable', 'config_id', 'times', 'reference')
    if any(getattr(first, name) != getattr(other, name) for name in fields):
        raise ValueError('Application, observable, target, frozen time grid and configuration must match')


@dataclass(frozen=True)
class FailurePolicy:
    application: str
    observable: str
    config_id: str
    times: tuple[float, ...]
    reference: float
    threshold: float
    alpha: float
    family_size: int
    n_calibration: int
    rank: int
    calibration_lineages: tuple[str, ...]
    calibration_seeds: tuple[int | None, ...]

    @property
    def attainable_bound(self) -> float:
        return (self.n_calibration + 1 - self.rank) / (self.n_calibration + 1)

    def _check(self, trajectory: Trajectory) -> None:
        for name in ('application', 'observable', 'config_id', 'times'):
            if getattr(trajectory, name) != getattr(self, name):
                raise ValueError(f'Frozen calibration identity mismatch: {name}')
        if trajectory.reference is not None and trajectory.reference != self.reference:
            raise ValueError('Frozen calibration target mismatch')

    def event(self, trajectory: Trajectory) -> dict:
        self._check(trajectory)
        if trajectory.reference is None:
            return {'status': 'ABSTAIN_UNKNOWN_REFERENCE', 'event': None,
                    'first_event_time': None, 'score': None}
        # Adjacent observations on the frozen full grid are required. An isolated
        # crossing and two nonconsecutive crossings are not a sustained event.
        above = [value > self.threshold for value in trajectory.errors]
        first = next((j for j in range(1, len(above)) if above[j - 1] and above[j]), None)
        return {'status': 'EVALUATED', 'event': first is not None,
                'first_event_time': None if first is None else trajectory.times[first],
                'first_event_index': first, 'score': trajectory.score(),
                'threshold': self.threshold, 'strict_comparison': True}

    def to_dict(self) -> dict:
        return {'schema': 'rnglab12.finite-horizon-policy.v1', **asdict(self),
                'attainable_marginal_bound': self.attainable_bound,
                'guarantee_scope': 'exchangeable whole-trajectory marginal rank bound; not conditional on realized calibration',
                'sustained_crossings': 2, 'comparison': 'strict_gt'}

    @classmethod
    def from_dict(cls, data: dict) -> 'FailurePolicy':
        if data.get('schema') != 'rnglab12.finite-horizon-policy.v1' or data.get('comparison') != 'strict_gt' or data.get('sustained_crossings') != 2:
            raise ValueError('Unknown calibration policy schema or event rule')
        fields = cls.__dataclass_fields__
        values = {name: data[name] for name in fields}
        for name in ('times', 'calibration_lineages', 'calibration_seeds'):
            values[name] = tuple(values[name])
        policy = cls(**values)
        _probability(policy.alpha, 'alpha')
        _family(policy.family_size)
        if type(policy.n_calibration) is not int or policy.n_calibration < 1 or policy.rank != _rank(policy.n_calibration, policy.alpha, policy.family_size) or not 1 <= policy.rank <= policy.n_calibration:
            raise ValueError('Invalid calibration count or rank')
        if len(policy.calibration_lineages) != policy.n_calibration or len(policy.calibration_seeds) != policy.n_calibration:
            raise ValueError('Calibration lineage count mismatch')
        if not math.isfinite(policy.threshold) or policy.threshold < 0 or not math.isfinite(policy.reference):
            raise ValueError('Invalid calibration threshold or reference')
        Trajectory([0] * len(policy.times), policy.times, policy.application,
                   policy.observable, policy.config_id, 'policy-validation', 'development', policy.reference)
        assert_disjoint_lineages([{'lineage_id': line, 'seed': seed, 'split': 'calibration'}
                                  for line, seed in zip(policy.calibration_lineages, policy.calibration_seeds)])
        return policy


def calibrate(trajectories: Sequence[Trajectory], alpha: float = .01, family_size: int = 1) -> FailurePolicy:
    """Freeze a Bonferroni-adjusted whole-trajectory sustained-crossing policy.

    A trajectory is the independent unit, including every planned monitoring
    look. Healthy calibration membership is an experimental-design input and
    cannot be established just by inspecting its observed trajectory.
    """
    alpha = _probability(alpha, 'alpha')
    family_size = _family(family_size)
    items = list(trajectories)
    if not items:
        raise ValueError('Independent healthy calibration trajectories are required')
    if any(item.split != 'calibration' for item in items):
        raise ValueError('Only the calibration partition may fit thresholds')
    assert_disjoint_lineages([item.lineage() for item in items])
    for item in items:
        _same_identity(items[0], item)
    scores = sorted(item.score() for item in items)
    rank = _rank(len(items), alpha, family_size)
    if rank > len(items):
        minimum = -(-(Fraction(family_size) / Fraction(str(alpha))).numerator //
                     (Fraction(family_size) / Fraction(str(alpha))).denominator) - 1
        raise ValueError(f'Insufficient independent trajectories: need at least {minimum} for requested finite rank bound')
    return FailurePolicy(items[0].application, items[0].observable, items[0].config_id,
                         items[0].times, float(items[0].reference), float(scores[rank - 1]),
                         alpha, family_size, len(items), rank,
                         tuple(item.lineage_id for item in items), tuple(item.seed for item in items))


def clopper_pearson_upper(events: int, count: int, confidence: float = .95, family_size: int = 1) -> float:
    """Exact one-sided binomial upper bound, with Bonferroni confidence control."""
    confidence = _probability(confidence, 'confidence')
    family_size = _family(family_size)
    if type(count) is not int or type(events) is not int or count < 1 or not 0 <= events <= count:
        raise ValueError('An integer event count between zero and the positive group count is required')
    tail = (1 - confidence) / family_size
    if events == count:
        return 1.0
    if events == 0:
        return -math.expm1(math.log(tail) / count)
    return float(beta.ppf(1 - tail, events + 1, count - events))


def validate_holdout(policy: FailurePolicy, trajectories: Sequence[Trajectory],
                     confidence: float = .95, family_size: int = 1,
                     target_rate: float | None = None) -> dict:
    """Evaluate a frozen policy once on new independent healthy seed groups.

    Repeatedly selecting a policy after seeing this validation invalidates its
    holdout interpretation; no function can enforce honest future data use.
    """
    confidence = _probability(confidence, 'confidence')
    family_size = _family(family_size)
    items = list(trajectories)
    if not items:
        raise ValueError('Independent healthy holdout trajectories are required')
    if any(item.split not in {'confirmation', 'final_holdout'} for item in items):
        raise ValueError('Use a distinct confirmation or final_holdout partition for validation')
    calibration = [{'lineage_id': line, 'seed': seed, 'split': 'calibration'}
                   for line, seed in zip(policy.calibration_lineages, policy.calibration_seeds)]
    assert_disjoint_lineages(calibration, [item.lineage() for item in items])
    results = [policy.event(item) for item in items]
    if any(item['event'] is None for item in results):
        return {'status': 'ABSTAIN_UNKNOWN_REFERENCE', 'groups': len(items),
                'events': None, 'upper_bound': None, 'passes': False}
    count = len(items)
    events = sum(item['event'] for item in results)
    upper = clopper_pearson_upper(events, count, confidence, family_size)
    target = policy.alpha / policy.family_size if target_rate is None else _probability(target_rate, 'target_rate')
    return {'status': 'EVALUATED', 'groups': count, 'events': events,
            'empirical_rate': events / count, 'upper_bound': upper,
            'simultaneous_confidence': confidence, 'confidence_family_size': family_size,
            'per_policy_target': target, 'calibration_per_policy_target': policy.alpha / policy.family_size,
            'validation_target_explicitly_overridden': target_rate is not None, 'passes': upper <= target,
            'criterion': 'one-sided exact binomial upper bound <= frozen per-policy false-alarm target',
            'assumptions': 'independent identically distributed healthy holdout groups; fixed policy; no holdout tuning',
            'first_event_times': [item['first_event_time'] for item in results]}
