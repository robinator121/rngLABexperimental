import json
import math

import numpy as np
import pytest

from rnglab12.calibration import (Trajectory, FailurePolicy, calibrate,
                                 clopper_pearson_upper, validate_holdout)
from rnglab12.seeds import derive_seed, seed_lineage, assert_disjoint_lineages
from rnglab12.uncertainty import summarize_series


def trajectory(group, errors=(1., 1.), split='calibration', **kwargs):
    lineage = seed_lineage(split, 'test-app', group)
    arguments = dict(errors=errors, times=range(len(errors)), application='test-app',
                     observable='mean', config_id='frozen-config-hash',
                     lineage_id=lineage['lineage_id'], split=split,
                     reference=0., seed=lineage['seed'])
    arguments.update(kwargs)
    return Trajectory(**arguments)


def test_seed_domains_and_relabelled_seed_leakage():
    values = [derive_seed(split, 'x', 1) for split in
              ('development', 'calibration', 'confirmation', 'final_holdout', 'pilot')]
    assert len(set(values)) == 5
    assert derive_seed('calibration', 'x', 1) == derive_seed('calibration', 'x', 1)
    assert derive_seed('calibration', 'x', 1, 'source') != derive_seed('calibration', 'x', 1, 'acceptance')
    with pytest.raises(ValueError):
        derive_seed('untracked', 'x', 1)
    original = seed_lineage('calibration', 'x', 1)
    relabelled = {**seed_lineage('final_holdout', 'x', 1), 'seed': original['seed']}
    with pytest.raises(ValueError, match='numeric seed'):
        assert_disjoint_lineages([original], [relabelled])


def test_exact_attainable_rank_and_minimum_sample_size():
    with pytest.raises(ValueError, match='at least 99'):
        calibrate([trajectory(j) for j in range(98)])
    policy = calibrate([trajectory(j, (j, j)) for j in range(99)])
    assert policy.rank == 99
    assert policy.threshold == 98
    assert policy.attainable_bound == .01
    assert 'not conditional' in policy.to_dict()['guarantee_scope']


def test_family_adjustment_and_serializable_frozen_identity():
    traces = [trajectory(j, (j, j)) for j in range(199)]
    single = calibrate(traces, alpha=.05)
    family = calibrate(traces, alpha=.05, family_size=2)
    assert family.rank > single.rank
    assert family.attainable_bound <= .025
    restored = FailurePolicy.from_dict(json.loads(json.dumps(family.to_dict())))
    assert restored == family
    with pytest.raises(ValueError, match='rank'):
        FailurePolicy.from_dict({**family.to_dict(), 'rank': 1})
    for field, value in [('application', 'different'), ('observable', 'other'),
                         ('config_id', 'other-config'), ('times', (0, 2)), ('reference', 1.)]:
        with pytest.raises(ValueError, match='mismatch'):
            restored.event(trajectory(5, split='confirmation', **{field: value}))


def test_calibration_uses_whole_monitoring_horizon_not_pointwise_quantiles():
    # Every healthy trajectory has one sustained spike, at a different pair of
    # looks. Each pointwise 85th percentile is zero, yet that policy alarms on
    # every trajectory. The whole-history score includes the monitoring search.
    traces = []
    for j in range(99):
        errors = np.zeros(20)
        errors[2 * (j % 10):2 * (j % 10) + 2] = 1
        traces.append(trajectory(j, errors))
    pointwise = np.quantile(np.asarray([t.errors for t in traces]), .85, axis=0)
    assert np.all(pointwise == 0)
    policy = calibrate(traces, alpha=.15)
    assert policy.threshold == 1
    assert sum(policy.event(t)['event'] for t in traces) == 0  # conservative strict ties
    assert all(t.score() > 0 for t in traces)  # pointwise zero threshold alarms


def test_strict_sustained_events_and_unknown_reference_abstention():
    policy = calibrate([trajectory(j, (1, 1, 1, 1)) for j in range(99)])
    assert policy.event(trajectory('tie', (1, 1, 1, 1), split='confirmation'))['event'] is False
    assert policy.event(trajectory('separate', (2, 0, 2, 0), split='confirmation'))['event'] is False
    event = policy.event(trajectory('consecutive', (0, 2, 2, 0), split='confirmation'))
    assert event['event'] is True and event['first_event_time'] == 2
    result = policy.event(trajectory('unknown', (2, 2, 2, 2), split='confirmation', reference=None))
    assert result['event'] is None and result['status'] == 'ABSTAIN_UNKNOWN_REFERENCE'
    with pytest.raises(ValueError, match='UNKNOWN_REFERENCE'):
        calibrate([trajectory(j, reference=None) for j in range(99)])


def test_invalid_trace_or_mixed_calibration_identity_rejected():
    for errors in [(math.nan, 0), (math.inf, 1), (-1, 0)]:
        with pytest.raises(ValueError):
            trajectory(1, errors)
    with pytest.raises(ValueError, match='times'):
        trajectory(1, times=(1, 1))
    with pytest.raises(ValueError, match='calibration partition'):
        calibrate([trajectory(j, split='development') for j in range(99)])
    with pytest.raises(ValueError, match='must match'):
        calibrate([trajectory(j, observable='other' if j == 98 else 'mean') for j in range(99)])
    same = trajectory('reused')
    with pytest.raises(ValueError, match='Reused lineage'):
        calibrate([same] * 99)


def test_holdout_independent_groups_family_confidence_and_explicit_target():
    policy = calibrate([trajectory(j) for j in range(99)])
    holdout = [trajectory(j, (0, 0), split='final_holdout') for j in range(256)]
    result = validate_holdout(policy, holdout, family_size=6, target_rate=.05)
    assert result['events'] == 0 and result['passes'] is True
    assert result['upper_bound'] == pytest.approx(1 - (.05 / 6) ** (1 / 256))
    assert result['validation_target_explicitly_overridden'] is True
    assert validate_holdout(policy, holdout, family_size=6)['passes'] is False
    assert clopper_pearson_upper(0, 256, family_size=6) > clopper_pearson_upper(0, 256)
    assert clopper_pearson_upper(256, 256) == 1
    with pytest.raises(ValueError, match='distinct'):
        validate_holdout(policy, [trajectory(1)])
    with pytest.raises(ValueError, match='lineage'):
        validate_holdout(policy, [trajectory(1, split='final_holdout', lineage_id=policy.calibration_lineages[0])])
    with pytest.raises(ValueError, match='numeric seed'):
        validate_holdout(policy, [trajectory(1, split='final_holdout', seed=policy.calibration_seeds[0])])
    with pytest.raises(ValueError, match='lineage'):
        validate_holdout(policy, [holdout[0], holdout[0]])


def test_geyer_correlation_and_batch_comparison():
    rng = np.random.default_rng(392015)
    innovations = rng.normal(size=16384)
    trace = np.empty_like(innovations)
    trace[0] = innovations[0] / math.sqrt(1 - .8**2)
    for j in range(1, len(trace)):
        trace[j] = .8 * trace[j - 1] + innovations[j]
    result = summarize_series(trace, reference=0.)
    assert result['geyer']['status'] == 'ESTIMATED'
    assert result['geyer']['mcse'] > 2 * result['naive_mcse']
    # Stationary AR(1) has tau=(1+phi)/(1-phi)=9. A broad fixed range tests
    # correlation handling without pretending this realization is a proof.
    assert 6 < result['geyer']['tau'] < 13
    assert result['geyer']['ess'] < len(trace) / 6
    assert set(result['batch_means']) == {'32', '128'}
    assert all(item['status'] == 'ESTIMATED' for item in result['batch_means'].values())
    assert result['stationarity_verified'] is False
    json.dumps(result, allow_nan=False)


def test_uncertainty_safeguards_unknown_target_and_batch_prefix():
    for bad in [np.arange(31), np.r_[np.arange(32), np.inf], np.ones((32, 2))]:
        with pytest.raises(ValueError):
            summarize_series(bad)
    constant = summarize_series(np.ones(1024), reference=1.)
    assert constant['status'] == 'ABSTAIN_CONSTANT_TRACE'
    assert constant['geyer']['mcse'] is None and constant['geyer']['covered'] is None
    result = summarize_series(np.random.default_rng(4).normal(size=1027))
    assert result['reference_status'] == 'ABSTAIN_UNKNOWN_REFERENCE'
    assert result['geyer']['covered'] is None
    assert result['batch_means']['32']['prefix_n'] == 1024
    assert result['batch_means']['32']['discarded_tail'] == 3
    assert result['batch_means']['32']['covered'] is None
    assert summarize_series(np.arange(32))['batch_means']['32']['status'] == 'INSUFFICIENT_BATCHES'
    with pytest.raises(ValueError):
        summarize_series(np.arange(32), batch_sizes=(32, 32))
