"""Independent window/event logic and small-kernel controls for the locked fault study."""
import math

import numpy as np
import pytest

import fault_study
import study


def test_rolling_windows_match_direct_slices_and_recorded_cross_products():
    r, n, window, stride = 2, 12, 4, 2
    raw = np.zeros((r, n, 3))
    raw[:, :, 1] = np.arange(r*n).reshape(r, n)/30
    raw[:, :, 2] = np.arange(r*n).reshape(r, n)/100
    source = np.arange(r*n*3, dtype=float).reshape(r, n, 3)/100
    ends, scores, failure = fault_study.windows(raw, source, window=window, stride=stride)
    assert np.array_equal(ends, [4, 6, 8, 10, 12])
    obs = study.observables(raw)
    direct = lambda a: np.column_stack([a[:, end-window:end].mean(1) for end in ends])
    up, ua = direct(source[:, :, 0]), direct(source[:, :, 1])
    assert scores['source_marginal'] == pytest.approx(np.maximum(abs(up-.5), abs(ua-.5)), abs=2e-15)
    # Cross-channel score averages the per-proposal cross product already
    # recorded by the kernel; it is not the product of window-averaged margins.
    assert scores['source_cross_channel'] == pytest.approx(abs(direct(source[:, :, 2])), abs=2e-15)
    assert scores['application_shape'] == pytest.approx(
        abs(direct(obs['shape'])-study.reference()['means']['shape']), abs=2e-15)
    errors = abs(direct(obs['CV'])-study.reference()['finite_energy'])
    expected = []
    for row in errors:
        events = [end for j, end in enumerate(ends[1:]) if row[j]>.05 and row[j+1]>.05]
        expected.append(events[0] if events else -1)
    assert failure.tolist() == expected


def _direct_cv_windows(monkeypatch, values):
    values = np.asarray(values, dtype=float)
    raw = np.zeros((len(values), values.shape[1], 3))
    raw[:, :, 0] = values
    # Isolate window/event bookkeeping from the independently tested energy
    # formula so exact equality to the strict .05 tolerance is testable.
    monkeypatch.setattr(study, 'observables', lambda a: {'CV':a[:, :, 0], 'shape':a[:, :, 2]})
    monkeypatch.setattr(study, 'reference', lambda: {'finite_energy':0., 'means':{'shape':0.}})
    return fault_study.windows(raw, np.full_like(raw, .5), window=1, stride=1)


def test_failure_requires_two_consecutive_windows_and_timestamps_second(monkeypatch):
    _, _, failure = _direct_cv_windows(monkeypatch, [
        [.06, .04, .06, .06, .02],  # first confirmed event at endpoint4
        [.06, .04, .06, .04, .06],  # isolated excursions never persist
        [.06, .06, .01, .07, .07],  # retain the first event, endpoint2
    ])
    assert failure.tolist() == [4, -1, 2]


def test_failure_tolerance_is_strict_and_one_window_cannot_confirm(monkeypatch):
    _, _, failure = _direct_cv_windows(monkeypatch, [[.05, .05], [.050001, .050001]])
    assert failure.tolist() == [-1, 2]
    _, _, short = _direct_cv_windows(monkeypatch, [[.9]])
    assert short.tolist() == [-1]


def test_first_hit_strict_threshold_censoring_and_onset_endpoint_rule():
    ends = np.array([1792, 2048, 2304, 2560])
    scores = np.array([[1., 1., 1., 1.], [2., 2., 1., 3.], [.5, .5, 2., .5], [.5, .5, .5, .5]])
    assert fault_study.first_hit(scores, ends, 1.).tolist() == [-1, 1792, 2304, -1]
    assert fault_study.first_hit(scores, ends, 1., start=2049).tolist() == [-1, 2560, 2304, -1]
    # Endpoint2304 intentionally includes preonset observations with window1024.
    assert 2304-1024+1 == 1281
    assert 2304-2049+1 == 256


def _python_fault_kernel(seed, initial, steps, beta, delta, fault, onset, severity):
    """Small direct action-difference kernel, separate from the compiled code."""
    rng = np.random.default_rng(seed)
    x = initial.copy()
    m = len(x)
    rows, recorded, accepted = [], [], 0
    def action(y):
        return m/(2*beta)*sum((y[(j+1)%m]-y[j])**2 for j in range(m)) + beta/(2*m)*sum(y*y)
    for sweep in range(steps):
        inputs = []
        for j in range(m):
            up, ua = rng.random(), rng.random()
            if sweep >= onset:
                if fault == 1:
                    ua = min(ua+severity, np.nextafter(1., 0.))
                elif fault == 2:
                    ua = up
                elif fault == 3:
                    up = min(up+severity, np.nextafter(1., 0.))
            candidate = x.copy()
            candidate[j] += delta*(2*up-1)
            ds = action(candidate)-action(x)
            if ds <= 0 or ua < math.exp(-ds):
                x = candidate
                accepted += 1
            inputs.append([up, ua, (up-.5)*(ua-.5)])
        c = x.mean()
        shape = np.mean((x-c)**2)
        td = m/(2*beta)-m/(2*beta**2)*sum((x[(j+1)%m]-x[j])**2 for j in range(m)) + np.mean(x*x)/2
        rows.append([td, c, shape])
        recorded.append(np.mean(inputs, axis=0))
    return np.asarray(rows), np.asarray(recorded), accepted/(steps*m)


@pytest.mark.parametrize('fault', [0, 1, 2, 3])
def test_small_fault_kernel_matches_independent_whole_action_implementation(fault):
    seed, steps, beta, delta, onset, severity = 13571, 12, 2., .4, 3, .08
    initial = np.array([.4, -.2, .1, .5])
    raw, source, acceptance = study.local_kernel(
        np.random.default_rng(seed), np.random.default_rng(82), np.random.default_rng(83),
        initial, steps, beta, delta, False, False, fault, onset, severity, True)
    expected = _python_fault_kernel(seed, initial, steps, beta, delta, fault, onset, severity)
    assert raw == pytest.approx(expected[0], rel=1e-12, abs=1e-12)
    assert source == pytest.approx(expected[1], abs=2e-15)
    assert acceptance == expected[2]


@pytest.mark.parametrize('fault', [1, 2, 3])
def test_source_transformations_use_raw_pair_and_clip_at_largest_float_below_one(fault):
    seed, m, steps, onset, severity = 33447, 4, 6, 2, .8
    # Delta zero isolates exact source transformations from transition effects.
    raw, recorded, acceptance = study.local_kernel(
        np.random.default_rng(seed), np.random.default_rng(82), np.random.default_rng(83),
        np.zeros(m), steps, 2., 0., False, False, fault, onset, severity, True)
    pairs = np.random.default_rng(seed).random((steps, m, 2))
    up, ua = pairs[:, :, 0].copy(), pairs[:, :, 1].copy()
    if fault == 1:
        ua[onset:] = np.minimum(ua[onset:]+severity, np.nextafter(1., 0.))
    elif fault == 2:
        ua[onset:] = up[onset:]
    else:
        up[onset:] = np.minimum(up[onset:]+severity, np.nextafter(1., 0.))
    expected = np.stack([up.mean(1), ua.mean(1), ((up-.5)*(ua-.5)).mean(1)], axis=1)
    assert recorded == pytest.approx(expected, abs=2e-15)
    assert acceptance == 1.
    assert np.array_equal(raw[:, 1:], np.zeros((steps, 2)))
    assert np.all((up>=0)&(up<1)&(ua>=0)&(ua<1))


@pytest.mark.parametrize('fault', [1, 2, 3])
def test_exact_2049_sweep_onset_preserves_every_prefault_raw_and_source_record(fault):
    def run(code):
        return study.local_kernel(
            np.random.default_rng(28901), np.random.default_rng(28902), np.random.default_rng(28903),
            np.array([.4, -.2, .1, .5]), 2052, 2., .4,
            True, False, code, 2048, .08, True)
    healthy_raw, healthy_source, _ = run(0)
    fault_raw, fault_source, _ = run(fault)
    assert np.array_equal(healthy_raw[:2048], fault_raw[:2048])
    assert np.array_equal(healthy_source[:2048], fault_source[:2048])
    assert not np.array_equal(healthy_source[2048], fault_source[2048])
    assert not np.array_equal(healthy_raw[2048:], fault_raw[2048:])
    if fault == 2:
        assert np.array_equal(fault_source[2048:, 0], fault_source[2048:, 1])
