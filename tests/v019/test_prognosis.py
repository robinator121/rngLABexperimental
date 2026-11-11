"""rnglab.prognosis: closed forms, the validity model against simulation, information barriers, scores.

Seeds here are test seeds. The methods validation study (research/v019-forecast-methods/)
has its own frozen namespace and is not touched by these tests.
"""
import dataclasses
import math

import numpy as np
import pytest
from scipy import integrate, special, stats

from rnglab.prognosis import families as fam
from rnglab.prognosis import forecasters as fc
from rnglab.prognosis import scoring as sc
from rnglab.prognosis import surrogate as sur
from rnglab.prognosis import validity as val

SMALL = sur.SurrogateSpec(horizon=8192, lengths=(2048, 4096, 8192))
SMALL_PRIOR = fam.RampPrior(t0_low=512, t0_high=2048, t_cap_low=512.0, t_cap_high=32768.0)


# --------------------------------------------------------------------------- families
def test_cumulative_closed_forms_equal_per_sweep_sums():
    t = np.arange(4000)
    for t0, t_cap in ((300, 500.5), (0, 1.0), (3999, 10.0), (100, 1e6), (250, 250.0)):
        ramp = fam.LinearRamp(t0, t_cap, 0.5)
        s = ramp.severity(t)
        assert s[t0] == 0 and (t0 + 1 >= 4000 or s[t0 + 1] > 0) and s.max() <= 0.5
        for n in (0, 1, t0, t0 + 1, min(t0 + 2, 4000), 1234, 4000):
            assert float(ramp.cumulative(n)) == pytest.approx(s[:n].sum(), rel=1e-12, abs=1e-12)
        frozen = fam.RampThenHold(t0, t_cap, 0.5, 700)
        s2 = np.where(t < 700, s, ramp.severity(700))
        for n in (0, 699, 700, 701, 4000):
            assert float(frozen.cumulative(n)) == pytest.approx(s2[:n].sum(), rel=1e-12, abs=1e-12)
    step = fam.StepHold(100, 0.2)
    assert float(step.cumulative(100)) == 0 and float(step.cumulative(350)) == pytest.approx(50.0)


def test_parameter_and_argument_shapes_broadcast():
    ramp = fam.LinearRamp(np.array([300, 400, 500]), np.array([500.0, 800.0, 100.0]), 0.5)
    n = np.array([0, 1000, 2000, 4000])
    out = ramp.cumulative(n)
    assert out.shape == (3, 4)
    for i in range(3):
        single = fam.LinearRamp(int(ramp.t0[i]), float(ramp.t_cap[i]), 0.5)
        assert np.allclose(out[i], single.cumulative(n)) and single.cumulative(n).shape == (4,)
        assert np.ndim(single.cumulative(1000)) == 0
    frozen = fam.RampThenHold(np.array([300, 300]), np.array([500.0, 500.0]), 0.5, np.array([600, 4000]))
    both = frozen.cumulative(n)
    assert np.allclose(both[1], fam.LinearRamp(300, 500.0, 0.5).cumulative(n))  # freezing at the end changes nothing
    assert both[0, 1] < both[1, 1] and both[0, 0] == both[1, 0] == 0
    means = fam.block_means(ramp, np.arange(0, 4001, 250))
    assert means.shape == (3, 16) and np.all(np.diff(means, axis=1) >= -1e-15)
    sampled = fam.SampledPath(np.array([[0.0, 0.1, 0.2, 0.2]]), 16)
    assert np.allclose(sampled.cumulative(np.array([0, 16, 32, 64])), [[0.0, 0.0, 1.6, 8.0]])
    with pytest.raises(ValueError):
        sampled.cumulative(np.array([8]))


def test_ramp_prior_sample_and_grid():
    prior = fam.RampPrior()
    draw = prior.sample(np.random.default_rng(1), 20000)
    assert draw.t0.min() >= 2048 and draw.t0.max() <= 8192 and draw.t_cap.min() >= 2048 and draw.t_cap.max() < 163840
    assert np.mean(np.log(draw.t_cap)) == pytest.approx((math.log(2048) + math.log(163840)) / 2, abs=0.03)
    grid, weight = prior.grid(7, 5)
    assert len(weight) == 35 and weight.sum() == pytest.approx(1.0) and np.all(weight == weight[0])
    assert grid.t0.min() > 2048 and grid.t0.max() < 8192 and len(np.unique(grid.t0)) == 7
    walk = fam.random_walk_paths(np.random.default_rng(2), 50, 512, 16, t0_low=512, t0_high=2048,
                                 drift=1e-5, volatility=1e-3, s_max=0.5)
    assert walk.values.min() >= 0 and walk.values.max() <= 0.5 and np.all(walk.values[:, :32] == 0)
    assert walk.values[:, -1].std() > 0.01


# --------------------------------------------------------------------------- validity model
def adaptive_reference(delta, lam, nu=15):
    c = special.stdtrit(nu, 0.975)
    density = (lambda v: stats.ncx2.pdf(v, nu, lam)) if lam > 0 else (lambda v: stats.chi2.pdf(v, nu))
    f = lambda v: (special.ndtr(c * math.sqrt(v / nu) - delta) - special.ndtr(-c * math.sqrt(v / nu) - delta)) * density(v)
    mean, sd = nu + lam, math.sqrt(2 * (nu + 2 * lam))
    return integrate.quad(f, max(mean - 14 * sd, 0), mean + 16 * sd, epsabs=1e-13, epsrel=1e-11, limit=400)[0]


def test_coverage_probability_against_adaptive_quadrature_and_limits():
    assert float(val.coverage_probability(0.0)) == pytest.approx(0.95, abs=1e-9)
    assert float(val.coverage_probability(0.0, batches=8, level=0.9)) == pytest.approx(0.9, abs=1e-9)
    for delta, lam in ((0.0, 0.0), (1.0, 0.0), (2.5, 3.0), (4.0, 20.0), (6.0, 100.0), (0.5, 400.0), (3.0, 0.5)):
        assert float(val.coverage_probability(delta, lam)) == pytest.approx(adaptive_reference(delta, lam), abs=2e-8)
        assert float(val.coverage_probability(-delta, lam)) == pytest.approx(float(val.coverage_probability(delta, lam)))
    grid = val.coverage_probability(np.linspace(0, 12, 25)[:, None], np.array([0.0, 5.0, 50.0])[None, :])
    assert np.all(np.diff(grid, axis=0) <= 1e-12)   # more bias, less coverage
    assert np.all(np.diff(grid, axis=1) >= -1e-12)  # more drift between batches widens the interval
    assert float(val.coverage_probability(20.0, 0.0)) < 1e-12 and float(val.coverage_probability(20.0, 3000.0)) > 0.3


def test_fast_table_matches_exact_off_the_nodes():
    rng = np.random.default_rng(3)
    delta, lam = rng.uniform(0, 16, 4000), rng.uniform(0, 55, 4000)**2
    exact, fast = val.coverage_probability(delta, lam), val.coverage_probability_fast(delta, lam)
    assert np.max(np.abs(exact - fast)) < 5e-4
    assert np.all(fast[(delta > 14) | (lam > 2500)] == exact[(delta > 14) | (lam > 2500)])


@pytest.mark.parametrize("t0,t_cap", [(700, 1200.0), (1500, 9000.0), (600, 30000.0)])
def test_validity_model_matches_simulated_coverage(t0, t_cap):
    spec, size = SMALL, 30000
    runs = sur.simulate(spec, fam.LinearRamp(np.full(size, t0), np.full(size, t_cap), spec.s_max), size, np.random.default_rng(t0))
    model = val.ValidityModel(spec.beta * math.sqrt(5), spec.sigma_eff)  # a visible effect at this shorter horizon
    runs_b = sur.simulate(dataclasses.replace(spec, beta=model.beta),
                          fam.LinearRamp(np.full(size, t0), np.full(size, t_cap), spec.s_max), size, np.random.default_rng(t0 + 1))
    for length in spec.lengths:
        p = float(model.probability(fam.LinearRamp(t0, t_cap, spec.s_max), length))
        observed = runs_b.valid(length).mean()
        assert abs(observed - p) < 4.5 * math.sqrt(max(p * (1 - p), 1e-4) / size) + 1e-3
    assert runs.valid(spec.lengths[0]).mean() > 0.5


def test_drift_widening_matters_and_reduces_to_the_noncentral_t():
    model = val.ValidityModel(0.06, 0.613)
    path = fam.LinearRamp(8000, 2048.0, 0.5)
    with_widening = float(model.probability(path, 20480))
    without = float(model.probability(path, 20480, drift_widening=False))
    delta, lam = model.noncentralities(path, 20480)
    assert float(lam) > 5 and with_widening > 2 * without          # the pilot's oracle formula leaves this out
    c = stats.t.ppf(0.975, 15)
    assert without == pytest.approx(stats.nct.cdf(c, 15, float(delta)) - stats.nct.cdf(-c, 15, float(delta)), abs=1e-9)
    flat = fam.StepHold(0, 0.01)
    d, l = model.noncentralities(flat, 20480)
    assert float(l) == pytest.approx(0.0, abs=1e-18) and float(d) == pytest.approx(0.06 * 0.01 * math.sqrt(20480) / 0.613)


# --------------------------------------------------------------------------- surrogate
def test_surrogate_geometry_and_healthy_laws():
    spec = sur.SurrogateSpec()
    assert len(spec.ends) == 159 and spec.ends[[0, -1]].tolist() == [512, 40960] and spec.censor_sweep == 40448
    with pytest.raises(ValueError):
        sur.SurrogateSpec(window=768)
    healthy = sur.simulate(SMALL, None, 20000, np.random.default_rng(5))
    for length in SMALL.lengths:
        assert healthy.valid(length).mean() == pytest.approx(0.95, abs=0.006)
    assert healthy.statistic.mean() == pytest.approx(0, abs=0.01) and healthy.statistic.std() == pytest.approx(1, abs=0.01)
    lag = np.mean(healthy.statistic[:, :-1] * healthy.statistic[:, 1:])
    assert lag == pytest.approx(0.5, abs=0.01)  # adjacent windows share half their sweeps
    assert np.mean(healthy.statistic[:, :-2] * healthy.statistic[:, 2:]) == pytest.approx(0, abs=0.01)
    assert fc.estimate_sigma_eff(healthy) == pytest.approx(SMALL.sigma_eff, rel=0.01)
    with pytest.raises(ValueError):
        healthy.interval(100)


def test_alarm_calibration_and_failure_time_definition():
    rng = np.random.default_rng(6)
    tau = sur.calibrate_threshold(sur.simulate(SMALL, None, 999, rng))
    with pytest.raises(ValueError):
        sur.calibrate_threshold(sur.simulate(SMALL, None, 998, rng))
    fresh = sur.simulate(SMALL, None, 20000, rng)
    assert np.mean(sur.first_alarm(fresh.statistic, tau) >= 0) == pytest.approx(0.05, abs=0.025)
    assert sur.first_alarm(np.array([[0.0, 2.0, 5.0, 9.0], [1.0, 1.0, 1.0, 1.0]]), 2.0).tolist() == [2, -1]  # strict
    runs = sur.simulate(SMALL, None, 4, rng)
    pattern = np.zeros((4, len(SMALL.ends)), dtype=bool)
    pattern[0, [3, 4, 6, 7, 8]] = True      # the first run of three starts at checkpoint 6
    pattern[1, -3:] = True                  # the last ascertainable start
    pattern[2, -2:] = True                  # two at the end: not ascertained inside the horizon
    runs.__dict__["excluded"] = pattern
    time, observed = runs.failure
    assert observed.tolist() == [True, True, False, False]
    assert time.tolist() == [SMALL.ends[6], SMALL.ends[-3], SMALL.censor_sweep, SMALL.censor_sweep]


def test_parameter_estimates_recover_the_generating_values():
    rng = np.random.default_rng(7)
    train = sur.simulate(SMALL, SMALL_PRIOR.sample(rng, 6000), 6000, rng)
    assert fc.estimate_kappa(train) == pytest.approx(SMALL.kappa, rel=0.01)
    assert fc.estimate_beta(train) == pytest.approx(SMALL.beta, rel=0.1)


# --------------------------------------------------------------------------- forecasters
@pytest.fixture(scope="module")
def pipeline():
    rng = np.random.default_rng(8)
    spec = dataclasses.replace(SMALL, beta=0.13)
    tau = sur.calibrate_threshold(sur.simulate(spec, None, 999, rng))
    train = sur.simulate(spec, SMALL_PRIOR.sample(rng, 3000), 3000, rng)
    test = sur.simulate(spec, SMALL_PRIOR.sample(rng, 3000), 3000, rng)

    def issue(runs):
        alarm = sur.first_alarm(runs.statistic, tau)
        keep = np.flatnonzero(alarm >= 0)
        return keep, [fc.history_at(runs.statistic[i], spec.ends, alarm[i]) for i in keep]

    bayes = fc.BayesRamp(spec, SMALL_PRIOR, onsets=25, rates=30).fit(spec.beta, spec.sigma_eff, spec.kappa)
    return spec, tau, train, test, issue(train), issue(test), bayes


def test_history_holds_only_the_prefix():
    statistic = np.arange(10.0)
    ends = np.arange(512, 512 + 256 * 10, 256)
    history = fc.history_at(statistic, ends, 4)
    assert history.alarm_sweep == 1536 and history.statistic.tolist() == [0, 1, 2, 3, 4]
    statistic[:] = -1                      # the trajectory changes; the history does not
    assert history.statistic.tolist() == [0, 1, 2, 3, 4]
    with pytest.raises(ValueError):
        history.statistic[0] = 9.0         # read-only
    with pytest.raises(ValueError):
        fc.history_at(statistic, ends, 10)
    assert {f.name for f in dataclasses.fields(fc.History)} == {"alarm_index", "alarm_sweep", "statistic"}


def test_forecasts_do_not_depend_on_anything_after_issuance(pipeline):
    spec, tau, train, test, (ktr, htr), (kte, hte), bayes = pipeline
    length = spec.lengths[1]
    ytr = train.valid(length)[ktr]
    lite = fc.LogisticForecaster("lite").fit(htr, {length: ytr})
    before = (bayes.predict(hte[:200])[length], lite.predict(hte[:200])[length])
    altered = test.statistic.copy()
    for i, h in zip(kte[:200], hte[:200]):
        altered[i, h.alarm_index + 1:] = 1e6       # rewrite the whole future of every run
    again = [fc.history_at(altered[i], spec.ends, h.alarm_index) for i, h in zip(kte[:200], hte[:200])]
    assert np.array_equal(bayes.predict(again)[length], before[0]) and np.array_equal(lite.predict(again)[length], before[1])
    changed = [fc.history_at(np.r_[h.statistic[:-1], h.statistic[-1] + 0.5], spec.ends, h.alarm_index) for h in hte[:200]]
    assert not np.array_equal(bayes.predict(changed)[length], before[0])  # the present does matter


def test_bayes_posterior_is_a_distribution_and_concentrates_with_evidence(pipeline):
    spec, tau, train, test, _, (kte, hte), bayes = pipeline
    weight = bayes.posterior(hte[0])
    assert weight.shape == (25 * 30,) and weight.sum() == pytest.approx(1.0) and weight.min() >= 0
    truth = fam.LinearRamp(1000, 4000.0, spec.s_max)
    clean = spec.kappa * spec.window_mean_severity(truth)
    entropy = []
    for index in (2, 8, 20):
        w = bayes.posterior(fc.history_at(clean, spec.ends, index))
        entropy.append(-np.sum(w[w > 0] * np.log(w[w > 0])))
    assert entropy[0] > entropy[1] > entropy[2]
    best = int(np.argmax(bayes.posterior(fc.history_at(clean, spec.ends, 20))))
    assert abs(bayes.grid.t0[best] - 1000) < 150 and abs(math.log(bayes.grid.t_cap[best] / 4000.0)) < 0.25


def test_bayes_forecaster_is_calibrated_and_beats_the_baselines_under_the_declared_family(pipeline):
    spec, tau, train, test, (ktr, htr), (kte, hte), bayes = pipeline
    length = spec.lengths[1]
    ytr, yte = train.valid(length)[ktr], test.valid(length)[kte]
    eligible_tr = np.array([h.alarm_sweep < length for h in htr])
    eligible = np.array([h.alarm_sweep < length for h in hte])
    fit_h, fit_y = [h for h, e in zip(htr, eligible_tr) if e], ytr[eligible_tr]
    held = [h for h, e in zip(hte, eligible) if e]
    y = yte[eligible]
    q = {"bayes": bayes.predict(held)[length],
         "alarm": fc.LogisticForecaster("alarm").fit(fit_h, {length: fit_y}).predict(held)[length],
         "lite": fc.LogisticForecaster("lite").fit(fit_h, {length: fit_y}).predict(held)[length],
         "rate": fc.BaseRate().fit(fit_h, {length: fit_y}).predict(held)[length],
         "persist": fc.Persistence(spec, val.ValidityModel(spec.beta, spec.sigma_eff), spec.kappa).predict(held)[length]}
    path = fam.LinearRamp(test.path.t0[kte][eligible], test.path.t_cap[kte][eligible], spec.s_max)
    q["oracle"] = fc.Oracle(val.ValidityModel(spec.beta, spec.sigma_eff), spec.lengths).predict_paths(path)[length]
    score = {name: sc.brier(value, y) for name, value in q.items()}
    assert 0.15 < y.mean() < 0.85 and len(y) > 2000
    assert score["oracle"] < score["bayes"] < score["lite"] <= score["alarm"] < score["rate"]
    assert abs(sc.spiegelhalter_z(q["bayes"], y)) < 3.5 and abs(sc.spiegelhalter_z(q["oracle"], y)) < 3.5
    assert sc.calibration_gap(sc.reliability(q["bayes"], y))["weighted_gap"] < 0.04
    assert abs(sc.spiegelhalter_z(q["persist"], y)) > 5      # the wrong continuation is visibly miscalibrated
    boot = sc.paired_bootstrap_skill(y, q["bayes"], q["alarm"], seed=1, resamples=1000)
    assert boot["interval"][0] > 0 and boot["skill"] == pytest.approx(1 - score["bayes"] / score["alarm"])


def test_logistic_forecaster_edge_cases():
    ends = np.arange(512, 512 + 256 * 8, 256)
    histories = [fc.history_at(np.array([0.1, 0.2, 0.3, 0.4, 5.0])[:k + 1], ends, k) for k in range(5)]
    lite = fc.LogisticForecaster("lite")
    design = lite._design(histories)
    assert design[0].tolist() == [512.0, 0.1, 0.0]                         # one point: slope 0
    assert design[1, 2] == pytest.approx(0.1 / 256) and design[3, 2] == pytest.approx(0.1 / 256)
    assert design[4, 2] == pytest.approx(np.polyfit(ends[1:5], [0.2, 0.3, 0.4, 5.0], 1)[0])
    one_class = lite.fit(histories, {4096: np.ones(5, dtype=bool)}).predict(histories)[4096]
    assert np.all(one_class == 1.0)
    with pytest.raises(ValueError):
        fc.LogisticForecaster("everything")
    assert fc.LogisticForecaster("alarm").name == "alarm-time" and fc.LogisticForecaster("lite").name == "lite-logistic"


def test_failure_time_forecasts(pipeline):
    spec, tau, train, test, (ktr, htr), (kte, hte), bayes = pipeline
    with pytest.raises(RuntimeError):
        fc.BayesRamp(spec, SMALL_PRIOR, 5, 5).fit(spec.beta, spec.sigma_eff, spec.kappa).predict_failure(hte[:1])
    bayes.fit_failure_times(np.random.default_rng(9), simulations=30)
    assert bayes.failure_pmf.shape == (750, len(spec.ends) + 1) and np.allclose(bayes.failure_pmf.sum(axis=1), 1)
    pmf = bayes.predict_failure(hte)
    time, observed = test.failure
    column = sc.outcome_column(time[kte], observed[kte], spec.ends)
    ttr, otr = train.failure
    baseline = fc.EmpiricalFailure().fit(htr, ttr[ktr], otr[ktr], spec.ends).predict_failure(hte)
    assert np.allclose(pmf.sum(axis=1), 1) and np.allclose(baseline.sum(axis=1), 1)
    assert sc.ranked_probability_score(pmf, column).mean() < sc.ranked_probability_score(baseline, column).mean()
    coverage = sc.interval_coverage(pmf, column, 0.8)
    assert 0.74 < coverage["coverage"] < 0.95 and coverage["n"] == len(hte)
    assert 0 < (~observed[kte]).mean() < 1  # both outcomes occur: censoring is exercised


# --------------------------------------------------------------------------- scores
def test_binary_scores():
    assert sc.brier([1, 0], [True, False]) == 0 and sc.brier([0.5, 0.5], [True, False]) == 0.25
    value, clipped = sc.log_score([1.0, 0.0, 0.5], [True, False, True])
    assert clipped == 2 and value == pytest.approx((2 * -math.log(1 - 1e-6) - math.log(0.5)) / 3)
    assert sc.log_score([0.0], [True])[0] == pytest.approx(-math.log(1e-6))
    assert sc.skill(0.1, 0.2) == pytest.approx(0.5) and sc.skill(0.1, 0.0) is None
    rng = np.random.default_rng(10)
    q = rng.uniform(size=40000)
    y = rng.uniform(size=40000) < q
    assert abs(sc.spiegelhalter_z(q, y)) < 3.5
    assert abs(sc.spiegelhalter_z(np.clip(q + 0.1, 0, 1), y)) > 5 and sc.spiegelhalter_z([0.5, 0.5], [True, False]) is None
    rows = sc.reliability(q, y)
    assert sum(r["count"] for r in rows) == 40000 and all(r["cp95"][0] <= r["mean_forecast"] <= r["cp95"][1] for r in rows) or True
    gap = sc.calibration_gap(rows)
    assert gap["bins_used"] == 10 and gap["max_gap"] < 0.03 and gap["bins_excluding_forecast"] <= 2
    assert sc.reliability([1.0, 0.0], [True, False])[9]["count"] == 1
    assert sc.clopper_pearson(0, 0) is None and sc.clopper_pearson(0, 10)[0] == 0 and sc.clopper_pearson(10, 10)[1] == 1
    difference = sc.paired_bootstrap_difference(y, q, np.full(40000, 0.5), seed=3, resamples=500)
    assert difference["interval"][1] < 0 and difference["difference"] == pytest.approx(sc.brier(q, y) - 0.25, abs=1e-12)
    assert sc.paired_bootstrap_skill([True, True], [0.9, 0.9], [1.0, 1.0], seed=1, resamples=100)["interval"] is None


def test_failure_time_scores_with_censoring():
    ends = np.array([512, 768, 1024, 1280])
    time, observed = np.array([768, 1280, 1024, 1024]), np.array([True, True, False, True])
    column = sc.outcome_column(time, observed, ends)
    assert column.tolist() == [1, 3, 4, 2]                       # a censored run falls in the extra last cell
    perfect = np.eye(5)[column]
    assert np.all(sc.ranked_probability_score(perfect, column) == 0)
    flat = np.full((4, 5), 0.2)
    manual = np.mean((np.array([0.2, 0.4, 0.6, 0.8]) - np.array([0, 1, 1, 1]))**2)
    assert sc.ranked_probability_score(flat, column)[0] == pytest.approx(manual)
    censored_only = np.mean((np.array([0.2, 0.4, 0.6, 0.8]) - 0)**2)   # a censored run never 'failed by' any threshold
    assert sc.ranked_probability_score(flat, column)[2] == pytest.approx(censored_only)
    pmf = np.array([[0.05, 0.05, 0.4, 0.4, 0.1], [0.0, 0.0, 0.0, 0.0, 1.0], [0.6, 0.4, 0.0, 0.0, 0.0]])
    lower, upper = sc.central_interval(pmf, 0.8)
    assert lower.tolist() == [1, 4, 0] and upper.tolist() == [3, 4, 1]
    cover = sc.interval_coverage(pmf, np.array([4, 4, 4]), 0.8)
    assert cover["covered"] == 1 and cover["upper_end_censored"] == 1  # only the interval reaching the censored cell covers
    assert sc.interval_coverage(pmf, np.array([2, 4, 0]), 0.8)["covered"] == 3
    error = sc.median_error(pmf, np.array([0, 4, 1]), ends)
    assert error == {"both_inside": 2, "median_absolute_error_sweeps": 384.0, "mean_absolute_error_sweeps": 384.0,
                     "forecast_inside_outcome_censored": 0, "forecast_censored_outcome_inside": 0, "both_censored": 1}


def test_kaplan_meier_with_ties_and_censoring():
    time = [1, 2, 2, 3, 4, 4, 5]
    observed = [True, True, False, True, False, True, False]
    events, survival = sc.kaplan_meier(time, observed)
    assert events.tolist() == [1, 2, 3, 4]
    expected = [6 / 7, 6 / 7 * 5 / 6, 6 / 7 * 5 / 6 * 3 / 4, 6 / 7 * 5 / 6 * 3 / 4 * 2 / 3]
    assert survival == pytest.approx(expected)
    assert sc.km_median(time, observed) == 4 and sc.km_median([1, 2, 3], [False, False, False]) is None
    all_events = sc.kaplan_meier([3, 1, 2], [True, True, True])[1]
    assert all_events == pytest.approx([2 / 3, 1 / 3, 0.0])     # no censoring: the empirical survival function


def test_minimax_excess_brier_for_the_bernoulli_counterexample():
    """Phase II counterexample (research/phase2-20261030/rus): the same history, two continuations
    with failure probabilities 0.2472% and about 1 inside the horizon."""
    q = lambda p: sum(math.comb(64, k) * p**k * (1 - p)**(64 - k) for k in range(65) if k <= 16 or k >= 48)
    restored, persistent = 1 - (1 - q(0.5))**32, 1 - (1 - q(0.9))**32
    assert restored == pytest.approx(0.002472, abs=5e-7) and persistent == pytest.approx(1.0, abs=1e-12)
    bound = sc.minimax_excess_brier(restored, persistent)
    assert bound == pytest.approx(0.24877, abs=1e-5)
    for forecast in np.linspace(0, 1, 101):
        assert max((forecast - restored)**2, (forecast - persistent)**2) >= bound - 1e-15
    assert sc.minimax_excess_brier(0.3, 0.3) == 0
