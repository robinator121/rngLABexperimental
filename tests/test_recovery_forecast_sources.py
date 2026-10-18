"""Independent regression controls for the audited source/forecast failures."""
import hashlib
import math
import numpy as np
import pytest
from rnglab14.sources.base import RandomSource,SourceManifest,SourceBudgetExceeded
from rnglab14.sources import PCG64Source,DegradedFaultSource,XORParityConditioner,GeigerMullerPulseSource,VonNeumannConditioner
from rnglab14.forecast import (NoChangeBaselinePredictor,LinearExtrapolationPredictor,
    FixedHorizonPredictor,evaluate_forecast_on_trajectory,summarize_forecasting_performance)
from rnglab14.forecast.survival import KaplanMeierBaseline,CensoringKaplanMeier,ipcw_brier


class Bytes(RandomSource):
    def __init__(self,data):
        super().__init__(SourceManifest('control','synthetic',True,False));self.data=bytes(int(v) for v in data);self.at=0
    def _read_bytes(self,n):
        data=bytes(self.data[(self.at+i)%len(self.data)] for i in range(n));self.at+=n;return data


def test_all_consumption_paths_hash_emitted_bytes_exactly_once():
    src=PCG64Source(seed=18);raw=src.read_bytes(64)
    assert src.stream_digest==hashlib.sha256(raw).hexdigest()
    other=PCG64Source(seed=18)
    for _ in range(8):other.random_float()
    assert other.stream_digest==src.stream_digest and other.consumed_bytes==64
    for bad in [-1,True,1.5]:
        with pytest.raises(ValueError):other.read_bytes(bad)
    limited=PCG64Source(seed=1,max_bytes=1)
    with pytest.raises(SourceBudgetExceeded):limited.read_bytes(2)
    assert limited.consumed_bytes==0
    for budget in [math.nan,math.inf,-1,True,1.5]:
        with pytest.raises(ValueError):PCG64Source(seed=1,max_bytes=budget)


def test_rho_fault_reaches_bytes_bits_and_floats():
    def source(rho):return DegradedFaultSource(PCG64Source(seed=71),autocorrelation_rho=rho)
    assert source(.8).read_bytes(128)!=source(0).read_bytes(128)
    assert [source(.8).random_float()]!=[source(0).random_float()]
    bits=np.asarray(source(.8).read_bits(32000))
    assert .47<float(bits.mean())<.53
    assert .75<float(np.corrcoef(bits[:-1],bits[1:])[0,1])<.85
    s=source(.8);raw=s.read_bytes(256)
    assert s.stream_digest==hashlib.sha256(raw).hexdigest()
    assert s.manifest.metadata['rho_mechanism']=='stationary_binary_markov'


def test_joint_dependent_xor_zero_counterexample_and_constant_inputs_rejected():
    rng=np.random.default_rng(0);a=rng.integers(0,256,8192,dtype=np.uint8);b=rng.integers(0,256,8192,dtype=np.uint8)
    with pytest.raises(ValueError,match='constant combined'):
        XORParityConditioner([Bytes(a),Bytes(b),Bytes(a^b)])
    with pytest.raises(ValueError,match='constant input'):
        XORParityConditioner([Bytes([0]),Bytes([0])])
    # Disabling diagnostics is explicit and NEVER changes the assurance label.
    unsafe=XORParityConditioner([Bytes(a),Bytes(b),Bytes(a^b)],check_independence=False)
    assert set(unsafe.read_bytes(128))=={0}
    assert unsafe.manifest.metadata['independence_established'] is False
    assert unsafe.manifest.metadata['entropy_certified'] is False


def test_nonlinear_uncorrelated_dependence_is_not_labeled_independent():
    x=np.column_stack([np.arange(64),255-np.arange(64)]).ravel();y=np.abs(2*x-255)
    assert abs(np.corrcoef(x,y)[0,1])<1e-12
    c=XORParityConditioner([Bytes(x),Bytes(y)])
    assert c.manifest.metadata['independence_established'] is False


def test_stuck_von_neumann_has_finite_budget_and_live_emulation_refused():
    with pytest.raises(ValueError,match='budget exhausted'):
        VonNeumannConditioner(Bytes([0]),max_raw_bits_per_output_bit=8).read_bytes(1)
    with pytest.raises(ValueError,match='not implemented'):
        GeigerMullerPulseSource(is_synthetic=False)


def test_von_neumann_buffers_contiguous_bits_and_reports_physical_consumption():
    source=Bytes([0b01010101]);vn=VonNeumannConditioner(source)
    assert vn.read_bytes(10)==bytes(10)
    assert vn.raw_bits_acquired==source.consumed_bytes*8==160
    assert vn.raw_bits_consumed==160
    assert vn.extraction_efficiency==.5  # deterministic01 is not IID Bernoulli


def test_xor_constant_run_guard_is_not_bypassed_by_float_sized_reads():
    rng=np.random.default_rng(0);a=rng.integers(0,256,256,dtype=np.uint8);b=rng.integers(0,256,256,dtype=np.uint8);c=rng.integers(0,256,256,dtype=np.uint8)
    c[128:]=a[128:]^b[128:]
    x=XORParityConditioner([Bytes(a),Bytes(b),Bytes(c)])
    assert x.random_float()==0
    with pytest.raises(ValueError,match='constant output run'):x.random_float()


def test_geiger_dead_time_modes_against_distinct_analytical_mean_laws():
    np_source=GeigerMullerPulseSource(mean_cps=1,dead_time_seconds=.5,seed=814)
    p_source=GeigerMullerPulseSource(mean_cps=1,dead_time_seconds=.5,seed=814,paralyzable=True)
    assert abs(np.mean([np_source._next_event_interval() for _ in range(8000)])-1.5)<.05
    assert abs(np.mean([p_source._next_event_interval() for _ in range(8000)])-math.exp(.5))<.05


def test_censored_bound_beyond_followup_is_unresolved_not_covered():
    f=evaluate_forecast_on_trajectory(NoChangeBaselinePredictor(.4),[.1]*50,.4,[10])[0]
    assert f.observed_followup_samples==40 and f.status_note=='CENSORED_UNRESOLVED'
    assert f.covered_by_lower_bound is None and f.is_false_reassurance is None
    report=summarize_forecasting_performance([f])
    assert report.unresolved_bound_count==1 and report.lower_bound_coverage_95 is None
    assert report.event_mae is None and report.is_well_calibrated is None


def test_abstention_and_post_event_never_score_as_forecast_success():
    fs=evaluate_forecast_on_trajectory(LinearExtrapolationPredictor(.4),[0.]*100,.4,[20,40,60])
    rep=summarize_forecasting_performance(fs)
    assert rep.abstention_count==3 and rep.emitted_count==0
    assert rep.event_mae is None and rep.lower_bound_coverage_95 is None and rep.is_well_calibrated is None
    fs=evaluate_forecast_on_trajectory(NoChangeBaselinePredictor(.4),[0.]*49+[.5]*51,.4,[20,50,60])
    assert [f.current_t for f in fs]==[20]


def test_confirmed_horizon_and_empty_population_semantics():
    fs=evaluate_forecast_on_trajectory(FixedHorizonPredictor(10,.4),[.1]*50,.4,[10],trajectory_id='one')
    assert fs[0].status_note=='CENSORED_BOUND_SUPPORTED'
    assert summarize_forecasting_performance(fs).unique_trajectory_ids_count==1
    with pytest.raises(ValueError):summarize_forecasting_performance([])
    with pytest.raises(ValueError):evaluate_forecast_on_trajectory(FixedHorizonPredictor(1,.4),[0.,math.nan],.4,[1])
    class Bad(FixedHorizonPredictor):
        def predict(self,history,current_t):return math.nan,math.nan,False
    with pytest.raises(ValueError):evaluate_forecast_on_trajectory(Bad(1,.4),[0.]*10,.4,[1])


def test_survival_baseline_and_brier_match_independent_arithmetic():
    model=KaplanMeierBaseline.fit([1,2,3,4],[True,False,True,False])
    assert model.survival(1)==.75 and model.survival(3)==.375
    # No training censoring before horizon3: IPCW reduces to ordinary Brier.
    g=CensoringKaplanMeier.fit([8,9,10,11],[True]*4)
    rep=ipcw_brier([1,4,5],[True,True,False],[.8,.2,.3],3,g,trajectory_ids=['a','b','c'])
    assert math.isclose(rep['brier_ipcw'],(.2**2+.2**2+.3**2)/3)
    with pytest.raises(ValueError):ipcw_brier([1,4],[True,False],[.5,.5],3,g,['same','same'])
    with pytest.raises(ValueError):model.failure_probability(10)
    # One censored training observation at1 => G(2-)=.75, G(3-)=.75.
    g2=CensoringKaplanMeier.fit([1,4,5,6],[False,True,True,True])
    rep=ipcw_brier([2,5,1],[True,False,False],[.7,.4,.2],3,g2)
    assert math.isclose(rep['brier_ipcw'],(.3**2/.75+.4**2/.75+0)/3)


def test_event_first_discrete_censoring_ties_against_independent_cartesian_law():
    # Exact independent Cartesian product: T in{1,3}, C in{1,2,4}.
    observed=[True,True,True,False,False,True];durations=[1,1,1,1,2,3]
    g=CensoringKaplanMeier.fit(durations,observed)
    assert math.isclose(g.survival(2,before=True),2/3)
    assert math.isclose(g.survival(3,before=True),1/3)
    for horizon in [2,3]:
        rep=ipcw_brier(durations,observed,[.5]*6,horizon,g)
        assert math.isclose(rep['brier_ipcw'],.25)
