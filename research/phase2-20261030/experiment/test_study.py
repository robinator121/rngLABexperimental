"""Independent mathematical and implementation checks for the Phase II study."""
import math
import sys
from pathlib import Path
import numpy as np
import pytest
from scipy.stats import t
import study


def test_finite_reference_via_precision_spectrum_and_logdet_derivative():
    for m,beta in [(4,.5),(8,2.),(16,2.),(32,4.)]:
        ref=study.reference(m,beta)
        np.testing.assert_allclose(np.linalg.eigvalsh(ref['precision']),np.sort(ref['normal_mode_precision']),rtol=1e-12,atol=1e-12)
        assert abs(ref['means']['TD']-ref['finite_energy'])<1e-12
        assert abs(ref['means']['CV']-ref['finite_energy'])<1e-12
        # Z is (M/beta)^(M/2) det(Q)^(-1/2), apart from beta-independent factors.
        def logz(b):
            q=study.reference(m,b)['precision']
            return .5*m*math.log(m/b)-.5*np.linalg.slogdet(q)[1]
        h=1e-5
        assert abs(-(logz(beta+h)-logz(beta-h))/(2*h)-ref['finite_energy'])<1e-8


def test_action_difference_matches_full_gaussian_log_density():
    g=np.random.default_rng(14091);ref=study.reference();m=study.M;dt=study.BETA/m
    for _ in range(30):
        x=g.normal(size=m);j=int(g.integers(m));new=x[j]+g.normal();old=x[j]
        left,right=x[(j-1)%m],x[(j+1)%m]
        local=((new-left)**2+(right-new)**2-(old-left)**2-(right-old)**2)/(2*dt)+dt*(new*new-old*old)/2
        xp=x.copy();xp[j]=new
        dense=.5*(xp@ref['precision']@xp-x@ref['precision']@x)
        assert abs(local-dense)<2e-12


@pytest.mark.parametrize('refresh,reverse',[(False,False),(True,False),(False,True)])
def test_compiled_kernel_matches_separately_written_numpy_kernel(refresh,reverse):
    initial=np.linspace(-.4,.6,16)
    compiled=study.local_kernel(np.random.default_rng(11),np.random.default_rng(12),np.random.default_rng(13),initial,40,2.,.4,refresh,reverse)
    scalar=study.numpy_kernel(np.random.default_rng(11),np.random.default_rng(12),np.random.default_rng(13),initial,40,2.,.4,refresh,reverse)
    np.testing.assert_allclose(compiled[0],scalar[0],rtol=0,atol=1e-13)
    assert compiled[2]==scalar[1]
    np.testing.assert_array_equal(initial,np.linspace(-.4,.6,16))


def test_energy_estimators_match_upstream_and_exact_gaussian_moments():
    upstream=Path(__file__).resolve().parents[3]/'work/rnglab-main-snapshot/src'
    sys.path.insert(0,str(upstream))
    from rnglab14.models.pimc_harmonic import harmonic_energy_estimators
    ref=study.reference();g=np.random.default_rng(39401);n=200000
    paths=g.normal(size=(n,16))@ref['cholesky'].T
    c=paths.mean(1);shape=((paths-c[:,None])**2).mean(1)
    td=16/4-16*((np.roll(paths,-1,axis=1)-paths)**2).sum(1)/8+.5*(shape+c*c)
    cv=.25+.5*c*c+shape
    up_td,up_cv=harmonic_energy_estimators(paths,2.)
    np.testing.assert_allclose(td,up_td,atol=1e-13,rtol=1e-13)
    np.testing.assert_allclose(cv,up_cv,atol=1e-13,rtol=1e-13)
    for name,values in [('TD',td),('CV',cv)]:
        assert abs(values.mean()-ref['means'][name])<6*math.sqrt(ref['variances'][name]/n)
        assert abs(values.var()-ref['variances'][name])<.025*ref['variances'][name]


def test_batch_variance_on_stationary_AR1_matches_exact_finite_sample_expectation():
    # The theoretical expectation includes covariance between batches; it is
    # below Var(mean), which is precisely the potential failure being studied.
    r,n,b,phi=2048,1024,32,.9
    g=np.random.default_rng(59131);y=np.empty((r,n));y[:,0]=g.normal(size=r)
    for j in range(1,n):y[:,j]=phi*y[:,j-1]+math.sqrt(1-phi**2)*g.normal(size=r)
    def variance_mean(length):
        k=np.arange(1,length)
        return (length+2*np.sum((length-k)*phi**k))/length**2
    expected=(variance_mean(b)-variance_mean(n))/(n//b-1)
    estimate=study.variance_methods(y,(b,))[f'BM{b}'][0].mean()
    assert abs(estimate/expected-1)<.04
    assert estimate<variance_mean(n)*.9


def test_batch_abstention_and_t_critical_match_declared_method():
    y=np.random.default_rng(401).normal(size=(12,3072))
    methods=study.variance_methods(y,(32,128,512))
    assert 'BM512' not in methods  # only six blocks; frozen min eight
    assert methods['BM128'][1]==float(t.ppf(.975,23))


def test_centroid_refresh_preserves_shape_action_and_target_factorization():
    g=np.random.default_rng(5001);ref=study.reference();x=g.normal(size=16);c=x.mean();target=.7
    xp=x+target-c
    np.testing.assert_allclose(x-c,xp-xp.mean(),atol=1e-15)
    np.testing.assert_allclose(np.diff(np.r_[x,x[0]]),np.diff(np.r_[xp,xp[0]]),atol=1e-15)
    change=.5*(xp@ref['precision']@xp-x@ref['precision']@x)
    assert abs(change-study.BETA*(target*target-c*c)/2)<1e-12
