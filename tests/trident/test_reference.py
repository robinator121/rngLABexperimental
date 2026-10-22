import numpy as np
import pytest
from rnglab16.core import stationary,two_state,variance_mean,enumeration,energy,augmented,ising,vmc,pimc
from rnglab14.models.pimc_harmonic import harmonic_energy_estimators

def test_stationary():
    for a,b in [(.1,.2),(.9,.9),(.001,.002)]:assert np.allclose(stationary(two_state(a,b)),[b/(a+b),a/(a+b)])
def test_reject_nonunique():
    with pytest.raises(ValueError):stationary(np.eye(2))
def test_periodic():assert np.allclose(stationary([[0,1],[1,0]]),[.5,.5])
def test_variance():assert variance_mean(.5,.5,100)==.25/100
def test_small_enumeration():
    # Independent hand-counted L2 spectrum: degeneracies E=-8:2, E=0:12, E=8:2.
    T=2.;z=4*np.cosh(8/T)+12
    assert np.isclose(enumeration(2,T)['energy'],-8*np.sinh(8/T)/z)
    assert energy(np.ones((3,3)))==-18
    assert np.isclose(enumeration(3,1e6)['energy'],0,atol=1e-5)
def test_augmented():
    for rho in [0,.5,.95]:
        P=augmented(rho=rho);pi=stationary(P)
        assert np.max(abs(pi@P-pi))<1e-12
        assert np.allclose(pi.reshape(2,8).sum(0),np.ones(8)/8)
def test_consumption():
    x,u,_=ising(3,2.3,3,100,411)
    y,v,_=ising(3,2.3,3,100,411,fault='bias')
    assert np.mean(v)>np.mean(u)+.1
    assert np.any(x!=y)
def test_masking():
    x=vmc(8,1000,53,'bias');assert np.all(x[:,:,0]==.5)
def test_pimc_derivative():
    x=np.array([.1,-.7,.3,1]);beta=2.;h=1e-5
    def S(b):return len(x)*np.sum((np.roll(x,1)-x)**2)/(2*b)+b*np.mean(x*x)/2
    independent=len(x)/(2*beta)+(S(beta+h)-S(beta-h))/(2*h)
    assert np.isclose(harmonic_energy_estimators(x,beta)[0],independent)
def test_paired_replay():
    a=ising(2,2,2,40,12)[0];b=ising(2,2,2,40,12)[0];assert np.array_equal(a,b)

def test_single_class_forecast_rejected():
    from rnglab16.core import fit_logistic
    with pytest.raises(ValueError):fit_logistic(np.ones((10,2)),np.zeros(10))
