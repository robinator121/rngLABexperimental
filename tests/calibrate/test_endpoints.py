import numpy as np
import pytest
from rnglab17.endpoints import *

def test_rank():assert rank_threshold(np.arange(399),.05)==379

def test_small_reject():
    with pytest.raises(ValueError):rank_threshold([1],.01)
def test_persistence():
    x=np.array([[[0.],[0.],[5.],[0.],[5.],[5.]]]);r=window_metrics(x,0,1,1,[0],[1]);assert first_events(r['persistent_error'],r['times'],[1])[0,0]==6

def test_strict_threshold():
    p=np.ones((3,4,1));assert np.all(first_events(p,np.arange(4),[1])==-1)
def test_rate_empty():assert rate(0,0)['rate'] is None

def test_batch_constant():assert batch_intervals(np.ones((3,64)),1,8)['covered'].all()
def test_marginal_not_simultaneous():assert rate(10,512,alpha=.05/12,one_sided=True)['interval'][1]>rate(10,512,one_sided=True)['interval'][1]

def test_paired_no_damage():assert causal_loss(np.zeros(64),np.zeros(64),0,.1)['paired_excess_squared_loss']==0

def test_calibration_independent():
    rng=np.random.default_rng(170001);threshold=rank_threshold(rng.normal(size=399));holdout=rng.normal(size=10000);assert .025<np.mean(holdout>threshold)<.08
