import pytest
from rnglab14.api import handle_request

def test_reference_api_and_unknown_operation():
    result=handle_request({'schema':'rnglab.headless.v1','operation':'harmonic_reference','parameters':{'beads':16,'beta':2.}})
    assert result['result']['finite_bead_energy']==pytest.approx(.655474152298)
    with pytest.raises(ValueError):handle_request({'schema':'rnglab.headless.v1','operation':'execute_code','parameters':{}})

def test_forecast_api_preserves_unresolved_not_success():
    response=handle_request({'schema':'rnglab.headless.v1','operation':'forecast_evaluate','parameters':{'trajectory':[.1]*10,'threshold':.4,'landmarks':[5],'horizon':50}})
    assert response['result']['summary']['unresolved_bound_count']==1
    assert response['result']['summary']['is_well_calibrated'] is None
