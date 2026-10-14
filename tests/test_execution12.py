import json
import numpy as np
import pytest
from rnglab12.execution import DESIGN, simulate, run_tasks

def test_healthy_baseline_has_multiple_genuine_pre_fault_looks():
    row=simulate({'application':'integration','split':'pilot','group':'timing'})
    assert sum(t<=4096 for t in row['times'])>=3
    assert len(row['errors'])==len(row['source_scores'])==len(row['times'])
    assert row['times'][-1]==DESIGN['steps']
    assert row['source_family']=='synthetic_fault'

def test_prospective_features_do_not_depend_on_future_intervention():
    base={'application':'integration','split':'pilot','group':'past-only','issued_at':3072}
    a=simulate({**base,'linear_scale':{'start_row':4096,'duration_rows':4096,'final_severity':.1}})
    b=simulate({**base,'linear_scale':{'start_row':4096,'duration_rows':4096,'final_severity':.8}})
    assert a['issuance']==b['issuance']
    assert a['uncertainty']['mean']!=b['uncertainty']['mean']

def test_resume_refuses_tampered_completed_result(tmp_path):
    task={'application':'integration','split':'pilot','group':'resume-test',
          'design':{**DESIGN,'steps':1024,'burn_in':128,'rolling':256}}
    rows=run_tasks(tmp_path,[task],workers=1)
    result=next((tmp_path/'results').glob('*.json'))
    value=json.loads(result.read_text());value['uncertainty']['mean']+=1
    result.write_text(json.dumps(value))
    with pytest.raises(ValueError,match='hash mismatch'):run_tasks(tmp_path,[task],workers=1)

def test_proposal_acceptance_endpoint_control_is_material():
    base={'application':'vmc','split':'pilot','group':'endpoint-control',
          'fault':{'mode':'zero_burst','scope':'all','start_row':4096,'duration_rows':512}}
    zero=simulate(base)
    positive=simulate({**base,'positive_acceptance_control':True})
    assert abs(zero['uncertainty']['mean']-zero['reference'])>1000*abs(positive['uncertainty']['mean']-positive['reference'])
