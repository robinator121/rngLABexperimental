import argparse,json
from rnglab12 import studies
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['calibration','healthy-validation','periodic','faults','sources','rus']);a=p.parse_args()
    fn={'calibration':studies.calibration_study,'healthy-validation':studies.healthy_validation,'periodic':studies.periodic_study,'faults':studies.fault_study,'sources':studies.source_study,'rus':studies.rus_study}[a.stage]
    result=fn();print(json.dumps(result,indent=2))
