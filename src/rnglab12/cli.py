"""Local, explicit finite-stage CLI; never starts services or model sessions."""
import argparse, hashlib, json
from pathlib import Path
from .execution import ROOT

def verify(root=ROOT):
    manifest=root/'SHA256SUMS'
    if not manifest.exists():raise ValueError('Release is not sealed yet')
    count=0
    for line in manifest.read_text().splitlines():
        expected,name=line.split('  ',1);p=(root/name).resolve()
        if not p.is_relative_to(root.resolve()) or not p.is_file():raise ValueError('Invalid or missing manifest path')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:raise ValueError('Integrity mismatch: '+name)
        count+=1
    return {'status':'PASS','files':count,'scope':'byte integrity only; not scientific correctness'}

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    sub.add_parser('status');sub.add_parser('verify')
    run=sub.add_parser('run');run.add_argument('stage',choices=['calibration','healthy-validation','periodic','faults','sources','rus'])
    args=p.parse_args(argv)
    if args.command=='verify':out=verify()
    elif args.command=='status':
        out={'live_process_status':'NOT_PROBED; saved receipts never imply a running service',
             'receipts':{str(x.relative_to(ROOT)):json.loads(x.read_text()).get('status','RECORDED') for x in sorted((ROOT/'runs').glob('*/receipt.json'))},
             'automatic_workers_or_provider_calls':False}
    else:
        from . import studies
        fn={'calibration':studies.calibration_study,'healthy-validation':studies.healthy_validation,
            'periodic':studies.periodic_study,'faults':studies.fault_study,'sources':studies.source_study,'rus':studies.rus_study}[args.stage]
        out=fn()
    print(json.dumps(out,indent=2));return 0

if __name__=='__main__':raise SystemExit(main())
