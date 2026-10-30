#!/usr/bin/env python3
"""Read-only integrity checks for the lightweight public Phase II payload."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

PHASE=Path(__file__).resolve().parents[1]


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(2**20),b''):
            h.update(block)
    return h.hexdigest()


def main():
    manifest=json.loads((PHASE/'MANIFEST.json').read_text())
    failures=[]
    for row in manifest['files']:
        path=PHASE/row['path']
        if not path.is_file() or path.stat().st_size!=row['bytes'] or sha(path)!=row['sha256']:
            failures.append('file integrity: '+row['path'])
    cp=PHASE/'experiment/confirmation-v1/protocol.json'
    fp=PHASE/'experiment/fault-study-v1/protocol.json'
    confirmation=json.loads(cp.read_text());fault=json.loads(fp.read_text())
    pairs=[(PHASE/'experiment/study.py',confirmation['script_sha256']),
           (PHASE/'experiment/study.py',fault['study_code_sha256']),
           (PHASE/'experiment/fault_study.py',fault['code_sha256']),
           (PHASE/'experiment/confirmation-v1/primary_decision.json',fault['healthy_gate_sha256'])]
    for relative,expected in manifest['frozen_protocols'].items():pairs.append((PHASE/relative,expected))
    for path,expected in pairs:
        if sha(path)!=expected:failures.append('frozen link: '+str(path.relative_to(PHASE)))
    if not json.loads((PHASE/'experiment/confirmation-v1/primary_decision.json').read_text())['fault_study_allowed']:
        failures.append('recorded healthy gate did not pass')
    if list(PHASE.rglob('*.npz')):failures.append('unexpected large NPZ in lightweight payload')
    lines=(PHASE/'SHA256SUMS').read_text().splitlines()
    if lines!=[r['sha256']+'  '+r['path'] for r in manifest['files']]:
        failures.append('public checksum inventory differs from manifest')
    print(json.dumps({'passed':not failures,'distributed_files_checked':len(manifest['files']),
                      'frozen_links_checked':len(pairs),'omitted_traces':len(manifest['omitted_traces']),
                      'simulations_executed':False,'historical_132_check_audit_rerun':False,
                      'failures':failures},indent=2))
    return 1 if failures else 0


if __name__=='__main__':raise SystemExit(main())
