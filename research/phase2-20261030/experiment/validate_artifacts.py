"""Read-only, independent Phase II artifact validation; never import study.py."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
import time

import numpy as np
from scipy.stats import t, beta as beta_dist

ROOT = Path(__file__).resolve().parent
PHASE = ROOT.parent
ROLES = ('initial', 'update', 'centroid', 'scan', 'independent_paths')
CONFIRM_PROTOCOL_HASH = '468e816c1de68bb8975303b4f725cbffbc3ab34bc34ceef505ab84046025e984'
FAULT_PROTOCOL_HASH = '72ede91d8c4e82caa0454e4d50351272ee26a9e9b161bf4f694bc2c7f70614a1'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(2**20), b''):
            h.update(chunk)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def seed_digest(stage, role, replica):
    message = f'rnglab-phase2-v1/{stage}/{role}/{replica:05d}'
    return hashlib.sha256(message.encode()).hexdigest()


def exact_energy(beads=16, inverse_temperature=2.):
    z = inverse_temperature/(2*beads)
    return .5/(math.tanh(beads*math.asinh(z))*math.sqrt(1+z*z))


def independent_bm512(raw, burn, retained, inverse_temperature=2.):
    x = raw[:, burn:burn+retained]
    if x.shape[1] != retained or retained % 512:
        raise ValueError('Unexpected retained window or incomplete batches')
    # Derive directly from recorded centroid and squared radius of gyration.
    cv = 1/(2*inverse_temperature) + .5*x[:, :, 1]**2 + x[:, :, 2]
    blocks = retained//512
    batch_means = cv.reshape(len(cv), blocks, 512).mean(axis=2)
    mcvar = batch_means.var(axis=1, ddof=1)/blocks
    means = cv.mean(axis=1)
    critical = float(t.ppf(.975, blocks-1))
    valid = np.isfinite(mcvar) & (mcvar>0)
    covered = valid & (abs(means-exact_energy()) <= critical*np.sqrt(np.maximum(mcvar, 0)))
    return means, mcvar, covered, valid, critical


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT/'artifact-validation.json')
    args = parser.parse_args()
    output = args.output
    log = output.with_suffix('.log')
    if output.exists() or log.exists():
        raise FileExistsError('Preserve prior validation outputs; use a new --output path')
    started = time.monotonic()
    checks, manifests, coverage_rows, prefix_rows = [], [], [], []

    def check(name, ok, **detail):
        row = {'check':name, 'passed':bool(ok), **detail}
        checks.append(row)
        print(('PASS ' if ok else 'FAIL ')+name, flush=True)

    try:
        baseline_path = PHASE/'baseline/baseline_verification.json'
        baseline = load(baseline_path)
        baseline_manifest = load(PHASE/'baseline/manifest.json')
        check('baseline audit JSON matches its manifest',
              digest(baseline_path)==baseline_manifest['files']['baseline_verification.json']['sha256'])
        check('baseline lists exactly ten original Phase I artifacts', len(baseline['source_files'])==10)
        for name, entry in baseline['source_files'].items():
            path = PHASE.parent/name
            actual = digest(path)
            check('original preserved: '+name, actual==entry['sha256'] and path.stat().st_size==entry['bytes'],
                  expected_sha256=entry['sha256'], actual_sha256=actual, bytes=path.stat().st_size)

        cp = ROOT/'confirmation-v1/protocol.json'
        fp = ROOT/'fault-study-v1/protocol.json'
        confirm, fault = load(cp), load(fp)
        check('frozen confirmation protocol hash', digest(cp)==CONFIRM_PROTOCOL_HASH, actual_sha256=digest(cp))
        check('frozen fault protocol hash', digest(fp)==FAULT_PROTOCOL_HASH, actual_sha256=digest(fp))
        study_hash, fault_hash = digest(ROOT/'study.py'), digest(ROOT/'fault_study.py')
        check('confirmation code hash', study_hash==confirm['script_sha256'])
        check('fault main study code hash', study_hash==fault['study_code_sha256'])
        check('fault adapter code hash', fault_hash==fault['code_sha256'])
        decision_path = ROOT/'confirmation-v1/primary_decision.json'
        decision = load(decision_path)
        check('fault prerequisite healthy decision hash', digest(decision_path)==fault['healthy_gate_sha256'])
        check('fault prerequisite healthy gate passed', decision['fault_study_allowed'] is True)
        check('confirmation model and primary design', confirm['M']==16 and confirm['beta']==2.
              and confirm['retained']==8192 and confirm['burn']==0 and confirm['replicas']==1024)
        check('fault model, counts and onset', fault['model']['M']==16 and fault['model']['beta']==2.
              and fault['calibration_replicas']==399 and fault['test_replicas']==512
              and fault['total_sweeps']==10240 and fault['fault_start_sweep_1based']==2049)

        expected = {
            'exploratory-v1': ['centroid_refresh', 'iid_gaussian', 'local_displaced_d04',
                'local_stationary_d015', 'local_stationary_d04', 'local_stationary_d10',
                'local_uniform_d04', 'local_zero_d04', 'reversible_scan'],
            'confirmation-v1': ['local_stationary', 'centroid_refresh', 'iid_gaussian'],
            'fault-study-v1': ['healthy_calibration', 'healthy', 'acceptance_bias', 'channel_coupling'],
        }
        global_seeds, stage_seeds, collisions = {}, {}, []
        total_seed_records = 0
        for folder, conditions in expected.items():
            for condition in conditions:
                mp = ROOT/folder/(condition+'.manifest.json')
                m = load(mp)
                stage = m['stage']
                wanted_stage = ('fault-calibration-v1' if condition=='healthy_calibration' else 'fault-test-v1') if folder=='fault-study-v1' else folder
                check('stage identity: '+folder+'/'+condition, stage==wanted_stage)
                params = m['parameters']
                path = ROOT/folder/m['output']
                sha = digest(path)
                check('NPZ checksum: '+folder+'/'+condition, sha==m['sha256'], actual_sha256=sha, bytes=path.stat().st_size)
                check('manifest code hash: '+folder+'/'+condition, m['script_sha256']==study_hash)
                check('manifest condition identity: '+folder+'/'+condition,
                      m['experiment_id']==condition and params['id']==condition and m['output']==condition+'.npz')
                seed_errors = []
                if len(m['seeds'])!=params['replicas']:
                    seed_errors.append('replica count')
                for i, row in enumerate(m['seeds']):
                    if row['replica']!=i:
                        seed_errors.append(f'replica index {i}')
                    if set(row)!={'replica', *ROLES}:
                        seed_errors.append(f'role keys {i}')
                    for role in ROLES:
                        actual = row[role]
                        key = (stage, role, i)
                        if actual!=seed_digest(*key):
                            seed_errors.append(f'{i}/{role}')
                        prior = global_seeds.get(actual)
                        if prior is not None and prior!=key:
                            collisions.append({'digest':actual, 'first':prior, 'second':key})
                        global_seeds[actual] = key
                        stage_seeds.setdefault(stage, set()).add(actual)
                        total_seed_records += 1
                check('seed rule and indices: '+folder+'/'+condition, not seed_errors, errors=seed_errors)
                with np.load(path, allow_pickle=False) as data:
                    raw, sources, acceptance = data['raw'], data['sources'], data['acceptance']
                    wanted_raw = (params['replicas'], params['steps'], 3)
                    wanted_sources = wanted_raw if params.get('source_record') else (0, 0, 0)
                    check('NPZ array domain: '+folder+'/'+condition,
                          raw.shape==wanted_raw and raw.dtype==np.float64 and np.isfinite(raw).all()
                          and np.all(raw[:, :, 2]>=0) and sources.shape==wanted_sources
                          and np.isfinite(sources).all() and acceptance.shape==(params['replicas'],)
                          and np.isfinite(acceptance).all() and np.all((acceptance>=0)&(acceptance<=1)))
                manifests.append({'folder':folder, 'condition':condition, 'stage':stage,
                                  'npz_sha256':sha, 'manifest_sha256':digest(mp), 'replicas':params['replicas'], 'steps':params['steps']})
        check('all seed domains collision free; intentional same-stage condition pairing allowed',
              not collisions, collisions=collisions, unique_seed_count=len(global_seeds),
              recorded_seed_count=total_seed_records, unique_per_stage={k:len(v) for k,v in stage_seeds.items()},
              meaning='Namespace separation and unique derived seeds; not proof of statistically independent entropy')

        for folder, conditions, burn, retained, label in [
            ('confirmation-v1', expected['confirmation-v1'], 0, 8192, 'primary'),
            ('fault-study-v1', ['healthy', 'acceptance_bias', 'channel_coupling'], 2048, 8192, 'postonset')]:
            analysis = load(ROOT/folder/'analysis.json')
            with (ROOT/folder/'replicate_endpoints.csv').open(newline='') as stream:
                endpoint = {(row['condition'], int(row['replica'])):row for row in csv.DictReader(stream)
                            if row['observable']=='CV' and row['method']=='BM512' and row['window']==label}
            for condition in conditions:
                with np.load(ROOT/folder/(condition+'.npz'), allow_pickle=False) as data:
                    means, mcvar, covered, valid, critical = independent_bm512(data['raw'], burn, retained)
                saved = [r for r in analysis if r['condition']==condition and r['window']==label
                         and r['observable']=='CV' and r['method']=='BM512']
                if len(saved)!=1:
                    raise ValueError(f'Expected one analysis row for {folder}/{condition}')
                row = saved[0]
                count = int(covered.sum())
                empvar = float(means.var(ddof=1))
                ratio = float(mcvar[valid].mean()/empvar)
                bias = float(means.mean()-exact_energy())
                ci = [0. if count==0 else float(beta_dist.ppf(.025, count, len(means)-count+1)),
                      1. if count==len(means) else float(beta_dist.ppf(.975, count+1, len(means)-count))]
                aggregate_ok = (row['replicas']==len(means) and row['covered']==count
                    and row['abstentions']==int((~valid).sum())
                    and abs(row['coverage']-count/len(means))<1e-15
                    and abs(row['target']-exact_energy())<5e-14
                    and abs(row['bias']-bias)<5e-14
                    and np.allclose(row['coverage_ci95'], ci, rtol=1e-13, atol=1e-15)
                    and math.isclose(row['estimated_to_empirical_variance_ratio'], ratio, rel_tol=1e-12, abs_tol=1e-14))
                check('independent aggregate CV BM512: '+folder+'/'+condition, aggregate_ok)
                errors = []
                for i in range(len(means)):
                    e = endpoint.get((condition, i))
                    if (e is None or (e['covered']=='True')!=bool(covered[i])
                        or not math.isclose(float(e['estimate']), float(means[i]), rel_tol=1e-13, abs_tol=1e-14)
                        or not math.isclose(float(e['estimated_variance']), float(mcvar[i]), rel_tol=1e-12, abs_tol=1e-16)):
                        errors.append(i)
                check('independent per-replica CV BM512: '+folder+'/'+condition, not errors, mismatched_replicas=errors)
                coverage_rows.append({'folder':folder, 'condition':condition, 'replicas':len(means),
                    'covered':count, 'coverage':count/len(means), 'coverage_CP95':ci,
                    'bias':bias, 'estimated_to_empirical_variance_ratio':ratio, 'critical_t15':critical})

        # Use bytes, rather than float equality, to catch even signed-zero changes.
        with np.load(ROOT/'fault-study-v1/healthy.npz', allow_pickle=False) as data:
            healthy_raw = np.ascontiguousarray(data['raw'][:, :2048]).view(np.uint8)
            healthy_sources = np.ascontiguousarray(data['sources'][:, :2048]).view(np.uint8)
        for condition in ['acceptance_bias', 'channel_coupling']:
            with np.load(ROOT/'fault-study-v1'/(condition+'.npz'), allow_pickle=False) as data:
                raw = np.ascontiguousarray(data['raw'][:, :2048]).view(np.uint8)
                sources = np.ascontiguousarray(data['sources'][:, :2048]).view(np.uint8)
                raw_same, source_same = np.array_equal(raw, healthy_raw), np.array_equal(sources, healthy_sources)
                row = {'condition':condition, 'replicas':512, 'prefix_sweeps':2048,
                       'raw_byte_equal':bool(raw_same), 'source_byte_equal':bool(source_same),
                       'raw_prefix_sha256':hashlib.sha256(raw).hexdigest(),
                       'source_prefix_sha256':hashlib.sha256(sources).hexdigest()}
                prefix_rows.append(row)
                check('all-replica prefault raw/source bytes identical: '+condition, raw_same and source_same, **row)
    except Exception as exc:
        check('validation completed without exception', False, exception_type=type(exc).__name__, message=str(exc))

    result = {'validator':'independent read-only Phase II artifact audit v1',
        'timestamp_utc':datetime.now(timezone.utc).isoformat(), 'validator_sha256':digest(__file__),
        'imports_study_module':False, 'source_and_data_mutations':False,
        'passed':bool(checks) and all(row['passed'] for row in checks), 'checks':checks,
        'condition_manifests':manifests, 'independently_recomputed_primary_coverage':coverage_rows,
        'full_replica_prefault_prefix_checks':prefix_rows, 'elapsed_seconds':time.monotonic()-started,
        'scope':'Artifact integrity, saved-result arithmetic and namespace/prefix checks; not proof of RNG independence or universal statistical validity'}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    log.write_text('\n'.join(('PASS ' if c['passed'] else 'FAIL ')+c['check'] for c in checks)
                   +f"\nOverall passed: {result['passed']}\nElapsed seconds: {result['elapsed_seconds']:.3f}\n")
    print(json.dumps({'passed':result['passed'], 'checks':len(checks), 'output':str(output),
                      'elapsed_seconds':result['elapsed_seconds']}, indent=2))
    return 0 if result['passed'] else 1


if __name__=='__main__':
    sys.exit(main())
