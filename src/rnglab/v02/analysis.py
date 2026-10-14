"""Derived exploratory diagnostics; primary preregistered outputs are preserved."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from rnglab.v02.runner import csv_write,engine_hash
from rnglab.v02.statistics import proportion_interval,sustained_events


def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);a=p.parse_args()
    root=a.results;manifest=json.loads((root/'manifest.json').read_text());config=manifest['identity']['config']
    if manifest['identity']['engine_sha256']!=engine_hash():raise ValueError('analysis engine differs from saved run')
    rows=[];alarms=[];source=[]
    for file in sorted((root/'trials').glob('*.json')):
        obj=json.loads(file.read_text())['payload'];r=obj['result'];onset=r['fault_onset_step']
        eligible=[w for w in obj['windows'] if w['absolute_error'] is not None and w['retained']>=config['failure_min_samples']]
        observed=[m for m in obj['monitor'] if onset is None or m['time']>onset]
        monobit=min((m['time'] for m in observed if m['p_frequency']<=m['cutoff']),default=None)
        source.append({'scenario':r['scenario'],'workload':r['workload'],'trial':r['trial'],
                       'output_alarm':r['first_alarm'] is not None,'source_alarm':r['source_first_alarm'] is not None,
                       'seed_material_bits':r['entropy_bits_consumed'],'reseeds':r['reseed_count']})
        alarms.append({'scenario':r['scenario'],'workload':r['workload'],'trial':r['trial'],'alarm':r['first_alarm'] is not None})
        for outcome in obj['outcomes']:
            post=[w for w in eligible if onset is None or w['time']>onset]
            failure,_=sustained_events([w['absolute_error'] for w in post],[w['time'] for w in post],outcome['threshold'],config['sustained_windows'])
            rows.append({'scenario':r['scenario'],'workload':r['workload'],'trial':r['trial'],'threshold':outcome['threshold'],
                         'fault_onset':onset,'primary_failure':outcome['failure_time'],'primary_detection':outcome['detection_time'],
                         'primary_failure_before_onset':bool(onset is not None and outcome['failure_time'] is not None and outcome['failure_time']<=onset),
                         'post_onset_failure_confirmation':failure,'monobit_component_detection':monobit,
                         'post_onset_failure_before_composite_detection':bool(failure is not None and (outcome['detection_time'] is None or failure<outcome['detection_time'])),
                         'post_onset_failure_before_monobit_component':bool(failure is not None and (monobit is None or failure<monobit)),
                         'endpoint_error':abs(r['error']),'endpoint_above_tolerance':abs(r['error'])>outcome['threshold']})
    csv_write(root/'exploratory_event_diagnostics.csv',rows);csv_write(root/'source_alarm_trials.csv',source)
    groups=defaultdict(list)
    for row in rows:groups[(row['scenario'],row['workload'],row['threshold'])].append(row)
    summary=[]
    for (scenario,workload,threshold),group in groups.items():
        summary.append({'scenario':scenario,'workload':workload,'threshold':threshold,'trials':len(group),
                        'pre_onset_primary_failures':sum(r['primary_failure_before_onset'] for r in group),
                        'post_onset_failures':sum(r['post_onset_failure_confirmation'] is not None for r in group),
                        'post_failure_before_composite':sum(r['post_onset_failure_before_composite_detection'] for r in group),
                        'post_failure_before_monobit':sum(r['post_onset_failure_before_monobit_component'] for r in group),
                        'monobit_detected':sum(r['monobit_component_detection'] is not None for r in group),
                        'endpoint_above_tolerance':sum(r['endpoint_above_tolerance'] for r in group)})
    csv_write(root/'exploratory_event_summary.csv',summary)
    healthy={'synthetic','prng','csprng','hybrid_healthy','hybrid_frozen','os_recorded','drand_seeded'}
    clean=[r for r in alarms if r['scenario'] in healthy]
    k=sum(r['alarm'] for r in clean);rate,lo,hi=proportion_interval(k,len(clean))
    grouped=[]
    for name in sorted({r['scenario'] for r in alarms}):
        for w in sorted({r['workload'] for r in alarms}):
            group=[r for r in alarms if r['scenario']==name and r['workload']==w]
            if group:
                n=sum(r['alarm'] for r in group);v,l,h=proportion_interval(n,len(group))
                grouped.append({'scenario':name,'workload':w,'alarms':n,'trials':len(group),'fraction':v,'ci_low':l,'ci_high':h})
    csv_write(root/'alarm_rates.csv',grouped)
    (root/'analysis_notes.json').write_text(json.dumps({'analysis_status':'exploratory diagnostics, computed after viewing primary outcomes; does not redefine primary failure',
                                                      'healthy_output_alarm_pool':{'alarms':k,'trials':len(clean),'fraction':rate,'ci95_low':lo,'ci95_high':hi,'note':'descriptive pooled interval; OS/beacon reuse across workloads and configuration pairing introduce cross-group dependence'},
                                                      'monobit_rule':'p_frequency compared with the prespecified composite per-test cutoff; a component diagnostic, not a separately optimized sequential detector',
                                                      'post_onset_rule':'reset three-window confirmation after onset; sensitivity analysis to preexisting application excursions'},indent=2))
    print('Wrote exploratory timing and source-alarm diagnostics')

if __name__=='__main__':main()
