"""Gated offline random-channel fault comparison after locked healthy validation."""
from __future__ import annotations
import argparse
import csv
import json
import math
import time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
from scipy.stats import t
import study

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'fault-study-v1'


def freeze():
    gate=json.loads((ROOT/'confirmation-v1/primary_decision.json').read_text())
    if not gate['fault_study_allowed']:raise RuntimeError('Frozen healthy baseline gate did not pass; no fault study permitted')
    path=OUT/'protocol.json'
    if path.exists():raise FileExistsError('Fault protocol already frozen')
    p={'experiment_id':'fault-study-v1','timestamp_utc':datetime.now(timezone.utc).isoformat(),
       'code_sha256':study.digest(__file__),'study_code_sha256':study.digest(ROOT/'study.py'),
       'healthy_gate_sha256':study.digest(ROOT/'confirmation-v1/primary_decision.json'),
       'classification':'frozen secondary fault experiment; healthy acceptance gate passed; no calibration/test streams generated yet',
       'model':{'M':16,'beta':2.,'delta':.4,'initial':'exact Gaussian target','kernel':'local Metropolis plus exact centroid heat bath every sweep'},
       'calibration_replicas':399,'test_replicas':512,'total_sweeps':10240,'fault_start_sweep_1based':2049,
       'postfault_analysis_sweeps_1based':[2049,10240],'primary_CI':'CV BM512 on8192 postonset sweeps; intended finiteM target',
       'faults':{'acceptance_bias':{'code':1,'severity':.08,'mapping':'u_accept -> min(u_accept+.08,nextafter(1,0))'},
                 'channel_coupling':{'code':2,'mapping':'u_accept=u_proposal (both margins uniform; joint law altered)'}},
       'fault_scope':'offline experimental adapter; source role identified; no physical data or live instruments used',
       'pairing':'healthy/fault conditions reuse initial/update/centroid/scan seeds within a replica; independent replica indices and disjoint cal/test/healthy-confirm namespaces',
       'monitors':{'window_sweeps':1024,'stride_sweeps':256,'nominal_whole_run_false_alarm':.05,
          'source_marginal':'max absolute deviation of proposal or acceptance window mean from .5',
          'source_cross_channel':'absolute window mean of (up-.5)*(ua-.5)',
          'application_shape':'absolute window mean shape r_g² deviation from finiteM target',
          'threshold':'ceil((calibration_n+1)*.95)-th order statistic of entire-path max; strict > threshold; no retuning',
          'healthy_validation_gate':'each monitor CP95 upper false-alarm bound<=.10 on512 healthy test paths; report actual rates, not exact equal FPR',
          'postonset_windows':'classified by endpoint>=2049; earliest endpoint2304 spans1281-2304 and intentionally includes prefault data',
          'multiple_monitors':'5% target is separate for each monitor, not a combined5% policy',
          'paired_comparison':'whole-run alert indicator fault minus healthy with replica bootstrap95; secondary descriptive'},
       'application_failure':{'observable':'CV energy','rolling_window':1024,'stride':256,'absolute_error_tolerance':.05,'persistence_windows':2,
          'meaning':'prespecified numerical tolerance excursion, not a physical safety limit; distinct from calibrated anomaly monitor',
          'censoring':'no event by10240 -> right censored; already failed before onset counted separately'},
       'primary_effects':['paired postonset CV mean difference with t CI across replicas','paired intended-target coverage change with replica bootstrap95'],
       'secondary':['TD/shape mean and uncertainty','source marginal/cross-channel summaries','whole-run alarm rates and delays','failure occurrence and alert-vs-failure delays'],
       'no_superiority_rule':'do not generalize application superiority; marginal-source and cross-channel monitors have different information; use paired comparisons and actual healthy false alarms',
       'seed_stages':['fault-calibration-v1','fault-test-v1']}
    study.save(path,p);print('Frozen fault protocol',path,study.digest(path))


def windows(raw,sources,window=1024,stride=256):
    obs=study.observables(raw);ends=np.arange(window,raw.shape[1]+1,stride)
    def rolling(a):
        cs=np.concatenate([np.zeros((len(a),1)),np.cumsum(a,axis=1)],axis=1)
        return (cs[:,ends]-cs[:,ends-window])/window
    up=rolling(sources[:,:,0]);ua=rolling(sources[:,:,1]);cross=rolling(sources[:,:,2]);shape=rolling(obs['shape']);cv=rolling(obs['CV'])
    scores={'source_marginal':np.maximum(abs(up-.5),abs(ua-.5)),
            'source_cross_channel':abs(cross),
            'application_shape':abs(shape-study.reference()['means']['shape'])}
    errors=abs(cv-study.reference()['finite_energy'])
    failure_hit=(errors[:,:-1]>.05)&(errors[:,1:]>.05)
    failures=np.full(len(raw),-1,dtype=int)
    for j,end in enumerate(ends[1:]):failures=np.where((failures<0)&failure_hit[:,j],end,failures)
    return ends,scores,failures


def first_hit(scores,ends,threshold,start=0):
    hits=(scores>threshold)&(ends[None,:]>=start)
    anyhit=hits.any(1);index=np.argmax(hits,axis=1)
    return np.where(anyhit,ends[index],-1)


def run():
    p=json.loads((OUT/'protocol.json').read_text())
    if p['code_sha256']!=study.digest(__file__) or p['study_code_sha256']!=study.digest(ROOT/'study.py'):raise RuntimeError('Code changed since freeze')
    gate_path=ROOT/'confirmation-v1/primary_decision.json'
    if study.digest(gate_path)!=p['healthy_gate_sha256'] or not json.loads(gate_path.read_text())['fault_study_allowed']:
        raise RuntimeError('Healthy gate changed or failed after fault protocol freeze')
    started=time.monotonic()
    base={'initial':'stationary','steps':10240,'delta':.4,'refresh':True,'source_record':True,'onset':2048}
    cal={**base,'id':'healthy_calibration','replicas':399}
    cal_path=study.run_condition(OUT,'fault-calibration-v1',cal)
    a=np.load(cal_path);ends,cal_scores,_=windows(a['raw'],a['sources']);rank=math.ceil((399+1)*.95)
    thresholds={k:float(np.sort(v.max(1))[rank-1]) for k,v in cal_scores.items()}
    study.save(OUT/'calibration.json',{'rank':rank,'n':399,'thresholds':thresholds,'full_path_maxima':{k:v.max(1).tolist() for k,v in cal_scores.items()},'test_not_generated_yet':True})
    del a
    configs=[{**base,'id':'healthy','replicas':512},
             {**base,'id':'acceptance_bias','replicas':512,'fault':1,'severity':.08},
             {**base,'id':'channel_coupling','replicas':512,'fault':2}]
    all_rows=[];rep_rows=[];monitor_rows=[];timing_rows=[];effects={};paired_alarm_effects=[];healthy_means=None;healthy_covered=None;healthymonitor={}
    for config in configs:
        path=study.run_condition(OUT,'fault-test-v1',config)
        data=np.load(path);raw=data['raw'];sources=data['sources'];obs=study.observables(raw)
        rows,rr,decomp,acf=study.analyze_window(path,2048,8192,'postonset',True)
        all_rows+=rows;rep_rows+=rr
        ends,scores,failure=windows(raw,sources)
        for name,score in scores.items():
            alert=first_hit(score,ends,thresholds[name]);post=first_hit(score,ends,thresholds[name],2049)
            k=int((alert>=0).sum());pre=int(((alert>=0)&(alert<=2048)).sum());postk=int((post>=0).sum())
            row={'condition':config['id'],'monitor':name,'events':k,'replicas':512,'whole_run_rate':k/512,
                 'whole_run_CI95':study.binomial_interval(k,512),'preonset_alerts':pre,'postonset_alerts':postk,
                 'threshold':thresholds[name],'median_postonset_delay':None if postk==0 else float(np.median(post[post>=0]-2049)),
                 'healthy_gate_pass':bool(study.binomial_interval(k,512)[1]<=.10) if config['id']=='healthy' else None}
            monitor_rows.append(row)
            if config['id']=='healthy':healthymonitor[name]=(alert>=0)
            else:
                delta=(alert>=0).astype(int)-healthymonitor[name].astype(int)
                boot=study.generator('analysis',f'alarm/{config["id"]}/{name}',0).integers(512,size=(10000,512))
                paired_alarm_effects.append({'condition':config['id'],'monitor':name,'paired_whole_run_alert_difference':float(delta.mean()),
                    'bootstrap_CI95':np.quantile(delta[boot].mean(1),[.025,.975]).tolist(),
                    'fault_only':int((delta==1).sum()),'healthy_only':int((delta==-1).sum()),'scope':'secondary; actual healthy false alarm rates differ across monitors'})
            for i in range(512):
                eligible=failure[i]>2048
                timing_rows.append({'condition':config['id'],'monitor':name,'replica':i,'first_alert_sweep':int(alert[i]),
                    'first_postonset_alert_sweep':int(post[i]),'failure_sweep':int(failure[i]),'failure_right_censored':bool(failure[i]<0),
                    'failed_before_onset':bool(0<=failure[i]<=2048),'postonset_failure':bool(eligible),
                    'alert_minus_failure':int(post[i]-failure[i]) if eligible and post[i]>=0 else None})
        cv=obs['CV'][:,2048:];means=cv.mean(1);bm=study.variance_methods(cv,(512,))['BM512'];covered=abs(means-study.reference()['finite_energy'])<=bm[1]*np.sqrt(bm[0])
        source_summary={'mean_proposal':float(sources[:,2048:,0].mean()),'mean_acceptance':float(sources[:,2048:,1].mean()),'mean_cross_product':float(sources[:,2048:,2].mean())}
        failed=int((failure>=0).sum());prefailed=int(((failure>=0)&(failure<=2048)).sum())
        failure_summary={'events':failed,'total':512,'before_onset':prefailed,'new_after_onset':failed-prefailed,'right_censored':512-failed,'CI95':study.binomial_interval(failed,512),
            'median_new_failure_sweep':None if failed==prefailed else float(np.median(failure[failure>2048]))}
        if config['id']=='healthy':healthy_means=means.copy();healthy_covered=covered.copy()
        else:
            difference=means-healthy_means;se=difference.std(ddof=1)/math.sqrt(512);q=float(t.ppf(.975,511))
            dc=covered.astype(int)-healthy_covered.astype(int);boot=study.generator('analysis',f'fault/{config["id"]}',0).integers(512,size=(10000,512))
            effects[config['id']]={'paired_mean_energy_change':float(difference.mean()),'paired_mean_change_CI95':[float(difference.mean()-q*se),float(difference.mean()+q*se)],
                'paired_coverage_change':float(dc.mean()),'paired_coverage_change_bootstrap_CI95':np.quantile(dc[boot].mean(1),[.025,.975]).tolist()}
        effects.setdefault(config['id'],{}).update({'source_summary':source_summary,'failure_summary':failure_summary})
        del data,raw,sources,obs
    study.save(OUT/'analysis.json',all_rows);study.save(OUT/'monitor_results.json',monitor_rows);study.save(OUT/'effects.json',effects)
    study.save(OUT/'paired_alarm_effects.json',paired_alarm_effects)
    for name,rows in [('replicate_endpoints.csv',rep_rows),('alarm_failure_times.csv',timing_rows)]:
        keys=list(dict.fromkeys(k for row in rows for k in row))
        with (OUT/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader();writer.writerows(rows)
    decision={'healthy_monitor_gates_pass':all(r['healthy_gate_pass'] for r in monitor_rows if r['condition']=='healthy'),
              'elapsed_seconds':time.monotonic()-started,'thresholds_frozen_before_test':True,'no_threshold_retuning':True,
              'comparison_scope':'same nominal5% whole-run calibration; actual healthy test rates reported; no claim of all-source diagnostic superiority'}
    study.save(OUT/'decision.json',decision);print(json.dumps({'decision':decision,'effects':effects,'monitors':monitor_rows},indent=2))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['freeze','run']);a=parser.parse_args()
    freeze() if a.command=='freeze' else run()


if __name__=='__main__':main()
