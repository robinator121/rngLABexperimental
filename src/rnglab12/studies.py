"""Predeclared bounded studies, immutable inputs, no automatic deployment."""
from pathlib import Path
from dataclasses import asdict
from collections import Counter
import json, math, time
import numpy as np
from scipy.stats import norm, t, beta
from .execution import ROOT,APPLICATIONS,DESIGN,CONFIG_ID,now,digest,write_new,run_tasks,source_identity
from .calibration import Trajectory,FailurePolicy,calibrate,validate_holdout
from .seeds import derive_seed

def freeze(name,body):
    path=ROOT/'protocols'/(name+'.json')
    if path.exists():
        prior=json.loads(path.read_text())
        if prior['body']!=body:raise ValueError('Frozen protocol differs; create a new version')
        return prior
    record={'frozen_at_utc':now(),'body':body,'body_sha256':digest(body),'source_identity':source_identity()}
    write_new(path,record);return record

def trajectory(row,kind='application',scales=None):
    errors=row['errors'];reference=row['reference'];observable=row['observable']
    if kind=='source':errors=row['source_scores'];reference=0.;observable='delivered_MSB_APT_style_score'
    if kind=='combined':
        errors=np.maximum(np.array(errors)/scales[0],np.array(row['source_scores'])/scales[1]).tolist()
        reference=0.;observable='joint_development_scaled_score'
    return Trajectory(errors,row['times'],row['application'],observable,row['config_id'],
                      row['lineage_id'],row['task']['split'],reference,row['seed'])

def calibration_study():
    protocol=freeze('healthy-calibration-v1',{
        'design':DESIGN,'applications':list(APPLICATIONS),'development_groups_per_app':64,
        'calibration_groups_per_app':511,'final_holdout_groups_per_app':512,
        'threshold_alpha_per_application':.01,'calibration_family_size':1,
        'validation':{'target_per_application_rate':.05,'confidence':.95,'confidence_family_size':6},
        'policies':['application','source','combined'],'primary_gate':'all six application policy CP upper bounds <= .05',
        'source_scope':'delivered uniform MSB rolling proportion; APT-style not NIST90B assessment',
        'uncertainty':'Geyer IMS primary comparison; batch32/128 retained; empirical coverage diagnostic, not selected on holdout',
        'no_holdout_tuning':True,'max_cpu_workers':4,'max_trajectories':6522})
    tasks=[{'application':a,'split':split,'group':g} for a in APPLICATIONS
           for split,n in [('development',64),('calibration',511)] for g in range(n)]
    rows=run_tasks(ROOT/'runs/healthy-calibration',tasks,classification='DEVELOPMENT_AND_CALIBRATION',protocol=protocol)
    output={}
    for app in APPLICATIONS:
        dev=[r for r in rows if r['application']==app and r['task']['split']=='development']
        cal=[r for r in rows if r['application']==app and r['task']['split']=='calibration']
        scales=[float(np.median([trajectory(r).score() for r in dev])),
                float(np.median([trajectory(r,'source').score() for r in dev]))]
        if min(scales)<=0:raise ValueError('Development scaling degenerate')
        output[app]={'scales_from_development':scales,'policies':{
            kind:calibrate([trajectory(r,kind,scales) for r in cal],alpha=.01).to_dict()
            for kind in ['application','source','combined']}}
    path=ROOT/'evidence/healthy-policies.json'
    if path.exists():
        if json.loads(path.read_text())!=output:raise ValueError('Policy recomputation mismatch')
    else:write_new(path,output)
    return output

def healthy_validation():
    policies=json.loads((ROOT/'evidence/healthy-policies.json').read_text())
    protocol=freeze('healthy-final-holdout-v1',{
        'calibration_protocol':digest(json.loads((ROOT/'protocols/healthy-calibration-v1.json').read_text())),
        'policies_sha256':digest(policies),'groups_per_application':512,'split':'final_holdout',
        'primary_gate':'six per-application one-sided CP upper bounds <= .05 at simultaneous confidence .95',
        'source_combined_bounds':'pointwise descriptive secondary results; not multiplicity-confirmed superiority',
        'stop_rule':'fixed sample count, failures retained, no policy refit'})
    rows=run_tasks(ROOT/'runs/healthy-final-holdout',[
        {'application':a,'split':'final_holdout','group':g} for a in APPLICATIONS for g in range(512)],
        classification='FROZEN_HEALTHY_VALIDATION',protocol=protocol)
    results={};coverage=[]
    for app in APPLICATIONS:
        selected=[r for r in rows if r['application']==app];p=policies[app];results[app]={}
        for kind in ['application','source','combined']:
            rule=FailurePolicy.from_dict(p['policies'][kind]);family=6 if kind=='application' else 1
            result=validate_holdout(rule,[trajectory(r,kind,p['scales_from_development']) for r in selected],family_size=family,target_rate=.05)
            old=max(.05,.1*abs(selected[0]['reference']))
            if kind=='application':result['historical_tolerance_crossings']=sum(trajectory(r).score()>old for r in selected)
            results[app][kind]=result
        for method in ['geyer','batch32','batch128']:
            summaries=[r['uncertainty']['geyer'] if method=='geyer' else r['uncertainty']['batch_means'][method[5:]] for r in selected]
            usable=[s for s in summaries if s['covered'] is not None];covered=sum(s['covered'] for s in usable)
            coverage.append({'application':app,'method':method,'groups':len(selected),'usable':len(usable),'covered':covered,
                             'coverage':covered/len(usable) if usable else None,'scope':'approximate finite-run interval diagnostic; no holdout method selection'})
    report={'by_application':results,'all_primary_pass':all(p['application']['passes'] for p in results.values()),
            'coverage':coverage,'at':now(),'claims':'per-app finite-horizon MT19937 healthy-law rates; not shifted/physical-source transfer'}
    write_new(ROOT/'tables/healthy-validation.json',report)
    return report

def fault_study():
    gate=json.loads((ROOT/'tables/healthy-validation.json').read_text())
    if not gate['all_primary_pass']:raise PermissionError('Healthy calibration gate failed; do not continue detector confirmation')
    policies=json.loads((ROOT/'evidence/healthy-policies.json').read_text())
    protocol=freeze('transient-faults-v1',{
        'design':DESIGN,'policy_sha256':digest(policies),'split':'confirmation','new_seed_domain':'fault-v012',
        'groups_per_app':64,'apps':list(APPLICATIONS),'onset_zero_based_row':4096,'duration_rows':512,
        'arms':['healthy','joint_zero','proposal_zero_or_selection_zero','acceptance_zero','joint_minpositive','joint_reset'],
        'primary_scope':'bounded mechanism/recovery replication; descriptive paired effects, not powered detector-superiority test',
        'source_comparison':{'applications':['integration','vmc'],'groups':32,'families':['synthetic_fault','deterministic_seeded_or_reseeded_csprng'],
            'arms':['healthy','output_zero','reseed_input_zero'],'source_population_calibration_transfer':'NOT_ASSUMED'},
        'endpoints':'paired retained-mean error, observed threshold crossing, pre-onset events, rolling/cumulative recovery and censoring',
        'no_nist_superiority_claim':True})
    tasks=[]
    for app in APPLICATIONS:
        for g in range(64):
            base={'application':app,'split':'confirmation','group':g,'seed_domain':'fault-v012'}
            tasks.append({**base,'arm':'healthy'})
            joint={'mode':'zero_burst','scope':'all','start_row':4096,'duration_rows':512}
            tasks.append({**base,'arm':'joint_zero','fault':joint})
            if app not in ['integration','pi']:
                proposal='selection' if app=='ising' else 'proposal'
                tasks.extend([{**base,'arm':'proposal_zero_or_selection_zero','fault':{**joint,'scope':proposal}},
                              {**base,'arm':'acceptance_zero','fault':{**joint,'scope':'acceptance'}},
                              {**base,'arm':'joint_reset','fault':joint,'reset_at':4608}])
            if app in ['mcmc','vmc','harmonic_pimc']:
                tasks.append({**base,'arm':'joint_minpositive','fault':joint,'positive_acceptance_control':True})
    rows=run_tasks(ROOT/'runs/transient-faults',tasks,classification='FROZEN_MECHANISM_REPLICATION',protocol=protocol)
    summary=[]
    for app in APPLICATIONS:
        p=policies[app];arms=sorted({r['task']['arm'] for r in rows if r['application']==app})
        healthy={r['task']['group']:r for r in rows if r['application']==app and r['task']['arm']=='healthy'}
        for arm in arms:
            selected=[r for r in rows if r['application']==app and r['task']['arm']==arm]
            events={kind:[FailurePolicy.from_dict(p['policies'][kind]).event(trajectory(r,kind,p['scales_from_development'])) for r in selected] for kind in ['application','source','combined']}
            rule=FailurePolicy.from_dict(p['policies']['application']);recover=[];cumulative=[]
            for r,e in zip(selected,events['application']):
                if not e['event'] or arm=='healthy':continue
                j=next((j for j in range(1,len(r['times'])) if r['times'][j-1]>=4608 and r['times'][j-1]>=e['first_event_time'] and max(r['errors'][j-1:j+1])<=rule.threshold),None)
                recover.append(None if j is None else r['times'][j]-4608)
                cumulative.append(abs(r['cumulative_estimates'][-1]-r['reference'])<=rule.threshold)
            errors=np.array([abs(r['uncertainty']['mean']-r['reference']) for r in selected])
            delta=errors-np.array([abs(healthy[r['task']['group']]['uncertainty']['mean']-r['reference']) for r in selected])
            summary.append({'application':app,'arm':arm,'groups':len(selected),'median_full_error':float(np.median(errors)),
                'mean_paired_error_difference':float(delta.mean()),'median_paired_error_difference':float(np.median(delta)),
                'application_events':sum(e['event'] for e in events['application']),
                'pre_onset_application_events':sum(e['event'] and e['first_event_time']<=4096 for e in events['application']),
                'source_events':sum(e['event'] for e in events['source']),'combined_events':sum(e['event'] for e in events['combined']),
                'application_event_without_source_event':sum(a['event'] and not s['event'] for a,s in zip(events['application'],events['source'])),
                'rct_reference_events':sum(r['rct_event'] is not None for r in selected),
                'recovery_observed':sum(d is not None for d in recover),'recovery_eligible':len(recover),
                'median_recovery_delay_among_observed':float(np.median([d for d in recover if d is not None])) if any(d is not None for d in recover) else None,
                'final_cumulative_below_tolerance_among_recovery_eligible':sum(cumulative),
                'scope':'descriptive paired replication; health policies separately calibrated at same target, no superiority claim'})
    write_new(ROOT/'tables/transient-faults.json',summary)
    return summary

def periodic_study():
    n=4096;groups=1536;cells=[('integration',16,1024),('integration',16,2048),('pi',8,1024),('pi',8,2048)]
    predictions=[]
    for app,q,m in cells:
        if app=='integration':
            raw=[1/(2*k+1) for k in range(5)];mu=raw[1];var=raw[2]-mu**2
            central4=raw[4]-4*mu*raw[3]+6*mu**2*raw[2]-3*mu**4
        else:
            prob=math.pi/4;var=16*prob*(1-prob);central4=256*prob*(1-prob)*((1-prob)**3+prob**3)
        v=var*((n-m)/n**2+m**2/(q*n**2));sumw4=(n-m)/n**4+q*(m/(q*n))**4
        fourth=3*v*v+(central4-3*var*var)*sumw4
        rse=math.sqrt(max(0,fourth/v**2-1)/groups)
        predictions.append({'application':app,'q':q,'m':m,'n':n,'variance':v,'inflation':v/(var/n),
                            'predicted_variance_relative_se':rse,'simultaneous_normal_approx_halfwidth':float(norm.ppf(1-.05/32)*rse)})
    protocol=freeze('periodic-variance-v1',{'groups_per_cell_per_cohort':groups,'cohorts':['confirmation','final_holdout'],
        'cells':predictions,'generator':'NumPy PCG64 explicitly noncryptographic IID simulation control',
        'paired_healthy_and_fault':'same seed and full native uniforms; first q native observations inside fault segment repeated',
        'fault_start':1024,'hypotheses':'variance/formula ratio within [.8,1.2]; ensemble mean within target±.02 integration/±.05 pi',
        'uncertainty':'independent-group bootstrap percentile, 10000 resamples; 8 variance and 8 mean claims one Bonferroni family16; approximate not exact finite-sample',
        'stop_rule':'1536 groups fixed; no extension or selective pooling; two fresh cohorts; every negative retained',
        'coverage':'fixed32 batch means nominal95; descriptive undercoverage, not repaired by ensemble mean agreement'})
    directory=ROOT/'runs/periodic';directory.mkdir(exist_ok=True)
    records=[];summaries=[];started=time.monotonic()
    for split in ['confirmation','final_holdout']:
        for prediction in predictions:
            app=prediction['application'];q=prediction['q'];m=prediction['m'];reference=1/3 if app=='integration' else math.pi
            cell=[]
            for g in range(groups):
                seed=derive_seed(split,f'periodic-v012/{app}/q{q}/m{m}',g)
                rng=np.random.Generator(np.random.PCG64(seed));u=rng.random((n,1 if app=='integration' else 2))
                y=u[:,0]**2 if app=='integration' else 4*(np.sum(u*u,axis=1)<1).astype(float)
                faulty=y.copy();faulty[1024:1024+m]=np.tile(y[1024:1024+q],m//q)
                batches=faulty.reshape(-1,32).mean(axis=1);mcse=float(batches.std(ddof=1)/math.sqrt(len(batches)))
                estimate=float(faulty.mean());half=float(t.ppf(.975,len(batches)-1))*mcse
                row={'application':app,'split':split,'group':g,'seed':seed,'q':q,'m':m,'n':n,'reference':reference,
                     'healthy_mean':float(y.mean()),'fault_mean':estimate,'batch32_mcse':mcse,
                     'covered':estimate-half<=reference<=estimate+half,'protocol_sha256':protocol['body_sha256'],
                     'source_type':'PCG64 noncryptographic synthetic control'}
                records.append(row);cell.append(row)
            values=np.array([r['fault_mean'] for r in cell]);healthy=np.array([r['healthy_mean'] for r in cell])
            rng=np.random.default_rng(derive_seed(split,'periodic-bootstrap/'+app+str(m),'fixed'))
            boot=[]
            for _ in range(100):
                b=rng.integers(0,groups,size=(100,groups));boot.extend(values[b].var(axis=1,ddof=1).tolist())
            quantiles=np.quantile(boot,[.05/32,1-.05/32])/prediction['variance']
            meanhalf=float(t.ppf(1-.05/32,groups-1))*values.std(ddof=1)/math.sqrt(groups)
            mean_ci=[float(values.mean()-meanhalf),float(values.mean()+meanhalf)];margin=.02 if app=='integration' else .05
            summaries.append({**prediction,'split':split,'groups':groups,'observed_variance':float(values.var(ddof=1)),
                'variance_ratio_to_theory':float(values.var(ddof=1)/prediction['variance']),
                'variance_ratio_simultaneous_approx_interval':quantiles.tolist(),'variance_equivalence_supported':bool(quantiles[0]>.8 and quantiles[1]<1.2),
                'observed_inflation':float(values.var(ddof=1)/healthy.var(ddof=1)),
                'ensemble_mean':float(values.mean()),'mean_simultaneous_approx_interval':mean_ci,
                'mean_equivalence_supported':mean_ci[0]>reference-margin and mean_ci[1]<reference+margin,
                'batch32_coverage':sum(r['covered'] for r in cell)/groups})
    write_new(directory/'records.json',records);write_new(ROOT/'tables/periodic-variance.json',summaries)
    write_new(directory/'receipt.json',{'at':now(),'protocol_sha256':protocol['body_sha256'],'records':len(records),
        'source_identity':source_identity(),'elapsed_seconds':time.monotonic()-started})
    return summaries

def source_study():
    protocol=freeze('source-comparison-v1',{'applications':['integration','vmc'],'groups_per_app':32,'split':'confirmation',
        'seed_domain':'source-v012','design':DESIGN,
        'configurations':['MT19937','HMAC_DRBG_R0','HMAC_DRBG_reseed_every2_fixed4096B_generates'],
        'arms':'healthy; post-uniform joint zero; bad reseed material for periodic reseed configuration',
        'reseed_fault':'zero first144 reseed-input bytes; initial seed untouched; underlying input consumed',
        'independence':'group seeds matched by label; streams differ across generator families; not identical variates',
        'scope':'engineering/descriptive comparison only, no physical measurements, entropy/security proof or MT healthy-calibration transfer'})
    tasks=[]
    for app in ['integration','vmc']:
        for g in range(32):
            for configuration in ['MT19937','HMAC_DRBG_R0','HMAC_DRBG_R2']:
                base={'application':app,'split':'confirmation','group':g,'seed_domain':'source-v012','configuration':configuration}
                if configuration!='MT19937':base['source_spec']={'family':'deterministic_seeded_or_reseeded_csprng','request_bytes':4096,**({'reseed_interval_generate_calls':2} if configuration.endswith('R2') else {})}
                tasks.append({**base,'arm':'healthy'})
                tasks.append({**base,'arm':'output_zero','fault':{'mode':'zero_burst','scope':'all','start_row':4096,'duration_rows':512}})
                if configuration.endswith('R2'):
                    tasks.append({**base,'arm':'reseed_input_zero','reseed_material_fault':{'mode':'zero_burst','placement':'reseed_input_bytes','start_byte':0,'duration_bytes':144}})
    rows=run_tasks(ROOT/'runs/source-comparison',tasks,classification='SOURCE_ENGINEERING_PILOT',protocol=protocol)
    summaries=[]
    for app in ['integration','vmc']:
        for configuration in ['MT19937','HMAC_DRBG_R0','HMAC_DRBG_R2']:
            for arm in ['healthy','output_zero','reseed_input_zero']:
                selected=[r for r in rows if r['application']==app and r['task']['configuration']==configuration and r['task']['arm']==arm]
                if not selected:continue
                means=np.array([r['uncertainty']['mean'] for r in selected]);ref=selected[0]['reference']
                summaries.append({'application':app,'configuration':configuration,'arm':arm,'groups':len(selected),
                    'mean_bias':float(means.mean()-ref),'variance_of_estimates':float(means.var(ddof=1)),
                    'median_absolute_error':float(np.median(abs(means-ref))),
                    'source_family':selected[0]['source_family'],'physical_samples_used':False,
                    'scope':'descriptive pilot; no calibrated monitor-transfer, reseed superiority, entropy or security claim'})
    write_new(ROOT/'tables/source-comparison.json',summaries)
    return summaries

def rus_study():
    from .reliability import Observation,train_predictor,validate,predict
    gate=json.loads((ROOT/'tables/healthy-validation.json').read_text())
    if not gate['by_application']['integration']['application']['passes']:raise PermissionError('Integration healthy gate failed')
    policies=json.loads((ROOT/'evidence/healthy-policies.json').read_text())
    policy=FailurePolicy.from_dict(policies['integration']['policies']['application'])
    config=digest({'design':DESIGN,'policy':policy.to_dict(),'fault':'uniform scale(t), start2048,end8192, final severity U(.04,.35)','issued':3072})
    protocol=freeze('prospective-rus-v1',{'population_config':config,'application':'integration','fault':'linear_scale',
        'seed_domain':'rus-v012','groups':{'development':128,'calibration':128,'final_holdout':256},
        'design':DESIGN,'fault_start':2048,'fault_end':8192,'severity_population':[.04,.35],'issued_at':3072,
        'features':'past512 uniform mean, mean change from preceding512, current rolling absolute error; no future fault parameters',
        'eligibility':'still valid at issuance; already-failed cases counted separately, never future-selected',
        'models':['application_aware','source_only','constant_conformal'],'alpha':.1,
        'censoring':'right censor retained with conservative conformal residual; unresolved validation UNKNOWN',
        'primary_gate':'application-aware worst-case pointwise95 unsafe upper<=.15 AND positive usable fraction>=.5',
        'baselines':'source-only and intercept-only conformal; fixed1024 lower bound and no-prediction also reported',
        'scope':'known stochastic defect population only; unknown types abstain; no conditional-risk guarantee or universal prediction'})
    def tasks_for(split,n):
        tasks=[]
        for g in range(n):
            rng=np.random.default_rng(derive_seed(split,'rus-v012/fault',g))
            tasks.append({'application':'integration','split':split,'group':g,'seed_domain':'rus-v012','issued_at':3072,
                          'linear_scale':{'start_row':2048,'duration_rows':6144,'final_severity':float(rng.uniform(.04,.35))}})
        return tasks
    def observations(rows,mode):
        kept=[];already=[]
        for row in rows:
            event=policy.event(trajectory(row));at=event['first_event_time'];issued=3072
            if at is not None and at<=issued:already.append(row['lineage_id']);continue
            features=row['issuance']['features']
            if mode=='source_only':features=features[:2]
            if mode=='constant_conformal':features=[1.]
            kept.append(Observation(row['lineage_id'],row['task']['split'],'integration','linear_scale',config,issued,tuple(features),issued,
                                    event_at=at,censored_at=8192 if at is None else None))
        return kept,already
    training=run_tasks(ROOT/'runs/rus-development-calibration',tasks_for('development',128)+tasks_for('calibration',128),classification='RUS_TRAIN_CALIBRATION',protocol=protocol)
    models={};pop={}
    for mode in ['application_aware','source_only','constant_conformal']:
        dev,ade=observations([r for r in training if r['task']['split']=='development'],mode)
        cal,aca=observations([r for r in training if r['task']['split']=='calibration'],mode)
        model=train_predictor(dev,cal,application='integration',fault='linear_scale',configuration=config)
        model.save(ROOT/'evidence'/('rus-model-'+mode+'.json'));models[mode]=model
        pop[mode]={'development_eligible':len(dev),'calibration_eligible':len(cal),'development_already_failed':len(ade),'calibration_already_failed':len(aca)}
    validation_protocol=freeze('rus-final-holdout-v1',{'parent_protocol_sha256':protocol['body_sha256'],'models':{m:p.model_id for m,p in models.items()},'groups':256,'no_refit':True})
    final=run_tasks(ROOT/'runs/rus-final-holdout',tasks_for('final_holdout',256),classification='FROZEN_PROSPECTIVE_RUS_VALIDATION',protocol=validation_protocol)
    reports={}
    for mode,model in models.items():
        obs,already=observations(final,mode);report=validate(model,obs)
        report['population']={**pop[mode],'final_groups_total':len(final),'final_eligible':len(obs),'already_failed_at_issuance':len(already)}
        report['model_id']=model.model_id
        report['event_only_absolute_error']=None
        errs=[abs(predict(model,r).lower_bound-r.remaining_observed) for r in obs if r.event_at is not None and predict(model,r).lower_bound is not None]
        if errs:report['event_only_absolute_error']={'n':len(errs),'median':float(np.median(errs)),'scope':'event-only descriptive, not censor-corrected accuracy'}
        reports[mode]=report
    obs,_=observations(final,'application_aware');fixed={'SAFE':0,'UNSAFE':0,'UNKNOWN':0}
    for row in obs:
        outcome='UNSAFE' if row.event_at is not None and row.remaining_observed<1024 else 'SAFE' if row.remaining_observed>=1024 else 'UNKNOWN';fixed[outcome]+=1
    reports['fixed1024']={'counts':fixed,'eligible':len(obs)}
    reports['no_prediction']={'abstentions':len(obs),'eligible':len(obs),'usable_fraction':0}
    primary=reports['application_aware'];envelope=primary['unsafe_risk_95pct_identification_and_sampling_envelope']
    usable=primary['positive_bound_count']/primary['independent_seed_groups'] if primary['independent_seed_groups'] else 0.
    reports['primary_gate']={'supported':bool(envelope is not None and envelope[1]<=.15 and usable>=.5),
        'unsafe_worst_case_upper':None if envelope is None else envelope[1],'positive_usable_fraction':usable,
        'criterion':'worst-case95pointwise unsafe upper<=.15 and positive usable fraction>=.5; known population only'}
    reports['scope']='Known population, discrete-grid event label, marginal construction. No superiority or conditional safety conclusion.'
    write_new(ROOT/'tables/rus-validation.json',reports)
    return reports
