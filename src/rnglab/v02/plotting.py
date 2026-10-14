"""Rebuild publication-style figures from saved trial outputs; no simulations rerun."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import SymmetricalLogLocator,LogLocator,LogFormatterMathtext,NullFormatter,FuncFormatter
import numpy as np
from rnglab.v02.statistics import mean_interval


def read(path):
    with Path(path).open() as f:
        rows=list(csv.DictReader(f))
    result=[]
    for row in rows:
        obj={}
        for k,v in row.items():
            if v=='':obj[k]=None
            elif v in {'True','False'}:obj[k]=v=='True'
            else:
                try:obj[k]=float(v)
                except ValueError:obj[k]=v
        result.append(obj)
    return result


def bootstrap_mean(values,seed=91221):
    x=np.asarray(values)
    rng=np.random.Generator(np.random.PCG64(seed))
    means=x[rng.integers(0,len(x),size=(2000,len(x)))].mean(axis=1)
    return float(x.mean()),float(np.quantile(means,.025)),float(np.quantile(means,.975))


def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();root,out=a.results,a.output;out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2,'savefig.dpi':180})
    files=[]
    def save(fig,name):
        fig.savefig(out/name,bbox_inches='tight');plt.close(fig);files.append(name)
    research=root/'research';config=json.loads((research/'manifest.json').read_text())['identity']['config']
    summary=read(research/'summary.csv');trials=read(research/'trial_summary.csv');outcomes=read(research/'outcomes.csv');detection=read(research/'detection_summary.csv')
    workloads=[w['name'] for w in config['workloads']];names=[s['name'] for s in config['scenarios']]
    tolerance={w['name']:w['thresholds'][0] for w in config['workloads']}
    fig,axes=plt.subplots(1,3,figsize=(16,6.5),sharey=True,layout='constrained')
    for ax,workload in zip(axes,workloads):
        for i,name in enumerate(names):
            r=next(r for r in summary if r['scenario']==name and r['workload']==workload)
            lo=r['mean_ci95_low']-r['target'];hi=r['mean_ci95_high']-r['target']
            ax.errorbar(r['bias'],i,xerr=[[r['bias']-lo],[hi-r['bias']]],fmt='o',capsize=3)
        ax.axvline(0,color='black',lw=1)
        ax.axvspan(-tolerance[workload],tolerance[workload],color='green',alpha=.08)
        ax.set_xscale('symlog',linthresh=tolerance[workload])
        locator=SymmetricalLogLocator(linthresh=tolerance[workload],base=10)
        locator.set_params(numticks=5)
        ax.xaxis.set_major_locator(locator)
        if workload=='vmc':
            ax.set_xticks([-.01,0,.1,10,1000,100000],labels=['−.01','0','.1','10','10³','10⁵'])
            ax.tick_params(axis='x',labelsize=9)
        elif workload=='ising':
            ax.set_xticks([-.1,0,.1],labels=['−.1','0','.1'])
        ax.set(title=workload,xlabel='Bias · symmetric log scale · pointwise 95% interval',yticks=np.arange(len(names)),yticklabels=names)
    axes[0].invert_yaxis();fig.suptitle('Synthetic faults and recorded-source comparisons · 100 trials each',fontsize=15)
    save(fig,'physical_comparison.png')

    def error_comparison(dataset,scenario_names,title,filename):
        rows=read(dataset/'trial_summary.csv');sums=read(dataset/'summary.csv')
        fig,axes=plt.subplots(1,3,figsize=(16,max(5,len(scenario_names)*.34)),sharey=True,layout='constrained')
        for ax,workload in zip(axes,workloads):
            for i,name in enumerate(scenario_names):
                group=[r for r in rows if r['scenario']==name and r['workload']==workload]
                mean,lo,hi=bootstrap_mean([abs(r['error'])/tolerance[workload] for r in group])
                fixed=next(r for r in sums if r['scenario']==name and r['workload']==workload)['unique_streams']<len(group)
                ax.plot(mean,i,'o',color='C1' if fixed else 'C0')
                if not fixed:ax.errorbar(mean,i,xerr=[[mean-lo],[hi-mean]],fmt='none',color='C0',capsize=3)
            ax.axvline(1,color='black',ls='--')
            ax.set(title=workload,xscale='log',xlabel='Mean absolute error / first tolerance',yticks=np.arange(len(scenario_names)),yticklabels=scenario_names)
            low,high=ax.get_xlim()
            if high/low<100:
                ax.xaxis.set_major_locator(LogLocator(base=10,subs=(1,2,5),numticks=6))
                ax.xaxis.set_major_formatter(FuncFormatter(lambda value,pos:f'{value:g}'))
            else:
                ax.xaxis.set_major_locator(LogLocator(base=10,numticks=6))
                ax.xaxis.set_major_formatter(LogFormatterMathtext(base=10))
            ax.xaxis.set_minor_formatter(NullFormatter())
        axes[0].invert_yaxis();fig.suptitle(title+'\n95% bootstrap intervals; orange points repeat a stream and have no interval',fontsize=13)
        save(fig,filename)
    error_comparison(research,names,'Application error normalized by workload-specific tolerances','normalized_error.png')

    fig,axes=plt.subplots(1,3,figsize=(16,6.5),sharey=True,layout='constrained')
    for ax,workload in zip(axes,workloads):
        group=[r for r in detection if r['workload']==workload and r['threshold']==tolerance[workload]]
        for i,name in enumerate(names):
            r=next(r for r in group if r['scenario']==name)
            for offset,prefix,color in [(-.13,'detection','C0'),(.13,'failure','C1')]:
                value=r[prefix+'_probability'];lo=r[prefix+'_ci_low'];hi=r[prefix+'_ci_high']
                ax.errorbar(value,i+offset,xerr=[[value-lo],[hi-value]],fmt='o',color=color,capsize=2,label=prefix if i==0 else None)
        ax.set(title=workload,xlim=(-.03,1.03),xlabel='Observed probability · exact 95% binomial interval',yticks=np.arange(len(names)),yticklabels=names)
        ax.legend()
    axes[0].invert_yaxis();fig.suptitle('Output detection and sustained application failure · clean crossing rates are controls',fontsize=14)
    save(fig,'detection_vs_failure.png')

    selected=[n for n in ['synthetic','bias_onset','correlation_ramp','freeze_burst','output_bias','partial_loss'] if n in names]
    fig,axes=plt.subplots(2,3,figsize=(17,9),layout='constrained')
    for col,workload in enumerate(workloads):
        vectors=[];positions=[]
        for i,name in enumerate(selected):
            group=[r for r in outcomes if r['scenario']==name and r['workload']==workload and r['threshold']==tolerance[workload]]
            delta=[r['delta_t'] for r in group if r['delta_t'] is not None]
            if delta:vectors.append(delta);positions.append(i)
            counts=[sum(r['status']==s for r in group) for s in ['both','failure_only','detection_only','neither']]
            left=0
            for j,(count,status) in enumerate(zip(counts,['both','failure only','detection only','neither'])):
                axes[1,col].barh(i,count,left=left,color=plt.get_cmap('tab10')(j),label=status if i==0 else None)
                if count>=3:axes[1,col].text(left+count/2,i,str(count),ha='center',va='center',fontsize=8)
                left+=count
            axes[0,col].text(.98,i,f'{len(delta)}/{len(group)} both',transform=axes[0,col].get_yaxis_transform(),ha='right',fontsize=8)
        if vectors:axes[0,col].boxplot(vectors,positions=positions,orientation='horizontal',showfliers=True,widths=.4,manage_ticks=False)
        for row in range(2):
            axes[row,col].set(yticks=np.arange(len(selected)),yticklabels=selected,title=workload)
            axes[row,col].invert_yaxis()
        axes[0,col].axvline(0,color='black',ls='--');axes[0,col].set_xlabel('ΔT = confirmed failure − detection (proposals)')
        axes[1,col].set_xlabel('Trials; all censoring/event categories retained')
    axes[1,0].legend(fontsize=8,loc='lower right')
    fig.suptitle('Lead-time distribution conditional on both events · counts show censoring',fontsize=15)
    save(fig,'detection_delays.png')

    rows=read(research/'windows.csv');chosen=['synthetic','bias_onset','correlation_ramp','freeze_burst','output_bias']
    fig,axes=plt.subplots(2,3,figsize=(17,8),layout='constrained')
    for col,workload in enumerate(workloads):
        for name in chosen:
            groups=defaultdict(list)
            for r in rows:
                if r['scenario']==name and r['workload']==workload:groups[r['time']].append(r)
            x=sorted(groups)
            axes[0,col].plot(x,[np.mean([r['bit_bias'] for r in groups[n]]) for n in x],label=name)
            mx=[];lo=[];hi=[]
            x2=[]
            for n in x:
                vals=[r['absolute_error']/tolerance[workload] for r in groups[n] if r['absolute_error'] is not None]
                if vals:
                    mean,l,h=mean_interval(vals);x2.append(n);mx.append(mean);lo.append(max(0,l));hi.append(h)
            line,=axes[1,col].plot(x2,mx,label=name)
            axes[1,col].fill_between(x2,lo,hi,color=line.get_color(),alpha=.1)
        for row in range(2):axes[row,col].axvline(6144,color='black',ls='--',lw=1)
        axes[0,col].set(title=workload,ylabel='Mean bias of actual output bits')
        axes[1,col].axhline(1,color='black',ls=':');axes[1,col].set_yscale('symlog',linthresh=1)
        axes[1,col].set(xlabel='Total proposals (including burn-in)',ylabel='Mean cumulative |error| / tolerance')
        axes[1,col].set_ylim(bottom=0)
    axes[0,0].legend(fontsize=8);fig.suptitle('Fault trajectories and application response · shading: pointwise 95% trial intervals',fontsize=14)
    save(fig,'degradation_trajectories.png')

    monitor_groups=defaultdict(list)
    for path in (research/'trials').glob('*.json'):
        if not any(path.name.startswith(name+'__') for name in chosen):continue
        payload=json.loads(path.read_text())['payload'];r=payload['result']
        for m in payload['monitor']:
            monitor_groups[(r['scenario'],r['workload'],m['time'])].append(-np.log10(max(m['min_p'],1e-300)))
    fig,axes=plt.subplots(1,3,figsize=(16,4.5),layout='constrained')
    for ax,w in zip(axes,workloads):
        for name in chosen:
            x=sorted(k[2] for k in monitor_groups if k[:2]==(name,w))
            ax.plot(x,[np.median(monitor_groups[(name,w,n)]) for n in x],label=name)
        cutoff=config['monitor_alpha']/(np.ceil((config['samples']+config['burn_in'])/config['window'])*8)
        ax.axhline(-np.log10(cutoff),color='black',ls='--',label='corrected alarm cutoff')
        ax.axvline(6144,color='gray',ls=':');ax.set(title=w,xlabel='Total proposals',ylabel='Median −log10(minimum test p); capped at 300')
    axes[0].legend(fontsize=8);fig.suptitle('Actual consumed-bit test trajectories · correlations can be visible with zero bit bias')
    save(fig,'monitor_trajectories.png')

    exploratory=read(research/'exploratory_event_summary.csv')
    fig,axes=plt.subplots(1,3,figsize=(13,4.5),layout='constrained')
    for ax,w in zip(axes,workloads):
        row=next(r for r in exploratory if r['scenario']=='correlation_ramp' and r['workload']==w and r['threshold']==tolerance[w])
        numbers=[row['monobit_detected'],100,row['post_failure_before_monobit'],row['post_failure_before_composite']]
        ax.bar(np.arange(4),np.array(numbers)/100,color=['C0','C1','C2','C3'])
        for i,n in enumerate(numbers):ax.text(i,n/100+.025,f'{int(n)}/100',ha='center',fontsize=9)
        ax.set(title=w,ylim=(0,1.13),xticks=np.arange(4),xticklabels=['frequency\ncomponent','composite\nmonitor','failure before\nfrequency','failure before\ncomposite'],ylabel='Observed fraction')
        ax.tick_params(axis='x',labelsize=8)
    fig.suptitle('Correlation-ramp blind spot · post-onset timing is an exploratory sensitivity analysis',fontsize=14)
    save(fig,'frequency_blindspot.png')

    mitigation=root/'mitigation';mcfg=json.loads((mitigation/'manifest.json').read_text())['identity']['config']
    comb=['synthetic','xor_one_frozen','xor_both_frozen','xor_identical','xor_attacker_copy','hash_one_frozen','hash_both_frozen','combined_seed','combined_reseed']
    error_comparison(mitigation,comb,'Source combinations · synthetic independent and adaptive-copy models','entropy_combinations.png')
    reseed=['synthetic','frozen_initial','frozen_initial_recover','frozen_initial_triggered','recovering_reseed','frequent_reseed','triggered_output_fault']
    error_comparison(mitigation,reseed,'Initialization and reseeding · exploratory 24-trial comparisons','reseeding_resilience.png')
    mt=read(mitigation/'trial_summary.csv');extra=[r for r in trials if r['scenario'] in ['csprng','hybrid_healthy','hybrid_frozen']]
    snames=['csprng','hybrid_healthy','hybrid_frozen','frozen_initial','frozen_initial_recover','frozen_initial_triggered','recovering_reseed']
    fig,axes=plt.subplots(1,3,figsize=(16,5),sharey=True,layout='constrained')
    for ax,w in zip(axes,workloads):
        for i,name in enumerate(snames):
            group=[r for r in mt+extra if r['scenario']==name and r['workload']==w]
            ax.barh(i-.16,np.mean([r['first_alarm'] is not None for r in group]),height=.28,color='C0',label='output alarm' if i==0 else None)
            ax.barh(i+.16,np.mean([r['source_first_alarm'] is not None for r in group]),height=.28,color='C1',label='actual seed/reseed alarm' if i==0 else None)
        ax.set(title=w,xlim=(0,1.05),xlabel='Descriptive alarm fraction',yticks=np.arange(len(snames)),yticklabels=snames)
    axes[0].invert_yaxis();axes[0].legend(fontsize=8)
    fig.suptitle('Source health alarms can coexist with ordinary-looking CSPRNG output\nFirst three cases: 100 trials; initialization/recovery cases: 24 trials',fontsize=13)
    save(fig,'source_vs_output_alarms.png')

    sensitivity=root/'sensitivity'
    if (sensitivity/'summary.csv').exists():
        scfg=json.loads((sensitivity/'manifest.json').read_text())['identity']['config'];sr=read(sensitivity/'trial_summary.csv')
        fig,axes=plt.subplots(1,2,figsize=(13,4.5),layout='constrained')
        for ax,kind in zip(axes,['bias','autocorrelation']):
            scenarios=[s for s in scfg['scenarios'] if s.get('mechanism')==kind]
            for w in workloads:
                points=[('synthetic',0)]+[(s['name'],s['severity']) for s in scenarios]
                means=[];lower=[];upper=[]
                for name,_ in points:
                    m,l,h=bootstrap_mean([abs(r['error'])/tolerance[w] for r in sr if r['scenario']==name and r['workload']==w])
                    means.append(m);lower.append(m-l);upper.append(h-m)
                ax.errorbar([s for _,s in points],means,yerr=[lower,upper],marker='o',capsize=3,label=w)
            ax.set(title=kind,yscale='log',xlabel='Mechanism-specific strength',ylabel='Mean absolute error / first tolerance');ax.legend()
        fig.suptitle('Exploratory physical-workload severity sweep · 24 trials and 95% bootstrap intervals')
        save(fig,'physical_severity.png')

    pimc=json.loads((root/'pimc-validation'/'summary.json').read_text())['cases']
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    x=np.array([1/r['beads']**2 for r in pimc]);means=np.array([r['mean'] for r in pimc])
    axes[0].errorbar(x,means,yerr=[[r['mean']-r['ci999_low'] for r in pimc],[r['ci999_high']-r['mean'] for r in pimc]],fmt='o',capsize=4,label='PIMC · 99.9% trial intervals')
    axes[0].plot(x,[r['finite_target'] for r in pimc],marker='.',label='exact primitive discretization')
    axes[0].axhline(pimc[0]['continuum_target'],color='black',ls='--',label='continuum quantum energy')
    axes[0].set(xlabel='1 / imaginary-time beads²',ylabel='Thermal energy',title='β=2 harmonic oscillator');axes[0].legend(fontsize=8)
    axes[1].loglog([r['beads'] for r in pimc],[abs(r['discretization_error']) for r in pimc],marker='o')
    axes[1].set(xlabel='Imaginary-time beads',ylabel='Exact |finite-discretization error|',title='Sampling error is separate from discretization')
    fig.suptitle('Genuine path-integral Quantum Monte Carlo validation · 48 trials total')
    save(fig,'pimc_discretization.png')
    validation=json.loads((root/'model-validation'/'summary.json').read_text())['cases']
    vmc=[r for r in validation if r['case'].startswith('vmc')]
    alpha=np.array([float(r['case'].split('_')[-1]) for r in vmc])
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    axes[0].errorbar(alpha,[r['mean'] for r in vmc],yerr=[[r['mean']-r['ci999_low'] for r in vmc],[r['ci999_high']-r['mean'] for r in vmc]],fmt='o',capsize=3,label='validated trials')
    xx=np.linspace(.45,1.5,100);axes[0].plot(xx,(xx+1/xx)/4,label='exact variational energy');axes[0].axhline(.5,ls='--',color='black',label='quantum ground state')
    axes[0].set(xlabel='Variational parameter α',ylabel='Quantum energy');axes[0].legend(fontsize=8)
    ising=[r for r in validation if r['case'].startswith('ising')]
    err=[r['mean']-r['target'] for r in ising]
    axes[1].errorbar(np.arange(len(ising)),err,yerr=[[r['mean']-r['ci999_low'] for r in ising],[r['ci999_high']-r['mean'] for r in ising]],fmt='o',capsize=3)
    axes[1].axhline(0,color='black',ls='--');axes[1].set(xticks=np.arange(len(ising)),xticklabels=[r['case'] for r in ising],ylabel='Energy error versus exact enumeration')
    fig.suptitle('Independent physical-model validation · 12 trials per case; 99.9% trial intervals')
    save(fig,'physical_validation.png')
    (out/'plot_metadata.json').write_text(json.dumps({'files':files,'matplotlib':matplotlib.__version__,'numpy':np.__version__,'bootstrap_seed':91221,'bootstrap_replicates':2000,'bootstrap_use':'analysis only; simulation draws use providers','uncertainty':'pointwise unless explicitly labeled familywise; fixed streams have no across-stream intervals'},indent=2))
    print(f'Saved {len(files)} figures to {out}')

if __name__=='__main__':main()
