"""Plots saved measurements; secondary fixed-width checks are explicitly post hoc."""
import csv
import json
import math
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import norm,t,beta as beta_dist
import analytic_dependence as a


def main():
    root=Path(__file__).resolve().parent
    data=json.loads((root/'results.json').read_text());rows=data['rows']
    colors=['#3465a4','#c66b16','#167f62'];labels=['IID centroid, φ=0','AR(1) centroid, φ=0.997','AR(1) centroid, φ=0.999']
    fig,axes=plt.subplots(1,2,figsize=(12,5),layout='constrained')
    methods=['BM32','BM128','BM512','oracle_known_variance_normal'];x=np.arange(4)
    for i,(phi,color,label) in enumerate(zip(a.PHIS,colors,labels)):
        selected=[next(r for r in rows if r['phi']==phi and r['method']==m) for m in methods]
        coverage=np.array([r['coverage'] for r in selected]);ci=np.array([r['coverage_ci95'] for r in selected])
        xx=x+(i-1)*.12
        axes[0].errorbar(xx,coverage,yerr=[coverage-ci[:,0],ci[:,1]-coverage],fmt='o-',color=color,capsize=3,label=label,lw=1.3)
        bm=selected[:3];emp=np.array([r['estimated_to_empirical_variance_ratio'] for r in bm]);ci=np.array([r['estimated_to_empirical_variance_ratio_bootstrap_ci95'] for r in bm])
        theory=np.array([r['exact_expected_BM_to_true_MCvariance_ratio'] for r in bm])
        axes[1].errorbar(np.arange(3)+(i-1)*.08,emp,yerr=[emp-ci[:,0],ci[:,1]-emp],fmt='o',color=color,capsize=3)
        axes[1].plot(np.arange(3),theory,'--',color=color,label=label+' (analytic)')
    axes[0].axhline(.95,color='black',ls='--',lw=1,label='Nominal 95%')
    axes[0].set(xticks=x,xticklabels=['Batch 32','Batch 128','Batch 512','Exact variance\nnormal interval'],ylabel='Observed coverage',ylim=(.15,1.01),title='Measured 95% interval coverage; CP 95% bars')
    axes[1].axhline(1,color='black',ls=':',lw=1)
    axes[1].set(xticks=np.arange(3),xticklabels=['Batch 32','Batch 128','Batch 512'],ylabel='Variance ratio',ylim=(0,1.22),title='Points: estimated / empirical + bootstrap 95%\nDashed: E[estimated] / exact, no fitted parameters')
    for ax in axes:ax.grid(alpha=.2)
    axes[0].legend(fontsize=8,loc='lower right');axes[1].legend(fontsize=7,loc='upper right')
    fig.suptitle('Harmonic target with intentionally correlated centroid updates — 1024 independent replicas per condition',fontsize=12)
    for ext in ('svg','png'):fig.savefig(root/f'analytic_dependence_coverage_variance.{ext}',dpi=180)
    plt.close(fig)
    # Post-hoc checks use existing data only. They replace the noisy estimated
    # variance with its exact expectation, retaining each BM method's t critical.
    secondary=[]
    fig,axs=plt.subplots(1,3,figsize=(12,4),layout='constrained')
    for ax,phi,color,label in zip(axs,a.PHIS,colors,labels):
        condition='phi'+str(phi).replace('.','p');y=np.load(root/f'{condition}-CV-traces.npz')['CV'];means=y.mean(1)
        p=a.parameters();exact=a.exact_mc_variance(a.N,phi,p);standardized=(means-p['energy_mean'])/math.sqrt(exact)
        xx=np.linspace(-3,5,400);ax.hist(standardized,bins=35,density=True,color=color,alpha=.5,label='Measured replica means')
        ax.plot(xx,norm.pdf(xx),'k--',lw=1.5,label='Standard normal comparison')
        ax.set(title=f'{label}\nExact skewness = {a.exact_mean_skewness(a.N,phi,p):.3f}',xlabel='(Mean CV − exact mean) / exact MCSE',ylabel='Density',xlim=(-3,5));ax.grid(alpha=.2)
        for b in a.BATCHES:
            blocks=a.N//b;critical=float(t.ppf(.975,blocks-1));expected=a.exact_expected_bm_mcvariance(a.N,b,phi,p)
            hits=int(np.count_nonzero(np.abs(means-p['energy_mean'])<=critical*math.sqrt(expected)))
            oracle_t_hits=int(np.count_nonzero(np.abs(means-p['energy_mean'])<=critical*math.sqrt(exact)))
            estimate=y.reshape(len(y),blocks,b).mean(2).var(1,ddof=1)/blocks
            actual_hits=np.abs(means-p['energy_mean'])<=critical*np.sqrt(estimate)
            below=int(np.count_nonzero(means+critical*np.sqrt(estimate)<p['energy_mean']))
            above=int(np.count_nonzero(means-critical*np.sqrt(estimate)>p['energy_mean']))
            secondary.append({'condition':condition,'phi':phi,'batch':b,'critical_t':critical,
                'classification':'post-hoc deterministic diagnostic on saved trajectories; not a newly held-out comparison',
                'actual_BM_coverage_count':int(actual_hits.sum()),
                'actual_BM_interval_entirely_below_target_count':below,
                'actual_BM_interval_entirely_above_target_count':above,
                'correlation_estimated_variance_with_energy_mean':float(np.corrcoef(estimate,means)[0,1]),
                'fixed_expected_BM_variance_interval':{'variance':expected,'coverage_count':hits,'replicas':len(y),'coverage':hits/len(y),'CP95':a.binomial_ci(hits,len(y))},
                'fixed_exact_MCvariance_same_t_interval':{'variance':exact,'coverage_count':oracle_t_hits,'replicas':len(y),'coverage':oracle_t_hits/len(y),'CP95':a.binomial_ci(oracle_t_hits,len(y))}})
        del y
    axs[0].legend(fontsize=8)
    fig.suptitle('Exact variance does not make finite-run energy averages Gaussian',fontsize=12)
    for ext in ('svg','png'):fig.savefig(root/f'analytic_dependence_mean_distributions.{ext}',dpi=180)
    plt.close(fig)
    (root/'secondary_fixed_width_diagnostics.json').write_text(json.dumps(secondary,indent=2)+'\n')


if __name__=='__main__':main()
