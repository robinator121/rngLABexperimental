"""Independent-seed physical-model validation before research experiments."""
import argparse,csv,json,time
from pathlib import Path
import numpy as np
from rnglab.providers import PCGSource,DirectProvider
from rnglab.v02.physics import VMC,Ising,ising_exact
from rnglab.v02.statistics import correlated_summary,mean_interval


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    a.output.mkdir(parents=True,exist_ok=False)
    rows=[];start=time.perf_counter()
    cases=[(f'vmc_alpha_{alpha}',lambda alpha=alpha:VMC(alpha,x0=4)) for alpha in [.5,.7,1.,1.4]]
    cases += [(f'ising_L{L}_T{T}',lambda L=L,T=T:Ising(L,T,initial='checkerboard')) for L,T in [(3,2.),(3,3.),(4,3.)]]
    for case,factory in cases:
        for trial in range(12):
            m=factory();u,_=DirectProvider(PCGSource([82291,trial])).sample(36864,2)
            e,aux=m.advance(u);stats=correlated_summary(e[4096:])
            rows.append({'case':case,'trial':trial,'target':m.target,'estimate':stats['mean'],
                         'mcse':stats['mcse'],'tau':stats['tau'],'ess':stats['ess'],
                         'acceptance':float(aux['accepted'][4096:].mean()),
                         'abs_magnetization':float(aux['abs_magnetization'][4096:].mean()) if 'abs_magnetization' in aux else None,
                         'abs_m_target':ising_exact(m.length,m.temperature)['abs_magnetization'] if isinstance(m,Ising) else None})
        print(case,'done',flush=True)
    with (a.output/'trials.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=[]
    for case,_ in cases:
        group=[r for r in rows if r['case']==case]
        mean,lo,hi=mean_interval([r['estimate'] for r in group],confidence=.999)
        target=group[0]['target']
        summary.append({'case':case,'target':target,'mean':mean,'ci999_low':lo,'ci999_high':hi,'target_in_ci':bool(lo-1e-12<=target<=hi+1e-12)})
    (a.output/'summary.json').write_text(json.dumps({'trials':84,'burn_in':4096,'retained':32768,'seconds':time.perf_counter()-start,'cases':summary},indent=2))
    print(json.dumps(summary,indent=2))
    if not all(r['target_in_ci'] for r in summary):
        raise SystemExit('Validation requires investigation: target outside 99.9% trial interval')

if __name__=='__main__':main()
