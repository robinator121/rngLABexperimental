"""Validated finite-temperature harmonic oscillator PIMC extension (separate budget)."""
import argparse,json,math,time
from pathlib import Path
import numpy as np
from rnglab.providers import DirectProvider,PCGSource
from rnglab.v02.statistics import correlated_summary,mean_interval
from rnglab.v02.runner import csv_write


def thermal_energy(beta):return .5/math.tanh(beta/2)


def finite_energy(beta,beads):
    eps=beta/beads
    return .5/(math.sqrt(1+eps*eps/4)*math.tanh(beads*math.asinh(eps/2)))


def action(path,beta):
    eps=beta/len(path)
    return float(np.sum((np.roll(path,-1)-path)**2)/(2*eps)+eps*np.sum(path**2)/2)


class PIMC:
    """Primitive ring-polymer action, random bead / centroid moves, centroid-virial energy."""
    def __init__(self,beta=2.,beads=16):
        if beta<=0 or beads<3:raise ValueError('beta>0 and beads>=3 required')
        self.beta,self.beads=beta,beads;self.path=np.zeros(beads);self.sum=0.;self.sumsq=0.;self.sumdiff=0.

    def advance(self,uniforms):
        if uniforms.shape[1]!=3:raise ValueError('PIMC consumes three uniforms per proposal')
        M,beta=self.beads,self.beta;eps=beta/M
        output=np.empty(len(uniforms));thermo=np.empty(len(uniforms));accepted=np.zeros(len(uniforms),dtype=bool)
        s,q,diff=self.sum,self.sumsq,self.sumdiff;p=self.path
        for k,(us,up,ua) in enumerate(uniforms):
            site=int(us*(M+1))
            width=1. if site==M else 2*math.sqrt(eps)
            delta=width*(2*up-1)
            if site==M:
                dq=2*delta*s+M*delta*delta;ddiff=0.;ds=M*delta
            else:
                old=p[site];new=old+delta
                dq=new*new-old*old;ds=delta
                left,right=p[(site-1)%M],p[(site+1)%M]
                ddiff=(new-left)**2+(right-new)**2-(old-left)**2-(right-old)**2
            delta_action=ddiff/(2*eps)+eps*dq/2
            if delta_action<=0 or ua==0 or math.log(ua)<-delta_action:
                if site==M:p+=delta
                else:p[site]+=delta
                s+=ds;q+=dq;diff+=ddiff;accepted[k]=True
            output[k]=1/(2*beta)+q/M-.5*(s/M)**2
            thermo[k]=M/(2*beta)-M*diff/(2*beta*beta)+q/(2*M)
        self.sum,self.sumsq,self.sumdiff=s,q,diff
        return output,thermo,accepted


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=False)
    rows=[];start=time.perf_counter();burn=16384;retained=131072
    for M in [4,8,16,32]:
        for trial in range(12):
            provider=DirectProvider(PCGSource([929882,M,trial]))
            model=PIMC(2.,M);parts=[];therm=[];acc=[]
            for offset in range(0,burn+retained,4096):
                u,_=provider.sample(4096,3);e,t,a=model.advance(u)
                if offset>=burn:parts.append(e);therm.append(t);acc.append(a)
            values=np.concatenate(parts);stats=correlated_summary(values)
            rows.append({'beads':M,'beta':2.,'trial':trial,'retained':retained,'burn_in':burn,'uniforms':3*(retained+burn),
                         'finite_target':finite_energy(2.,M),'continuum_target':thermal_energy(2.),
                         'discretization_error':finite_energy(2.,M)-thermal_energy(2.),'estimate':stats['mean'],
                         'thermodynamic_estimate':float(np.concatenate(therm).mean()),'acceptance':float(np.concatenate(acc).mean()),
                         'tau':stats['tau'],'ess':stats['ess'],'mcse':stats['mcse']})
        print(f'PIMC M={M}: 12 trials complete',flush=True)
    csv_write(args.output/'trials.csv',rows)
    summary=[]
    for M in [4,8,16,32]:
        group=[r for r in rows if r['beads']==M];mean,lo,hi=mean_interval([r['estimate'] for r in group],.999)
        thermo,tlo,thi=mean_interval([r['thermodynamic_estimate'] for r in group],.999)
        target=finite_energy(2.,M)
        summary.append({'beads':M,'mean':mean,'ci999_low':lo,'ci999_high':hi,'finite_target':target,'continuum_target':thermal_energy(2.),
                        'discretization_error':target-thermal_energy(2.),'target_in_ci':lo<=target<=hi,
                        'thermodynamic_mean':thermo,'thermodynamic_target_in_ci':tlo<=target<=thi})
    (args.output/'summary.json').write_text(json.dumps({'trials':48,'seconds':time.perf_counter()-start,'seed_scheme':'[929882,beads,trial] PCG64','cases':summary},indent=2))
    print(json.dumps(summary,indent=2))
    if not all(r['target_in_ci'] and r['thermodynamic_target_in_ci'] for r in summary):raise SystemExit('PIMC validation needs investigation')

if __name__=='__main__':main()
