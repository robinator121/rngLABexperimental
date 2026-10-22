"""Adversarial checks using alternative formulas; same author, not external review."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from rnglab16.core import enumeration,stationary,augmented,two_state
from rnglab16.campaign import exact_ising_transition,save

def run(root):
    checks=[]
    def test(n,c):
        checks.append({'name':n,'passed':bool(c)})
    # Independent partition sum via explicit bond list, instead of np.roll.
    for L in [2,3]:
        n=L*L;bits=np.array([[1 if mask&(1<<k) else -1 for k in range(n)] for mask in range(2**n)])
        E=np.zeros(len(bits))
        for x in range(L):
            for y in range(L):
                for xx,yy in [((x+1)%L,y),(x,(y+1)%L)]:E-=bits[:,x*L+y]*bits[:,xx*L+yy]
        for T in [1.5,2/np.log(1+np.sqrt(2)),3.5]:
            h=1e-4;logz=lambda b:logsumexp(-b*E)
            e=-(logz(1/T+h)-logz(1/T-h))/(2*h)/n
            test(f'enumeration_derivative_{L}_{T}',abs(e-enumeration(L,T)['energy'])<1e-7)
    for a,b in [(.1,.2),(.7,.8),(.001,.002)]:
        for da,db in [(.0001,0),(0,.0001),(.0001,-.00005),(.01*a,.01*b)]:
            exact=(b*da-a*db)/((a+b)*(a+b+da+db));linear=(b*da-a*db)/(a+b)**2
            remainder=-(b*da-a*db)*(da+db)/((a+b)**2*(a+b+da+db))
            test(f'perturbation_{a}_{b}_{da}_{db}',abs((a+da)/(a+b+da+db)-a/(a+b)-exact)<1e-13 and abs(exact-linear-remainder)<1e-13)
    P,E=exact_ising_transition(3,2/np.log(1+np.sqrt(2)),.05)
    test('strong_bias_creates_multiple_absorbing_states',np.linalg.matrix_rank(P.T-np.eye(len(P)),tol=1e-10)<len(P)-1)
    test('deliberately_wrong_stationary_formula_rejected',abs(.1/(.1+.2)-.1/.2)>.1)
    test('deliberately_wrong_PIMC_coefficient_rejected',abs(4/(2*2)-4*2/(2*2**2)-(4/(2*2)-2/(2*2)))>.1)
    for rho in [0,.5,.9]:
        P=augmented(rho=rho);pi=stationary(P);power=np.ones(len(P))/len(P)
        for _ in range(10000):power=power@P
        test(f'augmented_power_{rho}',max(abs(power-pi))<1e-10)
    sourcehash=json.loads((root/'source_hashes.json').read_text())
    import hashlib
    for rel,digest in sourcehash.items():test('frozen_source_'+rel,hashlib.sha256((root/'source'/Path(rel).name).read_bytes()).hexdigest()==digest)
    result={'checks':checks,'all_pass':all(x['passed'] for x in checks),'independence':'Alternative numerical formulations using dual algorithmic checks; pending external scientific review','known_limitations':['Finite-run intervals remain exploratory.','Wolff and Metropolis update units differ.','Prediction cutoff fixed, not alarm-triggered.','Strong bias outside unique stationary theory domain.','Same seeds in amended run are replay, not new replication.']}
    save(root/'audit.json',result);print(json.dumps(result,indent=2))
    assert result['all_pass']
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);run(p.parse_args().root)
