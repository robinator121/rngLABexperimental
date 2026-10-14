"""Supplied-uniform anharmonic VMC and independently written classical sampler.

The quantum trial is exp(-alpha*x²/2-trial_quartic*x⁴/4). Its sampling
distribution is independent of the physical quartic coupling. The classical
control deliberately matches that sampling density; it is not an independent
physical quantum calculation. Historical kernels remain untouched.
"""
from __future__ import annotations
import ctypes
from functools import lru_cache
import hashlib
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np
from scipy.integrate import quad

COMPILER_FLAGS = ['-O3','-std=c11','-shared','-fPIC']
C_SOURCE = r'''
#include <math.h>
void quantum(long n,const double *u,double a,double b,double width,
             double *state,double *xout,double *accepted) {
  double x=*state;
  for(long i=0;i<n;i++) {
    double y=x+width*(2*u[2*i]-1);
    double xx=x*x,yy=y*y;
    double logratio=-a*(yy-xx)-b/2*(yy*yy-xx*xx);
    if(logratio>=0 || u[2*i+1]==0 || log(u[2*i+1])<logratio) {
      x=y;accepted[i]=1;
    } else accepted[i]=0;
    xout[i]=x;
  }
  *state=x;
}
/* Written in classical potential differences rather than wavefunction ratios. */
void classical(long n,const double *u,double a,double b,double width,
               double *state,double *xout,double *accepted) {
  double current=*state;
  for(long i=0;i<n;i++) {
    double candidate=current+width*(2*u[2*i]-1);
    double s0=current*current,s1=candidate*candidate;
    double old_potential=a*s0+0.5*b*s0*s0;
    double new_potential=a*s1+0.5*b*s1*s1;
    double change=new_potential-old_potential;
    if(change<=0 || u[2*i+1]==0 || log(u[2*i+1]) < -change) {
      current=candidate;accepted[i]=1;
    } else accepted[i]=0;
    xout[i]=current;
  }
  *state=current;
}
'''
KERNEL_SHA256=hashlib.sha256(C_SOURCE.encode()).hexdigest()

@lru_cache(maxsize=1)
def compiled_kernels():
    compiler=shutil.which('cc') or shutil.which('gcc')
    if compiler is None:
        return None
    root=Path(__file__).resolve().parents[3]
    cache=root/'.cache/v04-kernels';cache.mkdir(parents=True,exist_ok=True)
    identity=hashlib.sha256((C_SOURCE+' '.join(COMPILER_FLAGS)).encode()).hexdigest()
    library=cache/(identity+'.so')
    if not library.exists():
        with tempfile.TemporaryDirectory(dir=cache) as tmp:
            source=Path(tmp)/'models.c';output=Path(tmp)/'models.so'
            source.write_text(C_SOURCE)
            subprocess.run([compiler,*COMPILER_FLAGS,str(source),'-o',str(output),'-lm'],check=True,capture_output=True)
            os.replace(output,library)
    lib=ctypes.CDLL(str(library));pointer=ctypes.POINTER(ctypes.c_double)
    for name in ['quantum','classical']:
        function=getattr(lib,name)
        function.argtypes=[ctypes.c_long,pointer,ctypes.c_double,ctypes.c_double,ctypes.c_double,pointer,pointer,pointer]
        function.restype=None
    return lib

@lru_cache(maxsize=32)
def intended_moments(alpha,trial_quartic):
    if trial_quartic==0:
        v=1/(2*alpha)
        return {'normalizer':math.sqrt(math.pi/alpha),'x':0.,'x2':v,'x3':0.,
                'x4':3*v*v,'x6':15*v**3,'sign':0.}
    def weighted(power):
        return 2*quad(lambda x:x**power*math.exp(-alpha*x*x-.5*trial_quartic*x**4),
                      0,np.inf,epsabs=1e-12,epsrel=1e-12)[0]
    z=weighted(0)
    return {'normalizer':z,'x':0.,'x2':weighted(2)/z,'x3':0.,
            'x4':weighted(4)/z,'x6':weighted(6)/z,'sign':0.}

def quantum_local_energy(x,alpha,coupling,trial_quartic=0.):
    x2=np.asarray(x)**2
    return alpha/2+(1-alpha*alpha+3*trial_quartic)*x2/2+(
        coupling-alpha*trial_quartic)*x2*x2-trial_quartic**2*x2*x2*x2/2

class Model:
    def __init__(self,spec):
        self.spec=dict(spec);self.kind=spec.get('kind','quantum')
        self.alpha=float(spec['alpha']);self.b=float(spec.get('trial_quartic',0.))
        self.coupling=float(spec.get('coupling',0.));self.proposal=float(spec.get('proposal',1.5))
        if self.kind not in {'quantum','classical'} or self.alpha<=0 or self.b<0 or self.coupling<0 or self.proposal<=0:
            raise ValueError('invalid finite normalizable model')
        self.x=0.;self.targets=dict(intended_moments(self.alpha,self.b))
        m=self.targets
        self.targets['energy']=(self.alpha/2+(1-self.alpha**2+3*self.b)*m['x2']/2+
                                (self.coupling-self.alpha*self.b)*m['x4']-self.b**2*m['x6']/2)
        self.targets['classical_potential']=self.alpha*m['x2']+.5*self.b*m['x4']
        backend=os.environ.get('RNGLAB_V04_BACKEND',spec.get('backend','auto'))
        self.library=None if backend=='python' else compiled_kernels()
        if backend=='compiled' and self.library is None:raise RuntimeError('compiled backend unavailable')
        self.backend='compiled' if self.library is not None else 'python'

    def advance(self,uniforms):
        u=np.asarray(uniforms,dtype=np.float64)
        if u.ndim!=2 or u.shape[1]!=2 or not np.isfinite(u).all() or np.any(u<0) or np.any(u>=1):
            raise ValueError('two finite uniforms in [0,1) required per proposal')
        u=np.ascontiguousarray(u);n=len(u);positions=np.empty(n);accepted=np.empty(n)
        if self.library is not None:
            pointer=ctypes.POINTER(ctypes.c_double);state=ctypes.c_double(self.x)
            getattr(self.library,self.kind)(n,u.ctypes.data_as(pointer),self.alpha,self.b,self.proposal,
                                          ctypes.byref(state),positions.ctypes.data_as(pointer),accepted.ctypes.data_as(pointer))
            self.x=state.value
        else:
            x=self.x
            for i,(up,ua) in enumerate(u):
                candidate=x+self.proposal*(2*up-1)
                if self.kind=='quantum':
                    xx,yy=x*x,candidate*candidate
                    ratio=-self.alpha*(yy-xx)-self.b/2*(yy*yy-xx*xx)
                else:
                    old=self.alpha*x*x+.5*self.b*(x*x)*(x*x)
                    new=self.alpha*candidate*candidate+.5*self.b*(candidate*candidate)*(candidate*candidate)
                    ratio=old-new
                take=ratio>=0 or ua==0 or math.log(ua)<ratio
                if take:x=candidate
                positions[i]=x;accepted[i]=take
            self.x=x
        x2=positions*positions
        return {'x':positions,'x2':x2,'x3':x2*positions,'x4':x2*x2,'x6':x2*x2*x2,
                'sign':np.sign(positions),'energy':quantum_local_energy(positions,self.alpha,self.coupling,self.b),
                'classical_potential':self.alpha*x2+.5*self.b*x2*x2,'accepted':accepted}

    def state_dict(self):return {'x':self.x}
    def load_state_dict(self,state):self.x=float(state['x'])
