"""New polynomial PIMC; the historical harmonic kernel remains untouched."""
import ctypes,hashlib,math,os,shutil,subprocess,tempfile
from functools import lru_cache
import numpy as np
from .common import ROOT
from rnglab.v03.physics import normal_modes,mode_probabilities

C_SOURCE=r'''
#include <math.h>
void polynomial(long n,const double *u,int M,double beta,double c4,double c2,double c1,double c0,double width,
 const double *lam,const double *basis,const double *cdf,double *path,double *coeff,double *out) {
 double eps=beta/M;
 for(long k=0;k<n;k++) {
  int mode=0;while(mode<M-1 && u[3*k]>=cdf[mode])mode++;
  double dc=width/sqrt(lam[mode])*(2*u[3*k+1]-1),newc=coeff[mode]+dc,dv=0;
  for(int j=0;j<M;j++) {
   double x=path[j],y=x+dc*basis[mode*M+j];
   dv+=c4*(y*y*y*y-x*x*x*x)+(c2-.5)*(y*y-x*x)+c1*(y-x);
  }
  double da=.5*lam[mode]*(newc*newc-coeff[mode]*coeff[mode])+eps*dv;
  int take=(da<=0 || u[3*k+2]==0 || log(u[3*k+2]) < -da);
  if(take){coeff[mode]=newc;for(int j=0;j<M;j++)path[j]+=dc*basis[mode*M+j];}
  double mean=0,square=0,v=0,f=0,xf=0,diff=0,right=0;
  for(int j=0;j<M;j++) {
   double x=path[j],force=4*c4*x*x*x+2*c2*x+c1,d=path[(j+1)%M]-x;
   mean+=x;square+=x*x;v+=c4*x*x*x*x+c2*x*x+c1*x+c0;
   f+=force;xf+=x*force;diff+=d*d;right+=(x>0 ? 1. : (x==0 ? .5 : 0.));
  }
  mean/=M;square/=M;v/=M;f/=M;xf/=M;
  out[7*k]=mean;out[7*k+1]=square;out[7*k+2]=v;
  out[7*k+3]=1/(2*beta)+v+.5*(xf-mean*f);
  out[7*k+4]=M/(2*beta)-M*diff/(2*beta*beta)+v;
  out[7*k+5]=right/M;out[7*k+6]=take;
 }
}
'''
NAMES=('centroid','x2','potential','cv_energy','th_energy','right_occupancy','accepted')

@lru_cache(maxsize=1)
def library():
    flags=['-O3','-std=c11','-shared','-fPIC']
    digest=hashlib.sha256((C_SOURCE+repr(flags)).encode()).hexdigest()
    folder=ROOT/'.cache/polynomial'/digest;folder.mkdir(parents=True,exist_ok=True);path=folder/'kernel.so'
    if not path.exists():
        cc=shutil.which('cc')or shutil.which('gcc')
        if cc is None:raise RuntimeError('Local C compiler required')
        with tempfile.TemporaryDirectory(dir=folder)as tmp:
            src=os.path.join(tmp,'kernel.c');open(src,'w').write(C_SOURCE)
            obj=os.path.join(tmp,'kernel.so')
            subprocess.run([cc,*flags,src,'-o',obj,'-lm'],check=True,capture_output=True)
            os.replace(obj,path)
    lib=ctypes.CDLL(str(path));ptr=ctypes.POINTER(ctypes.c_double)
    lib.polynomial.argtypes=[ctypes.c_long,ptr,ctypes.c_int]+[ctypes.c_double]*6+[ptr]*6
    lib.polynomial.restype=None
    return lib

def potential(x,c):return c[0]*x**4+c[1]*x*x+c[2]*x+c[3]
def force(x,c):return 4*c[0]*x**3+2*c[1]*x+c[2]
def double_well(lam=1.,a=1.,tilt=0.):return (lam,-2*lam*a*a,tilt,lam*a**4)

def action(x,beta,c):
    eps=beta/len(x)
    return float(np.sum((np.roll(x,-1)-x)**2)/(2*eps)+eps*np.sum(potential(x,c)))

def observables(x,beta,c):
    mean=x.mean();v=potential(x,c).mean()
    return np.array([mean,np.mean(x*x),v,1/(2*beta)+v+.5*np.mean((x-mean)*force(x,c)),
        len(x)/(2*beta)-len(x)*np.sum((np.roll(x,-1)-x)**2)/(2*beta*beta)+v,
        np.mean((x>0)+.5*(x==0)),0.])

class PolynomialPIMC:
    def __init__(self,c,beta=2.,beads=16,width=2.,path=None,backend='compiled'):
        self.c=tuple(c);self.beta=beta;self.beads=beads;self.width=width;self.backend=backend
        if len(c)!=4 or c[0]<0 or beta<=0 or beads<3:raise ValueError('Invalid confining polynomial setup')
        self.basis,self.lam,_=normal_modes(beta,beads)
        self.cdf=np.ascontiguousarray(np.cumsum(mode_probabilities(beta,beads,'cv_optimal')));self.cdf[-1]=1.
        self.path=np.ascontiguousarray(np.zeros(beads)if path is None else path,dtype=float)
        if self.path.shape!=(beads,)or not np.isfinite(self.path).all():raise ValueError('Invalid path')
        self.coeff=np.ascontiguousarray(self.basis@self.path)

    def advance(self,u):
        u=np.ascontiguousarray(u,dtype=float)
        if u.ndim!=2 or u.shape[1]!=3 or not np.isfinite(u).all()or np.any(u<0)or np.any(u>=1):raise ValueError('Nx3 finite uniforms in[0,1) required')
        out=np.empty((len(u),7),dtype=float)
        if self.backend=='compiled':
            ptr=lambda a:a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
            library().polynomial(len(u),ptr(u),self.beads,self.beta,*self.c,self.width,
                ptr(self.lam),ptr(self.basis),ptr(self.cdf),ptr(self.path),ptr(self.coeff),ptr(out))
        elif self.backend=='python':
            for k,(us,up,ua)in enumerate(u):
                mode=int(np.searchsorted(self.cdf,us,side='right'));dc=self.width/np.sqrt(self.lam[mode])*(2*up-1)
                candidate=self.path+dc*self.basis[mode]
                da=action(candidate,self.beta,self.c)-action(self.path,self.beta,self.c)
                take=da<=0 or ua==0 or math.log(ua)<-da
                if take:self.path=candidate;self.coeff[mode]+=dc
                out[k]=observables(self.path,self.beta,self.c);out[k,6]=take
        else:raise ValueError('Unknown backend')
        return out

    def state(self):return {'path':self.path.copy(),'coeff':self.coeff.copy()}
    def restore(self,state):
        self.path=np.ascontiguousarray(state['path']);self.coeff=np.ascontiguousarray(state['coeff'])
        if not np.allclose(self.coeff,self.basis@self.path,atol=1e-9,rtol=1e-10):raise ValueError('Inconsistent path checkpoint')
