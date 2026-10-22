import numpy as np
from scipy.special import expit
from scipy.optimize import minimize
from rnglab14.models.pimc_harmonic import harmonic_energy_estimators,exact_finite_bead_energy

TC=2/np.log(1+np.sqrt(2))

def stationary(P):
    P=np.asarray(P,float)
    if P.ndim!=2 or P.shape[0]!=P.shape[1] or np.any(P<0) or not np.allclose(P.sum(1),1):raise ValueError('stochastic square matrix required')
    if np.linalg.matrix_rank(P.T-np.eye(len(P)),tol=1e-10)!=len(P)-1:raise ValueError('nonunique stationary law')
    A=P.T-np.eye(len(P));A[-1]=1;b=np.zeros(len(P));b[-1]=1
    return np.linalg.solve(A,b)

def two_state(a,b):
    if not(0<a<1 and 0<b<1):raise ValueError('strictly interior ergodic parameters required')
    return np.array([[1-a,a],[b,1-b]])

def variance_mean(a,b,n):
    p=a/(a+b);r=1-a-b;k=np.arange(1,n)
    return p*(1-p)/n*(1+2*np.sum((1-k/n)*r**k))

def augmented(a=.25,b=.5,rho=.8,bins=8):
    # State=(algorithm X, previous driver U). Update driver, then algorithm.
    D=rho*np.eye(bins)+(1-rho)*np.ones((bins,bins))/bins
    P=np.zeros((2*bins,2*bins))
    for x in range(2):
        for old in range(bins):
            for new in range(bins):
                u=(new+.5)/bins;y=(1 if u<a else 0) if x==0 else (0 if u<b else 1)
                P[x*bins+old,y*bins+new]+=D[old,new]
    return P

def energy(spins):
    return -(spins*np.roll(spins,1,axis=-1)+spins*np.roll(spins,1,axis=-2)).sum(axis=(-1,-2))

def enumeration(L,T):
    if L not in (2,3,4):raise ValueError('bounded exact enumeration L2..4')
    n=L*L;s=((np.arange(2**n,dtype=np.uint32)[:,None]>>np.arange(n))&1).astype(np.int8)*2-1
    s=s.reshape(-1,L,L);E=energy(s).astype(float);M=s.sum((1,2)).astype(float)
    w=np.exp(-(E-E.min())/T);w/=w.sum()
    return dict(energy=float(w@E/n),magnetization=float(w@M/n),abs_magnetization=float(w@np.abs(M)/n),
        susceptibility=float((w@(M*M)-(w@M)**2)/(n*T)),heat_capacity=float((w@(E*E)-(w@E)**2)/(n*T*T)))

class Draws:
    """Each replica has an independent SeedSequence child, no shared random arrays."""
    def __init__(self,seed,reps):self.rngs=[np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(reps)];self.reps=reps
    def array(self,shape):return np.stack([g.random(shape) for g in self.rngs])

def ising(L,T,reps,steps,seed,fault='healthy',severity=.15,algorithm='metropolis'):
    draws=Draws(seed,reps);init=draws.array((L,L));s=np.where(init<.5,-1,1).astype(np.int8)
    # Preallocated uniforms permit exact paired healthy/fault comparisons without stream desynchronization.
    U=draws.array((steps,L*L,3));out=np.empty((reps,steps,5));source=np.empty((reps,steps));updates=np.zeros((reps,steps))
    prev=np.full(reps,.5);idx=np.arange(reps)
    for t in range(steps):
        consumed=[]
        for j in range(L*L):
            up=U[:,t,j,0];ua=U[:,t,j,1].copy()
            if fault=='bias':ua=np.clip(ua+severity,0,np.nextafter(1.,0))
            elif fault=='drift':ua=np.clip(ua+severity*max(0,t-steps//4)/(steps*3/4),0,np.nextafter(1.,0))
            elif fault=='abrupt' and t>=steps//3:ua=np.clip(ua+severity,0,np.nextafter(1.,0))
            elif fault=='periodic' and t%40<20:ua=np.clip(ua+severity,0,np.nextafter(1.,0))
            elif fault=='coupled':ua=up.copy()
            elif fault=='correlation':ua=np.where(U[:,t,j,2]<severity,prev,ua);prev=ua.copy()
            site=np.minimum((up*L*L).astype(int),L*L-1);x=site//L;y=site%L
            spin=s[idx,x,y];nn=s[idx,(x+1)%L,y]+s[idx,(x-1)%L,y]+s[idx,x,(y+1)%L]+s[idx,x,(y-1)%L]
            de=2*spin*nn;acc=(de<=0)|(ua<np.exp(-np.maximum(de,0)/T));s[idx,x,y]=np.where(acc,-spin,spin)
            consumed.append(ua);updates[:,t]+=acc
        e=energy(s)/(L*L);m=s.mean((1,2));out[:,t]=np.column_stack([e,m,abs(m),m*m,e*e]);source[:,t]=np.mean(consumed,axis=0)
    return out,source,updates/(L*L)

def wolff(L,T,steps,seed):
    g=np.random.default_rng(seed);s=g.choice([-1,1],size=(L,L));out=[];sizes=[];p=1-np.exp(-2/T)
    for _ in range(steps):
        start=(int(g.integers(L)),int(g.integers(L)));cluster={start};todo=[start];sign=s[start]
        while todo:
            x,y=todo.pop()
            for z in [((x+1)%L,y),((x-1)%L,y),(x,(y+1)%L),(x,(y-1)%L)]:
                if z not in cluster and s[z]==sign and g.random()<p:cluster.add(z);todo.append(z)
        for z in cluster:s[z]*=-1
        e=energy(s)/(L*L);m=s.mean();out.append([e,m,abs(m),m*m,e*e]);sizes.append(len(cluster))
    return np.asarray(out),np.asarray(sizes)

def vmc(reps,steps,seed,fault='healthy',alpha=.5):
    d=Draws(seed,reps);u=d.array((steps,2));x=np.zeros(reps);out=np.empty((reps,steps,3))
    for t in range(steps):
        up=u[:,t,0];ua=u[:,t,1]
        if fault=='bias':up=np.clip(up+.08,0,np.nextafter(1.,0))
        if fault=='coupled':ua=up
        xp=x+2*up-1;lr=-2*alpha*(xp*xp-x*x);x=np.where((lr>=0)|(ua<np.exp(np.minimum(0,lr))),xp,x)
        out[:,t]=np.column_stack([alpha+(.5-2*alpha**2)*x*x,x*x,x**4])
    return out

def pimc(reps,steps,M,seed,fault='healthy',beta=2.):
    d=Draws(seed,reps);x=d.array((M,))-.5;u=d.array((steps,M,2));out=np.empty((reps,steps,4));dt=beta/M
    for t in range(steps):
        for j in range(M):
            up=u[:,t,j,0];ua=u[:,t,j,1]
            if fault=='bias':up=np.clip(up+.08,0,np.nextafter(1.,0))
            if fault=='coupled':ua=up
            old=x[:,j].copy();new=old+.4*(2*up-1);left=x[:,(j-1)%M];right=x[:,(j+1)%M]
            ds=((new-left)**2+(right-new)**2-(old-left)**2-(right-old)**2)/(2*dt)+dt*(new*new-old*old)/2
            x[:,j]=np.where((ds<=0)|(ua<np.exp(np.minimum(0,-ds))),new,old)
        td,cv=harmonic_energy_estimators(x,beta);out[:,t]=np.column_stack([td,cv,(x*x).mean(1),x.mean(1)**2])
    return out

def fit_logistic(X,y):
    if len(np.unique(y))<2:raise ValueError('two label classes required; no forecasting performance claim')
    center=X.mean(0);scale=X.std(0);scale[scale<1e-9]=1
    Z=np.column_stack([np.ones(len(X)),(X-center)/scale])
    def objective(w):
        p=Z@w;return np.mean(np.logaddexp(0,p)-y*p)+.01*np.sum(w[1:]**2)
    opt=minimize(objective,np.zeros(Z.shape[1]),method='BFGS')
    return dict(center=center.tolist(),scale=scale.tolist(),coef=opt.x.tolist(),success=bool(opt.success))

def predict(model,X):
    return expit(np.column_stack([np.ones(len(X)),(X-np.array(model['center']))/np.array(model['scale'])])@model['coef'])
