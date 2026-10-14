"""Independent scalar Gaussian-copula transition prediction; no MCMC fitting."""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import ndtr,ndtri
from scipy.sparse import coo_matrix,eye
from scipy.sparse.linalg import LinearOperator,gmres
from rnglab.v03.physics import mode_probabilities

def prediction(rho,dx=.02,bound=8.,quadrature=40,beta=2.,beads=16):
    width=2.;x=np.linspace(-bound,bound,round(2*bound/dx)+1);n=len(x)
    offsets=np.arange(-round(width/dx),round(width/dx)+1)
    lo=np.maximum(0.,.5+.5*(offsets-.5)*dx/width);hi=np.minimum(1.,.5+.5*(offsets+.5)*dx/width)
    node,weight=leggauss(quadrature);rows=[];cols=[];data=[]
    for i in range(n):
        mask=(i+offsets>=0)&(i+offsets<n)&(offsets!=0);k=offsets[mask];left=lo[mask];right=hi[mask]
        delta=.5*(x[i+k]**2-x[i]**2);alpha=np.exp(-np.maximum(delta,0.))
        if rho==0:probability=(right-left)*alpha
        else:
            proposal=left[:,None]+(right-left)[:,None]*(node+1)[None,:]/2
            threshold=ndtri(np.minimum(alpha,np.nextafter(1.,0.)))[:,None]
            conditional=ndtr((threshold-rho*ndtri(proposal))/np.sqrt(1-rho*rho));conditional[delta<=0,:]=1.
            probability=(right-left)*(conditional@weight)/2
        rows.extend([i]*len(k));cols.extend((i+k).tolist());data.extend(probability.tolist())
    P=coo_matrix((data,(rows,cols)),shape=(n,n)).tocsr();diagonal=1-np.asarray(P.sum(axis=1)).ravel()
    assert diagonal.min()>-1e-12
    P+=coo_matrix((diagonal,(np.arange(n),np.arange(n))),shape=(n,n)).tocsr()
    pi=np.exp(-x*x/2);pi/=pi.sum()
    for iterations in range(10000):
        after=P.T@pi
        if np.sum(abs(after-pi))<2e-14:pi=after;break
        pi=after
    else:raise RuntimeError('Stationarity did not converge')
    p0=float(mode_probabilities(beta,beads,'cv_optimal')[0]);Q=(1-p0)*eye(n,format='csr')+p0*P
    mu=float(pi@x)/np.sqrt(beta);f=x/np.sqrt(beta)-mu
    A=LinearOperator((n,n),matvec=lambda v:v-Q@v+float(pi@v),dtype=float)
    g,info=gmres(A,f,rtol=1e-12,atol=1e-13,restart=100)
    if info:raise RuntimeError('Poisson solver failed:'+str(info))
    sigma2=float(2*np.dot(pi*f,g)-np.dot(pi,f*f))
    return {'rho':rho,'dx':dx,'quadrature':quadrature,'bound':bound,'centroid':mu,
        'standardized_second_moment':float(pi@(x*x)),'asymptotic_variance_centroid':sigma2,
        'zero_mode_selection_probability':p0,'stationarity_residual':float(np.sum(abs(P.T@pi-pi))),
        'poisson_residual':float(np.max(abs(A@g-f))),'row_sum_error':float(np.max(abs(np.asarray(P.sum(axis=1)).ravel()-1))),
        'edge_probability':float(pi[abs(x)>bound-.5].sum()),
        'role':'Unfitted stationary mean and Markov-CLT variance; spatial/quadrature convergence is numerical evidence, not a certified continuum bound'}
