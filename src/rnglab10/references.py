"""Independent Schrödinger and finite-M transfer references, without MCMC."""
import numpy as np
from scipy.linalg import eigh,eigh_tridiagonal
from .polynomial import potential,force

def hamiltonian_basis(c,size=96,beta=2.):
    n=size+6
    x=np.diag(np.sqrt(np.arange(1,n)/2),1)+np.diag(np.sqrt(np.arange(1,n)/2),-1)
    x2=x@x;x4=x2@x2
    H=np.diag(np.arange(size)+.5)+(c[1]-.5)*x2[:size,:size]+c[0]*x4[:size,:size]+c[2]*x[:size,:size]+c[3]*np.eye(size)
    energies,vec=eigh(H);w=np.exp(-beta*(energies-energies[0]));w/=w.sum()
    def avg(matrix):return float(w@np.diag(vec.T@matrix[:size,:size]@vec))
    return {'method':'Padded harmonic oscillator matrix basis','size':size,'energies':energies[:8].tolist(),
        'energy':float(w@energies),'centroid':avg(x),'x2':avg(x2),'x4':avg(x4),
        'potential':c[0]*avg(x4)+c[1]*avg(x2)+c[2]*avg(x)+c[3],
        'logZ':float(-beta*energies[0]+np.log(np.exp(-beta*(energies-energies[0])).sum())),
        'highest_thermal_weight':float(w[-1])}

def hamiltonian_fd(c,dx=.01,bound=6.,beta=2.,states=48):
    n=round(2*bound/dx);x=np.linspace(-bound,bound,n+1)[1:-1];dx=x[1]-x[0]
    diag=1/dx**2+potential(x,c);off=np.full(len(x)-1,-.5/dx**2)
    energies,vec=eigh_tridiagonal(diag,off,select='i',select_range=(0,states-1))
    w=np.exp(-beta*(energies-energies[0]));w/=w.sum();pi=(vec*vec)@w
    return {'method':'Independent second-order finite-difference Schrodinger operator','dx':float(dx),'bound':bound,
        'energies':energies[:8].tolist(),'energy':float(w@energies),'centroid':float(pi@x),
        'x2':float(pi@(x*x)),'potential':float(pi@potential(x,c)),
        'right_occupancy':float(pi@((x>0)+.5*(x==0))),
        'logZ':float(-beta*energies[0]+np.log(np.exp(-beta*(energies-energies[0])).sum())),
        'highest_thermal_weight':float(w[-1]),'edge_probability':float(pi[abs(x)>bound-.5].sum())}, {'grid':x,'density':pi/dx,'ground_density':vec[:,0]**2/dx}

def transfer(c,beads=16,beta=2.,dx=.04,bound=6.,kappa=1.):
    n=round(2*bound/dx);x=np.linspace(-bound,bound,n+1);dx=x[1]-x[0]
    weights=np.full(len(x),dx);weights[[0,-1]]/=2
    eps=beta/beads;v=potential(x,c);f=force(x,c);d2=(x[:,None]-x[None,:])**2
    K=np.sqrt(kappa*weights[:,None]*weights[None,:]/(2*np.pi*eps))*np.exp(-kappa*d2/(2*eps)-kappa*eps*(v[:,None]+v[None,:])/2)
    eig,vec=eigh(K);negative=float(abs(eig[eig<0]).sum());eig=np.maximum(eig,0);scale=eig[-1];eig/=scale
    w=eig**beads;Z=w.sum();pi=(vec*vec)@(w/Z)
    previous=(vec*eig[None,:]**(beads-1))@vec.T
    bond=previous*(K/scale)/Z
    thermodynamic=1/(2*eps)-float(np.sum(bond*d2))/(2*eps*eps)+pi@v
    X=vec.T@(x[:,None]*vec);F=vec.T@(f[:,None]*vec)
    cf=0.;centroid2=0.
    for lag in range(beads):
        factor=eig[:,None]**(beads-lag)*eig[None,:]**lag
        cf+=np.sum(X*F.T*factor)/Z/beads
        centroid2+=np.sum(X*X.T*factor)/Z/beads
    cv=1/(2*beta)+pi@v+.5*(pi@(x*f)-cf)
    predicted_gap=(beads-1)/(2*beta)*(1-1/kappa)
    result={'method':'Symmetric primitive Nystrom transfer matrix','beads':beads,'beta':beta,'dx':float(dx),'bound':bound,'kappa':kappa,
        'centroid':float(pi@x),'centroid2':float(centroid2),'x2':float(pi@(x*x)),
        'potential':float(pi@v),'cv_energy':float(cv),'th_energy':float(thermodynamic),
        'right_occupancy':float(pi@((x>0)+.5*(x==0))),
        'energy_gap':float(thermodynamic-cv),'predicted_gap':predicted_gap,
        'gap_identity_residual':float(thermodynamic-cv-predicted_gap),
        'logZ':float(beads*np.log(scale)+np.log(Z)),
        'negative_eigenvalue_mass_roundoff':negative,'edge_probability':float(pi[abs(x)>bound-.5].sum())}
    return result,{'grid':x,'marginal':pi}
