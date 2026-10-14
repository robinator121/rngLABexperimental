"""Synthetic postconditioning numerical-uniform defects, not raw TRNG data."""
import numpy as np
from scipy.special import ndtri,ndtr
from scipy.signal import lfilter

def inject(base,kind='healthy',value=0.):
    u=np.asarray(base,dtype=float).copy()
    if kind=='healthy':return u
    if kind=='copula':
        if not abs(value)<1:raise ValueError('abs(rho)<1 required')
        z=ndtri(np.clip(u[:,1:],2.**-54,np.nextafter(1.,0.)))
        u[:,2]=ndtr(value*z[:,0]+np.sqrt(1-value*value)*z[:,1])
    elif kind=='acceptance_power':
        if value<=0:raise ValueError('kappa>0 required')
        u[:,2]=u[:,2]**(1/value)
    elif kind=='temporal_acceptance':
        if not abs(value)<1:raise ValueError('abs(phi)<1 required')
        z=ndtri(np.clip(u[:,2],2.**-54,np.nextafter(1.,0.)))
        # Stationary first draw; subsequent innovations have variance1-phi².
        noise=z*np.sqrt(1-value*value);noise[0]=z[0]
        u[:,2]=ndtr(lfilter([1.],[1.,-value],noise))
    else:raise ValueError('Unknown defect:'+kind)
    u[:,2]=np.minimum(u[:,2],np.nextafter(1.,0.))
    return u
