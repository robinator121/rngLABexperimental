"""Natural-unit VMC and periodic 2D Ising. All randomness is supplied as uniforms."""
from functools import lru_cache
import math
import numpy as np


def vmc_local_energy(x, alpha):
    return alpha/2 + (1-alpha*alpha)*np.asarray(x)**2/2


def vmc_exact(alpha):
    if alpha <= 0:
        raise ValueError("alpha must be positive")
    return (alpha + 1/alpha)/4


def ising_energy(spins):
    return -int(np.sum(spins*np.roll(spins, -1, axis=0) + spins*np.roll(spins, -1, axis=1)))


@lru_cache(maxsize=32)
def ising_exact(length=4, temperature=3.0):
    if length not in (3,4) or temperature <= 0:
        raise ValueError("exact enumeration supports L=3 or 4 and T>0")
    states = ((np.arange(2**(length*length), dtype=np.uint32)[:,None] >> np.arange(length*length)) & 1).astype(np.int8)*2-1
    spins = states.reshape(-1,length,length)
    e = -np.sum(spins*np.roll(spins,-1,axis=1)+spins*np.roll(spins,-1,axis=2),axis=(1,2))
    m = np.mean(states,axis=1)
    weights = np.exp(-(e-e.min())/temperature)
    weights /= weights.sum()
    return {"energy":float(weights @ e / (length*length)), "magnetization":float(weights @ m),
            "abs_magnetization":float(weights @ np.abs(m)), "energy2":float(weights @ (e/(length*length))**2)}


class VMC:
    """psi_alpha(x) proportional to exp(-alpha*x²/2), symmetric box proposals."""
    def __init__(self, alpha=.7, proposal=1.5, x0=0.0):
        if alpha <= 0 or proposal <= 0:
            raise ValueError("alpha and proposal must be positive")
        self.alpha, self.proposal, self.x = alpha, proposal, x0
        self.target = vmc_exact(alpha)

    def advance(self, uniforms):
        n = len(uniforms)
        energy, position, accepted = np.empty(n), np.empty(n), np.zeros(n,dtype=bool)
        x, a, width = self.x, self.alpha, self.proposal
        for i, (up, ua) in enumerate(uniforms):
            proposed = x + width*(2*up-1)
            log_ratio = -a*(proposed*proposed-x*x)
            if log_ratio >= 0 or ua == 0 or math.log(ua) < log_ratio:
                x = proposed
                accepted[i] = True
            energy[i] = a/2+(1-a*a)*x*x/2
            position[i] = x
        self.x = x
        return energy, {"x":position,"x2":position**2,"accepted":accepted}


class Ising:
    """Random-site single-spin Metropolis; one step is ONE attempted flip, not a sweep."""
    def __init__(self, length=4, temperature=3., initial="ordered"):
        if length < 3 or temperature <= 0 or initial not in {"ordered","checkerboard"}:
            raise ValueError("Ising requires L>=3, T>0, ordered/checkerboard start")
        self.length, self.temperature = length, temperature
        self.spins = np.ones((length,length),dtype=np.int8)
        if initial == "checkerboard":
            self.spins = (1-2*(np.indices((length,length)).sum(axis=0)%2)).astype(np.int8)
        self.energy = ising_energy(self.spins)
        self.magnetization = int(self.spins.sum())
        self.target = ising_exact(length, temperature)["energy"] if length in (3,4) else None

    def advance(self, uniforms):
        n, L = len(uniforms), self.length
        energy, magnetization, accepted = np.empty(n), np.empty(n), np.zeros(n,dtype=bool)
        s, e, m = self.spins, self.energy, self.magnetization
        for k,(us,ua) in enumerate(uniforms):
            site = int(us*L*L)
            i,j = divmod(site,L)
            spin = int(s[i,j])
            delta = 2*spin*(int(s[(i-1)%L,j])+int(s[(i+1)%L,j])+int(s[i,(j-1)%L])+int(s[i,(j+1)%L]))
            if delta <= 0 or ua < math.exp(-delta/self.temperature):
                s[i,j] = -spin
                e += delta
                m -= 2*spin
                accepted[k] = True
            energy[k],magnetization[k] = e/(L*L),m/(L*L)
        self.energy,self.magnetization=e,m
        return energy,{"magnetization":magnetization,"abs_magnetization":np.abs(magnetization),"accepted":accepted}


class Independent:
    def __init__(self, kind="integral"):
        self.kind = kind
        self.target = 1/3 if kind == "integral" else math.pi

    def advance(self, uniforms):
        value = uniforms[:,0]**2 if self.kind == "integral" else 4*(np.sum(uniforms**2,axis=1)<=1)
        return value, {}


def make_model(config):
    kind = config["kind"]
    if kind == "vmc":
        return VMC(**{k:config[k] for k in ["alpha","proposal","x0"] if k in config})
    if kind == "ising":
        return Ising(**{k:config[k] for k in ["length","temperature","initial"] if k in config})
    if kind in {"integral","pi"}:
        return Independent(kind)
    raise ValueError(f"Unknown physical workload {kind}")
