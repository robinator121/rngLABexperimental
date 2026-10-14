"""Supplied-draw physical kernels for V0.3; historical V0.2 is unchanged.

One row is one attempted proposal (not an Ising or PIMC sweep). Compilation is
optional: the Python implementation is the executable reference. No kernel
obtains random numbers internally. Checkpoints preserve the whole chain state.
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
import warnings

import numpy as np
from scipy.special import ndtri

from rnglab.v02.physics import ising_energy, ising_exact, vmc_exact, vmc_local_energy
from rnglab.v02.pimc import action, finite_energy, thermal_energy

VERSION = "v03-physical-kernels-3"

# -ffast-math is intentionally excluded: acceptance and supplied-draw semantics
# must match the scalar Python reference, including log(0) handling.
_C_SOURCE = r"""
#include <math.h>
#include <stddef.h>

void vmc(long n, const double *u, double alpha, double width, double *state,
         double *energy, double *xtrace, double *accepted) {
    double x=state[0];
    for(long k=0;k<n;k++) {
        double proposed=x+width*(2*u[2*k]-1);
        double ratio=-alpha*(proposed*proposed-x*x);
        if(ratio>=0 || u[2*k+1]==0 || log(u[2*k+1])<ratio) {
            x=proposed; accepted[k]=1;
        } else accepted[k]=0;
        energy[k]=alpha/2+(1-alpha*alpha)*x*x/2;
        xtrace[k]=x;
    }
    state[0]=x;
}

void ising(long n, const double *u, int L, double T, int *s, long *state,
           double *energy, double *mag, double *accepted) {
    long e=state[0],m=state[1]; int sites=L*L;
    for(long k=0;k<n;k++) {
        int site=(int)(u[2*k]*sites),i=site/L,j=site%L,spin=s[site];
        int delta=2*spin*(s[((i+L-1)%L)*L+j]+s[((i+1)%L)*L+j]
                         +s[i*L+(j+L-1)%L]+s[i*L+(j+1)%L]);
        if(delta<=0 || u[2*k+1]<exp(-delta/T)) {
            s[site]=-spin; e+=delta;m-=2*spin;accepted[k]=1;
        } else accepted[k]=0;
        energy[k]=(double)e/sites;mag[k]=(double)m/sites;
    }
    state[0]=e;state[1]=m;
}

void pimc(long n, const double *u, int M, double beta, double cp,
          double bw, double cw, double *p, double *state, double *cv,
          double *thermo, double *qtrace, double *ctrace, double *accepted,
          double *centroid_move) {
    double eps=beta/M,s=state[0],q=state[1],diff=state[2];
    for(long k=0;k<n;k++) {
        double us=u[3*k],up=u[3*k+1],ua=u[3*k+2];
        int centroid=(cp==1 || (cp>0 && us>=1-cp));
        int site=centroid ? M : (int)(us/(1-cp)*M);
        if(site>=M && !centroid) site=M-1;
        double delta=(centroid ? cw:bw)*(2*up-1),dq,ddiff,ds;
        if(centroid) {
            dq=2*delta*s+M*delta*delta;ddiff=0;ds=M*delta;
        } else {
            double old=p[site],proposed=old+delta;
            double left=p[(site+M-1)%M],right=p[(site+1)%M];
            dq=proposed*proposed-old*old;ds=delta;
            ddiff=(proposed-left)*(proposed-left)+(right-proposed)*(right-proposed)
                 -(old-left)*(old-left)-(right-old)*(right-old);
        }
        double da=ddiff/(2*eps)+eps*dq/2;
        if(da<=0 || ua==0 || log(ua)<-da) {
            if(centroid) for(int j=0;j<M;j++) p[j]+=delta;
            else p[site]+=delta;
            s+=ds;q+=dq;diff+=ddiff;accepted[k]=1;
        } else accepted[k]=0;
        cv[k]=1/(2*beta)+q/M-.5*(s/M)*(s/M);
        thermo[k]=M/(2*beta)-M*diff/(2*beta*beta)+q/(2*M);
        qtrace[k]=q/M;ctrace[k]=s/M;centroid_move[k]=centroid;
    }
    state[0]=s;state[1]=q;state[2]=diff;
}

void pimc_modes(long n, const double *u, int M, double beta, double width,
                const double *lambda, const double *ringlambda, const double *cdf, double *c,
                double *state, double *cv, double *thermo, double *qtrace,
                double *ctrace, double *accepted, double *centroid_move) {
    double q=state[1],diff=state[2],s=state[0];
    for(long k=0;k<n;k++) {
        int mode=0;
        if(cdf==NULL) mode=(int)(u[3*k]*M);
        else while(mode<M-1 && u[3*k]>=cdf[mode]) mode++;
        double proposed=c[mode]+width/sqrt(lambda[mode])*(2*u[3*k+1]-1);
        double dq=proposed*proposed-c[mode]*c[mode];
        double da=lambda[mode]*dq/2;
        if(da<=0 || u[3*k+2]==0 || log(u[3*k+2])<-da) {
            c[mode]=proposed;q+=dq;diff+=ringlambda[mode]*dq;accepted[k]=1;
        } else accepted[k]=0;
        s=sqrt((double)M)*c[0];
        cv[k]=1/(2*beta)+q/M-.5*(s/M)*(s/M);
        thermo[k]=M/(2*beta)-M*diff/(2*beta*beta)+q/(2*M);
        qtrace[k]=q/M;ctrace[k]=s/M;centroid_move[k]=(mode==0);
    }
    state[0]=s;state[1]=q;state[2]=diff;
}
"""
KERNEL_SOURCE_SHA256 = hashlib.sha256(_C_SOURCE.encode()).hexdigest()


@lru_cache(maxsize=1)
def _compiled_library():
    """Compile a small local shared object once; use a hashed, atomic cache."""
    compiler = shutil.which("cc") or shutil.which("gcc")
    if compiler is None:
        return None
    cache_root = os.environ.get("RNGLAB_KERNEL_CACHE")
    cache = (Path(cache_root) if cache_root else Path(__file__).resolve().parents[3] / ".cache" / "kernels") / KERNEL_SOURCE_SHA256
    cache.mkdir(parents=True, exist_ok=True)
    library = cache / "kernels.so"
    if not library.exists():
        with tempfile.TemporaryDirectory(prefix="build-", dir=cache) as staging:
            src = Path(staging) / "kernels.c"
            so = Path(staging) / "kernels.so"
            src.write_text(_C_SOURCE)
            try:
                subprocess.run([compiler, "-O3", "-std=c99", "-fPIC", "-shared",
                                str(src), "-o", str(so), "-lm"], check=True,
                               capture_output=True, timeout=30)
                os.replace(so, library)
            except (OSError, subprocess.SubprocessError):
                return None
    lib = ctypes.CDLL(str(library))
    dp, ip, lp = (ctypes.POINTER(t) for t in
                  (ctypes.c_double, ctypes.c_int, ctypes.c_long))
    lib.vmc.argtypes = [ctypes.c_long, dp, ctypes.c_double, ctypes.c_double, dp, dp, dp, dp]
    lib.ising.argtypes = [ctypes.c_long, dp, ctypes.c_int, ctypes.c_double, ip, lp, dp, dp, dp]
    lib.pimc.argtypes = [ctypes.c_long, dp, ctypes.c_int, ctypes.c_double,
                        ctypes.c_double, ctypes.c_double, ctypes.c_double,
                        dp, dp, dp, dp, dp, dp, dp, dp]
    lib.pimc_modes.argtypes = [ctypes.c_long, dp, ctypes.c_int, ctypes.c_double,
                              ctypes.c_double, dp, dp, dp, dp, dp, dp, dp, dp, dp, dp, dp]
    lib.vmc.restype = lib.ising.restype = lib.pimc.restype = lib.pimc_modes.restype = None
    return lib


def _pointer(a, ctype=ctypes.c_double):
    return a.ctypes.data_as(ctypes.POINTER(ctype))


def _select_backend(backend):
    backend = os.environ.get("RNGLAB_KERNEL_BACKEND", backend)
    if backend not in {"auto", "python", "compiled"}:
        raise ValueError("backend must be auto, python, or compiled")
    if backend == "python":
        return None, "python"
    lib = _compiled_library()
    if lib is None and backend == "compiled":
        raise RuntimeError("compiled kernels requested but a local C compiler is unavailable")
    return lib, "compiled" if lib is not None else "python"


class Model:
    dimension: int
    kind: str
    initialization_uniform_count = 0
    mixing_note = "Independent integration; temporal input dependence can inflate variance."

    def _uniforms(self, u):
        a = np.ascontiguousarray(u, dtype=np.float64)
        if a.ndim != 2 or a.shape[1] != self.dimension:
            raise ValueError(f"{self.kind} requires an (n,{self.dimension}) uniform array")
        if not np.all(np.isfinite(a)) or np.any(a < 0) or np.any(a >= 1):
            raise ValueError("uniforms must be finite and in [0,1)")
        return a

    def state_dict(self):
        return {"version": VERSION, "kind": self.kind, "state": self._state()}

    def _initialization_uniforms(self, uniforms):
        u = np.asarray(uniforms, dtype=np.float64)
        if u.ndim != 1 or len(u) != self.initialization_uniform_count:
            raise ValueError(f"{self.kind} initialization requires {self.initialization_uniform_count} supplied uniforms")
        if not np.isfinite(u).all() or np.any(u < 0) or np.any(u >= 1):
            raise ValueError("initialization uniforms must be finite and in [0,1)")
        return u

    def initialize_from_uniforms(self, uniforms):
        self._initialization_uniforms(uniforms)
        if self.initialization_uniform_count:
            raise NotImplementedError("equilibrium initialization unavailable for this workload")
        return self.state_dict()

    def load_state_dict(self, checkpoint):
        if checkpoint["version"] != VERSION or checkpoint["kind"] != self.kind:
            raise ValueError("checkpoint workload or version differs")
        self._restore(checkpoint["state"])


class Independent(Model):
    dimension, kind = 1, "integral"

    def __init__(self, integrand="x2", backend="auto"):
        aliases = {"square": "x2", "indicator_half": "indicator",
                   "oscillatory": "cos", "cosine": "cos"}
        self.integrand = aliases.get(integrand, integrand)
        if self.integrand not in {"x2", "indicator", "cos"}:
            raise ValueError("integrand must be x2, indicator, or cos")
        self.target = {"x2": 1/3, "indicator": .5, "cos": 0.}[self.integrand]
        self.targets = {"value": self.target}
        self.backend = "numpy"

    def advance(self, uniforms):
        x = self._uniforms(uniforms)[:, 0]
        y = (x*x if self.integrand == "x2" else
             (x < .5).astype(float) if self.integrand == "indicator" else
             np.cos(2*np.pi*x))
        return y, {}

    def _state(self):
        return {"integrand": self.integrand}

    def _restore(self, state):
        if state["integrand"] != self.integrand:
            raise ValueError("integrand differs from checkpoint")


class VMC(Model):
    dimension, kind = 2, "vmc"
    initialization_uniform_count = 1
    mixing_note = "Stationary density exp(-alpha*x²); check burn-in, proposal scale, and ESS."

    def __init__(self, alpha=.7, proposal=1.5, x0=0., observable="energy", backend="auto"):
        if not np.isfinite([alpha, proposal, x0]).all() or alpha <= 0 or proposal <= 0:
            raise ValueError("finite alpha>0, proposal>0, and x0 required")
        self.alpha, self.proposal, self.x = float(alpha), float(proposal), float(x0)
        if observable not in {"energy", "x", "x2"}:
            raise ValueError("VMC observable must be energy, x, or x2")
        self.observable = observable
        self.target = vmc_exact(alpha)
        self.targets = {"energy": self.target, "x": 0., "x2": 1/(2*alpha),
                        "potential": 1/(4*alpha), "kinetic": alpha/4}
        self.target = self.targets[observable]
        self._lib, self.backend = _select_backend(backend)

    def advance(self, uniforms):
        u = self._uniforms(uniforms)
        e, position, accepted = (np.zeros(len(u)) for _ in range(3))
        if self._lib is not None:
            state = np.array([self.x])
            self._lib.vmc(len(u), _pointer(u), self.alpha, self.proposal,
                          _pointer(state), _pointer(e), _pointer(position), _pointer(accepted))
            self.x = float(state[0])
        else:
            x, a, width = self.x, self.alpha, self.proposal
            for k, (up, ua) in enumerate(u):
                proposed = x + width*(2*up-1)
                ratio = -a*(proposed*proposed-x*x)
                if ratio >= 0 or ua == 0 or math.log(ua) < ratio:
                    x = proposed
                    accepted[k] = 1
                e[k] = a/2+(1-a*a)*x*x/2
                position[k] = x
            self.x = float(x)
        aux = {"energy": e, "x": position, "x2": position**2, "accepted": accepted.astype(bool),
                   "potential": .5*position**2,
                   "kinetic": self.alpha/2-.5*self.alpha**2*position**2}
        return aux[self.observable], aux

    def initialize_from_uniforms(self, uniforms):
        u = self._initialization_uniforms(uniforms)
        # Zero is a continuous measure-zero endpoint; nextafter prevents an
        # infinite initial state for the finite uniform grid's endpoint draw.
        self.x = float(ndtri(max(u[0], np.nextafter(0., 1.)))/math.sqrt(2*self.alpha))
        return self.state_dict()

    def _state(self):
        return {"alpha": self.alpha, "proposal": self.proposal, "observable": self.observable, "x": self.x}

    def _restore(self, state):
        if (state["alpha"] != self.alpha or state["proposal"] != self.proposal
                or state["observable"] != self.observable):
            raise ValueError("VMC checkpoint parameters differ")
        if not np.isfinite(state["x"]):
            raise ValueError("nonfinite VMC checkpoint")
        self.x = float(state["x"])


@lru_cache(maxsize=16)
def _ising_equilibrium_distribution(length, temperature):
    """Exact finite Boltzmann categorical CDF, with ground-energy stabilization."""
    sites = length*length
    states = ((np.arange(2**sites, dtype=np.uint32)[:, None] >> np.arange(sites)) & 1).astype(np.int8)*2-1
    spins = states.reshape(-1, length, length)
    energy = -np.sum(spins*np.roll(spins, -1, axis=1)+spins*np.roll(spins, -1, axis=2), axis=(1, 2))
    probabilities = np.exp(-(energy-energy.min())/temperature)
    probabilities /= probabilities.sum()
    cdf = np.cumsum(probabilities)
    cdf[-1] = 1.
    states.setflags(write=False)
    cdf.setflags(write=False)
    return states, cdf


class Ising(Model):
    dimension, kind = 2, "ising"
    initialization_uniform_count = 1
    mixing_note = "One attempted spin flip per row; critical slowing and spin-sector mixing need separate checks."

    def __init__(self, length=4, temperature=3., initial="ordered", observable="energy", backend="auto"):
        if length not in {3, 4} or not np.isfinite(temperature) or temperature <= 0:
            raise ValueError("V0.3 exact-reference Ising supports L=3 or 4 and T>0")
        if initial not in {"ordered", "checkerboard"}:
            raise ValueError("Ising initial must be ordered or checkerboard")
        if observable not in {"energy", "magnetization", "abs_magnetization"}:
            raise ValueError("Ising observable must be energy, magnetization, or abs_magnetization")
        self.observable = observable
        self.length, self.temperature = int(length), float(temperature)
        self.spins = np.ones((length, length), dtype=np.int32)
        if initial == "checkerboard":
            self.spins = (1-2*(np.indices((length, length)).sum(axis=0) % 2)).astype(np.int32)
        self.energy, self.magnetization = ising_energy(self.spins), int(self.spins.sum())
        self.targets = ising_exact(length, temperature).copy()
        self.target = self.targets[observable]
        self._lib, self.backend = _select_backend(backend)

    def advance(self, uniforms):
        u = self._uniforms(uniforms)
        e, mag, accepted = (np.zeros(len(u)) for _ in range(3))
        if self._lib is not None:
            state = np.array([self.energy, self.magnetization], dtype=np.dtype(ctypes.c_long))
            self._lib.ising(len(u), _pointer(u), self.length, self.temperature,
                            _pointer(self.spins, ctypes.c_int), _pointer(state, ctypes.c_long),
                            _pointer(e), _pointer(mag), _pointer(accepted))
            self.energy, self.magnetization = (int(v) for v in state)
        else:
            s, total_e, total_m, L = self.spins, self.energy, self.magnetization, self.length
            for k, (us, ua) in enumerate(u):
                i, j = divmod(int(us*L*L), L)
                spin = int(s[i, j])
                delta = 2*spin*(int(s[(i-1) % L, j])+int(s[(i+1) % L, j])
                                +int(s[i, (j-1) % L])+int(s[i, (j+1) % L]))
                if delta <= 0 or ua < math.exp(-delta/self.temperature):
                    s[i, j] = -spin
                    total_e += delta
                    total_m -= 2*spin
                    accepted[k] = 1
                e[k], mag[k] = total_e/(L*L), total_m/(L*L)
            self.energy, self.magnetization = int(total_e), int(total_m)
        aux = {"energy": e, "magnetization": mag, "abs_magnetization": np.abs(mag),
                   "accepted": accepted.astype(bool)}
        return aux[self.observable], aux

    def initialize_from_uniforms(self, uniforms):
        u = self._initialization_uniforms(uniforms)
        states, cdf = _ising_equilibrium_distribution(self.length, self.temperature)
        selected = min(len(cdf)-1, int(np.searchsorted(cdf, u[0], side="right")))
        self.spins = states[selected].reshape(self.length, self.length).astype(np.int32)
        self.energy, self.magnetization = ising_energy(self.spins), int(self.spins.sum())
        return self.state_dict()

    def _state(self):
        return {"length": self.length, "temperature": self.temperature,
                "observable": self.observable,
                "spins": self.spins.tolist(), "energy": self.energy,
                "magnetization": self.magnetization}

    def _restore(self, state):
        if (state["length"] != self.length or state["temperature"] != self.temperature
                or state["observable"] != self.observable):
            raise ValueError("Ising checkpoint parameters differ")
        s = np.asarray(state["spins"], dtype=np.int32)
        if s.shape != (self.length, self.length) or not np.all(np.isin(s, [-1, 1])):
            raise ValueError("invalid spin checkpoint")
        if ising_energy(s) != state["energy"] or s.sum() != state["magnetization"]:
            raise ValueError("Ising checkpoint bookkeeping is inconsistent")
        self.spins = s.copy()
        self.energy, self.magnetization = int(state["energy"]), int(state["magnetization"])


@lru_cache(maxsize=64)
def finite_covariance(beta, beads):
    """Exact covariance of exp(-action), independent matrix reference."""
    if beta <= 0 or beads < 3:
        raise ValueError("beta>0 and beads>=3 required")
    eps = beta/beads
    identity = np.eye(beads)
    precision = ((eps+2/eps)*identity
                 -(np.roll(identity, 1, axis=0)+np.roll(identity, -1, axis=0))/eps)
    covariance = np.linalg.inv(precision)
    covariance.setflags(write=False)
    return covariance


@lru_cache(maxsize=64)
def normal_modes(beta, beads):
    """Real orthonormal Fourier basis and exact precision/ring eigenvalues.

    Rows are mode vectors. Zero, cosine/sine pairs, then the even-M Nyquist
    vector form a full real basis. A uniform random mode plus a symmetric box
    displacement defines a reversible Metropolis proposal with no Gaussian RNG.
    """
    j = np.arange(beads)
    vectors = [np.full(beads, 1/math.sqrt(beads))]
    frequencies = [0]
    for frequency in range(1, (beads+1)//2):
        phase = 2*np.pi*frequency*j/beads
        vectors.extend([math.sqrt(2/beads)*np.cos(phase), math.sqrt(2/beads)*np.sin(phase)])
        frequencies.extend([frequency, frequency])
    if beads % 2 == 0:
        vectors.append((-1.)**j/math.sqrt(beads))
        frequencies.append(beads//2)
    basis = np.ascontiguousarray(vectors)
    ring = 4*np.sin(np.pi*np.asarray(frequencies)/beads)**2
    precision = beta/beads+ring/(beta/beads)
    for a in (basis, ring, precision):
        a.setflags(write=False)
    return basis, precision, ring


def mode_probabilities(beta, beads, selection="uniform"):
    """State-independent random-scan probabilities, derived before measurement.

    For E_cv=constant+sum(w_j*c_j²), w_0=1/(2M), w_other=1/M.
    At the exact Gaussian target Var(w_j*c_j²)=2*w_j²/lambda_j².
    Every standardized box proposal has the same per-update transition kernel.
    Random-scan thinning therefore gives tau_j=(tau_update+1)/p_j-1.
    Since stationary cross-mode covariances vanish, the minimizing scan weights
    are proportional to sqrt(Var_j)=abs(w_j)/lambda_j. All remain positive.
    """
    if selection not in {"uniform", "cv_optimal"}:
        raise ValueError("mode_selection must be uniform or cv_optimal")
    if selection == "uniform":
        return np.full(beads, 1/beads)
    _, precision, _ = normal_modes(beta, beads)
    weights = np.full(beads, 1/beads)
    weights[0] /= 2
    probabilities = weights/precision
    probabilities /= probabilities.sum()
    return probabilities


class PIMC(Model):
    dimension, kind = 3, "pimc"
    mixing_note = "Local bead moves slow at large bead counts; centroid moves sample only the zero mode. Check long-wavelength ESS and initial-state sensitivity."

    def __init__(self, beta=2., beads=16, centroid_probability=.2,
                 bead_proposal=None, centroid_proposal=1., estimator="centroid_virial",
                 path0=0., sampler="local_centroid", mode_proposal=2., mode_selection="uniform", backend="auto"):
        if not np.isfinite([beta, centroid_probability, centroid_proposal]).all():
            raise ValueError("PIMC parameters must be finite")
        if beta <= 0 or int(beads) != beads or beads < 3 or not 0 <= centroid_probability <= 1:
            raise ValueError("beta>0, integer beads>=3, centroid_probability in [0,1] required")
        if estimator not in {"centroid_virial", "thermodynamic"}:
            raise ValueError("energy estimator must be centroid_virial or thermodynamic")
        if sampler not in {"local_centroid", "normal_mode"}:
            raise ValueError("sampler must be local_centroid or normal_mode")
        if not np.isfinite(mode_proposal) or mode_proposal <= 0:
            raise ValueError("mode_proposal must be finite and positive")
        self.beta, self.beads = float(beta), int(beads)
        self.initialization_uniform_count = self.beads
        self.centroid_probability = float(centroid_probability)
        self.bead_proposal = 2*math.sqrt(beta/beads) if bead_proposal is None else float(bead_proposal)
        self.centroid_proposal, self.estimator = float(centroid_proposal), estimator
        self.sampler, self.mode_proposal = sampler, float(mode_proposal)
        self.mode_selection = mode_selection
        if (not np.isfinite(self.bead_proposal) or self.bead_proposal <= 0
                or self.centroid_proposal <= 0):
            raise ValueError("positive finite proposal widths required")
        self.path = np.full(self.beads, path0, dtype=np.float64) if np.isscalar(path0) else np.asarray(path0, dtype=np.float64).copy()
        if self.path.shape != (self.beads,) or not np.isfinite(self.path).all():
            raise ValueError("finite scalar or beads-long initial path required")
        self.sum, self.sumsq, self.sumdiff = self._path_sums()
        self._basis, self._precision, self._ring = normal_modes(beta, self.beads)
        self.mode_probabilities = mode_probabilities(beta, self.beads, mode_selection)
        self._mode_cdf = np.cumsum(self.mode_probabilities)
        self._mode_cdf[-1] = 1.
        self.mode_coefficients = np.ascontiguousarray(self._basis @ self.path)
        self.target, self.continuum_target = finite_energy(beta, beads), thermal_energy(beta)
        # beta * <centroid²> = 1 by the zero precision eigenvalue being epsilon.
        qtarget = self.target
        # E_cv = 1/(2 beta) + <q/M> - <centroid²>/2 = <q/M>.
        self.targets = {"energy": self.target, "centroid_virial": self.target,
                        "thermodynamic": self.target, "x2": qtarget,
                        "centroid": 0., "centroid2": 1/beta,
                        "potential": qtarget/2, "kinetic": self.target-qtarget/2}
        self._lib, self.backend = _select_backend(backend)
        if centroid_probability == 1 and sampler == "local_centroid":
            warnings.warn("centroid-only moves are nonergodic for internal path modes", RuntimeWarning)

    def _path_sums(self):
        return (float(self.path.sum()), float(np.sum(self.path**2)),
                float(np.sum((np.roll(self.path, -1)-self.path)**2)))

    def initialize_from_uniforms(self, uniforms):
        """Exact finite-ring Gaussian initialization from M supplied clean draws."""
        u = self._initialization_uniforms(uniforms)
        self.mode_coefficients = np.ascontiguousarray(ndtri(np.maximum(u, np.nextafter(0., 1.)))/np.sqrt(self._precision))
        self.path = np.ascontiguousarray(self.mode_coefficients @ self._basis)
        self.sum, self.sumsq, self.sumdiff = self._path_sums()
        return self.state_dict()

    def advance(self, uniforms):
        u = self._uniforms(uniforms)
        cv, thermo, qtrace, ctrace, accepted, centroid_move = (np.zeros(len(u)) for _ in range(6))
        M, beta, cp = self.beads, self.beta, self.centroid_probability
        eps = beta/M
        if self.sampler == "normal_mode":
            state = np.array([self.sum, self.sumsq, self.sumdiff])
            if self._lib is not None:
                self._lib.pimc_modes(len(u), _pointer(u), M, beta, self.mode_proposal,
                                     _pointer(self._precision), _pointer(self._ring),
                                     None if self.mode_selection == "uniform" else _pointer(self._mode_cdf),
                                     _pointer(self.mode_coefficients), _pointer(state),
                                     _pointer(cv), _pointer(thermo), _pointer(qtrace),
                                     _pointer(ctrace), _pointer(accepted), _pointer(centroid_move))
            else:
                s, q, diff = state
                c = self.mode_coefficients
                for k, (us, up, ua) in enumerate(u):
                    mode = int(us*M) if self.mode_selection == "uniform" else int(np.searchsorted(self._mode_cdf, us, side="right"))
                    proposed = c[mode]+self.mode_proposal/math.sqrt(self._precision[mode])*(2*up-1)
                    dq = proposed*proposed-c[mode]*c[mode]
                    da = self._precision[mode]*dq/2
                    if da <= 0 or ua == 0 or math.log(ua) < -da:
                        c[mode] = proposed
                        q += dq
                        diff += self._ring[mode]*dq
                        accepted[k] = 1
                    s = math.sqrt(M)*c[0]
                    cv[k] = 1/(2*beta)+q/M-.5*(s/M)*(s/M)
                    thermo[k] = M/(2*beta)-M*diff/(2*beta*beta)+q/(2*M)
                    qtrace[k], ctrace[k], centroid_move[k] = q/M, s/M, mode == 0
                state = [s, q, diff]
            self.sum, self.sumsq, self.sumdiff = (float(v) for v in state)
            self.path = np.ascontiguousarray(self.mode_coefficients @ self._basis)
        elif self._lib is not None:
            state = np.array([self.sum, self.sumsq, self.sumdiff])
            self._lib.pimc(len(u), _pointer(u), M, beta, cp,
                           self.bead_proposal, self.centroid_proposal,
                           _pointer(self.path), _pointer(state), _pointer(cv), _pointer(thermo),
                           _pointer(qtrace), _pointer(ctrace), _pointer(accepted), _pointer(centroid_move))
            self.sum, self.sumsq, self.sumdiff = (float(v) for v in state)
        else:
            s, q, diff, p = self.sum, self.sumsq, self.sumdiff, self.path
            for k, (us, up, ua) in enumerate(u):
                centroid = cp == 1 or (cp > 0 and us >= 1-cp)
                site = M if centroid else min(M-1, int(us/(1-cp)*M))
                delta = (self.centroid_proposal if centroid else self.bead_proposal)*(2*up-1)
                if centroid:
                    dq, ddiff, ds = 2*delta*s+M*delta*delta, 0., M*delta
                else:
                    old, proposed = p[site], p[site]+delta
                    left, right = p[(site-1) % M], p[(site+1) % M]
                    dq, ds = proposed*proposed-old*old, delta
                    ddiff = ((proposed-left)*(proposed-left)+(right-proposed)*(right-proposed)
                             -(old-left)*(old-left)-(right-old)*(right-old))
                da = ddiff/(2*eps)+eps*dq/2
                if da <= 0 or ua == 0 or math.log(ua) < -da:
                    if centroid:
                        p += delta
                    else:
                        p[site] += delta
                    s += ds
                    q += dq
                    diff += ddiff
                    accepted[k] = 1
                cv[k] = 1/(2*beta)+q/M-.5*(s/M)*(s/M)
                thermo[k] = M/(2*beta)-M*diff/(2*beta*beta)+q/(2*M)
                qtrace[k], ctrace[k], centroid_move[k] = q/M, s/M, centroid
            self.sum, self.sumsq, self.sumdiff = float(s), float(q), float(diff)
        potential = .5*qtrace
        primary = cv if self.estimator == "centroid_virial" else thermo
        return primary, {"centroid_virial": cv, "thermodynamic": thermo,
                         "potential": potential, "kinetic": cv-potential,
                         "thermodynamic_kinetic": thermo-potential,
                         "x2": qtrace, "centroid": ctrace, "centroid2": ctrace**2,
                         "accepted": accepted.astype(bool), "centroid_move": centroid_move.astype(bool)}

    def _parameters(self):
        return {"beta": self.beta, "beads": self.beads,
                "centroid_probability": self.centroid_probability,
                "bead_proposal": self.bead_proposal, "centroid_proposal": self.centroid_proposal,
                "estimator": self.estimator, "sampler": self.sampler, "mode_proposal": self.mode_proposal,
                "mode_selection": self.mode_selection}

    def _state(self):
        return {**self._parameters(), "path": self.path.tolist(), "sum": self.sum,
                "sumsq": self.sumsq, "sumdiff": self.sumdiff,
                "mode_coefficients": self.mode_coefficients.tolist()}

    def _restore(self, state):
        if any(state[k] != v for k, v in self._parameters().items()):
            raise ValueError("PIMC checkpoint parameters differ")
        p = np.asarray(state["path"], dtype=np.float64)
        if p.shape != (self.beads,) or not np.isfinite(p).all():
            raise ValueError("invalid PIMC path checkpoint")
        self.path = p.copy()
        measured = np.array(self._path_sums())
        saved = np.array([state["sum"], state["sumsq"], state["sumdiff"]])
        if not np.allclose(measured, saved, rtol=1e-9, atol=1e-8):
            raise ValueError("PIMC checkpoint bookkeeping differs from path")
        self.sum, self.sumsq, self.sumdiff = (float(v) for v in saved)
        self.mode_coefficients = np.ascontiguousarray(state["mode_coefficients"], dtype=np.float64)
        if self.sampler == "normal_mode" and not np.allclose(self.mode_coefficients @ self._basis, p,
                                                            rtol=1e-9, atol=1e-8):
            raise ValueError("PIMC checkpoint Fourier modes differ from path")


def make_model(spec):
    """Construct a model from a compact JSON-compatible workload specification."""
    config = dict(spec)
    kind = config.pop("kind")
    if kind in {"integral", "integration"}:
        return Independent(**config)
    if kind in {"x2", "indicator", "cos"}:
        return Independent(integrand=kind, **config)
    constructors = {"vmc": VMC, "ising": Ising, "pimc": PIMC}
    if kind not in constructors:
        raise ValueError(f"unknown V0.3 workload {kind}")
    return constructors[kind](**config)
