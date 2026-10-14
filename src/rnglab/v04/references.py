"""Independent natural-unit oscillator references; no production kernels used.

For psi=Z^-1/2 exp(-alpha*x²/2-b*x⁴/4), log(pi)=-alpha*x²-b*x⁴/2.
Changing Hamiltonian coupling lambda at fixed (alpha,b) changes observables,
not this Metropolis target. The Gaussian family has b=0.
"""
from __future__ import annotations

from functools import lru_cache
import math
import numpy as np
from scipy.integrate import quad
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq

VERSION = "v04-independent-oscillator-reference-1"


def _validate(alpha, coupling=0., trial_quartic=0.):
    if not np.isfinite([alpha, coupling, trial_quartic]).all():
        raise ValueError("finite oscillator parameters required")
    if alpha < 0 or coupling < 0 or trial_quartic < 0 or (alpha == 0 and trial_quartic == 0):
        raise ValueError("normalizable alpha>=0, trial_quartic>=0 and stable coupling>=0 required")


def trial_log_density(x, alpha, trial_quartic=0.):
    _validate(alpha, trial_quartic=trial_quartic)
    x = np.asarray(x, dtype=float)
    return -alpha*x*x-trial_quartic*x**4/2


def trial_local_energy(x, alpha, coupling, trial_quartic=0.):
    _validate(alpha, coupling, trial_quartic)
    x = np.asarray(x, dtype=float)
    return (alpha/2+(1-alpha*alpha+3*trial_quartic)*x*x/2
            +(coupling-alpha*trial_quartic)*x**4-trial_quartic**2*x**6/2)


def gaussian_local_energy(x, alpha, coupling):
    return trial_local_energy(x, alpha, coupling, 0.)


def gaussian_targets(alpha, coupling):
    _validate(alpha, coupling)
    if alpha <= 0:
        raise ValueError("Gaussian alpha must be positive")
    moments = {"x": 0., "x2": 1/(2*alpha), "x3": 0., "x4": 3/(4*alpha**2),
               "x6": 15/(8*alpha**3), "x8": 105/(16*alpha**4)}
    energy = (alpha+1/alpha)/4+3*coupling/(4*alpha**2)
    A, B, C = alpha/2, (1-alpha**2)/2, coupling
    second = (A*A+2*A*B*moments["x2"]+(B*B+2*A*C)*moments["x4"]
              +2*B*C*moments["x6"]+C*C*moments["x8"])
    return {"family": "gaussian_trial", "alpha": alpha, "coupling": coupling,
            "trial_quartic": 0., **moments, "energy": energy,
            "energy_variance": float(max(0., second-energy*energy)),
            "density_normalization": math.sqrt(alpha/math.pi),
            "wavefunction_normalization": (alpha/math.pi)**.25,
            "kernel_caveat": "Hamiltonian coupling does not enter the sampling density at fixed alpha"}


def optimized_gaussian_alpha(coupling):
    if not np.isfinite(coupling) or coupling < 0:
        raise ValueError("stable finite coupling>=0 required")
    if coupling == 0:
        return 1.
    # dE/dalpha=0 iff alpha³-alpha-6lambda=0; unique root above1.
    upper = max(2., (6*coupling+1)**(1/3)+1)
    return float(brentq(lambda a: a*a*a-a-6*coupling, 1., upper, xtol=5e-15))


@lru_cache(maxsize=128)
def quadrature_targets(alpha, coupling, trial_quartic=0.):
    """Density quadrature and an independent |psi'|²/2 kinetic expectation."""
    _validate(alpha, coupling, trial_quartic)
    density = lambda x: math.exp(-alpha*x*x-trial_quartic*x**4/2)
    normalization, normalization_error = quad(density, -np.inf, np.inf, epsabs=2e-12, epsrel=2e-12)
    moments = {"x": 0., "x3": 0.}
    moment_errors = {}
    for power in [2, 4, 6, 8, 10, 12]:
        integral, error = quad(lambda x: x**power*density(x), -np.inf, np.inf, epsabs=2e-12, epsrel=2e-12)
        moments["x"+str(power)] = integral/normalization
        moment_errors["x"+str(power)] = error/normalization
    # g'=alpha*x+b*x³, independent of the expanded local-energy polynomial.
    kinetic, kinetic_error = quad(lambda x: .5*(alpha*x+trial_quartic*x**3)**2*density(x),
                                  -np.inf, np.inf, epsabs=2e-12, epsrel=2e-12)
    potential = .5*moments["x2"]+coupling*moments["x4"]
    energy = kinetic/normalization+potential
    local_mean, local_error = quad(lambda x: float(trial_local_energy(x, alpha, coupling, trial_quartic))*density(x),
                                   -np.inf, np.inf, epsabs=2e-12, epsrel=2e-12)
    local_second, second_error = quad(lambda x: float(trial_local_energy(x, alpha, coupling, trial_quartic))**2*density(x),
                                     -np.inf, np.inf, epsabs=2e-12, epsrel=2e-12)
    return {"alpha": alpha, "coupling": coupling, "trial_quartic": trial_quartic,
        "unnormalized_density_integral": normalization, "wavefunction_normalization": normalization**-.5,
        "density_normalization": 1/normalization, "normalization_quadrature_error": normalization_error,
        **moments, "moment_quadrature_errors": moment_errors,
        "kinetic_energy": kinetic/normalization, "potential_energy": potential,
        "energy": energy, "local_energy_quadrature": local_mean/normalization,
        "energy_variance": local_second/normalization-energy*energy,
        "energy_method_disagreement": abs(energy-local_mean/normalization),
        "kinetic_quadrature_error": kinetic_error/normalization,
        "local_energy_quadrature_error": local_error/normalization,
        "local_second_quadrature_error": second_error/normalization,
        "target_type": "Exact trial-family variational expectation, numerically integrated; not exact ground state"}


def classical_quartic_log_density(x, coupling, temperature=1.):
    if temperature <= 0 or coupling < 0 or not np.isfinite([coupling, temperature]).all():
        raise ValueError("stable coupling>=0 and finite temperature>0 required")
    x = np.asarray(x, float)
    return -(.5*x*x+coupling*x**4)/temperature


def classical_quartic_targets(coupling, temperature=1.):
    """Matched nonquantum canonical target, normalized by independent quadrature."""
    classical_quartic_log_density(0., coupling, temperature)
    reference = quadrature_targets(.5/temperature, 0., 2*coupling/temperature)
    return {"coupling": coupling, "temperature": temperature,
        **{name: reference[name] for name in ["x", "x2", "x3", "x4"]},
        "potential_energy": .5*reference["x2"]+coupling*reference["x4"],
        "density_normalization": reference["density_normalization"],
        "mapping_to_trial_density": {"alpha": .5/temperature, "trial_quartic": 2*coupling/temperature},
        "target_type": "Classical canonical density exp(-V/T), not a quantum energy reference"}


def schrodinger_reference(coupling, bound=8., dx=.02, states=3):
    """Independent Dirichlet finite-difference eigenproblem; dx/bound checks needed."""
    if coupling < 0 or bound <= 0 or dx <= 0 or not np.isfinite([coupling, bound, dx]).all():
        raise ValueError("stable finite coupling and positive bound/dx required")
    intervals = round(2*bound/dx)
    if abs(intervals*dx-2*bound) > 1e-10 or states < 1 or states >= intervals:
        raise ValueError("invalid Schrödinger grid/states")
    x = np.linspace(-bound, bound, intervals+1)[1:-1]
    diagonal = 1/dx**2+.5*x*x+coupling*x**4
    off = np.full(len(x)-1, -.5/dx**2)
    energies, eigenvectors = eigh_tridiagonal(diagonal, off, select="i", select_range=(0, states-1))
    ground_weights = eigenvectors[:, 0]**2
    return {"coupling": coupling, "bound": bound, "dx": dx,
        "energies": energies.tolist(), "ground_energy": float(energies[0]),
        "ground_x2": float(ground_weights @ (x*x)), "ground_x4": float(ground_weights @ x**4),
        "ground_probability_edge": float(ground_weights[abs(x) >= bound-.5].sum()),
        "reference_type": "Approximate physical eigenvalues; distinct from trial variational energy"}


def save_physics_validation(output):
    """Bounded deterministic evidence; no Monte Carlo trials or fitted outcomes."""
    import hashlib
    import json
    from pathlib import Path
    import time
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    config = {"version": VERSION, "fixed_gaussian_alpha": .7, "couplings": [0., .1, .2, .5],
        "nongaus_parameters": {"alpha": 1., "trial_quartic": .2, "coupling": .2},
        "schrodinger_grids": [[6., .04], [6., .02], [6., .01], [8., .01]],
        "reference_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "independence": "Deterministic formulas/quadrature/eigenproblems; zero random-trial identities"}
    path = output / "physics_config.json"
    if path.exists() and json.loads(path.read_text()) != config:
        raise ValueError("Reference configuration changed; use a separate output directory")
    path.write_text(json.dumps(config, indent=2))
    start = time.perf_counter(); trial_states = []
    for coupling in config["couplings"]:
        for alpha, name in [(.7, "fixed_gaussian"), (optimized_gaussian_alpha(coupling), "optimized_gaussian")]:
            analytic = gaussian_targets(alpha, coupling)
            numerical = quadrature_targets(alpha, coupling)
            trial_states.append({"name": name, "analytic": analytic, "quadrature": numerical,
                "energy_difference": numerical["energy"]-analytic["energy"]})
    nongaus = quadrature_targets(1., .2, .2)
    eigenproblems = [schrodinger_reference(coupling, bound=bound, dx=dx)
                    for coupling in config["couplings"] for bound, dx in config["schrodinger_grids"]]
    result = {"version": VERSION, "wall_seconds": time.perf_counter()-start,
        "gaussian_trial_states": trial_states, "nongaus_trial": nongaus,
        "matched_classical_target": {"log_density": "-x²-.1x⁴", "temperature": 1.,
            "x2": nongaus["x2"], "x4": nongaus["x4"],
            "classical_potential_energy": nongaus["x2"]+.1*nongaus["x4"]},
        "schrodinger_references": eigenproblems,
        "kernel_caveat": "Samealphaandtrial_quarticgiveidenticaltarget/kernelregardlessHamiltoniancoupling;optimizedGaussianswithscaledproposalsarealsoequivalentinstandardizedcoordinates",
        "independent_random_trials": 0}
    (output / "physics_validation.json").write_text(json.dumps(result, indent=2))
    return result


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[3] / "results/v0.4/references")
    args = parser.parse_args()
    results = save_physics_validation(args.output)
    print("Independent physics references saved:", args.output, "seconds", results["wall_seconds"])
