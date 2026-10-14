"""Independent general-log-density Metropolis/copula operator for V0.4.

No production sampler is imported. Symmetric proposal bins and exact grid-state
action ratios retain the intended *discrete* target at rho0. Copula acceptance
is integrated across bins; grid/domain/quadrature refinements remain necessary.
Hamiltonian/observable parameters are separate from the target/operator.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import math
import time
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import csr_matrix
from scipy.special import ndtr, ndtri

VERSION = "v04-general-copula-operator-1"


@dataclass
class Transition:
    grid: np.ndarray
    matrix: csr_matrix
    healthy_weights: np.ndarray
    log_target: np.ndarray
    proposal: float
    rho: float
    bound: float
    dx: float
    quadrature: int
    operator_sha256: str


def _evaluate(function, x):
    values = np.broadcast_to(np.asarray(function(x), float), x.shape)
    if not np.isfinite(values).all():
        raise ValueError("Target/observable callable must return finite grid values")
    return values


def build_transition(log_density, *, proposal=1.5, rho=0., bound=6., dx=.02, quadrature=16):
    if not np.isfinite([proposal, rho, bound, dx]).all():
        raise ValueError("finite transition parameters required")
    if proposal <= 0 or bound <= 0 or dx <= 0 or abs(rho) >= 1:
        raise ValueError("positive proposal/bound/dx and |rho|<1 required")
    intervals = round(2*bound/dx)
    if intervals < 2 or abs(intervals*dx-2*bound) > 1e-10:
        raise ValueError("domain must contain an integer number of dx intervals")
    if not isinstance(quadrature, int) or quadrature < 2:
        raise ValueError("integer quadrature>=2 required")
    grid = np.linspace(-bound, bound, intervals+1)
    target = _evaluate(log_density, grid)
    weights = np.exp(target-target.max()); weights /= weights.sum()
    offsets = np.arange(-math.ceil(proposal/dx+.5), math.ceil(proposal/dx+.5)+1)
    left = np.maximum(0., ((offsets-.5)*dx/proposal+1)/2)
    right = np.minimum(1., ((offsets+.5)*dx/proposal+1)/2)
    widths = np.maximum(0., right-left)
    keep = widths > 0
    offsets, left, right, widths = (a[keep] for a in (offsets, left, right, widths))
    n = len(grid)
    rows = np.broadcast_to(np.arange(n)[:, None], (n, len(offsets)))
    cols = rows+offsets[None, :]
    valid = (cols >= 0) & (cols < n) & (offsets[None, :] != 0)
    log_ratio = target[np.clip(cols, 0, n-1)]-target[:, None]
    accept = np.exp(np.minimum(log_ratio, 0.))
    if rho == 0:
        probability = accept*widths[None, :]
    else:
        z_accept = ndtri(np.minimum(accept, np.nextafter(1., 0.)))
        nodes, gauss_weights = leggauss(quadrature)
        mid, half = (left+right)/2, widths/2
        integrated = np.zeros_like(accept)
        for node, weight in zip(nodes, gauss_weights):
            z_proposal = ndtri(mid+node*half)
            conditional = ndtr((z_accept-rho*z_proposal[None, :])/math.sqrt(1-rho*rho))
            conditional[log_ratio >= 0] = 1.
            integrated += weight*conditional/2
        probability = integrated*widths[None, :]
    matrix = csr_matrix((probability[valid], (rows[valid], cols[valid])), shape=(n, n))
    off_mass = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(off_mass > 1+1e-12):
        raise ArithmeticError("Transition row exceeds unit probability")
    matrix += csr_matrix((np.maximum(0., 1-off_mass), (np.arange(n), np.arange(n))), shape=(n, n))
    matrix.sort_indices()
    digest = hashlib.sha256()
    for a in (grid.astype("<f8"), matrix.indptr.astype("<i8"), matrix.indices.astype("<i8"), matrix.data.astype("<f8")):
        digest.update(a.tobytes())
    return Transition(grid, matrix, weights, target, proposal, rho, bound, dx, quadrature, digest.hexdigest())


def stationary_distribution(operator, tolerance=2e-14, maximum_iterations=10000):
    probability = operator.healthy_weights.copy()
    transpose = operator.matrix.T.tocsr()
    for iteration in range(1, maximum_iterations+1):
        following = transpose @ probability
        following /= following.sum()
        change = np.sum(abs(following-probability)); probability = following
        if change < tolerance:
            residual = float(np.sum(abs(transpose @ probability-probability)))
            return probability, iteration, residual
    raise RuntimeError("Stationary sparse iteration did not converge")


def balance_diagnostics(operator, probability):
    flow = operator.matrix.multiply(probability[:, None]).tocsr()
    imbalance = flow-flow.T
    return {"maximum_flow_imbalance": float(np.max(abs(imbalance.data))) if imbalance.nnz else 0.,
        "total_absolute_probability_current": float(np.sum(abs(imbalance.data))/2),
        "stationarity_l1": float(np.sum(abs(operator.matrix.T @ probability-probability)))}


def observable_arrays(operator, observables=None):
    x = operator.grid
    arrays = {"x": x, "x2": x*x, "x3": x**3, "x4": x**4}
    for name, function in (observables or {}).items():
        arrays[name] = _evaluate(function, x) if callable(function) else np.asarray(function, float)
        if arrays[name].shape != x.shape or not np.isfinite(arrays[name]).all():
            raise ValueError("Observable arrays must match finite grid")
    return arrays


def stationary_moments(operator, observables=None):
    start = time.perf_counter()
    probability, iterations, residual = stationary_distribution(operator)
    arrays = observable_arrays(operator, observables)
    moments = {name: float(probability @ values) for name, values in arrays.items()}
    healthy = {name: float(operator.healthy_weights @ values) for name, values in arrays.items()}
    return {"version": VERSION, "operator_sha256": operator.operator_sha256,
        "rho": operator.rho, "proposal": operator.proposal, "bound": operator.bound,
        "dx": operator.dx, "quadrature": operator.quadrature, "grid_points": len(operator.grid),
        "moments": moments, "healthy_discrete_moments": healthy,
        "biases_against_healthy_grid": {name: moments[name]-healthy[name] for name in moments},
        "position_variance": moments["x2"]-moments["x"]**2,
        "position_third_central_moment": moments["x3"]-3*moments["x"]*moments["x2"]+2*moments["x"]**3,
        "row_probability_error": float(np.max(abs(np.asarray(operator.matrix.sum(axis=1)).ravel()-1))),
        "normalization_error": float(abs(probability.sum()-1)),
        "stationarity_l1_residual": residual, "power_iterations": iterations,
        "healthy_balance": balance_diagnostics(operator, operator.healthy_weights),
        "stationary_balance": balance_diagnostics(operator, probability),
        "edge_probability": float(probability[abs(operator.grid) >= operator.bound-.5].sum()),
        "seconds_for_stationary_solution": time.perf_counter()-start}


def phase_prediction(operator, *, healthy_rows=16384, faulty_rows=49152, recent_window=8192,
                     observables=None, tolerance=2e-14):
    if min(healthy_rows, faulty_rows, recent_window) < 1 or recent_window > faulty_rows:
        raise ValueError("Positive phase budgets and recent_window<=faulty_rows required")
    probability, _, _ = stationary_distribution(operator, tolerance)
    transpose = operator.matrix.T.tocsr()
    arrays = observable_arrays(operator, observables)
    difference = operator.healthy_weights-probability
    sums = {name: 0. for name in arrays}; recent_sums = sums.copy()
    for step in range(1, faulty_rows+1):
        difference = transpose @ difference
        for name, values in arrays.items():
            contribution = float(difference @ values)
            sums[name] += contribution
            if step > faulty_rows-recent_window:
                recent_sums[name] += contribution
        if np.sum(abs(difference)) < tolerance:
            break
    remaining = faulty_rows-step
    norm = float(np.sum(abs(difference)))
    results = {}
    for name, values in arrays.items():
        clean, corrupt = float(operator.healthy_weights @ values), float(probability @ values)
        mixture = (healthy_rows*clean+faulty_rows*corrupt)/(healthy_rows+faulty_rows)
        results[name] = {"healthy_mean": clean, "stationary_fault_mean": corrupt,
            "cumulative_prediction": mixture+sums[name]/(healthy_rows+faulty_rows),
            "recent_prediction": corrupt+recent_sums[name]/recent_window,
            "stationary_mixture": mixture, "transient_correction": sums[name]/(healthy_rows+faulty_rows),
            "truncation_absolute_bound": remaining*np.max(abs(values))*norm/(healthy_rows+faulty_rows)}
    return {"operator_sha256": operator.operator_sha256, "healthy_rows": healthy_rows,
        "faulty_rows": faulty_rows, "recent_window": recent_window,
        "transient_rows_propagated": step, "remaining_density_l1": norm, "observables": results}


def save_transition_predictions(output):
    """Small serial pre-confirmation operator/convergence campaign, never fitted."""
    import csv
    import json
    from pathlib import Path
    import resource
    from rnglab.v04.references import (optimized_gaussian_alpha, trial_local_energy,
        trial_log_density, gaussian_local_energy, quadrature_targets)
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    cases = [
        {"name": "gaussian_a07", "alpha": .7, "trial_quartic": 0., "coupling": .1,
         "proposal": 1.5, "rhos": [-.05, -.035, 0., .035, .05]},
        {"name": "gaussian_optimized_l01", "alpha": optimized_gaussian_alpha(.1),
         "trial_quartic": 0., "coupling": .1, "rhos": [-.05, 0., .05]},
        {"name": "gaussian_optimized_l05", "alpha": optimized_gaussian_alpha(.5),
         "trial_quartic": 0., "coupling": .5, "rhos": [-.05, 0., .05]},
        {"name": "nongaus_a1_b02_l02", "alpha": 1., "trial_quartic": .2, "coupling": .2,
         "proposal": 1.5, "rhos": [-.1, -.05, 0., .05, .1]},
    ]
    for case in cases:
        case.setdefault("proposal", 1.5*math.sqrt(.7/case["alpha"]))
    grids = [(6., .04, 16), (6., .02, 16), (6., .01, 16), (8., .01, 16), (8., .01, 32)]
    config = {"version": VERSION, "cases": cases, "grids": grids,
        "phase": {"burn": 8192, "retained_healthy": 16384, "retained_faulty": 49152, "recent_window": 8192},
        "status": "Independent operator predictions before V0.4 confirmation; no simulation outcomes used",
        "reference_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "physics_reference_sha256": hashlib.sha256((Path(__file__).parent / "references.py").read_bytes()).hexdigest(),
        "kernel_generalization": "NonGaussiantrialgenuinelychangeslogpi;lambdaatfixedGaussianalphachangesonlyenergyobservable",
        "classical_match": "Foralpha1,b.2,classicalV=x²+.1x⁴atT1hasexactlysametargetandoperator"}
    path = output / "transition_config.json"
    # JSON normalizes tuple grid specifications into lists on disk.
    if path.exists() and json.loads(path.read_text()) != json.loads(json.dumps(config)):
        raise ValueError("Pre-confirmation reference identity changed; use another output")
    path.write_text(json.dumps(config, indent=2))
    start = time.perf_counter(); results = []; fresh = reused = 0
    for case in cases:
        alpha, b, coupling = case["alpha"], case["trial_quartic"], case["coupling"]
        target = lambda x, a=alpha, b=b: trial_log_density(x, a, b)
        energy = lambda x, a=alpha, b=b, lam=coupling: trial_local_energy(x, a, lam, b)
        observables = {"energy": energy, "classical_potential": lambda x, a=alpha, b=b: a*x*x+b*x**4/2}
        if case["name"] == "gaussian_a07":
            observables.update({f"energy_lambda_{lam:g}": lambda x, lam=lam: gaussian_local_energy(x, .7, lam)
                                for lam in [0., .1, .5]})
        for rho in case["rhos"]:
            for bound, dx, quadrature in grids:
                identity = f"{case['name']}_rho{rho:g}_b{bound:g}_dx{dx:g}_q{quadrature}"
                checkpoint = output / "operators" / (identity+".json")
                if checkpoint.exists():
                    result = json.loads(checkpoint.read_text())
                    reused += 1
                else:
                    operator = build_transition(target, proposal=case["proposal"], rho=rho,
                                                bound=bound, dx=dx, quadrature=quadrature)
                    result = {"identity": identity, "case": case, **stationary_moments(operator, observables),
                              "healthy_continuum_quadrature": quadrature_targets(alpha, coupling, b)}
                    if (bound, dx, quadrature) == (8., .01, 32):
                        result["finite_phase_prediction"] = phase_prediction(operator, observables=observables)
                    checkpoint.parent.mkdir(parents=True, exist_ok=True)
                    temporary = checkpoint.with_suffix(".tmp"); temporary.write_text(json.dumps(result))
                    temporary.replace(checkpoint)
                    fresh += 1
                results.append(result)
        print("Operator references saved:", case["name"], "total", len(results), flush=True)
    finest = [r for r in results if r["bound"] == 8. and r["dx"] == .01 and r["quadrature"] == 32]
    convergence, symmetry = [], []
    for fine in finest:
        for result in results:
            if result["case"]["name"] != fine["case"]["name"] or result["rho"] != fine["rho"]:
                continue
            convergence.append({"identity": result["identity"], "finest_identity": fine["identity"],
                "moment_differences": {name: result["moments"][name]-fine["moments"][name]
                                       for name in fine["moments"]}})
        if fine["rho"] > 0:
            minus = next(r for r in finest if r["case"]["name"] == fine["case"]["name"] and r["rho"] == -fine["rho"])
            symmetry.append({"case": fine["case"]["name"], "rho": fine["rho"],
                "odd_x_error": fine["moments"]["x"]+minus["moments"]["x"],
                "odd_x3_error": fine["moments"]["x3"]+minus["moments"]["x3"],
                "even_energy_error": fine["moments"]["energy"]-minus["moments"]["energy"]})
    summary = {"deterministic_operator_calculations": len(results), "independent_random_trials": 0,
        "fresh_calculations_this_invocation": fresh, "reused_operator_checkpoints": reused,
        "wall_seconds": time.perf_counter()-start, "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "maximum_row_probability_error": max(r["row_probability_error"] for r in results),
        "maximum_stationarity_residual": max(r["stationarity_l1_residual"] for r in results),
        "finest_predictions": finest, "convergence": convergence, "symmetry": symmetry,
        "generality_limit": "SameGaussianparameters/samecorruptionshareoperatoracrossHamiltonianlambda;onlynonGaussianlogtargetchangesgenuineshape. ClassicalmatchedtargethasidenticalMCMCmechanism."}
    (output / "transition_summary.json").write_text(json.dumps(summary, indent=2))
    fields = ["case", "rho", "proposal", "bound", "dx", "quadrature", "operator_sha256", "x", "x2", "x3", "x4", "energy", "classical_potential", "stationarity_l1_residual"]
    with (output / "stationary_predictions.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for result in results:
            row = {k: result[k] for k in fields if k in result}; row["case"] = result["case"]["name"]
            row.update({k: result["moments"][k] for k in fields if k in result["moments"]}); writer.writerow(row)
    return summary


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[3] / "results/v0.4/references")
    args = parser.parse_args()
    result = save_transition_predictions(args.output)
    print("Independent operator calculations:", result["deterministic_operator_calculations"],
          "seconds", result["wall_seconds"], "maximum residual", result["maximum_stationarity_residual"])
