"""Independent sparse-grid reference for a copula-coupled VMC transition.

This module does not import or call the simulation kernels. A box displacement
d*(2u-1) is rounded into a symmetric position-grid bin; the accepted grid state
uses its exact Gaussian action difference. Conditional copula acceptance is
integrated over the whole proposal bin. Outside-grid proposals are rejected.
Thus rho=0 has exact discrete Gaussian detailed balance; refining dx approaches
the continuous-position process, and increasing the domain checks truncation.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import time

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import csr_matrix
from scipy.special import ndtr, ndtri

VERSION = "v03-copula-grid-reference-1"


@dataclass
class GridTransition:
    grid: np.ndarray
    matrix: csr_matrix
    alpha: float
    proposal: float
    rho: float
    dx: float
    bound: float
    quadrature: int
    gaussian: np.ndarray


def conditional_acceptance(proposal_uniform, threshold, rho):
    """P(U_accept <= threshold | U_proposal=u) for a Gaussian copula.

    The convention of strict vs nonstrict acceptance differs only at ideal
    continuous measure-zero endpoints. Full downhill acceptance is explicit.
    """
    if not np.isfinite(rho) or abs(rho) >= 1:
        raise ValueError("the reference requires |rho|<1")
    u, p = np.broadcast_arrays(np.asarray(proposal_uniform, float), np.asarray(threshold, float))
    if (not np.isfinite(u).all() or not np.isfinite(p).all()
            or np.any(u <= 0) or np.any(u >= 1) or np.any(p < 0) or np.any(p > 1)):
        raise ValueError("u in (0,1), threshold in [0,1] required")
    if rho == 0:
        return p.copy()
    z = ndtri(p)
    return ndtr((z-rho*ndtri(u))/math.sqrt(1-rho*rho))


def build_transition(alpha=.7, proposal=1.5, rho=0., bound=6., dx=.02, quadrature=8):
    if not np.isfinite([alpha, proposal, rho, bound, dx]).all():
        raise ValueError("finite reference parameters required")
    if alpha <= 0 or proposal <= 0 or bound <= 0 or dx <= 0 or abs(rho) >= 1:
        raise ValueError("positive alpha/proposal/bound/dx and |rho|<1 required")
    intervals = round(2*bound/dx)
    if abs(intervals*dx-2*bound) > 1e-10 or intervals < 2:
        raise ValueError("2*bound must be an integer number of grid intervals")
    if not isinstance(quadrature, int) or quadrature < 2:
        raise ValueError("integer quadrature>=2 required")
    grid = np.linspace(-bound, bound, intervals+1)
    n = len(grid)
    maximum_offset = math.ceil(proposal/dx+.5)
    offsets = np.arange(-maximum_offset, maximum_offset+1)
    bin_left = np.maximum(0., ((offsets-.5)*dx/proposal+1)/2)
    bin_right = np.minimum(1., ((offsets+.5)*dx/proposal+1)/2)
    widths = np.maximum(0., bin_right-bin_left)
    nonempty = widths > 0
    offsets, bin_left, bin_right, widths = (a[nonempty] for a in (offsets, bin_left, bin_right, widths))
    rows = np.broadcast_to(np.arange(n)[:, None], (n, len(offsets)))
    cols = rows+offsets[None, :]
    valid = (cols >= 0) & (cols < n) & (offsets[None, :] != 0)
    proposed = grid[np.clip(cols, 0, n-1)]
    uphill_action = alpha*(proposed**2-grid[:, None]**2)
    threshold = np.exp(-np.maximum(uphill_action, 0.))
    if rho == 0:
        # This special branch deliberately avoids a CDF/inverse roundtrip and
        # establishes exact discrete target balance up to floating arithmetic.
        probabilities = threshold*widths[None, :]
    else:
        nodes, weights = leggauss(quadrature)
        midpoints, halves = (bin_left+bin_right)/2, widths/2
        zthreshold = ndtri(np.minimum(threshold, np.nextafter(1., 0.)))
        integrated = np.zeros_like(threshold)
        for node, weight in zip(nodes, weights):
            zproposal = ndtri(midpoints+halves*node)
            value = ndtr((zthreshold-rho*zproposal[None, :])/math.sqrt(1-rho*rho))
            value[uphill_action <= 0] = 1.
            integrated += weight*value/2
        probabilities = integrated*widths[None, :]
    matrix = csr_matrix((probabilities[valid], (rows[valid], cols[valid])), shape=(n, n))
    row_offdiagonal = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(row_offdiagonal > 1+1e-12):
        raise ArithmeticError("transition row has probability above one")
    matrix += csr_matrix((np.maximum(0., 1-row_offdiagonal), (np.arange(n), np.arange(n))), shape=(n, n))
    gaussian = np.exp(-alpha*grid**2)
    gaussian /= gaussian.sum()
    return GridTransition(grid, matrix, alpha, proposal, rho, dx, bound, quadrature, gaussian)


def stationary_distribution(transition, tolerance=2e-14, maximum_iterations=10000):
    distribution = transition.gaussian.copy()
    transpose = transition.matrix.T.tocsr()
    for iteration in range(1, maximum_iterations+1):
        following = transpose @ distribution
        following /= following.sum()
        change = np.sum(abs(following-distribution))
        distribution = following
        if change < tolerance:
            residual = float(np.sum(abs(transpose @ distribution-distribution)))
            return distribution, iteration, residual
    raise RuntimeError("sparse power iteration did not converge")


def detailed_balance_diagnostics(matrix, distribution):
    flow = matrix.multiply(distribution[:, None]).tocsr()
    current = flow-flow.T
    return {"maximum_absolute_flow_difference": float(np.max(abs(current.data))) if current.nnz else 0.,
            "total_absolute_current": float(np.sum(abs(current.data))/2),
            "stationarity_l1": float(np.sum(abs(matrix.T @ distribution-distribution)))}


def solve_reference(**spec):
    start = time.perf_counter()
    transition = build_transition(**spec)
    distribution, iterations, residual = stationary_distribution(transition)
    x = transition.grid
    mean = float(distribution @ x)
    x2 = float(distribution @ (x*x))
    energy = transition.alpha/2+(1-transition.alpha**2)*x2/2
    gaussian_x2 = float(transition.gaussian @ (x*x))
    return {"version": VERSION, "spec": spec, "grid_points": len(x), "nonzero_transitions": transition.matrix.nnz,
        "mean_x": mean, "mean_x2": x2, "mean_local_energy": energy,
        "continuum_gaussian_x2": 1/(2*transition.alpha), "discrete_gaussian_x2": gaussian_x2,
        "continuum_variational_energy": (transition.alpha+1/transition.alpha)/4,
        "power_iterations": iterations, "stationarity_l1_residual": residual,
        "row_probability_error": float(np.max(abs(np.asarray(transition.matrix.sum(axis=1)).ravel()-1))),
        "gaussian_balance": detailed_balance_diagnostics(transition.matrix, transition.gaussian),
        "stationary_balance": detailed_balance_diagnostics(transition.matrix, distribution),
        "boundary_probability": float(distribution[abs(x) >= transition.bound-.5].sum()),
        "gaussian_boundary_probability": float(transition.gaussian[abs(x) >= transition.bound-.5].sum()),
        "seconds": time.perf_counter()-start}


def finite_horizon_prediction(transition, healthy_rows=16384, faulty_rows=49152,
                              recent_window=4096, tolerance=2e-14):
    """Expected finite averages, starting in the exact discrete healthy target.

    After onset, propagate the density difference from stationary and sum its
    observable contribution. Once its L1 norm is negligible, an explicit worst
    remaining-horizon bound accompanies truncation. No random draws are used.
    """
    if min(healthy_rows, faulty_rows, recent_window) < 1 or recent_window > faulty_rows:
        raise ValueError("positive horizons and recent_window<=faulty_rows required")
    stationary, _, _ = stationary_distribution(transition, tolerance=tolerance)
    transpose = transition.matrix.T.tocsr()
    difference = transition.gaussian-stationary
    x = transition.grid
    observables = {"x": x, "x2": x*x,
                   "energy": transition.alpha/2+(1-transition.alpha**2)*x*x/2}
    transient_sums = {name: 0. for name in observables}
    recent_transient = {name: 0. for name in observables}
    for attempted in range(1, faulty_rows+1):
        difference = transpose @ difference
        for name, observable in observables.items():
            contribution = float(difference @ observable)
            transient_sums[name] += contribution
            if attempted > faulty_rows-recent_window:
                recent_transient[name] += contribution
        if np.sum(abs(difference)) < tolerance:
            break
    remaining = faulty_rows-attempted
    output = {}
    for name, observable in observables.items():
        clean = float(transition.gaussian @ observable)
        corrupt = float(stationary @ observable)
        complete = ((healthy_rows*clean+faulty_rows*corrupt+transient_sums[name])
                    /(healthy_rows+faulty_rows))
        bound = float(remaining*np.max(abs(observable))*np.sum(abs(difference))/(healthy_rows+faulty_rows))
        output[name] = {"healthy_mean": clean, "faulty_stationary_mean": corrupt,
                       "predicted_cumulative_mean": complete,
                       "simple_stationary_mixture": (healthy_rows*clean+faulty_rows*corrupt)/(healthy_rows+faulty_rows),
                       "predicted_recent_mean": corrupt+recent_transient[name]/recent_window,
                       "stationary_mixture_transient_correction": transient_sums[name]/(healthy_rows+faulty_rows),
                       "remaining_horizon_absolute_bound": bound}
    return {"healthy_rows": healthy_rows, "faulty_rows": faulty_rows,
            "recent_window": recent_window, "density_difference_steps_computed": attempted,
            "remaining_density_l1": float(np.sum(abs(difference))), "observables": output}
