"""Unfitted copula response from differentiated stationarity at rho=0.

Derivatives integrate the proposal-bin inverse-normal moments analytically:
  dH/drho|0 = -phi(z_p)*z_u
  d²H/drho²|0 = z_p*phi(z_p)*(1-z_u²).
No simulation or nonzero-rho stationary outcomes enter these derivatives.
"""
from __future__ import annotations
import math
import numpy as np
from scipy.sparse import csr_matrix
from scipy.special import ndtri
from rnglab.v03.cross_reference import build_transition

VERSION = "v03-copula-response-1"


def transition_derivatives(alpha=.7, proposal=1.5, bound=8., dx=.01):
    reference = build_transition(alpha=alpha, proposal=proposal, rho=0., bound=bound, dx=dx)
    grid = reference.grid
    n = len(grid)
    maximum_offset = math.ceil(proposal/dx+.5)
    offsets = np.arange(-maximum_offset, maximum_offset+1)
    left = np.maximum(0., ((offsets-.5)*dx/proposal+1)/2)
    right = np.minimum(1., ((offsets+.5)*dx/proposal+1)/2)
    nonempty = right > left
    offsets, left, right = (a[nonempty] for a in (offsets, left, right))
    za, zb = ndtri(left), ndtri(right)
    phia, phib = np.exp(-za*za/2)/math.sqrt(2*math.pi), np.exp(-zb*zb/2)/math.sqrt(2*math.pi)
    integral_z = phia-phib
    integral_one_minus_z2 = np.where(np.isfinite(zb), zb, 0.)*phib-np.where(np.isfinite(za), za, 0.)*phia
    rows = np.broadcast_to(np.arange(n)[:, None], (n, len(offsets)))
    cols = rows+offsets[None, :]
    valid = (cols >= 0) & (cols < n) & (offsets[None, :] != 0)
    y = grid[np.clip(cols, 0, n-1)]
    uphill = alpha*(y*y-grid[:, None]**2)
    threshold = np.exp(-np.maximum(uphill, 0.))
    z = ndtri(np.minimum(threshold, np.nextafter(1., 0.)))
    phi = np.exp(-z*z/2)/math.sqrt(2*math.pi)
    first = -phi*integral_z[None, :]
    second = z*phi*integral_one_minus_z2[None, :]
    first[uphill <= 0], second[uphill <= 0] = 0., 0.
    matrices = []
    for values in (first, second):
        matrix = csr_matrix((values[valid], (rows[valid], cols[valid])), shape=(n, n))
        rowsum = np.asarray(matrix.sum(axis=1)).ravel()
        matrix += csr_matrix((-rowsum, (np.arange(n), np.arange(n))), shape=(n, n))
        matrices.append(matrix)
    return reference, matrices[0], matrices[1]


def _poisson_solution(transpose, gaussian, rhs, tolerance=2e-14):
    # Work in the normalized zero-mass subspace, where the transition contracts.
    rhs = rhs-gaussian*rhs.sum()
    term, solution = rhs.copy(), np.zeros_like(rhs)
    for iteration in range(1, 10001):
        solution += term
        term = transpose @ term
        term -= gaussian*term.sum()
        if np.sum(abs(term)) < tolerance:
            solution -= gaussian*solution.sum()
            residual = np.sum(abs(solution-transpose @ solution-rhs))
            return solution, iteration, float(residual)
    raise RuntimeError("differentiated stationary Poisson series did not converge")


def solve_response(alpha=.7, proposal=1.5, bound=8., dx=.01):
    reference, first, second = transition_derivatives(alpha, proposal, bound, dx)
    gaussian = reference.gaussian
    transpose = reference.matrix.T.tocsr()
    rhs1 = first.T @ gaussian
    density1, iterations1, residual1 = _poisson_solution(transpose, gaussian, rhs1)
    rhs2 = second.T @ gaussian+2*(first.T @ density1)
    density2, iterations2, residual2 = _poisson_solution(transpose, gaussian, rhs2)
    x = reference.grid
    energy = alpha/2+(1-alpha*alpha)*x*x/2
    slope = float(density1 @ x)
    qsecond = float(density2 @ (x*x)/2)
    return {"version": VERSION, "alpha": alpha, "proposal": proposal, "bound": bound, "dx": dx,
        "mean_x_linear_coefficient": slope, "mean_x_quadratic_coefficient": float(density2 @ x/2),
        "mean_x2_linear_coefficient": float(density1 @ (x*x)),
        "mean_x2_quadratic_coefficient": qsecond,
        "local_energy_linear_coefficient": float(density1 @ energy),
        "local_energy_quadratic_coefficient": float(density2 @ energy/2),
        "position_variance_quadratic_coefficient": qsecond-slope*slope,
        "first_poisson_residual_l1": residual1, "second_poisson_residual_l1": residual2,
        "first_poisson_terms": iterations1, "second_poisson_terms": iterations2,
        "first_density_mass": float(density1.sum()), "second_density_mass": float(density2.sum()),
        "first_density_odd_error_max": float(np.max(abs(density1+density1[::-1]))),
        "second_density_even_error_max": float(np.max(abs(density2-density2[::-1])))}
