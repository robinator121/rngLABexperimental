#!/usr/bin/env python3
"""Independent stdlib analytic checks for primitive harmonic PIMC, no sampling.

Outputs deterministic finite-M reference moments and checks closed-form energy
against eigenmode sums and a partition derivative. Does not import RNGLab.
"""
from __future__ import annotations

import json
import math
from pathlib import Path


def energy_closed(m: int, beta: float) -> float:
    z = beta / (2 * m)
    return 0.5 / (math.tanh(m * math.asinh(z)) * math.sqrt(1 + z*z))


def log_z(m: int, beta: float) -> float:
    return -math.log(2 * math.sinh(m * math.asinh(beta / (2*m))))


def reference(m: int, beta: float) -> dict:
    laplace = [4 * math.sin(math.pi*k/m)**2 for k in range(m)]
    q = [beta/m + m/beta * ll for ll in laplace]
    a_td = [1/(2*m) - m/(2*beta**2)*ll for ll in laplace]
    a_cv = [1/(2*m)] + [1/m]*(m-1)
    td_mean = m/(2*beta) + sum(a/qq for a, qq in zip(a_td, q))
    cv_mean = 1/(2*beta) + sum(a/qq for a, qq in zip(a_cv, q))
    td_var = 2 * sum((a/qq)**2 for a, qq in zip(a_td, q))
    cv_var = 2 * sum((a/qq)**2 for a, qq in zip(a_cv, q))
    cv_skew = 8*sum((a/qq)**3 for a, qq in zip(a_cv, q))/cv_var**1.5
    h = 1.e-5 * beta
    derivative = -(log_z(m, beta+h)-log_z(m, beta-h))/(2*h)
    value = energy_closed(m, beta)
    assert abs(td_mean-value) < 2.e-12
    assert abs(cv_mean-value) < 2.e-12
    assert abs(derivative-value) < 1.e-8
    # Spectral inverse and direct action form agree for each test model.
    covariance = [[sum(math.cos(2*math.pi*k*(i-j)/m)/q[k]
                       for k in range(m))/m for j in range(m)] for i in range(m)]
    precision = [[(beta/m if i==j else 0.0) for j in range(m)] for i in range(m)]
    for i in range(m):
        j = (i+1) % m
        precision[i][i] += m/beta
        precision[j][j] += m/beta
        precision[i][j] -= m/beta
        precision[j][i] -= m/beta
    residual = max(abs(sum(precision[i][k]*covariance[k][j] for k in range(m))
                       - (1.0 if i==j else 0.0)) for i in range(m) for j in range(m))
    assert residual < 2.e-11
    centroid_variance = 0.5/beta**2
    return {
        "beads":m, "beta":beta, "energy_finite_M":value,
        "energy_continuum":0.5/math.tanh(beta/2),
        "TD_mean_eigenmode_sum":td_mean, "CV_mean_eigenmode_sum":cv_mean,
        "energy_partition_derivative":derivative,
        "TD_marginal_variance":td_var, "CV_marginal_variance":cv_var,
        "CV_centroid_marginal_variance":centroid_variance,
        "CV_internal_marginal_variance":cv_var-centroid_variance,
        "CV_centroid_marginal_variance_fraction":centroid_variance/cv_var,
        "TD_CV_covariance":2*sum(a*b/qq**2 for a,b,qq in zip(a_td,a_cv,q)),
        "CV_marginal_skewness":cv_skew,
        "CV_iid_8192_mean_skewness":cv_skew/math.sqrt(8192),
        "centroid_mean_square":1/beta,
        "internal_mean_square":sum(1/q[k] for k in range(1,m))/m,
        "precision_inverse_max_residual":residual,
        "precision_eigenvalues":q,
    }


def main() -> None:
    models = [reference(m,beta) for m in (1,2,3,8,16,32) for beta in (0.5,2.,8.)]
    x = [math.sin(j+0.25) for j in range(16)]
    beta = 2.0
    c = sum(x)/len(x)
    r2 = sum((v-c)**2 for v in x)/len(x)
    cv_direct = 1/(2*beta) + sum(v*v for v in x)/len(x)-0.5*c*c
    cv_decomposed = 1/(2*beta) + 0.5*c*c+r2
    assert abs(cv_direct-cv_decomposed) < 1.e-14
    def action(path):
        return (len(path)/(2*beta)*sum((path[(j+1)%len(path)]-path[j])**2
                for j in range(len(path))) + beta/(2*len(path))*sum(v*v for v in path))
    c_new = 0.87
    x_shifted = [v+c_new-c for v in x]
    shift_error = abs(action(x_shifted)-action(x)-beta/2*(c_new*c_new-c*c))
    assert shift_error < 1.e-13
    result = {"protocol":"RNGLab_PhaseII_independent_analytic_reference_v1",
              "sampling":"none; deterministic analytic checks only",
              "models_checked":len(models), "all_checks_passed":True,
              "centroid_translation_action_identity_error":shift_error,
              "CV_decomposition_error":abs(cv_direct-cv_decomposed),
              "models":models}
    path=Path(__file__).with_name("harmonic_reference_values.json")
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"checks_passed":True,"models":len(models),"output":str(path)},indent=2))


if __name__ == "__main__":
    main()
