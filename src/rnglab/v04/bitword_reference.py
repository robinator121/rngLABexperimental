"""Independent fair-Markov-bit word CDF and idealized Metropolis reference.

Each stationary bit is Bernoulli(.5), P(same)=(1+rho)/2. Adjacent numeric
uniforms retain only the first53 of64 bits. Dependence *within* retained words
changes their numerical marginal, even though single bits stay fair. Kernel
predictions below explicitly idealize successive words as independent, with
separate word-gap and prefix-quantization bounds; they are not claimed to be
the exact invariant law of the full correlated64-bit source.
"""
from __future__ import annotations
from dataclasses import dataclass
import hashlib
import math
import numpy as np
from scipy.sparse import csr_matrix
from rnglab.v04.transition_reference import Transition, stationary_moments, phase_prediction

VERSION = "v04-markov-bitword-reference-1"


def _validate(rho, bits):
    if not np.isfinite(rho) or abs(rho) > 1:
        raise ValueError("symmetric Markov rho in[-1,1] required")
    if not isinstance(bits, int) or not 1 <= bits <= 53:
        raise ValueError("retained bits must be an integer in[1,53]")


def exact_word_moments(rho, bits=53):
    _validate(rho, bits)
    positions = np.arange(1, bits+1)
    weights = np.exp2(-positions.astype(float))
    covariance = .25*rho**abs(positions[:, None]-positions[None, :])
    mean = .5*weights.sum()
    variance = float(weights @ covariance @ weights)
    return {"rho": rho, "retained_bits": bits, "mean": float(mean), "variance": variance,
        "second_moment": float(variance+mean*mean),
        "infinite_bit_variance": (2+rho)/(12*(2-rho)),
        "zero_word_probability": .5*((1+rho)/2)**(bits-1),
        "meaning": "Numerical-word moments, not bit fairness or entropy estimates"}


def word_cdf_lt(threshold, rho, bits=53):
    """Exact finite-word P(U<threshold), via at mostbits prefix branches.

    Thresholds are the supplied float64 numbers. Scaling by a power of two is
    exact in their representable range; ceil implements the strict discrete
    inequality. No enumeration of2^53 words or Monte Carlo sampling is used.
    """
    _validate(rho, bits)
    t = np.asarray(threshold, dtype=float)
    if not np.isfinite(t).all():
        raise ValueError("CDF thresholds must be finite")
    scaled = np.ceil(np.clip(t, 0., 1.)*2.**bits).astype(np.uint64)
    probability = np.zeros(t.shape, float)
    equal_prefix = np.ones(t.shape, float)
    previous = np.zeros(t.shape, np.uint8)
    same = (1+rho)/2
    for index in range(bits):
        digit = ((scaled >> np.uint64(bits-index-1)) & np.uint64(1)).astype(np.uint8)
        pzero = np.full(t.shape, .5) if index == 0 else np.where(previous == 0, same, 1-same)
        probability += np.where(digit == 1, equal_prefix*pzero, 0.)
        equal_prefix *= np.where(digit == 0, pzero, 1-pzero)
        previous = digit
    return np.where(t <= 0, 0., np.where(scaled >= np.uint64(2**bits), 1., probability))


def acceptance_from_log_ratio(log_ratio, rho, bits=53):
    """Strict input acceptance, including the finite word's zero atom.

    Positive mathematical probabilities below the smallest nonzero word still
    accept word0. This avoids incorrectly losing that atom to exp underflow.
    """
    _validate(rho, bits)
    ratio = np.asarray(log_ratio, float)
    if not np.isfinite(ratio).all():
        raise ValueError("finite log ratios required")
    tiny = ratio <= -bits*math.log(2)
    p = np.exp(np.clip(ratio, -bits*math.log(2), 0.))
    values = word_cdf_lt(p, rho, bits)
    zero_mass = .5*((1+rho)/2)**(bits-1)
    return np.where(ratio >= 0, 1., np.where(tiny, zero_mass, values))


def source_approximation_bounds(rho, retained_bits=53, word_bits=64, source_retained_bits=53, uniforms=None):
    _validate(rho, retained_bits)
    if retained_bits > word_bits or source_retained_bits < retained_bits or source_retained_bits > 53:
        raise ValueError("invalid retained/prefix/word-bit relationship")
    gap = word_bits-retained_bits+1
    tv = abs(rho)**gap/2
    largest_prefix_mass = .5*((1+abs(rho))/2)**(retained_bits-1)
    return {"rho": rho, "prefix_bits": retained_bits, "actual_retained_bits": source_retained_bits,
        "word_bits": word_bits, "between_retained_prefix_transitions": gap,
        "conditional_next_word_TV_upper": tv,
        "full_path_independent_word_TV_union_upper": min(1., max(0, uniforms-1)*tv) if uniforms is not None else None,
        "numeric_prefix_truncation_upper": 2.**-retained_bits-2.**-source_retained_bits,
        "CDF_prefix_bracket_max_width": largest_prefix_mass,
        "source_precision": "Actual53-MSBsource; prefixtruncationandword-independenceidealizationareseparateapproximations",
        "caution": "TVboundsarelawbounds,notautomaticallyunbounded-energyexpectationboundsorcryptographicguarantees"}


def build_bitword_transition(log_density, *, proposal=1.5, rho=.2, bits=53,
                             proposal_rho=None, acceptance_rho=None, bound=8., dx=.02):
    """Generic operator with exact numerical-word CDFs, independent-word idealization."""
    _validate(rho, bits)
    rp = rho if proposal_rho is None else proposal_rho
    ra = rho if acceptance_rho is None else acceptance_rho
    _validate(rp, bits); _validate(ra, bits)
    intervals = round(2*bound/dx)
    if proposal <= 0 or bound <= 0 or dx <= 0 or abs(intervals*dx-2*bound) > 1e-10:
        raise ValueError("positive proposal/domain/dx and aligned domain required")
    grid = np.linspace(-bound, bound, intervals+1)
    target = np.broadcast_to(np.asarray(log_density(grid), float), grid.shape)
    if not np.isfinite(target).all():
        raise ValueError("finite target values required")
    healthy = np.exp(target-target.max()); healthy /= healthy.sum()
    offsets = np.arange(-math.ceil(proposal/dx+.5), math.ceil(proposal/dx+.5)+1)
    left = np.maximum(0., ((offsets-.5)*dx/proposal+1)/2)
    right = np.minimum(1., ((offsets+.5)*dx/proposal+1)/2)
    masses = word_cdf_lt(right, rp, bits)-word_cdf_lt(left, rp, bits)
    keep = masses > 0
    offsets, masses = offsets[keep], masses[keep]
    n = len(grid)
    rows = np.broadcast_to(np.arange(n)[:, None], (n, len(offsets)))
    cols = rows+offsets[None, :]
    valid = (cols >= 0) & (cols < n) & (offsets[None, :] != 0)
    ratio = target[np.clip(cols, 0, n-1)]-target[:, None]
    values = acceptance_from_log_ratio(ratio, ra, bits)*masses[None, :]
    matrix = csr_matrix((values[valid], (rows[valid], cols[valid])), shape=(n, n))
    row_mass = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(row_mass > 1+1e-12):
        raise ArithmeticError("Invalid row mass")
    matrix += csr_matrix((np.maximum(0., 1-row_mass), (np.arange(n), np.arange(n))), shape=(n, n))
    matrix.sort_indices(); digest = hashlib.sha256()
    for a in (grid.astype("<f8"), matrix.indptr.astype("<i8"), matrix.indices.astype("<i8"), matrix.data.astype("<f8")):
        digest.update(a.tobytes())
    # No numerical acceptance quadrature is used: finite-word CDF is explicit.
    return Transition(grid, matrix, healthy, target, proposal, rho, bound, dx, 0, digest.hexdigest())


def save_bitword_reference(output):
    import csv
    import json
    from pathlib import Path
    import resource
    import time
    from rnglab.v04.references import trial_local_energy, trial_log_density
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    roles = {"both": (.2, .2), "proposal_only": (.2, 0.), "acceptance_only": (0., .2), "healthy": (0., 0.)}
    config = {"version": VERSION, "target": {"alpha": 1., "trial_quartic": .2, "coupling": .2},
        "proposal": 1.5, "rho": .2, "word_bits": 64, "actual_retained_bits": 53,
        "prefix_checks": [16, 24, 53], "roles": roles,
        "grids": [[6., .04], [6., .02], [6., .01], [8., .01]],
        "scope": "Marginal-CDF exact for suppliedfiniteword; operator idealizeswords asindependent; no actualsource simulations, no fittedenergy outcomes",
        "pilot_awareness": "Root suppliedpilotmeanbeforethisreference; independentlyderivedpost-pilot/pre-confirmationprediction, notblindedpilotprediction",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    path = output / "config.json"
    if path.exists() and json.loads(path.read_text()) != json.loads(json.dumps(config)):
        raise ValueError("Frozen bitword reference identity changed")
    path.write_text(json.dumps(config, indent=2))
    start = time.perf_counter(); results = []; fresh = reused = 0
    target = lambda x: trial_log_density(x, 1., .2)
    observations = {"energy": lambda x: trial_local_energy(x, 1., .2, .2)}
    for role, (rp, ra) in roles.items():
        bits_list = [16, 24, 53] if role == "both" else [53]
        for bits in bits_list:
            for bound, dx in config["grids"]:
                identity = f"{role}_bits{bits}_b{bound:g}_dx{dx:g}"
                checkpoint = output / "operators" / (identity+".json")
                if checkpoint.exists():
                    result = json.loads(checkpoint.read_text()); reused += 1
                else:
                    operator = build_bitword_transition(target, proposal=1.5, rho=.2, bits=bits,
                        proposal_rho=rp, acceptance_rho=ra, bound=bound, dx=dx)
                    result = {"identity": identity, "role": role, "cdf_bits": bits,
                        "proposal_bit_rho": rp, "acceptance_bit_rho": ra,
                        **stationary_moments(operator, observations),
                        "word_independence_and_prefix_bounds": source_approximation_bounds(.2, bits, uniforms=147456)}
                    checkpoint.parent.mkdir(parents=True, exist_ok=True)
                    temporary = checkpoint.with_suffix(".tmp"); temporary.write_text(json.dumps(result)); temporary.replace(checkpoint)
                    fresh += 1
                results.append(result)
        print("Independent bitword references:", role, "completed", len(results), flush=True)
    primary = next(r for r in results if r["role"] == "both" and r["cdf_bits"] == 53 and r["bound"] == 8 and r["dx"] == .01)
    fine = [r for r in results if r["bound"] == 8 and r["dx"] == .01]
    primary_file = {"status": "Unfitted pre-fresh-confirmation independent-word prediction",
        "source_model": "StationaryfairbinaryMarkov,rho.2; MSBfirst64wordretains53bits; adjacentwordsproposal/acceptance",
        "sampler": "NonGaussianVMC alpha1,b.2,lambda.2,boxhalfwidth1.5",
        "source_active_schedule": "Rawcopybitfaultconstantfrominitialsource; burn8192retained65536; predictedstationarymomentafterburn",
        "reference_checkpoint": "operators/both_bits53_b8_dx0.01.json", "prediction": primary,
        "exact_word_moments": exact_word_moments(.2, 53),
        "exact_CDF_method": "BinaryprefixrecursionforstrictCDF,53branches; zeroacceptanceatomretainedfromlogratio",
        "independent_word_idealization": source_approximation_bounds(.2, 53, uniforms=147456),
        "finite_word_zero_atom": exact_word_moments(.2, 53)["zero_word_probability"],
        "caution": "Roleworddependenceisnotexactlyzero. TVunionboundisnotanunbounded-energyerrorbound. Freshactualpipelineexperimentfalsifiesidealizedpredictionwithinprospectivebroadband."}
    (output / "PRIMARY_PREDICTION.json").write_text(json.dumps(primary_file, indent=2))
    comparisons = [{"identity": r["identity"], "energy_minus_primary": r["moments"]["energy"]-primary["moments"]["energy"],
        "x2_minus_primary": r["moments"]["x2"]-primary["moments"]["x2"]} for r in results if r["role"] == "both"]
    summary = {"reference_calculations": len(results), "fresh_this_invocation": fresh, "reused_checkpoints": reused,
        "independent_random_trials": 0, "wall_seconds": time.perf_counter()-start,
        "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "primary_prediction": primary_file, "finest_by_role_prefix": fine,
        "prefix_grid_domain_comparisons": comparisons,
        "maximum_stationarity_residual": max(r["stationarity_l1_residual"] for r in results)}
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    fields = ["identity", "role", "cdf_bits", "proposal_bit_rho", "acceptance_bit_rho", "bound", "dx", "x", "x2", "x4", "energy", "stationarity_l1_residual"]
    with (output / "predictions.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields); writer.writeheader()
        for r in results:
            row = {key: r[key] for key in fields if key in r}
            row.update({key: r["moments"][key] for key in fields if key in r["moments"]}); writer.writerow(row)
    return summary


if __name__ == "__main__":
    import argparse
    from pathlib import Path
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[3] / "results/v0.4/references/bitword")
    args = parser.parse_args(); result = save_bitword_reference(args.output)
    print("Primaryindependentwordenergy", result["primary_prediction"]["prediction"]["moments"]["energy"], "seconds", result["wall_seconds"])
