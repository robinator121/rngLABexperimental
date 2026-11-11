"""Methods validation study for conditional forecasting on the Gaussian surrogate.

    PYTHONPATH=src python research/v019-forecast-methods/run_study.py            # the frozen study, once
    PYTHONPATH=src python research/v019-forecast-methods/run_study.py --smoke    # engineering run, own namespace

Everything is read from protocol.json (frozen; its SHA-256 is in protocol.sha256
and is checked before any draw). The smoke run uses the smoke namespace and
reduced sizes, writes under smoke/, and is excluded from every result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np

from rnglab.prognosis import families as fam
from rnglab.prognosis import forecasters as fc
from rnglab.prognosis import scoring as sc
from rnglab.prognosis import surrogate as sur
from rnglab.prognosis import validity as val
from rnglab.seeds import namespace_seed

HERE = Path(__file__).resolve().parent
SMOKE_SIZES = {"training": 600, "evaluation_per_scenario": 300, "onsets": 13, "rates": 12, "forward": 5, "resamples": 200}


def rng_for(namespace: str, *parts: str) -> np.random.Generator:
    return np.random.default_rng(namespace_seed(namespace, *parts))


def draw_paths(name: str, config: dict, spec: sur.SurrogateSpec, prior: fam.RampPrior, size: int,
               rng: np.random.Generator):
    family = config["family"]
    if family == "ramp" and config.get("prior") == "declared":
        return prior.sample(rng, size)
    low, high = config["t0"]
    t0 = rng.integers(low, high + 1, size=size)
    if family == "ramp":
        a, b = config["t_cap"]
        return fam.LinearRamp(t0, np.exp(rng.uniform(np.log(a), np.log(b), size=size)), spec.s_max)
    if family == "step_hold":
        a, b = config["level_log_uniform"]
        return fam.StepHold(t0, np.exp(rng.uniform(np.log(a), np.log(b), size=size)))
    if family == "random_walk":
        return fam.random_walk_paths(rng, size, spec.blocks, spec.block, t0_low=low, t0_high=high,
                                     drift=config["drift_per_sweep"], volatility=config["volatility_per_sqrt_sweep"],
                                     s_max=spec.s_max)
    raise ValueError(f"unknown family in scenario {name}")


def issue(runs: sur.Trajectories, threshold: float):
    alarm = sur.first_alarm(runs.statistic, threshold)
    keep = np.flatnonzero(alarm >= 0)
    return alarm, keep, [fc.history_at(runs.statistic[i], runs.spec.ends, alarm[i]) for i in keep]


def subset_path(path, rows: np.ndarray):
    """The severity paths of the selected runs."""
    if isinstance(path, fam.LinearRamp):
        return fam.LinearRamp(np.asarray(path.t0)[rows], np.asarray(path.t_cap)[rows], path.s_max)
    if isinstance(path, fam.RampThenHold):
        return fam.RampThenHold(np.asarray(path.t0)[rows], np.asarray(path.t_cap)[rows], path.s_max, np.asarray(path.freeze)[rows])
    if isinstance(path, fam.StepHold):
        return fam.StepHold(np.asarray(path.t0)[rows], np.asarray(path.level)[rows])
    return fam.SampledPath(path.values[rows], path.block)


def calibrated(q: np.ndarray, y: np.ndarray) -> dict:
    z = sc.spiegelhalter_z(q, y)
    gap = sc.calibration_gap(sc.reliability(q, y), minimum_count=50)
    ok = z is not None and abs(z) < 3 and gap["weighted_gap"] is not None and gap["weighted_gap"] < 0.05
    return {"spiegelhalter_z": z, **gap, "calibrated": bool(ok)}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(argv)
    protocol_bytes = (HERE / "protocol.json").read_bytes()
    protocol = json.loads(protocol_bytes)
    digest = hashlib.sha256(protocol_bytes).hexdigest()
    recorded = (HERE / "protocol.sha256").read_text().split()[0] if (HERE / "protocol.sha256").is_file() else None
    if not args.smoke and digest != recorded:
        raise SystemExit("protocol.json does not match protocol.sha256: nothing generated")
    namespace = protocol["smoke_namespace"] if args.smoke else protocol["namespace"]
    out_dir = HERE / ("smoke" if args.smoke else "results")
    if not args.smoke and (out_dir / "results.json").exists():
        raise SystemExit("results/results.json exists: the frozen study runs once")
    sizes = {"training": protocol["sets"]["training"], "evaluation_per_scenario": protocol["sets"]["evaluation_per_scenario"],
             "onsets": protocol["posterior_grid"]["onsets"], "rates": protocol["posterior_grid"]["rates"],
             "forward": protocol["forward_simulations_per_hypothesis"], "resamples": protocol["bootstrap"]["resamples"]}
    if args.smoke:
        sizes.update(SMOKE_SIZES)
    started = time.time()
    spec = sur.SurrogateSpec(**{**protocol["surrogate"], "lengths": tuple(protocol["surrogate"]["lengths"])})
    prior = fam.RampPrior(**protocol["declared_prior"])
    level = protocol["bootstrap"]["level"]

    # --- calibration and training: everything a forecaster learns is fixed here
    healthy = sur.simulate(spec, None, protocol["sets"]["calibration_healthy"], rng_for(namespace, "calibration", "noise"))
    threshold = sur.calibrate_threshold(healthy)
    sigma_hat = fc.estimate_sigma_eff(healthy)
    train = sur.simulate(spec, prior.sample(rng_for(namespace, "training", "paths"), sizes["training"]),
                         sizes["training"], rng_for(namespace, "training", "noise"))
    beta_hat, kappa_hat = fc.estimate_beta(train), fc.estimate_kappa(train)
    _, keep_train, hist_train = issue(train, threshold)
    model = val.ValidityModel(beta_hat, sigma_hat, spec.batches, spec.level)
    bayes = fc.BayesRamp(spec, prior, sizes["onsets"], sizes["rates"]).fit(beta_hat, sigma_hat, kappa_hat)
    bayes_plain = fc.BayesRamp(spec, prior, sizes["onsets"], sizes["rates"], name="bayes-ramp-no-widening",
                               drift_widening=False).fit(beta_hat, sigma_hat, kappa_hat)
    persistence = fc.Persistence(spec, model, kappa_hat)
    fitted = {}
    for length in spec.lengths:
        chosen = np.array([h.alarm_sweep < length for h in hist_train])
        histories = [h for h, c in zip(hist_train, chosen) if c]
        y = train.valid(length)[keep_train][chosen]
        fitted[length] = {"n": int(chosen.sum()), "valid_fraction": float(y.mean()),
                          "base-rate": fc.BaseRate().fit(histories, {length: y}),
                          "alarm-time": fc.LogisticForecaster("alarm").fit(histories, {length: y}),
                          "lite-logistic": fc.LogisticForecaster("lite").fit(histories, {length: y})}
    bayes.fit_failure_times(rng_for(namespace, "training", "forward-simulation"), simulations=sizes["forward"])
    time_train, observed_train = train.failure
    empirical = fc.EmpiricalFailure().fit(hist_train, time_train[keep_train], observed_train[keep_train], spec.ends)
    results = {"schema": "rnglab-v019-forecast-methods-results-v1", "namespace": namespace, "smoke": args.smoke,
               "protocol_sha256": digest, "sizes": sizes,
               "environment": {"python": platform.python_version(), "numpy": np.__version__,
                               "scipy": __import__("scipy").__version__},
               "calibration": {"threshold": threshold, "sigma_eff_hat": sigma_hat,
                               "healthy_valid_fraction": {str(L): float(healthy.valid(L).mean()) for L in spec.lengths}},
               "training": {"beta_hat": beta_hat, "kappa_hat": kappa_hat, "alarms": len(hist_train),
                            "eligible": {str(L): fitted[L]["n"] for L in spec.lengths},
                            "valid_fraction": {str(L): fitted[L]["valid_fraction"] for L in spec.lengths}},
               "scenarios": {}}

    # --- evaluation: each scenario is opened once
    stored_forecasts = {}
    for name, config in protocol["scenarios"].items():
        size = sizes["evaluation_per_scenario"]
        if name == "E_freeze":
            source = stored_forecasts["E_continue"]
            freeze = np.where(source["alarm"] >= 0, spec.ends[np.maximum(source["alarm"], 0)], spec.horizon)
            path = fam.RampThenHold(source["path"].t0, source["path"].t_cap, spec.s_max, freeze)
            runs = sur.simulate(spec, path, size, rng_for(namespace, "E_continue", "noise"))  # the same noise
        else:
            path = draw_paths(name, config, spec, prior, size, rng_for(namespace, name, "paths"))
            runs = sur.simulate(spec, path, size, rng_for(namespace, name, "noise"))
        alarm, keep, histories = issue(runs, threshold)
        onset = np.asarray(getattr(path, "t0", np.full(size, -1)))
        row = {"runs": size, "alarms": len(keep), "no_alarm": int(size - len(keep)),
               "pre_onset_alarms": (int(np.sum(spec.ends[alarm[keep]] <= onset[keep])) if hasattr(path, "t0") else None),
               "lengths": {}}
        q_bayes, q_plain, q_persist = bayes.predict(histories), bayes_plain.predict(histories), persistence.predict(histories)
        q_oracle = fc.Oracle(model, spec.lengths).predict_paths(subset_path(path, keep))
        true_model = val.ValidityModel(spec.beta, spec.sigma_eff, spec.batches, spec.level)
        p_true = {L: true_model.probability(subset_path(path, keep), L) for L in spec.lengths}
        kept = {"alarm": alarm, "path": path, "keep": keep, "forecasts": {}, "p_true": p_true}
        for length in spec.lengths:
            chosen = np.array([h.alarm_sweep < length for h in histories], dtype=bool)
            selected = [h for h, c in zip(histories, chosen) if c]
            y = runs.valid(length)[keep][chosen]
            forecasts = {"base-rate": fitted[length]["base-rate"].predict(selected)[length],
                         "alarm-time": fitted[length]["alarm-time"].predict(selected)[length],
                         "lite-logistic": fitted[length]["lite-logistic"].predict(selected)[length],
                         "persistence": q_persist[length][chosen], "bayes-ramp": q_bayes[length][chosen],
                         "bayes-ramp-no-widening": q_plain[length][chosen], "oracle": q_oracle[length][chosen]}
            kept["forecasts"][length] = {"chosen": chosen, **forecasts}
            cell = {"eligible": int(chosen.sum()), "late_alarm": int(len(chosen) - chosen.sum()),
                    "valid_fraction": float(y.mean()) if len(y) else None, "forecasters": {}}
            if len(y) >= 50:
                baseline = sc.brier(forecasts["alarm-time"], y)
                for label, q in forecasts.items():
                    seed = namespace_seed(namespace, name, "analysis", f"skill-{label}-L{length}")
                    cell["forecasters"][label] = {
                        "brier": sc.brier(q, y), "log_score": sc.log_score(q, y)[0], "log_score_clipped": sc.log_score(q, y)[1],
                        "skill_over_alarm_time": sc.paired_bootstrap_skill(y, q, forecasts["alarm-time"], seed,
                                                                         sizes["resamples"], level),
                        "calibration": calibrated(q, y), "reliability": sc.reliability(q, y)}
                for first, second in (("bayes-ramp", "lite-logistic"), ("bayes-ramp", "bayes-ramp-no-widening"),
                                      ("lite-logistic", "alarm-time")):
                    seed = namespace_seed(namespace, name, "analysis", f"difference-{first}-{second}-L{length}")
                    cell[f"brier_difference:{first}-minus-{second}"] = sc.paired_bootstrap_difference(
                        y, forecasts[first], forecasts[second], seed, sizes["resamples"], level)
                cell["alarm_time_brier"] = baseline
            row["lengths"][str(length)] = cell
        # failure time
        time_, observed = runs.failure
        column = sc.outcome_column(time_[keep], observed[keep], spec.ends)
        pmf = {"bayes-ramp": bayes.predict_failure(histories), "empirical": empirical.predict_failure(histories)}
        rps = {k: sc.ranked_probability_score(v, column) for k, v in pmf.items()}
        boot = np.random.default_rng(namespace_seed(namespace, name, "analysis", "rps-skill"))
        values = []
        for _ in range(sizes["resamples"]):
            index = boot.integers(0, len(column), size=len(column))
            values.append(1 - rps["bayes-ramp"][index].mean() / rps["empirical"][index].mean())
        alarm_sweep = spec.ends[alarm[keep]]
        remaining = np.maximum(time_[keep] - alarm_sweep, 0)
        row["failure_time"] = {
            "runs": len(column), "censored": int(np.sum(~observed[keep])),
            "already_failed_at_issuance": int(np.sum(observed[keep] & (time_[keep] <= alarm_sweep))),
            "km_median_failure_sweep": sc.km_median(time_[keep], observed[keep]),
            "km_median_remaining_sweeps": sc.km_median(remaining, observed[keep]),
            "rps": {k: float(v.mean()) for k, v in rps.items()},
            "rps_skill": 1 - float(rps["bayes-ramp"].mean() / rps["empirical"].mean()),
            "rps_skill_interval": np.quantile(values, [(1 - level) / 2, 1 - (1 - level) / 2]).tolist(),
            "interval_coverage": {k: sc.interval_coverage(v, column, protocol["failure_time"]["interval_level"]) for k, v in pmf.items()},
            "median_error": {k: sc.median_error(v, column, spec.ends) for k, v in pmf.items()}}
        results["scenarios"][name] = row
        if name.startswith("E_"):
            stored_forecasts[name] = kept

    # --- R6: the identical-history pair
    a, b = stored_forecasts["E_continue"], stored_forecasts["E_freeze"]
    pair = {"alarms_identical": bool(np.array_equal(a["alarm"], b["alarm"])), "lengths": {}}
    for length in spec.lengths:
        fa, fb = a["forecasts"][length], b["forecasts"][length]
        labels = [k for k in fa if k not in ("chosen", "oracle")]
        same = bool(np.array_equal(fa["chosen"], fb["chosen"])
                    and all(np.array_equal(fa[k], fb[k]) for k in labels))
        pa, pb = a["p_true"][length][fa["chosen"]], b["p_true"][length][fb["chosen"]]
        bound = float(np.mean((pa - pb)**2) / 4) if same else None
        excess = {k: {"continue": float(np.mean((fa[k] - pa)**2)), "freeze": float(np.mean((fb[k] - pb)**2))} for k in labels} if same else None
        pair["lengths"][str(length)] = {
            "forecasts_identical": same, "bound_worst_world_excess": bound, "excess_over_path_oracle": excess,
            "mean_p_continue": float(pa.mean()), "mean_p_freeze": float(pb.mean()),
            "bound_holds_for_every_forecaster": (None if not same else
                                                 bool(all(max(v.values()) >= bound - 1e-12 for v in excess.values())))}
    results["identifiability_pair"] = pair
    results["seconds"] = time.time() - started
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "results.json").write_text(json.dumps(results, indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(f"wrote {out_dir / 'results.json'} in {results['seconds']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
