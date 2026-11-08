#!/usr/bin/env python3
"""Post-hoc check of the small refresh-condition bias (NOT pre-registered).

RNGLab v0.18 (written by Claude Code, 2026-11-06), after the frozen analysis
showed CV and TD bias intervals excluding zero in the refresh condition. For
each condition and estimator it reports the mean bias, its t statistic, and the
correlation between the CV and TD replica means. This is descriptive and
exploratory; it changes no gate or verdict.

    PYTHONPATH=src python research/v018-api-validation/posthoc_bias_check.py

Reads results/<condition>_replicas.csv; writes results/posthoc_bias_check.json.
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from rnglab14.models.pimc_harmonic import exact_finite_bead_energy  # noqa: E402


def main():
    reference = exact_finite_bead_energy(16, 2.0)
    out = {"status": "post hoc, not pre-registered; descriptive only", "reference_E_M": reference, "conditions": {}}
    for condition in ("local", "refresh", "iid"):
        with (HERE / "results" / f"{condition}_replicas.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        entry = {}
        for name in ("cv", "td"):
            d = np.array([float(r[f"{name}_mean"]) for r in rows]) - reference
            entry[f"{name}_bias"] = {"mean": float(d.mean()), "t": float(d.mean() / (d.std(ddof=1) / math.sqrt(len(d))))}
        cv = np.array([float(r["cv_mean"]) for r in rows])
        td = np.array([float(r["td_mean"]) for r in rows])
        entry["corr_cv_td_means"] = float(np.corrcoef(cv, td)[0, 1])
        out["conditions"][condition] = entry
    out["bias_intervals_reported"] = 6
    out["bias_intervals_excluding_zero"] = sum(abs(e[f"{n}_bias"]["t"]) > 1.96
                                               for e in out["conditions"].values() for n in ("cv", "td"))
    (HERE / "results" / "posthoc_bias_check.json").write_text(json.dumps(out, indent=1) + "\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
