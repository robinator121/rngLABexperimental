#!/usr/bin/env python3
"""Print the VMC masking contrast quoted in README.md.

RNGLab v0.18 audit helper (written by Claude Code, 2026-11-05). Runs exactly
the configuration of tests/test_models14.py::
test_vmc_zero_variance_masking_contrast_with_imperfect_wavefunction (proposal
bias_delta=0.25 on PCG64 seed 42, acceptance seed 99, initialization seed 1,
delta=1.0, 3000 steps, 500 burn-in) and records the values the test only
bounds. This is a demonstration of the textbook zero-variance property of VMC
at the exact trial function, not a new hypothesis. One deterministic run; no
uncertainty is attached to these numbers.

Usage: PYTHONPATH=src python scripts/v018/vmc_masking_values.py > out.json
"""
import json

from rnglab14.models import HarmonicOscillatorVMC
from rnglab14.sources import DegradedFaultSource, PCG64Source, StochasticChannels


def run(alpha: float) -> dict:
    channels = StochasticChannels(
        proposal=DegradedFaultSource(PCG64Source(seed=42), bias_delta=0.25),
        acceptance=PCG64Source(seed=99),
        initialization=PCG64Source(seed=1),
    )
    result = HarmonicOscillatorVMC(alpha=alpha, delta=1.0).run_simulation(channels, steps=3000, burn_in=500)
    keys = ("energy_mean", "energy_std", "energy_bias", "x2_mean", "x2_bias", "acceptance_rate")
    return {key: getattr(result, key) for key in keys}


if __name__ == "__main__":
    print(json.dumps({"script": "scripts/v018/vmc_masking_values.py",
                      "configuration": "tests/test_models14.py masking contrast, single run",
                      "alpha_0.5_exact": run(0.5), "alpha_0.6_imperfect": run(0.6)}, indent=2))
