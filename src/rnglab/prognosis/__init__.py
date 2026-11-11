"""Conditional forecasting of the remaining useful run of a Monte Carlo estimate under a degrading input.

Every forecast here is conditional on a declared law of future degradation
(`families`). `validity` maps a severity path to the probability that a
batch-means interval covers the reference; `surrogate` is a Gaussian workload
with known laws for checking methods; `forecasters` holds baselines and the
declared-model posterior forecaster; `scoring` holds proper scores,
reliability, and failure-time scores with explicit right-censoring.

docs/v019/FORECASTING_METHODS.md states what can and cannot be forecast and
reports the methods validation study in research/v019-forecast-methods/.
"""
from . import families, forecasters, scoring, surrogate, validity

__all__ = ["families", "forecasters", "scoring", "surrogate", "validity"]
