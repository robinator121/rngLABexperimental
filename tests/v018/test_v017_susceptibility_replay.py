"""Deterministic replay of the v0.17 healthy susceptibility counterexample.

Regenerates the v0.17 CALIBRATE healthy confirmation cell ising_L4_T0 (L=4,
T=1.5, 512 chains, 2048 sweeps, burn-in 512) from its original seed and
compares with v0.17 SUSCEPTIBILITY_PROBE.json row ising_L4_T0 (off-repo file,
SHA-256 b2953cff893416557939885d0e9e7a48f021c7c4e0f27c72c82570d188d9d5db).
This is a replay of saved results, not new confirmation; protocol
research/v018-estimators/protocol.json, data.regression_replay. Uses about
0.8 GB of memory for a few seconds.
"""
import hashlib
import math

import numpy as np

from rnglab import susceptibility
from rnglab16.core import ising


def test_v017_ising_L4_T0_counterexample_replays_exactly():
    seed = int.from_bytes(hashlib.sha256(b"CALIBRATE-20261022-v1/ising_L4_T0/confirmation").digest()[:8], "big")
    trace, _, _ = ising(4, 1.5, 512, 2048, seed, "healthy", 0.02)
    m = trace[:, 512:, 1]
    plugin = susceptibility.plugin(m, 16, 1.5)
    symmetric = susceptibility.symmetric(m, 16, 1.5)["estimate"]
    assert math.isclose(plugin.mean(), 0.7458277684119011, rel_tol=1e-9)
    assert math.isclose(symmetric.mean(), 10.403679741753471, rel_tol=1e-9)
    assert int(np.sum(np.all(np.sign(m[:, 1:]) * np.sign(m[:, :-1]) >= 0, axis=1))) == 470
