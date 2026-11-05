"""Starting-path digest and caller-supplied provenance in run_simulation.

v0.18 Phase 2 addition (gap G1 in research/v018-api-validation/CODEX_WIP_REVIEW.md),
written by Claude Code on top of the Codex run options.
"""
import hashlib
import math

import numpy as np
import pytest

from rnglab14.models.pimc_harmonic import PathIntegralHarmonicOscillator
from rnglab14.sources import PCG64Source, StochasticChannels

PATH = [-.6, .1, .4, .9]


def sources(seed=18101):
    return StochasticChannels(PCG64Source(seed), PCG64Source(seed+1),
                              PCG64Source(seed+2), PCG64Source(seed+3))


def digest(path):
    return hashlib.sha256(np.asarray(path, dtype="<f8").tobytes()).hexdigest()


def run(**options):
    return PathIntegralHarmonicOscillator(beads=4).run_simulation(sources(), sweeps=64, burn_in=16, **options)


@pytest.mark.parametrize("initialization", ["uniform", "exact_gaussian"])
def test_every_mode_records_the_digest_of_the_path_it_started_from(initialization):
    result = run(initialization=initialization, record_trace=True)
    protocol = result.initialization_protocol
    assert protocol["initial_path_sha256"] == digest(result.trace["initial_path"])
    assert protocol["initial_path_digest_scheme"] == "SHA-256 of little-endian float64 bytes in bead order"
    assert protocol["initial_path_provenance"] is None


@pytest.mark.parametrize("supplied", [
    PATH, np.array(PATH), np.array(PATH, dtype=">f8"),
    np.array([-.6, 0., .1, 0., .4, 0., .9, 0.])[::2],
])
def test_explicit_path_digest_ignores_container_byte_order_and_strides(supplied):
    result = run(initial_path=supplied)
    assert result.initialization_protocol["initial_path_sha256"] == digest(PATH)
    assert result.initialization_protocol["initial_path_provenance"] is None


def test_digest_distinguishes_one_ulp():
    nudged = list(PATH)
    nudged[2] = math.nextafter(nudged[2], 1.)
    assert run(initial_path=nudged).initialization_protocol["initial_path_sha256"] != digest(PATH)


def test_provenance_is_recorded_and_isolated_from_later_caller_mutation():
    provenance = {"namespace": "rnglab-v018-api-validation/local/init/0", "draw": {"paths": 2, "kept": 0}}
    result = run(initial_path=PATH, initial_path_provenance=provenance)
    recorded = result.initialization_protocol["initial_path_provenance"]
    assert recorded == provenance and recorded is not provenance
    provenance["namespace"] = "changed"
    provenance["draw"]["kept"] = 1
    assert recorded == {"namespace": "rnglab-v018-api-validation/local/init/0", "draw": {"paths": 2, "kept": 0}}
    assert result.initialization_protocol["mode"] == "explicit_path"


@pytest.mark.parametrize("option", [
    {"initial_path_provenance": {"origin": "x"}},
    {"initial_path": PATH, "initial_path_provenance": ["origin"]},
    {"initial_path": PATH, "initial_path_provenance": {1: "origin"}},
    {"initial_path": PATH, "initial_path_provenance": {"origin": object()}},
    {"initial_path": PATH, "initial_path_provenance": {"origin": math.nan}},
    {"initial_path": PATH, "initial_path_provenance": {"origin": [math.inf]}},
])
def test_invalid_provenance_fails_before_any_role_consumption(option):
    channel = sources()
    centroid = PCG64Source(19100)
    model = PathIntegralHarmonicOscillator(beads=4)
    with pytest.raises(ValueError):
        model.run_simulation(channel, sweeps=64, burn_in=16, centroid_source=centroid, **option)
    assert [s.consumed_bytes for s in (channel.proposal, channel.acceptance,
        channel.initialization, channel.resampling, centroid)] == [0]*5
    assert np.array_equal(model.path, np.zeros(4))
