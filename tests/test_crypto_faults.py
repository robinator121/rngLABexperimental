import numpy as np
import pytest

from rnglab_crypto.sources import ReplaySource
from rnglab_crypto.faults import FaultAdapter, FaultSpec, EntropyUnavailable, PLACEMENTS


def adapter(kind, **kwargs):
    return FaultAdapter(ReplaySource(bytes(32)), FaultSpec(kind, **kwargs), lab_only=True)


@pytest.mark.parametrize("kind", list(PLACEMENTS))
def test_fault_schedule_replay_and_shadow_recovery(kind):
    a = adapter(kind, onset=1, duration=2)
    b = adapter(kind, onset=1, duration=2)
    healthy = ReplaySource(bytes(32))
    for op in range(5):
        expected = healthy.read_bytes(32)
        if kind == "starvation" and op in (1, 2):
            with pytest.raises(EntropyUnavailable):
                a.read_bytes(32)
            with pytest.raises(EntropyUnavailable):
                b.read_bytes(32)
        else:
            aa, bb = a.read_bytes(32), b.read_bytes(32)
            assert aa == bb
            if op not in (1, 2):
                assert aa == expected
    assert a.operation == b.operation == 5


@pytest.mark.parametrize("kwargs", [{"intensity":float("nan")}, {"intensity":1.01}, {"duration":0},
                                  {"onset":True}, {"exposed_bits":9}, {"recovery":"silent"}])
def test_invalid_fault(kwargs):
    with pytest.raises(ValueError):
        FaultSpec("bias", **kwargs)


def test_fault_requires_opt_in():
    with pytest.raises(ValueError):
        FaultAdapter(ReplaySource(bytes(32)), FaultSpec("bias"))


def test_bit_bias_and_markov_mechanisms():
    bits = np.unpackbits(np.frombuffer(adapter("bias", intensity=.6).read_bytes(10000), dtype=np.uint8))
    assert .79 < bits.mean() < .81
    bits = np.unpackbits(np.frombuffer(adapter("serial", intensity=.8).read_bytes(10000), dtype=np.uint8))
    assert .48 < bits.mean() < .52
    assert .89 < (bits[1:] == bits[:-1]).mean() < .91


def test_periodic_and_partial_exposure_are_explicit():
    assert adapter("periodic", period_bytes=2).read_bytes(20)[:2]*10 == adapter("periodic", period_bytes=2).read_bytes(20)
    source = adapter("partial_exposure", entropy_bits=12, exposed_bits=8)
    value = int.from_bytes(source.read_bytes(32), "big")
    assert source.leak == {"known_prefix": value >> 4, "unknown_suffix_bits": 4}


def test_rollback_repeats_initial_drbg_output():
    source = adapter("rollback", onset=2)
    first = source.read_bytes(32)
    source.read_bytes(32)
    assert source.read_bytes(32) == first


def test_zero_intensity_is_paired_healthy_control():
    fault = adapter("bias", intensity=0)
    healthy = ReplaySource(bytes(32))
    assert fault.read_bytes(100) == healthy.read_bytes(100)
