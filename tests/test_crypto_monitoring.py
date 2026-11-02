import math
import pytest

from rnglab_crypto.monitoring import diagnostics, source_scores, health_screen, EntropyGate, SeedEvidence
from rnglab_crypto.risk import survival_summary, birthday_probability, event_rate


def test_bias_and_periodic_diagnostics_have_distinct_roles():
    assert diagnostics(bytes(1024))["rejected"]
    assert diagnostics(b"\xaa"*1024)["rejected"]  # balanced marginal, excessive runs
    assert health_screen(bytes(32))["latched_failure"]
    assert source_scores([bytes(32)])[0] == 16


def test_entropy_gate_never_credits_length_or_health_alone():
    gate = EntropyGate()
    with pytest.raises(RuntimeError):
        gate.check(SeedEvidence(256, "unvalidated"))
    with pytest.raises(RuntimeError):
        gate.check(SeedEvidence(256, "synthetic_model"))
    gate.check(SeedEvidence(256, "synthetic_model"), recovery=True)
    assert not gate.quarantined
    with pytest.raises(RuntimeError):
        gate.check(SeedEvidence(256, "synthetic_model"), healthy=False)


def test_censoring_is_not_a_failure_or_a_tail_median():
    result = survival_summary([0, None, None, None], 10)
    assert result["censored"] == 3
    assert result["median_operations"] is None
    assert result["curve"][-1]["survival"] == .75
    assert event_rate([False]*20)["interval"][1] > 0
    with pytest.raises(ValueError):
        survival_summary([10], 10)


def test_exact_birthday_risk():
    assert birthday_probability(2, 8) == pytest.approx(1/256)
    assert birthday_probability(3, 1) == 1
    assert birthday_probability(0, 96) == 0
    assert birthday_probability(1000, 96) == pytest.approx(1000*999/(2*2**96))


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf"), True])
def test_invalid_entropy_evidence(value):
    with pytest.raises(ValueError):
        SeedEvidence(value, "synthetic_model")
