import json
import pytest

from rnglab_crypto.experiments.config import CampaignConfig
from rnglab_crypto.experiments.campaign import run_campaign
from rnglab_crypto.reporting import freeze_run, finish_run, verify


@pytest.fixture(scope="module")
def campaign():
    return run_campaign(CampaignConfig.smoke())[0]


def test_synthetic_campaign_replay(campaign):
    assert run_campaign(CampaignConfig.smoke())[0] == campaign


def test_timeline_semantics_and_downstream_blind_spot(campaign):
    for rep in range(4):
        rows = {r["condition"]: r for r in campaign["nonces"]["rows"] if r["replicate"] == rep and r["algorithm"] == "AES-128-GCM"}
        assert rows["fixed"]["first_policy_violation"] == 17
        assert rows["counter_reset"]["first_policy_violation"] == 16
        assert rows["rollback"]["first_policy_violation"] == 16
        assert rows["counter_persistent"]["first_policy_violation"] is None
        assert rows["fixed"]["first_source_alarm"] == rows["healthy"]["first_source_alarm"]
        for name in ("fixed", "counter_reset"):
            assert not rows[name]["raw_source_fault_present"]
            assert rows[name]["source_detection_delay_operations"] is None
            assert not rows[name]["source_gate_prevention_with_eligible_source_warning"]
            assert rows[name]["prefault_or_healthy_alarm"] == (rows[name]["first_source_alarm"] is not None)
        for row in rows.values():
            assert row["first_policy_violation"] == row["first_demonstrated_confidentiality_compromise"]
            assert row["first_policy_violation"] == row["first_demonstrated_integrity_forgery"]
            assert row["guard_repeated_emissions"] == 0
            assert row["guard_blocked_operations"] == row["repeated_operations"]


def test_known_state_reseed_and_historical_damage(campaign):
    for row in campaign["recovery"]:
        assert all(v for v in row.values() if type(v) is bool)
    assert all(r["recovered"] for r in campaign["weak_keys"])
    assert all(r["attacker_exact_prediction"] for r in campaign["conditioning"])


def test_manifest_detects_tampering(tmp_path, campaign):
    out = freeze_run(tmp_path/"new", CampaignConfig.smoke())
    finish_run(out, campaign, {"wall_seconds": 0}, make_plots=False)
    assert verify(out)["verified_artifacts"] == 5
    (out/"results.json").write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        verify(out)
    with pytest.raises(FileExistsError):
        freeze_run(out, CampaignConfig.smoke())


@pytest.mark.parametrize("kwargs", [{"synthetic_only":False}, {"repetitions":True}, {"operations":4},
                                   {"fault_onset":256}, {"weak_seed_bits":[32]},
                                   {"weak_seed_bits":[16], "repetitions":512}])
def test_invalid_or_unbounded_configs(kwargs):
    with pytest.raises(ValueError):
        CampaignConfig(**kwargs)
