"""Stored-output verifier and campaign recount (research/v019-audit/verify.py).

The campaign recount is exercised on hand-built per-run views with a reduced
plan; its numbers are compared with direct recomputations written in the test.
The per-run verifier is exercised on records assembled from the independent
replay, then tampered with: every tampering must be reported.
"""
import copy
import importlib.util as _ilu
import json
import math
import sys as _sys
from pathlib import Path as _Path

import numpy as np
import pytest

_spec = _ilu.spec_from_file_location("_audit_loader", _Path(__file__).with_name("_audit_loader.py"))
au = _sys.modules.get("_audit_loader") or _ilu.module_from_spec(_spec)
if "_audit_loader" not in _sys.modules:
    _sys.modules["_audit_loader"] = au
    _spec.loader.exec_module(au)

from rnglab12.uncertainty import summarize_series

vf = au.load("v019_audit_verify", au.REPO / "research/v019-audit/verify.py")
rc, rp = vf.rc, vf.rp
NS = au.AUDIT_NAMESPACE
PLAN = {"calibration": ("healthy", 39), "healthy": ("healthy", 60), "fault/SHIFT": ("SHIFT", 80),
        "fault/LAG1": ("LAG1", 80), "fault/DUP": ("DUP", 80)}
FIT, EVAL = range(0, 48), range(48, 80)


def make_view(condition, family, replica, rng):
    """A synthetic run: noise statistics, a ramp in the target after onset, a biased trace after onset."""
    reference = rc.reference_energy()
    statistic = rng.exponential(0.6, size=159)
    cv = reference + rng.normal(scale=0.4, size=rc.H)
    t0 = t_cap = s_max = target = None
    if family != "healthy":
        t0, t_cap, s_max = int(rng.integers(2048, 8193)), float(np.exp(rng.uniform(np.log(2048), np.log(163840)))), 0.5
        severity = rc.severity_path(t0, t_cap, s_max)
        window_mean = np.array([severity[e - 512:e].mean() for e in rc.BATTERY_ENDS])
        target = rng.normal(size=159) + 400.0 * window_mean
        statistic = np.maximum(statistic, 0.5 * target**2 / math.log(10) / 2)
        cv = cv + 0.08 * severity
    return vf.RunView(condition=condition, family=family, replica=replica, t0=t0, t_cap=t_cap, s_max=s_max,
                      statistic=statistic, statistic_prime=statistic * 0.9, target=target,
                      mean_d=rng.normal(scale=0.15, size=157), uphill_fraction=0.65 + rng.normal(scale=0.003, size=157), cv=cv)


@pytest.fixture(scope="module")
def campaign():
    rng = np.random.default_rng(2026)
    views = {c: [make_view(c, family, i, rng) for i in range(count)] for c, (family, count) in PLAN.items()}
    return views, vf.recount_campaign(views, NS, plan=PLAN, fitting=FIT, evaluation=EVAL, resamples=300)


def test_campaign_recount_thresholds_and_healthy_gates(campaign):
    views, out = campaign
    assert out["status"] == "RECOUNTED" and out["frozen_plan"] is False and not any(out["missing"].values())
    maxima = sorted(v.statistic.max() for v in views["calibration"])
    assert out["thresholds"]["order_statistic"] == 38 and out["thresholds"]["full"] == maxima[37]
    assert out["thresholds"]["prime"] == pytest.approx(0.9 * maxima[37])
    alarms = sum(v.statistic.max() > maxima[37] for v in views["healthy"])
    assert out["V2_D"]["alarms"] == alarms and out["V2_D"]["runs"] == 60 and out["V2_D"]["unknown"] == 0
    reference = rc.reference_energy()
    for length in rc.LENGTHS:
        covered = sum(summarize_series(v.cv[:length], reference=reference, batch_sizes=(length // 16,))
                      ["batch_means"][str(length // 16)]["covered"] for v in views["healthy"])
        assert out["criterion_1"]["per_length"][length]["covered"] == covered
        mcse = [summarize_series(v.cv[:length], reference=reference, batch_sizes=(length // 16,))
                ["batch_means"][str(length // 16)]["mcse"] for v in views["calibration"]]
        assert out["sigma_eff"][length] == pytest.approx(math.sqrt(length) * np.mean(mcse), rel=1e-12)
    pooled = np.concatenate([v.mean_d for v in views["calibration"]])
    assert out["application_constants"]["sigma_D"] == pytest.approx(pooled.std(ddof=1), rel=1e-13)
    assert set(out["short_interval_coverage"]) == set(rc.BATTERY_ENDS[::16].tolist())


def test_campaign_recount_issuance_and_scores(campaign):
    views, out = campaign
    tau = out["thresholds"]["full"]
    for condition, (family, count) in PLAN.items():
        if family == "healthy":
            continue
        row = out["families"][family]
        assert sum(row["classes"].values()) == count
        alarms = {v.replica: next((int(e) for e, s in zip(rc.BATTERY_ENDS, v.statistic) if s > tau), None)
                  for v in views[condition]}
        for length in rc.LENGTHS:
            cell = row["per_length"][length]
            eligible = [r for r, a in alarms.items() if a is not None and a < length]
            assert cell["eligible"] == len(eligible)
            assert cell["eligible_fitting"] == sum(r in FIT for r in eligible)
            assert cell["eligible_evaluation"] == sum(r in EVAL for r in eligible)
            if cell["eligible_fitting"] and cell["eligible_evaluation"]:
                reference = rc.reference_energy()
                valid = {v.replica: bool(rc.endpoint(v.cv, length, reference)["valid"]) for v in views[condition]}
                q0 = np.mean([valid[r] for r in eligible if r in FIT])
                assert cell["q0"] == pytest.approx(q0)
                y = np.array([valid[r] for r in sorted(eligible) if r in EVAL], dtype=float)
                assert cell["scores"]["models"]["F0"]["brier"] == pytest.approx(np.mean((q0 - y)**2))
                assert cell["scores"]["models"]["F0+"]["skill"] in (0.0, None)
                assert set(cell["scores"]["models"]) == set(rc.MODELS)
                assert sum(b["count"] for b in cell["calibration"]["F-B-lite"]) == len(y)
        assert row["oracle_beta"] == pytest.approx(0.08, rel=0.25)  # the synthetic bias slope
    # eighty runs per family cannot meet the frozen minimums of 300 and 240
    assert out["criteria_2_3"]["criterion_2"] == out["criteria_2_3"]["criterion_3"] == "INSUFFICIENT_INFORMATION"
    assert out["criterion_4"]["result"] == "INSUFFICIENT_INFORMATION"
    assert out["verdicts"]["criterion_5"] == "INSUFFICIENT_INFORMATION"


def test_campaign_recount_missing_and_unknown_inputs(campaign):
    views, _ = campaign
    short = {**views, "calibration": views["calibration"][:-1]}
    out = vf.recount_campaign(short, NS, plan=PLAN, fitting=FIT, evaluation=EVAL, resamples=50)
    assert out["status"] == "INSUFFICIENT_INFORMATION" and out["missing"]["calibration"] == [38]
    broken = copy.copy(views["healthy"][0])
    broken.statistic = broken.statistic.copy()
    broken.statistic[4] = np.nan
    out = vf.recount_campaign({**views, "healthy": [broken] + views["healthy"][1:]}, NS, plan=PLAN,
                              fitting=FIT, evaluation=EVAL, resamples=50)
    assert out["V2_D"]["unknown"] == 1 and not out["issuance_interpretable"]
    assert out["verdicts"]["criterion_2"] == out["verdicts"]["criterion_4"] == "UNINTERPRETABLE"
    assert out["verdicts"]["criterion_1"] in ("PASS", "FAIL", "INSUFFICIENT_INFORMATION")  # reported separately


def test_fit_and_predict_fallback_rules():
    constants = {"sigma_D": 0.15, "mu_f": 0.65, "sigma_f": 0.003}
    rng = np.random.default_rng(5)
    view = make_view("fault/LAG1", "LAG1", 0, rng)
    early = vf.features(view, 768, constants)
    assert early["application_index"] is None and early["slope_count"] == 2 and early["A1_over_sigma_D"] is None
    assert vf.features(view, 512, constants)["target_slope"] == 0.0
    later = vf.features(view, 4096, constants)
    assert later["application_index"] == 12 and later["slope_count"] == 4
    assert later["A1_over_sigma_D"] == pytest.approx(view.mean_d[12] / 0.15)
    assert rc.APPLICATION_ENDS[12] == 4096  # last application window ends at the alarm: nothing later is read
    rows = [vf.features(view, int(a), constants) for a in rng.choice(rc.BATTERY_ENDS[4:100], size=60)]
    y = [int(r["target_at_alarm"] > np.median([q["target_at_alarm"] for q in rows])) for r in rows]
    gap = dict(rows[0], target_at_alarm=None)
    no_app = dict(rows[1], application_index=None, A1_over_sigma_D=None, A2_standardised=None)
    half_app = dict(rows[2], A2_standardised=None)
    out = vf.fit_and_predict(rows, y, [rows[3], gap, no_app, half_app])
    q0 = sum(y) / len(y)
    assert out["predictions"]["F0"] == [q0] * 4
    assert out["predictions"]["F-B-lite"][1] == q0 and out["fallbacks"]["F-B-lite"] == {"MISSING_FEATURE": 1}
    assert out["predictions"]["F-AB-lite"][2] == out["predictions"]["F-B-lite"][2]
    assert out["predictions"]["F-AB-lite"][3] == out["predictions"]["F-B-lite"][3]
    assert out["fallbacks"]["F-AB-lite"] == {"MISSING_FEATURE": 1, "NO_APPLICATION_CHECKPOINT": 1,
                                             "MISSING_APPLICATION_FEATURE": 1}
    assert out["predictions"]["F0+"][1] != q0  # the alarm sweep is present: F0+ still predicts
    one_class = vf.fit_and_predict(rows, [1] * len(rows), rows[:2])
    assert one_class["predictions"]["F-B-lite"] == [1.0, 1.0] and one_class["models"]["F0+"]["fallback"] == "ONE_CLASS"
    assert vf.fit_and_predict([], [], rows[:1])["q0"] is None
    saved = {"F0+": {"theta": [0.0, 0.0], "mu": [0.0], "scale": [1.0]}, "F-B-lite": {"fallback": "NOT_CONVERGED"},
             "F-AB-lite": {"fallback": "ONE_CLASS"}}
    from_saved = vf.fit_and_predict(rows, y, rows[:1], saved)   # saved parameters are used; nothing is refitted
    assert from_saved["predictions"]["F0+"] == [0.5] and from_saved["predictions"]["F-B-lite"] == [q0]
    assert from_saved["reasons"]["F-B-lite"] == ["NOT_CONVERGED"] and from_saved["models"]["F0+"]["source"] == "saved"
    with pytest.raises(rc.RecountError):
        vf.fit_and_predict(rows, y, rows[:1], {"F0+": saved["F0+"]})


def assemble_record(namespace, condition, family, replica):
    """A stored-format record built from the independent replay; battery p-values are placeholders."""
    run = rp.replay(namespace, condition, family, replica)
    arrays = {"CV": run["CV"], "TD": run["TD"], "marks": run["marks"], "severity": run["severity"]}
    rng = np.random.default_rng(replica)
    battery = []
    for end in rc.BATTERY_ENDS:
        if end == 512:  # the one checkpoint the tests regenerate: real statistics of the emitted window
            lo, hi = rc.window_words(run["marks"], 0, 512)
            real = rc.battery_checkpoint(run["emitted_words"][lo:hi])
            battery.append({"end": 512, **{k: real[k] for k in ("p_values", "S", "S_prime", "available_p_count",
                                                                 "frequency_signed_z", "lag1_feature", "duplicate_rate",
                                                                 "original_bits", "was_cut", "pair_count")}})
            continue
        p = rng.random(10).tolist()
        s, count = rc.alarm_statistic(p)
        battery.append({"end": int(end), "p_values": p, "S": s, "S_prime": rc.alarm_statistic(p[:8])[0],
                        "available_p_count": count})
    difference, uphill = run["TD"] - run["CV"], (np.diff(run["marks"], prepend=0) - 128) // 8
    application = [{"end": int(end), "mean_D": float(np.mean(difference[end - 1024:end])),
                    "uphill_count": int(uphill[end - 1024:end].sum()),
                    "uphill_fraction": int(uphill[end - 1024:end].sum()) / 16384} for end in rc.APPLICATION_ENDS]
    reference = rc.reference_energy()
    endpoints = {}
    for length in rc.LENGTHS:
        e = rc.endpoint(run["CV"], length, reference)
        endpoints[str(length)] = {**{k: e[k] for k in ("valid", "estimate", "signed_error", "mcse", "half_width",
                                                       "interval_95", "covered", "batch_means")},
                                  "invalid_reason": e["reason"]}
    profile = run["profile"] or {}
    record = {"namespace": namespace, "condition": condition, "family": family, "replica": replica,
              "seeds": run["seeds"], "status": "SUCCESS", "local_bytes": run["local_bytes"], "B0": run["B0"],
              "profile": {k: (None if profile.get(k) is None else (int(profile[k]) if k == "t0" else float(profile[k])))
                          for k in ("t0", "T_cap", "r", "s_max")},
              "emitted_sha256": run["emitted_sha256"], "base_sha256": run["base_sha256"],
              "prefix_sha256": run["prefix_sha256"], "start_sha256": run["start_sha256"],
              "array_hashes": {k: {"sha256": vf.array_sha256(a, vf.ARRAYS[k])} for k, a in arrays.items()},
              "conventional": battery, "application": application, "endpoints": endpoints}
    return record, arrays


@pytest.fixture(scope="module")
def stored(tmp_path_factory):
    root = tmp_path_factory.mktemp("audit_store")
    for condition, (family, _) in vf.PLAN.items():
        record, arrays = assemble_record(NS, condition, family, 0)
        record_path, trace_path = vf.stored_paths(root, condition, 0)
        record_path.parent.mkdir(parents=True)
        trace_path.parent.mkdir(parents=True)
        record_path.write_text(json.dumps(record, sort_keys=True, allow_nan=False))
        np.savez_compressed(trace_path, **arrays)
    return root


def test_verifier_accepts_consistent_runs_and_reports_every_tampering(stored):
    record, arrays = vf.load_stored(stored, "fault/DUP", 0)
    clean = vf.verify_run(record, arrays, NS, deep=True, checkpoint_stride=1000)
    assert clean["issues"] == [] and clean["deep"] and clean["compared"] > 1000

    def fields(rec=record, arr=arrays, namespace=NS, deep=False):
        return {i["field"] for i in vf.verify_run(rec, arr, namespace, deep=deep, checkpoint_stride=1000)["issues"]}

    # A run filed under the wrong namespace. Not the smoke namespace: the check would draw its profile.
    assert "seeds" in fields(namespace=NS + "-other")
    tampered = copy.deepcopy(record)
    tampered["profile"]["T_cap"] = float(np.nextafter(tampered["profile"]["T_cap"], np.inf))
    assert "profile" in fields(tampered)
    tampered = copy.deepcopy(record)
    tampered["B0"] += 8
    assert {"B0"} <= fields(tampered) and "replay.B0" in fields(tampered, deep=True)
    tampered = copy.deepcopy(record)
    tampered["conventional"][7]["S"] = tampered["conventional"][7]["S"] + 1e-9
    assert "battery[2304].S_from_stored_p" in fields(tampered)
    tampered = copy.deepcopy(record)  # a p-value changed together with its S: only the regenerated stream shows it
    tampered["conventional"][0]["p_values"][0] = 0.5
    tampered["conventional"][0]["S"], tampered["conventional"][0]["available_p_count"] = rc.alarm_statistic(
        tampered["conventional"][0]["p_values"])
    tampered["conventional"][0]["S_prime"] = rc.alarm_statistic(tampered["conventional"][0]["p_values"][:8])[0]
    assert "battery[512].p.frequency" not in fields(tampered) and "battery[512].p.frequency" in fields(tampered, deep=True)
    tampered = copy.deepcopy(record)
    tampered["conventional"].pop()
    assert "battery.ends" in fields(tampered)
    tampered = copy.deepcopy(record)
    tampered["endpoints"]["10240"]["valid"] = not tampered["endpoints"]["10240"]["valid"]
    assert "endpoint[10240].valid" in fields(tampered)
    tampered = copy.deepcopy(record)
    tampered["application"][3]["uphill_count"] += 1
    assert "application[1792].uphill_count" in fields(tampered)
    for key, index in (("CV", 100), ("TD", 40959), ("marks", 5), ("severity", 30000)):
        changed = dict(arrays)
        changed[key] = arrays[key].copy()
        changed[key][index] = arrays[key][index] + 8 if key == "marks" else np.nextafter(arrays[key][index], 9.0)
        found = fields(arr=changed, deep=True)
        assert f"{key}.sha256" in found, key
        assert (f"replay.{key}.bitwise" in found) or (key == "severity" and "severity.bitwise" in found), key
    truncated = {k: (v[:-1] if k == "CV" else v) for k, v in arrays.items()}
    assert "CV.shape" in fields(arr=truncated)
    assert "CV.dtype" in fields(arr={**arrays, "CV": arrays["CV"].astype("<f4")})


def test_command_line_report(stored, tmp_path, capsys):
    out = tmp_path / "verification.json"
    code = vf.main(["--root", str(stored), "--namespace", NS, "--replicas", "1", "--deep", "all",
                    "--checkpoint-stride", "1000", "--no-campaign", "--out", str(out)])
    report = json.loads(out.read_text())
    assert code == 0 and report["agrees"] and report["discrepancy_count"] == 0
    assert report["runs_planned"] == 5 and report["runs_deep"] == 5 and report["runs_missing"] == []
    assert report["protocol_sha256"] == (au.REPO / "research/v019-d-pilot/protocol.sha256").read_text().split()[0]
    assert set(report["code_sha256"]) == {"recount.py", "replay.py", "verify.py", "integrity.py"}
    # The audit's own engineering store has no attempt log, and an audit namespace does not need one.
    assert report["attempt_log_check"] == {"status": "ABSENT", "required": False, "compared": 0, "discrepancies": [], "agrees": True}
    assert "5 planned runs, 0 missing, 5 regenerated" in capsys.readouterr().out
    vf.stored_paths(stored, "healthy", 0)[1].rename(tmp_path / "moved.npz")
    try:
        assert vf.main(["--root", str(stored), "--namespace", NS, "--replicas", "1", "--deep", "none",
                        "--no-campaign", "--out", str(out)]) == 1
        report = json.loads(out.read_text())
        assert report["runs_missing"] == ["healthy/00000"] and not report["agrees"]
        assert report["campaign"]["status"] == "NOT_RUN"
    finally:
        (tmp_path / "moved.npz").rename(vf.stored_paths(stored, "healthy", 0)[1])


def test_i4_sample_matches_frozen_recipe():
    sample = vf.deep_sample("rnglab-v019-d-pilot")
    assert set(sample) == set(vf.PLAN) and all(len(v) == 5 == len(set(v)) for v in sample.values())
    assert all(max(v) < vf.PLAN[c][1] for c, v in sample.items())
    assert sample == vf.deep_sample("rnglab-v019-d-pilot") != vf.deep_sample("rnglab-v019-d-pilot-smoke")
