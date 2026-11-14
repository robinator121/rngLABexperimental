"""Stage D pilot analysis stages on a synthetic campaign, checked against the independent recount.

The campaign is 207 synthetic runs (tests/v019/_pilot_fixture.py) on a reduced plan; it is not
Monte Carlo output and carries audit-only labels. Every stage of the pipeline runs once on it.
The independent recount (research/v019-audit) shares no code with the pipeline; agreement between
the two is the main check here, followed by tampering that each side must report.

Needs the measurement layer (RNGLAB_DPILOT_MEASURE or research/v019-d-pilot/measure.py).
"""
import copy
import dataclasses
import importlib.util as _ilu
import json
import math
import shutil
import sys as _sys
from pathlib import Path as _Path

import numpy as np
import pytest
import scipy.optimize
import scipy.stats

_spec = _ilu.spec_from_file_location("_pilot_fixture", _Path(__file__).with_name("_pilot_fixture.py"))
fx = _sys.modules.get("_pilot_fixture") or _ilu.module_from_spec(_spec)
if "_pilot_fixture" not in _sys.modules:
    _sys.modules["_pilot_fixture"] = fx
    _spec.loader.exec_module(fx)

au, pc, pg, pa, pf = fx.au, fx.pc, fx.pg, fx.pa, fx.pf
vf = au.load("v019_audit_verify", au.REPO / "research/v019-audit/verify.py")
rc = vf.rc
PROTOCOL = json.loads((au.REPO / "research/v019-d-pilot/protocol.json").read_text())
ORACLE_FIELDS = {"oracle.unavailable", "oracle.brier"}


class Recording:
    """The measurement layer, with every run it is asked to open written down."""

    def __init__(self, measure):
        self._measure, self.opened = measure, []

    def __getattr__(self, name):
        return getattr(self._measure, name)

    def load_run(self, root, condition, replica, **options):
        self.opened.append((condition, replica))
        return self._measure.load_run(root, condition, replica, **options)


@pytest.fixture(scope="module")
def campaign(tmp_path_factory):
    ctx = fx.context(fx.TEST_PLAN, tmp_path_factory.mktemp("campaign"))
    pg.generate(ctx, engine_spec=fx.ENGINE_SPEC, workers=0)
    pa.calibrate(ctx)
    recorder = Recording(ctx.measure)
    pf.fit(dataclasses.replace(ctx, measure=recorder))
    pa.issue(ctx)
    pf.evaluate(ctx)
    pa.healthy_summary(ctx)
    pa.fault_summary(ctx)
    pa.implementation_checks(ctx, pg.replay_check(ctx, fx.SyntheticEngine(ctx), {"healthy": [0], "fault/DUP": [1]}))
    pa.criteria(ctx)
    return ctx, recorder.opened


def recount(ctx, root=None):
    plan = {c: (f, n) for c, f, n in ctx.plan.conditions}
    views = {c: [vf.load_view(ctx.root, c, i) for i in range(n)] for c, (_, n) in plan.items()}
    saved = vf.saved_models_from_primary(pc.read_json(ctx.result("models.json")))
    independent = vf.recount_campaign(views, ctx.plan.namespace, plan=plan, fitting=ctx.plan.fitting,
                                      evaluation=ctx.plan.evaluation, saved_models=saved)
    return independent, vf.compare_with_primary(root or ctx.root, independent)


@pytest.fixture(scope="module")
def recounted(campaign):
    return recount(campaign[0])


def rows_of(ctx, name):
    return [json.loads(line) for line in ctx.result(name).read_text().splitlines()]


# --------------------------------------------------------------------------- outputs and the independent recount
def test_every_result_file_exists_and_is_canonical(campaign):
    ctx, _ = campaign
    for name in ("calibration.json", "models.json", "model_manifest.json", "healthy.json", "faults.json", "scores.json",
                 "implementation_checks.json", "criteria.json", "run_manifest.json"):
        raw = ctx.result(name).read_bytes()
        assert raw == pc.canonical_json(json.loads(raw)), name
        assert b"NaN" not in raw and b"Infinity" not in raw
    for name in ("issuance.jsonl", "predictions.jsonl", "failures.jsonl"):
        raw = ctx.result(name).read_bytes()
        assert raw == pc.canonical_jsonl([json.loads(line) for line in raw.decode().splitlines()]), name
    listed = {entry["path"].split("/")[-1] for entry in PROTOCOL["outputs"]["files"]}
    produced = {p.name for p in ctx.root.glob("results/*") if p.is_file()}
    assert {"calibration.json", "issuance.jsonl", "models.json", "model_manifest.json", "predictions.jsonl", "healthy.json",
            "faults.json", "scores.json", "implementation_checks.json", "criteria.json", "run_manifest.json",
            "failures.jsonl"} <= produced & listed


def test_independent_recount_agrees_with_every_result_file(recounted):
    independent, comparison = recounted
    assert independent["status"] == "RECOUNTED" and comparison["files_absent"] == [] and comparison["compared"] > 1900
    other = [d for d in comparison["discrepancies"] if d["field"] not in ORACLE_FIELDS]
    assert other == []
    # Under the frozen oracle rule the two sides may disagree on which runs have no oracle value: the
    # library returns NaN or a number depending on the last bits of the noncentrality (F-01).
    assert len(comparison["discrepancies"]) <= 9


@pytest.mark.parametrize("name,change,expected", [
    ("calibration.json", lambda d: d.update(threshold_full=d["threshold_full"] + 1e-9), "threshold_full"),
    ("calibration.json", lambda d: d["sigma_eff"].update({"20480": d["sigma_eff"]["20480"] * (1 + 1e-9)}), "sigma_eff[20480]"),
    ("healthy.json", lambda d: d["coverage"]["10240"].update(covered=d["coverage"]["10240"]["covered"] + 1), "coverage[10240].covered"),
    ("healthy.json", lambda d: d["V2_D"].update(result="PASS"), "V2_D.result"),
    ("models.json", lambda d: d["families"]["DUP"]["20480"].update(q0=d["families"]["DUP"]["20480"]["q0"] + 1e-9), "q0"),
    ("scores.json", lambda d: d["families"]["LAG1"]["20480"]["models"]["F-B-lite"].update(
        brier=d["families"]["LAG1"]["20480"]["models"]["F-B-lite"]["brier"] + 1e-9), "F-B-lite.brier"),
    ("scores.json", lambda d: d["families"]["SHIFT"]["10240"]["models"]["F0+"]["calibration"][2].update(count=99), "F0+.calibration"),
    ("criteria.json", lambda d: d["results_before_precedence"].update({"1": "PASS"}), "1"),
    ("criteria.json", lambda d: d["results"].update({"4": "PASS"}), "final 4"),
])
def test_recount_reports_a_changed_result_file(campaign, recounted, tmp_path, name, change, expected):
    ctx, _ = campaign
    shutil.copytree(ctx.root / "results", tmp_path / "results", ignore=shutil.ignore_patterns("runs"))
    data = pc.read_json(tmp_path / "results" / name)
    change(data)
    (tmp_path / "results" / name).write_bytes(pc.canonical_json(data))
    comparison = vf.compare_with_primary(tmp_path, recounted[0])
    assert expected in {d["field"] for d in comparison["discrepancies"]} and not comparison["agrees"]


def test_recount_reports_changed_rows_and_missing_files(campaign, recounted, tmp_path):
    ctx, _ = campaign
    shutil.copytree(ctx.root / "results", tmp_path / "results", ignore=shutil.ignore_patterns("runs"))
    issuance = rows_of(ctx, "issuance.jsonl")
    row = next(r for r in issuance if r["family"] == "LAG1" and r["full"]["status"] == "ALARM")
    row["full"]["alarm_sweep"] += 256
    (tmp_path / "results/issuance.jsonl").write_bytes(pc.canonical_jsonl(issuance))
    predictions = rows_of(ctx, "predictions.jsonl")
    target = next(p for p in predictions if p["eligible"] and p["family"] == "DUP" and p["L"] == 20480)
    target["predictions"]["F-AB-lite"]["q"] += 1e-9
    (tmp_path / "results/predictions.jsonl").write_bytes(pc.canonical_jsonl(predictions))
    (tmp_path / "results/scores.json").unlink()
    comparison = vf.compare_with_primary(tmp_path, recounted[0])
    found = {(d["id"], d["field"]) for d in comparison["discrepancies"]}
    assert (f"issuance LAG1/{row['replica']:05d}", "alarm_sweep") in found
    assert ("predictions DUP L=20480", "q[F-AB-lite]") in found
    assert comparison["files_absent"] == ["scores.json"] and not comparison["agrees"]


# --------------------------------------------------------------------------- calibration
def test_threshold_is_the_prescribed_order_statistic(campaign):
    ctx, _ = campaign
    calibration = pc.read_json(ctx.result("calibration.json"))
    plan = ctx.plan
    assert calibration["status"] == "COMPLETE" and calibration["calibration_runs"] == 39 == len(calibration["maxima_full"])
    ordered = sorted(calibration["maxima_full"])
    tau = calibration["threshold_full"]
    assert tau == ordered[plan.threshold_order - 1] == ctx.h1lib.threshold(calibration["maxima_full"], order_statistic=38)
    assert tau != ordered[plan.threshold_order - 2] and tau != ordered[plan.threshold_order]      # neither neighbour
    assert sum(m > tau for m in ordered) == 39 - 38 and calibration["ties_at_threshold_full"] == 1
    defining = calibration["maxima_full"].index(tau)
    record, _ = pa.load_run(ctx, "calibration", defining)
    assert ctx.measure.issuance(record["conventional"], tau)["status"] == "NO_ALARM"          # equality is not an alarm
    assert ctx.measure.issuance(record["conventional"], float(np.nextafter(tau, 0)))["status"] == "ALARM"
    records = [pa.load_run(ctx, "calibration", i)[0] for i in range(39)]
    assert calibration["application_constants"] == ctx.h1lib.application_constants(records)
    for length in pc.LENGTHS:
        mcse = [r["endpoints"][str(length)]["mcse"] for r in records]
        assert calibration["sigma_eff"][str(length)] == pytest.approx(math.sqrt(length) * np.mean(mcse), rel=1e-14)
    assert set(calibration["input_record_sha256"]) == {f"calibration/{i:05d}" for i in range(39)}
    with pytest.raises(pc.PilotError, match="not overwritten"):
        pa.calibrate(ctx)                                    # thresholds are fixed once


def test_incomplete_calibration_makes_everything_downstream_unavailable(tmp_path):
    plan = pc.Plan(fx.NAMESPACE, (("calibration", "healthy", 3), ("healthy", "healthy", 2), ("fault/SHIFT", "SHIFT", 4),
                                  ("fault/LAG1", "LAG1", 4), ("fault/DUP", "DUP", 4)), range(0, 2), range(2, 4), 3, 3, "test")
    ctx = fx.context(plan, tmp_path)
    out = pg.generate(ctx, engine_spec=fx.PARTLY_FAILING_SPEC, workers=0)
    assert out["counts"]["FAILED_RECORD"] == 1 and out["counts"]["COMPLETE"] == 16
    calibration = pa.calibrate(ctx)
    assert calibration["status"] == "INSUFFICIENT_INFORMATION" and calibration["threshold_full"] is None
    assert calibration["maxima_full"][1] is None and calibration["application_constants"] is None    # not pooled over fewer runs
    assert all(v is None for v in calibration["sigma_eff"].values())
    models = pf.fit(ctx)
    assert all(cell["q0"] is None and cell["eligible_fitting"] == 0 for cells in models["families"].values() for cell in cells.values())
    rows = pa.issue(ctx)
    assert {r["full"]["status"] for r in rows} == {"UNKNOWN"} and all(r["full"]["eligibility"]["20480"]["eligible"] is None for r in rows)
    scores = pf.evaluate(ctx)
    assert all(cell["status"] == "INSUFFICIENT_INFORMATION" and cell["counts"]["unknown"] == 2
               for cells in scores["families"].values() for cell in cells.values())
    healthy = pa.healthy_summary(ctx)
    assert healthy["V2_D"]["result"] == "INSUFFICIENT_INFORMATION" and healthy["V2_D"]["unknown"] == 2
    pa.fault_summary(ctx)
    pa.implementation_checks(ctx, None)
    criteria = pa.criteria(ctx)
    assert [criteria["results"][k] for k in "234"] == ["UNINTERPRETABLE"] * 3 and criteria["interpretation_blocked"]
    independent, comparison = recount(ctx)
    assert independent["status"] == "INSUFFICIENT_INFORMATION" and comparison["discrepancies"] == []


# --------------------------------------------------------------------------- evaluation barrier and leakage
def test_fit_opens_fitting_runs_only_and_evaluation_is_opened_once(campaign):
    ctx, opened = campaign
    plan = ctx.plan
    assert len(opened) == 3 * len(plan.fitting) and all(c.startswith("fault/") and r in plan.fitting for c, r in opened)
    manifest = pc.read_json(ctx.result("model_manifest.json"))
    paths = [entry["path"] for entry in manifest["fitting_inputs"]]
    assert paths[0] == "results/calibration.json"
    assert not [p for p in paths[1:] if int(_Path(p).stem) not in plan.fitting or "/fault/" not in p]
    assert manifest["models_sha256"] == pc.sha256_file(ctx.result("models.json"))
    models = pc.read_json(ctx.result("models.json"))
    for cells in models["families"].values():
        for cell in cells.values():
            assert set(cell["fitting_ids"]) <= set(plan.fitting)
            assert all(set(m["excluded_ids"]) <= set(cell["fitting_ids"]) for m in cell["models"].values())
    predictions = rows_of(ctx, "predictions.jsonl")
    assert {p["replica"] for p in predictions} == set(plan.evaluation) and len(predictions) == 3 * 3 * len(plan.evaluation)
    assert all(p["models_sha256"] == manifest["models_sha256"] for p in predictions)
    with pytest.raises(pc.PilotError, match="opened once"):
        pf.evaluate(ctx)
    with pytest.raises(pc.PilotError, match="not overwritten"):
        pf.fit(ctx)                                           # no refit after the model file exists
    original = ctx.result("models.json").read_bytes()
    try:
        ctx.result("models.json").write_bytes(original + b" ")
        with pytest.raises(pc.PilotError, match="does not match its manifest"):
            pf.evaluate(ctx)
    finally:
        ctx.result("models.json").write_bytes(original)


def test_fit_refuses_a_run_that_is_not_a_fitting_run(tmp_path):
    class Leaky(pc.Plan):
        def split(self, condition, replica):
            return "evaluation"

    plan = Leaky(fx.NAMESPACE, fx.SMALL_PLAN.conditions, range(0, 1), range(1, 2), 2, 2, "test")
    ctx = fx.context(plan, tmp_path)
    ctx.result("calibration.json").parent.mkdir(parents=True)
    ctx.result("calibration.json").write_bytes(pc.canonical_json({"application_constants": None, "threshold_full": 1.0,
                                                                  "sigma_eff": {}}))
    with pytest.raises(pc.PilotError, match="evaluation barrier"):
        pf.fit(ctx)
    assert not ctx.result("models.json").exists()


def test_analysis_seeds_are_distinct_from_each_other_and_from_run_seeds():
    from rnglab.seeds import namespace_seed
    bootstrap = {namespace_seed(pc.PILOT_NAMESPACE, "fault/" + family, "analysis", f"bootstrap-L{length:05d}")
                 for family in pc.FAMILIES for length in pc.LENGTHS}
    assert len(bootstrap) == 9
    dlib = au.dlib()
    run_seeds = {s for condition, _, replica in pc.FROZEN_PLAN.runs() for s in dlib.seeds_for(pc.PILOT_NAMESPACE, condition, replica).values()}
    assert len(run_seeds) == 5 * 5999 and not bootstrap & run_seeds
    assert PROTOCOL["data"]["seeds"]["bootstrap"].startswith("seed=namespace_seed(ns,'fault/'+family,'analysis',f'bootstrap-L{L:05d}')")


# --------------------------------------------------------------------------- timing and horizons
def test_issuance_rows_respect_the_boundaries(campaign):
    ctx, _ = campaign
    rows = [r for r in rows_of(ctx, "issuance.jsonl")]
    assert len(rows) == 24 + 3 * 48 and not [r for r in rows if r["condition"] == "calibration"]
    tau = pc.read_json(ctx.result("calibration.json"))["threshold_full"]
    for r in rows:
        full = r["full"]
        record, arrays = pa.load_run(ctx, r["condition"], r["replica"])
        statistic = [c["S"] for c in record["conventional"]]
        first = next((c["end"] for c in record["conventional"] if c["S"] > tau), None)       # strict, chronological
        assert full["alarm_sweep"] == first and (full["status"] == "ALARM") == (first is not None)
        if first is not None:
            assert full["S_at_alarm"] == statistic[(first - 512) // 256] > tau
            assert full["issuance_endpoint"]["n"] == 16 * (first // 16) == first
            assert full["issuance_endpoint"]["already_invalid"] == (not full["issuance_endpoint"]["valid"])
        for length in pc.LENGTHS:
            assert full["eligibility"][str(length)]["eligible"] == (first is not None and first < length)   # alarm == L is not eligible
        if r["family"] != "healthy" and first is not None:
            assert full["onset_classification"] == ("pre_onset" if first <= r["t0"] else "post_onset")
            assert full["alarm_minus_t0"] == first - r["t0"]
        assert r["valid"] == {str(length): record["endpoints"][str(length)]["valid"] for length in pc.LENGTHS}
    assert {r["full"]["onset_classification"] for r in rows if r["family"] == "healthy"} <= {"no_alarm", None}


def test_features_use_nothing_after_the_alarm(campaign):
    ctx, _ = campaign
    constants = pc.read_json(ctx.result("calibration.json"))["application_constants"]
    record, _ = pa.load_run(ctx, "fault/LAG1", 3)
    alarm = 512 + 256 * 20
    base = pf.features(ctx, record, alarm, constants)
    later = copy.deepcopy(record)
    for checkpoint in later["conventional"][21:]:
        checkpoint["lag1_feature"] = 1e9
    for checkpoint in later["application"]:
        if checkpoint["end"] > alarm:
            checkpoint["mean_D"], checkpoint["uphill_fraction"] = 1e9, 1.0
    assert pf.features(ctx, later, alarm, constants) == base                 # the future does not enter
    now = copy.deepcopy(record)
    now["conventional"][20]["lag1_feature"] += 1.0
    assert pf.features(ctx, now, alarm, constants)["target_at_alarm"] == base["target_at_alarm"] + 1.0
    ends = [c["end"] for c in record["conventional"][17:21]]
    values = [c["lag1_feature"] for c in record["conventional"][17:21]]
    assert base["slope_count"] == 4 and base["target_slope"] == pytest.approx(np.polyfit(ends, values, 1)[0], rel=1e-9)
    assert base["application_end"] == alarm                                  # an application window ending at the alarm counts
    assert base["A1_over_sigma_D"] == pytest.approx(record["application"][(alarm - 1024) // 256]["mean_D"] / constants["sigma_D"])
    early = [pf.features(ctx, record, end, constants) for end in (512, 768, 1024, 1280)]
    assert [f["slope_count"] for f in early] == [1, 2, 3, 4] and early[0]["target_slope"] == 0.0
    assert early[1]["target_slope"] == pytest.approx((record["conventional"][1]["lag1_feature"] - record["conventional"][0]["lag1_feature"]) / 256)
    assert [f["application_reason"] for f in early] == ["NO_APPLICATION_CHECKPOINT", "NO_APPLICATION_CHECKPOINT", None, None]
    assert [f["application_end"] for f in early] == [None, None, 1024, 1280]
    gap = copy.deepcopy(record)
    gap["conventional"][18]["lag1_feature"] = None
    assert pf.features(ctx, gap, alarm, constants)["target_slope"] is None and pf.features(ctx, gap, alarm + 1024, constants)["target_slope"] is not None
    missing = copy.deepcopy(record)
    missing["application"][(alarm - 1024) // 256]["uphill_fraction"] = None
    assert pf.features(ctx, missing, alarm, constants)["application_reason"] == "MISSING_APPLICATION_FEATURE"
    assert pf.features(ctx, record, alarm, None)["application_reason"] == "MISSING_APPLICATION_FEATURE"   # no constants, no feature


# --------------------------------------------------------------------------- forecasters
def synthetic_rows(n=200, seed=1):
    rng = np.random.default_rng(seed)
    rows, labels = [], []
    for _ in range(n):
        alarm = float(rng.choice(np.arange(2048, 12000, 256)))
        target, slope = rng.normal(3, 2), rng.normal(0, 0.01)
        a1, a2 = rng.normal(), rng.normal()
        p = 1 / (1 + math.exp(-(-0.0003 * (alarm - 6000) - 0.5 * (target - 3) + 0.3 * a1)))
        rows.append({"alarm_sweep": alarm, "target_at_alarm": target, "target_slope": slope, "A1_over_sigma_D": a1,
                     "A2_standardised": a2, "application_reason": None})
        labels.append(int(rng.random() < p))
    return rows, labels


def test_the_fit_is_the_frozen_call_and_reaches_the_optimum(monkeypatch):
    call = PROTOCOL["forecasters"]["logistic"]["call"]
    for text in ("x0=np.zeros(d+1,dtype=np.float64)", "method='L-BFGS-B'", "bounds=None",
                 "'gtol':1e-4,'ftol':64*np.finfo(np.float64).eps,'maxiter':100,'maxfun':15000,'maxls':50,'maxcor':10"):
        assert text in call
    assert pf.OPTIMISER == {"method": "L-BFGS-B", "options": {"gtol": 1e-4, "ftol": 64 * float(np.finfo(np.float64).eps),
                                                              "maxiter": 100, "maxfun": 15000, "maxls": 50, "maxcor": 10}}
    assert pf.FEATURES == PROTOCOL["forecasters"]["features_in_order"]
    seen = {}
    real = scipy.optimize.minimize

    def spy(**options):
        seen.update(options)
        return real(**options)

    monkeypatch.setattr(pf.scipy.optimize, "minimize", spy)
    rows, labels = synthetic_rows()
    model = pf.fit_model("F-AB-lite", rows, labels, list(range(len(rows))))
    assert seen["method"] == "L-BFGS-B" and seen["bounds"] is None and seen["options"] == pf.OPTIMISER["options"]
    assert np.array_equal(seen["x0"], np.zeros(6)) and seen["x0"].dtype == np.float64
    assert model["status"] == "FITTED" and model["optimiser"]["terminated_by"] in ("gradient", "relative_reduction")
    x = np.array([[r[c] for c in pf.FEATURES["F-AB-lite"]] for r in rows])
    assert model["mu"] == pytest.approx(x.mean(axis=0)) and model["scale"] == pytest.approx(x.std(axis=0, ddof=0))
    reference = rc.fit_logistic_newton((x - x.mean(axis=0)) / x.std(axis=0), np.array(labels, dtype=float))
    mine = pf.predict({"q0": 0.5, "models": {name: pf.fit_model(name, rows, labels, list(range(len(rows))))
                                             for name in ("F0+", "F-B-lite", "F-AB-lite")}}, rows[0])
    theirs = rc.logistic_predict(reference["theta"], x.mean(axis=0), x.std(axis=0), x[:1])[0]
    assert mine["F-AB-lite"]["q"] == pytest.approx(theirs, abs=2e-4)          # the independent Newton optimum


def test_fallbacks_follow_the_protocol(monkeypatch):
    rows, labels = synthetic_rows(60, seed=2)
    ids = list(range(60))
    assert pf.fit_model("F0+", rows, [1] * 60, ids)["fallback_reason"] == "ONE_CLASS"
    assert pf.fit_model("F0+", [], [], [])["fallback_reason"] == "EMPTY_SUBSET"
    without_application = [{**r, "A1_over_sigma_D": None} for r in rows]
    assert pf.fit_model("F-AB-lite", without_application, labels, ids)["fallback_reason"] == "EMPTY_SUBSET"
    some = [{**r, "target_slope": None} if i < 5 else r for i, r in enumerate(rows)]
    partial = pf.fit_model("F-B-lite", some, labels, ids)
    assert partial["training_rows"] == 55 and partial["excluded_ids"] == [0, 1, 2, 3, 4]      # reported, not hidden
    constant = pf.fit_model("F0+", [{**r, "alarm_sweep": 4096.0} for r in rows], labels, ids)
    assert constant["constant_feature"] == [True] and constant["scale"] == [1.0]

    def no_convergence(**options):
        return scipy.optimize.OptimizeResult(x=np.ones(2), fun=1.0, jac=np.ones(2), success=False, status=1,
                                             message="STOP: TOTAL NO. OF ITERATIONS REACHED LIMIT", nit=100, nfev=120)

    monkeypatch.setattr(pf.scipy.optimize, "minimize", no_convergence)
    failed = pf.fit_model("F0+", rows, labels, ids)
    assert failed["status"] == "FALLBACK_F0" and failed["fallback_reason"] == "NOT_CONVERGED" and failed["w"] is None
    assert failed["optimiser"]["message"].startswith("STOP") and failed["optimiser"]["accepted"] is False

    def raises(**options):
        raise FloatingPointError("overflow")

    monkeypatch.setattr(pf.scipy.optimize, "minimize", raises)
    assert pf.fit_model("F0+", rows, labels, ids)["optimiser"]["exception"]["class"] == "FloatingPointError"


def test_prediction_rules():
    rows, labels = synthetic_rows(120, seed=3)
    ids = list(range(120))
    cell = {"q0": 0.4, "models": {name: pf.fit_model(name, rows, labels, ids) for name in ("F0+", "F-B-lite", "F-AB-lite")}}
    full = pf.predict(cell, rows[0])
    assert full["F0"] == {"q": 0.4, "fallback": None} and all(0 < full[m]["q"] < 1 and full[m]["fallback"] is None for m in pc.MODELS[1:])
    early = pf.predict(cell, {**rows[0], "application_reason": "NO_APPLICATION_CHECKPOINT", "A1_over_sigma_D": None, "A2_standardised": None})
    assert early["F-AB-lite"] == {"q": full["F-B-lite"]["q"], "fallback": "NO_APPLICATION_CHECKPOINT"}      # not zero-filled
    unavailable = pf.predict(cell, {**rows[0], "application_reason": "MISSING_APPLICATION_FEATURE", "A2_standardised": None})
    assert unavailable["F-AB-lite"] == {"q": full["F-B-lite"]["q"], "fallback": "MISSING_APPLICATION_FEATURE"}
    no_target = pf.predict(cell, {**rows[0], "target_at_alarm": None})
    assert no_target["F-B-lite"] == {"q": 0.4, "fallback": "MISSING_FEATURE"} and no_target["F-AB-lite"]["q"] == 0.4
    assert no_target["F0+"]["q"] == full["F0+"]["q"]                         # the alarm sweep is still there
    fallen = {"q0": 0.4, "models": {**cell["models"], "F-B-lite": pf.fit_model("F-B-lite", rows, [0] * 120, ids)}}
    assert pf.predict(fallen, rows[0])["F-B-lite"] == {"q": 0.4, "fallback": "ONE_CLASS"}
    nothing = pf.predict({"q0": None, "models": cell["models"]}, rows[0])
    assert all(v == {"q": None, "fallback": "NO_ELIGIBLE_FITTING_RUNS"} for v in nothing.values())


# --------------------------------------------------------------------------- oracle
def test_oracle_rule_is_frozen_and_the_amendment_changes_only_unavailable_values():
    assert pf.oracle_rule() == "frozen"                      # no approved amendment exists
    missing = pf.oracle_probability(None, 0.4, 0.1, 20480)
    assert missing["available"] is False and missing["reason"] == "ORACLE_PARAMETER_UNAVAILABLE"
    c = scipy.stats.t.ppf(0.975, df=15)
    changed = filled = 0
    for mean_severity in np.linspace(0.0, 0.5, 400):
        frozen = pf.oracle_probability(0.07, 0.4, float(mean_severity), 20480, "frozen")
        amended = pf.oracle_probability(0.07, 0.4, float(mean_severity), 20480, "A01")
        delta = 0.07 * mean_severity / (0.4 / math.sqrt(20480))
        if frozen["available"]:
            expected = float(scipy.stats.nct.cdf(c, df=15, nc=delta) - scipy.stats.nct.cdf(-c, df=15, nc=delta))
            assert frozen["q"] == expected == amended["q"] and not amended["lower_term_replaced"]
        else:
            assert frozen["q"] is None and frozen["reason"] == "NONFINITE_OR_OUT_OF_RANGE" and delta > 8.5
            assert amended["available"] and amended["lower_term_replaced"] and 0 <= amended["q"] < 2e-9
            assert amended["q"] == pytest.approx(rc.validity_probability_stable(delta), abs=1e-17)
            filled += 1
        changed += frozen["q"] != amended["q"]
    assert changed == filled
    if filled == 0:
        pytest.skip("this SciPy does not return NaN here")
    negative = pf.oracle_probability(-0.07, 0.4, 0.3, 40960, "A01")       # a negative bias: both terms finite
    assert negative["available"] and not negative["lower_term_replaced"]


# --------------------------------------------------------------------------- scores
def test_scores_equal_a_direct_recomputation(campaign):
    ctx, _ = campaign
    from rnglab.seeds import namespace_seed
    predictions = rows_of(ctx, "predictions.jsonl")
    scores = pc.read_json(ctx.result("scores.json"))
    for family in pc.FAMILIES:
        for length in pc.LENGTHS:
            cell = scores["families"][family][str(length)]
            rows = [p for p in predictions if p["family"] == family and p["L"] == length]
            eligible = [p for p in rows if p["eligible"]]
            counts = cell["counts"]
            assert counts["planned"] == len(rows) == 18 == counts["unknown"] + counts["no_alarm"] + counts["late_alarm"] + counts["eligible"]
            assert counts["eligible"] == len(eligible) == cell["n"] and [p["replica"] for p in eligible] == sorted(p["replica"] for p in eligible)
            y = [p["y"] for p in eligible]
            for model in pc.MODELS:
                q = [p["predictions"][model]["q"] for p in eligible]
                assert cell["models"][model]["brier"] == pytest.approx(sum((a - b)**2 for a, b in zip(q, y)) / len(y), rel=1e-12)
                assert sum(b["count"] for b in cell["models"][model]["calibration"]) == len(y)
            base = cell["models"]["F0+"]["brier"]
            assert cell["models"]["F0+"]["skill"] == 0.0 and cell["models"]["F-B-lite"]["skill"] == pytest.approx(1 - cell["models"]["F-B-lite"]["brier"] / base)
            assert cell["contrasts"]["FAB_versus_FB"]["delta_BS"] == pytest.approx(cell["models"]["F-B-lite"]["brier"] - cell["models"]["F-AB-lite"]["brier"])
            assert cell["bootstrap"]["seed"] == namespace_seed(ctx.plan.namespace, "fault/" + family, "analysis", f"bootstrap-L{length:05d}")
            assert cell["bootstrap"]["resamples"] == 10000 == PROTOCOL["scores"]["bootstrap"]["resamples"]
            low, high = cell["models"]["F-B-lite"]["skill_interval_criterion_4"]
            narrow = cell["models"]["F-B-lite"]["skill_interval_95"]
            assert low <= narrow[0] <= narrow[1] <= high                  # the criterion level is the wider one
            oracle = cell["oracle_all_evaluation_runs"]
            assert oracle["n"] == 18 and (oracle["brier"] is None) == (oracle["unavailable"] > 0)
    assert scores["oracle_rule"] == "frozen"


def test_a_missing_forecast_makes_the_primary_score_unavailable(campaign):
    ctx, _ = campaign
    y = np.array([1, 0, 1, 1, 0], dtype=float)
    complete = {name: [0.6, 0.3, 0.7, 0.8, 0.2] for name in pc.MODELS}
    assert pf.score_cell(ctx, "LAG1", 20480, y, complete)["status"] == "SCORED"
    holed = {**complete, "F-B-lite": [0.6, None, 0.7, 0.8, 0.2]}
    cell = pf.score_cell(ctx, "LAG1", 20480, y, holed)
    assert cell["status"] == "INSUFFICIENT_INFORMATION" and cell["bootstrap"] is None
    assert cell["models"]["F-B-lite"]["missing"] == 1 and cell["models"]["F-B-lite"]["brier"] is None
    assert cell["models"]["F-B-lite"]["available_case_brier_descriptive"] is not None and cell["models"]["F0+"]["skill"] is None
    empty = pf.score_cell(ctx, "LAG1", 20480, np.zeros(0), {name: [] for name in pc.MODELS})
    assert empty["status"] == "INSUFFICIENT_INFORMATION" and empty["n"] == 0
    perfect = {name: [1.0, 0.0, 1.0, 1.0, 0.0] for name in pc.MODELS}
    undefined = pf.score_cell(ctx, "LAG1", 20480, y, perfect)
    assert undefined["models"]["F-B-lite"]["skill"] is None and undefined["models"]["F-B-lite"]["skill_interval_criterion_4"] is None
    assert undefined["bootstrap"]["undefined_skill_resamples"] == 10000    # a zero baseline is never discarded
    for q in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, np.nextafter(0.2, 0)):
        table = pf.calibration_table(ctx, np.array([q]), np.array([1.0]))
        assert [b["count"] for b in table].index(1) == min(int(np.searchsorted(pf.EDGES, q, side="right")) - 1, 4)
    oracle = pf.score_oracle(ctx, y, [{"q": 0.5, "lower_term_replaced": False}] * 4 + [{"q": None, "lower_term_replaced": False}])
    assert oracle["brier"] is None and oracle["unavailable"] == 1 and oracle["available_case_brier_descriptive"] == 0.25


# --------------------------------------------------------------------------- criteria
def test_criterion_1_agrees_with_the_independent_rule_at_every_count():
    for k in range(1780, 2001):
        lower, upper = scipy.stats.beta.ppf(0.025, k, 2000 - k + 1), (1.0 if k == 2000 else scipy.stats.beta.ppf(0.975, k + 1, 2000 - k))
        healthy = {"coverage": {str(length): {"covered": k, "runs": 2000, "cp95": [float(lower), float(upper)]} for length in pc.LENGTHS}}
        mine, theirs = pa.criterion_1(healthy), rc.criterion_1({length: k for length in pc.LENGTHS})
        assert mine["result"] == theirs["result"]
        assert mine["per_length"]["10240"]["label"] == theirs["per_length"][10240]["label"]
    mixed = {"coverage": {"10240": {"cp95": [0.93, 0.97]}, "20480": {"cp95": [0.915, 0.94]}, "40960": {"cp95": [0.93, 0.97]}}}
    assert pa.criterion_1(mixed)["result"] == "INSUFFICIENT_INFORMATION" and pa.criterion_1(mixed)["action"] is None
    mixed["coverage"]["20480"]["cp95"] = [0.89, 0.915]
    assert pa.criterion_1(mixed) == {**pa.criterion_1(mixed), "result": "FAIL", "classification": "MISCALIBRATION", "action": "stop D as designed"}
    mixed["coverage"]["20480"]["cp95"] = [0.92, 0.98]
    assert pa.criterion_1(mixed)["result"] == "PASS"       # equality with a band edge is inside


def issuance_rows(family, eligible, valid, already, alarm=12000, unknown=0, planned=1000):
    rows = []
    for i in range(planned):
        status = "UNKNOWN" if i < unknown else "ALARM" if i < unknown + eligible else "NO_ALARM"
        sweep = alarm if status == "ALARM" else None
        rows.append({"family": family, "valid": {"20480": i - unknown < valid},
                     "full": {"status": status, "alarm_sweep": sweep,
                              "eligibility": {"20480": {"eligible": None if status == "UNKNOWN" else (sweep is not None and sweep < 20480)}},
                              "issuance_endpoint": {"already_invalid": status == "ALARM" and i - unknown < already}}})
    return rows


@pytest.mark.parametrize("spec", [
    dict(SHIFT=(900, 450, 90), LAG1=(900, 450, 90), DUP=(900, 855, 90)),
    dict(SHIFT=(900, 450, 90), LAG1=(900, 855, 90), DUP=(900, 855, 90)),
    dict(SHIFT=(900, 450, 90), LAG1=(299, 150, 30), DUP=(900, 855, 90)),
    dict(SHIFT=(900, 180, 449), LAG1=(900, 720, 449), DUP=(900, 450, 450)),
    dict(SHIFT=(300, 60, 0), LAG1=(300, 240, 0), DUP=(300, 241, 0)),
])
def test_criteria_2_and_3_agree_with_the_independent_rule(spec):
    rows = sum((issuance_rows(family, *values) for family, values in spec.items()), [])
    mine = pa.criteria_2_3(rows)
    theirs = rc.criteria_2_3({family: {"alarm": [r["full"]["alarm_sweep"] for r in rows if r["family"] == family],
                                       "valid": [r["valid"]["20480"] for r in rows if r["family"] == family],
                                       "already_invalid": [r["full"]["issuance_endpoint"]["already_invalid"] for r in rows if r["family"] == family]}
                              for family in spec})
    assert (mine["criterion_2"], mine["criterion_3"]) == (theirs["criterion_2"], theirs["criterion_3"])
    for family in spec:
        for key in ("eligible", "fraction_valid", "fraction_already_invalid", "median_remaining", "criterion_2", "criterion_3"):
            assert mine["per_family"][family][key] == theirs["per_family"][family][key]


def test_criteria_2_and_3_edges():
    late = sum((issuance_rows(f, 900, 450, 90, alarm=15361) for f in pc.FAMILIES), [])
    assert pa.criteria_2_3(late)["criterion_3"] == "FAIL" and pa.criteria_2_3(late)["criterion_2"] == "PASS"
    edge = sum((issuance_rows(f, 900, 450, 90, alarm=15360) for f in pc.FAMILIES), [])
    assert pa.criteria_2_3(edge)["criterion_3"] == "PASS"                    # a median of exactly 5,120 passes
    at_length = sum((issuance_rows(f, 900, 450, 90, alarm=20480) for f in pc.FAMILIES), [])
    assert pa.criteria_2_3(at_length)["per_family"]["DUP"]["eligible"] == 0   # an alarm at L is not before L
    unknown = issuance_rows("SHIFT", 900, 450, 90, unknown=1) + issuance_rows("LAG1", 900, 450, 90) + issuance_rows("DUP", 900, 855, 90)
    out = pa.criteria_2_3(unknown)
    assert out["per_family"]["SHIFT"]["criterion_2"] is None and out["criterion_2"] == "INSUFFICIENT_INFORMATION"


def test_criterion_4_and_5():
    def scores(n, skill, lower, baseline=0.2):
        cell = {"counts": {"eligible": n}, "models": {"F0+": {"brier": baseline},
                                                      "F-B-lite": {"skill": skill, "skill_interval_criterion_4": None if lower is None else [lower, 0.5]}}}
        return {"families": {"SHIFT": {"20480": cell}}}

    assert pa.criterion_4(scores(400, 0.10, 1e-9))["result"] == "PASS"                       # one passing family is enough
    assert pa.criterion_4(scores(400, 0.05, 1e-9))["result"] == "INSUFFICIENT_INFORMATION"   # none passes, two have no scores
    assert pa.criterion_4({"families": {}})["result"] == "INSUFFICIENT_INFORMATION"
    full = lambda *cells: {"families": {f: {"20480": c["families"]["SHIFT"]["20480"]} for f, c in zip(pc.FAMILIES, cells)}}
    weak = scores(400, 0.0, -0.1)
    assert pa.criterion_4(full(scores(400, 0.10, 1e-9), weak, weak))["result"] == "PASS"     # skill of exactly 0.10 passes
    assert pa.criterion_4(full(scores(400, 0.0999, 0.05), weak, weak))["result"] == "FAIL"
    assert pa.criterion_4(full(scores(400, 0.3, 0.0), weak, weak))["result"] == "FAIL"       # the lower bound must be above 0
    assert pa.criterion_4(full(scores(239, 0.5, 0.2), weak, weak))["result"] == "INSUFFICIENT_INFORMATION"
    assert pa.criterion_4(full(scores(400, 0.5, None), weak, weak))["result"] == "INSUFFICIENT_INFORMATION"
    assert pa.criterion_4(full(scores(400, None, None, baseline=0.0), weak, weak))["result"] == "INSUFFICIENT_INFORMATION"
    assert pa.criterion_5(None, "2026-12-01", "2027-01-10")["result"] == "INSUFFICIENT_INFORMATION"
    assert pa.criterion_5("2027-04-30", None, "2027-01-10")["result"] == "INSUFFICIENT_INFORMATION"
    on_time = pa.criterion_5("2027-04-30", "2026-12-22", "2027-01-30")
    assert on_time["result"] == "PASS" and on_time["results_deadline"] == "2027-01-30"
    assert pa.criterion_5("2027-04-30", "2026-12-23", "2027-01-30")["result"] == "FAIL"      # engineering one day late
    assert pa.criterion_5("2027-04-30", "2026-12-22", "2027-01-31")["result"] == "FAIL"
    assert pa.criterion_5("2027-05-31", "2026-12-01", "2027-02-28")["results_deadline"] == "2027-02-28"   # clamped to the month's end
    assert pa.criterion_5("2028-05-31", "2026-12-01", "2028-02-29")["results_deadline"] == "2028-02-29"   # leap year
    assert pa.criterion_5("2027-02-15", "2026-12-01", "2026-11-15")["results_deadline"] == "2026-11-15"   # across the year


def test_precedence_rules_in_the_criteria_file(campaign):
    ctx, _ = campaign
    criteria = pc.read_json(ctx.result("criteria.json"))
    assert criteria["V2_D"]["result"] == "FAIL"                          # 24 healthy runs cannot bound the rate below 0.10
    assert [criteria["results"][k] for k in "234"] == ["UNINTERPRETABLE"] * 3
    assert [criteria["results_before_precedence"][k] for k in "234"] == ["INSUFFICIENT_INFORMATION"] * 3   # fewer than 300 and 240
    assert criteria["results"]["5"] == "INSUFFICIENT_INFORMATION" and criteria["criterion_5"]["tdk_deadline"] is None
    assert criteria["implementation_checks_passed"] is True and criteria["interpretation_blocked"] is False
    assert all(criteria["counts_as_failed_for_the_stop_rule"][k] for k in "2345")
    checks = pc.read_json(ctx.result("implementation_checks.json"))
    assert checks["checks"]["I-3-step"]["expected_runs"] == ["fault/SHIFT/00000", "fault/LAG1/00000", "fault/DUP/00000"]
    assert checks["checks"]["I-1"]["runs_checked"] == 207 and checks["checks"]["I-3"]["runs_checked"] == 144
    healthy = pc.read_json(ctx.result("healthy.json"))
    assert len(healthy["short_interval_fixed_checkpoints"]["per_checkpoint"]) == 159
    assert healthy["short_interval_fixed_checkpoints"]["runs_with_short_intervals"] == 24
    faults = pc.read_json(ctx.result("faults.json"))
    for family in pc.FAMILIES:
        block = faults["families"][family]
        assert block["all"]["runs"] == 48 == block["fitting"]["runs"] + block["evaluation"]["runs"]
        assert sum(block["all"]["validity_pattern"].values()) == 48 and sum(block["all"]["alarm_class"].values()) == 48


def test_a_smoke_is_never_analysed(tmp_path):
    ctx = pc.Context(pc.SMOKE_PLAN, tmp_path, None, au.dlib(), None)
    for stage in (pa.calibrate, pa.issue, pa.healthy_summary, pa.fault_summary, pf.fit, pf.evaluate, pa.criteria):
        with pytest.raises(pc.PilotError, match="engineering and timing only"):
            stage(ctx)
    assert not list(tmp_path.iterdir())


def test_a_damaged_store_stops_the_analysis(campaign):
    ctx, _ = campaign
    trace = ctx.trace_path("fault/SHIFT", 5)
    original = trace.read_bytes()
    try:
        trace.write_bytes(original[:-200])                               # truncated
        with pytest.raises(pc.PilotError, match="no valid record"):
            pg.require_complete(ctx)
        assert pg.classify(ctx, "fault/SHIFT", 5)["reason"] == "trace hash differs from the record"
        with pytest.raises(ValueError):
            pa.load_run(ctx, "fault/SHIFT", 5)                           # the loader checks the hashes itself
    finally:
        trace.write_bytes(original)
    record_path = ctx.record_path("healthy", 2)
    original = record_path.read_bytes()
    try:
        record = json.loads(original)
        record["endpoints"]["20480"]["valid"] = not record["endpoints"]["20480"]["valid"]   # a plausible, wrong record
        record_path.write_bytes(pc.canonical_json(record))
        with pytest.raises(pc.PilotError, match="differ from the run manifest"):
            pg.require_complete(ctx)
    finally:
        record_path.write_bytes(original)
    assert pg.require_complete(ctx) == {"complete": 207, "failed_records": 0}


def test_criterion_5_of_both_implementations_agree_on_every_day_of_two_years():
    import datetime
    day, checked = datetime.date(2027, 1, 1), 0
    while day <= datetime.date(2028, 12, 31):
        deadline = day.isoformat()
        boundary = datetime.date.fromisoformat(pa.criterion_5(deadline, "2026-12-01", "2026-12-01")["results_deadline"])
        for engineering in ("2026-12-22", "2026-12-23"):
            for shift in (-1, 0, 1):
                results = (boundary + datetime.timedelta(days=shift)).isoformat()
                primary, independent = pa.criterion_5(deadline, engineering, results), rc.criterion_5(deadline, engineering, results)
                assert primary["result"] == independent["result"] and primary["results_deadline"] == independent["results_deadline"]
                assert primary["engineering_deadline"] == independent["engineering_deadline"]
                checked += 1
        day += datetime.timedelta(days=1)
    assert checked == 731 * 6


def criterion_5_tree(ctx, tmp_path, engineering, results, *, deadline="2027-04-30", written="2026-12-05T10:00:00+0100"):
    """A copy of the result files with a manifest and the criterion 5 block the pipeline would write."""
    shutil.copytree(ctx.root / "results", tmp_path / "results", ignore=shutil.ignore_patterns("runs"))
    (tmp_path / "pre_generation_manifest.json").write_bytes(pc.canonical_json({"tdk_deadline": deadline, "written_at": written}))
    criteria = pc.read_json(tmp_path / "results/criteria.json")
    criteria["criterion_5"] = pa.criterion_5(deadline, engineering, results)
    criteria["results_before_precedence"]["5"] = criteria["results"]["5"] = criteria["criterion_5"]["result"]
    (tmp_path / "results/criteria.json").write_bytes(pc.canonical_json(criteria))
    return criteria


def test_criterion_5_check_on_the_campaign_and_with_dates(campaign, tmp_path):
    ctx, _ = campaign
    bare = vf.check_criterion_5(ctx.root)                               # a test plan: no manifest, no dates
    assert bare["agrees"] and bare["compared"] == 7 and bare["independent"]["result"] == "INSUFFICIENT_INFORMATION"
    assert bare["manifest_present"] is False
    criterion_5_tree(ctx, tmp_path, "2026-12-01", "2026-12-20")
    dated = vf.check_criterion_5(tmp_path)
    assert dated["agrees"] and dated["independent"]["result"] == "PASS" and dated["independent"]["manifest_date"] == "2026-12-05"
    absent = vf.check_criterion_5(tmp_path / "nowhere")
    assert absent["status"] == "NOT_CHECKED" and not absent["agrees"]


@pytest.mark.parametrize("engineering, results, change, expected", [
    ("2026-12-01", "2026-12-20", lambda c: c["results"].update({"5": "FAIL"}), "final 5"),
    ("2026-12-23", "2026-12-28", lambda c: c["results_before_precedence"].update({"5": "PASS"}), "5"),
    ("2026-12-01", "2026-12-20", lambda c: c["criterion_5"].update(tdk_deadline="2027-06-30"), "5 declared deadline"),
    ("2026-12-01", "2026-12-20", lambda c: c["criterion_5"].update(results_deadline="2027-02-28"), "5 results deadline"),
    ("2026-12-01", "2026-12-20", lambda c: c["criterion_5"].update(engineering_deadline="2027-01-31"), "5 engineering deadline"),
    ("2026-12-01", "2026-12-20", lambda c: c["criterion_5"].update(result="FAIL"), "5 result"),
    ("2026-12-01", "2026-12-20", lambda c: c["criterion_5"].update(results_available_on="2026-12-04"), "5 dates"),
    ("2026-12-01", "2026-12-20", lambda c: c["criterion_5"].update(engineering_complete_on="2026-12-06"), "5 dates"),
    ("2026-12-01", "2026-12-20", lambda c: c["criterion_5"].update(results_available_on="20 December"), "5 dates"),
    ("2026-12-01", "2026-12-20", lambda c: c["criterion_5"].update(results_available_on="2027-02-15"), "5"),
])
def test_criterion_5_check_reports_a_changed_or_inconsistent_date(campaign, tmp_path, engineering, results, change, expected):
    ctx, _ = campaign
    criteria = criterion_5_tree(ctx, tmp_path, engineering, results)
    change(criteria)
    (tmp_path / "results/criteria.json").write_bytes(pc.canonical_json(criteria))
    out = vf.check_criterion_5(tmp_path)
    assert expected in {d["field"] for d in out["discrepancies"]} and not out["agrees"]


# --------------------------------------------------------------------------- interrupted runs: the audit log and unresolved runs
def test_the_attempt_log_of_the_campaign_is_checked_independently(campaign, recounted):
    ctx, _ = campaign
    check = vf.check_attempt_log(ctx.root, ctx.plan.namespace, fx.audit_plan(ctx.plan))
    assert check["agrees"] and check["status"] == "CHECKED" and check["journal_entries"] == 2 * 207 + 2
    assert check["compared"] == 16 * 207 + 2 and check["interruptions"] == 0 and check["unresolved"] == []
    assert check["run_manifest_present"] and not check["pre_generation_manifest_present"] and len(check["environments"]) == 1
    criteria = pc.read_json(ctx.result("criteria.json"))
    assert criteria["unresolved_infrastructure"] == {
        "status": None, "runs": [], "per_condition": {c: 0 for c, _, _ in ctx.plan.conditions},
        "criteria_that_cannot_pass": [], "max_infrastructure_retries": 3}
    assert recounted[0]["unresolved_infrastructure"] == {
        "runs": [], "per_condition": {c: 0 for c, _, _ in ctx.plan.conditions}, "criteria_that_cannot_pass": []}


def test_criteria_that_use_a_condition_are_the_same_on_both_sides():
    assert rc.criteria_blocked_by([]) == [] and rc.criteria_blocked_by(["calibration"]) == [2, 3, 4]
    assert rc.criteria_blocked_by(["healthy"]) == [1, 2, 3, 4] and rc.criteria_blocked_by(["fault/DUP", "fault/SHIFT"]) == [2, 3, 4]
    assert rc.criteria_blocked_by(["fault/LAG1", "healthy"]) == [1, 2, 3, 4]
    assert {k: tuple(map(int, v)) for k, v in pa.CRITERIA_USING.items()} == rc.CRITERIA_USING
    assert all("5" not in keys for keys in pa.CRITERIA_USING.values())         # the dates do not depend on any run


def test_an_unresolved_run_keeps_the_criteria_that_use_it_from_passing(tmp_path, monkeypatch):
    plan = pc.Plan(fx.NAMESPACE, (("calibration", "healthy", 39), ("healthy", "healthy", 3), ("fault/SHIFT", "SHIFT", 4),
                                  ("fault/LAG1", "LAG1", 4), ("fault/DUP", "DUP", 4)), range(0, 2), range(2, 4), 39, 38, "test")
    root = tmp_path / "store"
    ctx = fx.context(plan, root)
    fx.set_interruptions(root, {"healthy/1": ["OSError"] * 4})
    for _ in range(4):
        with pytest.raises(pg.InfrastructureInterruption):
            pg.generate(ctx, engine_spec=fx.INTERRUPTED_SPEC, workers=0)
    out = pg.generate(ctx, engine_spec=fx.INTERRUPTED_SPEC, workers=0)
    assert out["unresolved_infrastructure"] == ["healthy/00001"]
    assert out["counts"]["COMPLETE"] == 53 and out["counts"]["FAILED_RECORD"] == 1
    pa.calibrate(ctx)
    pf.fit(ctx)
    pa.issue(ctx)
    pf.evaluate(ctx)
    healthy = pa.healthy_summary(ctx)
    assert healthy["V2_D"]["unknown"] == 1                              # counted as a failed run, not left out
    pa.fault_summary(ctx)
    pa.implementation_checks(ctx, None)
    primary_rule, independent_rule = pa.criterion_1, rc.criterion_1
    # The plan is too small for any criterion to pass by itself: let criterion 1 pass, on both sides.
    monkeypatch.setattr(pa, "criterion_1", lambda summary: {**primary_rule(summary), "result": "PASS"})
    monkeypatch.setattr(rc, "criterion_1", lambda covered, runs: {**independent_rule(covered, runs), "result": "PASS"})
    criteria = pa.criteria(ctx)
    assert criteria["results_before_precedence"]["1"] == "PASS" and criteria["results"]["1"] == "INSUFFICIENT_INFORMATION"
    assert criteria["counts_as_failed_for_the_stop_rule"]["1"] is True
    assert "PASS" not in [criteria["results"][key] for key in "1234"]
    assert criteria["unresolved_infrastructure"] == {
        "status": "UNRESOLVED_INFRASTRUCTURE", "runs": [{"condition": "healthy", "replica": 1}],
        "per_condition": {"calibration": 0, "healthy": 1, "fault/SHIFT": 0, "fault/LAG1": 0, "fault/DUP": 0},
        "criteria_that_cannot_pass": ["1", "2", "3", "4"], "max_infrastructure_retries": 3}
    independent, comparison = recount(ctx)
    assert independent["criterion_1"]["result"] == "PASS" and independent["verdicts"]["criterion_1"] == "INSUFFICIENT_INFORMATION"
    assert independent["unresolved_infrastructure"] == {
        "runs": [{"condition": "healthy", "replica": 1}],
        "per_condition": {"calibration": 0, "healthy": 1, "fault/SHIFT": 0, "fault/LAG1": 0, "fault/DUP": 0},
        "criteria_that_cannot_pass": [1, 2, 3, 4]}
    assert [d for d in comparison["discrepancies"] if d["field"] not in ORACLE_FIELDS] == []
    check = vf.check_attempt_log(root, fx.NAMESPACE, fx.audit_plan(plan))
    assert check["agrees"], check["discrepancies"]
    assert check["unresolved"] == [{"condition": "healthy", "replica": 1}] and check["criteria_that_cannot_pass"] == [1, 2, 3, 4]
    # A criteria file that lets the criterion pass, or does not name the run, is reported by both checks.
    copy_root = tmp_path / "copy"
    shutil.copytree(root, copy_root)
    changed = copy.deepcopy(criteria)
    changed["results"]["1"] = "PASS"
    changed["unresolved_infrastructure"]["runs"] = []
    changed["unresolved_infrastructure"]["criteria_that_cannot_pass"] = []
    (copy_root / "results/criteria.json").write_bytes(pc.canonical_json(changed))
    found = {(d["id"], d["field"]) for d in vf.check_attempt_log(copy_root, fx.NAMESPACE, fx.audit_plan(plan))["discrepancies"]}
    assert found == {("criteria", "unresolved runs"), ("criteria", "criteria that cannot pass"), ("criteria", "final 1")}
    assert "final 1" in {d["field"] for d in vf.compare_with_primary(copy_root, independent)["discrepancies"]}
    # Without the unresolved run nothing is blocked: the same patched rule passes.
    views = {c: [vf.load_view(root, c, i) for i in range(n)] for c, (_, n) in fx.audit_plan(plan).items()}
    for view in views["healthy"]:
        view.failure_phase = None
    free = vf.recount_campaign(views, fx.NAMESPACE, plan=fx.audit_plan(plan), fitting=plan.fitting, evaluation=plan.evaluation,
                               saved_models=vf.saved_models_from_primary(pc.read_json(ctx.result("models.json"))))
    assert free["verdicts"]["criterion_1"] == "PASS" and free["unresolved_infrastructure"]["criteria_that_cannot_pass"] == []
