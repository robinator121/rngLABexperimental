"""Synthetic fixtures for the Stage D pipeline tests. No sampler runs and no pilot or smoke draw.

`SyntheticEngine` returns records in the measurement layer's own schema: endpoints, application
checkpoints and short intervals are computed by the measurement layer's functions on synthetic
traces, and battery checkpoints carry synthetic p-values. The traces are random numbers from
generators seeded by run labels of the audit-only namespace; they are not Monte Carlo output.

The engine is loaded by the driver through the specification "<this file>:synthetic_engine".
"""
from __future__ import annotations

import importlib.util
import math
import sys
import zlib
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("_audit_loader", _HERE / "_audit_loader.py")
au = sys.modules.get("_audit_loader") or importlib.util.module_from_spec(_spec)
if "_audit_loader" not in sys.modules:
    sys.modules["_audit_loader"] = au
    _spec.loader.exec_module(au)

PILOT = au.REPO / "research/v019-d-pilot"
pc = au.load("pilot_common", PILOT / "pilot_common.py")
pg = au.load("pilot_generate", PILOT / "pilot_generate.py")
pa = au.load("pilot_analysis", PILOT / "pilot_analysis.py")
pf = au.load("pilot_forecast", PILOT / "pilot_forecast.py")
NAMESPACE = au.AUDIT_NAMESPACE
# 39 calibration runs with the 38th order statistic give the same nominal rate, 2/40 = 0.05.
TEST_PLAN = pc.Plan(NAMESPACE, (("calibration", "healthy", 39), ("healthy", "healthy", 24), ("fault/SHIFT", "SHIFT", 48),
                                ("fault/LAG1", "LAG1", 48), ("fault/DUP", "DUP", 48)), range(0, 30), range(30, 48), 39, 38, "test")
SMALL_PLAN = pc.Plan(NAMESPACE, (("calibration", "healthy", 2), ("healthy", "healthy", 2), ("fault/SHIFT", "SHIFT", 2),
                                 ("fault/LAG1", "LAG1", 2), ("fault/DUP", "DUP", 2)), range(0, 1), range(1, 2), 2, 2, "test")
ENGINE_SPEC = f"{Path(__file__).resolve()}:synthetic_engine"
FAILING_SPEC = f"{Path(__file__).resolve()}:failing_engine"
PARTLY_FAILING_SPEC = f"{Path(__file__).resolve()}:partly_failing_engine"
NON_FINITE_SPEC = f"{Path(__file__).resolve()}:non_finite_engine"


def context(plan, root, *, need_measure: bool = True):
    """A pipeline context. With need_measure the test is skipped when the measurement layer is absent."""
    measure = au.measure() if need_measure else (au.load("v019_audit_measure", au.measure_path()) if au.measure_path() else None)
    c_analyze = au.load("d_pilot_c_analyze_test", au.REPO / "research/v018-h1-detection/analyze.py")
    return pc.Context(plan, Path(root), measure, measure.dlib if measure is not None else au.dlib(), c_analyze)


class SyntheticEngine:
    """Deterministic synthetic runs: a ramp profile, a biased trace, a battery statistic that follows the ramp."""

    def __init__(self, ctx):
        self.context = ctx

    def execute(self, condition, family, replica):
        ctx = self.context
        measure, dlib, plan = ctx.measure, ctx.dlib, ctx.plan
        rng = np.random.default_rng([zlib.crc32(condition.encode()), replica, 2026])
        horizon = pc.H
        if family == "healthy":
            profile = dict.fromkeys(dlib.PROFILE_FIELDS)
            severity = np.zeros(horizon, dtype=np.float64)
            cap = 1.0
        else:
            cap = float(dlib.CAPS[family])
            t0 = int(rng.integers(2048, 8193))
            t_cap = float(np.exp(rng.uniform(np.log(2048.0), np.log(163840.0))))
            profile = {"t0": t0, "T_cap": t_cap, "r": cap / t_cap, "s_max": cap}
            sweep = np.arange(horizon, dtype=np.float64)
            severity = np.where(sweep >= t0, np.minimum(cap, cap / t_cap * (sweep - t0)), 0.0)
        cv = measure.E_M + np.round(rng.normal(0.0, 0.4, horizon), 2) + 0.035 / cap * severity
        td = cv + np.round(rng.normal(0.0, 0.15, horizon), 2)
        marks = np.cumsum(128 + 8 * rng.integers(8, 13, horizon)).astype(np.int64)
        battery = []
        for end in measure.BATTERY_ENDS:
            z = float(rng.normal() + 60.0 * severity[end - 512:end].mean() / cap)
            p = rng.uniform(0.001, 1.0, 10).tolist()
            p[0] = math.erfc(abs(z) / math.sqrt(2))
            s, count = measure.statistic(p)
            s_prime, count_prime = measure.statistic(p[:8])
            battery.append({"end": int(end), "p_values": p, "S": s, "S_prime": s_prime, "available_p_count": count,
                            "available_p_count_prime": count_prime, "complete": True, "complete_prime": True, "failures": [],
                            "frequency_signed_z": z, "lag1_feature": z, "duplicate_rate": max(z, 0.0) / 1000.0, "was_cut": False})
        record = {"namespace": plan.namespace, "condition": condition, "family": family, "replica": replica,
                  "split": plan.split(condition, replica), "seeds": dlib.seeds_for(plan.namespace, condition, replica),
                  "sweeps": horizon, "profile": profile, "profile_source": "synthetic fixture",
                  "B0": None if family == "healthy" else int(marks[profile["t0"] - 1]), "local_bytes": int(marks[-1]),
                  "conventional": battery,
                  "application": [measure.application_checkpoint(cv, td, marks, end) for end in measure.APPLICATION_ENDS],
                  "endpoints": {str(length): measure.endpoint(cv, length, measure.E_M) for length in measure.LENGTHS},
                  "failures": [], "failure": None, "status": "SUCCESS", "completed_sweeps": horizon,
                  "I-1": {"check": "I-1", "passed": True, "mismatches": []}, "I-2": {"check": "I-2", "passed": True, "mismatches": []}}
        if family != "healthy":
            record["I-3"] = {"check": "I-3", "passed": True, "mismatches": []}
            if replica == 0:
                record["I-3-step"] = {"check": "I-3-step", "passed": True, "mismatches": []}
        if condition == "healthy":
            record["short_intervals"] = measure.short_interval_check(cv)
        record["generation_seconds"], record["check_seconds"], record["measurement_seconds"] = 0.0, 0.0, 0.0
        arrays = {"CV": cv.astype("<f8"), "TD": td.astype("<f8"), "marks": marks.astype("<i8"), "severity": severity.astype("<f8")}
        return record, arrays


class FailingEngine:
    """Every run raises in 'generation': exercises failure records without the measurement layer."""

    def __init__(self, ctx):
        self.context = ctx

    def execute(self, condition, family, replica):
        plan = self.context.plan
        seeds = self.context.dlib.seeds_for(plan.namespace, condition, replica)
        return pg.failure_record(plan, condition, family, replica, "generation", ValueError("synthetic failure"), seeds,
                                 {"completed_sweeps": 7, "last_completed_sweep": 6, "local_bytes": 1024}), None


class PartlyFailingEngine(SyntheticEngine):
    """Synthetic runs, except that calibration replica 1 fails in generation."""

    def execute(self, condition, family, replica):
        if condition == "calibration" and replica == 1:
            return FailingEngine(self.context).execute(condition, family, replica)
        return super().execute(condition, family, replica)


class NonFiniteEngine(SyntheticEngine):
    """Synthetic runs whose trace of healthy replica 0 holds a NaN: the storage layer must refuse it."""

    def execute(self, condition, family, replica):
        record, arrays = super().execute(condition, family, replica)
        if condition == "healthy" and replica == 0:
            arrays["CV"][17] = np.nan
        return record, arrays


def synthetic_engine(ctx):
    return SyntheticEngine(ctx)


def partly_failing_engine(ctx):
    return PartlyFailingEngine(ctx)


def non_finite_engine(ctx):
    return NonFiniteEngine(ctx)


def failing_engine(ctx):
    return FailingEngine(ctx)
