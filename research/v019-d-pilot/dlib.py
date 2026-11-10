"""Stage D C2 only: PROTOCOL.md sections 2.1–2.2 and 7; protocol.json governs.

No campaigns, endpoints or durable output. Raw streams live only in memory.
Stage C is imported without mutation; replay deliberately duplicates its math.
"""
from __future__ import annotations

import hashlib
import importlib.util
import math
import struct
from pathlib import Path

import numpy as np

from rnglab.seeds import namespace_seed
from rnglab14.models.pimc_harmonic import (
    PathIntegralHarmonicOscillator, sample_exact_finite_bead_paths,
)
from rnglab14.sources import PCG64Source
from rnglab14.sources.base import RandomSource, SourceManifest

REPO = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location(
    "d_pilot_stage_c", REPO / "research/v018-h1-detection/h1lib.py")
h1lib = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(h1lib)
TOP, STD = h1lib.TOP, h1lib.STD
H = 40960
# Stage C results/gates.json $.cells['SHIFT-5'/'LAG1-5'/'DUP-5'].severity.
CAPS = {"SHIFT": 0.046415888336127774, "LAG1": 0.5, "DUP": 0.5}
PROFILE_FIELDS = ("t0", "T_cap", "r", "s_max")
Profile = dict[str, int | np.float64 | None]


def seeds_for(ns: str, condition: str, i: int) -> dict[str, int]:
    """protocol.json data.seeds.exact_reuse: discard gfsr, add profile."""
    base = h1lib.seeds_for(ns, condition, i)
    seeds = {role: base[role] for role in ("local", "centroid", "start", "fault")}
    seeds["profile"] = namespace_seed(ns, condition, "profile", f"{i:05d}")
    return seeds


def draw_profile(family: str, seeds: dict[str, int]) -> Profile:
    """protocol.json profile.draw_order; healthy consumes no profile variates."""
    if family == "healthy":
        return dict.fromkeys(PROFILE_FIELDS)
    s_max = np.float64(CAPS[family])
    profile_rng = np.random.default_rng(seeds["profile"])
    t0 = int(profile_rng.integers(2048, 8193))
    T_cap = np.float64(np.exp(profile_rng.uniform(np.log(2048.0), np.log(163840.0))))
    r = np.float64(s_max / T_cap)
    return {"t0": t0, "T_cap": T_cap, "r": r, "s_max": s_max}


def severity(t: int, profile: Profile) -> np.float64:
    """protocol.json profile.severity, directly from t (no accumulated slope)."""
    t0 = profile["t0"]
    if t0 is None or t < t0:
        return np.float64(0.0)
    return np.minimum(np.float64(profile["s_max"]),
                      np.float64(profile["r"] * np.float64(t - t0)))


def severity_path(profile: Profile, sweeps: int = H) -> np.ndarray:
    """protocol.json profile.record: one float64 severity per integer sweep."""
    return np.fromiter((severity(t, profile) for t in range(sweeps)),
                       dtype=np.float64, count=sweeps)


class FaultyWords(RandomSource):
    """Sweep-driven local adapter; protocol.json fault_transforms and raw_bytes.

    step=True is solely the separate I-3 step test. B0/prefix remain None
    until the first onset word is read (including for truncated test runs).
    """

    def __init__(self, base: RandomSource, family: str, profile: Profile,
                 fault_seed: int, *, step: bool = False):
        if family not in ("healthy", *CAPS):
            raise ValueError(f"unknown family {family!r}")
        super().__init__(SourceManifest(
            source_id=f"d_{family}", family="synthetic_fault", is_deterministic=True,
            is_physical=False, metadata={"base": base.manifest.source_id,
                                         "family": family, "profile": dict(profile), "step": step}))
        self.base, self.family, self.profile = base, family, dict(profile)
        self.step = step
        self.current_sweep = 0
        self.emitted, self.raw_base = bytearray(), bytearray()
        self.B0: int | None = None
        self.fault_rng = np.random.default_rng(fault_seed)
        self.prev: bytes | None = None
        self.z: float | None = None

    @property
    def prefix_sha256(self) -> str | None:
        return None if self.B0 is None else hashlib.sha256(self.emitted[:self.B0]).hexdigest()

    def read_bytes(self, n: int) -> bytes:
        # RandomSource skips _read_bytes for n=0; enforce the contract here too.
        if isinstance(n, bool) or not isinstance(n, int) or n != 8:
            raise ValueError("local reads must be single 8-byte words")
        return super().read_bytes(n)

    def _transform(self, raw: bytes, s: np.float64) -> bytes:
        if self.family == "DUP":
            duplicate = self.fault_rng.random() < s
            return self.prev if duplicate and self.prev is not None else raw
        word = struct.unpack(">Q", raw)[0]
        u = (word >> 11) * 2.0**-53
        if self.family == "SHIFT":
            u = min(u + s, TOP)
        else:
            e = STD.inv_cdf(u + 2.0**-54)
            self.z = e if self.z is None else s * self.z + math.sqrt(1 - s * s) * e
            u = min(STD.cdf(self.z), TOP)
        return struct.pack(">Q", (int(u * 2.0**53) << 11) | (word & 0x7FF))

    def _read_bytes(self, n: int) -> bytes:
        if n != 8:
            raise ValueError("local reads must be single 8-byte words")
        raw = self.base.read_bytes(8)
        self.raw_base += raw
        out = raw
        if self.family != "healthy" and self.current_sweep >= self.profile["t0"]:
            if self.B0 is None:
                self.B0 = self.consumed_bytes
            s = np.float64(self.profile["s_max"]) if self.step else severity(self.current_sweep, self.profile)
            out = self._transform(raw, s)
        self.prev = out
        self.emitted += out
        return out


class SweepClock(RandomSource):
    """protocol.json sweep_clock.callback: mark, advance, then centroid read."""

    def __init__(self, base: RandomSource, local: FaultyWords):
        super().__init__(SourceManifest(source_id=f"d_clock_{base.manifest.source_id}",
                                        family="reference_prng", is_deterministic=True, is_physical=False))
        self.base, self.local = base, local
        self.marks: list[int] = []

    def read_bytes(self, n: int) -> bytes:
        if isinstance(n, bool) or not isinstance(n, int) or n != 16:
            raise ValueError("centroid reads must be exactly 16 bytes")
        return super().read_bytes(n)

    def _read_bytes(self, n: int) -> bytes:
        if n != 16:
            raise ValueError("centroid reads must be exactly 16 bytes")
        self.marks.append(self.local.consumed_bytes)
        self.local.current_sweep += 1
        return self.base.read_bytes(n)


def reapply(family: str, profile: Profile, seeds: dict[str, int], marks: np.ndarray,
            local_bytes: int, *, step: bool = False, base_stream: bytes | None = None) -> bytes:
    """Independent I-3: fresh PCG64 by default, or the recorded raw buffer.

    protocol.json sweep_clock.replay_marks and fault_transforms. No production
    transform/severity calls. The stored T_cap and cap independently yield r.
    """
    if family not in ("healthy", *CAPS):
        raise ValueError(f"unknown family {family!r}")
    if local_bytes < 0 or local_bytes % 8:
        raise ValueError("byte count must be a nonnegative multiple of eight")
    marks = np.asarray(marks)
    if (marks.ndim != 1 or marks.dtype != np.dtype("int64")
            or np.any(np.diff(np.concatenate(([0], marks))) <= 0)
            or np.any(marks % 8) or (local_bytes and (not len(marks) or marks[-1] != local_bytes))):
        raise ValueError("replay requires complete, increasing int64 word-aligned marks")
    if base_stream is not None and len(base_stream) != local_bytes:
        raise ValueError("recorded base length differs from local_bytes")
    base = PCG64Source(seeds["local"])
    rng = np.random.default_rng(seeds["fault"])
    out, previous, state = bytearray(), None, None
    for b in range(0, local_bytes, 8):
        raw = base.read_bytes(8) if base_stream is None else base_stream[b:b + 8]
        word = raw
        t = int(np.searchsorted(marks, b, side="right"))
        if family != "healthy" and t >= profile["t0"]:
            cap = np.float64(profile["s_max"])
            slope = np.float64(cap / np.float64(profile["T_cap"]))
            s = cap if step else np.minimum(cap, np.float64(slope * np.float64(t - profile["t0"])))
            if family == "DUP":
                q = rng.random()
                word = previous if q < s and previous is not None else raw
            else:
                w = struct.unpack(">Q", raw)[0]
                u = (w >> 11) * 2.0**-53
                if family == "SHIFT":
                    mapped = min(u + s, TOP)
                else:
                    innovation = STD.inv_cdf(u + 2.0**-54)
                    state = innovation if state is None else s * state + math.sqrt(1 - s * s) * innovation
                    mapped = min(STD.cdf(state), TOP)
                word = struct.pack(">Q", (int(mapped * 2.0**53) << 11) | (w & 0x7FF))
        out += word
        previous = word
    return bytes(out)


def _result(check: str, mismatches: list[dict], **details) -> dict:
    return {"check": check, "passed": not mismatches, "mismatches": mismatches, **details}


def _compare(errors: list[dict], field: str, actual, expected) -> None:
    if actual != expected:
        errors.append({"field": field, "actual": actual, "expected": expected})


def _bytes_compare(errors: list[dict], field: str, actual: bytes, expected: bytes) -> None:
    if actual != expected:
        first = next((i for i, (a, b) in enumerate(zip(actual, expected)) if a != b),
                     min(len(actual), len(expected)))
        errors.append({"field": field, "first_byte": first, "actual_length": len(actual),
                       "expected_length": len(expected), "actual_hex": actual[first:first + 8].hex(),
                       "expected_hex": expected[first:first + 8].hex()})


def check_i1(marks: np.ndarray, local_bytes: int, H: int, failed: bool,
             *, current_sweep: int | None = None) -> dict:
    """protocol.json implementation_checks.I-1; optional counter is not inferable from marks.

    Failed runs may have unmarked bytes from the incomplete sweep. The caller
    supplies the observed counter to check it too; generate_run always does so.
    """
    errors: list[dict] = []
    marks = np.asarray(marks)
    _compare(errors, "marks.dtype", str(marks.dtype), "int64")
    if marks.ndim != 1:
        return _result("I-1", errors + [{"field": "marks.shape", "actual": marks.shape, "expected": "1-D"}])
    if not failed:
        _compare(errors, "marks.count", len(marks), H)
        _compare(errors, "local_bytes", local_bytes, int(marks[-1]) if len(marks) else 0)
    elif len(marks) > H or (len(marks) and marks[-1] > local_bytes):
        errors.append({"field": "completed_sweeps", "count": len(marks), "local_bytes": local_bytes})
    increments = np.diff(np.concatenate(([0], marks)))
    bad = np.flatnonzero((increments < 128) | (increments > 256) | ((increments - 128) % 8 != 0))
    if len(bad):
        errors.append({"field": "sweep_bytes", "indices": bad.tolist(), "actual": increments[bad].tolist()})
    if current_sweep is not None:
        _compare(errors, "current_sweep", current_sweep, len(marks) if failed else H)
    return _result("I-1", errors, counter_checked=current_sweep is not None)


def check_i2(run: dict, *, expected_profile: Profile | None = None) -> dict:
    """Independent I-2 draw/path/onset/prefix checks (protocol.json implementation_checks).

    expected_profile is only for hand-built engineering fixtures: explicitly
    reports supplied-profile validation instead of claiming a seeded redraw.
    Unobserved onset in a shortened run has null B0 and prefix hash.
    """
    errors: list[dict] = []
    family, p = run["family"], run["profile"]
    if family == "healthy":
        expected = dict.fromkeys(PROFILE_FIELDS)
    elif expected_profile is not None:
        expected = expected_profile
    else:
        rng = np.random.default_rng(run["seeds"]["profile"])
        t0 = int(rng.integers(2048, 8193))
        duration = np.float64(np.exp(rng.uniform(np.log(2048.0), np.log(163840.0))))
        cap = np.float64(CAPS[family])
        expected = {"t0": t0, "T_cap": duration, "r": np.float64(cap / duration), "s_max": cap}
    for field in PROFILE_FIELDS:
        if field == "t0" or p[field] is None or expected[field] is None:
            _compare(errors, f"profile.{field}", p[field], expected[field])
        else:
            _bytes_compare(errors, f"profile.{field}", struct.pack(">d", p[field]), struct.pack(">d", expected[field]))
    path = np.asarray(run["severity"])
    _compare(errors, "severity.dtype", str(path.dtype), "float64")
    t = np.arange(run["sweeps"], dtype=np.int64)
    expected_path = np.zeros(len(t), dtype=np.float64)
    if family != "healthy":
        active = t >= p["t0"]
        expected_path[active] = np.minimum(np.float64(p["s_max"]),
                                         np.float64(p["r"]) * (t[active] - p["t0"]).astype(np.float64))
    _compare(errors, "severity.shape", path.shape, expected_path.shape)
    _bytes_compare(errors, "severity", path.astype("<f8").tobytes(), expected_path.astype("<f8").tobytes())
    marks = run["marks"]
    observed = family != "healthy" and p["t0"] < len(marks)
    onset = (0 if p["t0"] == 0 else int(marks[p["t0"] - 1])) if observed else None
    _compare(errors, "B0", run["B0"], onset)
    prefix = None
    if onset is not None:
        prefix = hashlib.sha256(PCG64Source(run["seeds"]["local"]).read_bytes(onset)).hexdigest()
        _compare(errors, "emitted_prefix_sha256", hashlib.sha256(run["emitted"][:onset]).hexdigest(), prefix)
    _compare(errors, "prefix_sha256", run["prefix_sha256"], prefix)
    return _result("I-2", errors, profile_check="supplied" if expected_profile is not None else "seeded",
                   onset_observed=observed)


def check_i3(run: dict, *, step: bool = False) -> dict:
    """protocol.json implementation_checks.I-3: raw regeneration and independent map."""
    errors: list[dict] = []
    base = PCG64Source(run["seeds"]["local"])
    raw = b"".join(base.read_bytes(8) for _ in range(run["local_bytes"] // 8))
    _bytes_compare(errors, "raw_base", run["raw_base"], raw)
    try:
        expected = reapply(run["family"], run["profile"], run["seeds"], run["marks"],
                           run["local_bytes"], step=step, base_stream=run["raw_base"])
    except (ValueError, ArithmeticError) as error:
        return _result("I-3", errors + [{"field": "reapply", "exception": repr(error)}])
    _bytes_compare(errors, "emitted", run["emitted"], expected)
    _compare(errors, "base_sha256", run["base_sha256"], hashlib.sha256(raw).hexdigest())
    _compare(errors, "emitted_sha256", run["emitted_sha256"], hashlib.sha256(expected).hexdigest())
    if run["B0"] is not None:
        _compare(errors, "prefix_sha256", run["prefix_sha256"], hashlib.sha256(expected[:run["B0"]]).hexdigest())
    return _result("I-3", errors)


def generate_run(ns: str, condition: str, family: str, i: int, sweeps: int = H,
                 *, profile: Profile | None = None) -> dict:
    """Minimal protocol.json application run; profile override is for engineering tests.

    Returns traces, metadata and in-memory buffers for I-2/I-3. Exceptions
    propagate to the caller for later campaign failure recording; no retries.
    The frozen sampler requires at least 32 sweeps even in shortened tests.
    """
    seeds = seeds_for(ns, condition, i)
    p = draw_profile(family, seeds) if profile is None else dict(profile)
    local = FaultyWords(PCG64Source(seeds["local"]), family, p, seeds["fault"])
    clock = SweepClock(PCG64Source(seeds["centroid"]), local)
    start = sample_exact_finite_bead_paths(PCG64Source(seeds["start"]), 16, 2.0, 2)[0]
    model = PathIntegralHarmonicOscillator(beads=16, temperature=0.5, displacement_delta=0.4)
    result = model.run_simulation(
        local, sweeps=sweeps, burn_in=0, initial_path=start, centroid_source=clock,
        initial_path_provenance={"namespace": ns, "condition": condition, "role": "start",
                                 "replica": i, "seed": seeds["start"]},
        uncertainty_batch_sizes=(512,), record_trace=True)
    marks = np.asarray(clock.marks, dtype=np.int64)
    emitted, raw = bytes(local.emitted), bytes(local.raw_base)
    return {"namespace": ns, "condition": condition, "family": family, "replica": i,
            "seeds": seeds, "sweeps": sweeps, "profile": p, "profile_source": "seeded" if profile is None else "supplied",
            "CV": np.asarray(result.trace["CV"], dtype=np.float64),
            "TD": np.asarray(result.trace["TD"], dtype=np.float64),
            "marks": marks, "severity": severity_path(p, sweeps), "B0": local.B0,
            "current_sweep": local.current_sweep, "local_bytes": local.consumed_bytes,
            "emitted": emitted, "raw_base": raw, "emitted_sha256": hashlib.sha256(emitted).hexdigest(),
            "base_sha256": hashlib.sha256(raw).hexdigest(), "prefix_sha256": local.prefix_sha256,
            "start_sha256": result.initialization_protocol["initial_path_sha256"],
            "I-1": check_i1(marks, local.consumed_bytes, sweeps, False, current_sweep=local.current_sweep)}
