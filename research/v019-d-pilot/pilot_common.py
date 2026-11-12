"""Stage D pilot: plan, paths, encodings, atomic files, module loading. protocol.json governs.

Shared by the pilot's stages (pilot_generate, pilot_analysis, pilot_forecast, run_pilot).
Nothing here draws a random number. Nothing here imports research/v019-audit: the
independent verification reads this pipeline's outputs and shares no code with it.

Frozen constants are module constants and cannot be passed in. A `Plan` carries only
sizes, so that tests can run a reduced plan; the command line builds the frozen plan
or the smoke plan and nothing else.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, Mapping

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
PROTOCOL_PATH = HERE / "protocol.json"
PILOT_NAMESPACE = "rnglab-v019-d-pilot"
SMOKE_NAMESPACE = "rnglab-v019-d-pilot-smoke"
RESERVED_NAMESPACES = (PILOT_NAMESPACE, SMOKE_NAMESPACE)
H = 40960
LENGTHS = (10240, 20480, 40960)
FAMILIES = ("SHIFT", "LAG1", "DUP")
MODELS = ("F0", "F0+", "F-B-lite", "F-AB-lite")
TARGETS = {"SHIFT": "frequency_signed_z", "LAG1": "lag1_feature", "DUP": "duplicate_rate"}
# protocol.json storage.parallelism
THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
PROTOCOL_WORKERS = 8
# protocol.json stop_criteria (thresholds fixed before generation)
COVERAGE_BAND = (0.92, 0.98)
V2D_MAXIMUM = 0.10
CRITERION_LENGTH = 20480
BALANCE_BAND = (0.2, 0.8)
ALREADY_INVALID_MAXIMUM = 0.5
MINIMUM_ELIGIBLE = 300
MINIMUM_REMAINING_MEDIAN = 5120
MINIMUM_SKILL = 0.10
MINIMUM_ELIGIBLE_EVALUATION = 240
CRITERION_4_LEVEL = 1 - 0.05 / 3
OTHER_LEVEL = 0.95
APPROVAL_DATE = "2026-11-10"
ENGINEERING_DEADLINE = "2026-12-22"
RESULTS = ("PASS", "FAIL", "INSUFFICIENT_INFORMATION", "UNINTERPRETABLE")


class PilotError(RuntimeError):
    """A condition under which the pilot must stop instead of continuing or repairing."""


@dataclass(frozen=True)
class Plan:
    """Sizes of a run plan. `kind` is 'pilot', 'smoke' or 'test'."""
    namespace: str
    conditions: tuple[tuple[str, str, int], ...]  # (condition key, family, replicas), protocol order
    fitting: range
    evaluation: range
    calibration_required: int
    threshold_order: int
    kind: str

    def __post_init__(self) -> None:
        if set(self.fitting) & set(self.evaluation):
            raise ValueError("fitting and evaluation replicas overlap")
        if self.kind not in ("pilot", "smoke", "test"):
            raise ValueError("plan kind must be pilot, smoke or test")
        smallest = min((n for _, family, n in self.conditions if family != "healthy"), default=0)
        if any(r < 0 or r >= smallest for r in (*self.fitting, *self.evaluation)):
            raise ValueError("a fitting or evaluation replica is outside the planned fault runs")

    def runs(self) -> Iterator[tuple[str, str, int]]:
        for condition, family, count in self.conditions:
            for replica in range(count):
                yield condition, family, replica

    def count(self, condition: str) -> int:
        return next(n for key, _, n in self.conditions if key == condition)

    def family(self, condition: str) -> str:
        return next(f for key, f, _ in self.conditions if key == condition)

    def total(self) -> int:
        return sum(n for _, _, n in self.conditions)

    def split(self, condition: str, replica: int) -> str:
        if condition == "calibration":
            return "calibration"
        if condition == "healthy":
            return "validation"
        return "fitting" if replica in self.fitting else "evaluation" if replica in self.evaluation else "unassigned"


_CONDITIONS = (("calibration", "healthy", 999), ("healthy", "healthy", 2000), ("fault/SHIFT", "SHIFT", 1000),
               ("fault/LAG1", "LAG1", 1000), ("fault/DUP", "DUP", 1000))
FROZEN_PLAN = Plan(PILOT_NAMESPACE, _CONDITIONS, range(0, 600), range(600, 1000), 999, 950, "pilot")
# protocol.json smoke: two runs per condition, engineering and timing only
SMOKE_PLAN = Plan(SMOKE_NAMESPACE, tuple((key, family, 2) for key, family, _ in _CONDITIONS),
                  range(0, 1), range(1, 2), 999, 950, "smoke")


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def protocol_sha256(check: bool = True) -> str:
    """Hash of protocol.json; refuses to continue when it differs from protocol.sha256."""
    actual = sha256_file(PROTOCOL_PATH)
    if check and actual != (HERE / "protocol.sha256").read_text().split()[0]:
        raise PilotError("protocol.json does not match protocol.sha256")
    return actual


def canonical_json(value: Any) -> bytes:
    """protocol.json storage.json: sorted keys, indent 2, UTF-8, no NaN or Infinity, final newline."""
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def canonical_jsonl(rows: list[Any]) -> bytes:
    return "".join(json.dumps(row, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n" for row in rows).encode("utf-8")


def fsync_directory(directory: Path) -> None:
    descriptor = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write(path: Path, data: bytes, *, overwrite: bool = False) -> str:
    """Write through a temporary file in the same directory, flush, rename. Returns the SHA-256.

    An existing target is an error unless overwrite=True: derived files are produced once.
    """
    path = Path(path)
    if path.exists() and not overwrite:
        raise PilotError(f"{path} already exists; it is not overwritten")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with open(temporary, "wb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    fsync_directory(path.parent)
    return sha256_bytes(data)


def read_json(path: Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def run_paths(condition: str, replica: int) -> tuple[Path, Path]:
    """(trace path, record path), relative to the output root (protocol.json storage)."""
    return Path("traces") / condition / f"{replica:05d}.npz", Path("results/runs") / condition / f"{replica:05d}.json"


def load_module(name: str, path: Path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@dataclass
class Context:
    """Where a stage works and which modules it calls.

    measure is the pilot's measurement layer (research/v019-d-pilot/measure.py). It is passed
    in, never searched for: the command line loads the file in this directory and stops if it
    is missing; tests pass the module they were given.
    """
    plan: Plan
    root: Path
    measure: Any
    dlib: Any
    c_analyze: Any

    @property
    def h1lib(self) -> Any:
        return self.dlib.h1lib

    def record_path(self, condition: str, replica: int) -> Path:
        return self.root / run_paths(condition, replica)[1]

    def trace_path(self, condition: str, replica: int) -> Path:
        return self.root / run_paths(condition, replica)[0]

    def result(self, name: str) -> Path:
        return self.root / "results" / name


def production_context(plan: Plan, root: Path) -> Context:
    """Context with the modules of this repository. Stops if the measurement layer is not in the tree."""
    measure_path = HERE / "measure.py"
    if not measure_path.is_file():
        raise PilotError("research/v019-d-pilot/measure.py is not in this tree: the pilot cannot run without "
                         "its measurement layer")
    measure = load_module("d_pilot_measure_production", measure_path)
    return Context(plan, Path(root), measure, measure.dlib,
                   load_module("d_pilot_c_analyze", REPO / "research/v018-h1-detection/analyze.py"))


def approved_amendments(directory: Path = HERE) -> dict[str, Any]:
    """Amendments the author has approved, from amendments.json next to the protocol.

    The file does not exist until the author writes it. Each entry names its document and
    that document's SHA-256; a missing or altered document is an error, not a silent fallback.
    """
    path = directory / "amendments.json"
    if not path.is_file():
        return {}
    amendments = read_json(path)
    for key, entry in amendments.items():
        document = directory / entry["document"]
        if not document.is_file() or sha256_file(document) != entry["sha256"]:
            raise PilotError(f"amendment {key}: document missing or altered")
        if not entry.get("approved_by") or not entry.get("approved_on"):
            raise PilotError(f"amendment {key}: approval is not recorded")
    return amendments


def software_versions() -> dict[str, str]:
    import platform

    import numpy
    import scipy
    return {"python": platform.python_version(), "numpy": numpy.__version__, "scipy": scipy.__version__,
            "libc": " ".join(platform.libc_ver()), "machine": platform.machine()}


def finite_or_none(value: Any) -> float | None:
    import math
    return float(value) if value is not None and isinstance(value, (int, float)) and math.isfinite(value) else None
