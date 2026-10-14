"""Read-only resource checks, bounded admission, and reproducibility receipts.

There is deliberately no cloud provisioning, package installation, remote upload,
or automatic scientific campaign execution in this module.
"""
from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import re
import resource
import shutil
import subprocess
import threading
from typing import Any, Callable, Mapping

import numpy as np


MIB = 1024**2
GIB = 1024**3


def _utc_now():
    return datetime.now(timezone.utc).isoformat()


def _canonical(value) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def canonical_digest(value) -> str:
    return hashlib.sha256(_canonical(value)).hexdigest()


def _small_text(path: str | Path, limit: int = 65_536):
    try:
        with Path(path).open("r") as stream:
            text = stream.read(limit + 1)
        return text if len(text) <= limit else None
    except (OSError, UnicodeError):
        return None


def _positive_int(text):
    try:
        value = int(text)
        return value if value > 0 and value < 2**60 else None
    except (TypeError, ValueError):
        return None


def _cgroup_directory():
    text = _small_text("/proc/self/cgroup") or ""
    for line in text.splitlines():
        parts = line.split(":", 2)
        if len(parts) == 3 and parts[0] == "0":
            root = Path("/sys/fs/cgroup")
            candidate = (root / parts[2].lstrip("/")).resolve()
            if candidate.is_relative_to(root) and candidate.exists():
                return candidate
    return Path("/sys/fs/cgroup")


def inventory_resources(*, disk_path: str | Path = ".", gpu_timeout_seconds: float = 3.0,
                        run: Callable = subprocess.run) -> dict[str, Any]:
    """Bounded local inspection; no network requests or credential discovery.

    A snapshot is not a reservation. Memory availability and GPU free memory can
    change when unrelated applications run; refresh before admitting real work.
    """
    if not 0 < gpu_timeout_seconds <= 5:
        raise ValueError("GPU inventory timeout must be in (0,5] seconds")
    host_cores = os.cpu_count()
    try:
        affinity_cores = len(os.sched_getaffinity(0))
    except (AttributeError, OSError):
        affinity_cores = host_cores
    cpu_model = platform.processor() or None
    for line in (_small_text("/proc/cpuinfo") or "").splitlines():
        if line.lower().startswith("model name"):
            cpu_model = line.split(":", 1)[-1].strip()
            break
    cgroup = _cgroup_directory()
    quota_cores = None
    cpu_max = (_small_text(cgroup / "cpu.max") or "").split()
    if len(cpu_max) == 2 and cpu_max[0] != "max":
        quota, period = _positive_int(cpu_max[0]), _positive_int(cpu_max[1])
        if quota and period:
            quota_cores = quota / period
    candidates = [float(v) for v in (host_cores, affinity_cores, quota_cores) if v is not None]
    effective_cores = min(candidates) if candidates else None
    memory = {}
    for line in (_small_text("/proc/meminfo") or "").splitlines():
        key, _, value = line.partition(":")
        if key in ("MemTotal", "MemAvailable"):
            pieces = value.split()
            if pieces and _positive_int(pieces[0]):
                memory[key] = int(pieces[0]) * 1024
    memory_limit = _positive_int((_small_text(cgroup / "memory.max") or "").strip())
    memory_current_text = (_small_text(cgroup / "memory.current") or "").strip()
    try:
        memory_current = int(memory_current_text)
    except ValueError:
        memory_current = None
    cgroup_available = max(0, memory_limit - memory_current) if memory_limit is not None and memory_current is not None else None
    available_candidates = [v for v in (memory.get("MemAvailable"), cgroup_available) if v is not None]
    disk = shutil.disk_usage(disk_path)
    gpu = {"status": "UNAVAILABLE", "devices": [], "probe": "nvidia-smi"}
    if shutil.which("nvidia-smi"):
        try:
            result = run(["nvidia-smi", "--query-gpu=index,name,memory.total,memory.free",
                          "--format=csv,noheader,nounits"], capture_output=True, text=True,
                         timeout=gpu_timeout_seconds, check=False)
            if result.returncode == 0:
                devices = []
                for row in csv.reader(result.stdout.splitlines()):
                    if len(row) != 4:
                        raise ValueError("unexpected GPU inventory format")
                    index, name, total, free = (cell.strip() for cell in row)
                    devices.append({"index": index, "name": name,
                                    "total_bytes": int(total) * MIB, "available_bytes": int(free) * MIB})
                gpu.update(status="AVAILABLE" if devices else "UNAVAILABLE", devices=devices)
            else:
                gpu["status"] = "PROBE_FAILED"  # Do not copy arbitrary process stderr.
        except subprocess.TimeoutExpired:
            gpu["status"] = "PROBE_TIMEOUT"
        except (OSError, ValueError):
            gpu["status"] = "PROBE_FAILED"
    limits = {}
    for label, constant in (("open_files", "RLIMIT_NOFILE"), ("processes", "RLIMIT_NPROC"),
                            ("address_space_bytes", "RLIMIT_AS")):
        if hasattr(resource, constant):
            soft, hard = resource.getrlimit(getattr(resource, constant))
            limits[label] = {"soft": None if soft == resource.RLIM_INFINITY else soft,
                             "hard": None if hard == resource.RLIM_INFINITY else hard}
    return {"schema_version": 1, "observed_at_utc": _utc_now(),
            "operating_system": {"system": platform.system(), "release": platform.release(),
                                 "machine": platform.machine(), "python": platform.python_version()},
            "cpu": {"model": cpu_model, "logical_cores": host_cores, "affinity_cores": affinity_cores,
                    "cgroup_quota_cores": quota_cores, "effective_cores": effective_cores},
            "memory": {"physical_total_bytes": memory.get("MemTotal"),
                       "host_available_bytes": memory.get("MemAvailable"),
                       "cgroup_limit_bytes": memory_limit, "cgroup_current_bytes": memory_current,
                       "effective_available_bytes": min(available_candidates) if available_candidates else None},
            "gpu": gpu, "disk": {"total_bytes": disk.total, "used_bytes": disk.used, "free_bytes": disk.free},
            "process_limits": limits, "network": {"status": "NOT_PROBED"},
            "limitations": ["Read-only local snapshot, not a reservation or process-lifetime guarantee.",
                            "GPU visibility is not scientific CPU/GPU execution parity.",
                            "Sandbox/network/account permissions are not inferred from local hardware."]}


def backend_capabilities():
    gpu_packages = {}
    for name in ("numba", "torch", "cupy"):
        try:
            gpu_packages[name] = importlib.util.find_spec(name) is not None
        except (ImportError, ValueError):
            gpu_packages[name] = False
    return {"python_cpu": "AVAILABLE", "cpu_compiler": shutil.which("cc") or shutil.which("gcc"),
            "compiled_cpu_parity": "NOT_EXECUTED", "gpu_packages_present": gpu_packages,
            "gpu_scientific_backend": "NOT_EXECUTED" if any(gpu_packages.values()) else "UNAVAILABLE",
            "gpu_scientific_parity": "NOT_EXECUTED",
            "prior_local_model_inference": "HISTORICAL_REPORT_ONLY_SEPARATE_FROM_SCIENTIFIC_BACKEND"}


@dataclass(frozen=True)
class JobRequest:
    job_id: str
    cpu_cores: float
    ram_bytes: int
    gpu_bytes: int = 0
    gpu_index: str | None = None

    def __post_init__(self):
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,100}", self.job_id):
            raise ValueError("invalid bounded job identifier")
        if isinstance(self.cpu_cores, bool) or not math.isfinite(self.cpu_cores) or self.cpu_cores <= 0:
            raise ValueError("CPU cores must be finite and positive")
        for label in ("ram_bytes", "gpu_bytes"):
            value = getattr(self, label)
            if type(value) is not int or value < 0 or label == "ram_bytes" and value == 0:
                raise ValueError("invalid memory reservation")
        if self.gpu_index is not None and not re.fullmatch(r"[0-9]+", str(self.gpu_index)):
            raise ValueError("GPU index must be a numeric local index")


class JobAdmission:
    """Thread-safe reservation ledger; the caller owns execution and supervision.

    Admissions neither grant scientific authorization nor enforce OS memory/CPU
    limits. Combine with process/cgroup limits in an approved execution service.
    """

    def __init__(self, inventory: Mapping[str, Any], *, max_concurrent: int = 2,
                 reserve_cpu: float = 1.0, reserve_ram_bytes: int = 512 * MIB,
                 reserve_gpu_bytes: int = 512 * MIB,
                 inventory_provider: Callable | None = None):
        if type(max_concurrent) is not int or max_concurrent < 1:
            raise ValueError("positive concurrency bound required")
        if not math.isfinite(reserve_cpu) or reserve_cpu < 0 or type(reserve_ram_bytes) is not int or reserve_ram_bytes < 0:
            raise ValueError("invalid headroom reservation")
        if type(reserve_gpu_bytes) is not int or reserve_gpu_bytes < 0:
            raise ValueError("invalid GPU headroom reservation")
        self.max_concurrent, self.reserve_cpu, self.reserve_ram_bytes = max_concurrent, reserve_cpu, reserve_ram_bytes
        self.reserve_gpu_bytes = reserve_gpu_bytes
        self._inventory = json.loads(_canonical(dict(inventory)))
        self._inventory_provider = inventory_provider
        self._running: dict[str, JobRequest] = {}
        self._lock = threading.RLock()

    def _capacity(self):
        cores = self._inventory.get("cpu", {}).get("effective_cores")
        ram = self._inventory.get("memory", {}).get("effective_available_bytes")
        return (max(0.0, cores - self.reserve_cpu) if cores is not None else None,
                max(0, ram - self.reserve_ram_bytes) if ram is not None else None)

    def refresh(self, inventory: Mapping[str, Any]):
        with self._lock:
            self._inventory = json.loads(_canonical(dict(inventory)))

    def admit(self, request: JobRequest):
        if not isinstance(request, JobRequest):
            raise TypeError("JobRequest required")
        with self._lock:
            if self._inventory_provider:
                self.refresh(self._inventory_provider())
            if request.job_id in self._running:
                prior = self._running[request.job_id]
                matches = (request.cpu_cores, request.ram_bytes, request.gpu_bytes) == (prior.cpu_cores, prior.ram_bytes, prior.gpu_bytes)
                matches = matches and (request.gpu_index is None or str(request.gpu_index) == str(prior.gpu_index))
                return {"status": "ALREADY_ADMITTED" if matches else "REJECTED_JOB_ID_CONFLICT",
                        "job_id": request.job_id, "reservation": asdict(prior)}
            cores, ram = self._capacity()
            if cores is None or ram is None:
                return {"status": "DEFERRED", "job_id": request.job_id, "reason": "RESOURCE_INVENTORY_INCOMPLETE"}
            gpu = self._inventory.get("gpu", {})
            devices = {str(d["index"]): {**d, "available_bytes": max(0, d["available_bytes"] - self.reserve_gpu_bytes)}
                       for d in gpu.get("devices", [])}
            running = list(self._running.values())
            assigned_gpu = str(request.gpu_index) if request.gpu_index is not None else None
            if request.gpu_bytes:
                if assigned_gpu is None:
                    free = {index: device["available_bytes"] - sum(r.gpu_bytes for r in running if str(r.gpu_index) == index)
                            for index, device in devices.items()}
                    eligible = [index for index, value in free.items() if value >= request.gpu_bytes]
                    assigned_gpu = max(eligible, key=free.get) if eligible else None
                if assigned_gpu is None or assigned_gpu not in devices:
                    return {"status": "DEFERRED", "job_id": request.job_id, "reason": "NO_ELIGIBLE_GPU"}
                if request.gpu_bytes > devices[assigned_gpu]["available_bytes"]:
                    return {"status": "DEFERRED", "job_id": request.job_id, "reason": "GPU_MEMORY_LIMIT"}
            if request.cpu_cores > cores or request.ram_bytes > ram:
                return {"status": "DEFERRED", "job_id": request.job_id, "reason": "HOST_CAPACITY_LIMIT"}
            if len(running) >= self.max_concurrent:
                return {"status": "DEFERRED", "job_id": request.job_id, "reason": "CONCURRENCY_LIMIT"}
            if sum(r.cpu_cores for r in running) + request.cpu_cores > cores:
                return {"status": "DEFERRED", "job_id": request.job_id, "reason": "CPU_RESERVED"}
            if sum(r.ram_bytes for r in running) + request.ram_bytes > ram:
                return {"status": "DEFERRED", "job_id": request.job_id, "reason": "RAM_RESERVED"}
            if request.gpu_bytes and (sum(r.gpu_bytes for r in running if str(r.gpu_index) == assigned_gpu)
                                      + request.gpu_bytes > devices[assigned_gpu]["available_bytes"]):
                return {"status": "DEFERRED", "job_id": request.job_id, "reason": "GPU_RESERVED"}
            admitted = JobRequest(request.job_id, request.cpu_cores, request.ram_bytes,
                                  request.gpu_bytes, assigned_gpu)
            self._running[request.job_id] = admitted
            return {"status": "ADMITTED", "job_id": request.job_id, "reservation": asdict(admitted),
                    "scientific_authorization": "NOT_GRANTED_BY_RESOURCE_ADMISSION"}

    def release(self, job_id: str):
        with self._lock:
            return self._running.pop(job_id, None) is not None

    def select_pending(self, requests):
        """One finite pass through at most 1,000 queued requests; no job launch."""
        if not isinstance(requests, (list, tuple)) or len(requests) > 1000:
            raise ValueError("pending selection must be a finite list of at most 1000 jobs")
        return [self.admit(request) for request in requests]

    def snapshot(self):
        with self._lock:
            cores, ram = self._capacity()
            return {"schema_version": 1, "max_concurrent": self.max_concurrent,
                    "cpu_capacity_after_headroom": cores, "ram_capacity_after_headroom": ram,
                    "gpu_headroom_bytes_per_device": self.reserve_gpu_bytes,
                    "active": [asdict(value) for value in self._running.values()],
                    "note": "Reconcile live PIDs after process restart; stale reservations are not running jobs."}


def file_digest(path: str | Path, *, max_bytes: int = 64 * MIB):
    """Streaming digest with an explicit per-file cap and mutation detection."""
    if type(max_bytes) is not int or max_bytes < 0:
        raise ValueError("invalid artifact hashing budget")
    path = Path(path)
    with path.open("rb") as stream:
        before = os.fstat(stream.fileno())
        if before.st_size > max_bytes:
            raise ValueError("artifact exceeds checkpoint hashing budget")
        digest = hashlib.sha256()
        consumed = 0
        while True:
            chunk = stream.read(min(65_536, max_bytes - consumed + 1))
            if not chunk:
                break
            consumed += len(chunk)
            if consumed > max_bytes:
                raise ValueError("artifact exceeded checkpoint hashing budget")
            digest.update(chunk)
        after = os.fstat(stream.fileno())
    if (before.st_size, before.st_mtime_ns, before.st_ino) != (after.st_size, after.st_mtime_ns, after.st_ino):
        raise ValueError("artifact changed while hashing")
    return {"sha256": digest.hexdigest(), "bytes": consumed}


_SECRET_FIELDS = {"api_key", "apikey", "password", "private_key", "private_key_bytes",
                  "secret", "access_token", "refresh_token", "key_hex", "value_hex"}


def _check_public_json(value, *, depth=0):
    if depth > 32:
        raise ValueError("checkpoint nesting exceeds limit")
    if isinstance(value, dict):
        for key, item in value.items():
            if not isinstance(key, str) or key.lower() in _SECRET_FIELDS:
                raise ValueError("checkpoint contains a prohibited credential/private-state field")
            _check_public_json(item, depth=depth + 1)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _check_public_json(item, depth=depth + 1)
    elif value is not None and not isinstance(value, (str, int, float, bool)):
        raise ValueError("checkpoint must contain plain JSON research state")
    elif isinstance(value, float) and not math.isfinite(value):
        raise ValueError("checkpoint values must be finite")


class CheckpointStore:
    """Atomic bounded public research checkpoints; no pickle or secret state.

    Resume verifies exact source/protocol/code manifests. Sampler continuation
    remains the caller's responsibility; this receipt does not assert parity.
    """

    def __init__(self, *, max_checkpoint_bytes: int = MIB, max_artifact_bytes: int = 64 * MIB):
        if type(max_checkpoint_bytes) is not int or not 1 <= max_checkpoint_bytes <= 16 * MIB:
            raise ValueError("checkpoint bound must be in [1 byte,16 MiB]")
        if type(max_artifact_bytes) is not int or max_artifact_bytes < 0:
            raise ValueError("invalid artifact hashing budget")
        self.max_checkpoint_bytes, self.max_artifact_bytes = max_checkpoint_bytes, max_artifact_bytes

    def save(self, path: str | Path, *, state: Mapping, protocol: Mapping,
             source_manifest: Mapping, code_manifest: Mapping, backend: str,
             completed_units: int, artifacts: Mapping[str, str | Path] | None = None):
        if type(completed_units) is not int or completed_units < 0:
            raise ValueError("completed_units must be a nonnegative integer")
        if artifacts is not None and len(artifacts) > 100:
            raise ValueError("at most 100 checkpoint artifacts allowed")
        artifact_manifest = {}
        remaining = self.max_artifact_bytes
        for label, file in (artifacts or {}).items():
            artifact_manifest[label] = file_digest(file, max_bytes=remaining)
            remaining -= artifact_manifest[label]["bytes"]
        payload = {"schema_version": 1, "classification": "public_research_checkpoint",
                   "created_at_utc": _utc_now(), "backend": backend, "completed_units": completed_units,
                   "state": dict(state), "protocol_sha256": canonical_digest(dict(protocol)),
                   "source_manifest_sha256": canonical_digest(dict(source_manifest)),
                   "code_manifest_sha256": canonical_digest(dict(code_manifest)),
                   "artifacts": artifact_manifest}
        _check_public_json(payload)
        encoded = _canonical(payload)
        if len(encoded) > self.max_checkpoint_bytes:
            raise ValueError("checkpoint exceeds serialization budget")
        path = Path(path)
        temporary = path.with_name(path.name + ".new")
        # Exclusive creation prevents accidental concurrent writer overwrite.
        try:
            fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            with os.fdopen(fd, "wb") as stream:
                stream.write(encoded)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, path)
        except BaseException:
            # Only unlink our explicitly named temporary if this call created it.
            if "fd" in locals():
                temporary.unlink(missing_ok=True)
            raise
        return {"status": "CHECKPOINT_SAVED", "checkpoint_sha256": hashlib.sha256(encoded).hexdigest(),
                "checkpoint_bytes": len(encoded), "completed_units": completed_units}

    def load(self, path: str | Path, *, protocol: Mapping, source_manifest: Mapping,
             code_manifest: Mapping, artifacts: Mapping[str, str | Path] | None = None,
             expected_checkpoint_sha256: str | None = None):
        path = Path(path)
        with path.open("rb") as stream:
            encoded = stream.read(self.max_checkpoint_bytes + 1)
        if len(encoded) > self.max_checkpoint_bytes:
            raise ValueError("checkpoint exceeds input budget")
        checkpoint_digest = hashlib.sha256(encoded).hexdigest()
        if expected_checkpoint_sha256 is not None and checkpoint_digest != expected_checkpoint_sha256:
            raise ValueError("checkpoint trusted receipt digest mismatch")
        payload = json.loads(encoded)
        _check_public_json(payload)
        if payload.get("schema_version") != 1 or payload.get("classification") != "public_research_checkpoint":
            raise ValueError("unsupported checkpoint schema")
        if type(payload.get("completed_units")) is not int or payload["completed_units"] < 0 or not isinstance(payload.get("state"), dict):
            raise ValueError("invalid checkpoint state/counters")
        for field, expected in (("protocol_sha256", protocol), ("source_manifest_sha256", source_manifest),
                                ("code_manifest_sha256", code_manifest)):
            if payload.get(field) != canonical_digest(dict(expected)):
                raise ValueError("checkpoint " + field + " mismatch")
        provided = artifacts or {}
        if set(provided) != set(payload.get("artifacts", {})):
            raise ValueError("checkpoint artifact set mismatch")
        remaining = self.max_artifact_bytes
        for label, file in provided.items():
            actual = file_digest(file, max_bytes=remaining)
            remaining -= actual["bytes"]
            if actual != payload["artifacts"][label]:
                raise ValueError("checkpoint artifact digest mismatch")
        return {"status": "RESUME_VERIFIED" if expected_checkpoint_sha256 is not None else "RESUME_MANIFESTS_VERIFIED",
                "payload": payload, "trusted_receipt_verified": expected_checkpoint_sha256 is not None,
                "checkpoint_sha256": checkpoint_digest,
                "scientific_confirmation": False}


def compare_execution_parity(reference_output, candidate_output, *, reference_state: Mapping,
                             candidate_state: Mapping, protocol: Mapping, source_manifest: Mapping,
                             reference_backend: str = "python_cpu", candidate_backend: str = "compiled_cpu",
                             atol: float = 1e-12, rtol: float = 1e-12):
    """Numeric output and final-state parity for one explicitly named fixture."""
    if not math.isfinite(atol) or not math.isfinite(rtol) or min(atol, rtol) < 0:
        raise ValueError("finite nonnegative parity tolerances required")
    a, b = np.asarray(reference_output, dtype=float), np.asarray(candidate_output, dtype=float)
    if a.size > 1_000_000 or b.size > 1_000_000:
        raise ValueError("parity fixture exceeds element budget")
    finite = bool(np.isfinite(a).all() and np.isfinite(b).all())
    shape_match = a.shape == b.shape
    outputs_match = bool(shape_match and finite and np.allclose(a, b, atol=atol, rtol=rtol))
    state_match = _state_close(reference_state, candidate_state, atol=atol, rtol=rtol)
    return {"schema_version": 1, "status": "PARITY_PASSED" if outputs_match and state_match else "PARITY_FAILED",
            "reference_backend": reference_backend, "candidate_backend": candidate_backend,
            "output_shape": list(a.shape), "output_shape_match": shape_match, "finite_outputs": finite,
            "outputs_match": outputs_match, "state_match": state_match,
            "maximum_absolute_output_difference": float(np.max(np.abs(a - b))) if shape_match and finite and a.size else None,
            "absolute_tolerance": atol, "relative_tolerance": rtol,
            "protocol_sha256": canonical_digest(dict(protocol)),
            "source_manifest_sha256": canonical_digest(dict(source_manifest)),
            "scientific_confirmation": False,
            "scope": "Named finite software fixture; no general backend or scientific validation claim"}


def _state_close(left, right, *, atol, rtol):
    if isinstance(left, Mapping) and isinstance(right, Mapping):
        return set(left) == set(right) and all(_state_close(left[k], right[k], atol=atol, rtol=rtol) for k in left)
    if isinstance(left, (list, tuple, np.ndarray)) and isinstance(right, (list, tuple, np.ndarray)):
        a, b = np.asarray(left), np.asarray(right)
        if a.shape != b.shape or a.size > 1_000_000:
            return False
        try:
            return bool(np.isfinite(a).all() and np.isfinite(b).all() and np.allclose(a, b, atol=atol, rtol=rtol))
        except (TypeError, ValueError):
            return False
    if isinstance(left, (int, float, np.number)) and isinstance(right, (int, float, np.number)):
        return bool(np.isfinite(left) and np.isfinite(right) and np.isclose(left, right, atol=atol, rtol=rtol))
    return type(left) is type(right) and left == right


@dataclass(frozen=True)
class SSHConfig:
    host: str
    user: str
    identity_file: str
    known_hosts_file: str
    port: int = 22
    timeout_seconds: int = 5

    def __post_init__(self):
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.:-]{0,252}", self.host) or self.host.startswith("-"):
            raise ValueError("explicit hostname/IP required")
        if not re.fullmatch(r"[a-z_][a-z0-9_-]{0,31}", self.user):
            raise ValueError("explicit SSH username required")
        if type(self.port) is not int or not 1 <= self.port <= 65_535:
            raise ValueError("invalid SSH port")
        if type(self.timeout_seconds) is not int or not 1 <= self.timeout_seconds <= 5:
            raise ValueError("SSH timeout must be 1–5 seconds")
        for value in (self.identity_file, self.known_hosts_file):
            if not isinstance(value, str) or not value or "\0" in value or "\n" in value:
                raise ValueError("explicit local key and pinned known-hosts paths required")


def read_only_ssh_preflight(config: SSHConfig, *, run: Callable = subprocess.run):
    """One requested host, pinned identity, fixed read-only command; no upload.

    Key contents are never read by Python or included in returned receipts. No
    interactive login, password fallback, ssh-keyscan, sudo or installation occurs.
    """
    if not isinstance(config, SSHConfig):
        raise TypeError("SSHConfig required")
    public = {"host": config.host, "user": config.user, "port": config.port,
              "timeout_seconds": config.timeout_seconds, "operation": "READ_ONLY_PREFLIGHT"}
    key_path, host_path = Path(config.identity_file).expanduser(), Path(config.known_hosts_file).expanduser()
    if not key_path.is_file():
        return {**public, "status": "IDENTITY_FILE_UNAVAILABLE"}
    if not host_path.is_file() or host_path.stat().st_size == 0:
        return {**public, "status": "AWAITING_PINNED_HOST_KEY"}
    if os.name == "posix" and key_path.stat().st_mode & 0o077:
        return {**public, "status": "IDENTITY_PERMISSIONS_TOO_OPEN"}
    if os.name == "posix" and host_path.stat().st_mode & 0o022:
        return {**public, "status": "PINNED_HOSTS_PERMISSIONS_TOO_OPEN"}
    ssh = shutil.which("ssh")
    if ssh is None:
        return {**public, "status": "SSH_UNAVAILABLE"}
    command = [ssh, "-F", "/dev/null", "-o", "BatchMode=yes", "-o", "StrictHostKeyChecking=yes",
               "-o", "UserKnownHostsFile=" + str(host_path), "-o", "IdentitiesOnly=yes",
               "-o", "GlobalKnownHostsFile=/dev/null", "-o", "UpdateHostKeys=no", "-o", "VerifyHostKeyDNS=no",
               "-o", "ForwardAgent=no", "-o", "ClearAllForwardings=yes", "-o", "PermitLocalCommand=no",
               "-o", "RequestTTY=no",
               "-o", "PasswordAuthentication=no", "-o", "KbdInteractiveAuthentication=no",
               "-o", "ConnectTimeout=" + str(config.timeout_seconds), "-o", "ConnectionAttempts=1",
               "-i", str(key_path), "-p", str(config.port), config.user + "@" + config.host,
               "LC_ALL=C uname -srm; getconf _NPROCESSORS_ONLN; getconf PAGESIZE; getconf _PHYS_PAGES; command -v python3"]
    try:
        result = run(command, capture_output=True, text=True, timeout=config.timeout_seconds, check=False)
    except subprocess.TimeoutExpired:
        return {**public, "status": "UNAVAILABLE_TIMEOUT", "provisioned": False}
    except OSError:
        return {**public, "status": "SSH_EXECUTION_FAILED", "provisioned": False}
    if result.returncode != 0:
        return {**public, "status": "UNAVAILABLE_OR_AUTHENTICATION_FAILED", "returncode": result.returncode,
                "provisioned": False, "note": "Host availability/authentication not established; stderr intentionally omitted"}
    lines = result.stdout.splitlines()
    # Never echo arbitrary remote output: parse only the fixed expected fields.
    if (len(lines) != 5 or len(lines[0]) > 256 or not all(line.isdigit() for line in lines[1:4])
            or not re.fullmatch(r"/[A-Za-z0-9_./+-]{1,200}", lines[4])):
        return {**public, "status": "UNEXPECTED_PREFLIGHT_OUTPUT", "provisioned": False}
    return {**public, "status": "READ_ONLY_PREFLIGHT_PASSED", "provisioned": False,
            "remote_system": lines[0], "logical_cpu_cores": int(lines[1]),
            "physical_ram_bytes": int(lines[2]) * int(lines[3]), "python_executable": lines[4],
            "scientific_execution_parity": "NOT_EXECUTED",
            "note": "Reachability and inventory only; no RNGLab deployment or parity established"}


def deployment_plan(*, provider: str | None = None, approved_budget: Mapping | None = None,
                    local_only: bool = True):
    """Reviewable plan only, explicitly refusing to invent a provisioned host."""
    if type(local_only) is not bool:
        raise ValueError("local_only policy must be explicit boolean")
    configured = bool(provider and approved_budget)
    return {"schema_version": 1,
            "status": "LOCAL_ONLY_USER_REQUESTED" if local_only else
            ("AWAITING_EXPLICIT_PROVISIONING_AUTHORIZATION" if configured else "AWAITING_USER_PROVIDER_AND_BUDGET"),
            "execution_policy": "local_only" if local_only else "future_remote_proposal",
            "provider": provider, "approved_budget": dict(approved_budget) if approved_budget else None,
            "provisioned_hosts": [], "provisioning_executed": False,
            "remote_work_status": "DEFERRED_BY_USER_LOCAL_ONLY_DECISION" if local_only else "AWAITING_CONFIGURATION_AND_APPROVAL",
            "planned_steps": ["Keep all current execution local; remote steps below are documentation for future explicit authorization.",
                              "Select provider, region, CPU/RAM/disk and maximum billed lifetime with the user.",
                              "Obtain explicit provisioning and spending authorization.",
                              "Pin host identity through a trusted channel and run read-only inventory.",
                              "Review minimal source/dependency transfer and credentials policy.",
                              "Run bounded CPU fixture parity, checkpoint/resume and resource-admission checks.",
                              "Require separate approval before scientific confirmation or persistent services."],
            "historical_host": {"host": "80.158.39.137", "user": "ubuntu",
                                "last_observation": "Previously attempted read-only SSH timed out after 5 seconds; current status unverified",
                                "v011_deployment_receipt": "NOT_FOUND"},
            "cloud_api_calls": 0, "api_cost": 0,
            "gpu_scientific_parity": "NOT_EXECUTED"}
