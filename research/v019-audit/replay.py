"""Independent end-to-end replay of one v0.19 Stage D pilot run from its seed labels.

Nothing here calls the pilot driver, the fault adapter, the sweep clock, the
random-source classes or the sampler class. The run is rebuilt from the frozen
definitions (research/v019-d-pilot/protocol.json, keys application, sweep_clock,
profile, fault_transforms):

- base, centroid and start streams come from NumPy's PCG64 (recount.pcg64_words);
- the sweep counter is the loop index of this file's own sweep loop, so the
  severity seen by a word and the byte marks are derived, not copied;
- the local Metropolis update and the centroid refresh are written out directly.

Two parts are mirrored rather than independent, because bit-identical output
needs the same floating-point operations: the start-path construction (the
Gaussian reference sampler's matrix algebra) and the elementary operations of
the sampler update (including `** 2`, which is libm `pow` and is not always
equal to a correctly rounded product). The comparison target is the stored
run, so any plumbing error between the definitions and the stored bytes shows.
"""
from __future__ import annotations

import hashlib
import importlib.util
import math
import sys
from pathlib import Path
from statistics import NormalDist
from typing import Any

import numpy as np

_HERE = Path(__file__).resolve().parent


def _load(name: str, path: Path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


rc = _load("v019_audit_recount", _HERE / "recount.py")
BEADS, BETA, DELTA = 16, 2.0, 0.4
_STD = NormalDist()


def start_path(start_seed: int) -> np.ndarray:
    """First of two exact finite-bead Gaussian paths (P application.start_call), mirrored algebra."""
    identity = np.eye(BEADS)
    difference = np.roll(identity, -1, axis=0) - identity
    precision = BETA / BEADS * identity + BEADS / BETA * (difference.T @ difference)
    cholesky = np.linalg.cholesky(np.linalg.inv(precision))
    count = 2 * BEADS
    pairs = (count + 1) // 2
    uniforms = (rc.pcg64_words(start_seed, 2 * pairs) >> np.uint64(11)).astype(float).reshape(pairs, 2) * 2.0**-53
    radius = np.sqrt(-2.0 * np.log1p(-uniforms[:, 0]))
    angles = 2.0 * math.pi * uniforms[:, 1]
    normals = np.column_stack((radius * np.cos(angles), radius * np.sin(angles))).ravel()[:count]
    return (normals.reshape(2, BEADS) @ cholesky.T)[0]


def replay(namespace: str, condition: str, family: str, replica: int, sweeps: int = rc.H,
           profile: dict[str, Any] | None = None, caps: dict[str, float] | None = None) -> dict[str, Any]:
    """Regenerate marks, traces, severity and stream hashes of one run.

    `profile` overrides the seeded draw (engineering fixtures only). A LAG1
    inverse-CDF domain error propagates: under the protocol it is a run failure.
    """
    seeds = rc.run_seeds(namespace, condition, replica)
    if family == "healthy":
        t0 = t_cap = s_max = None
    else:
        if profile is None:
            caps = rc.protocol_caps() if caps is None else caps
            profile = rc.draw_profile(seeds["profile"], caps[family])
        t0, t_cap, s_max = int(profile["t0"]), float(profile["T_cap"]), float(profile["s_max"])
    severity = rc.severity_path(t0, t_cap, s_max, sweeps)
    base = rc.pcg64_words(seeds["local"], 2 * BEADS * sweeps).tolist()  # at most 32 words per sweep
    refresh = (rc.pcg64_words(seeds["centroid"], 2 * sweeps) >> np.uint64(11)).astype(float) * 2.0**-53
    fault_rng = np.random.default_rng(seeds["fault"])
    start = start_path(seeds["start"])
    path = start.tolist()
    dtau = BETA / BEADS
    k_spring = 1.0 / dtau
    emitted: list[int] = []
    marks = np.empty(sweeps, dtype=np.int64)
    paths = np.empty((sweeps, BEADS), dtype=np.float64)
    pointer, previous, state, onset_word = 0, None, None, None
    root_beta = math.sqrt(BETA)
    active, s = False, 0.0

    def next_uniform() -> float:
        """Consume one local word under the current sweep's severity; return its 53-bit uniform."""
        nonlocal pointer, previous, state
        word = base[pointer]
        pointer += 1
        if active:
            if family == "DUP":
                if fault_rng.random() < s and previous is not None:
                    word = previous
            else:
                u = (word >> 11) * 2.0**-53
                if family == "SHIFT":
                    mapped = min(u + s, rc.TOP)
                else:
                    e = _STD.inv_cdf(u + 2.0**-54)
                    state = e if state is None else s * state + math.sqrt(1 - s * s) * e
                    mapped = min(_STD.cdf(state), rc.TOP)
                word = (int(mapped * 2.0**53) << 11) | (word & 0x7FF)
        previous = word
        emitted.append(word)
        return (word >> 11) * 2.0**-53

    for t in range(sweeps):
        active = t0 is not None and t >= t0
        s = float(severity[t])
        if active and onset_word is None:
            onset_word = pointer
        for j in range(BEADS):
            x_old = path[j]
            x_new = x_old + (next_uniform() * 2.0 - 1.0) * DELTA
            left, right = path[j - 1], path[(j + 1) % BEADS]
            d_spring = 0.5 * k_spring * ((x_new - left) ** 2 + (right - x_new) ** 2
                                         - (x_old - left) ** 2 - (right - x_old) ** 2)
            d_action = d_spring + dtau * 0.5 * (x_new ** 2 - x_old ** 2)
            if d_action <= 0.0:
                path[j] = x_new
            elif next_uniform() < math.exp(-d_action):
                path[j] = x_new
        marks[t] = 8 * pointer
        centre = math.fsum(value / BEADS for value in path)
        u, v = float(refresh[2 * t]), float(refresh[2 * t + 1])
        new_centre = math.sqrt(-2.0 * math.log1p(-u)) * math.cos(2.0 * math.pi * v) / root_beta
        path = [(value - centre) + new_centre for value in path]
        paths[t] = path
    springs = np.sum((np.roll(paths, -1, axis=-1) - paths)**2, axis=-1)
    x2 = np.mean(paths * paths, axis=-1)
    centroid = np.mean(paths, axis=-1)
    td = BEADS / (2.0 * BETA) - BEADS * springs / (2.0 * BETA**2) + 0.5 * x2
    # The sampler squares the scalar centroid with libm pow; an array square is a rounded product and
    # differs from it by one unit in the last place in about 0.09% of sweeps. Mirror the scalar power.
    centroid_squared = np.array([value ** 2 for value in centroid.tolist()], dtype=np.float64)
    cv = 1.0 / (2.0 * BETA) + x2 - 0.5 * centroid_squared
    words = np.array(emitted, dtype=np.uint64)
    stream = rc.words_to_bytes(words)
    raw = rc.words_to_bytes(np.array(base[:pointer], dtype=np.uint64))
    b0 = None if onset_word is None else 8 * onset_word
    return {"seeds": seeds, "profile": profile, "t0": t0, "severity": severity, "marks": marks, "CV": cv, "TD": td,
            "emitted_words": words, "local_bytes": 8 * pointer, "B0": b0,
            "emitted_sha256": hashlib.sha256(stream).hexdigest(), "base_sha256": hashlib.sha256(raw).hexdigest(),
            "prefix_sha256": None if b0 is None else hashlib.sha256(stream[:b0]).hexdigest(),
            "start_sha256": hashlib.sha256(start.astype("<f8").tobytes()).hexdigest()}
