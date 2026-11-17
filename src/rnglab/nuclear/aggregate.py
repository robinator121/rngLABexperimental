"""RNGLab v0.20 nuclear module, milestone B2: aggregated (population-level) decay.

Populations are counts per nuclide; no nucleus is allocated. Over an interval
of length delta every nucleus follows one decay-path category of its start
nuclide (``paths``), independently of all others, so the category counts of
the n_i nuclei in nuclide i are Multinomial(n_i, q_i(delta)) and independent
across i. Sampling them gives exactly the law of the counts at the grid times
and of the decays per (nuclide, channel) in each interval. Photons are then
Binomial(decays through the channel, photon yield), and neutrons equal the
beta-delayed-neutron decays.

The multinomials are sampled as conditional binomials, p_k = q_k / T_k with
T_k = sum_{j >= k} q_j (decimal, ``IntervalLaw.conditional``); each binomial by
inversion at the midpoint of its uniform's 2**-53 cell, one uniform per part of
at most 2**50 trials (``binomial`` states which law that realises). A binomial
whose probability exceeds 1/2 is sampled through its complement, computed in
decimal.

Conventions and limits
----------------------
- The length of interval k is the float64 difference t[k+1] - t[k]; the laws
  are evaluated at that value (a relative difference of at most 2**-53 from
  the exact difference of the two float times).
- Total population at most 2**63 - 1; the neutron count of one interval must
  also fit (schemes with several beta-delayed-neutron steps on one path lower
  the limit accordingly). Both are checked before any random number is read.
- No state is shared with callers: the time grid is copied, plan arrays and
  returned arrays are read-only, a sink receives detached copies of every
  interval's outputs, and a checkpoint is built from fresh Python objects.
- A checkpoint carries the SHA-256 of its own contents and the identity of
  the stream (generator, seed, absolute word position). A resume refuses a
  checkpoint whose checksum, types, layout, word count or conservation is
  wrong, and a stream with another seed or at another position. The checksum
  is an integrity check, not an authentication.

What aggregation loses (compared with Layer A)
----------------------------------------------
Nucleus identities and lineage, event timestamps within an interval (only the
interval of each decay is known), the order of decays within an interval, and
per-nucleus histories. It also changes how random numbers enter: a fixed
number of uniforms per interval, independent of the population, instead of
3*D words per nucleus. Agreement with Layer A is statistical (tests that do not
reject, with a stated sensitivity), for healthy generators; it does not carry
over to how a defective generator distorts each engine.

Random-number layout "nuclear-aggregate-layout-v1"
--------------------------------------------------
Per interval, in this order: for every nuclide in scheme order with K >= 2
categories, one slot per category k < K - 1; then one slot per (nuclide,
channel) in scheme order whose photon yield lies strictly between 0 and 1.
Every slot holds ``parts`` words (parts = ceil(total population / part_max)),
read even when its count is 0. Consumption is therefore a fixed
``words_per_interval`` per interval, so interval k starts at word k times that,
which is what checkpoint and resume rely on.
"""
from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from .binomial import PART_MAX, binomial_draw, binomial_inverse
from .decay import MODE_SPECS, DecayScheme, Mode
from .manifest import environment
from .paths import IntervalLaw, interval_law, path_categories, tree_probabilities
from .streams import WordStream

LAYOUT = "nuclear-aggregate-layout-v1"
CHECKPOINT_SCHEMA = "rnglab-nuclear-aggregate-checkpoint-v1"
MANIFEST_SCHEMA = "rnglab-nuclear-aggregate-manifest-v1"
INT64_MAX = 2 ** 63 - 1


def _check_initial(scheme: DecayScheme, initial: Mapping[str, int]) -> np.ndarray:
    if not isinstance(initial, Mapping):
        raise TypeError("initial populations must be a mapping name -> count")
    unknown = set(initial) - set(scheme.index)
    if unknown:
        raise ValueError(f"unknown nuclides in initial populations: {sorted(unknown)}")
    counts = []
    for nuclide in scheme.nuclides:
        count = initial.get(nuclide.name, 0)
        if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 0:
            raise ValueError(f"initial count of {nuclide.name} must be a nonnegative integer")
        counts.append(int(count))
    if sum(counts) > INT64_MAX:
        raise ValueError("total population exceeds 2**63 - 1; not supported in exact mode")
    return np.array(counts, dtype=np.int64)


def _check_times(times: Sequence[float]) -> np.ndarray:
    grid = np.array(times, dtype=np.float64)  # a copy: the caller's array is never aliased
    if grid.ndim != 1 or grid.size < 2 or not np.all(np.isfinite(grid)) or grid[0] < 0 or np.any(np.diff(grid) <= 0):
        raise ValueError("times must be a finite, strictly increasing grid of at least two points from t >= 0")
    grid.setflags(write=False)
    return grid


@dataclass(frozen=True)
class SourcePlan:
    """Sampling plan of one start nuclide for one interval length."""

    start: int
    law: IntervalLaw
    p_small: np.ndarray        # (K - 1,) probability actually inverted, <= 1/2
    flip: np.ndarray           # (K - 1,) True when the conditional probability exceeds 1/2
    first_category: int        # global index of the source's first category


class AggregatePlan:
    """Everything fixed before sampling: categories, interval laws, slots and word layout."""

    def __init__(self, scheme: DecayScheme, times: Sequence[float], total: int, *, part_max: int = PART_MAX):
        if not isinstance(scheme, DecayScheme):
            raise TypeError("scheme must be a DecayScheme")
        if isinstance(part_max, bool) or not isinstance(part_max, int) or not 1 <= part_max <= PART_MAX:
            raise ValueError(f"part_max must be an integer in [1, {PART_MAX}]")
        self.scheme = scheme
        self.times = _check_times(times)
        self.total = int(total)
        self.part_max = part_max
        self.parts = max(1, -(-self.total // part_max))
        self.categories = [path_categories(scheme, i) for i in range(len(scheme))]
        self.first_category = np.cumsum([0] + [len(c) for c in self.categories[:-1]]).astype(np.int64)
        self.n_categories = int(sum(len(c) for c in self.categories))
        width = max(scheme.max_channels, 1)
        self.n_channels = width
        self.category_end = np.zeros(self.n_categories, dtype=np.int64)
        self.category_steps = np.zeros((self.n_categories, len(scheme) * width), dtype=np.int64)
        for i, cats in enumerate(self.categories):
            for k, category in enumerate(cats):
                g = self.first_category[i] + k
                self.category_end[g] = category.end
                for parent, channel in category.steps:
                    self.category_steps[g, parent * width + channel] += 1
        self.sources = [i for i, cats in enumerate(self.categories) if len(cats) >= 2]
        self.yield_channels = [(j, k, c.photon_yield) for j, n in enumerate(scheme.nuclides)
                               for k, c in enumerate(n.channels) if 0.0 < c.photon_yield < 1.0]
        self.full_yield_channels = [(j, k) for j, n in enumerate(scheme.nuclides)
                                    for k, c in enumerate(n.channels) if c.photon_yield == 1.0]
        self.neutron_channels = [(j, k) for j, n in enumerate(scheme.nuclides)
                                 for k, c in enumerate(n.channels) if c.mode is Mode.BETA_DELAYED_NEUTRON]
        neutron_columns = [j * width + k for j, k in self.neutron_channels]
        self.max_neutrons_per_nucleus = (int(self.category_steps[:, neutron_columns].sum(axis=1).max())
                                         if neutron_columns else 0)
        if self.max_neutrons_per_nucleus * self.total > INT64_MAX:
            raise ValueError("the neutron count of one interval could exceed 2**63 - 1 for this scheme and "
                             "population; not supported")
        self.slot_offset = {}
        slot = 0
        for i in self.sources:
            self.slot_offset[i] = slot
            slot += len(self.categories[i]) - 1
        self.photon_slot_offset = slot
        self.slots = slot + len(self.yield_channels)
        self.words_per_interval = self.slots * self.parts
        self.deltas = np.diff(self.times)
        for array in (self.first_category, self.category_end, self.category_steps, self.deltas):
            array.setflags(write=False)
        self.interval_plans: dict[float, dict[int, SourcePlan]] = {}
        # The path-tree matrix exponential is a cross-check only; where it is not finite (rates near the float
        # limits) the comparison is counted as unavailable, never as agreement.
        self.diagnostics: dict[str, Any] = {"max_abs_decimal_minus_tree": 0.0, "tree_comparisons": 0,
                                            "tree_comparisons_unavailable": 0, "max_digits": 0}
        for delta in sorted(set(self.deltas.tolist())):
            plans = {}
            for i in self.sources:
                law = interval_law(scheme, i, delta)
                p, r = law.conditional()
                flip = p > 0.5
                p_small = np.where(flip, r, p)
                p_small.setflags(write=False)
                flip.setflags(write=False)
                plans[i] = SourcePlan(i, law, p_small, flip, int(self.first_category[i]))
                tree = tree_probabilities(scheme, i, delta)
                if np.all(np.isfinite(tree)):
                    self.diagnostics["tree_comparisons"] += 1
                    self.diagnostics["max_abs_decimal_minus_tree"] = max(
                        self.diagnostics["max_abs_decimal_minus_tree"], float(np.max(np.abs(law.as_float() - tree))))
                else:
                    self.diagnostics["tree_comparisons_unavailable"] += 1
                self.diagnostics["max_digits"] = max(self.diagnostics["max_digits"], law.digits_used)
            self.interval_plans[delta] = plans

    def layout(self) -> dict[str, Any]:
        return {"layout": LAYOUT, "parts": self.parts, "part_max": self.part_max, "slots": self.slots,
                "words_per_interval": self.words_per_interval, "categories": self.n_categories,
                "multinomial_slots": self.photon_slot_offset, "photon_slots": len(self.yield_channels)}


def _split_sizes(n: np.ndarray, parts: int, part_max: int) -> np.ndarray:
    offsets = np.arange(parts, dtype=np.int64) * part_max
    return np.clip(n[:, None] - offsets[None, :], 0, part_max)


def _vector_binomials(n: np.ndarray, p_small: np.ndarray, flip: np.ndarray, uniforms: np.ndarray,
                      parts: int, part_max: int) -> np.ndarray:
    """Binomial(n, p) for many slots at once; uniforms has shape (slots, parts)."""
    sizes = _split_sizes(n, parts, part_max)
    draws = binomial_inverse(sizes.ravel(), np.repeat(p_small, parts), uniforms.ravel()).reshape(sizes.shape)
    total = draws.sum(axis=1)
    return np.where(flip, n - total, total)


def _interval_record(plan: AggregatePlan, populations: np.ndarray, counts: np.ndarray,
                     photons: np.ndarray) -> dict[str, np.ndarray]:
    """One interval's outputs, as read-only arrays (a sink cannot alter the run through them)."""
    new_populations = np.zeros(len(plan.scheme), dtype=np.int64)
    np.add.at(new_populations, plan.category_end, counts)
    decays = (counts @ plan.category_steps).reshape(len(plan.scheme), plan.n_channels)
    neutrons = np.int64(sum(int(decays[j, k]) for j, k in plan.neutron_channels))
    record = {"populations": new_populations, "path_counts": np.array(counts, dtype=np.int64),
              "decays": decays, "photons": np.array(photons, dtype=np.int64), "neutrons": neutrons}
    for key in ("populations", "path_counts", "decays", "photons"):
        record[key].setflags(write=False)
    return record


def _photons_from(plan: AggregatePlan, decays: np.ndarray, sampled: np.ndarray) -> np.ndarray:
    photons = np.zeros_like(decays)
    for (j, k) in plan.full_yield_channels:
        photons[j, k] = decays[j, k]
    for slot, (j, k, _) in enumerate(plan.yield_channels):
        photons[j, k] = sampled[slot]
    return photons


def _chain_digest(previous: bytes, interval: int, record: Mapping[str, np.ndarray]) -> bytes:
    digest = hashlib.sha256(previous)
    digest.update(f"interval:{interval}".encode())
    for name in ("populations", "path_counts", "decays", "photons"):
        digest.update(name.encode() + np.ascontiguousarray(np.asarray(record[name], dtype="<i8")).tobytes())
    digest.update(np.asarray([record["neutrons"]], dtype="<i8").tobytes())
    return digest.digest()


def _initial_digest(plan: AggregatePlan, populations: np.ndarray) -> bytes:
    header = json.dumps({"scheme_sha256": plan.scheme.sha256(), "times": [repr(t) for t in plan.times.tolist()],
                         "layout": plan.layout(), "initial": populations.tolist()}, sort_keys=True)
    return hashlib.sha256(header.encode()).digest()


@dataclass(frozen=True, eq=False)
class AggregateRun:
    """Result of an aggregated run (or of the part of it after a resume).

    ``populations`` has one row per grid time from ``first_interval`` on;
    ``path_counts``, ``decays``, ``photons`` and ``neutrons`` one row per
    interval. Arrays are empty when the run was made with keep_outputs=False.
    """

    plan: AggregatePlan
    initial: dict[str, int]
    first_interval: int
    first_word: int
    engine: str
    populations: np.ndarray
    path_counts: np.ndarray
    decays: np.ndarray
    photons: np.ndarray
    neutrons: np.ndarray
    final_populations: np.ndarray
    digest: str
    words_read: int
    layout: str = field(default=LAYOUT)

    def summary(self) -> dict[str, Any]:
        names = [n.name for n in self.plan.scheme.nuclides]
        return {"engine": self.engine, "layout": self.plan.layout(), "first_interval": self.first_interval,
                "intervals": int(self.plan.deltas.size - self.first_interval), "first_word": self.first_word,
                "words_read": self.words_read,
                "final_populations": dict(zip(names, self.final_populations.tolist())), "digest": self.digest,
                "plan_diagnostics": dict(self.plan.diagnostics)}


def stream_identity(stream: WordStream) -> dict[str, Any]:
    """Which generator a stream reads and where it stands, in words from the generator's origin.

    The position is the source's declared start_word (0 if it declares none) plus the words it
    has emitted. A checkpoint stores this, and a resume requires the same identity, so a stream
    with another seed or at another position is refused.
    """
    manifest = stream.source.manifest
    metadata = manifest.metadata if isinstance(manifest.metadata, Mapping) else {}
    start = metadata.get("start_word", 0)
    if isinstance(start, bool) or not isinstance(start, int) or start < 0:
        start = 0
    consumed = int(stream.source.consumed_bytes)
    return {"family": manifest.family, "generator": metadata.get("generator"), "seed": manifest.seed,
            "word_position": start + consumed // 8, "byte_aligned": consumed % 8 == 0}


def _checkpoint_checksum(body: Mapping[str, Any]) -> str:
    try:
        text = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (TypeError, ValueError) as error:
        raise ValueError("checkpoint is not canonical JSON data") from error
    return hashlib.sha256(text.encode()).hexdigest()


def make_checkpoint(plan: AggregatePlan, initial: Mapping[str, int], next_interval: int, populations: np.ndarray,
                    running: bytes, words_consumed: int, stream_state: Mapping[str, Any]) -> dict[str, Any]:
    """A JSON-serialisable state from which the run continues bit for bit.

    It shares no object with the engine. ``checkpoint_sha256`` is the SHA-256 of the canonical JSON
    of every other field: it detects corruption and inconsistent edits (it is an integrity check,
    not an authentication: whoever can rewrite the checksum can forge a checkpoint).
    """
    body = {"schema": CHECKPOINT_SCHEMA, "scheme_sha256": plan.scheme.sha256(),
            "times": [repr(t) for t in plan.times.tolist()], "layout": plan.layout(),
            "initial": {k: int(v) for k, v in initial.items()}, "next_interval": int(next_interval),
            "populations": [int(x) for x in populations], "running_sha256": running.hex(),
            "words_consumed": int(words_consumed), "stream": dict(stream_state)}
    return {**body, "checkpoint_sha256": _checkpoint_checksum(body)}


def _strict_int(value: object, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"checkpoint {what} must be an integer")
    return value


def _resume_state(plan: AggregatePlan, initial: Mapping[str, int], checkpoint: Mapping[str, Any],
                  stream: WordStream):
    """Validate a checkpoint and return the state it holds.

    Refused: wrong keys; a checksum that does not match the other fields; another scheme, grid,
    layout or initial population; mistyped or fractional fields (nothing is coerced); a word count
    that contradicts the layout; populations that are negative or do not conserve the total; a
    malformed running digest; a stream whose generator, seed or word position differs.
    """
    if not isinstance(checkpoint, Mapping):
        raise ValueError("checkpoint must be a mapping")
    expected = make_checkpoint(plan, initial, 0, np.zeros(len(plan.scheme), dtype=np.int64), b"", 0, {})
    if set(checkpoint) != set(expected):
        raise ValueError("checkpoint keys do not match the schema")
    body = {key: value for key, value in checkpoint.items() if key != "checkpoint_sha256"}
    if checkpoint["checkpoint_sha256"] != _checkpoint_checksum(body):
        raise ValueError("checkpoint checksum does not match its contents (corrupted or edited)")
    for key in ("schema", "scheme_sha256", "times", "layout", "initial"):
        if checkpoint[key] != expected[key]:
            raise ValueError(f"checkpoint does not match this run: {key}")
    next_interval = _strict_int(checkpoint["next_interval"], "next_interval")
    if not 0 <= next_interval <= plan.deltas.size:
        raise ValueError("checkpoint interval out of range")
    if _strict_int(checkpoint["words_consumed"], "words_consumed") != next_interval * plan.words_per_interval:
        raise ValueError("checkpoint word count does not match the layout")
    raw = checkpoint["populations"]
    if not isinstance(raw, list) or len(raw) != len(plan.scheme):
        raise ValueError("checkpoint populations must be a list with one count per nuclide")
    counts = [_strict_int(x, "population") for x in raw]
    if any(c < 0 for c in counts) or sum(counts) != plan.total:
        raise ValueError("checkpoint populations are negative or do not conserve the total")
    running = checkpoint["running_sha256"]
    if not isinstance(running, str) or len(running) != 64 or any(c not in "0123456789abcdef" for c in running):
        raise ValueError("checkpoint running_sha256 must be 64 lowercase hex digits")
    here = stream_identity(stream)
    if checkpoint["stream"] != here or not here["byte_aligned"]:
        raise ValueError("the stream is not the checkpoint's generator at the checkpoint's word position "
                         f"(checkpoint {checkpoint['stream']!r}, stream {here!r})")
    return next_interval, np.array(counts, dtype=np.int64), bytes.fromhex(running)


def simulate_aggregate(scheme: DecayScheme, initial: Mapping[str, int], times: Sequence[float], stream: WordStream,
                       *, part_max: int | None = None, resume_from: Mapping[str, Any] | None = None,
                       checkpoint_every: int = 0, on_checkpoint: Callable[[dict], None] | None = None,
                       sink: Callable[[int, dict], None] | None = None, keep_outputs: bool = True,
                       plan: AggregatePlan | None = None) -> AggregateRun:
    """Aggregated decay over the grid ``times`` (vectorised engine).

    With ``resume_from`` the stream must already be positioned at the
    checkpoint's word (for NumpyPCG64Source: start_word = first word of the run
    + words_consumed). ``sink(interval, record)`` receives every interval's
    record; ``keep_outputs=False`` keeps only the final state (bounded memory).
    """
    if not isinstance(stream, WordStream):
        raise TypeError("stream must be a WordStream")
    populations = _check_initial(scheme, initial)
    if plan is None:
        plan = AggregatePlan(scheme, times, int(populations.sum()), part_max=PART_MAX if part_max is None else part_max)
    elif not isinstance(plan, AggregatePlan):
        raise TypeError("plan must be an AggregatePlan")
    if (plan.scheme is not scheme or plan.total != int(populations.sum())
            or not np.array_equal(plan.times, _check_times(times)) or part_max not in (None, plan.part_max)):
        raise ValueError("plan does not match the scheme, population, times and part_max")
    if isinstance(checkpoint_every, bool) or not isinstance(checkpoint_every, int) or checkpoint_every < 0:
        raise ValueError("checkpoint_every must be a nonnegative integer")
    if bool(checkpoint_every) != (on_checkpoint is not None):
        raise ValueError("checkpoint_every and on_checkpoint must be given together")
    initial = {n.name: int(c) for n, c in zip(scheme.nuclides, populations)}
    running = _initial_digest(plan, populations)
    first_interval = 0
    if resume_from is not None:
        first_interval, populations, running = _resume_state(plan, initial, resume_from, stream)
    first_word = stream.words_read
    kept: dict[str, list] = {"populations": [populations.copy()], "path_counts": [], "decays": [], "photons": [],
                             "neutrons": []}
    for interval in range(first_interval, plan.deltas.size):
        sources = plan.interval_plans[float(plan.deltas[interval])]
        u = stream.uniforms(plan.words_per_interval).reshape(plan.slots, plan.parts)
        counts = np.zeros(plan.n_categories, dtype=np.int64)
        for i, n in enumerate(plan.categories):
            if len(n) == 1:
                counts[plan.first_category[i]] = populations[i]
        remaining = {i: populations[i] for i in plan.sources}
        stages = max((len(plan.categories[i]) - 1 for i in plan.sources), default=0)
        for stage in range(stages):
            active = [i for i in plan.sources if len(plan.categories[i]) - 1 > stage]
            n = np.array([remaining[i] for i in active], dtype=np.int64)
            p_small = np.array([sources[i].p_small[stage] for i in active])
            flip = np.array([sources[i].flip[stage] for i in active])
            slots = [plan.slot_offset[i] + stage for i in active]
            x = _vector_binomials(n, p_small, flip, u[slots], plan.parts, plan.part_max)
            for i, value in zip(active, x):
                counts[plan.first_category[i] + stage] = value
                remaining[i] -= value
        for i in plan.sources:
            counts[plan.first_category[i] + len(plan.categories[i]) - 1] = remaining[i]
        decays = (counts @ plan.category_steps).reshape(len(scheme), plan.n_channels)
        if plan.yield_channels:
            n = np.array([decays[j, k] for j, k, _ in plan.yield_channels], dtype=np.int64)
            y = np.array([yield_ for _, _, yield_ in plan.yield_channels])
            flip = y > 0.5
            sampled = _vector_binomials(n, np.where(flip, 1.0 - y, y), flip, u[plan.photon_slot_offset:],
                                        plan.parts, plan.part_max)
        else:
            sampled = np.zeros(0, dtype=np.int64)
        record = _interval_record(plan, populations, counts, _photons_from(plan, decays, sampled))
        running = _chain_digest(running, interval, record)
        populations = record["populations"]
        if sink is not None:  # detached copies: nothing a sink does to them can reach the engine or its outputs
            sink(interval, {key: np.array(value, copy=True) for key, value in record.items()})
        if keep_outputs:
            for key in kept:
                kept[key].append(record[key])
        if checkpoint_every and (interval + 1) % checkpoint_every == 0:
            on_checkpoint(make_checkpoint(plan, initial, interval + 1, populations, running,
                                          (interval + 1) * plan.words_per_interval, stream_identity(stream)))
    return _finish(plan, initial, first_interval, first_word, "vectorised", kept, populations, running, stream,
                   keep_outputs)


def _finish(plan, initial, first_interval, first_word, engine, kept, populations, running, stream, keep_outputs):
    width = (len(plan.scheme), plan.n_channels)
    if keep_outputs:
        arrays = {"populations": np.array(kept["populations"], dtype=np.int64),
                  "path_counts": np.array(kept["path_counts"], dtype=np.int64).reshape(-1, plan.n_categories),
                  "decays": np.array(kept["decays"], dtype=np.int64).reshape(-1, *width),
                  "photons": np.array(kept["photons"], dtype=np.int64).reshape(-1, *width),
                  "neutrons": np.array(kept["neutrons"], dtype=np.int64)}
    else:
        arrays = {"populations": np.zeros((0, len(plan.scheme)), dtype=np.int64),
                  "path_counts": np.zeros((0, plan.n_categories), dtype=np.int64),
                  "decays": np.zeros((0, *width), dtype=np.int64), "photons": np.zeros((0, *width), dtype=np.int64),
                  "neutrons": np.zeros(0, dtype=np.int64)}
    for value in arrays.values():
        value.setflags(write=False)
    final = populations.copy()
    final.setflags(write=False)
    return AggregateRun(plan, dict(initial), first_interval, first_word, engine, arrays["populations"],
                        arrays["path_counts"], arrays["decays"], arrays["photons"], arrays["neutrons"], final,
                        running.hex(), stream.words_read - first_word)


def simulate_aggregate_reference(scheme: DecayScheme, initial: Mapping[str, int], times: Sequence[float],
                                 stream: WordStream, *, part_max: int = PART_MAX) -> AggregateRun:
    """Pure-Python reference engine: one slot at a time, same layout; for cross-checks only.

    It shares the plan (laws and conditional probabilities) and the inversion
    routine with ``simulate_aggregate``; the bookkeeping (slot order, splitting,
    flips, photon and neutron accounting, digests) is written out again here.
    """
    populations = _check_initial(scheme, initial)
    plan = AggregatePlan(scheme, times, int(populations.sum()), part_max=part_max)
    initial = {n.name: int(c) for n, c in zip(scheme.nuclides, populations)}
    running = _initial_digest(plan, populations)
    first_word = stream.words_read
    kept: dict[str, list] = {"populations": [populations.copy()], "path_counts": [], "decays": [], "photons": [],
                             "neutrons": []}
    width = max(scheme.max_channels, 1)
    for interval in range(plan.deltas.size):
        sources = plan.interval_plans[float(plan.deltas[interval])]
        counts = [0] * plan.n_categories
        for i, cats in enumerate(plan.categories):
            start = int(plan.first_category[i])
            if len(cats) == 1:
                counts[start] = int(populations[i])
                continue
            remaining = int(populations[i])
            for k in range(len(cats) - 1):
                words = stream.uniforms(plan.parts)
                value = binomial_draw(remaining, float(sources[i].p_small[k]), bool(sources[i].flip[k]), words,
                                      part_max=plan.part_max)
                counts[start + k] = value
                remaining -= value
            counts[start + len(cats) - 1] = remaining
        decays = [[0] * width for _ in scheme.nuclides]
        new_populations = [0] * len(scheme)
        for i, cats in enumerate(plan.categories):
            for k, category in enumerate(cats):
                c = counts[int(plan.first_category[i]) + k]
                new_populations[category.end] += c
                for parent, channel in category.steps:
                    decays[parent][channel] += c
        photons = [[0] * width for _ in scheme.nuclides]
        for j, nuclide in enumerate(scheme.nuclides):
            for k, channel in enumerate(nuclide.channels):
                if channel.photon_yield == 1.0:
                    photons[j][k] = decays[j][k]
        for j, k, y in plan.yield_channels:
            words = stream.uniforms(plan.parts)
            flip = y > 0.5
            photons[j][k] = binomial_draw(decays[j][k], 1.0 - y if flip else y, flip, words, part_max=plan.part_max)
        neutrons = sum(decays[j][k] for j, n in enumerate(scheme.nuclides) for k, c in enumerate(n.channels)
                       if c.mode is Mode.BETA_DELAYED_NEUTRON)
        record = {"populations": np.array(new_populations, dtype=np.int64),
                  "path_counts": np.array(counts, dtype=np.int64), "decays": np.array(decays, dtype=np.int64),
                  "photons": np.array(photons, dtype=np.int64), "neutrons": np.int64(neutrons)}
        running = _chain_digest(running, interval, record)
        populations = record["populations"]
        for key in kept:
            kept[key].append(record[key])
    return _finish(plan, initial, 0, first_word, "reference", kept, populations, running, stream, True)


def aggregate_manifest(*, experiment: str, run_id: str, run: AggregateRun, stream: WordStream,
                       code: Mapping[str, Any], extra: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Manifest of one aggregated run: configuration, stream provenance, output digest."""
    return {"schema": MANIFEST_SCHEMA, "experiment": experiment, "run_id": run_id, "code": dict(code),
            "environment": environment(),
            "rng": {"role": "physics", **run.plan.layout(), "first_word_of_run": run.first_word,
                    "words_read_by_run": run.words_read, **stream.provenance()},
            "configuration": {"scheme": run.plan.scheme.to_dict(), "scheme_sha256": run.plan.scheme.sha256(),
                              "initial": dict(run.initial), "times": [repr(t) for t in run.plan.times.tolist()],
                              "engine": run.engine},
            "outputs": run.summary(), "extra": dict(extra or {})}
