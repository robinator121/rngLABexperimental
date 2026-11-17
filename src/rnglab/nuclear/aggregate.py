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
- No state is shared with callers. A plan is immutable once built: its
  attributes cannot be rebound, its containers are tuples and read-only
  mappings, and its arrays are read-only views of immutable bytes (their
  WRITEABLE flag cannot be set again). It keeps its own snapshot of the
  scheme's nuclides and the scheme's SHA-256; a run refuses a scheme that no
  longer matches that digest and reads nothing else of the caller's while it
  runs. Results are immutable in the same way; a sink receives detached,
  writable copies of every interval's outputs; a checkpoint is built from
  fresh Python objects. Conservation of the total population is checked after
  every interval.
- A run owns its stream's source. After every interval and every callback the
  engine checks that the source has emitted exactly the layout's words, all
  through this stream; a read by anyone else (another WordStream on the same
  source, a sink drawing numbers) stops the run with RuntimeError. Reads of a
  source that this source wraps cannot be seen and are the caller's to avoid.
  The run records its stream's provenance when it ends; manifests are built
  from that record, never from a live stream.
- Checkpoint and resume need a replayable generator: exactly one of
  NumpyPCG64Source, rnglab14's PCG64Source or Xoshiro256PlusPlusSource, whose
  output is fixed by the parameters a checkpoint records (generator, seed and,
  for PCG64Source, the stream id) and a word position. Any other source (fault
  adapters, physical or OS entropy, conditioned sources) is refused before any
  random number is read, because a checkpoint could not identify its sequence.
- A checkpoint carries the SHA-256 of its own contents, of the scheme and of
  the time grid, and the stream's replay identity. A resume refuses a
  checkpoint whose checksum, types, layout, word count or conservation is
  wrong, and a stream with another generator, seed, stream id or position.
  The checksum is an integrity check, not an authentication.

Memory and cost
---------------
With keep_outputs=False the outputs a run retains do not grow with the number
of intervals I. The plan does: it holds the grid (O(I)) and one interval law
per distinct interval length and start nuclide, O(L * C) decimals for L
distinct lengths and C categories in all (a uniform grid has L = 1). Building
it also runs the independent tree check, a dense K x K matrix exponential per
start nuclide with K categories (8 K**2 bytes, O(K**3) time), for every
distinct length. K grows with branching, exponentially with the depth of a
chain whose nuclides have two channels; a plan counts the categories exactly
before enumerating them and refuses more than ``max_categories`` (default
1024) for one start nuclide. A checkpoint holds the grid's SHA-256, not the
grid, so each costs O(nuclides).

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
from types import MappingProxyType
from typing import Any, Callable, Mapping, Sequence

import numpy as np

from rnglab14.sources.reference import PCG64Source, Xoshiro256PlusPlusSource

from .binomial import PART_MAX, binomial_draw, binomial_inverse
from .decay import MODE_SPECS, DecayScheme, Mode
from .manifest import environment
from .paths import IntervalLaw, interval_law, path_categories, tree_probabilities
from .streams import NumpyPCG64Source, WordStream

LAYOUT = "nuclear-aggregate-layout-v1"
CHECKPOINT_SCHEMA = "rnglab-nuclear-aggregate-checkpoint-v2"
MANIFEST_SCHEMA = "rnglab-nuclear-aggregate-manifest-v2"
INT64_MAX = 2 ** 63 - 1
MAX_CATEGORIES = 1024


def _frozen(values: Any, dtype: Any = None) -> np.ndarray:
    """A read-only copy backed by an immutable bytes object, so its WRITEABLE flag cannot be set again."""
    array = np.ascontiguousarray(values if dtype is None else np.asarray(values, dtype=dtype))
    return np.frombuffer(array.tobytes(), dtype=array.dtype).reshape(array.shape)


def category_counts(scheme: DecayScheme) -> list[int]:
    """Number of path categories of every start nuclide, counted exactly without enumerating them."""
    memo: dict[int, int] = {}

    def count(i: int) -> int:
        if i not in memo:
            memo[i] = 1 + sum(count(scheme.index[c.daughter]) for c in scheme.nuclides[i].channels)
        return memo[i]

    return [count(i) for i in range(len(scheme))]


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
    return _frozen(grid)


def _grid_sha256(grid: np.ndarray) -> str:
    return hashlib.sha256(json.dumps([repr(t) for t in grid.tolist()]).encode()).hexdigest()


@dataclass(frozen=True)
class SourcePlan:
    """Sampling plan of one start nuclide for one interval length."""

    start: int
    law: IntervalLaw
    p_small: np.ndarray        # (K - 1,) probability actually inverted, <= 1/2
    flip: np.ndarray           # (K - 1,) True when the conditional probability exceeds 1/2
    first_category: int        # global index of the source's first category


class AggregatePlan:
    """Everything fixed before sampling: categories, interval laws, slots and word layout.

    Immutable once built (see the module docstring): attributes cannot be rebound or deleted,
    containers are tuples and read-only mappings, arrays are read-only views of immutable bytes.
    ``scheme`` is the caller's object, kept to check identity; the engines use ``nuclides`` (the
    scheme's tuple of frozen nuclides, taken when the plan is built) and ``scheme_sha256``.
    """

    _sealed = False

    def __setattr__(self, name: str, value: Any) -> None:
        if self._sealed:
            raise AttributeError(f"AggregatePlan is immutable; cannot set {name!r}")
        object.__setattr__(self, name, value)

    def __delattr__(self, name: str) -> None:
        raise AttributeError(f"AggregatePlan is immutable; cannot delete {name!r}")

    def __init__(self, scheme: DecayScheme, times: Sequence[float], total: int, *, part_max: int = PART_MAX,
                 max_categories: int = MAX_CATEGORIES):
        if not isinstance(scheme, DecayScheme):
            raise TypeError("scheme must be a DecayScheme")
        if isinstance(part_max, bool) or not isinstance(part_max, int) or not 1 <= part_max <= PART_MAX:
            raise ValueError(f"part_max must be an integer in [1, {PART_MAX}]")
        if isinstance(total, bool) or not isinstance(total, (int, np.integer)) or not 0 <= total <= INT64_MAX:
            raise ValueError("total must be an integer in [0, 2**63 - 1]")
        if isinstance(max_categories, bool) or not isinstance(max_categories, int) or max_categories < 1:
            raise ValueError("max_categories must be a positive integer")
        counted = category_counts(scheme)
        if max(counted) > max_categories:
            raise ValueError(f"nuclide {scheme.nuclides[counted.index(max(counted))].name} has {max(counted)} path "
                             f"categories, more than max_categories = {max_categories}; the plan would need "
                             f"O(K**2) memory and O(K**3) time per start nuclide")
        self.scheme = scheme
        self.scheme_sha256 = scheme.sha256()
        self.nuclides = tuple(scheme.nuclides)
        self.n_nuclides = len(self.nuclides)
        self.times = _check_times(times)
        self.grid_sha256 = _grid_sha256(self.times)
        self.total = int(total)
        self.part_max = part_max
        self.parts = max(1, -(-self.total // part_max))
        self.categories = tuple(tuple(path_categories(scheme, i)) for i in range(len(scheme)))
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
        self.sources = tuple(i for i, cats in enumerate(self.categories) if len(cats) >= 2)
        self.yield_channels = tuple((j, k, c.photon_yield) for j, n in enumerate(self.nuclides)
                                    for k, c in enumerate(n.channels) if 0.0 < c.photon_yield < 1.0)
        self.full_yield_channels = tuple((j, k) for j, n in enumerate(self.nuclides)
                                         for k, c in enumerate(n.channels) if c.photon_yield == 1.0)
        self.neutron_channels = tuple((j, k) for j, n in enumerate(self.nuclides)
                                      for k, c in enumerate(n.channels) if c.mode is Mode.BETA_DELAYED_NEUTRON)
        neutron_columns = [j * width + k for j, k in self.neutron_channels]
        self.max_neutrons_per_nucleus = (int(self.category_steps[:, neutron_columns].sum(axis=1).max())
                                         if neutron_columns else 0)
        if self.max_neutrons_per_nucleus * self.total > INT64_MAX:
            raise ValueError("the neutron count of one interval could exceed 2**63 - 1 for this scheme and "
                             "population; not supported")
        slot_offset = {}
        slot = 0
        for i in self.sources:
            slot_offset[i] = slot
            slot += len(self.categories[i]) - 1
        self.slot_offset = MappingProxyType(slot_offset)
        self.photon_slot_offset = slot
        self.slots = slot + len(self.yield_channels)
        self.words_per_interval = self.slots * self.parts
        self.deltas = _frozen(np.diff(self.times))
        self.first_category = _frozen(self.first_category)
        self.category_end = _frozen(self.category_end)
        self.category_steps = _frozen(self.category_steps)
        interval_plans: dict[float, Mapping[int, SourcePlan]] = {}
        # The path-tree matrix exponential is a cross-check only; where it is not finite (rates near the float
        # limits) the comparison is counted as unavailable, never as agreement.
        diagnostics: dict[str, Any] = {"max_abs_decimal_minus_tree": 0.0, "tree_comparisons": 0,
                                       "tree_comparisons_unavailable": 0, "max_digits": 0}
        for delta in sorted(set(self.deltas.tolist())):
            plans = {}
            for i in self.sources:
                law = interval_law(scheme, i, delta)
                p, r = law.conditional()
                flip = p > 0.5
                plans[i] = SourcePlan(i, law, _frozen(np.where(flip, r, p)), _frozen(flip),
                                      int(self.first_category[i]))
                tree = tree_probabilities(scheme, i, delta)
                if np.all(np.isfinite(tree)):
                    diagnostics["tree_comparisons"] += 1
                    diagnostics["max_abs_decimal_minus_tree"] = max(
                        diagnostics["max_abs_decimal_minus_tree"], float(np.max(np.abs(law.as_float() - tree))))
                else:
                    diagnostics["tree_comparisons_unavailable"] += 1
                diagnostics["max_digits"] = max(diagnostics["max_digits"], law.digits_used)
            interval_plans[delta] = MappingProxyType(plans)
        self.interval_plans = MappingProxyType(interval_plans)
        self.diagnostics = MappingProxyType(diagnostics)
        self._sealed = True

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
    new_populations = np.zeros(plan.n_nuclides, dtype=np.int64)
    np.add.at(new_populations, plan.category_end, counts)
    decays = (counts @ plan.category_steps).reshape(plan.n_nuclides, plan.n_channels)
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
    header = json.dumps({"scheme_sha256": plan.scheme_sha256, "times": [repr(t) for t in plan.times.tolist()],
                         "layout": plan.layout(), "initial": populations.tolist()}, sort_keys=True)
    return hashlib.sha256(header.encode()).digest()


@dataclass(frozen=True, eq=False)
class AggregateRun:
    """Result of an aggregated run (or of the part of it after a resume).

    ``populations`` has one row per grid time from ``first_interval`` on;
    ``path_counts``, ``decays``, ``photons`` and ``neutrons`` one row per
    interval. Arrays are empty when the run was made with keep_outputs=False.
    Immutable: ``initial`` is a read-only mapping, arrays are read-only views of
    immutable bytes, and ``provenance`` is the canonical JSON of the stream's
    provenance when the run ended (``stream_provenance()`` returns a fresh copy).
    """

    plan: AggregatePlan
    initial: Mapping[str, int]
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
    provenance: str
    layout: str = field(default=LAYOUT)

    def stream_provenance(self) -> dict[str, Any]:
        return json.loads(self.provenance)

    def summary(self) -> dict[str, Any]:
        names = [n.name for n in self.plan.nuclides]
        return {"engine": self.engine, "layout": self.plan.layout(), "first_interval": self.first_interval,
                "intervals": int(self.plan.deltas.size - self.first_interval), "first_word": self.first_word,
                "words_read": self.words_read,
                "final_populations": dict(zip(names, self.final_populations.tolist())), "digest": self.digest,
                "plan_diagnostics": dict(self.plan.diagnostics)}


REPLAYABLE = ("rnglab.nuclear.streams.NumpyPCG64Source", "rnglab14.sources.reference.PCG64Source",
              "rnglab14.sources.reference.Xoshiro256PlusPlusSource")


def replay_identity(source: Any) -> dict[str, Any] | None:
    """The parameters that fix a source's output sequence, and its origin, or None if it is not replayable.

    Only the exact classes in REPLAYABLE qualify (a subclass or wrapper could change the sequence):
    NumpyPCG64Source (seed; its declared start_word is where it begins), rnglab14's PCG64Source
    (seed and stream id) and Xoshiro256PlusPlusSource (seed). Returns {"generator", "parameters",
    "start_word"}; a seed or stream id that is not a nonnegative integer makes the source unidentifiable.
    """
    kind = type(source)
    metadata = source.manifest.metadata if hasattr(source, "manifest") else None
    seed = getattr(getattr(source, "manifest", None), "seed", None)
    if not isinstance(metadata, Mapping) or isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        return None
    if kind is NumpyPCG64Source:
        start = metadata.get("start_word")
        if isinstance(start, bool) or not isinstance(start, int) or start < 0:
            return None
        return {"generator": REPLAYABLE[0], "parameters": {"seed": seed}, "start_word": start}
    if kind is PCG64Source:
        stream_id = metadata.get("stream_id")
        if isinstance(stream_id, bool) or not isinstance(stream_id, int) or stream_id < 0:
            return None
        return {"generator": REPLAYABLE[1], "parameters": {"seed": seed, "stream_id": stream_id}, "start_word": 0}
    if kind is Xoshiro256PlusPlusSource:
        return {"generator": REPLAYABLE[2], "parameters": {"seed": seed}, "start_word": 0}
    return None


def stream_identity(stream: WordStream) -> dict[str, Any]:
    """Which sequence a stream reads and where it stands, in words from the generator's origin.

    The position is the source's start word plus the words it has emitted. A checkpoint stores
    this, and a resume requires the same identity, so another generator, seed, stream id or
    position is refused. Raises ValueError for a source that is not replayable (see replay_identity).
    """
    identity = replay_identity(stream.source)
    if identity is None:
        raise ValueError(f"{type(stream.source).__module__}.{type(stream.source).__qualname__} is not a replayable "
                         f"generator; checkpoint and resume support exactly {', '.join(REPLAYABLE)}")
    consumed = int(stream.source.consumed_bytes)
    return {"generator": identity["generator"], "parameters": identity["parameters"],
            "word_position": identity["start_word"] + consumed // 8, "byte_aligned": consumed % 8 == 0}


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
    body = {"schema": CHECKPOINT_SCHEMA, "scheme_sha256": plan.scheme_sha256,
            "grid_sha256": plan.grid_sha256, "intervals": int(plan.deltas.size), "layout": plan.layout(),
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
    malformed running digest; a stream whose generator, seed, stream id or word position differs, or
    a source that is not replayable.
    """
    if not isinstance(checkpoint, Mapping):
        raise ValueError("checkpoint must be a mapping")
    expected = make_checkpoint(plan, initial, 0, np.zeros(plan.n_nuclides, dtype=np.int64), b"", 0, {})
    if set(checkpoint) != set(expected):
        raise ValueError("checkpoint keys do not match the schema")
    body = {key: value for key, value in checkpoint.items() if key != "checkpoint_sha256"}
    if checkpoint["checkpoint_sha256"] != _checkpoint_checksum(body):
        raise ValueError("checkpoint checksum does not match its contents (corrupted or edited)")
    for key in ("schema", "scheme_sha256", "grid_sha256", "intervals", "layout", "initial"):
        if checkpoint[key] != expected[key]:
            raise ValueError(f"checkpoint does not match this run: {key}")
    next_interval = _strict_int(checkpoint["next_interval"], "next_interval")
    if not 0 <= next_interval <= plan.deltas.size:
        raise ValueError("checkpoint interval out of range")
    if _strict_int(checkpoint["words_consumed"], "words_consumed") != next_interval * plan.words_per_interval:
        raise ValueError("checkpoint word count does not match the layout")
    raw = checkpoint["populations"]
    if not isinstance(raw, list) or len(raw) != plan.n_nuclides:
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
    The run owns the stream's source: a read of it by anyone else while the run
    lasts (a sink or checkpoint callback included) raises RuntimeError.
    Checkpoint and resume need a replayable generator (``replay_identity``).
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
    if scheme.sha256() != plan.scheme_sha256:
        raise ValueError("the scheme has changed since the plan was built")
    if isinstance(checkpoint_every, bool) or not isinstance(checkpoint_every, int) or checkpoint_every < 0:
        raise ValueError("checkpoint_every must be a nonnegative integer")
    if bool(checkpoint_every) != (on_checkpoint is not None):
        raise ValueError("checkpoint_every and on_checkpoint must be given together")
    if checkpoint_every or resume_from is not None:
        stream_identity(stream)  # refuses a source that is not replayable, before any random number is read
    initial = MappingProxyType({n.name: int(c) for n, c in zip(plan.nuclides, populations)})
    running = _initial_digest(plan, populations)
    first_interval = 0
    if resume_from is not None:
        first_interval, populations, running = _resume_state(plan, initial, resume_from, stream)
    first_word = stream.words_read
    ledger = _Ledger(stream)
    kept: dict[str, list] = {"populations": [populations.copy()], "path_counts": [], "decays": [], "photons": [],
                             "neutrons": []}
    for interval in range(first_interval, plan.deltas.size):
        expected_words = (interval + 1 - first_interval) * plan.words_per_interval
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
        decays = (counts @ plan.category_steps).reshape(plan.n_nuclides, plan.n_channels)
        if plan.yield_channels:
            n = np.array([decays[j, k] for j, k, _ in plan.yield_channels], dtype=np.int64)
            y = np.array([yield_ for _, _, yield_ in plan.yield_channels])
            flip = y > 0.5
            sampled = _vector_binomials(n, np.where(flip, 1.0 - y, y), flip, u[plan.photon_slot_offset:],
                                        plan.parts, plan.part_max)
        else:
            sampled = np.zeros(0, dtype=np.int64)
        record = _interval_record(plan, populations, counts, _photons_from(plan, decays, sampled))
        _check_interval(plan, interval, record)
        ledger.check(expected_words, f"during interval {interval}")
        running = _chain_digest(running, interval, record)
        populations = record["populations"]
        if sink is not None:  # detached copies: nothing a sink does to them can reach the engine or its outputs
            sink(interval, {key: np.array(value, copy=True) for key, value in record.items()})
            ledger.check(expected_words, f"by the sink after interval {interval}")
        if keep_outputs:
            for key in kept:
                kept[key].append(record[key])
        if checkpoint_every and (interval + 1) % checkpoint_every == 0:
            on_checkpoint(make_checkpoint(plan, initial, interval + 1, populations, running,
                                          (interval + 1) * plan.words_per_interval, stream_identity(stream)))
            ledger.check(expected_words, f"by the checkpoint callback after interval {interval}")
    return _finish(plan, initial, first_interval, first_word, "vectorised", kept, populations, running, stream,
                   keep_outputs)


class _Ledger:
    """The words a run has read against the bytes its source has emitted since the run began."""

    def __init__(self, stream: WordStream):
        self.stream = stream
        self.words = stream.words_read
        self.emitted = int(stream.source.consumed_bytes)

    def check(self, expected_words: int, where: str) -> None:
        words = self.stream.words_read - self.words
        emitted = int(self.stream.source.consumed_bytes) - self.emitted
        if words != expected_words or emitted != 8 * expected_words:
            raise RuntimeError(f"the run's source was read outside its layout {where}: the layout accounts for "
                               f"{expected_words} words, the run's stream has read {words} and the source has "
                               f"emitted {emitted} bytes; the run's random numbers are no longer the layout's")


def _check_interval(plan: AggregatePlan, interval: int, record: Mapping[str, np.ndarray]) -> None:
    """Conservation, signs and photons within decays after every interval (exact, in Python integers)."""
    populations = record["populations"].tolist()
    if (sum(populations) != plan.total or min(populations) < 0 or int(record["path_counts"].min(initial=0)) < 0
            or int(record["photons"].min(initial=0)) < 0 or bool(np.any(record["photons"] > record["decays"]))):
        raise RuntimeError(f"interval {interval} violates conservation, yields a negative count or more photons "
                           f"than decays (populations {populations}, total {plan.total})")


def _finish(plan, initial, first_interval, first_word, engine, kept, populations, running, stream, keep_outputs):
    width = (plan.n_nuclides, plan.n_channels)
    if keep_outputs:
        arrays = {"populations": np.array(kept["populations"], dtype=np.int64),
                  "path_counts": np.array(kept["path_counts"], dtype=np.int64).reshape(-1, plan.n_categories),
                  "decays": np.array(kept["decays"], dtype=np.int64).reshape(-1, *width),
                  "photons": np.array(kept["photons"], dtype=np.int64).reshape(-1, *width),
                  "neutrons": np.array(kept["neutrons"], dtype=np.int64)}
    else:
        arrays = {"populations": np.zeros((0, plan.n_nuclides), dtype=np.int64),
                  "path_counts": np.zeros((0, plan.n_categories), dtype=np.int64),
                  "decays": np.zeros((0, *width), dtype=np.int64), "photons": np.zeros((0, *width), dtype=np.int64),
                  "neutrons": np.zeros(0, dtype=np.int64)}
    arrays = {key: _frozen(value) for key, value in arrays.items()}
    provenance = {**stream.provenance(), "replay_identity": replay_identity(stream.source),
                  "first_word_of_run": first_word, "words_read_by_run": stream.words_read - first_word}
    return AggregateRun(plan, MappingProxyType(dict(initial)), first_interval, first_word, engine,
                        arrays["populations"], arrays["path_counts"], arrays["decays"], arrays["photons"],
                        arrays["neutrons"], _frozen(populations), running.hex(), stream.words_read - first_word,
                        json.dumps(provenance, sort_keys=True, default=str))


def simulate_aggregate_reference(scheme: DecayScheme, initial: Mapping[str, int], times: Sequence[float],
                                 stream: WordStream, *, part_max: int = PART_MAX) -> AggregateRun:
    """Pure-Python reference engine: one slot at a time, same layout; for cross-checks only.

    It shares the plan (laws and conditional probabilities) and the inversion
    routine with ``simulate_aggregate``; the bookkeeping (slot order, splitting,
    flips, photon and neutron accounting, digests) is written out again here.
    Like the vectorised engine it reads each interval's words before sampling,
    so a failure inside an interval leaves the stream at that interval's end.
    """
    populations = _check_initial(scheme, initial)
    plan = AggregatePlan(scheme, times, int(populations.sum()), part_max=part_max)
    nuclides = plan.nuclides
    initial = {n.name: int(c) for n, c in zip(nuclides, populations)}
    running = _initial_digest(plan, populations)
    first_word = stream.words_read
    ledger = _Ledger(stream)
    kept: dict[str, list] = {"populations": [populations.copy()], "path_counts": [], "decays": [], "photons": [],
                             "neutrons": []}
    width = plan.n_channels
    for interval in range(plan.deltas.size):
        sources = plan.interval_plans[float(plan.deltas[interval])]
        words = stream.uniforms(plan.words_per_interval)
        slot = 0
        counts = [0] * plan.n_categories
        for i, cats in enumerate(plan.categories):
            start = int(plan.first_category[i])
            if len(cats) == 1:
                counts[start] = int(populations[i])
                continue
            remaining = int(populations[i])
            for k in range(len(cats) - 1):
                value = binomial_draw(remaining, float(sources[i].p_small[k]), bool(sources[i].flip[k]),
                                      words[slot * plan.parts:(slot + 1) * plan.parts], part_max=plan.part_max)
                slot += 1
                counts[start + k] = value
                remaining -= value
            counts[start + len(cats) - 1] = remaining
        decays = [[0] * width for _ in nuclides]
        new_populations = [0] * len(nuclides)
        for i, cats in enumerate(plan.categories):
            for k, category in enumerate(cats):
                c = counts[int(plan.first_category[i]) + k]
                new_populations[category.end] += c
                for parent, channel in category.steps:
                    decays[parent][channel] += c
        photons = [[0] * width for _ in nuclides]
        for j, nuclide in enumerate(nuclides):
            for k, channel in enumerate(nuclide.channels):
                if channel.photon_yield == 1.0:
                    photons[j][k] = decays[j][k]
        for j, k, y in plan.yield_channels:
            flip = y > 0.5
            photons[j][k] = binomial_draw(decays[j][k], 1.0 - y if flip else y, flip,
                                          words[slot * plan.parts:(slot + 1) * plan.parts], part_max=plan.part_max)
            slot += 1
        if slot != plan.slots:
            raise RuntimeError(f"reference engine used {slot} slots, the layout has {plan.slots}")
        neutrons = sum(decays[j][k] for j, n in enumerate(nuclides) for k, c in enumerate(n.channels)
                       if c.mode is Mode.BETA_DELAYED_NEUTRON)
        record = {"populations": np.array(new_populations, dtype=np.int64),
                  "path_counts": np.array(counts, dtype=np.int64), "decays": np.array(decays, dtype=np.int64),
                  "photons": np.array(photons, dtype=np.int64), "neutrons": np.int64(neutrons)}
        _check_interval(plan, interval, record)
        ledger.check((interval + 1) * plan.words_per_interval, f"during interval {interval}")
        running = _chain_digest(running, interval, record)
        populations = record["populations"]
        for key in kept:
            kept[key].append(record[key])
    return _finish(plan, initial, 0, first_word, "reference", kept, populations, running, stream, True)


def aggregate_manifest(*, experiment: str, run_id: str, run: AggregateRun, code: Mapping[str, Any],
                       extra: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Manifest of one aggregated run: configuration, stream provenance, output digest.

    The stream provenance is the record the run took when it ended (source manifest, replay
    identity, words and bytes read, the source's SHA-256 of every byte it emitted), so no other
    stream, and no later state of the same stream, can be attributed to the run.
    """
    if not isinstance(run, AggregateRun):
        raise TypeError("run must be an AggregateRun")
    if run.plan.scheme.sha256() != run.plan.scheme_sha256:
        raise ValueError("the scheme has changed since the run; its manifest would misreport the configuration")
    return {"schema": MANIFEST_SCHEMA, "experiment": experiment, "run_id": run_id, "code": dict(code),
            "environment": environment(),
            "rng": {"role": "physics", **run.plan.layout(), **run.stream_provenance()},
            "configuration": {"scheme": run.plan.scheme.to_dict(), "scheme_sha256": run.plan.scheme_sha256,
                              "initial": dict(run.initial), "times": [repr(t) for t in run.plan.times.tolist()],
                              "engine": run.engine},
            "outputs": run.summary(), "extra": dict(extra or {})}
