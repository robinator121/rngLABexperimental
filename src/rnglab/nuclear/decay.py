"""RNGLab v0.20 nuclear module, Layer A: event-level radioactive decay.

Synthetic nuclides only. Half-lives are in abstract time units. Z, A and the
isomer level are bookkeeping labels; the schemes used in RNGLab put Z at 140 or
above, outside the chart of known nuclides, so that no scheme can be mistaken
for nuclear data.

Physics model and its assumptions
---------------------------------
- Nuclei decay independently of each other and of the environment.
- In a state with half-life T the waiting time is exponential with rate
  lambda = ln 2 / T. The channel is then chosen with its branching ratio,
  independently of the waiting time. This is equivalent to competing
  independent exponential clocks with partial rates b_k * lambda.
- A decay is an instantaneous event. A photon that accompanies a decay
  (prompt gamma, the gamma of an isomeric transition, the X-ray after electron
  capture) carries the decay's timestamp and is emitted with the channel's
  photon yield.
- Only the beta-delayed neutron channel (``Mode.BETA_DELAYED_NEUTRON``) emits a
  neutron. No other mode does.
- Schemes are acyclic, so a nucleus passes through each nuclide at most once
  and makes at most D decays, D being the longest decay path.
- Finite half-lives must lie in HALF_LIFE_RANGE. Branching ratios must sum to 1
  within 1e-12 and are never renormalised; the law used by the samplers and by
  every analytic reference is ``Nuclide.channel_probabilities()``.
- Not modelled: internal conversion and Auger electrons (photon yield below 1
  stands in for them), bremsstrahlung, neutrinos, positron annihilation
  photons, recoil, energies, geometry. Detector effects belong to Layer C.

Random-number layout "nuclear-decay-layout-v1"
----------------------------------------------
- Nuclei are numbered 0..N-1 in scheme order of their initial nuclide.
- Nucleus i owns words [i*K, (i+1)*K) of its stream, K = 3*D.
- Decay step r of a nucleus uses its word 3r for the waiting time
  (tau = -log1p(-u) / lambda), word 3r+1 for the channel (the first channel k
  whose cumulative branching ratio exceeds u) and word 3r+2 for the photon
  (emitted iff u < photon yield). Unused words are read and discarded.

The layout makes consumption independent of the outcome (N*K words for every
history), makes chunked and unchunked runs agree bit for bit, and maps every
stream position to one (nucleus, step, purpose). Uniforms are multiples of
2**-53 below 1, so a sampled waiting time is at most 53 ln 2 / lambda (about
36.7 mean lives); the sampled law differs from the exponential by at most
2**-53 in distribution function.
"""
from __future__ import annotations

import enum
import hashlib
import json
import math
from dataclasses import dataclass, field
from typing import Any, Iterator, Mapping, Sequence

import numpy as np

from .streams import WordStream

LAYOUT = "nuclear-decay-layout-v1"
WORDS_PER_STEP = 3
LN2 = math.log(2.0)
BRANCHING_TOLERANCE = 1e-12
HALF_LIFE_RANGE = (1e-300, 1e300)  # finite half-lives outside it are refused (rates and waiting times stay finite)


class Mode(enum.Enum):
    ALPHA = "alpha"
    BETA_MINUS = "beta-"
    BETA_PLUS = "beta+"
    ELECTRON_CAPTURE = "ec"
    ISOMERIC_TRANSITION = "it"
    BETA_DELAYED_NEUTRON = "beta-n"


@dataclass(frozen=True)
class ModeSpec:
    code: int
    delta_z: int
    delta_a: int
    particles: tuple[str, ...]
    photon: str


MODE_SPECS: dict[Mode, ModeSpec] = {
    Mode.ALPHA: ModeSpec(1, -2, -4, ("alpha",), "gamma"),
    Mode.BETA_MINUS: ModeSpec(2, 1, 0, ("electron",), "gamma"),
    Mode.BETA_PLUS: ModeSpec(3, -1, 0, ("positron",), "gamma"),
    Mode.ELECTRON_CAPTURE: ModeSpec(4, -1, 0, (), "xray"),
    Mode.ISOMERIC_TRANSITION: ModeSpec(5, 0, 0, (), "gamma"),
    Mode.BETA_DELAYED_NEUTRON: ModeSpec(6, 1, -1, ("electron", "neutron"), "gamma"),
}
MODE_BY_CODE = {spec.code: mode for mode, spec in MODE_SPECS.items()}
PARTICLES = ("alpha", "electron", "positron", "gamma", "xray", "neutron")
PARTICLE_CODES = {name: code for code, name in enumerate(PARTICLES, start=1)}
NEUTRON_MODE_CODE = MODE_SPECS[Mode.BETA_DELAYED_NEUTRON].code

EVENT_COLUMNS = ("event_id", "nucleus", "step", "time", "parent", "daughter",
                 "channel", "mode", "photon", "neutron", "previous_event")
EVENT_DTYPES = {"event_id": "<i8", "nucleus": "<i8", "step": "<i2", "time": "<f8", "parent": "<i2",
                "daughter": "<i2", "channel": "<i2", "mode": "|i1", "photon": "|b1", "neutron": "|b1",
                "previous_event": "<i8"}


def _real(value: object, what: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.integer, np.floating)):
        raise TypeError(f"{what} must be a real number")
    value = float(value)
    if math.isnan(value):
        raise ValueError(f"{what} must not be NaN")
    return value


def _integer(value: object, what: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise TypeError(f"{what} must be an integer")
    return int(value)


@dataclass(frozen=True)
class Channel:
    """One decay channel: daughter nuclide, mode, branching ratio, photon yield."""

    daughter: str
    mode: Mode
    branching: float
    photon_yield: float = 0.0

    def __post_init__(self):
        if not isinstance(self.daughter, str) or not self.daughter:
            raise ValueError("daughter must be a nonempty name")
        if not isinstance(self.mode, Mode):
            raise TypeError("mode must be a Mode")
        branching = _real(self.branching, "branching ratio")
        photon_yield = _real(self.photon_yield, "photon yield")
        if not 0.0 < branching <= 1.0:
            raise ValueError("branching ratio must lie in (0, 1]")
        if not 0.0 <= photon_yield <= 1.0:
            raise ValueError("photon yield must lie in [0, 1]")
        object.__setattr__(self, "branching", branching)
        object.__setattr__(self, "photon_yield", photon_yield)


@dataclass(frozen=True)
class Nuclide:
    """A synthetic nuclide state. ``half_life = math.inf`` means stable."""

    name: str
    z: int
    a: int
    half_life: float
    channels: tuple[Channel, ...] = ()
    level: int = 0

    def __post_init__(self):
        if not isinstance(self.name, str) or not self.name:
            raise ValueError("nuclide name must be a nonempty string")
        z, a, level = (_integer(self.z, "Z"), _integer(self.a, "A"), _integer(self.level, "level"))
        if not 0 <= z <= a or level < 0:
            raise ValueError(f"{self.name}: need 0 <= Z <= A and level >= 0")
        half_life = _real(self.half_life, "half-life")
        if not half_life > 0.0:
            raise ValueError(f"{self.name}: half-life must be positive (math.inf for stable)")
        if math.isfinite(half_life) and not HALF_LIFE_RANGE[0] <= half_life <= HALF_LIFE_RANGE[1]:
            raise ValueError(f"{self.name}: finite half-lives must lie in {HALF_LIFE_RANGE}")
        channels = tuple(self.channels)
        if not all(isinstance(channel, Channel) for channel in channels):
            raise TypeError(f"{self.name}: channels must be Channel objects")
        if math.isinf(half_life) != (not channels):
            raise ValueError(f"{self.name}: a nuclide is stable if and only if it has no channels")
        if channels:
            total = math.fsum(channel.branching for channel in channels)
            if abs(total - 1.0) > BRANCHING_TOLERANCE:
                raise ValueError(f"{self.name}: branching ratios sum to {total!r}, not 1 "
                                 f"(tolerance {BRANCHING_TOLERANCE}); they are not renormalised")
            keys = [(channel.daughter, channel.mode) for channel in channels]
            if len(set(keys)) != len(keys):
                raise ValueError(f"{self.name}: duplicate (daughter, mode) channel")
        object.__setattr__(self, "z", z)
        object.__setattr__(self, "a", a)
        object.__setattr__(self, "level", level)
        object.__setattr__(self, "half_life", half_life)
        object.__setattr__(self, "channels", channels)
        if channels and not all(p > 0.0 for p in self.channel_probabilities()):
            raise ValueError(f"{self.name}: a channel has no probability left under the cumulative thresholds "
                             "(the declared ratios round it away)")

    @property
    def stable(self) -> bool:
        return math.isinf(self.half_life)

    def channel_probabilities(self) -> tuple[float, ...]:
        """The branching law the sampler uses, and the one every analytic reference uses.

        Channel k is chosen when c_{k-1} <= u < c_k, with c_k the running sum of
        the declared ratios and the last threshold set to 1. Its probability is
        c_k - c_{k-1} (to within 2**-53, the uniform resolution); it differs from
        the declared ratio by at most the normalisation tolerance. A stable
        nuclide has no channels and an empty tuple.
        """
        if not self.channels:
            return ()
        running, edges = 0.0, [0.0]
        for channel in self.channels[:-1]:
            running += channel.branching
            edges.append(running)
        edges.append(1.0)
        return tuple(high - low for low, high in zip(edges, edges[1:]))


def _read_only(array: np.ndarray) -> np.ndarray:
    array.setflags(write=False)
    return array


class DecayScheme:
    """A validated, acyclic set of synthetic nuclides and their decay channels."""

    def __init__(self, name: str, nuclides: Sequence[Nuclide]):
        if not isinstance(name, str) or not name:
            raise ValueError("scheme name must be a nonempty string")
        nuclides = tuple(nuclides)
        if not nuclides or not all(isinstance(n, Nuclide) for n in nuclides):
            raise ValueError("a scheme needs at least one Nuclide")
        names = [n.name for n in nuclides]
        if len(set(names)) != len(names):
            raise ValueError("nuclide names must be unique")
        states = [(n.z, n.a, n.level) for n in nuclides]
        if len(set(states)) != len(states):
            raise ValueError("two nuclides share (Z, A, level)")
        index = {n.name: i for i, n in enumerate(nuclides)}
        for parent in nuclides:
            for channel in parent.channels:
                if channel.daughter not in index:
                    raise ValueError(f"{parent.name}: unknown daughter {channel.daughter!r}")
                daughter = nuclides[index[channel.daughter]]
                spec = MODE_SPECS[channel.mode]
                if (daughter.z - parent.z, daughter.a - parent.a) != (spec.delta_z, spec.delta_a):
                    raise ValueError(f"{parent.name} -> {daughter.name}: (dZ, dA) = "
                                     f"({daughter.z - parent.z}, {daughter.a - parent.a}) contradicts "
                                     f"mode {channel.mode.value} ({spec.delta_z}, {spec.delta_a})")
                if channel.mode is Mode.ISOMERIC_TRANSITION and not daughter.level < parent.level:
                    raise ValueError(f"{parent.name}: an isomeric transition must lower the level")
        self.name = name
        self.nuclides = nuclides
        self.index = index
        self.depth = self._longest_path()

        size = len(nuclides)
        width = max((len(n.channels) for n in nuclides), default=0)
        self.max_channels = width
        thresholds = np.full((size, max(width - 1, 1)), np.inf)
        daughters = np.full((size, max(width, 1)), -1, dtype=np.int16)
        modes = np.zeros((size, max(width, 1)), dtype=np.int8)
        yields = np.zeros((size, max(width, 1)))
        for i, nuclide in enumerate(nuclides):
            cumulative = 0.0
            for k, channel in enumerate(nuclide.channels):
                daughters[i, k] = index[channel.daughter]
                modes[i, k] = MODE_SPECS[channel.mode].code
                yields[i, k] = channel.photon_yield
                if k < len(nuclide.channels) - 1:
                    cumulative += channel.branching
                    thresholds[i, k] = cumulative
        self.decay_constants = _read_only(np.array([0.0 if n.stable else LN2 / n.half_life for n in nuclides]))
        self.cum_thresholds = _read_only(thresholds)
        self.daughter_index = _read_only(daughters)
        self.mode_code = _read_only(modes)
        self.photon_yield = _read_only(yields)

    def _longest_path(self) -> int:
        depth: dict[int, int] = {}
        visiting: set[int] = set()

        def visit(i: int) -> int:
            if i in depth:
                return depth[i]
            if i in visiting:
                raise ValueError(f"decay scheme has a cycle through {self.nuclides[i].name}")
            visiting.add(i)
            best = 0
            for channel in self.nuclides[i].channels:
                best = max(best, 1 + visit(self.index[channel.daughter]))
            visiting.discard(i)
            depth[i] = best
            return best

        return max(visit(i) for i in range(len(self.nuclides)))

    def __len__(self) -> int:
        return len(self.nuclides)

    @property
    def words_per_nucleus(self) -> int:
        return WORDS_PER_STEP * self.depth

    def initial_states(self, initial: Mapping[str, int]) -> np.ndarray:
        """Nucleus states in layout order: scheme order of the initial nuclide."""
        if not isinstance(initial, Mapping):
            raise TypeError("initial populations must be a mapping name -> count")
        unknown = set(initial) - set(self.index)
        if unknown:
            raise ValueError(f"unknown nuclides in initial populations: {sorted(unknown)}")
        counts = []
        for nuclide in self.nuclides:
            count = initial.get(nuclide.name, 0)
            if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 0:
                raise ValueError(f"initial count of {nuclide.name} must be a nonnegative integer")
            counts.append(int(count))
        return np.repeat(np.arange(len(self.nuclides), dtype=np.int16), counts)

    def paths(self, start: str) -> Iterator[tuple[tuple[int, ...], float]]:
        """Every decay-path prefix from ``start``: (nuclide indices, product of branching ratios).

        The empty prefix (no decay) is included. Two channels with the same
        daughter give two prefixes with the same nuclide sequence.
        """
        def walk(path: tuple[int, ...], weight: float):
            yield path, weight
            nuclide = self.nuclides[path[-1]]
            for channel, probability in zip(nuclide.channels, nuclide.channel_probabilities()):
                yield from walk(path + (self.index[channel.daughter],), weight * probability)

        yield from walk((self.index[start],), 1.0)

    def precursors(self) -> dict[str, float]:
        """Delayed-neutron precursors of the scheme: name -> neutron emission probability per decay."""
        return {n.name: math.fsum(p for c, p in zip(n.channels, n.channel_probabilities())
                                  if c.mode is Mode.BETA_DELAYED_NEUTRON)
                for n in self.nuclides
                if any(c.mode is Mode.BETA_DELAYED_NEUTRON for c in n.channels)}

    def to_dict(self) -> dict[str, Any]:
        return {"name": self.name, "nuclides": [
            {"name": n.name, "z": n.z, "a": n.a, "level": n.level,
             "half_life": None if n.stable else n.half_life,
             "channels": [{"daughter": c.daughter, "mode": c.mode.value, "branching": c.branching,
                           "photon_yield": c.photon_yield} for c in n.channels]}
            for n in self.nuclides]}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> "DecayScheme":
        nuclides = []
        for item in data["nuclides"]:
            channels = tuple(Channel(c["daughter"], Mode(c["mode"]), c["branching"], c.get("photon_yield", 0.0))
                             for c in item.get("channels", ()))
            half_life = math.inf if item["half_life"] is None else item["half_life"]
            nuclides.append(Nuclide(item["name"], item["z"], item["a"], half_life, channels, item.get("level", 0)))
        return cls(data["name"], nuclides)

    def sha256(self) -> str:
        text = json.dumps(self.to_dict(), sort_keys=True, separators=(",", ":"), allow_nan=False)
        return hashlib.sha256(text.encode()).hexdigest()


def _check_horizon(horizon: object) -> float:
    horizon = _real(horizon, "horizon")
    if horizon < 0.0:
        raise ValueError("horizon must be nonnegative (math.inf allowed)")
    return horizon


def _canonical_digest(parts: Sequence[tuple[str, np.ndarray, str]], header: Mapping[str, Any]) -> str:
    digest = hashlib.sha256(json.dumps(header, sort_keys=True, separators=(",", ":")).encode())
    for name, array, dtype in parts:
        array = np.ascontiguousarray(np.asarray(array).astype(dtype, copy=False))
        digest.update(f"\0{name}\0{dtype}\0{array.size}\0".encode())
        digest.update(array.tobytes())
    return digest.hexdigest()


@dataclass(frozen=True, eq=False)
class DecayRun:
    """Result of one decay simulation. ``events`` is sorted by event_id.

    event_id = nucleus * D + step, so identifiers are fixed by the layout and
    the same in every engine; ``previous_event`` is the event that created the
    nucleus's current state (-1 for its first decay).
    """

    scheme: DecayScheme
    horizon: float
    initial_state: np.ndarray
    final_state: np.ndarray
    events: Mapping[str, np.ndarray]
    words_per_nucleus: int
    first_word: int
    engine: str
    layout: str = field(default=LAYOUT)

    @property
    def n_nuclei(self) -> int:
        return int(self.initial_state.size)

    @property
    def n_events(self) -> int:
        return int(self.events["event_id"].size)

    @property
    def words_used(self) -> int:
        return self.n_nuclei * self.words_per_nucleus

    def populations(self, times: Sequence[float]) -> np.ndarray:
        """Counts per nuclide at each time (state after every event at or before t)."""
        times = np.atleast_1d(np.asarray(times, dtype=np.float64))
        if np.any(np.isnan(times)) or np.any(times < 0) or np.any(times > self.horizon):
            raise ValueError("population times must lie in [0, horizon]")
        size = len(self.scheme)
        out = np.empty((times.size, size), dtype=np.int64)
        start = np.bincount(self.initial_state, minlength=size)
        when, parent, daughter = self.events["time"], self.events["parent"], self.events["daughter"]
        for j in range(size):
            entries = np.sort(when[daughter == j])
            exits = np.sort(when[parent == j])
            out[:, j] = (start[j] + np.searchsorted(entries, times, side="right")
                         - np.searchsorted(exits, times, side="right"))
        return out

    def decay_counts(self) -> np.ndarray:
        """Decays per parent nuclide within [0, horizon]."""
        return np.bincount(self.events["parent"], minlength=len(self.scheme)).astype(np.int64)

    def channel_counts(self) -> np.ndarray:
        """Decays per (parent, channel), shape (nuclides, max channels)."""
        width = max(self.scheme.max_channels, 1)
        flat = self.events["parent"].astype(np.int64) * width + self.events["channel"]
        return np.bincount(flat, minlength=len(self.scheme) * width).reshape(len(self.scheme), width)

    def photon_counts(self) -> np.ndarray:
        """Photons per (parent, channel), shape (nuclides, max channels)."""
        width = max(self.scheme.max_channels, 1)
        flat = self.events["parent"].astype(np.int64) * width + self.events["channel"]
        return np.bincount(flat[self.events["photon"]], minlength=len(self.scheme) * width).reshape(
            len(self.scheme), width)

    def decay_times(self, nuclide: str) -> np.ndarray:
        """Times of the decays of one nuclide, in event_id order."""
        return self.events["time"][self.events["parent"] == self.scheme.index[nuclide]]

    def emissions(self) -> dict[str, np.ndarray]:
        """Physical emissions (no detector effects), sorted by (time, event_id, particle).

        Columns: time, particle (code in PARTICLE_CODES), event_id, nucleus.
        """
        mode, photon = self.events["mode"], self.events["photon"]
        select = {
            "alpha": mode == MODE_SPECS[Mode.ALPHA].code,
            "electron": (mode == MODE_SPECS[Mode.BETA_MINUS].code) | (mode == NEUTRON_MODE_CODE),
            "positron": mode == MODE_SPECS[Mode.BETA_PLUS].code,
            "gamma": photon & (mode != MODE_SPECS[Mode.ELECTRON_CAPTURE].code),
            "xray": photon & (mode == MODE_SPECS[Mode.ELECTRON_CAPTURE].code),
            "neutron": mode == NEUTRON_MODE_CODE,
        }
        time = np.concatenate([self.events["time"][m] for m in select.values()])
        particle = np.concatenate([np.full(int(m.sum()), PARTICLE_CODES[name], dtype=np.int8)
                                   for name, m in select.items()])
        event_id = np.concatenate([self.events["event_id"][m] for m in select.values()])
        nucleus = np.concatenate([self.events["nucleus"][m] for m in select.values()])
        order = np.lexsort((particle, event_id, time))
        return {"time": time[order], "particle": particle[order], "event_id": event_id[order],
                "nucleus": nucleus[order]}

    def digest(self) -> str:
        """SHA-256 of the scheme, horizon, layout, initial and final states and every event column."""
        header = {"scheme_sha256": self.scheme.sha256(), "horizon": repr(self.horizon), "layout": self.layout,
                  "words_per_nucleus": self.words_per_nucleus, "n_nuclei": self.n_nuclei}
        parts = [("initial_state", self.initial_state, "<i2"), ("final_state", self.final_state, "<i2")]
        parts += [(name, self.events[name], EVENT_DTYPES[name]) for name in EVENT_COLUMNS]
        return _canonical_digest(parts, header)

    def summary(self) -> dict[str, Any]:
        names = [n.name for n in self.scheme.nuclides]
        return {"engine": self.engine, "layout": self.layout, "n_nuclei": self.n_nuclei,
                "n_events": self.n_events, "words_per_nucleus": self.words_per_nucleus,
                "first_word": self.first_word, "words_used": self.words_used,
                "final_counts": dict(zip(names, np.bincount(self.final_state, minlength=len(names)).tolist())),
                "decay_counts": dict(zip(names, self.decay_counts().tolist())),
                "events_sha256": self.digest()}


def _assemble(columns: dict[str, list[np.ndarray]], depth: int) -> dict[str, np.ndarray]:
    events = {}
    for name in ("event_id", "nucleus", "step", "time", "parent", "daughter", "channel", "mode", "photon"):
        parts = columns[name]
        events[name] = (np.concatenate(parts) if parts else np.empty(0)).astype(EVENT_DTYPES[name])
    order = np.argsort(events["event_id"], kind="stable")
    events = {name: values[order] for name, values in events.items()}
    events["neutron"] = events["mode"] == NEUTRON_MODE_CODE
    events["previous_event"] = np.where(events["step"] > 0, events["event_id"] - 1, -1).astype(np.int64)
    for values in events.values():
        values.setflags(write=False)
    return {name: events[name] for name in EVENT_COLUMNS}


def simulate_decays(scheme: DecayScheme, initial: Mapping[str, int], horizon: float, stream: WordStream,
                    *, chunk: int = 65536) -> DecayRun:
    """Event-level decay simulation, vectorised over nuclei in chunks (memory bounded)."""
    if not isinstance(scheme, DecayScheme):
        raise TypeError("scheme must be a DecayScheme")
    if not isinstance(stream, WordStream):
        raise TypeError("stream must be a WordStream")
    horizon = _check_horizon(horizon)
    if isinstance(chunk, bool) or not isinstance(chunk, int) or chunk < 1:
        raise ValueError("chunk must be a positive integer")
    initial_state = scheme.initial_states(initial)
    final_state = initial_state.copy()
    depth, width = scheme.depth, scheme.words_per_nucleus
    first_word = stream.words_read
    rates = scheme.decay_constants
    columns: dict[str, list[np.ndarray]] = {name: [] for name in EVENT_COLUMNS}
    total = initial_state.size
    for start in range(0, total if width else 0, chunk):
        stop = min(total, start + chunk)
        u = stream.uniforms((stop - start) * width).reshape(stop - start, width)
        state = final_state[start:stop]
        clock = np.zeros(stop - start)
        active = np.ones(stop - start, dtype=bool)
        for step in range(depth):
            rate = rates[state]
            candidates = np.flatnonzero(active & (rate > 0.0))
            if candidates.size == 0:
                break
            new_time = clock[candidates] + (-np.log1p(-u[candidates, WORDS_PER_STEP * step]) / rate[candidates])
            inside = new_time <= horizon
            active[candidates[~inside]] = False
            decaying = candidates[inside]
            if decaying.size == 0:
                break
            parent = state[decaying]
            channel = np.count_nonzero(
                u[decaying, WORDS_PER_STEP * step + 1][:, None] >= scheme.cum_thresholds[parent], axis=1)
            daughter = scheme.daughter_index[parent, channel]
            nucleus = start + decaying.astype(np.int64)
            columns["event_id"].append(nucleus * depth + step)
            columns["nucleus"].append(nucleus)
            columns["step"].append(np.full(decaying.size, step))
            columns["time"].append(new_time[inside])
            columns["parent"].append(parent)
            columns["daughter"].append(daughter)
            columns["channel"].append(channel)
            columns["mode"].append(scheme.mode_code[parent, channel])
            columns["photon"].append(u[decaying, WORDS_PER_STEP * step + 2] < scheme.photon_yield[parent, channel])
            state[decaying] = daughter
            clock[decaying] = new_time[inside]
    initial_state.setflags(write=False)
    final_state.setflags(write=False)
    return DecayRun(scheme, horizon, initial_state, final_state, _assemble(columns, depth), width, first_word,
                    engine="vectorised")


def simulate_decays_reference(scheme: DecayScheme, initial: Mapping[str, int], horizon: float,
                              stream: WordStream) -> DecayRun:
    """Pure-Python reference engine: one nucleus at a time, for cross-checks only (slow).

    It shares the scheme definition and the layout with simulate_decays, and
    nothing else: rates, cumulative thresholds, the word-to-uniform conversion
    and the sampling loop are recomputed here.
    """
    horizon = _check_horizon(horizon)
    starts = scheme.initial_states(initial).tolist()
    nuclides = scheme.nuclides
    rates = [0.0 if math.isinf(n.half_life) else math.log(2.0) / n.half_life for n in nuclides]
    thresholds = []
    for nuclide in nuclides:
        running, cuts = 0.0, []
        for channel in nuclide.channels[:-1]:
            running += channel.branching
            cuts.append(running)
        thresholds.append(cuts)
    depth = scheme.depth
    width = 3 * depth
    first_word = stream.words_read
    rows, finals = [], []
    for nucleus, state in enumerate(starts):
        raw = stream.raw(width)
        u = [(int.from_bytes(raw[8 * j:8 * j + 8], "big") >> 11) * 2.0 ** -53 for j in range(width)]
        clock = 0.0
        for step in range(depth):
            if rates[state] == 0.0:
                break
            new_time = clock + (-math.log1p(-u[3 * step]) / rates[state])
            if not new_time <= horizon:
                break
            k = 0
            while k < len(thresholds[state]) and u[3 * step + 1] >= thresholds[state][k]:
                k += 1
            channel = nuclides[state].channels[k]
            daughter = scheme.index[channel.daughter]
            rows.append((nucleus * depth + step, nucleus, step, new_time, state, daughter, k,
                         MODE_SPECS[channel.mode].code, u[3 * step + 2] < channel.photon_yield))
            state, clock = daughter, new_time
        finals.append(state)
    names = ("event_id", "nucleus", "step", "time", "parent", "daughter", "channel", "mode", "photon")
    columns = {name: [np.array([row[i] for row in rows])] if rows else [] for i, name in enumerate(names)}
    initial_state = np.array(starts, dtype=np.int16)
    final_state = np.array(finals, dtype=np.int16)
    initial_state.setflags(write=False)
    final_state.setflags(write=False)
    return DecayRun(scheme, horizon, initial_state, final_state, _assemble(columns, depth), width, first_word,
                    engine="reference")
