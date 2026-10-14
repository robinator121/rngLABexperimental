"""Uniform-marginal dependence processes in *proposal-row* time.

Mathematical models (before floating-point discretization):

* copy: U_t = U_{t-1} with independent probability r_t, otherwise V_t,
  where V_t and U_-1 are independent Uniform(0,1). Induction proves
  uniform marginals for any deterministic r_t schedule, including ramps.
  At constant r, Cov(f(U_t),f(U_{t+k})) = Var(f(U))*r**k. The finite-N
  variance of the average is sigma²/N times
  [1 + 2*sum((1-k/N)*r**k, k=1,...,N-1)], and the infinite-N inflation
  is (1+r)/(1-r), r<1. The mean remains E[f(U)], even at r=1.
* gaussian_ar1: Z_t = r_t*Z_{t-1} + sqrt(1-r_t²)*epsilon_t,
  with independent stationary N(0,1) initialization and innovations.
  Each Z_t remains normal for deterministic schedules; U_t=Phi(Z_t).
  At constant r, Corr(U_t,U_{t+k}) = (6/pi)*asin(r**k/2), which is
  different from the latent-normal correlation r**k.
* lag_copy: the copy process on each residue class modulo lag, initialized
  with independent uniforms. Corr(U_t,U_{t+k}) is r**(k/lag) for positive
  multiples of lag, otherwise zero (constant r).
* block_reuse: each deterministic block starts with a fresh uniform anchor;
  later rows independently reuse it with scheduled probability r_t.
  Marginals are uniform, but the covariance depends on position in a block.
* cross_correlation: Gaussian copula between two columns within a row.
  It has no temporal dependence; Corr(U_a,U_b)=6/pi*asin(r_t/2).

Uniform marginals alone do not establish correct MCMC sampling: conditioned
on the chain history, proposal/acceptance inputs may be dependent and the
usual transition-kernel/detailed-balance argument may fail.

Finite implementation: PCG64 supplies a 53-bit dyadic uniform grid. Copula
models use inverse/CDF transforms in IEEE float64, clipping the inverse at
the smallest grid midpoint and the result below one. Their continuous
marginal proof has this explicit floating-point approximation. Optional
quantize_bits projects every model to the same dyadic grid. These are
reproducible statistical fault models, not entropy or security proofs.
"""
from __future__ import annotations

import copy
import hashlib
import json
from numbers import Integral

import numpy as np
from scipy.signal import lfilter
from scipy.special import ndtr, ndtri


_ALIASES = {"clean": "iid", "copy_previous": "copy", "copy_previous_value": "copy",
            "gaussian": "gaussian_ar1", "ar1": "gaussian_ar1", "copula": "gaussian_ar1",
            "periodic_reuse": "lag_copy", "cross": "cross_correlation"}
_KINDS = {"iid", "copy", "gaussian_ar1", "lag_copy", "block_reuse", "cross_correlation"}
_MIN_INVERSE = 2.0 ** -54
_MAX_UNIFORM = np.nextafter(1.0, 0.0)


def _integer(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def _normal(values):
    return ndtri(np.clip(values, _MIN_INVERSE, _MAX_UNIFORM))


def _schedule_dict(value):
    value = dict(value)
    aliases = {"onset_step": "onset", "onset_steps": "onset", "onset_proposal": "onset",
               "ramp_steps": "ramp", "end_step": "end", "end_steps": "end",
               "recovery_steps": "recovery", "period_steps": "period", "active_steps": "active"}
    for old, new in aliases.items():
        if old in value:
            if new in value and value[new] != value[old]:
                raise ValueError(f"conflicting {old}/{new}")
            value[new] = value.pop(old)
    for key in ("onset", "ramp", "recovery"):
        value[key] = _integer(value.get(key, 0), key)
    if value.get("end") is not None:
        value["end"] = _integer(value["end"], "end")
        if value["end"] <= value["onset"]:
            raise ValueError("end must exceed onset")
    if value.get("period") is not None:
        value["period"] = _integer(value["period"], "period", 1)
        value["active"] = _integer(value.get("active", value["period"]), "active", 1)
        if value["active"] > value["period"]:
            raise ValueError("active must be <= period")
    if value.get("trajectory", "linear") not in {"linear", "smoothstep", "cosine"}:
        raise ValueError("unknown schedule trajectory")
    return value


def schedule_envelope(positions, spec):
    """Deterministic envelope, onset-inclusive, end-exclusive without recovery.

    For ramp R>0, the value is zero at onset and reaches one at onset+R.
    For recovery D>0, end starts recovery and end+D has envelope zero.
    Multiple episodes combine by maximum, never by summing probabilities.
    """
    positions = np.asarray(positions, dtype=np.int64)
    episodes = spec.get("episodes")
    if episodes is not None:
        if not isinstance(episodes, (list, tuple)) or not episodes:
            raise ValueError("episodes must be a nonempty list")
        return np.maximum.reduce([schedule_envelope(positions, s) for s in episodes])
    spec = _schedule_dict(spec)
    onset, ramp = spec["onset"], spec["ramp"]
    envelope = np.clip((positions - onset) / ramp, 0., 1.) if ramp else (positions >= onset).astype(float)
    trajectory = spec.get("trajectory", "linear")
    if trajectory == "smoothstep":
        envelope = envelope * envelope * (3 - 2 * envelope)
    elif trajectory == "cosine":
        envelope = 0.5 * (1 - np.cos(np.pi * envelope))
    if spec.get("end") is not None:
        end, recovery = spec["end"], spec["recovery"]
        down = np.clip((end + recovery - positions) / recovery, 0., 1.) if recovery else (positions < end)
        if trajectory == "smoothstep":
            down = down * down * (3 - 2 * down)
        elif trajectory == "cosine":
            down = 0.5 * (1 - np.cos(np.pi * down))
        envelope *= down
    if spec.get("period") is not None:
        envelope *= (positions - onset) % spec["period"] < spec["active"]
    return envelope


def copy_variance_inflation(rho, n):
    """Exact finite-N inflation of Var(mean(f(U))) for stationary copy.

    Returns N at rho=1; asymptotic (1+rho)/(1-rho) applies only when rho<1.
    This concerns a single-uniform integrand, not an MCMC chain's output.
    """
    n = _integer(n, "n", 1)
    rho = float(rho)
    if not np.isfinite(rho) or not 0 <= rho <= 1:
        raise ValueError("copy rho must be in [0,1]")
    lags = np.arange(1, n, dtype=float)
    return float(1 + 2 * np.dot(1 - lags / n, np.power(rho, lags)))


def gaussian_uniform_correlation(latent_rho):
    """Pearson correlation of two uniform Gaussian-copula coordinates."""
    rho = np.asarray(latent_rho, dtype=float)
    if np.any(~np.isfinite(rho)) or np.any(np.abs(rho) > 1):
        raise ValueError("latent correlation outside [-1,1]")
    return 6 / np.pi * np.arcsin(rho / 2)


def uniform_diagnostics(values, lags=(1, 2, 8)):
    """Descriptive diagnostics only; p-values would require dependence-aware nulls."""
    values = np.asarray(values, dtype=float)
    if values.ndim != 2 or len(values) < 2:
        raise ValueError("values must have shape (rows>=2, columns)")
    correlations = {}
    for lag in lags:
        lag = _integer(lag, "lag", 1)
        if lag >= len(values):
            continue
        row = []
        for column in range(values.shape[1]):
            x, y = values[:-lag, column], values[lag:, column]
            row.append(float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else None)
        correlations[str(lag)] = row
    cross = np.eye(values.shape[1]).tolist()
    for i in range(values.shape[1]):
        for j in range(i):
            a, b = values[:, i], values[:, j]
            corr = float(np.corrcoef(a, b)[0, 1]) if np.std(a) > 0 and np.std(b) > 0 else None
            cross[i][j] = cross[j][i] = corr
    return {"rows": len(values), "mean": values.mean(axis=0).tolist(),
            "variance": values.var(axis=0, ddof=1).tolist(), "lag_correlations": correlations,
            "cross_correlations": cross,
            "adjacent_equal_fraction": (values[1:] == values[:-1]).mean(axis=0).tolist(),
            "conditional_next_means_previous_quartile": [
                [float(values[1:, c][mask].mean()) if np.any(mask) else None
                 for mask in [(values[:-1, c] >= q/4) & (values[:-1, c] < (q+1)/4) for q in range(4)]]
                for c in range(values.shape[1])]}


class UniformProcess:
    """Generate fixed-dimensional rows with independent channel substreams.

    Every call consumes exactly count base and decision draws per column,
    regardless of fault, channel selection, or downstream acceptance outcome.
    A fixed seed/dimension/channels gives identical base draws across specs.
    Sampling in arbitrary chunks produces identical arrays. The seed is a
    reproducibility identity; independent experiments must use disjoint seeds.
    """

    def __init__(self, spec, seed, dimension, channels=None):
        self.dimension = _integer(dimension, "dimension", 1)
        self.spec = copy.deepcopy(dict(spec or {}))
        self.kind = _ALIASES.get(self.spec.get("kind", self.spec.get("process", "iid")),
                                 self.spec.get("kind", self.spec.get("process", "iid")))
        if self.kind not in _KINDS:
            raise ValueError(f"unknown uniform process {self.kind}")
        self.rho = float(self.spec.get("rho", self.spec.get("strength", 0.0)))
        lo = -1 if self.kind in {"gaussian_ar1", "cross_correlation"} else 0
        if not np.isfinite(self.rho) or not lo <= self.rho <= 1:
            raise ValueError(f"rho must be in [{lo},1]")
        defaults = {1: ("proposal",), 2: ("proposal", "acceptance"), 3: ("site", "proposal", "acceptance")}
        self.channels = tuple(channels or defaults.get(self.dimension, tuple(f"column_{i}" for i in range(self.dimension))))
        if len(self.channels) != self.dimension or len(set(self.channels)) != self.dimension:
            raise ValueError("channels must uniquely name every column")
        selected = self.spec.get("channels", self.spec.get("channel", "all"))
        if selected in ("all", None):
            selected = list(range(self.dimension))
        elif isinstance(selected, (str, Integral)):
            selected = [selected]
        if not isinstance(selected, (list, tuple)):
            raise ValueError("channels selection must be names or indices")
        self.selected = tuple(dict.fromkeys(self._column(x) for x in selected))
        self.lag = _integer(self.spec.get("lag", 8), "lag", 1)
        self.block_size = _integer(self.spec.get("block_size", self.spec.get("block", 8)), "block_size", 1)
        bits = self.spec.get("quantize_bits")
        self.quantize_bits = None if bits is None else _integer(bits, "quantize_bits", 1)
        if self.quantize_bits is not None and self.quantize_bits > 53:
            raise ValueError("quantize_bits cannot exceed 53")
        self.schedule = copy.deepcopy(self.spec.get("schedule", self.spec))
        if self.schedule.get("episodes") is not None:
            self.schedule["episodes"] = [_schedule_dict(s) for s in self.schedule["episodes"]]
        else:
            self.schedule = _schedule_dict(self.schedule)
        schedule_envelope(np.empty(0, dtype=np.int64), self.schedule)
        sequence = (np.random.SeedSequence(seed.entropy, spawn_key=seed.spawn_key, pool_size=seed.pool_size)
                    if isinstance(seed, np.random.SeedSequence) else np.random.SeedSequence(seed))
        self._seed_identity = {"entropy": np.asarray(sequence.entropy, dtype=object).tolist(),
                               "spawn_key": list(sequence.spawn_key), "pool_size": sequence.pool_size}
        children = sequence.spawn(3 * self.dimension)
        self._base = [np.random.Generator(np.random.PCG64(s)) for s in children[:self.dimension]]
        self._decisions = [np.random.Generator(np.random.PCG64(s)) for s in children[self.dimension:2*self.dimension]]
        self._initial = [np.random.Generator(np.random.PCG64(s)) for s in children[2*self.dimension:]]
        self._previous = np.array([r.random() for r in self._initial])
        self._latent = _normal(self._previous)
        self._normal_transforms = np.ones(self.dimension, dtype=np.int64)
        self._history = np.column_stack([r.random(self.lag) for r in self._initial])
        self._anchor = self._previous.copy()
        self.position = 0
        self._applied = np.zeros(self.dimension, dtype=np.int64)
        pair = self.spec.get("pair", ("proposal", "acceptance") if {"proposal", "acceptance"} <= set(self.channels) else (0, 1))
        self.pair = tuple(self._column(x) for x in pair) if self.kind == "cross_correlation" else ()
        if self.kind == "cross_correlation" and (len(self.pair) != 2 or self.pair[0] == self.pair[1]):
            raise ValueError("cross_correlation needs two distinct columns")

    def _column(self, value):
        if isinstance(value, str):
            if value not in self.channels:
                raise ValueError(f"unknown channel {value}")
            return self.channels.index(value)
        return _integer(value, "column") if isinstance(value, Integral) and 0 <= value < self.dimension else self._bad_column(value)

    @staticmethod
    def _bad_column(value):
        raise ValueError(f"invalid column {value}")

    def sample(self, count):
        count = _integer(count, "count")
        if count == 0:
            return np.empty((0, self.dimension))
        base = np.column_stack([r.random(count) for r in self._base])
        decisions = np.column_stack([r.random(count) for r in self._decisions])
        positions = np.arange(self.position, self.position + count, dtype=np.int64)
        coefficient = self.rho * schedule_envelope(positions, self.schedule)
        output = base.copy()
        if self.kind in {"copy", "lag_copy"}:
            gates = decisions < coefficient[:, None]
            for column in self.selected:
                if self.kind == "copy":
                    origin = np.maximum.accumulate(np.where(~gates[:, column], np.arange(1, count + 1), 0))
                    output[:, column] = np.concatenate(([self._previous[column]], base[:, column]))[origin]
                else:
                    for residue in range(self.lag):
                        offsets = np.flatnonzero(positions % self.lag == residue)
                        if not len(offsets):
                            continue
                        origin = np.maximum.accumulate(np.where(~gates[offsets, column], np.arange(1, len(offsets) + 1), 0))
                        history = self._history[residue, column]
                        output[offsets, column] = np.concatenate(([history], base[offsets, column]))[origin]
                        self._history[residue, column] = output[offsets[-1], column]
                self._applied[column] += np.count_nonzero(gates[:, column])
        elif self.kind == "gaussian_ar1":
            # Constant-coefficient runs use the validated linear recurrence in
            # scipy.signal; ramp rows are singletons, preserving exact semantics.
            cuts = np.r_[0, np.flatnonzero(np.diff(coefficient) != 0) + 1, count]
            for column in self.selected:
                innovations = _normal(base[:, column])
                self._normal_transforms[column] += count
                latent = np.empty(count)
                previous = self._latent[column]
                for begin, end in zip(cuts[:-1], cuts[1:]):
                    r = coefficient[begin]
                    latent[begin:end], _ = lfilter([np.sqrt(max(0., 1-r*r))], [1., -r],
                                                   innovations[begin:end], zi=[r*previous])
                    previous = latent[end - 1]
                self._latent[column] = previous
                output[:, column] = np.minimum(ndtr(latent), _MAX_UNIFORM)
                # Guarantees rho=0 controls match the IID base bit-for-bit.
                output[coefficient == 0, column] = base[coefficient == 0, column]
                self._applied[column] += np.count_nonzero(coefficient)
        elif self.kind == "block_reuse":
            gates = decisions < coefficient[:, None]
            starts = np.flatnonzero(positions % self.block_size == 0)
            for column in self.selected:
                origins = np.maximum.accumulate(np.where(positions % self.block_size == 0, np.arange(1, count + 1), 0))
                anchors = np.concatenate(([self._anchor[column]], base[:, column]))[origins]
                output[gates[:, column], column] = anchors[gates[:, column]]
                if len(starts):
                    self._anchor[column] = base[starts[-1], column]
                self._applied[column] += np.count_nonzero(gates[:, column])
        elif self.kind == "cross_correlation":
            first, second = self.pair
            z = coefficient * _normal(base[:, first]) + np.sqrt(np.maximum(0., 1-coefficient**2)) * _normal(base[:, second])
            self._normal_transforms[first] += count
            self._normal_transforms[second] += count
            output[:, second] = np.minimum(ndtr(z), _MAX_UNIFORM)
            output[coefficient == 0, second] = base[coefficient == 0, second]
            self._applied[second] += np.count_nonzero(coefficient)
        self._previous = output[-1].copy()
        self.position += count
        if self.quantize_bits is not None:
            scale = 2. ** self.quantize_bits
            output = np.minimum(np.floor(output * scale) / scale, (scale-1)/scale)
        return output

    @property
    def provenance(self):
        return {"version": "rnglab-v03-uniforms-1", "spec": copy.deepcopy(self.spec),
                "seed_identity": copy.deepcopy(self._seed_identity),
                "dimension": self.dimension, "channels": list(self.channels),
                "selected_columns": list(self.selected), "cross_pair": list(self.pair),
                "clock": "zero-based proposal row",
                "generator": "NumPy PCG64; independent base/decision/init substreams per column",
                "precision": "53-bit dyadic base; float64 Gaussian CDF/inverse approximation",
                "quantize_bits": self.quantize_bits,
                "security": "synthetic dependence model; no physical-entropy or security claim"}

    @property
    def diagnostics(self):
        return {"position": self.position, "applied_by_column": self._applied.tolist(),
                "kind": self.kind, "rho": self.rho,
                "draw_counters": self.draw_counters,
                "stationary_copy_asymptotic_inflation": (1+self.rho)/(1-self.rho) if self.kind == "copy" and self.rho < 1 else None,
                "gaussian_uniform_lag1": float(gaussian_uniform_correlation(self.rho)) if self.kind == "gaussian_ar1" else None,
                "marginal_uniformity": "analytic in ideal arithmetic; see provenance precision",
                "conditional_on_chain_history": "uniformity/independence is not guaranteed under dependence"}

    @property
    def emitted_uniforms(self):
        return self.position * self.dimension

    @property
    def draw_counters(self):
        """Counts of consumed output and latent draws, never entropy estimates.

        Each PCG random() variate consumes one raw 64-bit generator output and
        exposes 53 precision bits. Gaussian variates are inverse-CDF transforms
        of those same uniforms, so they consume no additional random draws.
        """
        return {"emitted_uniforms": self.emitted_uniforms,
                "proposal_rows": self.position,
                "base_uniform_draws": self.emitted_uniforms,
                "decision_uniform_draws": self.emitted_uniforms,
                "initial_uniform_draws": self.dimension * (1+self.lag),
                "normal_random_draws": 0,
                "inverse_normal_transforms": int(self._normal_transforms.sum()),
                "base_uniform_draws_by_column": [self.position] * self.dimension,
                "decision_uniform_draws_by_column": [self.position] * self.dimension,
                "inverse_normal_transforms_by_column": self._normal_transforms.tolist(),
                "raw_pcg64_outputs": 2*self.emitted_uniforms + self.dimension*(1+self.lag),
                "meaning": "algorithmic draws; no entropy/security interpretation"}

    def _fingerprint(self):
        provenance = self.provenance
        # A checkpoint carries all generator states and the original seed identity;
        # its restored identity must not depend on a throwaway constructor seed.
        provenance.pop("seed_identity")
        return hashlib.sha256(json.dumps(provenance, sort_keys=True).encode()).hexdigest()

    def state_dict(self):
        return {"fingerprint": self._fingerprint(), "position": self.position,
                "seed_identity": copy.deepcopy(self._seed_identity),
                "base": [copy.deepcopy(r.bit_generator.state) for r in self._base],
                "decisions": [copy.deepcopy(r.bit_generator.state) for r in self._decisions],
                "initial": [copy.deepcopy(r.bit_generator.state) for r in self._initial],
                "previous": self._previous.tolist(), "latent": self._latent.tolist(),
                "history": self._history.tolist(), "anchor": self._anchor.tolist(),
                "applied": self._applied.tolist(), "normal_transforms": self._normal_transforms.tolist()}

    def load_state_dict(self, state):
        if state["fingerprint"] != self._fingerprint():
            raise ValueError("state belongs to a different specification/dimension/channels")
        for name in ("base", "decisions", "initial"):
            for rng, value in zip(getattr(self, "_" + name), state[name], strict=True):
                rng.bit_generator.state = copy.deepcopy(value)
        self.position = _integer(state["position"], "position")
        for name in ("previous", "latent", "history", "anchor"):
            setattr(self, "_" + name, np.asarray(state[name], dtype=float).copy())
        self._applied = np.asarray(state["applied"], dtype=np.int64).copy()
        self._normal_transforms = np.asarray(state["normal_transforms"], dtype=np.int64).copy()
        self._seed_identity = copy.deepcopy(state["seed_identity"])
