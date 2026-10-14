"""Bounded synthetic digitized-source -> conditioning -> generator pipelines.

These models never claim hardware entropy or cryptographic security. All keys
are public experimental artifacts. A clean 256-bit initial key models adequate
initial material; reproducible PCG-derived keys do not demonstrate real entropy.
Raw-bit time and proposal-row time are distinct. Optional diagnostic probes are
observational copies and never feed a no-reseed generator.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import hmac
import math
import time

import numpy as np
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms


VERSION = "v04-source-pipelines-1"
CHUNK_BITS = 131072


class PipelineExhausted(RuntimeError):
    """A bounded extractor cannot supply the requested output; no fallback."""
    def __init__(self, message, metadata):
        super().__init__(message)
        self.metadata = metadata


def _int(value, name, minimum=0):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be integer >= {minimum}")
    return int(value)


def _seed(seed, *suffix):
    entropy = list(seed) if isinstance(seed, (list, tuple, np.ndarray)) else [int(seed)]
    return np.random.SeedSequence([*entropy, *suffix])


def bits_to_uniform64(bits, rows):
    """Consecutive MSB-first 64-bit words, top53 bits, exact [0,1) dyadic grid."""
    bits = np.asarray(bits, dtype=np.uint8)
    if len(bits) != rows*128 or np.any(bits > 1):
        raise ValueError("exactly rows*2*64 binary digits required")
    packed = np.packbits(bits).tobytes()
    words = np.frombuffer(packed, dtype=">u8")
    return ((words >> 11).astype(np.float64)*2.**-53).reshape(rows, 2)


def bytes_to_uniform64(data, rows):
    if len(data) != rows*16:
        raise ValueError("exactly16 bytes per two-uniform row required")
    words = np.frombuffer(data, dtype=">u8")
    return ((words >> 11).astype(np.float64)*2.**-53).reshape(rows, 2)


def _features(bits):
    n = len(bits);ones = int(bits.sum())
    variance = ones/n*(1-ones/n) if n else None
    if n > 1 and variance and np.std(bits[:-1]) > 0 and np.std(bits[1:]) > 0:
        correlation = float(np.corrcoef(bits[:-1], bits[1:])[0, 1])
    else:
        correlation = None
    changes = np.flatnonzero(bits[1:] != bits[:-1])+1 if n > 1 else np.empty(0, int)
    longest = int(np.diff(np.r_[0, changes, n]).max()) if n else 0
    return {"bits": n, "ones_fraction": ones/n if n else None,
            "lag1": correlation, "transition_fraction": len(changes)/(n-1) if n > 1 else None,
            "monobit_score": abs(2*ones-n)/math.sqrt(n) if n else None,
            "longest_run": longest, "constant": bool(n and (ones == 0 or ones == n))}


def _planned_windows(faults, extent, width, maximum):
    windows = set()
    for fault in faults:
        points = [fault.get("onset_bit", 0)]
        if fault.get("ramp_bits"):
            points.append(fault.get("onset_bit", 0)+fault["ramp_bits"])
        if fault.get("end_bit") is not None:
            points.extend([fault["end_bit"], fault["end_bit"]+fault.get("recovery_bits", 0)])
        for point in points:
            windows.update(max(0, point//width+d) for d in (-1, 0, 1))
    windows.update(range(min(4, extent//width)))
    windows = {w for w in windows if (w+1)*width <= extent}
    if len(windows) > maximum:
        windows = set(sorted(windows)[:maximum])
    remaining = maximum-len(windows)
    if remaining:
        windows.update(int(v) for v in np.linspace(0, max(0, extent//width-1), remaining))
    return sorted(windows)[:maximum]


class BitAudit:
    def __init__(self, starts=(), width=512):
        self.position = self.ones = self.same = self.pairs = 0
        self.first = self.last = None
        self.hash = hashlib.sha256()
        self.width = width
        self.pending = {int(start): [] for start in starts}
        self.windows = []

    def update(self, bits):
        n = len(bits)
        if not n:
            return
        start, end = self.position, self.position+n
        self.hash.update(np.asarray(bits, dtype=np.uint8).tobytes())
        self.ones += int(bits.sum())
        self.same += int(np.count_nonzero(bits[1:] == bits[:-1]))
        self.pairs += n-1
        if self.last is not None:
            self.same += int(self.last == bits[0]);self.pairs += 1
        else:
            self.first = int(bits[0])
        self.last = int(bits[-1])
        for target in list(self.pending):
            a, b = max(start, target), min(end, target+self.width)
            if a < b:
                self.pending[target].append(bits[a-start:b-start].copy())
                if b == target+self.width:
                    values = np.concatenate(self.pending.pop(target))
                    self.windows.append({"raw_bit_start": target, "raw_bit_end": b, **_features(values)})
        self.position = end

    def summary(self):
        n = self.position
        mean = self.ones/n if n else None
        return {"bits_observed": n, "sha256_uint8_bits": self.hash.hexdigest(),
                "ones_fraction": mean, "same_adjacent_fraction": self.same/self.pairs if self.pairs else None,
                "constant": bool(n and self.ones in (0, n)),
                "windows": sorted(self.windows, key=lambda w:w["raw_bit_start"]),
                "diagnostic_window_bits": self.width,
                "scope": "descriptive source diagnostics; no entropy-rate or security inference"}


def _faults(spec):
    values = spec.get("faults")
    if values is None:
        values = [] if spec.get("kind", "iid") in {"iid", "clean"} else [spec]
    aliases = {"markov": "copy", "autocorrelation": "copy", "frozen_block": "freeze",
               "intermittent_freeze": "freeze", "ramp_bias": "bias"}
    result = []
    for value in values:
        fault = dict(value);fault["kind"] = aliases.get(fault["kind"], fault["kind"])
        if fault["kind"] not in {"bias", "copy", "freeze", "periodic"}:
            raise ValueError("unknown raw-bit fault")
        for key in ("onset_bit", "ramp_bits", "recovery_bits"):
            fault[key] = _int(fault.get(key, 0), key)
        if fault.get("end_bit") is not None:
            fault["end_bit"] = _int(fault["end_bit"], "end_bit")
            if fault["end_bit"] <= fault["onset_bit"]:
                raise ValueError("end_bit must exceed onset_bit")
        if fault.get("period_bits"):
            period = _int(fault["period_bits"], "period_bits", 1)
            active = _int(fault.get("active_bits", period), "active_bits", 1)
            if active > period:
                raise ValueError("active_bits exceeds period_bits")
            fault["active_bits"] = active
        if fault["kind"] == "bias" and ("p1" in fault or "probability_one" in fault):
            p = float(fault.get("p1", fault.get("probability_one")))
            if not 0 <= p <= 1:
                raise ValueError("p1 outside[0,1]")
            fault["strength"] = 2*abs(p-.5);fault["value"] = int(p >= .5)
        else:
            fault["strength"] = float(fault.get("strength", fault.get("rho", 1.)))
        if not np.isfinite(fault["strength"]) or not 0 <= fault["strength"] <= 1:
            raise ValueError("fault strength outside[0,1]")
        if fault.get("value", 0) not in (0, 1):
            raise ValueError("frozen/bias value must be0 or1")
        pattern = fault.get("pattern", "0011")
        if not pattern or set(pattern)-{"0", "1"}:
            raise ValueError("periodic pattern must be binary")
        fault["pattern"] = pattern
        result.append(fault)
    return result


def _envelope(positions, fault):
    onset, ramp = fault["onset_bit"], fault["ramp_bits"]
    envelope = np.clip((positions-onset)/ramp, 0, 1) if ramp else (positions >= onset).astype(float)
    if fault.get("end_bit") is not None:
        end, recovery = fault["end_bit"], fault["recovery_bits"]
        envelope *= np.clip((end+recovery-positions)/recovery, 0, 1) if recovery else (positions < end)
    if fault.get("period_bits"):
        envelope *= (positions-onset)%fault["period_bits"] < fault["active_bits"]
    return envelope


class RawSource:
    def __init__(self, spec, seed, maximum_bits, diagnostic_extent=None, max_windows=32):
        self.faults = _faults(spec);self.spec = dict(spec)
        self.maximum_bits = maximum_bits;self.position = 0
        self.base = np.random.Generator(np.random.PCG64(_seed(seed, 401)))
        self.decisions = [np.random.Generator(np.random.PCG64(_seed(seed, 402, i))) for i in range(len(self.faults))]
        self.previous = [int(np.random.default_rng(_seed(seed, 403, i)).integers(2)) for i in range(len(self.faults))]
        self.pending = np.empty(0, np.uint8)
        starts = _planned_windows(self.faults, diagnostic_extent or maximum_bits, 512, max_windows)
        self.audit = BitAudit([w*512 for w in starts])
        self.base_raw64_words = 0
        self.decision_uniform_draws = 0

    def read(self, count):
        count = _int(count, "raw count")
        if self.position+count > self.maximum_bits:
            raise PipelineExhausted("bounded raw source exhausted", self.metadata())
        output = np.empty(count, np.uint8)
        for offset in range(0, count, CHUNK_BITS):
            take = min(CHUNK_BITS, count-offset)
            missing = max(0, take-len(self.pending))
            if missing:
                words = self.base.bit_generator.random_raw((missing+63)//64)
                self.base_raw64_words += len(words)
                self.pending = np.concatenate((self.pending, np.unpackbits(words.astype(">u8").view(np.uint8))))
            bits = self.pending[:take].copy();self.pending = self.pending[take:]
            positions = np.arange(self.position, self.position+take)
            for index, fault in enumerate(self.faults):
                gates = self.decisions[index].random(take) < fault["strength"]*_envelope(positions, fault)
                self.decision_uniform_draws += take
                if fault["kind"] == "copy":
                    origins = np.maximum.accumulate(np.where(~gates, np.arange(1, take+1), 0))
                    bits = np.concatenate(([self.previous[index]], bits))[origins].astype(np.uint8)
                elif fault["kind"] == "periodic":
                    pattern = np.array(list(fault["pattern"]), dtype=np.uint8)
                    bits[gates] = pattern[positions[gates]%len(pattern)]
                else:
                    bits[gates] = fault.get("value", 1 if fault["kind"] == "bias" else 0)
                self.previous[index] = int(bits[-1])
            output[offset:offset+take] = bits
            self.audit.update(bits);self.position += take
        return output

    def metadata(self):
        return {"model": self.spec, "normalized_faults": self.faults,
                "raw_bits_consumed": self.position, "base_raw64_words": self.base_raw64_words,
                "raw_fault_initial_bernoulli_draws": len(self.faults),
                "decision_uniform_draws": self.decision_uniform_draws,
                "audit": self.audit.summary(), "audit_layer": "digitized_source_bits",
                "classification": "synthetic digitized source; reproducible public PCG seed; no measured physical entropy"}


class MultiSource:
    """XOR inputs with explicit independence or shared-input relationships."""
    def __init__(self, specs, seed, maximum_bits, relationship, correlation, extent, max_windows):
        per_source_windows = max_windows//len(specs) if specs else 0
        self.sources = [RawSource(s, [*seed, 911, i], maximum_bits, extent, per_source_windows) for i,s in enumerate(specs)]
        self.relationship = relationship;self.correlation = correlation
        self.rng = np.random.default_rng(_seed(seed, 912))
        self.position = 0;self.maximum_bits = maximum_bits
        self.audit = BitAudit();self.effective = [BitAudit() for _ in specs]
        self.equal = self.observed = self.both_one = 0
        if len(specs) < 2 or relationship not in {"independent", "identical", "complement", "correlated_copy"}:
            raise ValueError("XOR requires>=2 inputs and explicit valid relationship")
        if not 0 <= correlation <= 1:
            raise ValueError("source correlation outside[0,1]")

    def read(self, count):
        values = [s.read(count) for s in self.sources]
        if self.relationship == "identical":
            values[1] = values[0].copy()
        elif self.relationship == "complement":
            values[1] = 1-values[0]
        elif self.relationship == "correlated_copy":
            mask = self.rng.random(count) < self.correlation
            values[1][mask] = values[0][mask]
        for audit, value in zip(self.effective, values):
            audit.update(value)
        self.equal += int(np.count_nonzero(values[0] == values[1]));self.observed += count
        self.both_one += int(np.count_nonzero(values[0] & values[1]))
        output = np.bitwise_xor.reduce(values)
        self.audit.update(output);self.position += count
        return output

    def metadata(self):
        n = self.observed
        means = [audit.ones/n if n else None for audit in self.effective[:2]]
        denominator = math.sqrt(means[0]*(1-means[0])*means[1]*(1-means[1])) if n else 0
        cross = (self.both_one/n-means[0]*means[1])/denominator if denominator else None
        return {"raw_bits_consumed": self.position, "raw_input_bits_total": sum(s.position for s in self.sources),
                "relationship": self.relationship, "copy_probability": self.correlation,
                "input_equal_fraction": self.equal/self.observed if self.observed else None,
                "raw_input_cross_correlation": cross,
                "pre_relationship_inputs": [s.metadata() for s in self.sources],
                "effective_inputs": [a.summary() for a in self.effective], "audit": self.audit.summary(),
                "audit_layer": "XOR_conditioning_output_not_raw_inputs",
                "effective_input_layer": "digitized_bits_after_explicit_synthetic_source_relationship_before_XOR",
                "assumption": "XOR needs genuinely independent input quality; identical inputs cancel. A fresh independent fair input masks arbitrary independent other inputs."}


class Conditioner:
    def __init__(self, source, kind):
        if kind not in {"none", "von_neumann", "sha256"}:
            raise ValueError("conditioner must be none/von_neumann/sha256")
        self.source, self.kind = source, kind
        self.pending = np.empty(0, np.uint8);self.audit = BitAudit();self.hash_blocks = 0
        self.returned_bits = 0

    def read(self, count):
        out = np.empty(count, np.uint8);cursor = 0
        while cursor < count:
            if len(self.pending):
                take = min(len(self.pending), count-cursor)
                out[cursor:cursor+take] = self.pending[:take]
                self.pending = self.pending[take:];cursor += take
                continue
            missing = count-cursor
            if self.kind == "none":
                transformed = self.source.read(min(CHUNK_BITS, missing))
            elif self.kind == "sha256":
                blocks = min(CHUNK_BITS//512, (missing+255)//256)
                raw = self.source.read(blocks*512).reshape(blocks, 512)
                packed = np.packbits(raw, axis=1)
                # No counter: identical source blocks MUST yield identical hashes.
                digests = b"".join(hashlib.sha256(b"rnglab-v04/condition/"+row.tobytes()).digest() for row in packed)
                transformed = np.unpackbits(np.frombuffer(digests, np.uint8));self.hash_blocks += blocks
            else:
                available = self.source.maximum_bits-self.source.position
                take = min(CHUNK_BITS, max(512, missing*4), available//2*2)
                if take < 2:
                    raise PipelineExhausted("von Neumann stalled/exhausted; no fabricated output", self.metadata())
                raw = self.source.read(take).reshape(-1, 2)
                transformed = raw[raw[:,0] != raw[:,1],0]
            self.audit.update(transformed);self.pending = transformed
        self.returned_bits += count
        return out

    def metadata(self):
        return {"kind": self.kind, "returned_bits": self.returned_bits,
                "generated_conditioned_bits": self.audit.position, "buffered_bits": len(self.pending),
                "raw_bits_consumed": self.source.position, "hash_blocks": self.hash_blocks,
                "audit": self.audit.summary(),
                "assumptions": {"none":"no conditioning",
                    "von_neumann":"fair independent extracted bits only under iid stationary Bernoulli assumptions; dependence/drift can invalidate them; freezes yield no output",
                    "sha256":"deterministic whitening, not a proven seeded extractor or entropy guarantee; identical blocks repeat"}[self.kind]}


def _uniform_metadata(u):
    n = len(u)
    if not n:
        return {"rows": 0}
    correlation = np.corrcoef(u.T)[0,1] if n > 1 and np.all(u.std(axis=0)>0) else None
    lag = [float(np.corrcoef(u[:-1,j],u[1:,j])[0,1]) if n>2 and np.std(u[:-1,j])>0 and np.std(u[1:,j])>0 else None for j in range(2)]
    return {"rows": n, "mean": u.mean(axis=0).tolist(), "variance": u.var(axis=0).tolist(),
            "lag1": lag, "proposal_acceptance_correlation": float(correlation) if correlation is not None else None,
            "adjacent_equal_fraction": (u[1:] == u[:-1]).mean(axis=0).tolist() if n>1 else None,
            "histogram16": [np.histogram(u[:,j], np.linspace(0,1,17))[0].tolist() for j in range(2)],
            "scope":"descriptive marginal/joint tests; numerical plausibility is not cryptographic security"}


def generate_pipeline(spec, seed, rows):
    """Return Nx2 uniforms plus explicit layer/probe/key/draw provenance.

    generator: direct, pcg64, philox, chacha20. conditioner: none, sha256,
    von_neumann, xor (requires xor_sources). initialization: clean_independent
    (default) or raw_prefix. Periodic reseed defaults to HMAC mixing old state;
    replace_key_reset_counter is an explicitly unsafe separate processing fault.
    """
    tick = time.perf_counter();spec = dict(spec);seed = list(seed) if isinstance(seed,(list,tuple,np.ndarray)) else [int(seed)]
    rows = _int(rows, "rows")
    if rows > spec.get("max_rows", 2_000_000):
        raise ValueError("bounded max_rows exceeded")
    generator = spec.get("generator", "direct")
    if generator not in {"direct", "pcg64", "philox", "chacha20"}:
        raise ValueError("unknown generator")
    conditioner = spec.get("conditioner", "none")
    source_spec = spec.get("source", {"kind":"iid"})
    interval = spec.get("reseed_interval_rows")
    if interval is not None:
        interval = _int(interval, "reseed_interval_rows", 1)
        if generator == "direct":
            raise ValueError("direct generation has no reseed state")
    initial = spec.get("initialization", "clean_independent")
    if initial not in {"clean_independent", "raw_prefix"}:
        raise ValueError("unknown initialization source")
    required = rows*128
    expected_raw = required*(4 if conditioner == "von_neumann" else 2 if conditioner == "sha256" else 1) if generator == "direct" else 256*((rows-1)//interval+2 if rows and interval else 1)
    max_raw = _int(spec.get("max_raw_bits", max(65536, required*8) if generator == "direct" and conditioner == "von_neumann" else max(65536, expected_raw*2)), "max_raw_bits", 256)
    max_windows = _int(spec.get("diagnostic_max_windows", 32), "diagnostic_max_windows", 0)
    if max_windows > 32:
        raise ValueError("at most32 diagnostic windows per source")
    def factory(extent, cap):
        if conditioner == "xor" or spec.get("xor_sources"):
            sources = spec.get("xor_sources", [source_spec, {"kind":"iid"}])
            return MultiSource(sources, seed, cap, spec.get("relationship", "independent"),
                float(spec.get("source_correlation", 0.)), extent, max_windows)
        return RawSource(source_spec, seed, cap, extent, max_windows)
    raw = factory(expected_raw, max_raw)
    conditional = Conditioner(raw, "none" if conditioner == "xor" else conditioner)
    metadata = {"version": VERSION, "spec": spec, "seed_identity": seed,
        "classification":"synthetic source-to-application experiment; public deterministic seeds/keys, never production secrets",
        "layers":["modeled digitized source bits", "optional XOR/extraction/hash conditioning", "optional deterministic generator", "64-bit word to53-bit uniform", "adjacent words assigned proposal/acceptance"],
        "fault_layer":"digitized source bits; replace_key_reset_counter if selected is separately a generator-processing fault",
        "uniform_mapping":"MSB-first64-bit word, discard11LSBs, divide top53 by2^53; exact[0,1)",
        "clock":"raw bits and proposal rows are distinct; variable extraction and sparse reseeding do not imply equal clocks",
        "adequate_seed_scope":"clean material models the adequate-initial-entropy assumption; public PCG seed realizes zero secret entropy conditional on published experiment identity"}
    keys=[];reseed_events=[]
    if generator == "direct":
        try:
            bits = conditional.read(required)
            u = bits_to_uniform64(bits, rows)
        except PipelineExhausted as error:
            error.metadata = {**metadata,"status":"conditioner_exhausted", "raw_source":raw.metadata(),"conditioner":conditional.metadata()}
            raise
        key = None;initial_consumed = 0
    else:
        before = raw.position
        if initial == "raw_prefix":
            material = np.packbits(conditional.read(256)).tobytes()
        else:
            words = np.random.PCG64(_seed(seed, 12021)).random_raw(4)
            material = words.astype(">u8").tobytes()
        initial_consumed = raw.position-before
        key = material
        policy = spec.get("reseed_policy", "mix_old_state")
        if policy not in {"mix_old_state", "replace_key_reset_counter"}:
            raise ValueError("unknown reseed policy")
        u = np.empty((rows,2));epoch = 0
        for start in range(0, max(rows,1), interval or max(rows,1)):
            count = min(interval or rows, rows-start)
            if start and interval:
                offset = raw.position;fresh = np.packbits(conditional.read(256)).tobytes();epoch += 1
                key = (hmac.new(key, b"rnglab-v04/reseed/"+epoch.to_bytes(8,"big")+fresh, hashlib.sha256).digest()
                       if policy == "mix_old_state" else fresh)
                reseed_events.append({"proposal_row":start,"raw_bit_start":offset,"raw_bit_end":raw.position,
                    "conditioned_seed_bits":256,"fresh_material_sha256":hashlib.sha256(fresh).hexdigest(),"policy":policy})
            nonce_epoch = 0 if policy == "replace_key_reset_counter" else epoch
            keys.append({"epoch":epoch,"proposal_row":start,"key_hex_public_experimental":key.hex(),
                         "key_sha256":hashlib.sha256(key).hexdigest(),"nonce_epoch":nonce_epoch})
            if generator == "chacha20":
                cipher = Cipher(algorithms.ChaCha20(key, bytes(8)+nonce_epoch.to_bytes(8,"little")),mode=None).encryptor()
                u[start:start+count] = bytes_to_uniform64(cipher.update(bytes(count*16)),count)
            else:
                entropy = np.frombuffer(key,dtype=">u4").astype(np.uint32).tolist()
                bit_generator = np.random.PCG64(entropy) if generator == "pcg64" else np.random.Philox(entropy)
                u[start:start+count] = np.random.Generator(bit_generator).random((count,2))
    metadata.update(status="complete",raw_source=raw.metadata(),conditioner=conditional.metadata(),
        generator=generator,output_rows=rows,emitted_uniforms=rows*2,
        raw_source_bits_consumed_for_generation=raw.position,initial_raw_bits_consumed=initial_consumed,
        post_initial_raw_bits_consumed=raw.position-initial_consumed,
        clean_initial_material_bits=256 if generator != "direct" and initial == "clean_independent" else 0,
        output_raw64_words=rows*2,discarded_word_bits=rows*2*11,
        generator_initialization=initial if generator != "direct" else None,
        initial_key_sha256=keys[0]["key_sha256"] if keys else None,
        experimental_public_keys=keys,reseed_events=reseed_events,
        repeated_key_nonce_pairs=len(keys)-len({(k["key_sha256"],k["nonce_epoch"]) for k in keys}),
        uniform_sha256=hashlib.sha256(u.astype("<f8").tobytes()).hexdigest(),uniform_diagnostics=_uniform_metadata(u))
    metadata["raw_input_bits_total_for_generation"] = sum(s.position for s in raw.sources) if isinstance(raw,MultiSource) else raw.position
    metadata["raw_input_diagnostics"] = [a.summary() for a in raw.effective] if isinstance(raw,MultiSource) else [raw.audit.summary()]
    metadata["raw_input_diagnostic_scope"] = "before conditioning; effective source inputs after any explicitly modeled shared-input relationship"
    metadata["conditioning_stages"] = (["xor"] if isinstance(raw,MultiSource) else [])+([conditional.kind] if conditional.kind != "none" else [])
    metadata["clean_initial_raw64_words"] = 4 if generator != "direct" and initial == "clean_independent" else 0
    if generator != "direct" and spec.get("diagnostic_probe_bits",150000):
        probe_bits = _int(spec.get("diagnostic_probe_bits",150000),"diagnostic_probe_bits")
        probe = factory(probe_bits,probe_bits)
        for start in range(0,probe_bits,CHUNK_BITS):
            probe.read(min(CHUNK_BITS,probe_bits-start))
        metadata["observational_raw_probe"]={"feeds_generator":False,"bits":probe_bits,"source":probe.metadata(),
            "scope":"separate copy of hypothetical raw sampling clock; never counted as consumed generator seed/reseed bits"}
    metadata["elapsed_seconds"] = time.perf_counter()-tick
    return u,metadata
