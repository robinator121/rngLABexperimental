"""Bit sources, uniform adapters, and a demonstration reseeding construction."""
from typing import Protocol
import hashlib
import hmac
import numpy as np
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms
from rnglab.degradation import Fault, FaultSpec

BITS_PER_UNIFORM = 32
WEIGHTS = np.exp2(-np.arange(1, BITS_PER_UNIFORM + 1, dtype=float))


class BitSource(Protocol):
    def read_bits(self, count: int) -> np.ndarray: ...


class RandomnessProvider(Protocol):
    def sample(self, count: int, dimension: int) -> tuple[np.ndarray, np.ndarray]: ...


class SyntheticSource:
    """Reproducible Bernoulli(1/2) model, NOT a physical TRNG."""
    def __init__(self, seed, faults=()):
        children = np.random.SeedSequence(seed).spawn(1 + len(faults))
        self.rng = np.random.Generator(np.random.PCG64(children[0]))
        self.faults = [Fault(FaultSpec(**f), child) for f, child in zip(faults, children[1:])]
        self.position = 0

    def read_bits(self, count):
        if count < 0:
            raise ValueError("count must be nonnegative")
        out = (self.rng.random(count) < 0.5).astype(np.uint8)
        for fault in self.faults:
            out = fault.apply(out, self.position)
        self.position += count
        return out


class PCGSource:
    """Standard seeded NumPy PCG64 with buffered raw output."""
    def __init__(self, seed):
        self.rng = np.random.PCG64(seed)
        self.pending = np.empty(0, dtype=np.uint8)

    def read_bits(self, count):
        if count < 0:
            raise ValueError("count must be nonnegative")
        missing = max(0, count - len(self.pending))
        if missing:
            raw = self.rng.random_raw((missing + 63) // 64).astype("<u8").tobytes()
            self.pending = np.concatenate((self.pending, np.unpackbits(np.frombuffer(raw, dtype=np.uint8))))
        result, self.pending = self.pending[:count].copy(), self.pending[count:]
        return result


class ChaChaSource:
    """Established ChaCha20 keystream; deterministic experimental seeds only."""
    def __init__(self, key: bytes, epoch=0):
        if len(key) != 32:
            raise ValueError("ChaCha20 requires a 32-byte key")
        # Zero initial block counter; epoch occupies final eight nonce bytes.
        # Works with both historical 64/64 and newer 32/96 API interpretations.
        self.cipher = Cipher(algorithms.ChaCha20(key, bytes(8) + epoch.to_bytes(8, "little")), mode=None).encryptor()
        self.pending = np.empty(0, dtype=np.uint8)
        self.bytes_emitted = 0

    def read_bits(self, count):
        if count < 0:
            raise ValueError("count must be nonnegative")
        missing = max(0, count - len(self.pending))
        if missing:
            size = (missing + 7) // 8
            if self.bytes_emitted + size > 2**38:
                raise ValueError("ChaCha20 per-epoch byte limit exceeded")
            raw = self.cipher.update(bytes(size))
            self.bytes_emitted += size
            self.pending = np.concatenate((self.pending, np.unpackbits(np.frombuffer(raw, dtype=np.uint8))))
        result, self.pending = self.pending[:count].copy(), self.pending[count:]
        return result


def bits_to_uniform(bits, count, dimension):
    """MSB first, exact dyadic grid {0,...,2^32-1}/2^32; never returns 1."""
    return (bits.reshape(count, dimension, BITS_PER_UNIFORM) * WEIGHTS).sum(axis=2)


class DirectProvider:
    def __init__(self, source):
        self.source = source
        self.position = 0
        self.events = []

    def sample(self, count, dimension):
        if count < 0 or dimension < 1:
            raise ValueError("invalid sample shape")
        bits = self.source.read_bits(count * dimension * BITS_PER_UNIFORM)
        self.position += count
        return bits_to_uniform(bits, count, dimension), bits


def entropy_bytes(source):
    return np.packbits(source.read_bits(256)).tobytes()


class HybridProvider:
    """Demonstration HMAC state mixing + ChaCha20, NOT a validated DRBG."""
    def __init__(self, initial_source, reseed_source=None, interval=None):
        if interval is not None and (not isinstance(interval, int) or interval < 1 or reseed_source is None):
            raise ValueError("reseed interval must be positive and have a source")
        self.key = entropy_bytes(initial_source)
        self.stream = ChaChaSource(self.key)
        self.reseed_source = reseed_source
        self.interval = interval
        self.position = 0
        self.epoch = 0
        self.events = []

    def _reseed(self):
        entropy_offset = self.reseed_source.position
        fresh = entropy_bytes(self.reseed_source)
        self.epoch += 1
        self.key = hmac.new(self.key, b"rnglab-v0.1/reseed/" + self.epoch.to_bytes(8, "big") + fresh, hashlib.sha256).digest()
        self.stream = ChaChaSource(self.key, self.epoch)
        self.events.append({"sample_index": self.position, "epoch": self.epoch, "source_bit_offset": entropy_offset})

    def sample(self, count, dimension):
        if count < 0 or dimension < 1:
            raise ValueError("invalid sample shape")
        remaining, parts = count, []
        while remaining:
            if self.interval and self.position and self.position % self.interval == 0:
                self._reseed()
            take = remaining if self.interval is None else min(remaining, self.interval - self.position % self.interval)
            parts.append(self.stream.read_bits(take * dimension * BITS_PER_UNIFORM))
            self.position += take
            remaining -= take
        bits = np.concatenate(parts) if parts else np.empty(0, dtype=np.uint8)
        return bits_to_uniform(bits, count, dimension), bits


def make_provider(config, seed):
    """Keep clean initialization identical between hybrid/control configurations."""
    kind = config["provider"]
    if kind == "synthetic":
        return DirectProvider(SyntheticSource(seed, config.get("faults", [])))
    if kind == "prng":
        return DirectProvider(PCGSource(seed))
    if kind in {"csprng", "hybrid"}:
        initial = SyntheticSource(seed, config.get("initial_faults", []))
        if kind == "csprng":
            return HybridProvider(initial)
        reseed = SyntheticSource([*seed, 9271], config.get("faults", []))
        return HybridProvider(initial, reseed, config["reseed_interval"])
    raise ValueError(f"Unknown provider: {kind}")
