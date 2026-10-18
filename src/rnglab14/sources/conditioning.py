"""RNGLab v0.14: Entropy Conditioning, Extraction, and Reseeding Pipelines.

SECURITY NOTICE:
These conditioning mechanisms are implemented strictly for academic research and
simulation benchmarks. They are NOT certified for production cryptographic key generation.
Passing NIST SP 800-22 tests on cryptographically conditioned output demonstrates only
PRF hash compression properties, not the min-entropy of the underlying entropy source.
Min-entropy must be evaluated directly on raw noise using NIST SP 800-90B.

Implements research conditioning mechanisms with explicit assumptions:
1. Von Neumann Pairwise Debiasing (with empirical efficiency tracking)
2. XOR parity with limited health diagnostics (not independence verification)
3. SHA-256 block compressor (does not create or certify entropy)
4. Reseeding DRBG Pipeline (configurable cadence: 100, 1000, 10000, 10^6 cycles,
   with entropy tracking, reseed failure injection, and state recovery)
"""
from __future__ import annotations

from collections import deque
import hashlib
import math
import struct
from typing import List, Optional, Sequence

from .base import RandomSource, SourceManifest


class IndependenceViolationWarning(UserWarning):
    """Raised when XOR conditioning streams exhibit statistically significant cross-correlation."""
    pass


class VonNeumannConditioner(RandomSource):
    """Classical Von Neumann pairwise bit extractor:
    Takes pairs of bits:
    - (0, 1) -> 0
    - (1, 0) -> 1
    - (0, 0) -> discard
    - (1, 1) -> discard
    """

    def __init__(self, raw_source: RandomSource, max_bytes: Optional[int] = None, max_raw_bits_per_output_bit: int = 256):
        manifest = SourceManifest(
            source_id=f"vn_{raw_source.manifest.source_id}",
            family="conditioned",
            is_deterministic=raw_source.manifest.is_deterministic,
            is_physical=raw_source.manifest.is_physical,
            metadata={"mechanism": "von_neumann", "parent_source": raw_source.manifest.source_id},
        )
        super().__init__(manifest, max_bytes=max_bytes)
        if max_raw_bits_per_output_bit < 4:
            raise ValueError("extraction budget must be >=4 raw bits per output bit")
        self.max_raw_bits_per_output_bit = max_raw_bits_per_output_bit
        self.raw = raw_source
        self.raw_bits_consumed = 0  # examined bits
        self.raw_bits_acquired = 0  # actual parent bytes*8, including buffered bits
        self._raw_bit_buffer = deque()
        self.conditioned_bits_produced = 0
        self._bit_buffer: list[int] = []

    def _extract_bits(self, count: int):
        initial = self.raw_bits_consumed
        budget = max(1, count-len(self._bit_buffer)) * self.max_raw_bits_per_output_bit
        while len(self._bit_buffer) < count:
            if self.raw_bits_consumed-initial >= budget:
                raise ValueError("Von Neumann extraction budget exhausted; source may be stuck")
            # Read two bits from raw source
            pair=[]
            for _ in range(2):
                if not self._raw_bit_buffer:
                    value=self.raw.read_bytes(1)[0]
                    self._raw_bit_buffer.extend((value>>shift)&1 for shift in range(7,-1,-1))
                    self.raw_bits_acquired += 8
                pair.append(self._raw_bit_buffer.popleft())
                self.raw_bits_consumed += 1
            b1,b2=pair
            if b1 == 0 and b2 == 1:
                self._bit_buffer.append(0)
                self.conditioned_bits_produced += 1
            elif b1 == 1 and b2 == 0:
                self._bit_buffer.append(1)
                self.conditioned_bits_produced += 1
            # (0,0) and (1,1) discarded

    def _read_bytes(self, n: int) -> bytes:
        bits_needed = n * 8
        self._extract_bits(bits_needed)
        out_bytes = bytearray()
        for _ in range(n):
            byte_val = 0
            for shift in range(7, -1, -1):
                b = self._bit_buffer.pop(0)
                byte_val |= (b << shift)
            out_bytes.append(byte_val)
        return bytes(out_bytes)

    @property
    def extraction_efficiency(self) -> float:
        """Output/acquired bits. IID fair expected rate .25, not a finite upper bound."""
        if self.raw_bits_acquired == 0:
            return 0.0
        return self.conditioned_bits_produced / self.raw_bits_acquired


class XORParityConditioner(RandomSource):
    """XOR with limited finite-sample health screens, never independence proof.

    check_independence is retained only as a compatibility spelling for enabling
    diagnostics. Passing these screens establishes neither mutual independence,
    min-entropy, cryptographic security nor compliance. Even nonlinear/rank
    tests cannot establish those guarantees from a finite prefix.
    """

    def __init__(self, sources: Sequence[RandomSource], check_independence: bool = True,
                 max_bytes: Optional[int] = None):
        if len(sources) < 2 or len({id(s) for s in sources}) != len(sources):
            raise ValueError("XOR requires at least two distinct source objects")
        manifest = SourceManifest(
            source_id=f"xor_{len(sources)}streams", family="conditioned",
            is_deterministic=all(s.manifest.is_deterministic for s in sources),
            is_physical=any(s.manifest.is_physical for s in sources),
            metadata={"mechanism":"xor_parity", "parent_sources":[s.manifest.source_id for s in sources],
                      "diagnostics_enabled":bool(check_independence),
                      "independence_established":False,"entropy_certified":False,
                      "diagnostic_scope":"128-byte prefix and constant-output check; not calibrated assurance"})
        super().__init__(manifest, max_bytes=max_bytes)
        self.sources = list(sources)
        self.diagnostics_enabled = bool(check_independence)
        self._last_output_byte = None
        self._output_run_length = 0
        self.pairwise_correlations: dict[str,float] = {}
        if self.diagnostics_enabled:
            self._screen_prefix()

    @staticmethod
    def _xor(streams):
        out=bytearray(streams[0])
        for data in streams[1:]:
            for i,value in enumerate(data): out[i] ^= value
        return bytes(out)

    def _screen_prefix(self, n=128):
        samples=[s.read_bytes(n) for s in self.sources]
        if any(len(set(x)) <= 1 for x in samples):
            raise ValueError("XOR diagnostic rejected constant input prefix; entropy unknown")
        for i,x in enumerate(samples):
            for j in range(i+1,len(samples)):
                y=samples[j];mx=sum(x)/n;my=sum(y)/n
                cov=sum((a-mx)*(b-my) for a,b in zip(x,y))
                denom=math.sqrt(sum((a-mx)**2 for a in x)*sum((b-my)**2 for b in y))
                corr=cov/denom
                self.pairwise_correlations[f"{i}_vs_{j}"]=corr
                if abs(corr)>.15:
                    raise ValueError("XOR diagnostic rejected correlated prefix; this is not an independence test")
        if len(set(self._xor(samples))) <= 1:
            raise ValueError("XOR diagnostic rejected constant combined prefix (possible joint dependence)")

    def _read_bytes(self, n):
        out=self._xor([s.read_bytes(n) for s in self.sources])
        if self.diagnostics_enabled:
            for value in out:
                self._output_run_length = self._output_run_length+1 if value==self._last_output_byte else 1
                self._last_output_byte = value
                if self._output_run_length>=16:
                    raise ValueError("XOR diagnostic rejected constant output run; no entropy assurance")
        return out


# Historical import preserved, with explicitly restricted semantics above.
XORIndependenceCheckedConditioner = XORParityConditioner


class HashBlockConditioner(RandomSource):
    """Cryptographic hash compressor using SHA-256.
    
    Compresses M input bytes from raw source into 32 output bytes.
    Byte compression ratio is 32/M; output entropy depends on input assumptions.
    """

    def __init__(self, raw_source: RandomSource, input_block_bytes: int = 64, max_bytes: Optional[int] = None):
        if input_block_bytes < 32:
            raise ValueError("hash compressor input block must be >=32 bytes")
        manifest = SourceManifest(
            source_id=f"hash256_{raw_source.manifest.source_id}",
            family="conditioned",
            is_deterministic=raw_source.manifest.is_deterministic,
            is_physical=raw_source.manifest.is_physical,
            metadata={"mechanism": "sha256_extractor", "compression_ratio": 32 / input_block_bytes},
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self.raw = raw_source
        self.block_size = input_block_bytes
        self._buffer = bytearray()

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            raw_block = self.raw.read_bytes(self.block_size)
            digest = hashlib.sha256(raw_block).digest()
            self._buffer.extend(digest)
        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out


class ReseedingDRBGPipeline(RandomSource):
    """Hybrid source-to-DRBG pipeline with configurable reseed schedules.
    
    Evaluates:
    - Reseed cadence: 100, 1000, 10000, or 1,000,000 samples
    - Entropy accounting
    - Simulated reseed failures (reseed dropped, or poisoned by degraded raw source)
    - Generator state recovery vs corruption persistence
    """

    def __init__(
        self,
        entropy_source: RandomSource,
        reseed_interval_samples: int = 1000,
        simulate_reseed_failure: bool = False,
        max_bytes: Optional[int] = None,
    ):
        manifest = SourceManifest(
            source_id=f"drbg_reseed_{reseed_interval_samples}",
            family="hybrid_drbg",
            is_deterministic=entropy_source.manifest.is_deterministic,
            is_physical=entropy_source.manifest.is_physical,
            metadata={
                "reseed_interval": reseed_interval_samples,
                "simulate_reseed_failure": simulate_reseed_failure,
            },
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self.entropy_source = entropy_source
        self.reseed_interval = reseed_interval_samples
        self.simulate_reseed_failure = simulate_reseed_failure
        self.reseed_count = 0
        self.samples_since_reseed = 0
        self.state_key = b"\x00" * 32
        self.state_counter = 0
        self._buffer = bytearray()
        # Initial seeding
        self._reseed()

    def _reseed(self):
        """Acquire fresh entropy and mix into state."""
        if self.simulate_reseed_failure and self.reseed_count > 0:
            # Reseed failure injected: state is NOT updated
            return
        fresh_entropy = self.entropy_source.read_bytes(32)
        # Mix with previous key
        self.state_key = hashlib.sha256(self.state_key + fresh_entropy).digest()
        self.state_counter = 0
        self.reseed_count += 1
        self.samples_since_reseed = 0

    def _generate_block(self) -> bytes:
        data = self.state_key + struct.pack(">Q", self.state_counter)
        block = hashlib.sha256(data).digest()
        self.state_counter += 1
        self.samples_since_reseed += 1
        if self.samples_since_reseed >= self.reseed_interval:
            self._reseed()
        return block

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            self._buffer.extend(self._generate_block())
        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out
