"""RNGLab v0.20 nuclear module: from RNGLab random sources to numbers.

Every random number the nuclear simulations use is read through
``RandomSource.read_bytes`` (rnglab14.sources.base), the audited consumption
boundary that counts and SHA-256 hashes the emitted bytes. This module fixes
how those bytes become numbers:

- a *word* is 8 consecutive bytes read big-endian;
- a *uniform* is ``(word >> 11) * 2**-53``, a multiple of 2**-53 in [0, 1).

Both conventions are those of ``RandomSource.random_float``, so a uniform here
equals the uniform any other RNGLab consumer would read from the same bytes.
``WordStream`` never buffers: a request for n words reads exactly 8n bytes, so
the source's byte count and stream digest cover exactly the words handed out.
A fault adapter (any RandomSource wrapping another) placed between a source and
a WordStream therefore acts on exactly the words the physics sees.

``NumpyPCG64Source`` exposes numpy's PCG64 raw 64-bit outputs as a RandomSource,
seeded as rnglab.seeds does (``SeedSequence(seed)``). Word k of the source is
the k-th raw output. Distinct seed labels prevent accidental stream reuse; they
do not prove statistical independence between streams.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from rnglab.seeds import namespace_seed
from rnglab14.sources.base import RandomSource, SourceManifest

WORD_BYTES = 8
UNIFORM_SCALE = 2.0 ** -53


def words_to_uniforms(words: np.ndarray) -> np.ndarray:
    """Uniforms (word >> 11) * 2**-53 for an array of uint64 words."""
    words = np.asarray(words, dtype=np.uint64)
    return (words >> np.uint64(11)).astype(np.float64) * UNIFORM_SCALE


def _count(n: object, what: str) -> int:
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 0:
        raise ValueError(f"{what} must be a nonnegative integer")
    return int(n)


class WordStream:
    """Hands out 64-bit words read from one RandomSource, counting them.

    ``first_source_byte`` is the source's byte count when the stream was
    created: an offset into the source that the stream did not read itself.
    """

    def __init__(self, source: RandomSource):
        if not isinstance(source, RandomSource):
            raise TypeError("source must be an rnglab RandomSource")
        self.source = source
        self.words_read = 0
        self.first_source_byte = int(source.consumed_bytes)

    def raw(self, n: int) -> bytes:
        """Exactly 8n bytes from the source (n words)."""
        n = _count(n, "word count")
        data = self.source.read_bytes(WORD_BYTES * n)
        self.words_read += n
        return data

    def words(self, n: int) -> np.ndarray:
        """n words as uint64 (big-endian decoding)."""
        return np.frombuffer(self.raw(n), dtype=">u8").astype(np.uint64)

    def uniforms(self, n: int) -> np.ndarray:
        """n uniforms in [0, 1), one per word."""
        return words_to_uniforms(self.words(n))

    @property
    def bytes_read(self) -> int:
        return WORD_BYTES * self.words_read

    def provenance(self) -> dict[str, Any]:
        """Identity and consumption of the stream, for run manifests."""
        manifest = self.source.manifest.to_dict()
        return {
            "source": manifest,
            "word_convention": "8-byte big-endian words; uniform = (word >> 11) * 2**-53",
            "first_source_byte": self.first_source_byte,
            "words_read": self.words_read,
            "bytes_read": self.bytes_read,
            "source_consumed_bytes": int(self.source.consumed_bytes),
            "source_stream_sha256": self.source.stream_digest,
        }


class NumpyPCG64Source(RandomSource):
    """RandomSource over numpy's PCG64 raw 64-bit outputs (big-endian words).

    ``start_word`` jumps the generator ahead by that many outputs with
    ``PCG64.advance``; skipped words are not read, so they are neither counted
    nor hashed. The jump is recorded in the manifest.
    """

    def __init__(self, seed: int, *, label: str = "", start_word: int = 0, max_bytes: int | None = None):
        if isinstance(seed, bool) or not isinstance(seed, (int, np.integer)) or seed < 0:
            raise ValueError("seed must be a nonnegative integer")
        start_word = _count(start_word, "start_word")
        seed = int(seed)
        bit_generator = np.random.PCG64(np.random.SeedSequence(seed))
        if start_word:
            bit_generator.advance(start_word)
        manifest = SourceManifest(
            source_id=f"numpy-pcg64:{seed:016x}+{start_word}",
            family="reference_prng",
            is_deterministic=True,
            is_physical=False,
            seed=seed,
            metadata={
                "generator": "numpy PCG64, raw 64-bit outputs",
                "seeding": "numpy.random.SeedSequence(seed)",
                "seed_label": label,
                "start_word": start_word,
                "numpy": np.__version__,
            },
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self._bit_generator = bit_generator
        self._pending = b""

    def _read_bytes(self, n: int) -> bytes:
        missing = n - len(self._pending)
        if missing > 0:
            count = (missing + WORD_BYTES - 1) // WORD_BYTES
            fresh = self._bit_generator.random_raw(count).astype(">u8").tobytes()
            buffer = self._pending + fresh
        else:
            buffer = self._pending
        out, self._pending = buffer[:n], buffer[n:]
        return out


def labelled_source(namespace: str, *parts: object, start_word: int = 0) -> NumpyPCG64Source:
    """NumpyPCG64Source seeded with namespace_seed(namespace, *parts)."""
    label = "/".join([namespace, *map(str, parts)])
    return NumpyPCG64Source(namespace_seed(namespace, *parts), label=label, start_word=start_word)
