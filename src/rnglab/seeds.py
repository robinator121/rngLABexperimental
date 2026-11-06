"""RNGLab v0.18: namespace-derived seeds.

A seed is the first 8 bytes (big-endian) of SHA-256 of
"<namespace>/<part>/<part>/...". Replicas inside one seeded unit are
numpy.random.SeedSequence(seed).spawn(count) children. Distinct strings give
distinct seeds with overwhelming probability; that prevents accidental stream
reuse but is not a proof of statistical independence between streams.
"""
from __future__ import annotations

import hashlib

import numpy as np


def namespace_seed(namespace: str, *parts: object) -> int:
    """Return the 64-bit seed for namespace/part/part/..."""
    if not isinstance(namespace, str) or not namespace:
        raise ValueError("namespace must be a nonempty string")
    label = "/".join([namespace, *map(str, parts)])
    return int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")


def replica_generators(seed: int, count: int) -> list[np.random.Generator]:
    """PCG64 generators for the SeedSequence(seed).spawn(count) children."""
    if type(count) is not int or count < 1:
        raise ValueError("count must be a positive integer")
    return [np.random.default_rng(child) for child in np.random.SeedSequence(seed).spawn(count)]
