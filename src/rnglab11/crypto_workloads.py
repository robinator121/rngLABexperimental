"""Small, test-only cryptographic consumers returning aggregate evidence.

Absence of observed collisions is not a claim of cryptographic security. Never
use experimental/degraded sources here to generate actual operational keys.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import math

from .sources import ByteSource, SourceBudgetExceeded


def _limits(count: int, width: int, max_count: int, max_bytes: int):
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise ValueError("count must be a positive integer")
    if not isinstance(width, int) or isinstance(width, bool) or not 1 <= width <= 4096:
        raise ValueError("record width must be in [1,4096]")
    if count > max_count or count * width > max_bytes:
        raise SourceBudgetExceeded("cryptographic workload exceeds bounded test limits")


def nonce_collision_experiment(source: ByteSource, *, count: int = 256, nonce_bytes: int = 12,
                               max_count: int = 100_000, max_bytes: int = 1_048_576):
    """Exact finite collision counts; values/nonces are not exported."""
    _limits(count, nonce_bytes, max_count, max_bytes)
    consumed_before = source.bytes_consumed
    counts = Counter(source.read(nonce_bytes) for _ in range(count))
    pair_collisions = sum(n * (n - 1) // 2 for n in counts.values())
    return {"workload": "test_only_nonce_uniqueness", "count": count,
            "nonce_bytes": nonce_bytes, "unique_values": len(counts),
            "duplicate_values": count - len(counts), "collision_pairs": pair_collisions,
            "largest_multiplicity": max(counts.values()),
            "expected_collision_pairs_if_iid_uniform": math.ldexp(count * (count - 1) / 2, -8 * nonce_bytes),
            "bytes_consumed": source.bytes_consumed - consumed_before,
            "security_claim": "none; IID uniform expectation is only a comparison model"}


def stream_repetition_experiment(source: ByteSource, *, count: int = 256, block_bytes: int = 16,
                                 max_count: int = 100_000, max_bytes: int = 1_048_576):
    """Count aligned repeated blocks and adjacent repetition, with finite memory."""
    _limits(count, block_bytes, max_count, max_bytes)
    consumed_before = source.bytes_consumed
    counts = Counter()
    previous = None
    adjacent = 0
    for _ in range(count):
        block = source.read(block_bytes)
        adjacent += int(previous == block)
        counts[block] += 1
        previous = block
    return {"workload": "test_only_aligned_stream_repetition", "count": count,
            "block_bytes": block_bytes, "unique_blocks": len(counts),
            "repeated_blocks": count - len(counts), "adjacent_repetitions": adjacent,
            "largest_multiplicity": max(counts.values()),
            "bytes_consumed": source.bytes_consumed - consumed_before,
            "security_claim": "none; alignment/finite sample can miss other repetitions"}


def ed25519_key_experiment(source: ByteSource, *, count: int = 16, max_count: int = 1024,
                            max_bytes: int = 32_768):
    """Derive/sign/verify with pyca cryptography, exporting no private key bytes.

    from_private_bytes accepts a 32-byte seed; it does not verify that the seed
    has adequate entropy. Only public-key digests and aggregate counts survive.
    """
    _limits(count, 32, max_count, max_bytes)
    try:
        from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey
        from cryptography.hazmat.primitives import serialization
    except ImportError as exc:
        raise RuntimeError("ed25519 test requires the optional cryptography dependency") from exc
    consumed_before = source.bytes_consumed
    public_digests = Counter()
    successful_verifications = 0
    message = b"RNGLab v0.11 bounded test-only key experiment"
    for _ in range(count):
        private_key = Ed25519PrivateKey.from_private_bytes(source.read(32))
        public_key = private_key.public_key()
        public_key.verify(private_key.sign(message), message)
        successful_verifications += 1
        public = public_key.public_bytes(serialization.Encoding.Raw, serialization.PublicFormat.Raw)
        public_digests[hashlib.sha256(public).hexdigest()] += 1
    return {"workload": "test_only_ed25519_from_32byte_seed", "count": count,
            "unique_public_keys": len(public_digests),
            "duplicate_public_keys": count - len(public_digests),
            "successful_sign_verify": successful_verifications,
            "bytes_consumed": source.bytes_consumed - consumed_before,
            "private_material_exported": False,
            "security_claim": "none; signing success and uniqueness do not establish seed security"}
