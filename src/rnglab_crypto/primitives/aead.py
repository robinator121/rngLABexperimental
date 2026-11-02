"""AEAD consumers and fixed one-block synthetic misuse experiments."""
from __future__ import annotations

import hashlib
import sqlite3

from cryptography.exceptions import InvalidTag
from cryptography.hazmat.primitives.ciphers.aead import AESGCM, ChaCha20Poly1305

ALGORITHMS = ("AES-128-GCM", "AES-256-GCM", "ChaCha20-Poly1305")


def cipher(algorithm, key):
    if algorithm not in ALGORITHMS:
        raise ValueError("unsupported AEAD")
    size = 16 if algorithm == "AES-128-GCM" else 32
    if len(key) != size:
        raise ValueError("wrong key size")
    return ChaCha20Poly1305(key) if algorithm == "ChaCha20-Poly1305" else AESGCM(key)


class NonceReuse(RuntimeError):
    pass


class NonceGuard:
    """Finite, per-key in-process duplicate guard; no persistent/replay guarantees."""
    def __init__(self, limit=100_000):
        self.seen = set()
        self.limit = limit

    def reserve(self, nonce):
        if len(nonce) != 12:
            raise ValueError("96-bit nonce required")
        if nonce in self.seen:
            raise NonceReuse("duplicate candidate blocked before encryption")
        if len(self.seen) >= self.limit:
            raise RuntimeError("nonce guard full; fail closed")
        self.seen.add(nonce)


class NonceLedger:
    """Transactional per-key counter; commit reservation BEFORE returning nonce.

    All writers for a key MUST use the same durable database. File loss, database
    snapshot rollback, unrelated ledgers and malicious storage are outside this
    guarantee. Rotate keys if durable state cannot be trusted. Single-host
    SQLite transactions serialize concurrent writers; crashes can skip values.
    """
    def __init__(self, path, key: bytes):
        if len(key) not in (16, 32):
            raise ValueError("AEAD key required")
        self.key_id = hashlib.sha256(key).hexdigest()
        self.db = sqlite3.connect(str(path), timeout=30, isolation_level=None)
        self.db.execute("PRAGMA synchronous=FULL")
        self.db.execute("CREATE TABLE IF NOT EXISTS counters (key_id TEXT PRIMARY KEY, next INTEGER NOT NULL)")
        self.db.execute("INSERT OR IGNORE INTO counters VALUES (?, 0)", (self.key_id,))

    def allocate(self):
        self.db.execute("BEGIN IMMEDIATE")
        try:
            n, = self.db.execute("SELECT next FROM counters WHERE key_id=?", (self.key_id,)).fetchone()
            if not 0 <= n < 2**63-1:
                raise OverflowError("counter exhausted or corrupt; rotate key")
            self.db.execute("UPDATE counters SET next=? WHERE key_id=?", (n+1, self.key_id))
            self.db.execute("COMMIT")
            return n.to_bytes(12, "big")
        except BaseException:
            self.db.execute("ROLLBACK")
            raise

    def close(self):
        self.db.close()


class GuardedAEAD:
    """Opt-in local consumer: reject reused nonce before invoking encryption."""
    def __init__(self, algorithm, key, limit=100_000):
        self.aead = cipher(algorithm, key)
        self.guard = NonceGuard(limit)

    def encrypt(self, nonce, plaintext, aad=None):
        self.guard.reserve(nonce)
        return self.aead.encrypt(nonce, plaintext, aad)


def _mul(x, y):
    """GCM field multiplication, SP800-38D algorithm 1. Analysis only."""
    z, v = 0, y
    for i in range(128):
        if x & (1 << (127-i)):
            z ^= v
        v = (v >> 1) ^ (0xE1000000000000000000000000000000 if v & 1 else 0)
    return z


def _inverse(x):
    if not x:
        raise ValueError("zero has no field inverse")
    result, power = 1 << 127, 2**128-2
    while power:
        if power & 1:
            result = _mul(result, x)
        x = _mul(x, x)
        power >>= 1
    return result


def _xor(a, b):
    if len(a) != len(b):
        raise ValueError("equal-length operands required")
    return bytes(x ^ y for x, y in zip(a, b))


def _gcm_one_block_candidate(first, second):
    """Internal algebra for two equal-length, no-AAD, one-block lab records."""
    c1, c2 = int.from_bytes(first[:16], "big"), int.from_bytes(second[:16], "big")
    t1, t2 = int.from_bytes(first[16:], "big"), int.from_bytes(second[16:], "big")
    h_squared = _mul(t1 ^ t2, _inverse(c1 ^ c2))
    changed = c1 ^ (1 << 127)
    tag = t1 ^ _mul(changed ^ c1, h_squared)
    return changed.to_bytes(16, "big") + tag.to_bytes(16, "big")


def nonce_misuse_demo(algorithm: str, synthetic_seed: bytes, *, reuse=True):
    """Generates its own local records; takes no external ciphertexts/targets.

    Reused nonce: exact known-plaintext XOR recovery for all three AEADs;
    accepted new-message forgery for GCM. ChaCha integrity misuse is analytical
    only in this initial suite; a bit flip is still rejected and is reported.
    """
    key = hashlib.sha256(b"BLUE demo"+synthetic_seed).digest()
    key = key[:16] if algorithm == "AES-128-GCM" else key
    aead = cipher(algorithm, key)
    n1, n2 = bytes(12), bytes(12) if reuse else (1).to_bytes(12, "big")
    p1, p2 = b"lab message one!", b"lab message two!"
    c1, c2 = aead.encrypt(n1, p1, None), aead.encrypt(n2, p2, None)
    recovered = _xor(_xor(c1[:16], c2[:16]), p1)
    tampered = bytes([c1[0] ^ 128]) + c1[1:]
    try:
        aead.decrypt(n1, tampered, None)
        tamper_rejected = False
    except InvalidTag:
        tamper_rejected = True
    forgery_accepted = None
    if algorithm.startswith("AES"):
        forged = _gcm_one_block_candidate(c1, c2)
        try:
            plaintext = aead.decrypt(n1, forged, None)
            forgery_accepted = plaintext != p1 and plaintext == bytes([p1[0] ^ 128])+p1[1:]
        except InvalidTag:
            forgery_accepted = False
    return {"xor_plaintext_recovered": recovered == p2,
            "new_message_forgery_accepted": forgery_accepted,
            "unaided_bit_flip_rejected": tamper_rejected,
            "integrity_scope": "one_block_GCM_measured" if algorithm.startswith("AES") else "forgery_not_attempted"}
