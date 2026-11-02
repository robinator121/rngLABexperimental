"""Conditioning is not an entropy certificate. Expansion cannot add seed entropy."""
import hashlib

from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms
from cryptography.hazmat.primitives.kdf.hkdf import HKDF
from rnglab11.drbg import HMACDRBG


def condition(raw: bytes, method: str = "sha256") -> bytes:
    if method == "sha256":
        return hashlib.sha256(raw).digest()
    if method == "hkdf":
        # Public fixed salt provides domain separation, not additional entropy.
        return HKDF(hashes.SHA256(), 32, b"RNGLab BLUE public salt", b"research seed").derive(raw)
    raise ValueError("unknown conditioning method")


def expand(seed: bytes, n: int, mechanism: str = "hmac") -> bytes:
    if len(seed) != 32 or type(n) is not int or not 1 <= n <= 65536:
        raise ValueError("32-byte seed and bounded positive output required")
    if mechanism == "hmac":
        return HMACDRBG(seed, b"BLUE public nonce", security_strength=128).read(n)
    if mechanism == "chacha20":
        # pyca's raw ChaCha20: 64-bit counter || 64-bit nonce, all zero here.
        # Each experiment instantiates a new stream. This is not a NIST DRBG.
        return Cipher(algorithms.ChaCha20(seed, bytes(16)), mode=None).encryptor().update(bytes(n))
    raise ValueError("unknown expansion mechanism")
