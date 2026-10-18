"""NIST SP 800-90A Rev. 1 Deterministic Random Bit Generators (DRBGs).

Implements:
1. HMAC_DRBG (Section 10.1.2) with SHA-256
2. Known Answer Test (KAT) vectors for automated verification
"""
from __future__ import annotations

import hashlib
import hmac
import struct
from typing import Any


from rnglab11.drbg import HMACDRBG as _HMACDRBG


class HMAC_DRBG(_HMACDRBG):
    """NIST SP 800-90A Rev. 1 HMAC_DRBG using SHA-256."""

    def __init__(
        self,
        entropy_input: bytes,
        nonce: bytes,
        personalization_string: bytes = b"",
        security_strength: int = 256,
    ):
        super().__init__(
            entropy_input=bytes(entropy_input),
            nonce=bytes(nonce),
            personalization=bytes(personalization_string),
            security_strength=security_strength,
        )

    def reseed(self, entropy_input: bytes, additional_input: bytes = b"") -> None:
        super().reseed(bytes(entropy_input), additional_input=bytes(additional_input))

    def generate(self, num_bytes: int, additional_input: bytes = b"") -> bytes:
        return self.read(num_bytes, additional_input=bytes(additional_input))


# Official NIST SP 800-90A Worked Example Vectors (SHA-256, Strength 128)
# https://csrc.nist.gov/CSRC/media/Projects/Cryptographic-Standards-and-Guidelines/documents/examples/HMAC_DRBG.pdf
NIST_WORKED_EXAMPLE = {
    "entropy_input": bytes(range(0x37)),
    "nonce": bytes(range(0x20, 0x28)),
    "first_block_hex": (
        "D67B8C1734F46FA3F763CF57C6F9F4F2DC1089BD8BC1F6F023950BFC56176352"
        "08C8501238AD7A4400DEFEE46C640B61AF77C2D1A3BFAA90EDE5D207406E5403"
    ),
    "second_block_hex": (
        "8FDAEC20F8B421407059E3588920DA7EDA9DCE3CF8274DFA1C59C108C1D0AA9B"
        "0FA38DA5C792037C4D33CD070CA7CD0C5608DBA8B885654639DE2187B74CB263"
    ),
}


def verify_kat_hmac_drbg() -> bool:
    """Run NIST SP 800-90A Known Answer Test verification."""
    v = NIST_WORKED_EXAMPLE
    drbg = HMAC_DRBG(
        entropy_input=v["entropy_input"],
        nonce=v["nonce"],
        security_strength=128,
    )
    first = drbg.generate(64)
    second = drbg.generate(64)
    return (
        first == bytes.fromhex(v["first_block_hex"])
        and second == bytes.fromhex(v["second_block_hex"])
    )
