"""
Finite research treatments and a conditional universal-hash reference.

SECURITY NOTICE:
These conditioning routines (von Neumann, SHA-256 block extraction, Toeplitz hashing)
are implemented strictly for experimental research, academic demonstration, and
stochastic benchmark treatments. They are NOT certified or audited for production
cryptography, cryptographic key generation, or security-sensitive parameters.

Passing statistical test batteries (such as NIST SP 800-22) after cryptographic
hashing or PRF expansion demonstrates only the diffusion properties of the hash,
NOT the entropy of the underlying physical source. True physical entropy must be
assessed directly on raw unconditioned noise using NIST SP 800-90B min-entropy
estimators and continuous health tests.
"""
from __future__ import annotations
import hashlib
import math
from pathlib import Path
from .integrity import digest, write_json


def pack_bits(bits):
    value = count = 0
    for bit in bits:
        if bit not in (0, 1):
            raise ValueError("expected binary symbols")
        value = (value << 1) | bit
        count += 1
        if count == 8:
            yield value
            value = count = 0


def treatments(symbol_path, out):
    """One symbol per byte; MSB-first packed outputs. No expansion or replay."""
    out = Path(out)
    out.mkdir(exist_ok=False)
    count = vn_count = equal = 0
    pending = None
    pack_value = pack_count = vn_value = vn_pending = 0
    hash_block = bytearray()
    with Path(symbol_path).open("rb") as inp, (out/"raw.bin").open("xb") as raw, \
         (out/"vn.bin").open("xb") as vn, (out/"sha256.bin").open("xb") as hashed:
        for chunk in iter(lambda: inp.read(65536), b""):
            for bit in chunk:
                if bit not in (0, 1):
                    raise ValueError("expected one binary symbol per byte")
                count += 1
                pack_value = (pack_value << 1) | bit
                pack_count += 1
                if pack_count == 8:
                    raw.write(bytes([pack_value]))
                    hash_block.append(pack_value)
                    pack_value = pack_count = 0
                    if len(hash_block) == 64:
                        hashed.write(hashlib.sha256(b"RNGLab15 research\0" + hash_block).digest())
                        hash_block.clear()
                if pending is None:
                    pending = bit
                else:
                    if pending != bit:
                        vn_value = (vn_value << 1) | pending
                        vn_pending += 1
                        vn_count += 1
                        if vn_pending == 8:
                            vn.write(bytes([vn_value]))
                            vn_value = vn_pending = 0
                    else:
                        equal += 1
                    pending = None
    result = {"input_symbol_sha256": digest(symbol_path), "input_bits": count,
              "raw_tail_bits_discarded": pack_count, "vn_equal_pairs_discarded": equal,
              "vn_output_bits_before_byte_truncation": vn_count, "vn_tail_bits_discarded": vn_pending,
              "vn_unpaired_input_bits": int(pending is not None),
              "hash_tail_input_bytes_discarded": len(hash_block),
              "validated_entropy_credit_bits": 0,
              "qualification": "research candidate only; VN needs pair assumptions; fixed hashing has no extraction bound",
              "outputs": {name: {"bytes": (out/f"{name}.bin").stat().st_size,
                                  "sha256": digest(out/f"{name}.bin")} for name in ("raw", "vn", "sha256")}}
    write_json(out/"treatments.json", result)
    return result


def toeplitz_reference(bits, seed_bits, output_bits, *, conditional_min_entropy, epsilon):
    """Universal Toeplitz hashing, ONLY conditional on external source/seed evidence.

    Leftover-hash sufficient budget m <= floor(k - 2 log2(1/epsilon)). This
    reference computes a transform; it cannot establish k or seed independence.
    Seed length n+m-1; rows use seed[n-1-i+j]. Explicitly capped quadratic work.
    """
    n = len(bits)
    if not 1 <= n <= 8192 or not 1 <= output_bits <= n:
        raise ValueError("bounded positive dimensions required")
    if not (math.isfinite(conditional_min_entropy) and 0 <= conditional_min_entropy <= n
            and math.isfinite(epsilon) and 0 < epsilon < 1):
        raise ValueError("invalid conditional extraction budget")
    if output_bits > math.floor(conditional_min_entropy - 2*math.log2(1/epsilon)):
        raise ValueError("insufficient entropy budget for requested security distance")
    if len(seed_bits) != n+output_bits-1 or any(b not in (0, 1) for b in (*bits, *seed_bits)):
        raise ValueError("invalid binary input or independent seed dimensions")
    return [sum(bits[i]*seed_bits[n-1-i+j] for i in range(n)) % 2 for j in range(output_bits)]
