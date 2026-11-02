import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.asymmetric import ec, ed25519, utils
from cryptography.hazmat.primitives.kdf.hkdf import HKDF

from rnglab11.drbg import HMACDRBG
from rnglab_crypto.primitives.aead import (ALGORITHMS, cipher, nonce_misuse_demo, NonceLedger,
                                         NonceReuse, GuardedAEAD, _mul, _inverse)
from rnglab_crypto.primitives.signatures import signature_demo
from rnglab_crypto.primitives.conditioning import condition, expand


@pytest.mark.parametrize("filename", ["hmac_drbg.json", "hmac_drbg_reseed.json"])
def test_nist_cavp_hmac_vectors(filename):
    v = json.loads((Path(__file__).parent/"crypto_vectors"/filename).read_text())
    drbg = HMACDRBG(bytes.fromhex(v["EntropyInput"][0]), bytes.fromhex(v["Nonce"][0]))
    if "EntropyInputReseed" in v:
        drbg.reseed(bytes.fromhex(v["EntropyInputReseed"][0]))
    drbg.read(128)
    assert drbg.read(128).hex() == v["ReturnedBits"][0]


@pytest.mark.parametrize("algorithm,size,expected", [
    ("AES-128-GCM", 16, "0388dace60b6a392f328c2b971b2fe78ab6e47d42cec13bdf53a67b21257bddf"),
    ("AES-256-GCM", 32, "cea7403d4d606b6e074ec5d3baf39d18d0d1c8a799996bf0265b98b5d48ab919")])
def test_gcm_fixed_zero_vectors(algorithm, size, expected):
    # GCM specification test cases: 96-bit zero IV, 128-bit zero plaintext, empty AAD.
    aead = cipher(algorithm, bytes(size))
    assert aead.encrypt(bytes(12), bytes(16), None).hex() == expected
    assert aead.decrypt(bytes(12), bytes.fromhex(expected), None) == bytes(16)


def test_rfc8439_aead_vector():
    key = bytes(range(0x80, 0xa0))
    nonce = bytes.fromhex("070000004041424344454647")
    aad = bytes.fromhex("50515253c0c1c2c3c4c5c6c7")
    plaintext = b"Ladies and Gentlemen of the class of '99: If I could offer you only one tip for the future, sunscreen would be it."
    expected = bytes.fromhex(
        "d31a8d34648e60db7b86afbc53ef7ec2a4aded51296e08fea9e2b5a736ee62d6"
        "3dbea45e8ca9671282fafb69da92728b1a71de0a9e060b2905d6a5b67ecd3b36"
        "92ddbd7f2d778b8c9803aee328091b58fab324e4fad675945585808b4831d7bc"
        "3ff4def08e4b7a9de576d26586cec64b61161ae10b594f09e26a7e902ecbd0600691")
    assert cipher("ChaCha20-Poly1305", key).encrypt(nonce, plaintext, aad) == expected


def test_rfc6979_p256_sha256_sample():
    key = ec.derive_private_key(int("C9AFA9D845BA75166B5C215767B1D6934E50C3DB36E89B127B8A622B120F6721", 16), ec.SECP256R1())
    signature = key.sign(b"sample", ec.ECDSA(hashes.SHA256(), deterministic_signing=True))
    r, s = utils.decode_dss_signature(signature)
    assert r == int("EFD48B2AACB6A8FD1140DD9CD45E81D69D2C877B56AAF991C34D0EA84EAF3716", 16)
    assert s == int("F7CB1C942D657C41D436C7A1B6E29F65F3E900DBB9AFF4064DC4AB2F843ACDA8", 16)


def test_rfc8032_empty_message():
    key = ed25519.Ed25519PrivateKey.from_private_bytes(bytes.fromhex("9d61b19deffd5a60ba844af492ec2cc44449c5697b326919703bac031cae7f60"))
    assert key.public_key().public_bytes_raw().hex() == "d75a980182b10ab7d54bfed3c964073a0ee172f3daa62325af021a68f707511a"
    assert key.sign(b"").hex() == ("e5564300c360ac729086e2cc806e828a84877f1eb8e5d974d873e065224901555f"
                                  "b8821590a33bacc61e39701cf9b46bd25bf5f0595bbe24655141438e7a100b")


def test_rfc5869_case1():
    result = HKDF(hashes.SHA256(), 42, bytes(range(13)), bytes(range(0xf0, 0xfa))).derive(b"\x0b"*22)
    assert result.hex() == "3cb25f25faacd57a90434f64d0362f2a2d2d0a90cf1a5a4c5db02d56ecc4c5bf34007208d5b887185865"


@pytest.mark.parametrize("algorithm", ALGORITHMS)
def test_aead_misuse_and_fresh_nonce_negative_control(algorithm):
    bad = nonce_misuse_demo(algorithm, b"test")
    good = nonce_misuse_demo(algorithm, b"test", reuse=False)
    assert bad["xor_plaintext_recovered"] and not good["xor_plaintext_recovered"]
    assert bad["unaided_bit_flip_rejected"] and good["unaided_bit_flip_rejected"]
    if algorithm.startswith("AES"):
        assert bad["new_message_forgery_accepted"] and not good["new_message_forgery_accepted"]
    else:
        assert bad["new_message_forgery_accepted"] is None


def test_gcm_field_identity_and_inverse():
    x = 0x123456789ABCDEF
    assert _mul(x, 1 << 127) == x
    assert _mul(x, _inverse(x)) == 1 << 127
    with pytest.raises(ValueError):
        _inverse(0)


def test_signature_demonstration():
    assert all(v for v in signature_demo(b"fixture").values() if type(v) is bool)


def test_nonce_guard_blocks_before_library_encrypt(monkeypatch):
    guarded = GuardedAEAD("AES-128-GCM", bytes(16))
    guarded.encrypt(bytes(12), b"synthetic")
    class MustNotEncrypt:
        def encrypt(self, *args):
            pytest.fail("duplicate reached primitive")
    guarded.aead = MustNotEncrypt()
    with pytest.raises(NonceReuse):
        guarded.encrypt(bytes(12), b"synthetic")


def test_nonce_ledger_restart_concurrency_and_exhaustion(tmp_path):
    path = tmp_path/"nonce.sqlite"
    ledger = NonceLedger(path, bytes(32))
    assert ledger.allocate() == bytes(12)
    ledger.close()
    def writer(_):
        source = NonceLedger(path, bytes(32))
        try:
            return [source.allocate() for _ in range(8)]
        finally:
            source.close()
    with ThreadPoolExecutor(max_workers=3) as pool:
        nonces = [n for batch in pool.map(writer, range(3)) for n in batch]
    assert len(set(nonces)) == 24 and bytes(12) not in nonces
    ledger = NonceLedger(path, bytes(32))
    ledger.db.execute("UPDATE counters SET next=?", (2**63-1,))
    with pytest.raises(OverflowError):
        ledger.allocate()
    ledger.close()


@pytest.mark.parametrize("method", ["sha256", "hkdf"])
@pytest.mark.parametrize("mechanism", ["hmac", "chacha20"])
def test_known_conditioned_seed_exactly_predicts(method, mechanism):
    seed = condition(bytes(32), method)
    assert expand(seed, 1024, mechanism) == expand(seed, 1024, mechanism)
    assert expand(seed, 1024, mechanism) != bytes(1024)
