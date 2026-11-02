"""Contained P-256 nonce-reuse algebra, verified by maintained library code."""
import hashlib

from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.asymmetric import ec, ed25519, utils

P256_ORDER = 0xFFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551


def signature_demo(synthetic_seed: bytes):
    q = P256_ORDER
    d = int.from_bytes(hashlib.sha256(b"private"+synthetic_seed).digest(), "big") % (q-1)+1
    k = int.from_bytes(hashlib.sha256(b"nonce"+synthetic_seed).digest(), "big") % (q-1)+1
    private = ec.derive_private_key(d, ec.SECP256R1())
    public = private.public_key()
    messages = (b"BLUE synthetic signature A", b"BLUE synthetic signature B")
    z = [int.from_bytes(hashlib.sha256(m).digest(), "big") for m in messages]
    # Point multiplication comes from the library, only signature scalar algebra
    # is exposed here to construct deliberate misuse unavailable in its safe API.
    r = ec.derive_private_key(k, ec.SECP256R1()).public_key().public_numbers().x % q
    s = [(pow(k, -1, q)*(h+r*d)) % q for h in z]
    if r == 0 or 0 in s or s[0] == s[1]:
        raise RuntimeError("degenerate synthetic demonstration")
    for msg, scalar in zip(messages, s):
        public.verify(utils.encode_dss_signature(r, scalar), msg, ec.ECDSA(hashes.SHA256()))
    recovered_k = (z[0]-z[1])*pow(s[0]-s[1], -1, q) % q
    recovered_d = (s[0]*recovered_k-z[0])*pow(r, -1, q) % q
    forged_message = b"BLUE locally forged synthetic message"
    forged = ec.derive_private_key(recovered_d, ec.SECP256R1()).sign(
        forged_message, ec.ECDSA(hashes.SHA256(), deterministic_signing=True))
    public.verify(forged, forged_message, ec.ECDSA(hashes.SHA256()))
    deterministic = ec.ECDSA(hashes.SHA256(), deterministic_signing=True)
    a, b = [private.sign(m, deterministic) for m in messages]
    for msg, sig in zip(messages, (a, b)):
        public.verify(sig, msg, ec.ECDSA(hashes.SHA256()))
    ed = ed25519.Ed25519PrivateKey.from_private_bytes(hashlib.sha256(synthetic_seed).digest())
    ed_sig = ed.sign(messages[0])
    ed.public_key().verify(ed_sig, messages[0])
    return {"ecdsa_repeated_nonce_private_key_recovered": recovered_d == d,
            "ecdsa_recovered_key_signature_verified": True,
            "rfc6979_same_message_same_signature": a == private.sign(messages[0], deterministic),
            "rfc6979_distinct_messages_distinct_r": utils.decode_dss_signature(a)[0] != utils.decode_dss_signature(b)[0],
            "ed25519_same_message_same_signature": ed_sig == ed.sign(messages[0]),
            "ed25519_distinct_messages_distinct_signature": ed_sig != ed.sign(messages[1]),
            "signing_rng_dependency": "RFC6979 and Ed25519 need no fresh per-message entropy; key generation still does"}
