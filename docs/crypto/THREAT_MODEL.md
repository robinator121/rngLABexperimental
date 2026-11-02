# Cryptographic threat model

This package is an academic, local laboratory. All plaintexts, keys, seeds, nonce histories and signatures in attack demonstrations are generated inside the harness. There are no network-target, credential, imported ciphertext, or external-key recovery interfaces. Replay seeds are intentionally public. The experiment's “healthy secret” is a **modeled information restriction**, never a claim that a published seed is secret.

The causal layers are raw source → SHA256/HKDF conditioning → HMAC-DRBG or ChaCha20 deterministic expansion → key/nonce/signature consumer → observed security outcome → defensive policy. The nonce experiment also includes direct downstream faults that bypass the source. This distinction matters: it would be a confounder to make a downstream reset artificially visible as a raw-source defect.

| Experiment | Attacker knowledge/capability | Claimed and excluded outcomes |
|---|---|---|
| Bounded weak keys | Known KDF, 4–16-bit deliberately finite seed dictionary, one synthetic AEAD record with known plaintext/AAD | Enumerated seed recovery only; no search of production-strength key space |
| Nonce reuse | Two same-key same-nonce synthetic records; one known plaintext; public tags, nonce, algorithm | Exact XOR disclosure. GCM one-block, equal-length, empty-AAD tag algebra produces an accepted changed message. No universal forgery API |
| ECDSA | Two signatures on distinct known messages under one deliberately repeated P-256 scalar nonce | Recover synthetic private scalar and verify a new signature. Not a claim about properly used ECDSA or Ed25519 |
| Conditioning | Known short seed, public salt/nonce, request schedule | Exact prediction after expansion despite statistical nonrejection; no proof of randomness from finite diagnostics |
| State compromise | Exact modeled DRBG snapshot and reseed inputs that are either known or withheld | Measure clone prediction. A withheld independent input is a model assumption; divergent clones alone are not a proof of full security restoration |
| Partial exposure | Prefix of a bounded 4–16-bit seed/state identifier | Search remaining suffix; not an attack on partially exposed native HMAC/ChaCha internal state |
| Defensive monitoring | Raw blocks, seed provenance/budget assumptions, per-key nonce allocation history | Detect only visible deviations. Attacker knowledge and future reset regimes are not identifiable from byte distributions alone |

Security events are separate fields:

1. Fault onset: laboratory injection time.
2. Model insecurity onset: a declared reduction of attacker uncertainty or loss of a required invariant; not necessarily observed misuse. For small-support/serial/drift faults this flags a departure from the intended source model, not a proof that every resulting ciphertext is compromised.
3. First source alarm: observed score threshold crossing, possibly pre-fault or a healthy false alarm.
4. First policy violation: attempted reuse of a nonce under one key.
5. First demonstrated compromise: exact plaintext recovery; GCM forgery has a distinct field.
6. Mitigation time: source quarantine or pre-encryption nonce rejection; successful prevention requires no violating emission.

Indices are zero-based operations. `null` means not observed within the horizon, not “secure forever.” Compromise is immediate at reuse in this chosen-known-plaintext harness; no post-reuse computational attack delay is estimated. For known seeds, theoretical insecurity can precede any alarm or collision. Negative lead time is retained, not clipped.

The experiment does not model side channels, malicious primitives, protocol replay protection, distributed key lifecycle, malicious OS/hardware, backups containing keys, or remote attacks. A local duplicate-nonce guard has finite memory and no crash persistence. The SQLite counter is per key and safe across tested reopen/concurrent-writer scenarios only when every writer shares trustworthy durable storage. Snapshot rollback and lost databases require key rotation or an external monotonic authority; SQLite does not solve those threats.

Paired arms are isolated counterfactual service worlds. They deliberately share replay key/source material to control nuisance variation; their ciphertexts are never pooled into one real key epoch. Shared material across algorithm comparisons is likewise experimental pairing, not a recommendation to reuse production keys across algorithms. Monitoring is offline sequential replay with causal per-block scores and operation-indexed event times, not a deployed real-time service or measured wall-clock reaction deadline.

Maintained `cryptography` supplies AES-GCM, ChaCha20-Poly1305, P-256 point operations/deterministic signatures, Ed25519 and HKDF. The small GCM field and ECDSA scalar calculations construct test misuse; they are not production cipher/signature implementations.
