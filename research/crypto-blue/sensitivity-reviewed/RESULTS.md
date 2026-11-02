# RNGLab BLUE: measured synthetic campaign

Seed namespace: `blue-sensitivity-v1`. Independent paired replicates: 64. Operations per trajectory: 512; injected onset: 128 (zero-based).

These are local synthetic experiments. Statistical nonrejection does not establish security. All replay material is public; ‘healthy’ means the declared secret-unexposed simulation model.

## Weak-key experiment

| Uniform seed support | Successful bounded recovery | Mean guesses | Mean duplicate keys / 64 |
|---|---|---|---|
| 8 bits | 64/64 (94.4–100.0% exact 95% CI) | 125.3 | 7.22 |
| 12 bits | 64/64 (94.4–100.0% exact 95% CI) | 2105.9 | 0.64 |

Only the intentionally bounded seed dictionaries were searched. No production-strength AES key search was attempted.

## Conditioning blind spot

| Conditioner / expander | Known-seed output: three diagnostics do not reject | Reference output: do not reject | Exact attacker prediction |
|---|---|---|---|
| sha256/hmac | 63/64 (91.6–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) |
| sha256/chacha20 | 63/64 (91.6–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) |
| hkdf/hmac | 64/64 (94.4–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) |
| hkdf/chacha20 | 64/64 (94.4–100.0% exact 95% CI) | 63/64 (91.6–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) |

Monobit, block-frequency and runs diagnostics use Holm correction within each stream at 1%. This is a diagnostic subset, not the complete NIST STS battery. Intervals across rows are marginal, not simultaneous.

## Nonce failure and detection

| Condition | Observed collision trajectories | Source alarm trajectories | Repeated emissions with source gate | Repeated emissions with nonce guard |
|---|---|---|---|---|
| healthy | 0/64 (0.0–5.6% exact 95% CI) | 0/64 (0.0–5.6% exact 95% CI) | 0 | 0 |
| low_8 | 64/64 (94.4–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) | 0 | 0 |
| low_12 | 64/64 (94.4–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) | 0 | 0 |
| low_16 | 45/64 (57.6–81.1% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) | 0 | 0 |
| fixed | 64/64 (94.4–100.0% exact 95% CI) | 0/64 (0.0–5.6% exact 95% CI) | 24512 | 0 |
| serial | 10/64 (7.8–26.9% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) | 0 | 0 |
| drift | 63/64 (91.6–100.0% exact 95% CI) | 64/64 (94.4–100.0% exact 95% CI) | 0 | 0 |
| rollback | 64/64 (94.4–100.0% exact 95% CI) | 0/64 (0.0–5.6% exact 95% CI) | 8192 | 0 |
| counter_reset | 64/64 (94.4–100.0% exact 95% CI) | 0/64 (0.0–5.6% exact 95% CI) | 8192 | 0 |
| counter_persistent | 0/64 (0.0–5.6% exact 95% CI) | 0/64 (0.0–5.6% exact 95% CI) | 0 | 0 |

Source alarms during downstream fixed-nonce/counter-reset faults are alarms on an unchanged healthy raw source, not successful source-fault detections. Incidental/pre-fault quarantine can prevent misuse but is excluded from credited early warnings.

The same nonce paths are paired across all three AEADs; they are not three independent samples. For each first collision, the pipeline tests exact known-plaintext recovery; GCM additionally tests a new one-block forgery. ChaCha20-Poly1305 forgery is not implemented or claimed as measured.

| Algorithm | Collision paths | Demonstrated plaintext recovery | Accepted new-message forgery |
|---|---|---|---|
| AES-128-GCM | 438 | 438 | 438 |
| AES-256-GCM | 438 | 438 | 438 |
| ChaCha20-Poly1305 | 438 | 438 | not attempted |

Counts pooled over fault conditions are descriptive, not estimates of real-world prevalence.

Frozen healthy source-monitor calibration check: **PASS** (CP95 upper false-alarm bound 0.0560; criterion <=0.10).

Four-cell conditioning criterion (one-sided exact binomial versus 90% nonrejection, Holm family alpha .05): sha256/hmac: supported, sha256/chacha20: supported, hkdf/hmac: supported, hkdf/chacha20: supported.

## Comparison matrix

| RNG condition | Statistical tests | Cryptographic outcome | Security implication | Detection | Mitigation |
|---|---|---|---|---|---|
| Known 8-bit seed, SHA256/HKDF then HMAC/ChaCha expansion | Measured nonrejection rates above | Exact output prediction measured | Unpredictability absent when seed is known | Raw low-support health screen alarms; provenance identifies known seed | Require independently justified startup entropy; public salt adds no entropy |
| Bounded weak key seed | Output appearance is not the recovery criterion | AES-256-GCM known-plaintext dictionary recovery and duplicate Ed25519 keys measured | Key length exceeds effective seed uncertainty | Model entropy-budget gate blocks | OS/library key generation with trustworthy entropy infrastructure |
| Reduced nonce seed support | Raw monitor rates above | First collision, XOR disclosure and GCM forgery measured on colliding paths | Nonce uniqueness failure is distinct from key weakness | Source monitor plus per-key reservation guard | Durable counter or bounded guard |
| Downstream fixed nonce / counter reset | Raw source remains the paired reference | Misuse and its consequences measured | Healthy raw randomness cannot prevent an API/state-management fault | Nonce policy monitor at attempted reuse | Guard before encryption; shared persistent ledger |
| Unique persistent counter across simulated process restart | Counter output is intentionally structured; source reference unchanged | No observed reuse in this finite campaign | Unpredictable nonces are unnecessary for these AEADs if unique per key | Counter reservation invariants | Commit allocation before use; preserve state or rotate key |
| Known DRBG state plus known reseed | Output tests cannot establish secrecy | Clone continues to predict, measured | Unvalidated reseed cannot restore unknown state | State provenance / policy; no statistical proof | Withheld model input plus rekey changes clone outcome; past ciphertext remains exposed |
| ECDSA repeated signing nonce | No statistical battery needed | Synthetic private scalar recovered; new signature verifies | Message-specific nonce misuse compromises ECDSA key | Nonce policy; repeated r across distinct messages is evidence | RFC6979 deterministic signing; adequate key entropy |
| RFC6979 / Ed25519 deterministic signing | Repeated signatures of same message expected | Same-message determinism and distinct-message behavior measured | No fresh signing entropy needed by these constructions | Library known-answer checks | Protect key generation and private state |

## Signatures, recovery and controls

**Signatures:**

- `ecdsa_repeated_nonce_private_key_recovered`: 64/64.
- `ecdsa_recovered_key_signature_verified`: 64/64.
- `rfc6979_same_message_same_signature`: 64/64.
- `rfc6979_distinct_messages_distinct_r`: 64/64.
- `ed25519_same_message_same_signature`: 64/64.
- `ed25519_distinct_messages_distinct_signature`: 64/64.

**Recovery:**

- `known_state_prediction`: 64/64.
- `known_reseed_still_predictable`: 64/64.
- `withheld_reseed_breaks_old_clone_prediction`: 64/64.
- `old_state_clone_rejected_after_rekey`: 64/64.
- `past_ciphertext_still_decryptable`: 64/64.
- `retained_unknown_state_resists_known_reseed_clone`: 64/64.
- `unvalidated_startup_and_reseed_blocked`: 64/64.
- `restored_model_budget_releases_quarantine`: 64/64.
- `starvation_blocks_output`: 64/64.

Fresh-nonce negative controls: 0/192 successful XOR recoveries; 0 accepted GCM forgeries. Bit-flip rejection is recorded separately from forged-tag acceptance.

## Stability window and limits

The source score cutoff is frozen from separate healthy paths using a maximum over all monitored times and score components. Whole-path false alarms are measured on the evaluation controls. A nonce guard acts at the attempted policy violation; its prevention is not forecasting. Null events are administratively right-censored. Kaplan–Meier curves and exact common-horizon binomial intervals are recorded without extrapolating censored tails.

Alarm-to-compromise differences are conditional descriptive windows. A known seed is insecure immediately under the threat model, even before an observable collision. Identical byte histories can accompany different attacker knowledge or different future resets: no output-only universal warning time is identifiable. Birthday-risk curves are analytical results conditional on IID uniform finite support.

The initial suite demonstrates infrastructure failures and scoped prevention. It provides no break of AES or ChaCha20, no certified entropy estimate, no validated physical TRNG, and no deployment-wide forecast guarantee. See the frozen protocol, per-replicate JSON, source hashes, audit, and security limitations.
