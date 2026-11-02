# RNGLab BLUE: measured synthetic campaign

Seed namespace: `blue-confirm-v1`. Independent paired replicates: 128. Operations per trajectory: 256; injected onset: 64 (zero-based).

These are local synthetic experiments. Statistical nonrejection does not establish security. All replay material is public; ‘healthy’ means the declared secret-unexposed simulation model.

## Weak-key experiment

| Uniform seed support | Successful bounded recovery | Mean guesses | Mean duplicate keys / 64 |
|---|---|---|---|
| 8 bits | 128/128 (97.2–100.0% exact 95% CI) | 126.8 | 7.06 |
| 12 bits | 128/128 (97.2–100.0% exact 95% CI) | 1980.1 | 0.59 |

Only the intentionally bounded seed dictionaries were searched. No production-strength AES key search was attempted.

## Conditioning blind spot

| Conditioner / expander | Known-seed output: three diagnostics do not reject | Reference output: do not reject | Exact attacker prediction |
|---|---|---|---|
| sha256/hmac | 128/128 (97.2–100.0% exact 95% CI) | 127/128 (95.7–100.0% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) |
| sha256/chacha20 | 121/128 (89.1–97.8% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) |
| hkdf/hmac | 128/128 (97.2–100.0% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) |
| hkdf/chacha20 | 127/128 (95.7–100.0% exact 95% CI) | 126/128 (94.5–99.8% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) |

Monobit, block-frequency and runs diagnostics use Holm correction within each stream at 1%. This is a diagnostic subset, not the complete NIST STS battery. Intervals across rows are marginal, not simultaneous.

## Nonce failure and detection

| Condition | Observed collision trajectories | Source alarm trajectories | Repeated emissions with source gate | Repeated emissions with nonce guard |
|---|---|---|---|---|
| healthy | 0/128 (0.0–2.8% exact 95% CI) | 8/128 (2.7–11.9% exact 95% CI) | 0 | 0 |
| low_8 | 128/128 (97.2–100.0% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) | 0 | 0 |
| low_12 | 126/128 (94.5–99.8% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) | 0 | 0 |
| low_16 | 35/128 (19.8–35.9% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) | 0 | 0 |
| fixed | 128/128 (97.2–100.0% exact 95% CI) | 8/128 (2.7–11.9% exact 95% CI) | 23456 | 0 |
| serial | 4/128 (0.9–7.8% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) | 0 | 0 |
| drift | 92/128 (63.2–79.5% exact 95% CI) | 128/128 (97.2–100.0% exact 95% CI) | 0 | 0 |
| rollback | 128/128 (97.2–100.0% exact 95% CI) | 6/128 (1.7–9.9% exact 95% CI) | 8000 | 0 |
| counter_reset | 128/128 (97.2–100.0% exact 95% CI) | 8/128 (2.7–11.9% exact 95% CI) | 7937 | 0 |
| counter_persistent | 0/128 (0.0–2.8% exact 95% CI) | 8/128 (2.7–11.9% exact 95% CI) | 0 | 0 |

The same nonce paths are paired across all three AEADs; they are not three independent samples. For each first collision, the pipeline tests exact known-plaintext recovery; GCM additionally tests a new one-block forgery. ChaCha20-Poly1305 forgery is not implemented or claimed as measured.

| Algorithm | Collision paths | Demonstrated plaintext recovery | Accepted new-message forgery |
|---|---|---|---|
| AES-128-GCM | 769 | 769 | 769 |
| AES-256-GCM | 769 | 769 | 769 |
| ChaCha20-Poly1305 | 769 | 769 | not attempted |

Counts pooled over fault conditions are descriptive, not estimates of real-world prevalence.

Frozen healthy source-monitor calibration check: **INSUFFICIENT EVIDENCE** (CP95 upper false-alarm bound 0.1194; criterion <=0.10).

Four-cell conditioning criterion (one-sided exact binomial versus 90% nonrejection, Holm family alpha .05): sha256/hmac: supported, sha256/chacha20: insufficient evidence, hkdf/hmac: supported, hkdf/chacha20: supported.

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

- `ecdsa_repeated_nonce_private_key_recovered`: 128/128.
- `ecdsa_recovered_key_signature_verified`: 128/128.
- `rfc6979_same_message_same_signature`: 128/128.
- `rfc6979_distinct_messages_distinct_r`: 128/128.
- `ed25519_same_message_same_signature`: 128/128.
- `ed25519_distinct_messages_distinct_signature`: 128/128.

**Recovery:**

- `known_state_prediction`: 128/128.
- `known_reseed_still_predictable`: 128/128.
- `withheld_reseed_breaks_old_clone_prediction`: 128/128.
- `old_state_clone_rejected_after_rekey`: 128/128.
- `past_ciphertext_still_decryptable`: 128/128.
- `retained_unknown_state_resists_known_reseed_clone`: 128/128.
- `unvalidated_startup_and_reseed_blocked`: 128/128.
- `restored_model_budget_releases_quarantine`: 128/128.
- `starvation_blocks_output`: 128/128.

Fresh-nonce negative controls: 0/384 successful XOR recoveries; 0 accepted GCM forgeries. Bit-flip rejection is recorded separately from forged-tag acceptance.

## Stability window and limits

The source score cutoff is frozen from separate healthy paths using a maximum over all monitored times and score components. Whole-path false alarms are measured on the evaluation controls. A nonce guard acts at the attempted policy violation; its prevention is not forecasting. Null events are administratively right-censored. Kaplan–Meier curves and exact common-horizon binomial intervals are recorded without extrapolating censored tails.

Alarm-to-compromise differences are conditional descriptive windows. A known seed is insecure immediately under the threat model, even before an observable collision. Identical byte histories can accompany different attacker knowledge or different future resets: no output-only universal warning time is identifiable. Birthday-risk curves are analytical results conditional on IID uniform finite support.

The initial suite demonstrates infrastructure failures and scoped prevention. It provides no break of AES or ChaCha20, no certified entropy estimate, no validated physical TRNG, and no deployment-wide forecast guarantee. See the frozen protocol, per-replicate JSON, source hashes, audit, and security limitations.
