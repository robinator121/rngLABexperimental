# RNGLab BLUE: measured synthetic campaign

Seed namespace: `blue-smoke-v1`. Independent paired replicates: 4. Operations per trajectory: 64; injected onset: 16 (zero-based).

These are local synthetic experiments. Statistical nonrejection does not establish security. All replay material is public; ‘healthy’ means the declared secret-unexposed simulation model.

## Weak-key experiment

| Uniform seed support | Successful bounded recovery | Mean guesses | Mean duplicate keys / 64 |
|---|---|---|---|
| 4 bits | 4/4 (39.8–100.0% exact 95% CI) | 8.0 | 48.00 |

Only the intentionally bounded seed dictionaries were searched. No production-strength AES key search was attempted.

## Conditioning blind spot

| Conditioner / expander | Known-seed output: three diagnostics do not reject | Reference output: do not reject | Exact attacker prediction |
|---|---|---|---|
| sha256/hmac | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) |
| sha256/chacha20 | 3/4 (19.4–99.4% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) |
| hkdf/hmac | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) |
| hkdf/chacha20 | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) |

Monobit, block-frequency and runs diagnostics use Holm correction within each stream at 1%. This is a diagnostic subset, not the complete NIST STS battery. Intervals across rows are marginal, not simultaneous.

## Nonce failure and detection

| Condition | Observed collision trajectories | Source alarm trajectories | Repeated emissions with source gate | Repeated emissions with nonce guard |
|---|---|---|---|---|
| healthy | 0/4 (0.0–60.2% exact 95% CI) | 1/4 (0.6–80.6% exact 95% CI) | 0 | 0 |
| low_4 | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) | 0 | 0 |
| low_8 | 4/4 (39.8–100.0% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) | 0 | 0 |
| fixed | 4/4 (39.8–100.0% exact 95% CI) | 1/4 (0.6–80.6% exact 95% CI) | 166 | 0 |
| serial | 0/4 (0.0–60.2% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) | 0 | 0 |
| drift | 1/4 (0.6–80.6% exact 95% CI) | 4/4 (39.8–100.0% exact 95% CI) | 0 | 0 |
| rollback | 4/4 (39.8–100.0% exact 95% CI) | 1/4 (0.6–80.6% exact 95% CI) | 64 | 0 |
| counter_reset | 4/4 (39.8–100.0% exact 95% CI) | 1/4 (0.6–80.6% exact 95% CI) | 64 | 0 |
| counter_persistent | 0/4 (0.0–60.2% exact 95% CI) | 1/4 (0.6–80.6% exact 95% CI) | 0 | 0 |

The same nonce paths are paired across all three AEADs; they are not three independent samples. For each first collision, the pipeline tests exact known-plaintext recovery; GCM additionally tests a new one-block forgery. ChaCha20-Poly1305 forgery is not implemented or claimed as measured.

| Algorithm | Collision paths | Demonstrated plaintext recovery | Accepted new-message forgery |
|---|---|---|---|
| AES-128-GCM | 21 | 21 | 21 |
| AES-256-GCM | 21 | 21 | 21 |
| ChaCha20-Poly1305 | 21 | 21 | not attempted |

Counts pooled over fault conditions are descriptive, not estimates of real-world prevalence.

Frozen healthy source-monitor calibration check: **INSUFFICIENT EVIDENCE** (CP95 upper false-alarm bound 0.8059; criterion <=0.10).

Four-cell conditioning criterion (one-sided exact binomial versus 90% nonrejection, Holm family alpha .05): sha256/hmac: insufficient evidence, sha256/chacha20: insufficient evidence, hkdf/hmac: insufficient evidence, hkdf/chacha20: insufficient evidence.

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

- `ecdsa_repeated_nonce_private_key_recovered`: 4/4.
- `ecdsa_recovered_key_signature_verified`: 4/4.
- `rfc6979_same_message_same_signature`: 4/4.
- `rfc6979_distinct_messages_distinct_r`: 4/4.
- `ed25519_same_message_same_signature`: 4/4.
- `ed25519_distinct_messages_distinct_signature`: 4/4.

**Recovery:**

- `known_state_prediction`: 4/4.
- `known_reseed_still_predictable`: 4/4.
- `withheld_reseed_breaks_old_clone_prediction`: 4/4.
- `old_state_clone_rejected_after_rekey`: 4/4.
- `past_ciphertext_still_decryptable`: 4/4.
- `retained_unknown_state_resists_known_reseed_clone`: 4/4.
- `unvalidated_startup_and_reseed_blocked`: 4/4.
- `restored_model_budget_releases_quarantine`: 4/4.
- `starvation_blocks_output`: 4/4.

Fresh-nonce negative controls: 0/12 successful XOR recoveries; 0 accepted GCM forgeries. Bit-flip rejection is recorded separately from forged-tag acceptance.

## Stability window and limits

The source score cutoff is frozen from separate healthy paths using a maximum over all monitored times and score components. Whole-path false alarms are measured on the evaluation controls. A nonce guard acts at the attempted policy violation; its prevention is not forecasting. Null events are administratively right-censored. Kaplan–Meier curves and exact common-horizon binomial intervals are recorded without extrapolating censored tails.

Alarm-to-compromise differences are conditional descriptive windows. A known seed is insecure immediately under the threat model, even before an observable collision. Identical byte histories can accompany different attacker knowledge or different future resets: no output-only universal warning time is identifiable. Birthday-risk curves are analytical results conditional on IID uniform finite support.

The initial suite demonstrates infrastructure failures and scoped prevention. It provides no break of AES or ChaCha20, no certified entropy estimate, no validated physical TRNG, and no deployment-wide forecast guarantee. See the frozen protocol, per-replicate JSON, source hashes, audit, and security limitations.
