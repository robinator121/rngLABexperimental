# Scientific review and correction record

The initial implementation/campaign checkpoint is commit `41e7c4c`. This review was performed by the same implementing agent, with independent arithmetic and published vectors; it is not an external cryptographic audit.

## Corrections before final delivery

1. **Fault attribution:** initial reporting could label a post-onset alarm during a downstream fixed-nonce/reset fault as a source “detection.” The raw bytes in those arms are exactly the healthy paired bytes. Reviewed reporting marks every such alarm as a healthy-source false alarm, gives no source detection latency, and credits no eligible warning. Incidental quarantine prevention is retained separately from warning effectiveness.
2. **Actual source-gated execution:** the initial source-control outcome was computed from the deterministic emitted prefix. The reviewed code also executes an independent AEAD consumer only before quarantine and derives its repeated emissions from that emitted list. Counts did not change; this verifies the operational gate rather than relying only on its counterfactual arithmetic.
3. **Layer provenance:** the bounded-support adapter's default placement describes DRBG initialization. The nonce experiment now explicitly overrides the target layer to raw nonce seed material; each event row records its actual fault layer.

No fault intensities, sample sizes, seed namespaces, thresholds, hypothesis tests or acceptance criteria changed. Reviewed campaigns replay the same independent units. Archived initial runs remain intact and are reproducible from the checkpoint; **do not add their sample counts** to reviewed runs.

## Falsification and limits

- Published vectors test HMAC-DRBG before/after reseed, both AES-GCM key sizes, ChaCha20-Poly1305, RFC6979, Ed25519 and HKDF. They do not constitute exhaustive implementation or standards validation.
- Fresh-nonce GCM controls reject the constructed forgery. Unaided bit flips are separately rejected. This prevents confusing normal corruption rejection with integrity protection under nonce misuse.
- Changed-message ECDSA signatures and a recovered-key new signature are verified by the maintained library. Same-message deterministic signatures are not misclassified as misuse.
- Raw hashes/alarms are equal for healthy and downstream-fault arms; no source detector can distinguish these paths from those inputs alone. Known-seed and hidden-seed output bytes can likewise be identical while attacker uncertainty differs.
- Kaplan–Meier values are independently recomputed as empirical survival under common administrative censoring; missing events are not treated as failures or infinite safe lifetime. Guard and source-gate operation counts are checked for consistency.
- Confirmation's healthy false-alarm gate failed; it was neither suppressed nor retuned. The single conditioner/expander confirmation failure (Holm-adjusted p≈.05096) was retained. Sensitivity outcomes do not override either decision.
- A post-hoc full 256-seed diagnostic enumeration and full 8/12/16-bit nonce-map enumeration are labeled secondary exact finite-domain checks. They clarify sampling/support assumptions and are not substituted for preregistered evidence.
- Serial/drift cases without collision remain negative finite-horizon observations. Their lack of measured compromise is not a general proof of security. Counter-state assumptions, known versus unknown reseed inputs, and need for application rekeying are explicit.
- Physical audit results from the earlier detector report are attributed, not reproduced or promoted to entropy validation. No real-data security claim is made.

The independent result-accounting script passed **28,764 checks** on reviewed confirmation and **14,463** on sensitivity. These are structural/arithmetic checks, not thousands of independent cryptographic experiments. Exact scientific-result SHA256 values are:

- Confirmation: `898f92db249eca38c3faf1585a3b38980b332bbed619780fa2d7ca7d28d07df3`
- Sensitivity: `1a90e4705f6ff43151e20351f8a0dc5c16f6cc789891797ac9eb7e87af71214d`

The core suite passed **322 tests**: 257 unchanged baseline tests and 65 new crypto checks. A clean installed-wheel smoke run outside the source tree also completed. CI adds CPU-only smoke/replay verification on its existing Python 3.11/3.12 matrix; remote CI execution is not claimed.

An initial attempt to collect the separate Phase II research tests under the root pytest configuration failed on standalone import paths. This is retained in `evidence/legacy-research.log`; it is an invocation error, not counted as a passing test run. The project's documented isolated replay preparation/environment is used for the subsequent research regression check, recorded separately.

That documented isolated run subsequently passed **35/35 Phase II checks** in 2.99 s, using its pinned upstream comparison source and dependency environment. These are separate from the 322 core tests. Confirmation and sensitivity reviewed scientific JSON both reproduced byte-for-byte. The clean installed-wheel smoke JSON matched the current source-checkout smoke JSON, and `pip check` reported no broken requirements.
