# Proposed research paper outline

Working title: **From Randomness Diagnostics to Cryptographic Preconditions: Measured Infrastructure Failures and Detection Limits**.

1. Motivation: distinguish simulation reliability, distributional randomness and attacker-conditional unpredictability. State an infrastructure question rather than a primitive-break claim.
2. Threat model and causal layers: source, conditioner, deterministic generator, consumer and defensive controls; synthetic-only scope.
3. Repository/assurance audit: reusable RNGLab components, unvalidated entropy heuristics, corrected interpretation and fixed-vector checks.
4. Experimental methods: frozen design, paired trajectories, independent calibration, bounded attacks, uncertainty, censoring and evaluation splits.
5. Measured results: weak-seed search, diagnostic blind spots, AEAD misuse, ECDSA algebra/deterministic controls and recovery.
6. Security stability windows: distinguish insecurity, alarm, policy breach and demonstrated compromise; explain selection/censoring and the identical-observation impossibility examples.
7. Defensive effectiveness: source monitoring versus precondition enforcement, durable counters, budget policies and quarantine; false alarms and overhead.
8. Physical-source integration: historical report attribution, missing provenance, read-only boundary and proposed measurement program.
9. Limitations and falsification: model assumptions, finite horizons, negative results, no certification or practical production-key recovery; independent validation still needed.
10. Reproducibility and availability: public seeds, environment, source/artifact hashes, complete per-replica evidence, tests and commands.

The initial campaign supports a bounded systems demonstration and measured detection study. A broad forecasting or hardware entropy-security paper requires additional independently reviewed experiments; the current results must not be stretched into those claims.
