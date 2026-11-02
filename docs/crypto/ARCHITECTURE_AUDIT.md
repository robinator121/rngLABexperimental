# Repository audit for RNGLab BLUE

Audit base: `8501f9a850eb96b667dfcbcf99f479943a7cbfbc` (remote default branch retrieved 2026-11-02).
Work occurs in a fresh clone on `feature/crypto-blue`. An earlier local checkout had unrelated staged/unstaged edits and was left untouched. Synced project sources are read-only. No physics, finance or acquisition code is changed.

The actual baseline was **257 passed in 4.58 s**, Python 3.14.7. The host default interpreter initially lacked NumPy; a new `.venv` was installed from the repository's declared dependencies plus pytest and matplotlib. Exact pins and the unmodified baseline JUnit/log are in `research/crypto-blue/evidence/`. This local run does not establish Python 3.11/3.12 compatibility; CI covers those versions.

| Area inspected | Implemented behavior and limitation | BLUE decision |
|---|---|---|
| `rnglab11/sources.py` | Bounded reads, consumed-output hashes, finite files, explicit exhaustion; metadata is not entropy evidence | Preserve conventions for budgets, provenance and explicit errors |
| `rnglab11/drbg.py` | HMAC-DRBG SHA256 with update, reseed, request limits, opted-in state export; length checks do not assess entropy | Reuse this implementation; independently exercise official CAVP no-reseed and reseed vectors |
| `rnglab11/crypto_workloads.py` | Counts nonce/block duplicates; derives Ed25519 seeds, signs and verifies; no AEAD confidentiality/integrity experiment | Add actual AEAD consumers and scoped misuse demonstrations |
| `rnglab11/monitoring.py` | Prospective monitors with configuration/digest-bound calibration and separate histories | Preserve frozen whole-path calibration principle |
| `rnglab13/sources`, `standards` | Distinct source interface; custom HMAC counter reseeding is not the same mechanism as the HMAC-DRBG wrapper | Do not label the custom counter generator as SP800-90A HMAC-DRBG |
| `rnglab13/standards/entropy_90b.py` | Despite names/docstrings, compression is empirical block Shannon entropy multiplied by 0.85; Markov is maximum empirical one-step transition; collision uses an approximate recurrence formula; MCV ignores its alpha argument | No cryptographic entropy credit from these heuristics. Their minimum is not an established conservative bound on attacker uncertainty |
| `rnglab13/math/impossibility.py` | Identical-history/different-future construction states identifiability limitations | Extend to unknown attacker state and downstream restart regimes; abstain from universal forecasting |
| `rnglab14/sources/base.py` | `read_bytes` is one audited accounting boundary; `SourceManifest` has timestamps; byte/bit/float interfaces | Subclass `RandomSource`; fix replay timestamp label; explicitly document call schedule |
| `rnglab14/sources/faults.py` | Persistent binary Markov path when rho is nonzero; otherwise bias shifts/clips bytes, not Bernoulli bit bias | Implement explicitly parameterized bit faults in isolated BLUE module, with no modification to legacy fault semantics |
| `rnglab14/statistics/multitest.py`, `forecast/survival.py` | Holm correction; KM with explicit censoring support; IPCW assumptions documented | Reuse Holm and KM; no IPCW needed for complete common-horizon follow-up |
| `rnglab14/monitoring/dual_monitor.py` | Source/application distinction; trace helper requires fault and failure onset together | BLUE event schema permits fault without observed failure; source alarm, policy violation and compromise remain separate |
| `rnglab15/ingest.py`, `assessment.py` | Explicit schema, acquisition order retained, streaming timing, RCT/APT; descriptive screens explicitly disclaim entropy validation | Reuse ingestion/timing and binary health screens with stated assumptions |
| `rnglab15/conditioning.py` | Finite raw/VN/hash treatments; conditional Toeplitz budget; entropy credit zero without external assumptions | Follow zero-credit physical replay policy; no physical key generation |
| `rnglab17/endpoints.py` | Rank-based frozen threshold and exact binomial rates; paired-trajectory analysis patterns | Reuse rank threshold and two-sided Clopper–Pearson rates with input validation |
| `research/phase2-20261030` | Separate calibration/evaluation namespaces, paired replicas, frozen protocols, healthy gates, logs, hashes, and explicit data omissions | Adopt these conventions; distinguish replay from independent replication |
| CI/packaging | Python >=3.11, NumPy/SciPy/cryptography already declared; CI installs scientific dependencies and runs `tests` on 3.11/3.12 | Add isolated crypto extra, CLI and artifact smoke/replay checks; preserve existing test job |

The Phase II RUS report is particularly useful: an identical alarm history has different conditional future failure laws. BLUE cannot infer a secret seed's exposure from identical output bytes either. An observable alarm-to-collision interval is therefore not a guaranteed interval of continued cryptographic security.

New package: `sources/` (replay/OS/capture/physical audit), `faults/` (scheduled injection), `primitives/` (library consumers and contained algebra), `experiments/` (validated configs and paired campaign), `monitoring/` (diagnostics/policy), `risk/` (censoring and analytical birthday risk), `reporting/` (immutable run artifacts). Existing public APIs remain unchanged.

Physical data availability: the repository includes an assessment report but not the original full detector export. No matching raw waveform file was found in the current project's `sources/` or output inventory. Historical report measurements are attributed, not presented as newly reproduced. See `PHYSICAL_INTEGRATION.md`.
