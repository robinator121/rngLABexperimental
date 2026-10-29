# v0.14.1 correction candidate

This private remediation branch corrects implementation and evidence semantics. It is not a confirmed universal research platform, standards implementation or production monitoring release. Original v014 results remain unchanged in their original worktree and Git history. Superseded v014 reports in this checkout carry warning banners.

- **PIMC:** `M/(2β)−MΣΔx²/(2β²)+mean(V)` replaces the incorrect spring coefficient. CV is correct. Finite-bead and continuum targets are distinct. Existing v012 IAT/ESS/batch diagnostics are reused; short Metropolis chains still under-cover. A separately labeled exact Gaussian reference sampler does not conceal that failure.
- **Forecasting:** evaluate pre-event landmarks only; censored bounds beyond observed followup remain unresolved; no-prediction errors/calibration are null. Heuristic bounds are not calibrated95% claims. A development-fit survival baseline and IPCW Brier scorer explicitly require independent trajectories/censoring/support. Dedicated reverse KM handles event-first discrete ties. Repeated landmarks are not independent trials.
- **Sources:** one byte-consumption/hash boundary, provenance format2. Adapters implement `_read_bytes`; direct/float/bit/int APIs hash emitted bytes once. The synthetic Markov fault now reaches emitted bytes consumed by floats, with rho explicitly scoped to contiguous bits. Short bit calls discard padding; this is documented rather than labeled contiguous sampling. VN buffers contiguous bits and reports acquired versus examined bits, with a finite extraction budget.
- **Conditioning:** `XORParityConditioner` is the honest API name; the old import is a compatibility alias. Prefix and rolling constant-output diagnostics reject known counterexamples but do not establish independence, entropy or security. Explicitly disabling diagnostics never changes that assurance status.
- **Physical claims:** PRNG emulators cannot be marked live hardware. Poisson dead-time modes now implement distinct mathematical laws, still synthetic only. Standards catalog entries are limited diagnostics/gaps; false full-compliance reports are superseded.
- **Monitoring:** frozen whole-chain/trajectory rank calibration; no uncalibrated VMC flag presented as reliable detection. Censored lead times are undefined; active gaps, pre-fault false alarms and post-onset detections are separate. Held-out combined monitoring was not universally superior.

Same original-source/forecast executable contracts: **0/9 before →9/9 after**. Physics exact-reference and held-out validation, preserved false positives, short-chain undercoverage, and all protocols/source hashes are under `../../evidence/physics`. Forecast v1's incorrect censor-tie fitting is preserved and invalidated explicitly; corrected v2 uses entirely fresh seed namespaces. No outcome was relabeled confirmation.

Headless API: `python -m rnglab14.api` reads one bounded `rnglab.headless.v1` JSON request from stdin, performs a known deterministic operation, and returns JSON. No provider/server/daemon starts. For example:

```json
{"schema":"rnglab.headless.v1","operation":"harmonic_reference","parameters":{"beads":16,"beta":2.0}}
```

Historical scientific v012 code/results are preserved, not retroactively repaired. Finance uses its separate namespace and independently corrected uncertainty/contracts. See the recovery root for gate status, finite local verification proof, accepted candidate review queue and unresolved questions. Candidate acceptance is not deployment or scientific confirmation.
