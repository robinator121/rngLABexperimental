# Reproduce RNGLab BLUE

Use Python 3.11+ and a clean checkout of `feature/crypto-blue`. The recorded local execution used Python 3.14.7; CI is configured for 3.11 and 3.12. Tests and ordinary campaigns require only CPU resources.

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install '.[crypto]'
pytest tests -q

# Fast development execution, not scientific confirmation
python -m rnglab_crypto run --smoke --output /tmp/blue-smoke-new
python -m rnglab_crypto verify /tmp/blue-smoke-new

# Registered confirmation and sensitivity configurations
python -m rnglab_crypto run \
  --config research/crypto-blue/configs/confirmation.json \
  --output /tmp/blue-confirm-new
python -m rnglab_crypto run \
  --config research/crypto-blue/configs/sensitivity.json \
  --output /tmp/blue-sensitivity-new
python -m rnglab_crypto verify /tmp/blue-confirm-new
python -m rnglab_crypto verify /tmp/blue-sensitivity-new

# Exact replay requires the same code and compatible numerical dependencies
cmp research/crypto-blue/confirmation-reviewed/results.json /tmp/blue-confirm-new/results.json
cmp research/crypto-blue/sensitivity-reviewed/results.json /tmp/blue-sensitivity-new/results.json

# Independent result-accounting and finite-support review
python scripts/crypto_review.py /tmp/blue-confirm-new /tmp/blue-sensitivity-new \
  --output /tmp/blue-independent-review.json

# Optional observational outputs; never part of deterministic results.json
python scripts/crypto_overhead.py /tmp/blue-overhead-new.json
python -m rnglab_crypto os-reference --output /tmp/blue-os-new
```

Every output directory must be new: the CLI refuses to overwrite a previous campaign. A failed run can leave its frozen protocol; choose a new directory for another run. `protocol.json` freezes the configuration, seed-derivation convention, methodology hash where available and source hashes before generation. `results.json` includes all replicate outcomes, source score trajectories, confidence intervals, decisions and counterfactual comparisons. `RESULTS.md` and four figure files are generated from those data. `manifest.json` hashes the exact artifact inventory; `verify` rejects altered/missing/extra files. This is integrity checking, not a signed provenance certificate.

The installed entry point `rnglab-crypto` is equivalent to `python -m rnglab_crypto`. An installed wheel may not include the repository methodology document; that field is then null, while package code hashes remain available. Use the checkout for complete scientific reproduction.

Replay material derives from SHA256 of a compact JSON list `["BLUE-v1", namespace, replicate, role]`. HMAC-DRBG advances per Generate/read call: changing read sizes/schedule changes streams. Within-replica sharing and separate calibration namespaces are deliberate. Reusing the same namespace is replay, not an independent experiment. Register new namespaces/configurations before new scientific evaluation.

The environment manifest records NumPy, SciPy, cryptography, OpenSSL, Python and platform. `research/crypto-blue/evidence/requirements-lock.txt` records the original installed versions; omit its editable local-path line if using it to recreate dependencies, then install the current checkout. Elapsed wall time, platform strings and OS capture metadata are intentionally not expected to replay. Exact JSON floating-point identity is verified in the recorded environment; different supported versions may differ numerically. Such differences require review, not automatic acceptance based only on rounded summaries.

Archived `confirmation-v1`/`sensitivity-v1` were produced by checkpoint commit `41e7c4c`. The reviewed runs retain the same experimental design and seeds, with corrected alarm attribution and explicit execution of the source-gated encryption branch. They are **replays of the same experiments, not additional independent samples**. The smoke run predates final reporting refinements and is development evidence only. See `SCIENTIFIC_REVIEW.md` for the exact change rationale.

Physical audits require the unavailable original file and a declared schema/provenance sidecar; see `PHYSICAL_INTEGRATION.md`. No hardware is required for the synthetic test suite.

The additional 35 historical Phase II checks use their own prepared scratch layout and pinned environment: follow `research/phase2-20261030/REPRODUCE.md`. Collecting those standalone files under the root project's importlib pytest configuration produces import errors; it is not the documented invocation. BLUE verification used the documented preparer without modifying the original research payload.
