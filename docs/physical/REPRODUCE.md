# Reproduce the recorded v0.15 experiments

The Python 3.14.7 version pins are in `requirements-physical.lock`. Use a fresh
virtual environment; install those pins and the project with no build isolation.
The registry records exact dependency versions, code SHA-256 values, Git parent,
input hashes, the complete configuration, and every output hash. The source bundle
and patch contain no measured raw dataset. You must have the listed local inputs.
The content-addressed private vault is not a public data archive or tamper-proof
storage; mode 0444 prevents accidental editing, and use-time hashes detect changes.

Run from the repository root. Choose a fresh result directory (here `/tmp/v015-repro`):

```sh
mkdir /tmp/v015-repro
python -m pytest -q -p no:cacheprovider
python -m rnglab15 inventory /home/robi/rnglab-data /tmp/v015-repro/inventory.json
python -m rnglab15 run examples/physical/local-20261019.json /tmp/v015-repro/measured-run-01 --vault /tmp/v015-private-vault
# Copy the delivered archived OS reference file; do not resample it for exact replay.
python -m rnglab15 simulate /tmp/v015-repro/synthetic-campaign-01 --trials 32 --reference /path/to/delivered/outputs/synthetic-campaign-01/os-reference.bin
python scripts/run_physical15_replays.py /tmp/v015-repro
```

This reproduces feature values, treatment bytes and application outcomes. Registry
paths, elapsed timing, platform and Git parent reflect the new run and therefore
need not be byte-identical. The synthetic result payload is byte-reproducible in
the same dependency/runtime environment from its archived OS bytes. New OS bytes
create a new experiment, not the same seed-based replay. Physical segmentation
uses finite disjoint file ranges; role seeds are paired across treatments.

Optional report plots use `scripts/report_physical15.py`. The report expects a
resource-usage JSON with measured `elapsed_seconds` and Linux `max_rss_kib_linux`;
copying the old timing is not a new benchmark. The delivered plots are standalone
PNG files. RCT/APT alarms and hashing are diagnostics, never a secure-output gate
in this research lane. Failed input validation retains a failed registry and
partial evidence. Pick a fresh run directory after resolving a schema mismatch.

The historical v0.14 baseline had 120 passing tests before any edit. The v0.15
suite additionally covers schema failures, timestamps/resets, missing sensor
values, raw hash preservation, interrupted bounded capture, immutable snapshots,
source mutation, conditioner vectors/tails, exact APT cutoffs and RCT boundaries,
latched alarms, finite-source exhaustion, entropy-budget refusal, NIST-tool
blockers/invocation bounds, all four model bridges and full OS-fixture replay.

The old project distribution name `rnglab-research-v012` is retained for compatibility;
its version is now `0.15.0.dev1`, with an additive `rnglab-physical` command. No
legacy research results are recomputed or overwritten by these commands.
