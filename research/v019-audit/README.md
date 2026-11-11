# v0.19 audit: independent verification of the Stage D pilot

Code in this directory shares nothing with the pilot's own pipeline (`research/v019-d-pilot/`), the Stage C library or the package's uncertainty code. It reads only stored outputs and the frozen `protocol.json`. It generates no pilot data; its engineering runs use the namespace `rnglab-v019-audit-smoke`.

| File | Purpose |
|---|---|
| `recount.py` | independent formulas for every pilot quantity (imports the standard library, NumPy and SciPy only) |
| `replay.py` | regenerates a whole run from its seed labels with its own sweep loop |
| `verify.py` | compares stored runs with their recomputation; recounts the campaign; writes a JSON report of discrepancies |
| `integrity.py` | atomic run commits, classification of partial and corrupt runs, manifests |
| `preflight.py` | checks before generation; `golden.json` holds the digests it compares |
| `bench.py`, `audit_numbers.py` | timings, memory, and the numbers quoted in `docs/v019/IMPLEMENTATION_AUDIT.md` |
| `results/` | `benchmark.json`, `audit_numbers.json` |

## Commands

```bash
PYTHONPATH=src .venv/bin/python -m pytest tests/v019 -q
```

```bash
PYTHONPATH=src .venv/bin/python research/v019-audit/preflight.py
```

```bash
PYTHONPATH=src .venv/bin/python research/v019-audit/verify.py --root research/v019-d-pilot/smoke --namespace rnglab-v019-d-pilot-smoke --replicas 2 --deep all --no-campaign --out verification_smoke.json
```

```bash
PYTHONPATH=src .venv/bin/python research/v019-audit/verify.py --root research/v019-d-pilot --namespace rnglab-v019-d-pilot --deep sample --workers 6 --out results/verification.json
```

`--deep all` regenerates every run (about 2 s each); `--deep sample` regenerates the protocol's I-4 sample of 25.

Tests that need the measurement layer or the Stage C traces skip unless `RNGLAB_DPILOT_MEASURE` (path of `measure.py`) or `RNGLAB_STAGEC_TRACES` (the `traces` directory) is set. `RNGLAB_AUDIT_FULL=1` adds a full-length storage round trip.

Findings: `docs/v019/IMPLEMENTATION_AUDIT.md`.
