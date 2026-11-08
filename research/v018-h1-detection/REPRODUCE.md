# Reproducing v0.18 Stage C

Environment: the Phase II lock (`research/phase2-20261030/environment/requirements-lock.txt`). The run used Python 3.14.9, NumPy 2.5.3 and SciPy 1.18.1 (`results/run_manifest.json`). CPU only. The scripts that generate or replay runs refuse to start unless `OPENBLAS_NUM_THREADS=1` and `OMP_NUM_THREADS=1`.

## 1. Tests

```bash
PYTHONPATH=src .venv/bin/pytest tests -q
```

`tests/v018/test_h1_lib.py` checks the added SP 800-22 tests against the standard's worked examples, and checks the fault adapters, the sweep clock and the GFSR. `tests/v018/test_h1_campaign.py` runs the whole pipeline on a tiny campaign under the smoke namespace.

## 2. Pilot and grid (already frozen; `pilot/summary.json`, `grid.json`)

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-h1-detection/run_pilot.py run --workers 8
```

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-h1-detection/run_pilot.py summarize
```

```bash
.venv/bin/python research/v018-h1-detection/make_grid.py
```

## 3. Generation (8 workers; the recorded run took 22:30:59–00:21:25 UTC per `results/run.log`)

Check out `0aca45a` or later. `run` refuses to start unless `protocol.json`, `grid.json` and every file listed in `manifest_pregeneration.json` match their recorded SHA-256. From `456ebcd` on, `analyze.py` and `recount.py` differ by Amendment 01, so generate from `0aca45a` to pass that check, then return to the current commit for step 4. It never overwrites a finished block; a restart skips completed blocks. To reproduce into a fresh folder, move the existing `traces/` and `results/` aside first. Every cell file should reproduce the SHA-256 recorded in `results/run_manifest.json`.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-h1-detection/run_campaign.py run --workers 8
```

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-h1-detection/run_campaign.py replay
```

## 4. Analysis, recount, Amendment 01 comparison, tables, figures

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-h1-detection/analyze.py
```

```bash
PYTHONPATH=src .venv/bin/python research/v018-h1-detection/recount.py
```

```bash
PYTHONPATH=src .venv/bin/python research/v018-h1-detection/amendment01_compare.py
```

```bash
.venv/bin/python research/v018-h1-detection/report_tables.py
```

```bash
.venv/bin/python research/v018-h1-detection/report_numbers.py
```

```bash
MPLBACKEND=Agg .venv/bin/python research/v018-h1-detection/plot.py
```

- **`analyze.py`** writes `results/gates.json` and `results/secondary.json`. The secondaries read `traces/`.
- **`recount.py`** writes `results/recount.json`, which should report 319 of 319 comparisons in agreement.
- **`amendment01_compare.py`** runs the superseded `analyze.py` from `0aca45a` and writes `results/amendment01_comparison.json`.
- **`amendment01_record.py`** wrote the pre-analysis record `results/amendment01_preanalysis.json`. It is meaningful only before an analysis has run, so do not re-run it after one.
