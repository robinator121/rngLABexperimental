# Reproducing v0.18 Phase 2

Environment: the Phase II lock (`research/phase2-20261030/environment/requirements-lock.txt`). The run used Python 3.14.9, NumPy 2.5.3 and SciPy 1.18.1. CPU only.

## 1. Tests

```bash
PYTHONPATH=src .venv/bin/pytest -q
```

The `tests/v018/test_pimc_default_golden.py` test pins the sampler's pre-v0.18 default behaviour bit for bit. `tests/v018/test_api_validation_scripts.py` runs the whole pipeline on 4 replicas per condition, under the smoke namespace, in a temporary folder.

## 2. Generation (4 cores; the recorded run took 09:07:17–09:55:49 UTC per `results/run.log`; traces about 3 GiB per `design.json`)

Check out `ff948e2` or later. `run` refuses to start unless `protocol.json` and every file listed in `manifest_pregeneration.json` match their recorded SHA-256. It never overwrites a finished block; a restart skips completed blocks.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-api-validation/run_campaign.py run --workers 4
```

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-api-validation/run_campaign.py replay
```

To reproduce into a fresh folder, move the existing `traces/` and `results/` aside first. The run is deterministic: every block should reproduce the SHA-256 recorded in `results/run_manifest.json`, and every per-replica CSV the hash recorded there.

## 3. Analysis, recount, tables, figures

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-api-validation/analyze.py
```

```bash
PYTHONPATH=src .venv/bin/python research/v018-api-validation/recount.py
```

```bash
.venv/bin/python research/v018-api-validation/report_tables.py
```

```bash
MPLBACKEND=Agg .venv/bin/python research/v018-api-validation/plot.py
```

```bash
PYTHONPATH=src .venv/bin/python research/v018-api-validation/posthoc_bias_check.py
```

`analyze.py` writes `results/gates.json` and `results/secondary.json`. `recount.py` writes `results/recount.json`, which should report all 21 comparisons in agreement. `posthoc_bias_check.py` is a post-hoc check that was not pre-registered.

## 4. Design numbers and smoke (no campaign data)

```bash
PYTHONPATH=src .venv/bin/python research/v018-api-validation/design_calculations.py > research/v018-api-validation/design.json
```

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-api-validation/smoke/smoke_timing.py > out.json
```

The smoke uses the separate namespace `rnglab-v018-api-validation-smoke`, and its timings vary from machine to machine.
