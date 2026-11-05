# Reproduce the v0.18 Phase 1 estimator campaign

CPU only. The original run used Python 3.14.9 with the Phase II lock (`research/phase2-20261030/environment/requirements-lock.txt`: NumPy 2.5.3, SciPy 1.18.1), 3 worker processes and one BLAS thread each. Disk: about 2 GB for the Ising traces (gitignored).

## 1. Environment

```bash
python3 -m venv .venv
.venv/bin/pip install -r research/phase2-20261030/environment/requirements-lock.txt
PYTHONPATH=src .venv/bin/pytest tests/v018 -q
```

`tests/v018/test_v017_susceptibility_replay.py` regenerates the v0.17 counterexample from its original seed and must match it to relative 1e-9.

## 2. Check the frozen protocol and manifest

```bash
cd research/v018-estimators && sha256sum -c protocol.sha256 && cd ../..
```

`manifest_pregeneration.json` lists the SHA-256 of every file the campaign depends on and all 1,512 seed units. `run_campaign.py run` refuses to start if any hash differs.

## 3. Regenerate into a fresh copy

The runner never overwrites existing outputs. To replay without touching the published results, clone the repository at the commit that added `results/`, delete `research/v018-estimators/results/ar1`, `results/ising` and `traces/` in that clone, then:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-estimators/run_campaign.py run --workers 3
PYTHONPATH=src .venv/bin/python research/v018-estimators/analyze.py
PYTHONPATH=src .venv/bin/python research/v018-estimators/recount.py
MPLBACKEND=Agg PYTHONPATH=src .venv/bin/python research/v018-estimators/plot.py
```

Results are deterministic given the seeds and library versions. Compare `results/gates.json` and the per-unit NPZ files with the published ones; `results/ising/*.npz` record the SHA-256 of each trace file (`trace_sha256`).

## 4. Smoke and design numbers

```bash
OPENBLAS_NUM_THREADS=1 PYTHONPATH=src .venv/bin/python research/v018-estimators/smoke/smoke_timing.py > /tmp/smoke.json
.venv/bin/python research/v018-estimators/design_calculations.py > /tmp/design.json
```

The smoke uses the separate namespace `rnglab-v018-estimators-smoke`; timings vary by machine.

## Files

| File | Role |
|---|---|
| `PROTOCOL.md`, `protocol.json`, `protocol.sha256` | frozen design |
| `manifest_pregeneration.json` | code hashes and seeds, committed before the first draw |
| `run_campaign.py` | generation and per-unit estimators |
| `analyze.py` | frozen gate evaluation → `results/gates.json`, `results/*_table.csv` |
| `recount.py` | independent recount → `results/recount.json` |
| `plot.py` | figures → `figures/*.png` |
| `results/ar1/*.npz`, `results/ising/*.npz` | per-replica, per-chain and per-ensemble endpoints |
| `traces/*.npz` (not in Git) | Ising integer total magnetization and energy for all chains and sweeps |
