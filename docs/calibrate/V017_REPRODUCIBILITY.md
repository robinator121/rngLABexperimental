# Reproduction
Private isolated branch feat/v017-calibrate, cloned from sealed v0.16 commit21f93873db149cc2939d44ad6a344deea4cbbbb9. No historical result/source modified, no remote push, no service restart. Existing historical process was inspected read-only and left running.

Environment: Python >= 3.11 with NumPy >= 2.0, SciPy, Matplotlib, and pytest. CPU-only execution, two workers, one numerical thread each. Runtime/memory receipts are in execution.json and every cell. Fresh seed namespaces are SHA256-derived, and replicas use SeedSequence children. This prevents identified namespace reuse, not a proof of PRNG statistical independence.

From repo/, using that interpreter:
```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q tests/calibrate tests/trident tests/test_recovery_physics.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m rnglab17.campaign ../new-confirmation
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/calibrate_analyze.py ../new-confirmation
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/calibrate_forecast.py ../new-confirmation
```
The forecast command checks the calibrated L3/Tc gate. Cell JSON+NPZ checkpoints are preserved; source/config mismatch rejects resume. Hard wall/memory limits within a cell are not guaranteed by the runner; stated stage budget is checked at checkpoints.

Canonical read-only verification: scripts/calibrate_audit.py uses the delivered confirmation-v2/ directory; adjust RUN to verify a reproduction. scripts/calibrate_forensic.py additionally needs the separately preserved sibling v0.16 package. Main CALIBRATE forecasting uses the portable reference/v016_theory_grid.json; provenance preserves its original and executed source identity.

Exploratory controls and plots:
```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/calibrate_physics_controls.py
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/calibrate_sensitivity.py
MPLCONFIGDIR=/tmp/calibrate-mpl PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/calibrate_report.py
```
The latter scripts use canonical package paths; update those paths for a separate reproduction directory. Do not write a rerun over sealed original results. Main configuration/protocol is frozen; pilot/controls/sensitivity results are distinct classifications. Seven figure PNGs use black Noto Serif text and colored data. No PDF was requested.
