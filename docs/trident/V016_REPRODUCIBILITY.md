# Reproduction

Branch `feat/v016-trident`, isolated private clone of sealed recovery source705a77e. No remote push. Existing v0.12/v0.13/v0.14 archives and outside worktrees not modified.

Python environment: standard Python 3 virtual environment (`python -m venv .venv`). Dependencies: NumPy, SciPy, Matplotlib, pytest plus package dependencies.

From `repo/`:
```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m pytest -q tests/trident
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -m rnglab16.campaign ../new-reproduction
```
Copy the current `src/rnglab16/*.py` into `../new-reproduction/source/` before running the adversarial script. Then:
```sh
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/trident_analyze.py ../new-reproduction
PYTHONPATH=src OPENBLAS_NUM_THREADS=1 python scripts/trident_audit.py ../new-reproduction
```
Report generator references the canonical confirmation-run004 path; adjust RUN for a different reproduction output. Config/source hashes must match for resume. Existing checkpoint JSON is reused; results are never appended to a different-source run. SeedSequence version/NumPy version recorded; cross-version bit-identical replay not promised. Failed and partial attempts are not pooled.

Measured execution:
{
  "elapsed_seconds": 63.33819212200001,
  "peak_rss_KiB": 217508,
  "python": "3.14.7",
  "numpy": "2.5.3",
  "platform": "Linux-7.2.3-arch1-3-x86_64-with-glibc2.44",
  "config_sha256": "95b928faa4f04b16998c9965cf07c2e48efe1c0526d9c21da89f7bbb06cf9bc0",
  "completed": true
}

Each compressed NPZ contains complete sequential observable traces; JSON records include unit counts, base seed, endpoints and artifact paths. SOURCE and SHA256SUMS identify all delivered artifacts. RSS is process peak, not aggregate host or energy consumption.


The extra saved-trace audit can be repeated with `python repo/scripts/trident_trace_checks.py` from the package root. It reads existing traces and generates bootstrap reductions, not new Monte Carlo trajectories. All five figure PNGs were visually inspected; text is black Noto Serif and data are colored.
