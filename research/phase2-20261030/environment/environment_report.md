> Publication note: this is the completed local experiment report. Statements about files not yet pushed or an isolated patch describe that run before publication. The repository now includes the optional patch. Full raw traces are not distributed; see [data availability](../DATA_AVAILABILITY.md).

# RNGLab Phase II environment and repository verification

Date: 2026-10-30 (Europe/Budapest). Source commit: `dd96adba23b8100aba9147d71fffc63932fc6b0d`.

The previous missing-dependency blocker is resolved in a dedicated virtual environment. The complete unchanged repository test suite passes: **238 passed in 9.63 seconds**. No repository defect was demonstrated by that baseline run. A subsequently authorized optional harmonic centroid update was added only to the isolated copy; its final complete suite passes **257 tests in 9.87 seconds**. The remote CI failure has a separate, explicit account-level cause: GitHub reports that the jobs were never started because the account is locked due to a billing issue.

## Isolation and provenance

- Original preserved snapshot: `LOCAL_WORKSPACE/work/rnglab-main-snapshot`.
- Isolated test copy: `LOCAL_WORKSPACE/work/phase2-env/repo`.
- Virtual environment interpreter: `LOCAL_WORKSPACE/work/phase2-env/venv/bin/python`.
- Original snapshot and previous outputs were not edited. Git status in the isolated copy was clean before and after baseline testing. The later optional patch changes only the two model files and adds one regression-test file listed in the manifest.
- No applicable `AGENTS.md` was found in the workspace ancestors or inside the source snapshot.
- Tests import directly from the isolated copy's `src` directory. The package has not been installed in editable mode, so experiment commands should set `PYTHONPATH` or explicitly add the intended source directory.

## Installed environment

The repository's `docs/calibrate/requirements-reproduction.txt` was installed without altering its five direct pins.

| Component | Version |
|---|---:|
| CPython | 3.14.7 |
| NumPy | 2.5.3 |
| SciPy | 1.18.1 |
| Matplotlib | 3.11.2 |
| cryptography | 50.0.2 |
| pytest | 9.1.1 |
| Numba, optional accelerator | 0.68.1 |
| llvmlite | 0.50.0 |
| pip | 26.2.1 |

Numba and llvmlite were installed after the full suite, without changing NumPy or any other pinned scientific package. A compiled float64 loop produced the expected exact numeric result (285.0), confirming that compilation works on this interpreter. This is a compiler smoke check; any accelerated experimental kernel still needs its own numerical comparison against an independent implementation. `pip check` reports no broken requirements. There was an optional NumPy configuration display warning about PyYAML; it does not affect imports, numerical execution, or the test suite.

The host reports `Linux-7.2.3-arch1-3-x86_64-with-glibc2.44` and 16 logical CPUs. The complete installed distribution list and interpreter details are recorded in `environment-manifest.json`, and all versions, including transitive packages, are frozen in `requirements-lock.txt`. NumPy's build and BLAS configuration is recorded in `numpy-build-config.txt`.

## Test result

Command, executed from the isolated repository:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg \
PYTHONPATH=LOCAL_WORKSPACE/work/phase2-env/repo/src \
LOCAL_WORKSPACE/work/phase2-env/venv/bin/python \
-m pytest tests -v --durations=20 \
--junitxml=LOCAL_WORKSPACE/outputs/phase2-20261030/environment/pytest-initial.xml
```

- Result: exit code 0, 238 passed, no failures or skipped tests reported.
- Pytest's test duration: 9.63 seconds; outer process wall time: 10.478 seconds.
- Scope includes physical-source parsers, source consumption/provenance, harmonic-oscillator references, MCMC uncertainty, VMC masking, forecasting/censoring, standards adapters, and financial prototypes.
- This is the current repository's automated regression suite. Passing it does not establish nominal confidence-interval coverage or stationarity across independent simulation campaigns; Phase II scientific experiments address those separately.
- This local run uses CPython 3.14.7, not the remote workflow's 3.11/3.12 matrix. Compatibility with those two interpreters was not newly exercised here.

The full log is `pytest-initial.log`; machine-readable results are `pytest-initial.xml` and `pytest-initial-run.json`.

## Remote CI diagnosis

The [latest workflow run](https://github.com/robinator121/rngLABexperimental/actions/runs/40931415718) matches the audited commit. The public GitHub API lists two failed jobs, `test (3.11)` and `test (3.12)`, created on 2026-10-29 at 18:30:14 UTC and completed at 18:30:17 UTC. Both jobs have no executed steps and runner id 0. Their check annotations explicitly state:

> The job was not started because your account is locked due to a billing issue.

This establishes an account infrastructure blocker rather than evidence that tests failed. The repository owner can resolve the account issue in GitHub billing/account settings and rerun the workflow. No GitHub account, workflow, repository, or billing setting was changed during this work.

The companion notice says the `ubuntu-latest` label begins migration to Ubuntu 26 on October 19, 2026. This is informational and is not the annotated failure cause. Raw responses are retained in `github-actions-runs.json`, `github-actions-jobs.json`, `github-check-runs.json`, and the two `github-check-*-annotations.json` files.

## Recreate the environment

Use Python 3.14.7 on a compatible Linux host, create an isolated environment, and install the saved full lockfile:

```sh
python3.14 -m venv work/recreated-rnglab-env
work/recreated-rnglab-env/bin/python -m pip install \
-r outputs/phase2-20261030/environment/requirements-lock.txt
work/recreated-rnglab-env/bin/python -m pip check
```

The lockfile fixes versions but does not pin downloaded wheel hashes. The saved installation logs identify the wheel filenames used. The source commit and interpreter version must also be retained when comparing scientific results.

## Deliverable index

- `environment-manifest.json`: interpreter, platform, package versions, Git provenance, compiler smoke result.
- `requirements-lock.txt`: full exact installed distribution versions.
- `dependency-install.log`, `numba-feasibility.log`, `numba-install.log`, `venv-create.log`: installation evidence.
- `pip-check.log`, `numpy-build-config.txt`: dependency and numerical-library configuration evidence.
- `pytest-initial.log`, `pytest-initial.xml`, `pytest-initial-run.json`: complete test evidence.
- `github-*.json`: remote CI metadata and causal annotations.
- `pimc-centroid-heatbath.patch`, `centroid_update_patch.md`: optional harmonic centroid update and its review note.
- `pytest-centroid-targeted.*`, `pytest-centroid-full.*` and their run JSON files: added-update verification.
- `patch-apply-check.log`: read-only application check against the pristine commit.
- `SHA256SUMS`: artifact integrity inventory, excluding itself.

The optional research improvement is supplied as `pimc-centroid-heatbath.patch`, with design, API, assumptions and verification details in `centroid_update_patch.md`. It passed 51 targeted physics tests, then the final full suite passed 257 tests. Its patch passes `git apply --check` against the pristine audited commit. The original snapshot remains unchanged. This patch does not address the independent GitHub billing blocker.


## Fault experiment review and saved-artifact audit

The frozen fault driver passed 14 independent small-kernel/window/event checks. The completed campaign then passed 132 read-only artifact checks, including original Phase I preservation, all 16 archive checksums, seed derivations, independent primary CV BM512 arithmetic and full-replica prefault byte equality. Details are in `fault_study_review.md` and `reproduction_resources.md`; tests and audit evidence are indexed by the environment manifest.

The final consolidated Phase II suite passed **35 tests in 7.57 seconds**: standalone 9, fault logic 14, baseline checker 3, RUS 5, and analytic dependence 4. The upstream suite was excluded from that consolidated run. Full evidence is saved under `pytest-phase2-consolidated.*` and its run metadata JSON.
