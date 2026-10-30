> Publication note: this is the completed local experiment report. Statements about files not yet pushed or an isolated patch describe that run before publication. The repository now includes the optional patch. Full raw traces are not distributed; see [data availability](../DATA_AVAILABILITY.md).

# Phase II environment and reproduction resources

The dedicated interpreter is `LOCAL_WORKSPACE/work/phase2-env/venv/bin/python`. `requirements-lock.txt` records every installed distribution. Principal versions are CPython 3.14.7, NumPy 2.5.3, SciPy 1.18.1, Matplotlib 3.11.2, cryptography 50.0.2, pytest 9.1.1, Numba 0.68.1 and llvmlite 0.50.0. The host has 16 logical CPUs; the scientific driver uses one process and a Numba CPU kernel, with no GPU or MPI execution. BLAS/OpenMP thread limits are 1.

Verification totals: unchanged repository 238 passing tests; optional centroid patch 257 passing full-suite tests; consolidated Phase II suite 35 passing tests (9 standalone, 14 fault, 3 baseline, 5 RUS, 4 analytic); independent saved-artifact audit 132 passing checks. These are distinct checks and are not added together as independent statistical evidence.

The 16 saved condition archives occupy1.661 GiB compressed. Summed per-condition generation time is137.19 seconds, including archive writing but excluding analysis. The complete fault campaign reports102.38 seconds. Timings are observations on this host, not performance guarantees. Per-condition timing and size values are in `experiment-resource-snapshot.json`.

A confirmation condition's raw 1024×8192×3 float64 array occupies 192 MiB before compression. A fault condition's raw 512×10240×3 array and equally sized source array occupy 120 MiB each. FFTs, observable transforms and analysis copies add memory; peak resident memory was not measured. The full experiment output tree currently occupies1.686 GiB. Preserve enough space for both trace archives and expanded analysis arrays.

To inspect saved results without rerunning simulations, run the independent validator from the distributed experiment folder with this interpreter and choose a fresh output path:

```sh
LOCAL_WORKSPACE/work/phase2-env/venv/bin/python \
LOCAL_WORKSPACE/outputs/phase2-20261030/experiment/validate_artifacts.py \
--output /absolute/path/to/new-validation.json
```

The validator pins the distributed confirmation/fault protocol hashes. Its first completed result is already stored in `experiment/artifact-validation.json`. It refuses to overwrite its chosen output or log. It neither imports the generation driver nor modifies source, protocol or trace files.

For a fresh simulation replay, use a separate empty result tree with copies of the frozen generation scripts and protocol JSONs. The drivers refuse to overwrite existing traces. Run the frozen healthy confirmation before the fault campaign; the fault driver verifies that its prerequisite healthy-decision hash and status match. Do not refreeze or edit an existing protocol to make a failed check pass. Saved manifests record all seed domains, numerical parameters, versions, code hashes, trace checksums and generation times. A new research design should receive a new protocol, stage namespace and result directory.

The consolidated Phase II run passed all 35 checks in 7.57 seconds after all campaign scripts were frozen. Its complete log, JUnit XML and command metadata are `pytest-phase2-consolidated.log`, `pytest-phase2-consolidated.xml` and `pytest-phase2-consolidated-run.json`. The upstream repository suite was excluded from this consolidated run because it had already been executed.
