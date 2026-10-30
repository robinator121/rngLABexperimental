# Reproduce the public Phase II study

The public payload contains frozen generation scripts, protocols, seed manifests, analysis tables and figures. Large simulation trace archives are omitted. Seven lightweight Phase I originals are published under `research/phase1-20261030`; the full set of ten original Phase I inputs is not distributed because its raw trace and two physical inventories are omitted. A full replay regenerates roughly 2 GiB of Phase II traces on a CPU; additional analysis arrays require RAM and temporary space. The historical execution used one process, NumPy/SciPy and a Numba CPU kernel, with BLAS and OpenMP limited to one thread.

The published outcomes and their retained seed manifests are observations from the original campaign. Regenerating a campaign is a separate action; do it in a new scratch directory. The preparer below copies the frozen Python files and confirmation/fault protocol bytes unchanged and materializes upstream commit `dd96adba23b8100aba9147d71fffc63932fc6b0d` using the local Git object database. It does not need the Phase I inputs or run simulations.

## Prepare the environment and an empty replay tree

From the repository root, use CPython 3.14.7 for the closest replay of the recorded environment. The complete package versions are in `environment/requirements-lock.txt`; other Python versions or numerical builds have not been validated for byte-exact replay.

```sh
python3.14 -m venv .venv-phase2
.venv-phase2/bin/python -m pip install -r research/phase2-20261030/environment/requirements-lock.txt
.venv-phase2/bin/python -m pip check

.venv-phase2/bin/python research/phase2-20261030/reproduction/prepare_public_replay.py \
    --destination ../rnglab-phase2-replay
```

The destination must not exist. The preparer reads the pinned commit through `git archive`, so working-tree edits and a newer repository HEAD cannot change the upstream reference. A shallow clone may lack that commit; fetch it with `git fetch origin dd96adba23b8100aba9147d71fffc63932fc6b0d` and retry using a fresh destination. No source tarball is required.

The resulting layout includes:

```text
rnglab-phase2-replay/
  pytest.ini
  public-replay-preparation.json
  work/rnglab-main-snapshot/src/
  outputs/phase2-20261030/
    experiment/confirmation-v1/protocol.json
    experiment/fault-study-v1/protocol.json
    experiment/*.py
    analytic-dependence/*.py
    baseline/*.py
    rus/*.py
    theory/*.py
```

The manifest records the pinned source revision and checksums of every copied script and protocol. No published analysis outputs are copied into the empty result folders.

## Run the standalone 35 checks

Set the interpreter and scratch paths from the repository root:

```sh
PHASE2_PYTHON="$PWD/.venv-phase2/bin/python"
PHASE2_SCRATCH="$(cd .. && pwd)/rnglab-phase2-replay"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MPLBACKEND=Agg

(cd "$PHASE2_SCRATCH" && "$PHASE2_PYTHON" -m pytest -c pytest.ini -v \
    --junitxml=phase2-tests.xml)
```

The scratch configuration explicitly selects the standalone study (9), fault logic (14), baseline arithmetic checker (3), RUS construction (5) and analytic dependence checks (4). It avoids inheriting the repository's upstream pytest configuration. The baseline arithmetic tests operate on deterministic formulas and do not read the omitted Phase I files. The upstream comparison test resolves the pinned source through the prepared directory layout. The original consolidated run passed all 35 checks; running them again verifies the new environment.

## Regenerate the scientific campaigns

Use the same interpreter, scratch variables and thread limits defined above. These commands write only to the prepared tree. Generation drivers preserve existing main trace archives; reruns should use another new scratch tree.

```sh
"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/experiment/study.py" explore
"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/experiment/study.py" confirm
"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/experiment/fault_study.py" run

"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/analytic-dependence/analytic_dependence.py" --run
"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/analytic-dependence/analyze_and_plot.py"
"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/rus/rus_identifiability.py"
"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/theory/harmonic_reference_checks.py"

"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/experiment/plot_results.py"
"$PHASE2_PYTHON" "$PHASE2_SCRATCH/outputs/phase2-20261030/experiment/plot_fault_results.py"
```

Run confirmation before the fault study. Its required healthy conditions must pass their prespecified gates. The fault protocol also pins the healthy-decision file's checksum, so a numerical build that changes those bytes will stop at that check. Do not edit a frozen protocol to bypass a mismatch; inspect the environment and result differences, or create a separately versioned new study. The `freeze` commands are intentionally absent because the original confirmation/fault protocol bytes are already supplied.

The analytic dependence command constructs a target-preserving synthetic transition. It is distinct from the source faults in the fault experiment. The RUS command constructs two possible future laws conditional on one fixed observed prefix; it does not fit or validate a general forecast. The theory command performs deterministic reference checks.

## Published-payload and full-trace validation

The lightweight publication retains the main and fault analysis outputs, original seed manifests, protocol hashes and trace checksums. The supplied public payload validator rehashes the files that are actually distributed. From the repository root:

```sh
.venv-phase2/bin/python research/phase2-20261030/reproduction/verify_public_payload.py
```

The historical **132-check artifact audit** is a separate result. Its `experiment/validate_artifacts.py` checks the full set of ten Phase I inputs and the original large Phase II archives as well as independently recomputing primary intervals. It cannot complete from the lightweight public payload alone. Regenerating Phase II traces does not recreate the omitted Phase I trace or physical inventories. The retained audit JSON/log document what passed on the original full local dataset; their presence is not a claim that all 132 checks were rerun by a public reader.

A byte-for-byte file hash confirms integrity of a supplied artifact. Seed namespace separation confirms the recorded seed derivation; neither establishes independent physical entropy or a general guarantee of interval coverage. Interpret the regenerated results within the workload, protocols and limitations described in the research reports.
