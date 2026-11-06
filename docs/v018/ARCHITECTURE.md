# v0.18 architecture: one namespace over the validated components

Written by Claude Code (Opus 5.5), 2026-11-05, Phase 0. The author decides which components count as validated; the table below states the evidence for each.

## Layout

```text
src/
├── rnglab/                 # v0.18 namespace (new modules) + legacy v0.1–v0.8 subpackages
│   ├── __init__.py         # legacy, unchanged ("Experimental framework V0.1", __version__ 0.1.0)
│   ├── sources.py          # v0.18: source hierarchy and channels        (re-export)
│   ├── faults.py           # v0.18: fault adapters                       (re-export; Phase 3 adds)
│   ├── models.py           # v0.18: PIMC, VMC, Ising, decision support   (re-export; Phase 2 extends PIMC)
│   ├── references.py       # v0.18: exact references                     (re-export)
│   ├── monitors.py         # v0.18: monitor calibration, health tests    (re-export; Phase 3 adds)
│   ├── stats.py            # v0.18: uncertainty and multiplicity         (re-export)
│   ├── longrun.py          # v0.18 Phase 1 (new): BM, OBM, SV, Geyer wrapper, batch-length diagnostic
│   ├── convergence.py      # v0.18 Phase 1 (new): rank-normalized split-R-hat
│   ├── susceptibility.py   # v0.18 Phase 1 (new): plug-in P, symmetry-aware S, ensemble jackknife G
│   ├── seeds.py            # v0.18 Phase 1 (new): namespace-derived seeds
│   └── v02/ … v08/, cli.py, degradation/, diagnostics/, experiments/, providers/, simulations/   # legacy, unchanged
├── rnglab10/ … rnglab17/   # historical snapshots, unchanged
├── rnglab_crypto/          # BLUE, unchanged
└── rnglab_finance/         # finance PoC, unchanged
```

`src/rnglab/` already held the legacy v0.1–v0.8 framework, so the v0.18 modules sit beside it as top-level modules. No legacy file was moved or edited, and no new name collides with a legacy subpackage. The v0.18 modules do not import the legacy subpackages.

**Open decision for the author.** `import rnglab` still reports `__version__ = "0.1.0"` and the legacy docstring, because changing `src/rnglab/__init__.py` is outside new v0.18 code. Options: leave it (current state), update the docstring and version, or move the legacy subpackages to another name (this needs your approval under rule 8).

## Phase 0 rule: re-export, do not copy

Every v0.18 name is the original object (`tests/v018/test_v018_namespace.py` checks identity), so the historical snapshots and their tests keep working, and a fix in v0.18 cannot silently diverge from the code that produced past results. New behaviour (Phases 1–3) will be added as new functions or new optional arguments. Phase 2 must change `rnglab14.models.pimc_harmonic.run_simulation`, and the milestone requires the defaults to stay bit-identical, which regression tests will pin.

## Mapping

| v0.18 name | Origin | Validation evidence | Status |
|---|---|---|---|
| `sources.RandomSource`, `SourceManifest`, `SourceBudgetExceeded` | `rnglab14.sources.base` | `tests/test_sources14.py`, `test_recovery_forecast_sources.py` (byte accounting, digest boundary) | tested |
| `sources.StochasticChannels` (the "channels") | `rnglab14.sources.base` | used by every model test; distinct objects do not prove independence | tested |
| `sources.PCG64Source`, `Xoshiro256PlusPlusSource`, `MersenneTwisterSource`, `ChaCha20Source`, `OSEntropySource` | `rnglab14.sources.reference` | `tests/test_sources14.py`; throughput measured in `docs/v018/evidence/source_throughput.json` | tested |
| `faults.DegradedFaultSource` | `rnglab14.sources.faults` | `tests/test_recovery_forecast_sources.py` | tested; not the Phase II channel faults |
| `models.PathIntegralHarmonicOscillator`, `PIMCResult`, `harmonic_energy_estimators`, `sample_exact_finite_bead_paths` | `rnglab14.models.pimc_harmonic` | `tests/test_pimc_centroid_heatbath.py` (19 tests), `tests/test_recovery_physics.py`; Phase II independent formula checks | formulas validated; interval coverage through this API **not yet validated** (Phase 2) |
| `models.harmonic_centroid_heatbath` | `rnglab14.models.pimc_harmonic` | `tests/test_pimc_centroid_heatbath.py`; `research/phase2-20261030/environment/centroid_update_patch.md` | tested |
| `models.HarmonicOscillatorVMC`, `VMCResult` | `rnglab14.models.vmc_harmonic` | `tests/test_models14.py` | tested |
| `models.Ising2D`, `IsingObservableResult` | `rnglab14.models.ising_2d` | `tests/test_models14.py` | tested |
| `models.ising_vectorized` | `rnglab16.core.ising` | `tests/trident/test_reference.py`; v0.16/v0.17 campaigns (outputs off-repo) | tested; plug-in susceptibility invalid at low T (Phase 1) |
| `models.wolff_vectorized` | `rnglab16.core.wolff` | v0.16/v0.17 exploratory cluster runs only (outputs off-repo); **no unit test** | untested |
| `models.MonteCarloDecisionSupport`, `DecisionOption`, `DecisionEvaluationResult` | `rnglab14.models.decision_support` | `tests/test_models14.py` | tested; Phase 4 pilot only |
| `references.exact_finite_bead_energy` | `rnglab14.models.pimc_harmonic` | equals the Phase II target 0.6554741522981742 to 1e-14 (namespace test) | validated |
| `references.ising_exact_enumeration` | `rnglab16.core.enumeration` | L4, T=1.5 susceptibility equals v0.17's 10.403226593514693 to 1e-12 (namespace test); `tests/trident/test_reference.py` | validated, L ≤ 4 only |
| `references.stationary`, `two_state`, `two_state_mean_variance`, `augmented_driver_kernel`, `ISING_TC` | `rnglab16.core` | `tests/trident/test_reference.py` | tested |
| `monitors.FrozenRankCalibration` | `rnglab14.monitoring.calibration` | `tests/test_recovery_physics.py` | tested |
| `monitors.rank_threshold` | `rnglab17.endpoints` | `tests/calibrate/test_endpoints.py` | tested |
| `monitors.RepetitionCountTest`, `AdaptiveProportionTest`, `HealthTestResult` | `rnglab14.standards.health_90b` | `tests/test_standards14.py` | tested |
| `stats.summarize_series` | `rnglab12.uncertainty` | `tests/test_statistics12.py` | tested; BM and Geyer only |
| `stats.batch_intervals`, `rate`, `causal_loss`, `window_metrics` | `rnglab17.endpoints` | `tests/calibrate/test_endpoints.py` | tested |
| `stats.summarize_independent_reference` | `rnglab14.models.pimc_harmonic` | `tests/test_recovery_physics.py` (through `run_exact_reference`) | tested |
| `stats.bonferroni_correction`, `holm_bonferroni_correction`, `benjamini_hochberg_fdr`, `MultiTestingResult` | `rnglab14.statistics.multitest` | `tests/test_forecast_monitoring14.py` | tested |

"Tested" means unit or regression tests exist and pass; it does not mean a held-out coverage or calibration study. "Validated" is reserved for agreement with an independent exact reference.

## Deliberately not in the v0.18 namespace

| Component | Where it stays | Reason |
|---|---|---|
| Physical emulators (GM, photodiode, serial) and conditioners (von Neumann, hash, DRBG, XOR) | `rnglab14.sources` | Synthetic research treatments; physical claims superseded (`docs/archive/V0141_CORRECTIONS.md`). |
| `DualIntegrityMonitor` | `rnglab14.monitoring` | Its "0.0% hazardous gap" result is superseded (`AUDIT.md` §2.2). |
| RUS predictors and trajectories | `rnglab13.forecast`, `rnglab14.forecast` | Out of v0.18 scope; Phase II RUS result is a non-identifiability counterexample. |
| SP 800-22 subset (3 tests), 90B entropy estimators, HMAC-DRBG | `rnglab14.standards` | Not used by v0.18 campaigns. |
| Phase II standalone sampler, channel faults, monitors | `research/phase2-20261030/experiment/*.py` | Frozen bytes of a completed study; must not become an import dependency. Phase 3 reimplements the needed adapters and monitors in `rnglab.faults` / `rnglab.monitors` and checks them against these definitions. |
| Crypto BLUE, finance | `rnglab_crypto`, `rnglab_finance` | Separate research lines. |

## Where v0.18 work will land

| Phase | Code | Tests | Research package |
|---|---|---|---|
| 1 estimator repair | new modules `rnglab/longrun.py`, `rnglab/convergence.py`, `rnglab/susceptibility.py`, `rnglab/seeds.py` | `tests/v018/test_longrun.py`, `test_convergence_susceptibility.py`, `test_v017_susceptibility_replay.py` | `research/v018-estimators/` |
| 2 API validation | optional arguments on `rnglab14.models.pimc_harmonic.PathIntegralHarmonicOscillator.run_simulation`, defaults unchanged | `tests/v018/` (defaults, options, counters, digests) | `research/v018-api-validation/` |
| 3 severity and monitors | `rnglab/faults.py`, `rnglab/monitors.py` | `tests/v018/` | `research/v018-severity/` |
| 4 decision pilot (exploratory) | script only | — | `research/v018-decision-pilot/` |

## Test counts

| Scope | Count | Source |
|---|---:|---|
| Baseline suite on `main` 697c943 | 322 passed | `docs/v018/evidence/pytest_baseline.txt` |
| New v0.18 namespace tests | 9 | `tests/v018/test_v018_namespace.py` |
| Full suite after Phase 0 | 331 passed | `docs/v018/evidence/pytest_phase0.txt` |
| Phase II research checks (separate, not re-run) | 35 | `research/phase2-20261030/environment/pytest-phase2-consolidated.xml` |
