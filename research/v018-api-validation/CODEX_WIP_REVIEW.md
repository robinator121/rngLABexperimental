# Review of the uncommitted Codex work `feat/pimc-run-options` (for Phase 2)

Reviewed 2026-11-05 by Claude Code (Opus 5.5) for the author, who decides whether to port it. Source: an uncommitted working tree in a local Codex output folder (`~/.codex/.chatgpt-projects/…/portfolio-audit-20261101/repos/rngLABexperimental`), dated 1 November 2026, based on commit `8501f9a`. `src/rnglab14/models/pimc_harmonic.py` is identical on `8501f9a` and the current `main`, so the patch applies cleanly. Nothing from it has been copied into this repository yet.

## What it changes

| File | Change | Relevant to Phase 2 |
|---|---|---|
| `src/rnglab14/models/pimc_harmonic.py` (+95 lines) | New keyword-only options on `run_simulation`: `initialization` (`"uniform"` default, or `"exact_gaussian"`), `initial_path`, `uncertainty_batch_sizes` (default `(32, 128)`), `record_trace` (default False). New result fields `initialization_protocol`, `source_provenance` (per-role byte counters, digests before and after, aliases) and `trace` (initial path, all retained paths, TD, CV, centroid, rg²) | yes |
| `tests/test_pimc_run_options.py` (181 lines, 23 tests) | Invalid options fail before any source is read; exact-Gaussian start matches the reference sampler and its covariance; an explicit path is copied and consumes nothing; the trace re-derives every estimator; recording a trace changes nothing; the role ledger counts bytes and aliases correctly | yes |
| `.github/workflows/ci.yml`, `pyproject.toml`, `scripts/smoke_installed.py`, `README.md`, `docs/CURRENT_REPRODUCTION.md`, `docs/REPRODUCIBILITY.md` | Packaging and CI changes and documentation edits | **no**; the README edit would conflict with the v0.18 rewrite |

## Test results on the current code (scratch copy of HEAD `59a448c` plus the patch)

- The patch applies with `git apply --check` and no conflicts.
- Its own 23 tests pass.
- Full suite: 373 passed, 2 failed (`test_physical15`). Both failures also occur **without** the patch in the same scratch copy, so they come from the copy's location (`rnglab15` rejects paths that are not regular files), not from the patch. In the real repository those tests pass.

## Against the milestone specification

| Milestone requirement (`research/phase2-20261030/reports/Research_advancement_and_next_milestone.md`) | Status in the Codex patch |
|---|---|
| Optional exact initial-path argument, shape and finiteness validated before sampling | **Met** (validated before any source read; tested) |
| …copied | **Met** (tested: no shared memory; caller's array unchanged) |
| …with its provenance recorded | **Partial.** Only `mode = "explicit_path"` is recorded. There is no digest of the supplied path and no way to say where it came from (gap G1) |
| Batch sizes exposed, or an audited raw-trace output, so BM512 runs through the package machinery | **Met, both ways** (`uncertainty_batch_sizes`, `record_trace`) |
| Full precision preserved in saved estimates and traces | **Met** through the trace and `summarize_series`; the legacy scalar fields keep their rounding by design |
| Default behaviour and default source consumption unchanged when the options are absent | **Plausible but not pinned.** No test compares every pre-existing result field, bit for bit, against the pre-patch code (gap G2). The new fields are additive, and `source_provenance` is now filled even in default calls |
| Regression tests for defaults, options, source counters and digests | **Mostly met**; defaults only partly (see G2) |

Other observations:
- The `"exact_gaussian"` mode draws two reference paths and discards one, because `sample_exact_finite_bead_paths` requires at least two. The waste is recorded in `initialization_protocol`.
- Invalid role objects now raise `TypeError` before any source is read. That changes behaviour only for invalid input.
- The style matches the surrounding module.

## Recommendation

**Port the two PIMC files only, credited to Codex, with two additions before Phase 2's protocol is frozen:**

- **G1 – path provenance:** record the SHA-256 of the supplied path's float64 bytes, plus an optional caller-supplied provenance mapping, in `initialization_protocol`.
- **G2 – golden default test:** before applying the patch, capture every legacy result field of a few fixed-seed default calls (with and without a centroid source) from the current code. After the patch, assert that they are bit-identical.

Do not port the CI, packaging or documentation changes.

The alternative is to reimplement from the specification. That would cost more for little gain: the patch meets the specification apart from G1 and G2, and its tests are thorough. Either way, the AI-disclosure section must record the authorship: Codex wrote the base patch and tests, Claude Code did the review and the G1/G2 additions.

## Port record (2026-11-05)

The author chose to port with both additions. It was done in this order, so each step can be checked on its own:

1. **G2 – golden default test, captured before the patch.** `tests/v018/test_pimc_default_golden.py` runs six fixed-seed default calls: a single source, the `rnglab15` applications call, the Phase 2 model with and without a separate centroid source, an aliased centroid source, and faulted mixed sources. It pins every pre-existing `PIMCResult` field, the final path, and each role source's byte counter and digest, bit for bit. The values in `tests/v018/data/pimc_default_golden.json` were captured from the committed, unpatched module; that file records the commit and the module's SHA-256. The capture script refuses to run if the module has uncommitted changes.
   - **Sensitivity check (manual, not saved as a script).** In two scratch copies of `src/`, the test failed on a one-ulp change to the proposal step and on a single extra byte read from the initialization source.
2. **Codex files applied unchanged.** `src/rnglab14/models/pimc_harmonic.py` and `tests/test_pimc_run_options.py` are byte-identical to the Codex working tree. The golden test still passes, which shows the defaults are unchanged.
3. **G1 – start-path provenance.** Every mode now records `initial_path_sha256`: the SHA-256 of the starting path as little-endian float64 bytes in bead order. A new keyword `initial_path_provenance`, allowed only together with `initial_path`, stores a JSON-serialisable mapping as its JSON round trip, recorded as given and not verified. Invalid provenance fails before any source is read. Tests: `tests/v018/test_pimc_initial_path_provenance.py`.

Full suite after the port: `docs/v018/evidence/pytest_phase2_api.txt`. The CI, packaging and documentation changes in the Codex tree were not ported.
