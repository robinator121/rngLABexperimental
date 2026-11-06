# Amendment 01 to protocol rnglab-v018-estimators-v1 (pre-analysis)

**Status: approved by the author on 2026-11-05 ("Approve and analyze") and committed before `analyze.py` was first run on campaign outputs.** Drafted by Claude Code (Opus 5.5). Generation had finished (19 of 19 settings, 10:39 UTC) before approval; no output had been analyzed. `protocol.json` (SHA-256 `ab2f16a6…0ec66`), its gates, endpoints, sample sizes and namespaces are **unchanged**. No generated file is modified.

## Why

The automated Codex review on PR #2 reported three implementation defects, and the author's ChatGPT review raised two interpretation questions:

1. `rnglab.convergence.split_rhat` used normal scores Φ⁻¹((r − 3/8)/(S − 1/4)). The cited definition (Vehtari et al. 2021, Blom scores) is Φ⁻¹((r − 3/8)/(S + 1/4)). The wrong denominator makes complementary ranks asymmetric, so R̂(x) ≠ R̂(−x); reproduced: 1.02347 vs 1.02420 on a 4×20 example.
2. Zero within-chain variance returned NaN even when the chains sat at different constant values, and the analysis counted NaN as "not flagged", so total non-mixing could be read as converged.
3. `analyze.py` gave verdicts on whatever gates existed, so a partial set of passing cells could yield "supported", contrary to the protocol's "all nine cells".
4. The control gate L-C requires 20 separate 95% intervals (10 AR(1) settings × 2 criteria) to all contain their targets. Under a perfect generator and independence, at least one misses with probability 1 − 0.95²⁰ ≈ 0.64. The frozen wording "a failure means a generator defect" is therefore too strong. This is a design flaw in the frozen protocol.
5. Estimators S and G both assume that m² is stationary after the 512-sweep burn-in. Symmetric initial spins make the sign of m exactly symmetric at every sweep, but they do not make the |m| distribution stationary. The protocol had no diagnostic that separates an estimator failure from a burn-in failure.

## What had been seen when this was written

Only the runner's progress log (`results/run.log`: file names and wall seconds) and process status. No endpoint, gate value or per-unit number from the campaign outputs had been read. The pipeline checks used for this amendment ran on smoke-namespace data in a scratch folder.

## Changes

| ID | Change | Files | Effect on endpoints |
|---|---|---|---|
| A1 | Blom denominator S + 1/4; +∞ for separated constant chains; NaN only when all draws are tied | `src/rnglab/convergence.py`, `tests/v018/test_convergence_susceptibility.py` | Implementation corrected to the cited definition. The generation run stored R̂ computed with the old code in `results/ising/*.npz` (`rhat_m`, `rhat_m2`, `rhat_e`); those files stay **byte-identical**. `amend01_from_traces.py` recomputes R̂ from the saved integer traces (after checking each trace's SHA-256 against the generation record) into `results/amend01/*.npz`. Gate D-F1 and the R̂ secondary tables use the recomputed values. `results/gates.json` reports, per cell and series, the maximum difference and the old and new flagged counts. |
| A2 | Verdicts require every registered cell, setting, unit count and gate ID; otherwise `NOT_EVALUABLE_INCOMPLETE`, and no gates are evaluated on incomplete inputs | `research/v018-estimators/analyze.py` | Enforces the protocol's verdict text; no gate changes |
| A3 | The frozen L-C gate is evaluated and reported unchanged. Added secondary view: Bonferroni over the 20 control intervals (per-interval level 0.0025), using the exact χ² interval for the variance ratio and a t interval for the mean (the AR(1) replica mean is exactly Gaussian) | `analyze.py` | Interpretation only. **Pre-specified reading:** an L-C FAIL is a *flag for investigation*, not a confirmed generator defect. If the adjusted view also flags a setting, that setting is regenerated under the new namespace `rnglab-v018-estimators-amend01` as a separately labelled check, and both results are reported. The frozen "controls_valid" verdict is reported as computed, whatever happens. |
| A4 | Split-half burn-in diagnostic per Ising cell: second-half minus first-half chain means of m² (relative to χ) and of energy (relative to the exact energy), 95% t-intervals across chains | `amend01_from_traces.py`, `analyze.py` | Secondary, descriptive. **Pre-specified reading:** if an S-S or S-G bias criterion fails in a cell whose split-half interval excludes 0, the failure is reported as "estimator or burn-in, not separable"; without that drift signal, as "estimator". |

## Not changed

`protocol.json`, `protocol.sha256`, every gate criterion and unit count, the seed namespace, `run_campaign.py`, `src/rnglab/longrun.py`, `src/rnglab/susceptibility.py`, and every file written by the generation run.

## Provenance and reproduction

- Generation used the code pinned in `manifest_pregeneration.json` (commit `f3be97a`). `run_campaign.py` checks those hashes, so a byte-exact replay of generation must be run at that commit: after this amendment, `analyze.py` and `convergence.py` no longer match the manifest.
- Analysis uses the amended files. Their SHA-256 values are in `amendment_01.json`.
- Order: generation (at `f3be97a`) → `amend01_from_traces.py` → `analyze.py` → `recount.py` → `plot.py`.

## Option for the author (not part of this amendment)

D-F1 rests on R̂, whose implementation changed after the freeze. If that is considered compromised, a separately pre-registered replication of the L4 T1.5 cell under a fresh namespace takes about 10 minutes of CPU.
