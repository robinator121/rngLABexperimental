# AI contributions log (v0.18)

One line per contribution, newest last. Claude Code (Anthropic) did everything not listed here: design, protocols, campaign code, analysis, reports, and review of the work below. The author made every decision.

| Date | AI | What it did or said |
|---|---|---|
| 2026-11-01 | OpenAI Codex | Wrote the PIMC run options and their 23 tests; ported unchanged in Phase 2 |
| 2026-11-05 | OpenAI Codex | Found 3 Phase 1 analysis bugs (R̂ denominator, NaN for constant chains, verdicts on partial gate sets); fixed in Amendment 01 |
| 2026-11-05 | ChatGPT (PR #2 review) | Flagged the multiplicity problem in the Phase 1 control gates; documented in Amendment 01 |
| 2026-11-07 | Codex (gpt-6-astra) | CLI test prompt: answered correctly |
| 2026-11-07 | Antigravity (Gemini) | CLI test prompt: answered correctly |
| 2026-11-07 | Local qwen2.5:3b | Asked for H2 degradation models; gave generic reliability metrics instead (mean time to failure, bit error rate, state repetition). Weak; one idea kept: a shrinking generator cycle as a degradation mode |
| 2026-11-07 | Codex (gpt-6-astra, high) | Wrote the Stage C `report_tables.py` and `plot.py` (255 lines) against sample outputs. Claude reviewed: format-only, no recomputed statistics, runs cleanly; accepted unchanged |
| 2026-11-07 | Antigravity (Gemini) | Code review of Stage C analysis: first two runs produced nothing (headless mode auto-denied file reads) |
| 2026-11-07 | Antigravity (Gemini) | Literature search from memory, 15 references; Crossref check (`scripts/v018/verify_dois.py` → `docs/v018/evidence/agy_literature_check.json`): 7 fully match, 6 DOIs point to other papers, 2 DOIs do not exist. New verified: Fiala et al. 2012 and Snir et al. 2014 (silent data corruption in HPC), Vattulainen et al. 1995, James 1990, Bauke & Mertens 2007 |
| 2026-11-07 | Antigravity (Gemini) | Code review, third try with code pasted inline: produced nothing (tried a shell command; headless mode denied it). Review reassigned to Codex |
| 2026-11-07 | Codex (gpt-6-astra, high, read-only) | Reviewed the Stage C analysis against the protocol: 7 valid findings (implementation failures not blocking verdicts; incomplete replay passing; Bonferroni divisor not frozen at 24; recount ignoring zero-variance intervals; possible divide-by-zero; recount not checking rates; missing 95% paired intervals). None changes this campaign; fixes proposed as pre-analysis Amendment 01 |
