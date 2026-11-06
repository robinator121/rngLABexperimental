# Commit hash map for the `feature/v018` history rewrite (2026-11-06)

On 2026-11-06, at the author's request, the 18 commits of `feature/v018` that are not on `main` were rewritten to remove the `Co-Authored-By: Claude …` trailer, which 11 of them carried. This was done with `git filter-branch --msg-filter`, and the branch was force-pushed to GitHub.

**What changed and what did not:**
- **Changed:** the commit hashes and the commit messages, which lost only the trailer line.
- **Unchanged,** checked for every commit:
  - the file trees, so every file's content and SHA-256 is unchanged;
  - the author, committer and dates;
  - the rest of each message. The plain-text AI disclosure lines in 8 commit messages are kept.
- **Unchanged by design:** the in-repository AI-disclosure statements (README, reports, protocols).

**Why this map exists.** Several committed files cite the old hashes. Some are frozen protocols, manifests and past results, which are never edited (research rule 3); they are listed below. Use this table to resolve an old hash to the commit that now carries the same tree.

Generated with:

```bash
paste <(git rev-list --reverse origin/main..<old tip>) <(git rev-list --reverse origin/main..feature/v018)
```

The old tip is `5a2a83b5f867f820c88698c925ad02aaafef388e`. A local backup branch (`backup/feature-v018-before-trailer-rewrite-20261106`) keeps the old commits on the author's machine; it was not pushed.

| # | Old commit | New commit | Subject | Files citing the old hash |
|---:|---|---|---|---|
| 1 | `7548438b752b78c63febee06615182b89ca5c073` | `96f773b0524b6cfbd5789962a9fc5fdb9fcd7271` | v0.18 Phase 0: claim audit, consolidated namespace, README rewrite | — |
| 2 | `44676eb4034c71f418af5d07d71712f1297439fc` | `ec280a359df012969053743d230883037a668b59` | v0.18 Phase 1: freeze estimator-repair protocol | `research/v018-estimators/README.md`, `research/v018-estimators/REPORT.md` |
| 3 | `7cfa3d7f1c0b110ac2a1c8105622a7d3b5a4f206` | `4781a86e83add96f69bd01f185cb6d437001e2d9` | v0.18 Phase 1: estimators, diagnostic and campaign code (pre-generation) | `research/v018-estimators/manifest_pregeneration.json` |
| 4 | `f3be97a3698f2ddd9e6316e956411938ba672810` | `056d6410495b97283868bfa55a969d801258d5be` | v0.18 Phase 1: pre-generation manifest | `research/v018-estimators/AMENDMENT_01.md`, `research/v018-estimators/README.md`, `research/v018-estimators/REPORT.md`, `research/v018-estimators/REPRODUCE.md`, `research/v018-estimators/amendment_01.json` |
| 5 | `cd6eb793e9521d65858d3331a8dd4ce40599beda` | `b54fbf6d9d5915a6b29a8f3428e57b33b5aea9d4` | v0.18 Phase 1: recount, figures and reproduction guide | `docs/v018/evidence/ci_steps_local.txt` |
| 6 | `657b666742c89e8b1498dcdb6925840a372e6514` | `c15e398ff866be702a517c7f933ddf6a24091b1a` | docs: CI has never run; record the local CI-step run | — |
| 7 | `e9f102260a52809940ec30ad115394faa1ca4ef6` | `802a8c826f0e87220b2bbe010e0725b0253adcc0` | v0.18 Phase 1: Amendment 01 (pre-analysis) and adopted roadmap | `research/v018-estimators/REPORT.md` |
| 8 | `59a448c70710166b1e3af82c4c9da527b0f0d4fe` | `6f49513eac76450f0168a2a7f3cbb39b12c45745` | v0.18 Phase 1: results, independent recount, figures and report | `docs/v018/evidence/pytest_phase2_api.txt`, `research/v018-api-validation/CODEX_WIP_REVIEW.md`, `tests/v018/data/pimc_default_golden.json` |
| 9 | `caaf784b21d23ec80938560417bf6cdc79d5b45c` | `592acb0230232034264552ba11a27439680f5d95` | v0.18 Phase 2: golden pin of PIMC default behaviour | `research/v018-api-validation/REPORT.md` |
| 10 | `c47e44f178015dd95448a021acf89a97ff7152f0` | `8b531c8c41b992e522017b6899d0cb7e7484d003` | v0.18 Phase 2: PIMC run options from Codex WIP (unchanged) | `research/v018-api-validation/REPORT.md` |
| 11 | `4cf55a879325ecac2b848c4bd0f159d57fb3957d` | `5322e0863f07a1510a63981c71e7e5de0801333e` | v0.18 Phase 2: start-path digest and provenance (G1) | `research/v018-api-validation/PROTOCOL.md`, `research/v018-api-validation/REPORT.md`, `research/v018-api-validation/protocol.json` |
| 12 | `7585b1fca71b32ae72882a7f8f5def0a46e45039` | `5ef2e4284de6e71e1aee59edc1d81e5316e1dd8a` | tests: use a test-only provenance label | — |
| 13 | `0cbefb466e743640e88fc44a921b2cacc8fbbfb5` | `bbb94dc88e7e9d2682179995b346961afc6030ad` | v0.18 Phase 2: freeze API-validation protocol | `research/v018-api-validation/README.md`, `research/v018-api-validation/REPORT.md` |
| 14 | `f9fa0f220178f40e9280eb999eb41d8451f80360` | `9cc4b1966ca2243232030b689794b35a08766eba` | v0.18 Phase 2: campaign, analysis and recount scripts | `research/v018-api-validation/REPORT.md`, `research/v018-api-validation/manifest_pregeneration.json`, `research/v018-api-validation/results/run_manifest.json` |
| 15 | `ff948e2e63b7eebab0e39177e6568f114b95d458` | `8fa78114a774de4057cc7e16afaf3dd4ebc5f116` | v0.18 Phase 2: pre-generation manifest | `docs/v018/evidence/pytest_phase2.txt`, `research/v018-api-validation/README.md`, `research/v018-api-validation/REPORT.md`, `research/v018-api-validation/REPRODUCE.md`, `research/v018-api-validation/results/report_tables.md`, `research/v018-api-validation/results/run_manifest.json` |
| 16 | `a2800497759c6e75a7b85bedd2e2292c01522135` | `1a2eb9f9317f1e604b45fd07ffe1bc5dea3b8705` | v0.18 Phase 2: results, independent recount, figures and report | — |
| 17 | `8523d32b715b3e78a653fff3bef6a69ae7231bcb` | `9bd267c4a98e3681c948bbfd10d6fa4172e3c725` | v0.18 Stage B: novelty review (widened Phase 3.5), adopted as GO | — |
| 18 | `5a2a83b5f867f820c88698c925ad02aaafef388e` | `015bed62d6934c9196daf25bc0d30b21e6cb1b77` | v0.18 Stage C: freeze H1 detection protocol | — |

Produced with Claude Code (Anthropic) under the author's direction.
