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

## Second rewrite (2026-11-08)

On 2026-11-08, at the author's request, the commit messages of `feature/v018` were rewritten again to remove AI-attribution text. This was done with `git filter-branch --msg-filter`, and the branch was force-pushed. The author will add the AI disclosure once the project is finished.

**What changed and what did not:**
- **Changed:** the hashes of 22 of the 28 branch commits, and the messages of 20 of them. Commits 1 to 6 had no such text, so they keep their hashes.
- **Unchanged,** checked for every commit:
  - the file trees, so every file's content and SHA-256 is unchanged;
  - the author, committer and dates.
- **Unchanged by design:** the in-repository statements (READMEs, reports, protocols, `docs/v018/AI_LOG.md`), and file names such as `CODEX_WIP_REVIEW.md`, which commit messages still cite.
- **Earlier citations:** the old hashes in the table above that point to the first rewrite now resolve through this table as well.

The old tip is `19f8b045fcca11c1fbf08f967f18b35c37bfb653`. A local backup branch (`backup/feature-v018-before-disclosure-rewrite-20261108`) keeps the old commits on the author's machine; it was not pushed.

**Scripts updated in the same commit as this section:** `amendment01_compare.py` and `amendment01_record.py` in `research/v018-h1-detection/` read the superseded scripts with `git show <commit>:…`. Their commit constant now names the rewritten manifest commit, whose tree is identical. The frozen records they produced are not edited.

| # | Old commit | New commit | Subject | Files citing the old hash |
|---:|---|---|---|---|
| 7 | `802a8c826f0e87220b2bbe010e0725b0253adcc0` | `99330afbe949b1fa6dca628b2ec55817eb4e4b84` | v0.18 Phase 1: Amendment 01 (pre-analysis) and adopted roadmap | — |
| 8 | `6f49513eac76450f0168a2a7f3cbb39b12c45745` | `a20c6bd9ffa09b468660af12705d3cb35d7dbbdf` | v0.18 Phase 1: results, independent recount, figures and report | — |
| 9 | `592acb0230232034264552ba11a27439680f5d95` | `0f9c0da8f6eb11d0f674ec51abc8f2e2a94acff7` | v0.18 Phase 2: golden pin of PIMC default behaviour | — |
| 10 | `8b531c8c41b992e522017b6899d0cb7e7484d003` | `b7753b3453c183ec363ada14cea163df32412130` | v0.18 Phase 2: PIMC run options from WIP (unchanged) | — |
| 11 | `5322e0863f07a1510a63981c71e7e5de0801333e` | `4ee2c4bef47e095f1f0560d9b2c509150c3b7bff` | v0.18 Phase 2: start-path digest and provenance (G1) | — |
| 12 | `5ef2e4284de6e71e1aee59edc1d81e5316e1dd8a` | `c5084d81c4ac8bd647befe51e26a1e8dc388ed33` | tests: use a test-only provenance label | — |
| 13 | `bbb94dc88e7e9d2682179995b346961afc6030ad` | `1a9cea7ec0732952e531ec67f6f6399328e3a450` | v0.18 Phase 2: freeze API-validation protocol | — |
| 14 | `9cc4b1966ca2243232030b689794b35a08766eba` | `fd7c1eebf3fd232bc909e2a5607fceb6f388f326` | v0.18 Phase 2: campaign, analysis and recount scripts | — |
| 15 | `8fa78114a774de4057cc7e16afaf3dd4ebc5f116` | `0375aa6803ee190b257e0799f92549cafb5e569c` | v0.18 Phase 2: pre-generation manifest | — |
| 16 | `1a2eb9f9317f1e604b45fd07ffe1bc5dea3b8705` | `6217889f1409dd1276750ef8d018cf506e963335` | v0.18 Phase 2: results, independent recount, figures and report | — |
| 17 | `9bd267c4a98e3681c948bbfd10d6fa4172e3c725` | `88e1b3b48b125937a4d438e8dbad6e565468100a` | v0.18 Stage B: novelty review (widened Phase 3.5), adopted as GO | — |
| 18 | `015bed62d6934c9196daf25bc0d30b21e6cb1b77` | `33d63f7480b3762ccd1b82f6f6a442b1ea909632` | v0.18 Stage C: freeze H1 detection protocol | `research/v018-h1-detection/README.md`, `research/v018-h1-detection/REPORT.md` |
| 19 | `ec26a161070f1dfea603fef3c1fd8a489a7d4c93` | `e040528efe797cee3b48cb4dba952d8516801a0b` | docs: commit hash map for the feature/v018 trailer rewrite | — |
| 20 | `676d102ca7d456badd67289935a6f44a9528e0a4` | `e74ce8e4f986b7a08b410116853fb1c24d281427` | v0.18 Stage C: fault adapters, conventional battery, monitors, pilot | — |
| 21 | `4d1410c0dd44a0ccef0caf98efa16d9503fed8f7` | `7b2e8f4abd95db43c7bb78df748d1ab167121e84` | v0.18 Stage C: pilot summary and campaign grids | `research/v018-h1-detection/REPORT.md` |
| 22 | `95f5f1e655cdf40eed3441e62afcd672715ae31b` | `796f913b46573e8c9e8010a26b78d2c86a6c3861` | v0.18 Stage C: campaign, analysis and recount scripts | `research/v018-h1-detection/REPORT.md`, `research/v018-h1-detection/manifest_pregeneration.json`, `research/v018-h1-detection/results/run_manifest.json` |
| 23 | `0aca45ae88ab76f369a0118578f37139c0adca94` | `e5c014022dc7a11be6dbc2b8a73b36890b62f101` | v0.18 Stage C: pre-generation manifest | `docs/v018/evidence/pytest_stageC_amendment01.txt`, `research/v018-h1-detection/AMENDMENT_01.md`, `research/v018-h1-detection/README.md`, `research/v018-h1-detection/REPORT.md`, `research/v018-h1-detection/REPRODUCE.md`, `research/v018-h1-detection/results/amendment01_preanalysis.json`, `research/v018-h1-detection/results/report_tables.md`, `research/v018-h1-detection/results/run_manifest.json` |
| 24 | `7a188d9887d293eb8648967f69aa7fa8f43a3438` | `5079bdbbce61f4b3bb65efb550afd1f6fde64ba6` | v0.18 Stage C: campaign data | `research/v018-h1-detection/REPORT.md` |
| 25 | `456ebcd83a147ee87c29274688c0fc5854e1d8a1` | `3f495685df8b5d904e54fca858910fdbf006b4c8` | v0.18 Stage C: Amendment 01, approved before analysis | `research/v018-h1-detection/REPORT.md`, `research/v018-h1-detection/REPRODUCE.md` |
| 26 | `11d3d94b5f8c41fd4844e91c498cf90756c2554a` | `5e7088f8cbd7470560286e09162dabd99643e0be` | v0.18 Stage C: report scripts, contribution log and citation check | — |
| 27 | `2656da5d2cbdc99333d351fcbdf2d0a63db708af` | `edafbdaa50bb67734e8ac5f1d6f2b2f7fe1ac7a3` | v0.18 Stage C: analysis results | — |
| 28 | `19f8b045fcca11c1fbf08f967f18b35c37bfb653` | `bed49942136b0a8c312a2d858230ce0963f92f24` | v0.18 Stage C: report | — |
