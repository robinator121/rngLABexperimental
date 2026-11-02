# RNGLab BLUE evidence package

Start with the [scientific findings](../../docs/crypto/RESULTS.md) and [reproduction instructions](../../docs/crypto/REPRODUCE.md).

- `confirmation-reviewed/`: final 128-replica paired campaign, immutable protocol/results/report/plots/hash manifest.
- `sensitivity-reviewed/`: planned 64-replica longer-horizon sensitivity study.
- `configs/`: unchanged registered designs for those two studies.
- `evidence/`: original baseline, crypto and final tests, independent review, replay/build checks, dependency pins and separate overhead observations.
- `os-reference/`: one captured OS CSPRNG statistical sample, deliberately public and never used as a real key.
- `confirmation-v1/`, `sensitivity-v1/`: preserved initial checkpoint runs; reviewed versions correct attribution and execute the source-gated branch explicitly. These are the **same experimental units**, not extra samples.
- `smoke/`: development evidence only, excluded from scientific inference.

All adversarial demonstrations use locally generated synthetic material. No unvalidated physical samples are used to provision keys. Source and result hashes, raw score trajectories and explicit censoring support review. Current statistical limitations—including failed confirmation criteria—are retained rather than selected away.
