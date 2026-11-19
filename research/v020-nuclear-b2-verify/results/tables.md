# B2 verification: tables from `results/verification.json`

Namespace `rnglab-v020-nuclear-b2-verify`; code `e0d3b09af706d3a13a33479e29045f27268154d8` (committed tree; files of earlier attempts: 0; ignored files allowed (bytecode caches, smoke/, egg-info): 2; other uncommitted or ignored files: none); started 2026-11-19T17:43:04.829519+00:00; 39.4 s in total.

## Verdicts

- R1: **PASS**
- R2: **MET**
- R0: NC1: C1 holds: the false-FAIL rate is below 0.02 (95% upper bound 0.0097)
- R0: NC2: C1 holds: the false-FAIL rate is below 0.02 (95% upper bound 0.0088)
- R0: R1 (a FAIL in NC1 or NC2): a correct engine FAILs in 835 of 50000 null batteries (rate 0.0167, 95% interval 0.0156 to 0.0179); reported beside C1, which is per family

| Deterministic check | Holds |
|---|---|
| D8 interval laws | True |
| D7 ensemble streams | True |
| D9 big populations | True |
| D3 xref | True |
| D4 replay | True |
| D5 checkpoint | True |
| D6 sink | True |
| D1 consumption | True |
| D2 bookkeeping and Layer A invariants | True |
| D10 planned tests ran | True |

D8: 113 interval laws; largest |decimal - path-tree expm| 1.57e-13; largest |end states - Layer A expm| 1.57e-13; largest |sum - 1| 5.91e-47; brackets hold: True; most digits used: 50; non-finite references: none.

| XREF case | Parts | Digests equal | Arrays equal | Streams equal | Reference engine seconds |
|---|---:|---|---|---|---:|
| S2, part_max 1125899906842624 | 1 | True | True | True | 0.06 |
| S4, part_max 997 | 101 | True | True | True | 0.02 |

## Families (Holm, family-wise alpha 0.01)

| Family | Tests | Smallest p | Smallest Holm-adjusted p | Rejections |
|---|---:|---:|---:|---|
| NC1 | 222 | 0.00379 | 0.84 | none |
| NC2 | 85 | 0.0303 | 1 | none |

### NC1: every test

| Test | Statistic (z for binomial and moment tests) | p | Holm-adjusted p |
|---|---:|---:|---:|
| NC1 S1 BT-STATE Sa-400 t=0.5 | 0.153 | 0.696 | 1 |
| NC1 S1 BT-STATE Sa-400 t=1.0 | 0.762 | 0.383 | 1 |
| NC1 S1 BT-STATE Sa-400 t=2.0 | 0.833 | 0.361 | 1 |
| NC1 S1 BT-STATE Sa-400 t=4.0 | 0.535 | 0.464 | 1 |
| NC1 S1 BT-STATE Sa-400 t=6.0 | 2.15 | 0.143 | 1 |
| NC1 S1 BT-DECAY Sa-400 | -1.47 | 0.143 | 1 |
| NC1 S1 BT-BINS Sa-400 | 5.8 | 0.326 | 1 |
| NC1 S1 BT-PHOTON Sa-400/beta-->Sb-400 | 0.013 | 0.992 | 1 |
| NC1 S1 BT-PHOTON-T Sa-400/beta-->Sb-400 | 3.23 | 0.664 | 1 |
| NC1 S2 BT-STATE Xa-400 t=0.5 | 4.6 | 0.596 | 1 |
| NC1 S2 BT-STATE Xa-400 t=1.0 | 11.8 | 0.107 | 1 |
| NC1 S2 BT-STATE Xa-400 t=2.0 | 9.82 | 0.278 | 1 |
| NC1 S2 BT-STATE Xa-400 t=4.0 | 6.07 | 0.64 | 1 |
| NC1 S2 BT-STATE Xa-400 t=8.0 | 4.87 | 0.771 | 1 |
| NC1 S2 BT-STATE Xa-400 t=12.0 | 5.88 | 0.66 | 1 |
| NC1 S2 BT-DECAY Xa-400 | -0.497 | 0.619 | 1 |
| NC1 S2 BT-DECAY Xb-396 | 0.917 | 0.36 | 1 |
| NC1 S2 BT-DECAY Xc-400m | -0.936 | 0.35 | 1 |
| NC1 S2 BT-DECAY Xc-400 | -0.812 | 0.417 | 1 |
| NC1 S2 BT-DECAY Xd-396 | 1.34 | 0.182 | 1 |
| NC1 S2 BT-DECAY Xe-396 | -0.571 | 0.57 | 1 |
| NC1 S2 BT-BINS Xa-400 | 3.7 | 0.717 | 1 |
| NC1 S2 BT-CHANNEL Xa-400 | 1.05 | 0.294 | 1 |
| NC1 S2 BT-PHOTON Xa-400/alpha->Xb-396 | -0.648 | 0.52 | 1 |
| NC1 S2 BT-PHOTON Xa-400/beta-->Xc-400m | 1.03 | 0.304 | 1 |
| NC1 S2 BT-PHOTON Xb-396/ec->Xd-396 | 0.543 | 0.59 | 1 |
| NC1 S2 BT-PHOTON Xc-400m/it->Xc-400 | -1.33 | 0.184 | 1 |
| NC1 S2 BT-CHANNEL Xe-396 | 0.64 | 0.531 | 1 |
| NC1 S2 BT-PHOTON Xe-396/beta-->Xg-396 | 0.888 | 0.378 | 1 |
| NC1 S2 BT-PHOTON-T Xa-400/alpha->Xb-396 | 8.33 | 0.215 | 1 |
| NC1 S2 BT-PHOTON-T Xa-400/beta-->Xc-400m | 3 | 0.809 | 1 |
| NC1 S2 BT-PHOTON-T Xb-396/ec->Xd-396 | 6.88 | 0.332 | 1 |
| NC1 S2 BT-PHOTON-T Xc-400m/it->Xc-400 | 8.61 | 0.197 | 1 |
| NC1 S2 BT-PHOTON-T Xe-396/beta-->Xg-396 | 5.02 | 0.414 | 1 |
| NC1 S2 BT-PATH Xa-400 | 37.5 | 0.707 | 1 |
| NC1 S2 BT-PATH Xb-396 | 7.93 | 0.541 | 1 |
| NC1 S2 BT-PATH Xc-400m | 11.1 | 0.852 | 1 |
| NC1 S2 BT-PATH Xc-400 | 16.4 | 0.29 | 1 |
| NC1 S2 BT-PATH Xd-396 | 2.96 | 0.705 | 1 |
| NC1 S2 BT-PATH Xe-396 | 4.96 | 0.665 | 1 |
| NC1 S3 BT-STATE Ya-410 t=1.0 | 1.25 | 0.534 | 1 |
| NC1 S3 BT-STATE Ya-410 t=5.0 | 0.649 | 0.723 | 1 |
| NC1 S3 BT-STATE Ya-410 t=10.0 | 2.19 | 0.334 | 1 |
| NC1 S3 BT-STATE Ya-410 t=20.0 | 6.09 | 0.0476 | 1 |
| NC1 S3 BT-STATE Ya-410 t=40.0 | 1.57 | 0.457 | 1 |
| NC1 S3 BT-DECAY Ya-410 | -0.523 | 0.601 | 1 |
| NC1 S3 BT-DECAY Yb-410 | -0.865 | 0.387 | 1 |
| NC1 S3 BT-BINS Ya-410 | 7.08 | 0.215 | 1 |
| NC1 S4 BT-STATE G1 t=0.5 | 0.0652 | 0.968 | 1 |
| NC1 S4 BT-STATE G1 t=1.0 | 0.409 | 0.815 | 1 |
| NC1 S4 BT-STATE G1 t=2.0 | 1.57 | 0.455 | 1 |
| NC1 S4 BT-STATE G1 t=5.0 | 2.53 | 0.282 | 1 |
| NC1 S4 BT-STATE G1 t=13.0 | 2.33 | 0.127 | 1 |
| NC1 S4 BT-DECAY G1 | 0.11 | 1 | 1 |
| NC1 S4 BT-BINS G1 | 1.2 | 0.879 | 1 |
| NC1 S4 BT-NEUTRON G1 | 8.67 | 0.123 | 1 |
| NC1 S4 BT-STATE G2 t=0.5 | 0.674 | 0.714 | 1 |
| NC1 S4 BT-STATE G2 t=1.0 | 1.69 | 0.429 | 1 |
| NC1 S4 BT-STATE G2 t=2.0 | 2.85 | 0.241 | 1 |
| NC1 S4 BT-STATE G2 t=5.0 | 0.86 | 0.65 | 1 |
| NC1 S4 BT-STATE G2 t=13.0 | 1.11 | 0.573 | 1 |
| NC1 S4 BT-DECAY G2 | -0.302 | 0.759 | 1 |
| NC1 S4 BT-BINS G2 | 2.92 | 0.712 | 1 |
| NC1 S4 BT-NEUTRON G2 | 5.44 | 0.365 | 1 |
| NC1 S4 BT-STATE G3 t=0.5 | 0.125 | 0.724 | 1 |
| NC1 S4 BT-STATE G3 t=1.0 | 0.00163 | 0.968 | 1 |
| NC1 S4 BT-STATE G3 t=2.0 | 0.279 | 0.598 | 1 |
| NC1 S4 BT-STATE G3 t=5.0 | 1.86e-04 | 0.989 | 1 |
| NC1 S4 BT-STATE G3 t=13.0 | 0.225 | 0.635 | 1 |
| NC1 S4 BT-DECAY G3 | -0.474 | 0.639 | 1 |
| NC1 S4 BT-BINS G3 | 1.37 | 0.928 | 1 |
| NC1 S4 BT-NEUTRON G3 | 1.37 | 0.928 | 1 |
| NC1 S4 BT-CHANNEL G1 | -1.53 | 0.128 | 1 |
| NC1 S4 BT-CHANNEL G2 | 1.01 | 0.311 | 1 |
| NC1 S5 BT-STATE Za-450 t=1.0 | 0.00727 | 0.932 | 1 |
| NC1 S5 BT-STATE Za-450 t=1000.0 | 0.541 | 0.969 | 1 |
| NC1 S5 BT-STATE Za-450 t=3000.0 | 2.18 | 0.703 | 1 |
| NC1 S5 BT-DECAY Za-450 | -0.574 | 0.566 | 1 |
| NC1 S5 BT-DECAY Zb-450 | -0.574 | 0.566 | 1 |
| NC1 S5 BT-DECAY Zc-446 | -0.464 | 0.641 | 1 |
| NC1 S5 BT-DECAY Zd-446 | -0.56 | 0.576 | 1 |
| NC1 S5 BT-BINS Za-450 | 0.351 | 0.95 | 1 |
| NC1 S5 BT-PHOTON Zb-450/alpha->Zc-446 | 1.83 | 0.0673 | 1 |
| NC1 S5 BT-CHANNEL Zd-446 | -0.0806 | 0.935 | 1 |
| NC1 S5 BT-PHOTON-T Zb-450/alpha->Zc-446 | 3.39 | 0.335 | 1 |
| ENS ET-MEAN Xa-400 t=2.0 | -0.188 | 0.852 | 1 |
| ENS ET-VAR Xa-400 t=2.0 | 417 | 0.515 | 1 |
| ENS ET-MEAN Xb-396 t=2.0 | 0.147 | 0.882 | 1 |
| ENS ET-VAR Xb-396 t=2.0 | 377 | 0.45 | 1 |
| ENS ET-MEAN Xc-400m t=2.0 | 1.44 | 0.149 | 1 |
| ENS ET-VAR Xc-400m t=2.0 | 447 | 0.0909 | 1 |
| ENS ET-MEAN Xc-400 t=2.0 | 0.842 | 0.4 | 1 |
| ENS ET-VAR Xc-400 t=2.0 | 422 | 0.399 | 1 |
| ENS ET-MEAN Xd-396 t=2.0 | -0.222 | 0.826 | 1 |
| ENS ET-VAR Xd-396 t=2.0 | 383 | 0.588 | 1 |
| ENS ET-MEAN Xe-396 t=2.0 | 0.502 | 0.615 | 1 |
| ENS ET-VAR Xe-396 t=2.0 | 385 | 0.829 | 1 |
| ENS ET-MEAN Xf-392 t=2.0 | -1.67 | 0.0955 | 1 |
| ENS ET-VAR Xf-392 t=2.0 | 372 | 0.362 | 1 |
| ENS ET-MEAN Xa-400 t=8.0 | -0.637 | 0.525 | 1 |
| ENS ET-VAR Xa-400 t=8.0 | 440 | 0.14 | 1 |
| ENS ET-MEAN Xb-396 t=8.0 | -0.654 | 0.517 | 1 |
| ENS ET-VAR Xb-396 t=8.0 | 378 | 0.621 | 1 |
| ENS ET-MEAN Xc-400m t=8.0 | 0.598 | 0.55 | 1 |
| ENS ET-VAR Xc-400m t=8.0 | 424 | 0.267 | 1 |
| ENS ET-MEAN Xc-400 t=8.0 | -0.598 | 0.551 | 1 |
| ENS ET-VAR Xc-400 t=8.0 | 413 | 0.601 | 1 |
| ENS ET-MEAN Xd-396 t=8.0 | 2.05 | 0.0406 | 1 |
| ENS ET-VAR Xd-396 t=8.0 | 383 | 0.612 | 1 |
| ENS ET-MEAN Xe-396 t=8.0 | 0.561 | 0.574 | 1 |
| ENS ET-VAR Xe-396 t=8.0 | 411 | 0.642 | 1 |
| ENS ET-MEAN Xf-392 t=8.0 | -1.14 | 0.253 | 1 |
| ENS ET-VAR Xf-392 t=8.0 | 363 | 0.198 | 1 |
| ENS ET-MEAN Xg-396 t=8.0 | 0.729 | 0.466 | 1 |
| ENS ET-VAR Xg-396 t=8.0 | 402 | 0.892 | 1 |
| ENS ET-MEAN Xh-395 t=8.0 | 0.241 | 0.809 | 1 |
| ENS ET-VAR Xh-395 t=8.0 | 355 | 0.179 | 1 |
| ENS ET-VARSUM Xa-400+Xb-396 t=2.0 | 409 | 0.705 | 1 |
| ENS ET-VARSUM Xa-400+Xc-400m t=2.0 | 420 | 0.45 | 1 |
| ENS ET-VARSUM Xa-400+Xc-400 t=2.0 | 363 | 0.191 | 1 |
| ENS ET-VARSUM Xa-400+Xd-396 t=2.0 | 381 | 0.521 | 1 |
| ENS ET-VARSUM Xa-400+Xe-396 t=2.0 | 432 | 0.243 | 1 |
| ENS ET-VARSUM Xa-400+Xf-392 t=2.0 | 460 | 0.0388 | 1 |
| ENS ET-VARSUM Xb-396+Xc-400m t=2.0 | 429 | 0.288 | 1 |
| ENS ET-VARSUM Xb-396+Xc-400 t=2.0 | 428 | 0.308 | 1 |
| ENS ET-VARSUM Xb-396+Xd-396 t=2.0 | 395 | 0.912 | 1 |
| ENS ET-VARSUM Xb-396+Xe-396 t=2.0 | 374 | 0.384 | 1 |
| ENS ET-VARSUM Xb-396+Xf-392 t=2.0 | 341 | 0.0313 | 1 |
| ENS ET-VARSUM Xc-400m+Xc-400 t=2.0 | 397 | 0.963 | 1 |
| ENS ET-VARSUM Xc-400m+Xd-396 t=2.0 | 415 | 0.563 | 1 |
| ENS ET-VARSUM Xc-400m+Xe-396 t=2.0 | 459 | 0.0396 | 1 |
| ENS ET-VARSUM Xc-400m+Xf-392 t=2.0 | 378 | 0.47 | 1 |
| ENS ET-VARSUM Xc-400+Xd-396 t=2.0 | 456 | 0.0516 | 1 |
| ENS ET-VARSUM Xc-400+Xe-396 t=2.0 | 398 | 0.981 | 1 |
| ENS ET-VARSUM Xc-400+Xf-392 t=2.0 | 427 | 0.324 | 1 |
| ENS ET-VARSUM Xd-396+Xe-396 t=2.0 | 369 | 0.281 | 1 |
| ENS ET-VARSUM Xd-396+Xf-392 t=2.0 | 367 | 0.255 | 1 |
| ENS ET-VARSUM Xe-396+Xf-392 t=2.0 | 361 | 0.176 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xa-400@8 | 445 | 0.11 | 1 |
| ENS ET-VARSUM2 Xb-396@2+Xb-396@8 | 392 | 0.832 | 1 |
| ENS ET-VARSUM2 Xc-400m@2+Xc-400m@8 | 485 | 0.00379 | 0.84 |
| ENS ET-VARSUM2 Xc-400@2+Xc-400@8 | 392 | 0.82 | 1 |
| ENS ET-VARSUM2 Xd-396@2+Xd-396@8 | 422 | 0.403 | 1 |
| ENS ET-VARSUM2 Xe-396@2+Xe-396@8 | 414 | 0.568 | 1 |
| ENS ET-VARSUM2 Xf-392@2+Xf-392@8 | 358 | 0.135 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xb-396@8 | 424 | 0.38 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xc-400m@8 | 427 | 0.322 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xc-400@8 | 396 | 0.923 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xd-396@8 | 394 | 0.888 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xe-396@8 | 404 | 0.833 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xf-392@8 | 411 | 0.655 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xg-396@8 | 390 | 0.776 | 1 |
| ENS ET-VARSUM2 Xa-400@2+Xh-395@8 | 410 | 0.684 | 1 |
| COMP CT-HOMO t=0.5 | 7.93 | 0.339 | 1 |
| COMP CT-HOMO t=1.0 | 4.48 | 0.812 | 1 |
| COMP CT-HOMO t=2.0 | 3.62 | 0.89 | 1 |
| COMP CT-HOMO t=4.0 | 6.83 | 0.555 | 1 |
| COMP CT-HOMO t=8.0 | 7.33 | 0.502 | 1 |
| COMP CT-HOMO t=12.0 | 3.15 | 0.924 | 1 |
| COMP CT-VAR Xa-400 t=2.0 | 0.891 | 0.248 | 1 |
| COMP CT-VAR Xb-396 t=2.0 | 1.11 | 0.288 | 1 |
| COMP CT-VAR Xc-400m t=2.0 | 0.789 | 0.0184 | 1 |
| COMP CT-VAR Xc-400 t=2.0 | 0.938 | 0.522 | 1 |
| COMP CT-VAR Xd-396 t=2.0 | 0.93 | 0.468 | 1 |
| COMP CT-VAR Xe-396 t=2.0 | 0.99 | 0.918 | 1 |
| COMP CT-VAR Xf-392 t=2.0 | 1.11 | 0.306 | 1 |
| COMP CT-VARSUM Xa-400+Xb-396 t=2.0 | 0.896 | 0.272 | 1 |
| COMP CT-VARSUM Xa-400+Xc-400m t=2.0 | 0.936 | 0.51 | 1 |
| COMP CT-VARSUM Xa-400+Xc-400 t=2.0 | 0.979 | 0.832 | 1 |
| COMP CT-VARSUM Xa-400+Xd-396 t=2.0 | 1.09 | 0.37 | 1 |
| COMP CT-VARSUM Xa-400+Xe-396 t=2.0 | 0.839 | 0.079 | 1 |
| COMP CT-VARSUM Xa-400+Xf-392 t=2.0 | 0.814 | 0.0398 | 1 |
| COMP CT-VARSUM Xb-396+Xc-400m t=2.0 | 0.931 | 0.477 | 1 |
| COMP CT-VARSUM Xb-396+Xc-400 t=2.0 | 0.96 | 0.681 | 1 |
| COMP CT-VARSUM Xb-396+Xd-396 t=2.0 | 0.901 | 0.296 | 1 |
| COMP CT-VARSUM Xb-396+Xe-396 t=2.0 | 1.16 | 0.128 | 1 |
| COMP CT-VARSUM Xb-396+Xf-392 t=2.0 | 1.24 | 0.0342 | 1 |
| COMP CT-VARSUM Xc-400m+Xc-400 t=2.0 | 0.892 | 0.255 | 1 |
| COMP CT-VARSUM Xc-400m+Xd-396 t=2.0 | 0.841 | 0.0844 | 1 |
| COMP CT-VARSUM Xc-400m+Xe-396 t=2.0 | 0.752 | 0.00458 | 1 |
| COMP CT-VARSUM Xc-400m+Xf-392 t=2.0 | 1 | 0.961 | 1 |
| COMP CT-VARSUM Xc-400+Xd-396 t=2.0 | 0.831 | 0.0649 | 1 |
| COMP CT-VARSUM Xc-400+Xe-396 t=2.0 | 1.07 | 0.486 | 1 |
| COMP CT-VARSUM Xc-400+Xf-392 t=2.0 | 0.983 | 0.862 | 1 |
| COMP CT-VARSUM Xd-396+Xe-396 t=2.0 | 0.958 | 0.668 | 1 |
| COMP CT-VARSUM Xd-396+Xf-392 t=2.0 | 0.957 | 0.659 | 1 |
| COMP CT-VARSUM Xe-396+Xf-392 t=2.0 | 1.11 | 0.3 | 1 |
| COMP CT-DECAY Xa-400 | 8.77 | 0.187 | 1 |
| COMP CT-DECAY Xb-396 | 4.22 | 0.647 | 1 |
| COMP CT-DECAY Xc-400m | 1.11 | 0.981 | 1 |
| COMP CT-DECAY Xc-400 | 3.38 | 0.759 | 1 |
| COMP CT-DECAY Xd-396 | 6.91 | 0.33 | 1 |
| COMP CT-DECAY Xe-396 | 5.62 | 0.468 | 1 |
| COMP CT-PHOTON Xa-400/alpha->Xb-396 | 0.11 | 0.74 | 1 |
| COMP CT-PHOTON Xa-400/beta-->Xc-400m | 1.47 | 0.226 | 1 |
| COMP CT-PHOTON Xb-396/ec->Xd-396 | 0.375 | 0.54 | 1 |
| COMP CT-PHOTON Xc-400m/it->Xc-400 | 3.76 | 0.0524 | 1 |
| COMP CT-PHOTON Xe-396/beta-->Xg-396 | 1.07 | 0.3 | 1 |
| BIG 1000000000000 BIG-STATE t=0.5 | 9.26 | 0.321 | 1 |
| BIG 1000000000000 BIG-STATE t=1.0 | 10.1 | 0.255 | 1 |
| BIG 1000000000000 BIG-STATE t=2.0 | 12.9 | 0.117 | 1 |
| BIG 1000000000000 BIG-STATE t=4.0 | 4.46 | 0.813 | 1 |
| BIG 1000000000000 BIG-STATE t=8.0 | 9.15 | 0.33 | 1 |
| BIG 1000000000000 BIG-STATE t=12.0 | 6.24 | 0.62 | 1 |
| BIG 1000000000000 BIG-DECAY Xa-400 | 0.485 | 0.628 | 1 |
| BIG 1000000000000 BIG-DECAY Xb-396 | -1.37 | 0.172 | 1 |
| BIG 1000000000000 BIG-DECAY Xc-400m | 1.56 | 0.119 | 1 |
| BIG 1000000000000 BIG-DECAY Xc-400 | 0.958 | 0.338 | 1 |
| BIG 1000000000000 BIG-DECAY Xd-396 | -1.33 | 0.183 | 1 |
| BIG 1000000000000 BIG-DECAY Xe-396 | 1.63 | 0.102 | 1 |
| BIG 1152921504606859321 BIG-STATE t=0.5 | 2.91 | 0.94 | 1 |
| BIG 1152921504606859321 BIG-STATE t=1.0 | 2.98 | 0.936 | 1 |
| BIG 1152921504606859321 BIG-STATE t=2.0 | 7.51 | 0.483 | 1 |
| BIG 1152921504606859321 BIG-STATE t=4.0 | 8.7 | 0.369 | 1 |
| BIG 1152921504606859321 BIG-STATE t=8.0 | 8.75 | 0.363 | 1 |
| BIG 1152921504606859321 BIG-STATE t=12.0 | 5.73 | 0.678 | 1 |
| BIG 1152921504606859321 BIG-DECAY Xa-400 | 0.945 | 0.345 | 1 |
| BIG 1152921504606859321 BIG-DECAY Xb-396 | -1.11 | 0.266 | 1 |
| BIG 1152921504606859321 BIG-DECAY Xc-400m | 1.46 | 0.144 | 1 |
| BIG 1152921504606859321 BIG-DECAY Xc-400 | 1.74 | 0.0815 | 1 |
| BIG 1152921504606859321 BIG-DECAY Xd-396 | -1.25 | 0.211 | 1 |
| BIG 1152921504606859321 BIG-DECAY Xe-396 | 1.99 | 0.0461 | 1 |

### NC2: every test

| Test | Statistic (z for binomial and moment tests) | p | Holm-adjusted p |
|---|---:|---:|---:|
| NC2 S1 BT-STATE Sa-400 t=0.5 | 0.194 | 0.66 | 1 |
| NC2 S1 BT-STATE Sa-400 t=1.0 | 0.104 | 0.747 | 1 |
| NC2 S1 BT-STATE Sa-400 t=2.0 | 0.888 | 0.346 | 1 |
| NC2 S1 BT-STATE Sa-400 t=4.0 | 3.16 | 0.0756 | 1 |
| NC2 S1 BT-STATE Sa-400 t=6.0 | 0.086 | 0.769 | 1 |
| NC2 S1 BT-DECAY Sa-400 | -0.293 | 0.769 | 1 |
| NC2 S1 BT-BINS Sa-400 | 8.58 | 0.127 | 1 |
| NC2 S1 BT-PHOTON Sa-400/beta-->Sb-400 | -1.64 | 0.1 | 1 |
| NC2 S1 BT-PHOTON-T Sa-400/beta-->Sb-400 | 5.85 | 0.321 | 1 |
| NC2 S2 BT-STATE Xa-400 t=0.5 | 6.83 | 0.337 | 1 |
| NC2 S2 BT-STATE Xa-400 t=1.0 | 6.31 | 0.504 | 1 |
| NC2 S2 BT-STATE Xa-400 t=2.0 | 6.84 | 0.554 | 1 |
| NC2 S2 BT-STATE Xa-400 t=4.0 | 12.2 | 0.143 | 1 |
| NC2 S2 BT-STATE Xa-400 t=8.0 | 16 | 0.0418 | 1 |
| NC2 S2 BT-STATE Xa-400 t=12.0 | 10.3 | 0.242 | 1 |
| NC2 S2 BT-DECAY Xa-400 | 1.67 | 0.0974 | 1 |
| NC2 S2 BT-DECAY Xb-396 | -0.292 | 0.771 | 1 |
| NC2 S2 BT-DECAY Xc-400m | 0.535 | 0.594 | 1 |
| NC2 S2 BT-DECAY Xc-400 | -0.386 | 0.703 | 1 |
| NC2 S2 BT-DECAY Xd-396 | -0.437 | 0.662 | 1 |
| NC2 S2 BT-DECAY Xe-396 | -0.237 | 0.815 | 1 |
| NC2 S2 BT-BINS Xa-400 | 9.71 | 0.138 | 1 |
| NC2 S2 BT-CHANNEL Xa-400 | -0.61 | 0.541 | 1 |
| NC2 S2 BT-PHOTON Xa-400/alpha->Xb-396 | -0.612 | 0.544 | 1 |
| NC2 S2 BT-PHOTON Xa-400/beta-->Xc-400m | -0.983 | 0.328 | 1 |
| NC2 S2 BT-PHOTON Xb-396/ec->Xd-396 | 1.69 | 0.0913 | 1 |
| NC2 S2 BT-PHOTON Xc-400m/it->Xc-400 | 0.469 | 0.644 | 1 |
| NC2 S2 BT-CHANNEL Xe-396 | -0.949 | 0.342 | 1 |
| NC2 S2 BT-PHOTON Xe-396/beta-->Xg-396 | -0.923 | 0.359 | 1 |
| NC2 S2 BT-PHOTON-T Xa-400/alpha->Xb-396 | 4.54 | 0.604 | 1 |
| NC2 S2 BT-PHOTON-T Xa-400/beta-->Xc-400m | 4.95 | 0.551 | 1 |
| NC2 S2 BT-PHOTON-T Xb-396/ec->Xd-396 | 11 | 0.0885 | 1 |
| NC2 S2 BT-PHOTON-T Xc-400m/it->Xc-400 | 2.14 | 0.906 | 1 |
| NC2 S2 BT-PHOTON-T Xe-396/beta-->Xg-396 | 2.1 | 0.836 | 1 |
| NC2 S2 BT-PATH Xa-400 | 36.9 | 0.732 | 1 |
| NC2 S2 BT-PATH Xb-396 | 9.18 | 0.421 | 1 |
| NC2 S2 BT-PATH Xc-400m | 20.9 | 0.232 | 1 |
| NC2 S2 BT-PATH Xc-400 | 12.7 | 0.55 | 1 |
| NC2 S2 BT-PATH Xd-396 | 1.66 | 0.893 | 1 |
| NC2 S2 BT-PATH Xe-396 | 15.5 | 0.0304 | 1 |
| NC2 S3 BT-STATE Ya-410 t=1.0 | 1.07 | 0.586 | 1 |
| NC2 S3 BT-STATE Ya-410 t=5.0 | 0.759 | 0.684 | 1 |
| NC2 S3 BT-STATE Ya-410 t=10.0 | 1.06 | 0.588 | 1 |
| NC2 S3 BT-STATE Ya-410 t=20.0 | 5.06 | 0.0797 | 1 |
| NC2 S3 BT-STATE Ya-410 t=40.0 | 1.58 | 0.454 | 1 |
| NC2 S3 BT-DECAY Ya-410 | 1.25 | 0.212 | 1 |
| NC2 S3 BT-DECAY Yb-410 | 1.21 | 0.228 | 1 |
| NC2 S3 BT-BINS Ya-410 | 8.06 | 0.153 | 1 |
| NC2 S4 BT-STATE G1 t=0.5 | 0.538 | 0.764 | 1 |
| NC2 S4 BT-STATE G1 t=1.0 | 0.25 | 0.883 | 1 |
| NC2 S4 BT-STATE G1 t=2.0 | 2.2 | 0.332 | 1 |
| NC2 S4 BT-STATE G1 t=5.0 | 1.61 | 0.446 | 1 |
| NC2 S4 BT-STATE G1 t=13.0 | 0.0476 | 0.827 | 1 |
| NC2 S4 BT-DECAY G1 | 0.11 | 1 | 1 |
| NC2 S4 BT-BINS G1 | 4.72 | 0.318 | 1 |
| NC2 S4 BT-NEUTRON G1 | 7.8 | 0.167 | 1 |
| NC2 S4 BT-STATE G2 t=0.5 | 0.667 | 0.717 | 1 |
| NC2 S4 BT-STATE G2 t=1.0 | 1.02 | 0.6 | 1 |
| NC2 S4 BT-STATE G2 t=2.0 | 0.362 | 0.835 | 1 |
| NC2 S4 BT-STATE G2 t=5.0 | 6.99 | 0.0303 | 1 |
| NC2 S4 BT-STATE G2 t=13.0 | 0.478 | 0.788 | 1 |
| NC2 S4 BT-DECAY G2 | -0.302 | 0.759 | 1 |
| NC2 S4 BT-BINS G2 | 8.01 | 0.156 | 1 |
| NC2 S4 BT-NEUTRON G2 | 2.26 | 0.812 | 1 |
| NC2 S4 BT-STATE G3 t=0.5 | 0.0281 | 0.867 | 1 |
| NC2 S4 BT-STATE G3 t=1.0 | 0.0307 | 0.861 | 1 |
| NC2 S4 BT-STATE G3 t=2.0 | 0.00398 | 0.95 | 1 |
| NC2 S4 BT-STATE G3 t=5.0 | 0.37 | 0.543 | 1 |
| NC2 S4 BT-STATE G3 t=13.0 | 0.635 | 0.426 | 1 |
| NC2 S4 BT-DECAY G3 | -0.797 | 0.428 | 1 |
| NC2 S4 BT-BINS G3 | 1.08 | 0.956 | 1 |
| NC2 S4 BT-NEUTRON G3 | 1.08 | 0.956 | 1 |
| NC2 S4 BT-CHANNEL G1 | 0.218 | 0.827 | 1 |
| NC2 S4 BT-CHANNEL G2 | 0.622 | 0.532 | 1 |
| NC2 S5 BT-STATE Za-450 t=1.0 | 0.0241 | 0.877 | 1 |
| NC2 S5 BT-STATE Za-450 t=1000.0 | 0.785 | 0.94 | 1 |
| NC2 S5 BT-STATE Za-450 t=3000.0 | 2.67 | 0.614 | 1 |
| NC2 S5 BT-DECAY Za-450 | -1.14 | 0.255 | 1 |
| NC2 S5 BT-DECAY Zb-450 | -1.14 | 0.255 | 1 |
| NC2 S5 BT-DECAY Zc-446 | -1.03 | 0.306 | 1 |
| NC2 S5 BT-DECAY Zd-446 | -1.07 | 0.284 | 1 |
| NC2 S5 BT-BINS Za-450 | 2.85 | 0.415 | 1 |
| NC2 S5 BT-PHOTON Zb-450/alpha->Zc-446 | -1.01 | 0.313 | 1 |
| NC2 S5 BT-CHANNEL Zd-446 | 0.406 | 0.687 | 1 |
| NC2 S5 BT-PHOTON-T Zb-450/alpha->Zc-446 | 3.81 | 0.282 | 1 |

Tests not applicable (fewer than two cells after pooling): NC1 S5 BT-STATE Za-450 t=0.001, NC2 S5 BT-STATE Za-450 t=0.001.

Planned tests that could not be computed (failed items): none.

## R0: calibration (complete null batteries)

`results/calibration.json`: 50000 batteries, namespace `rnglab-v020-nuclear-b2-calibration-record`, code `25ff5516b2dd88518e478ce4064d3e50ae0f00b6` (committed tree; files of earlier attempts: 0; ignored files allowed (bytecode caches, smoke/, egg-info): 2; other uncommitted or ignored files: none); attempt 1, deviation: none; replay check: 0 chunks of earlier attempts compared, 0 mismatches. C1 (per family) holds: **True**.

| Family | Tests | Batteries that FAIL under the null | Rate | 95% interval | C1 |
|---|---:|---:|---:|---|---|
| NC1 | 222 | 440 | 0.0088 | 0.0080 to 0.0097 | True |
| NC2 | 85 | 399 | 0.0080 | 0.0072 to 0.0088 | True |
| R1 (NC1 or NC2) | | 835 | 0.0167 | 0.0156 to 0.0179 | reported, not a criterion |

### NC1: test types at the Holm level (4.50e-05)

| Type | Tests | Batteries with a rejection | Share of the union bound | 95% interval |
|---|---:|---:|---:|---|
| BIG-DECAY | 12 | 15 | 0.55 | 0.31 to 0.92 |
| BIG-STATE | 12 | 36 | 1.33 | 0.93 to 1.84 |
| BT-BINS | 7 | 21 | 1.33 | 0.82 to 2.04 |
| BT-CHANNEL | 5 | 13 | 1.15 | 0.61 to 1.97 |
| BT-DECAY | 16 | 24 | 0.67 | 0.43 to 0.99 |
| BT-NEUTRON | 3 | 17 | 2.52 | 1.47 to 4.03 |
| BT-PATH | 6 | 19 | 1.41 | 0.85 to 2.20 |
| BT-PHOTON | 7 | 20 | 1.27 | 0.77 to 1.96 |
| BT-PHOTON-T | 7 | 19 | 1.21 | 0.73 to 1.88 |
| BT-STATE | 34 | 80 | 1.04 | 0.83 to 1.30 |
| CT-DECAY | 6 | 14 | 1.04 | 0.57 to 1.74 |
| CT-HOMO | 6 | 13 | 0.96 | 0.51 to 1.64 |
| CT-PHOTON | 5 | 9 | 0.80 | 0.37 to 1.52 |
| CT-VAR | 7 | 12 | 0.76 | 0.39 to 1.33 |
| CT-VARSUM | 21 | 44 | 0.93 | 0.68 to 1.25 |
| ET-MEAN | 16 | 33 | 0.92 | 0.63 to 1.29 |
| ET-VAR | 16 | 41 | 1.14 | 0.82 to 1.54 |
| ET-VARSUM | 21 | 46 | 0.97 | 0.71 to 1.30 |
| ET-VARSUM2 | 15 | 17 | 0.50 | 0.29 to 0.81 |

### NC2: test types at the Holm level (1.18e-04)

| Type | Tests | Batteries with a rejection | Share of the union bound | 95% interval |
|---|---:|---:|---:|---|
| BT-BINS | 7 | 37 | 0.90 | 0.63 to 1.24 |
| BT-CHANNEL | 5 | 33 | 1.12 | 0.77 to 1.58 |
| BT-DECAY | 16 | 66 | 0.70 | 0.54 to 0.89 |
| BT-NEUTRON | 3 | 17 | 0.96 | 0.56 to 1.54 |
| BT-PATH | 6 | 35 | 0.99 | 0.69 to 1.38 |
| BT-PHOTON | 7 | 31 | 0.75 | 0.51 to 1.07 |
| BT-PHOTON-T | 7 | 42 | 1.02 | 0.74 to 1.38 |
| BT-STATE | 34 | 201 | 1.01 | 0.87 to 1.15 |

## Positive controls

Level: p < 4.50e-05.

| Control | Target | p | Rejected |
|---|---|---:|---|
| PC-B1 | BT-STATE Sa-400 t=2.0 | 2.19e-20 | True |
| PC-B2 single run | BT-STATE Xa-400 t=12.0 | 0 | True |
| PC-B2 comparison | CT-HOMO t=12.0, ENS-A against ENS-NAIVE | 0 | True |

PC-B1, every test (descriptive, outside the families):

| Test | p |
|---|---:|
| PC-B1 S1 BT-STATE Sa-400 t=0.5 | 1.12e-11 |
| PC-B1 S1 BT-STATE Sa-400 t=1.0 | 9.21e-20 |
| PC-B1 S1 BT-STATE Sa-400 t=2.0 | 2.19e-20 |
| PC-B1 S1 BT-STATE Sa-400 t=4.0 | 1.69e-15 |
| PC-B1 S1 BT-STATE Sa-400 t=6.0 | 1.03e-12 |
| PC-B1 S1 BT-DECAY Sa-400 | 2.13e-12 |
| PC-B1 S1 BT-BINS Sa-400 | 1.52e-28 |
| PC-B1 S1 BT-PHOTON Sa-400/beta-->Sb-400 | 0.688 |
| PC-B1 S1 BT-PHOTON-T Sa-400/beta-->Sb-400 | 0.0694 |

PC-B2 single run, every test (descriptive, outside the families):

| Test | p |
|---|---:|
| PC-B2 S2 BT-STATE Xa-400 t=0.5 | 0 |
| PC-B2 S2 BT-STATE Xa-400 t=1.0 | 0 |
| PC-B2 S2 BT-STATE Xa-400 t=2.0 | 0 |
| PC-B2 S2 BT-STATE Xa-400 t=4.0 | 0 |
| PC-B2 S2 BT-STATE Xa-400 t=8.0 | 0 |
| PC-B2 S2 BT-STATE Xa-400 t=12.0 | 0 |

## Huge populations

| Initial nuclei | Parts | Seconds | Final populations |
|---:|---:|---:|---|
| 1000000000000 | 1 | 0.10 | [15624939867, 3333306787, 3742666453, 63930753258, 19931420780, 48371961057, 606734597529, 214497391409, 23832962860] |
| 1152921504606859321 | 1025 | 21.35 | [18014398383629835, 3843057020713582, 4315011348429696, 73706781418711248, 22979375576363809, 55769377545493112, 699518114726262191, 247297849646569626, 27477538940686222] |

## Single-machine timings (seconds, wall clock; not HPC measurements)

| Run | Seconds |
|---|---:|
| NC1 S1 | 0.007 |
| NC1 S2 | 0.041 |
| NC1 S3 | 0.009 |
| NC1 S4 | 0.011 |
| NC1 S5 | 0.019 |
| NC2 S1 | 0.006 |
| NC2 S2 | 0.036 |
| NC2 S3 | 0.010 |
| NC2 S4 | 0.009 |
| NC2 S5 | 0.020 |
| ensembles | 17.213 |
| BIG 1000000000000 | 0.098 |
| BIG 1152921504606859321 | 21.345 |
