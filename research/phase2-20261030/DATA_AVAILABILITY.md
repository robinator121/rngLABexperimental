# Data availability and provenance

## Distributed evidence

This repository contains the Phase II generation scripts and their tests, the two frozen confirmation/fault protocols, the exact healthy-decision file, complete seed manifests, per-condition results and replica endpoint tables, alarm/failure times, figures, mathematical checks, the RUS prefix and future outcomes, dependency pins, and recorded local verification evidence.

Generation code and locally frozen protocols are preserved byte for byte. The pre-run statements in the protocols are point-in-time design records, not claims that the experiments are still unexecuted. The fault protocol's prerequisite healthy-decision hash is preserved as well.

## Omitted large and local artifacts

Nineteen compressed Phase II trace archives, approximately 2 GiB in total with supporting material, remain outside Git. The public `MANIFEST.json` lists their relative filenames, original byte sizes and SHA-256 values as **not distributed**. There is no promised release, object-store mirror or download URL. Do not infer availability from a checksum.

The original Phase I compressed trace and two detector-file inventories are also omitted. Seven smaller historical Phase I files are published separately. Physical detector files and an entropy-source validation were not part of these new computational experiments. Virtual environments, interpreter caches, installation scratch logs and the redundant pristine-source tarball are not distributed.

The public preparer uses the local Git object for the pinned upstream commit instead of a tarball. [The reproduction guide](REPRODUCE.md) explains how to regenerate Phase II traces without needing the original Phase I inputs. Reusing the same seeds is reproducibility, not independent replication. A separate independent study needs new registered seed namespaces.

## Which checks can a reader run?

- `reproduction/verify_public_payload.py` verifies checksums of files actually included and frozen script/protocol/decision links. It does not run simulations or make a new coverage claim.
- The 35 standalone numerical and regression checks can run in the fresh prepared workspace described in the guide.
- Replaying the campaigns reconstructs new Phase II trace files under the fixed designs, subject to the recorded environment and numerical reproducibility limitations.
- The historical `experiment/validate_artifacts.py` requires the complete original Phase I inputs plus full Phase II archives. The recorded 132-pass result documents the original local audit; that command is not fully runnable from this lightweight publication alone. Regenerating Phase II does not recreate the missing Phase I history.

## Historical and public inventories

The root `MANIFEST.json` and `SHA256SUMS` describe this public Phase II directory. The files in `provenance/local-full-*` preserve the full original local inventory, including omitted traces and other local-only artifacts. Original component manifests and component `SHA256SUMS` files are historical full-experiment records and can refer to omitted files or the original unnormalized reports. Their existence does not mean every listed item is in Git.

Publication notes were prepended to public report copies. Private workspace prefixes in copies of documentation, logs and metadata were replaced with `LOCAL_WORKSPACE` or `LOCAL_HOME`. Those changes do not alter experiment measurements. The transformation record contains original and published SHA-256 values for each copied Phase II file. Historical source/generation hashes should be interpreted against their unchanged code or the original hash shown in that record, not silently substituted with a new protocol hash.

No new numerical campaign or package-API coverage study was executed merely to publish the results. The optional repository update is supported by the previously executed regression suite; its equivalent package-API coverage campaign remains the proposed next milestone.
