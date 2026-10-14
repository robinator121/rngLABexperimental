"""Domain-separated reproducible seeds and checks against partition leakage.

Distinct derived seeds support a reproducible design; hashing alone does not
prove physical independence or turn a deterministic generator into a CSPRNG.
"""
from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping

SPLITS = frozenset({'development', 'calibration', 'confirmation', 'final_holdout', 'pilot'})
_NAMESPACE = 'RNGLab-v0.12-seed-v1'


def _identity(split: str, application: str, group: str | int, stream: str) -> dict:
    if split not in SPLITS:
        raise ValueError(f'Unknown seed partition: {split!r}')
    if not isinstance(application, str) or not application.strip():
        raise ValueError('A nonempty application identity is required')
    if not isinstance(stream, str) or not stream.strip():
        raise ValueError('A nonempty stream identity is required')
    if isinstance(group, bool) or not isinstance(group, (str, int)) or not str(group):
        raise ValueError('A nonempty string or integer group identity is required')
    return {'namespace': _NAMESPACE, 'split': split, 'application': application,
            'group': str(group), 'stream': stream}


def _digest(identity: dict) -> str:
    return hashlib.sha256(json.dumps(identity, sort_keys=True, separators=(',', ':'),
                                     ensure_ascii=True).encode('ascii')).hexdigest()


def derive_seed(split: str, application: str, group: str | int, stream: str = 'source') -> int:
    """Return a 128-bit seed with explicit partition/application/group/stream domains."""
    return int(_digest(_identity(split, application, group, stream))[:32], 16)


def seed_lineage(split: str, application: str, group: str | int, stream: str = 'source') -> dict:
    identity = _identity(split, application, group, stream)
    digest = _digest(identity)
    return {**identity, 'seed': int(digest[:32], 16), 'lineage_id': digest}


def assert_disjoint_lineages(*partitions: Iterable[Mapping]) -> None:
    """Reject duplicate group identities or numeric seeds, even if IDs were relabeled.

    Callers must retain original lineage metadata when transforming a trace.
    An intentionally paired condition uses one group once, not separate entries
    purporting to be independent calibration or validation groups.
    """
    seen_ids: set[str] = set()
    seen_seeds: set[int] = set()
    for partition in partitions:
        for record in partition:
            if record.get('split') not in SPLITS:
                raise ValueError('Every lineage needs an authorized seed partition')
            lineage = record.get('lineage_id')
            if not isinstance(lineage, str) or not lineage:
                raise ValueError('Every lineage needs a nonempty lineage_id')
            if lineage in seen_ids:
                raise ValueError('Reused lineage: independent groups/partitions must be disjoint')
            seen_ids.add(lineage)
            seed = record.get('seed')
            if seed is not None:
                if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
                    raise ValueError('Seeds must be nonnegative integers')
                if seed in seen_seeds:
                    raise ValueError('Reused numeric seed, including across relabeled partitions')
                seen_seeds.add(seed)
