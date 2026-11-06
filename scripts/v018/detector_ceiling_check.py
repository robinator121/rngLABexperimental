#!/usr/bin/env python3
"""Per-record ceiling and trigger-alignment counts for the detector export.

RNGLab v0.18 audit helper (written by Claude Code, 2026-11-05). It recomputes
the per-record figures quoted in docs/physical/MUON_DETECTOR_ASSESSMENT.md
section 2, which scripts/inspect_detector_export.py does not report. Same
assumed layout: field 0, field 1, then waveform codes. Standard library only.
No entropy, unit or particle-identity inference.

Usage: python3 scripts/v018/detector_ceiling_check.py EXPORT.txt > out.json
"""
import argparse
import collections
import hashlib
import json
from pathlib import Path

CEILING = 65535
TRIGGER_THRESHOLD = 10000
SEARCH = range(100, 160)
MAIN_WINDOW = range(110, 160)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("export", type=Path)
    path = parser.parse_args().export
    digest = hashlib.sha256()
    records = any_ceiling = window_ceiling = 0
    first_crossing = collections.Counter()
    ceiling_per_record = []
    for line in path.open("rb"):
        digest.update(line)
        tokens = line.split()
        if not tokens:
            continue
        wave = [int(token) for token in tokens[2:]]
        records += 1
        count = wave.count(CEILING)
        ceiling_per_record.append(count)
        any_ceiling += count > 0
        window_ceiling += any(wave[i] == CEILING for i in MAIN_WINDOW)
        crossing = next((i for i in SEARCH if wave[i] > TRIGGER_THRESHOLD), None)
        first_crossing[crossing] += 1
    ceiling_per_record.sort()
    result = {
        "script": "scripts/v018/detector_ceiling_check.py",
        "sha256": digest.hexdigest(),
        "records": records,
        "records_reaching_ceiling_anywhere": any_ceiling,
        "fraction_reaching_ceiling_anywhere": any_ceiling / records,
        "records_reaching_ceiling_in_samples_110_159": window_ceiling,
        "fraction_reaching_ceiling_in_samples_110_159": window_ceiling / records,
        "median_ceiling_samples_per_record": ceiling_per_record[records // 2],
        "first_index_above_10000_in_samples_100_159": {str(k): v for k, v in sorted(
            first_crossing.items(), key=lambda item: (item[0] is None, item[0] or 0))},
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
