#!/usr/bin/env python3
"""Stream an export in the assumed field0, field1, waveform... text layout.

Python standard library only. Does not validate entropy or identify detector fields.
Usage: python3 inspect_export.py '/path/to/export.txt' > inventory.json
"""
import argparse
import collections
import hashlib
import json
import math
from pathlib import Path


def inspect(path):
    before = path.stat()
    digest = hashlib.sha256()
    widths = collections.Counter()
    leading = collections.Counter()
    endings = collections.Counter()
    errors = []
    noninteger_count = too_short_count = blank_count = valid = samples = ceiling = floor = 0
    minimum = maximum = first_time = last_time = previous = gap_min = gap_max = None
    nonincreasing = gap_gcd = time_gcd = gap_sum = byte_count = 0
    for lineno, line in enumerate(path.open('rb'), 1):
        byte_count += len(line)
        digest.update(line)
        endings['CRLF' if line.endswith(b'\r\n') else 'LF' if line.endswith(b'\n') else 'none'] += 1
        tokens = line.split()
        widths[len(tokens)] += 1
        if not tokens:
            blank_count += 1
            continue
        try:
            # Reject decimal points, exponent notation and nonnumeric suffixes.
            values = [int(token) for token in tokens]
        except ValueError:
            noninteger_count += 1
            if len(errors) < 20:
                errors.append({'line': lineno, 'problem': 'not all integer tokens'})
            continue
        if len(values) < 3:
            too_short_count += 1
            if len(errors) < 20:
                errors.append({'line': lineno, 'problem': 'fewer than three fields'})
            continue
        valid += 1
        leading[values[0]] += 1
        timestamp = values[1]
        if first_time is None:
            first_time = timestamp
        time_gcd = math.gcd(time_gcd, abs(timestamp))
        if previous is not None:
            gap = timestamp - previous
            gap_gcd = math.gcd(gap_gcd, abs(gap))
            gap_sum += gap
            gap_min = gap if gap_min is None else min(gap_min, gap)
            gap_max = gap if gap_max is None else max(gap_max, gap)
            nonincreasing += gap <= 0
        previous = last_time = timestamp
        waveform = values[2:]
        minimum = min(waveform) if minimum is None else min(minimum, min(waveform))
        maximum = max(waveform) if maximum is None else max(maximum, max(waveform))
        samples += len(waveform)
        ceiling += waveform.count(65535)
        floor += waveform.count(0)
    after = path.stat()
    return {
        'path': str(path.resolve()), 'bytes_stat': before.st_size, 'bytes_read': byte_count,
        'decimal_MB': before.st_size / 1e6, 'binary_MiB': before.st_size / 2**20,
        'sha256': digest.hexdigest(), 'valid_numeric_records': valid,
        'fields_per_line': dict(widths), 'assumed_waveform_lengths': {k-2:v for k,v in widths.items() if k>=3},
        'blank_lines': blank_count, 'noninteger_lines': noninteger_count,
        'too_short_lines': too_short_count, 'first_20_errors': errors,
        'line_endings': dict(endings), 'field0_counts': dict(leading),
        'first_field1': first_time, 'last_field1': last_time,
        'field1_span': None if not valid else last_time-first_time,
        'field1_gcd': time_gcd, 'gap_gcd': gap_gcd,
        'gap_min': gap_min, 'gap_max': gap_max,
        'gap_mean': gap_sum/(valid-1) if valid>1 else None,
        'nonincreasing_steps': nonincreasing, 'waveform_values': samples,
        'waveform_min': minimum, 'waveform_max': maximum,
        'values_equal_65535': ceiling, 'values_equal_0': floor,
        'source_unchanged_during_read': before.st_size==after.st_size and before.st_mtime_ns==after.st_mtime_ns,
        'interpretation': 'Assumes two leading fields. No unit inference, entropy estimate or TRNG validation.'
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('file', type=Path)
    args = parser.parse_args()
    print(json.dumps(inspect(args.file), indent=2))
