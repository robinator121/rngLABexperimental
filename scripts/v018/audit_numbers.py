#!/usr/bin/env python3
"""Index the numbers quoted in documentation against saved machine outputs.

RNGLab v0.18 audit helper (written by Claude Code, 2026-11-05). For every
numeric token in README.md and docs/**/*.md, search the text outputs (JSON,
CSV, TXT, LOG, XML) of the campaign that the document describes for a value
that rounds to the quoted value at the quoted precision. Percentages are also
matched as fractions.

Each document is searched only in its own evidence group (SCOPES below), so a
number cannot be "backed" by an unrelated table that happens to contain it.
A decoy control repeats every lookup with the value multiplied by 1.0371; the
decoy match rate estimates how often a number of that precision matches the
same evidence by chance. A real match rate close to the decoy rate is
uninformative.

A match shows only that the number occurs in a saved output of the right
campaign. It does not show that the output supports the sentence around it.

Usage (from the repository root):
    python scripts/v018/audit_numbers.py --out docs/v018/evidence/audit_numbers.csv \
        [--external GROUP=DIR ...]

External directories are outside Git; their matches are reported in a separate
column because they cannot be reproduced from the repository.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import json
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
TEXT_SUFFIXES = {".json", ".csv", ".txt", ".log", ".xml", ".tsv"}
MAX_EVIDENCE_BYTES = 64 * 1024 * 1024
DECOY_FACTOR = 1.0371

# Evidence group -> repository paths (files or directories) holding its outputs.
GROUPS = {
    "v014": ["runs/v014-integrated-experiment"],
    "v016": [],
    "v017": ["docs/calibrate"],
    "phase2": ["research/phase2-20261030", "research/phase1-20261030"],
    "crypto": ["research/crypto-blue"],
    "detector": ["docs/v018/evidence/detector_inventory.json",
                 "docs/v018/evidence/detector_ceiling.json"],
    "benchmarks": [],
    "tests": ["docs/v018/evidence/pytest_baseline.txt"],
}
# Document path prefix -> evidence groups searched (first matching prefix wins).
SCOPES = [
    ("README.md", ["v014", "v016", "v017", "phase2", "crypto", "detector", "benchmarks", "tests"]),
    ("docs/calibrate/", ["v017", "v016"]),
    ("docs/trident/", ["v016"]),
    ("docs/crypto/", ["crypto", "tests"]),
    ("docs/physical/", ["detector"]),
    ("docs/benchmarks/", ["benchmarks"]),
    ("docs/", ["v014"]),
]

NUMBER = re.compile(r"(?<![\w.])[-−]?(?:\d{1,3}(?:,\d{3})+|\d+)?(?:\.\d+)?(?:[eE][-+]?\d+)?%?")
SKIP_LINE = re.compile(r"Document Identifier|Provenance SHA-256|^\*\*Date|SHA-256:|sha256|^\s*\|\s*:?-{3}")
SKIP_SPAN = re.compile(r"`[^`]*`|\]\([^)]*\)|https?://\S+|\b[0-9a-f]{12,}\b|"
                       r"\b(?:19|20)\d\d-\d\d-\d\d\b|\bv\d+(?:\.\d+)+\b|\b\d\d:\d\d(?::\d\d)?\b")


def doc_files() -> list[Path]:
    files = [REPO / "README.md"]
    files += sorted(p for p in (REPO / "docs").rglob("*.md") if "v018" not in p.parts)
    return files


def scope_for(rel: str) -> list[str]:
    return next(groups for prefix, groups in SCOPES if rel.startswith(prefix))


def text_files(root: Path) -> list[Path]:
    paths = [root] if root.is_file() else list(root.rglob("*")) if root.is_dir() else []
    return sorted(p for p in paths if p.is_file() and p.suffix.lower() in TEXT_SUFFIXES
                  and p.stat().st_size <= MAX_EVIDENCE_BYTES)


def parse_number(token: str) -> tuple[float, int, bool] | None:
    raw = token.replace("−", "-").replace(",", "")
    percent = raw.endswith("%")
    raw = raw.rstrip("%")
    if not raw or raw in {"-", "."} or not any(c.isdigit() for c in raw):
        return None
    try:
        value = float(raw)
    except ValueError:
        return None
    mantissa = raw.split("e")[0].split("E")[0]
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    return value, decimals, percent


class EvidenceIndex:
    def __init__(self, files: list[Path], base: Path):
        values: list[tuple[float, str]] = []
        for path in files:
            try:
                text = path.read_text(errors="ignore")
            except OSError:
                continue
            label = str(path.relative_to(base)) if path.is_relative_to(base) else str(path)
            for match in NUMBER.finditer(text):
                parsed = parse_number(match.group(0))
                if parsed is not None:
                    values.append((parsed[0], label))
        values.sort(key=lambda item: item[0])
        self.keys = [v for v, _ in values]
        self.labels = [label for _, label in values]
        self.files = len(files)

    def lookup(self, value: float, decimals: int, percent: bool) -> list[str]:
        found = set()
        for v, d in [(value, decimals)] + ([(value / 100.0, decimals + 2)] if percent else []):
            half = 0.5 * 10.0 ** -d
            lo = bisect.bisect_left(self.keys, v - half - 1e-12)
            hi = bisect.bisect_right(self.keys, v + half + 1e-12)
            found.update(self.labels[lo:hi])
        return sorted(found)


def discriminating(value: float, decimals: int, percent: bool) -> bool:
    """Small whole numbers and round values match almost any output by chance."""
    if decimals == 0 and not percent:
        return abs(value) >= 100 and value % 100 != 0
    significant = len(f"{abs(value):.{decimals}f}".replace(".", "").lstrip("0"))
    return significant >= 3


def tokens(line: str):
    cleaned = SKIP_SPAN.sub(" ", line)
    for match in NUMBER.finditer(cleaned):
        parsed = parse_number(match.group(0))
        if parsed is None:
            continue
        value, decimals, percent = parsed
        if decimals == 0 and not percent and 1900 <= value <= 2100:
            continue
        yield match.group(0), value, decimals, percent


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--external", action="append", default=[], metavar="GROUP=DIR")
    args = parser.parse_args(argv)

    repo_index = {name: EvidenceIndex([f for rel in paths for f in text_files(REPO / rel)], REPO)
                  for name, paths in GROUPS.items()}
    external_roots: dict[str, list[Path]] = {}
    for item in args.external:
        group, _, directory = item.partition("=")
        if group not in GROUPS or not directory:
            parser.error(f"--external needs GROUP=DIR with GROUP in {sorted(GROUPS)}")
        external_roots.setdefault(group, []).append(Path(directory))
    external_index = {group: [(root, EvidenceIndex(text_files(root), root)) for root in roots]
                      for group, roots in external_roots.items()}

    rows, summary = [], {}
    for doc in doc_files():
        rel = str(doc.relative_to(REPO))
        groups = scope_for(rel)
        stats = summary.setdefault(rel, {"evidence_groups": groups, "tokens": 0, "discriminating": 0,
                                         "repo": 0, "external_only": 0, "unmatched": 0,
                                         "decoy_repo_or_external": 0})
        for lineno, line in enumerate(doc.read_text(errors="ignore").splitlines(), 1):
            if SKIP_LINE.search(line):
                continue
            for text, value, decimals, percent in tokens(line):
                def search(v):
                    repo_hits = [f"{g}:{h}" for g in groups for h in repo_index[g].lookup(v, decimals, percent)]
                    ext_hits = [f"{g}:{root.name}/{h}" for g in groups
                                for root, index in external_index.get(g, [])
                                for h in index.lookup(v, decimals, percent)]
                    return repo_hits, ext_hits
                hits, ext_hits = search(value)
                decoy_hits, decoy_ext = search(value * DECOY_FACTOR)
                disc = discriminating(value, decimals, percent)
                stats["tokens"] += 1
                if disc:
                    stats["discriminating"] += 1
                    stats["repo"] += bool(hits)
                    stats["external_only"] += bool(ext_hits) and not hits
                    stats["unmatched"] += not hits and not ext_hits
                    stats["decoy_repo_or_external"] += bool(decoy_hits or decoy_ext)
                rows.append({"doc": rel, "line": lineno, "token": text, "decimals": decimals,
                             "discriminating": disc, "repo_match": ";".join(hits[:3]),
                             "external_match": ";".join(ext_hits[:3]),
                             "decoy_match": bool(decoy_hits or decoy_ext)})
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    args.out.with_suffix(".summary.json").write_text(json.dumps({
        "script": "scripts/v018/audit_numbers.py",
        "decoy_factor": DECOY_FACTOR,
        "repo_evidence_files": {name: index.files for name, index in repo_index.items()},
        "external_evidence_files": {group: {str(root): index.files for root, index in pairs}
                                    for group, pairs in external_index.items()},
        "definition": ("discriminating = whole number >=100 and not a multiple of 100, or a value "
                       "with >=3 significant digits, or any percentage. repo / external_only / "
                       "unmatched partition the discriminating tokens. decoy counts discriminating "
                       "tokens whose value*decoy_factor also matches the same evidence groups."),
        "documents": summary}, indent=2) + "\n")
    print(f"{len(rows)} tokens from {len(summary)} documents -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
