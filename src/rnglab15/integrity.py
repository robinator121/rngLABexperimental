"""Bounded raw capture and content-addressed, read-only evidence snapshots."""
from __future__ import annotations
import hashlib
import json
import os
from pathlib import Path
import stat
import tempfile

CHUNK = 1024 * 1024
MAX_LINE = 128 * 1024


def canonical(value):
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"


def write_json(path, value):
    with Path(path).open("x", encoding="utf-8") as f:
        f.write(canonical(value))
        f.flush()
        os.fsync(f.fileno())


def signature(s):
    return (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)


def digest(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError("only regular, non-symlink files are accepted")
    before = path.stat()
    h = hashlib.sha256()
    with path.open("rb") as f:
        if signature(os.fstat(f.fileno())) != signature(before):
            raise ValueError("input replaced before read")
        for block in iter(lambda: f.read(CHUNK), b""):
            h.update(block)
        if signature(os.fstat(f.fileno())) != signature(before):
            raise ValueError("input changed during read")
    if signature(path.stat()) != signature(before):
        raise ValueError("input replaced during read")
    return h.hexdigest()


def lines(path, max_line=MAX_LINE):
    """Reject oversized records before allocating an unbounded line."""
    with Path(path).open("rb") as f:
        while True:
            line = f.readline(max_line + 1)
            if not line:
                return
            if len(line) > max_line:
                raise ValueError("record exceeds configured byte limit")
            yield line.decode("utf-8", errors="strict")


def snapshot(path, vault):
    """Copy original bytes, verify stability, fsync, then publish without overwrite.

    Mode 0444 guards against accidental modification, not an adversarial owner.
    Hash verification is required on every subsequent use.
    """
    path, vault = Path(path), Path(vault)
    vault.mkdir(parents=True, exist_ok=True)
    original = digest(path)
    before = signature(path.stat())
    fd, tmp = tempfile.mkstemp(prefix=".pending-", dir=vault)
    try:
        h = hashlib.sha256()
        with os.fdopen(fd, "wb") as out, path.open("rb") as inp:
            for block in iter(lambda: inp.read(CHUNK), b""):
                h.update(block)
                out.write(block)
            out.flush()
            os.fsync(out.fileno())
        if h.hexdigest() != original or signature(path.stat()) != before:
            raise ValueError("source changed while snapshotting")
        target = vault / original
        os.chmod(tmp, 0o444)
        try:
            os.link(tmp, target)
        except FileExistsError:
            if digest(target) != original:
                raise ValueError("corrupt vault object")
        dirfd = os.open(vault, os.O_DIRECTORY)
        try:
            os.fsync(dirfd)
        finally:
            os.close(dirfd)
        return target, original
    finally:
        os.unlink(tmp)


def capture(chunks, target, max_bytes):
    """Finite transport-neutral acquisition. Interrupted captures remain .partial.

    Feed raw serial/socket read chunks here; do not decode, repair or drop bytes.
    A caller must supply a finite timeout for each transport read.
    """
    if not isinstance(max_bytes, int) or max_bytes <= 0:
        raise ValueError("positive finite byte limit required")
    target = Path(target)
    if target.exists():
        raise FileExistsError(target)
    partial = target.with_name(target.name + ".partial")
    h, total = hashlib.sha256(), 0
    with partial.open("xb") as out:
        for chunk in chunks:
            if not isinstance(chunk, bytes):
                raise TypeError("capture chunks must be bytes")
            if total + len(chunk) > max_bytes:
                raise ValueError("capture byte budget exceeded; partial preserved")
            out.write(chunk)
            h.update(chunk)
            total += len(chunk)
        out.flush()
        os.fsync(out.fileno())
    os.chmod(partial, 0o444)
    os.link(partial, target)
    partial.unlink()
    return {"bytes": total, "sha256": h.hexdigest(), "complete": True}


def inventory(root):
    """No content preview; explicit privacy review is needed before ingest."""
    root = Path(root).resolve(strict=True)
    result = []
    for p in sorted(root.rglob("*")):
        s = p.lstat()
        if p.is_dir() and not p.is_symlink():
            continue
        item = {"path": str(p.relative_to(root)), "size_bytes": s.st_size,
                "permissions": stat.filemode(s.st_mode), "mode": oct(stat.S_IMODE(s.st_mode)),
                "uid": s.st_uid, "gid": s.st_gid, "mtime_ns": s.st_mtime_ns,
                "symlink": p.is_symlink(), "format": "unclassified",
                "provenance": "unverified", "sample_rate_hz": None}
        if stat.S_ISREG(s.st_mode):
            item["sha256"] = digest(p)
            if p.suffix == ".csv":
                import csv
                header = next(csv.reader(lines(p)), [])
                item["columns"] = header
                if "time_us" in header and "ax" in header:
                    item["format"] = "sensor_csv"
                elif "timestamp" in header:
                    item["format"] = "event_csv_candidate"
                else:
                    item["format"] = "csv_unclassified"
            elif p.suffix == ".txt":
                item["format"] = "text_unclassified"
            elif p.suffix == ".py":
                item["format"] = "capture_program"
        result.append(item)
    return result
