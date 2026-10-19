"""Auditable adapter to external official NIST executables, never a substitute."""
from __future__ import annotations
from pathlib import Path
import shutil
import subprocess
from .integrity import digest,write_json


def assess(symbol_file, out, *, mode="non_iid", executable=None, bits_per_symbol=1, timeout=300):
    if mode not in ("iid","non_iid") or not 1<=bits_per_symbol<=8 or not 0<timeout<=3600:
        raise ValueError("invalid NIST invocation")
    source=Path(symbol_file)
    source_hash=digest(source)
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    metadata={"mode":mode,"symbol_sha256":source_hash,"symbols":source.stat().st_size,
              "bits_per_symbol":bits_per_symbol,"validated_entropy_credit_bits":0,
              "qualification":"external output requires source-model and restart review; no automatic entropy approval"}
    try:
        if source.stat().st_size<1_000_000:
            metadata.update(status="blocked",reason="fewer than one million sequential raw symbols")
            return metadata
        found=shutil.which(executable or f"ea_{mode}")
        if not found:
            metadata.update(status="blocked",reason=f"official ea_{mode} executable unavailable")
            return metadata
        maximum=(1<<bits_per_symbol)-1
        with source.open("rb") as f:
            for block in iter(lambda:f.read(1024*1024),b""):
                if any(x>maximum for x in block):
                    raise ValueError("symbol value exceeds declared bit width")
        command=[found,"-i","-a",str(source.resolve()),str(bits_per_symbol)]
        metadata.update(command=command,executable_sha256=digest(found))
        with (out/"stdout.txt").open("x") as stdout,(out/"stderr.txt").open("x") as stderr:
            try:
                result=subprocess.run(command,stdout=stdout,stderr=stderr,timeout=timeout,check=False)
                metadata.update(status="executed" if result.returncode==0 else "tool_failed",returncode=result.returncode)
            except subprocess.TimeoutExpired:
                metadata.update(status="timed_out",timeout_seconds=timeout)
        if digest(source)!=source_hash:
            raise ValueError("NIST input changed")
        return metadata
    except Exception as exc:
        metadata.update(status="failed",reason=f"{type(exc).__name__}: {exc}")
        raise
    finally:
        write_json(out/"assessment.json",metadata)
