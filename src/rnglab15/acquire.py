"""Optional bounded serial capture, preserving startup messages and bad records."""
from __future__ import annotations
from pathlib import Path
import math
import time
from .integrity import capture,write_json,digest


def serial_capture(port, target, *, seconds, max_bytes, baud=115200, provenance):
    if not math.isfinite(seconds) or seconds<=0 or max_bytes<=0 or baud<=0 or not provenance:
        raise ValueError("finite capture bounds and hardware provenance required")
    import serial  # optional pyserial, no port discovery or opening on import
    target=Path(target)
    # Opening UNO ports may reset the MCU; preserve startup and time discontinuities.
    started=time.time_ns()
    with serial.Serial(port,baud,timeout=0.25) as device:
        deadline=time.monotonic()+seconds
        def chunks():
            total=0
            while time.monotonic()<deadline and total<max_bytes:
                block=device.read(min(4096,max_bytes-total))
                if block:
                    total+=len(block)
                    yield block
        result=capture(chunks(),target,max_bytes)
    result.update(started_wall_ns=started,ended_wall_ns=time.time_ns(),port=port,baud=baud,
                  provenance=provenance,max_bytes=max_bytes,seconds_limit=seconds,
                  transport_note="host receive times are not hardware event timestamps; no dropped-byte guarantee")
    write_json(target.with_name(target.name+'.manifest.json'),result)
    return result
