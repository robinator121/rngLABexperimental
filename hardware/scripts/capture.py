
import argparse
import csv
import re
import time
from datetime import datetime, timezone
from pathlib import Path

import serial

parser = argparse.ArgumentParser()
parser.add_argument("experiment")
parser.add_argument("--seconds", type=int, default=300)
args = parser.parse_args()

if not re.fullmatch(r"[A-Za-z0-9_-]+", args.experiment):
    parser.error("Use letters, digits, underscores or hyphens")

if args.seconds <= 0:
    parser.error("--seconds must be positive")

root = Path.home() / "rnglab-data"
root.mkdir(parents=True, exist_ok=True)

session = datetime.now(timezone.utc).strftime(
    "%Y%m%dT%H%M%SZ"
)

filename = root / f"{args.experiment}_{session}.csv"

sensor_columns = [
    "time_us",
    "ax", "ay", "az",
    "gx", "gy", "gz",
    "pressure_raw",
    "temp_raw",
    "humidity_raw",
    "ir_pulses",
    "ir_age_us",
]

columns = [
    "experiment",
    "sample_id",
    "host_timestamp_utc",
    *sensor_columns,
]

count = 0

print(f"Experiment: {args.experiment}")
print(f"Duration: {args.seconds} seconds")
print(f"Output: {filename}")

try:
    with serial.Serial(
        "/dev/ttyACM0", 115200, timeout=1
    ) as device:

        # Opening the port resets the Arduino Uno.
        time.sleep(2)
        start = time.monotonic()

        with filename.open(
            "w", newline=""
        ) as file:

            writer = csv.writer(file)
            writer.writerow(columns)

            try:
                while time.monotonic() - start < args.seconds:
                    raw = device.readline().decode(
                        "utf-8", errors="replace"
                    ).strip()

                    fields = raw.split(",")

                    if len(fields) != len(sensor_columns):
                        continue

                    if not fields[0].isdigit():
                        continue

                    received = datetime.now(
                        timezone.utc
                    ).isoformat()

                    writer.writerow([
                        args.experiment,
                        count,
                        received,
                        *fields,
                    ])

                    count += 1

                    if count % 1000 == 0:
                        file.flush()
                        print(f"Captured {count:,} rows")

            except KeyboardInterrupt:
                print("\nStopped manually.")

except serial.SerialException as exc:
    print(f"Serial error: {exc}")

print(f"Saved {count:,} rows to {filename}")

