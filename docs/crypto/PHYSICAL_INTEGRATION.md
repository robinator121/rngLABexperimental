# Physical entropy integration audit

The historical [`MUON_DETECTOR_ASSESSMENT.md`](../physical/MUON_DETECTOR_ASSESSMENT.md) documents a 117,962-record, 1,023-sample-per-waveform export and reports extensive digitizer clipping and an 80-unit timestamp lattice. Those numbers are **historical reported measurements**, not reproduced by BLUE: the original full raw export is absent from this repository and the current project inventory. No detector identity, timestamp-unit calibration or original acquisition metadata has been independently verified in this task. The device must not be described as a validated muon TRNG.

The report's declared acquisition interval is 2026-10-19 14:05:19–14:27:23 CEST. Host/acquisition timestamp labels do not establish the physical clock's units. Timing GCD, counter reversals, sample order, clipping and dependencies must be checked on original bytes before choosing a source symbol. Nominal 16-bit storage does not establish 16 effective bits or 16 bits of entropy. Triggered waveform shape, pulse saturation, clock quantization, conditional selection of quiet windows and correlated samples can severely reduce useful uncertainty.

BLUE supplies a functioning **read-only** audit command:

```bash
python -m rnglab_crypto physical-audit /path/to/meas_data.txt \
  --metadata /path/to/acquisition.json --output /new/path/physical-audit.json
```

Minimal sidecar:

```json
{
  "samples_per_record": 1023,
  "provenance": "original acquisition identifier and instrument/export references",
  "acquisition_start": null,
  "acquisition_end": null,
  "native_adc_bits": null
}
```

Optional `seconds_per_tick` is accepted only with a nonempty `clock_evidence` reference. The command preserves record order, rejects malformed dimensions/out-of-range codes, computes a full-file SHA256, timestamp reversals/repeats/GCD, rail fractions, distinct codes and within-record raw lag-one correlation. It never joins the end of one waveform to the next. Correlation includes pulse shape, is descriptive and is not a noise-innovation estimate. Effective resolution remains null without calibration. It assigns zero validated cryptographic entropy credit.

`CapturedSource` also provides finite read-only, content-hashed replay with explicit physical/provenance labels. Exhaustion fails rather than looping. Capturing a real sample makes it available for statistical diagnostics or simulated downstream treatments; it does not make a public fixed file a fresh secret source. The default crypto campaigns use no physical data and generate no real keys from recordings.

Required next measurements before an entropy-security claim:

1. Instrument and sensor identification; ADC encoding, true resolution, voltage calibration, bandwidth, sample clock, preprocessing and trigger settings.
2. Raw, sequential, free-running noise with control runs (blocked/disconnected sensor where meaningful), lower-gain unclipped acquisition, timestamps and dropped-event counters.
3. Physical separation and environmental logs for temperature, bias, supply, coupling, periodic pickup, dead time and restart transients. Evidence that any claimed radiation/muon events have the stated physical origin.
4. A prespecified source symbol and conditioning pipeline; conditional predictors and dependence evaluation; worst-case operating regimes and attacker side-information model.
5. SP800-90B-appropriate sequential/restart collections and qualified assessment, using independently reviewed tools. A million correlated ADC codes are not automatically a million independent entropy outputs.
6. Assessed startup/reseed budgets and calibrated online health tests on raw symbols, with fail-closed handling and independent lifecycle review.

Physical integration is therefore **implemented as an audit/replay boundary and tested on synthetic fixtures**, while the real-data campaign and entropy qualification remain blocked by unavailable original data and measurement evidence. This limitation does not block the synthetic crypto campaign.
