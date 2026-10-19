# Detector export assessment for RNGLab

Analyzed 19 October 2026 · Full local file · Descriptive assessment, no TRNG validation

## What the data supports

The local desktop export is **606.6 MB**, with **117,962 records and 120.7 million waveform values**. It is far larger than the preliminary 4-record sample. It is useful for waveform analysis, instrument monitoring, and selecting candidate physical entropy measurements. It does **not** establish that the pulses are muons, that every pulse is radiation, or that the exported numbers are suitable random bits.

The strongest findings are an almost fixed trigger position, extensive ceiling clipping, many smaller pulses in the nominal baseline, pronounced pulse recovery, and event gaps that are neither exponentially distributed nor fully uncorrelated. No validated entropy per sample, entropy per event, or secure output bit rate is assigned.

## 1. Verified file and layout

| Check | Result |
|---|---|
| Input | `meas_data_2026-10-19 14.05.19.txt` on the desktop |
| Exact size | 606,632,259 bytes = 606.632259 MB = 578.529605 MiB |
| Nonempty data records | 117,962 |
| Fields on every row | 1,025 whitespace-separated decimal integers |
| Working interpretation | Field 0, field 1, then 1,023 waveform codes |
| Waveform sample indices | 0–1,022 (all indices in this report are zero based) |
| Total waveform codes | 120,675,126 |
| Observed code range | 0–65,535 |
| Line endings / malformed records | All CRLF / none found |
| Identical full waveforms | None, using 128-bit waveform digests; independent of timestamp |
| Source changed during scan | No size or modification-time change |

SHA-256 of the entire original file:

`feba450d95b074d3b0fff50df8357a682e90ea62614a5948e3e4e106149b4e10`

The earlier attachment is **20,550 bytes and four records**. Numeric comparison confirms that these are exactly the first four records of this desktop file; whitespace/line endings differ. Conclusions from only that attachment would have been very limited. Here, integrity counts, basic features, code histograms, timing and feature correlations use **all records**. Ensemble plots and secondary width measurements use clearly specified subsets.

The file has no header, instrument settings, column labels, or units. The two-field interpretation is supported by a constant first value, monotonic second value, and pulse-like sequence afterward. It remains an inferred schema. There are **1,023 stored waveform values**, not an assumed 1,024; the exporter manual/source is needed to establish whether that is intentional. A valid final row does not prove the acquisition was lossless or complete before export.

| Field | Observation | Interpretation and remaining uncertainty |
|---|---|---|
| 0 | Zero in all 117,962 rows | Could be channel, board, or status; no label can be recovered from a constant |
| 1 | Strictly increasing from 6,708,800 to 1,323,954,661,040; all multiples of 80 | Strong timestamp/counter candidate; epoch, units, clock and conversion remain unknown |
| 2 onward | Fast positive peaks, ceiling plateaus and recovery in codes 0–65,535 | Digitized waveform is strongly supported; unsigned 16-bit storage is consistent, but native ADC resolution and any rescaling are unknown |

## 2. Trigger alignment and saturation

Within samples 100–159, the first code above 10,000 is at **124 in 117,464 records (99.58%)** and at **123 in the remaining 498**. Every record crosses that threshold. The nominal main-pulse maximum is at least 10,748. This is strong evidence of acquisition or export alignment to a large pulse, rather than free-running, unselected waveform capture. The true hardware trigger threshold is not known: 10,000 is an analysis threshold chosen here.

**108,473 records (91.96%)** reach 65,535 somewhere. The nominal main-pulse window (110–159) reaches it in **108,393 records (91.89%)**. Across all waveforms, **627,833 samples (0.520%)** equal 65,535. The median record contains six such samples. These repeated flat maxima strongly suggest an export/ADC ceiling. Whether clipping occurs in the ADC, analog chain, or software needs documentation.

Peak height is therefore censored in most events. Neither the clipped maximum nor a simple integral of these waveforms can be treated as calibrated deposited energy. Clipping duration can track signal size but is not an independent measurement of amplitude or entropy. The first maximum index is also biased toward the start of a plateau.

![Waveform overview](01_waveforms.png)

![Alignment and clipping](02_alignment_clipping.png)

## 3. Baseline, small pulses and recovery

Using the median of samples 0–99 for each record, the overall median pre-pulse level is **3,341.25 codes**; 90% of record medians lie between **3,310.5 and 3,397**. This window is often pulse contaminated: the median raw standard deviation is **118.3 codes**, while the median robust adjacent-difference scale is **25.2 codes**. The latter is a noise-scale proxy, not a calibrated white-noise RMS.

Only **254 records (0.215%)** have a pre-window maximum less than 100 codes above its median. In that selected quiet subset, the median raw standard deviation is **19.4 codes** and the difference scale **17.3 codes**. This gives a more useful view of background fluctuations, but the selection itself truncates large excursions and must not be mistaken for an unbiased source characterization.

After linear detrending each pre-window, pooled lag-1 autocorrelation is **0.855** across all records and **0.146** in the quiet subset. Adjacent samples are therefore dependent; counting each ADC word or low bit as fresh independent entropy is unjustified. The quiet estimate comes from only 254 short, selected windows, not a long free-running noise run.

With a fixed peak-finding rule (prominence at least 150 codes, separation at least three samples), records contain on average **4.02 pre-window peaks** and **34.03 peaks in samples 220–1,022**. These are algorithmic local maxima, not certified photon or particle counts. In 6,000 evenly spaced records, 203,971 post-window peaks have median half-prominence width **4.40 samples** and median prominence **275 codes**. Isolated small peaks show a rapid rise and slower decay; this is waveform morphology, not an identified physical mechanism.

There is a large undershoot: the median per-record difference between samples 160–219 and the pre-window median is **−1,121.5 codes**. Even at samples 900–1,022 it is **−393.5 codes**. The post-pulse region has not usually returned to the pre-pulse level by the end of the stored waveform. Post-pulse “noise” includes this deterministic recovery and cannot be pooled with the pre-window without a model.

**Hypotheses to test:** the small peaks could include sensor dark pulses, photoelectron signals, correlated avalanches, pickup or other electronic effects. The undershoot could reflect AC coupling, shaping, overload recovery or digital processing. The export alone does not separate these explanations.

![Baseline diagnostics](03_baseline_noise.png)

![Pulse shapes](08_pulse_shapes.png)

## 4. Timing and stability

There are **117,961 positive gaps**, with no repeats or backward steps in field 1. Their greatest common divisor is **80 raw units**, also the common divisor of all timestamp values. This is an observed lattice, not proof that the hardware clock resolution is exactly 80 units. All four lowest bits of the stored timestamps are consequently zero; extracting timestamp LSBs directly would fail immediately.

| Gap measurement | Raw field-1 units |
|---|---:|
| Minimum | 1,674,720 |
| Median | 9,276,800 |
| Mean | 11,223,607.4 |
| 95th percentile | 24,821,440 |
| Maximum | 79,862,880 |
| Total first-to-last span | 1,323,947,952,240 |

If the raw unit is a nanosecond, the span would be **22 min 3.948 s**, the mean gap **11.224 ms**, and the interval-based event rate **89.10/s**. This interpretation is consistent with the filename time (14:05:19) and file modification time (14:27:23.891, Budapest), but it is **a plausible unit inference, not a calibration**. No waveform time scale can be inferred from that agreement.

The gap coefficient of variation is **0.615**, substantially below the value 1 of an exponential waiting-time model. The maximum empirical-CDF difference from an exponential with the fitted mean is **0.251**, a large shape difference. This is a descriptive comparison; a standard fixed-parameter KS p-value would not apply to the fitted reference. Consecutive gaps have Pearson autocorrelation **−0.0432**; most later lags are much smaller. For orientation only, an IID reference band would be approximately ±0.0057. This is evidence against treating the recorded gap sequence as independent exponential arrivals.

Across 50 equal timestamp-span bins, counts range from **2,248 to 2,427**, with count variance/mean **0.555**. Rate appears broadly stable over this run, while arrivals are more regular than the simplest unfiltered Poisson model predicts. In 100 equal-record blocks, pre-baseline medians range only **3,337–3,345 codes**. This single acquisition cannot establish stability across restarts, temperature changes or different days.

**Hypotheses to test:** trigger holdoff/dead time, filtering, coincidence selection, acquisition buffering, event merging, clock/export processing, or a source with structured timing. A minimum observed gap alone does not measure dead time. The present data cannot identify which mechanism causes the non-exponential distribution.

![Timing diagnostics](04_timing.png)

## 5. Correlations and unusual records

Spearman rank correlations across all events are strong between ceiling occupancy and the main-window area (**+0.945**), ceiling occupancy and early post level (**−0.952**), and early and late post levels (**+0.973**). The main-window area and early post level are correlated **−0.984**. These features mostly describe different aspects of the same pulse and recovery, so their individual marginal variations cannot be added as independent entropy contributions. The ceiling/early-post relationship also persists in consecutive feature differences (rank correlation **−0.934**), so it is not explained solely by a slow run-wide drift.

Pre-baseline and pre-window peak count correlate **+0.619**; the difference-scale proxy and pre-window peak count correlate **+0.665**. This supports pulse contamination as a contributor to apparent “noise.” The preceding gap versus pre-baseline correlation is small (**+0.0143**); that does not rule out nonlinear or conditional dependence. None of these correlations establishes causality.

**1,249 records (1.06%)** contain more than one contiguous excursion above 10,000; five contain three. **1,198** have a new excursion starting at sample 180 or later. **28 records** contain a total of **3,490 zero codes**; one has 438 zeros. Those floor hits and extra pulses deserve inspection before selecting entropy windows. They could reflect overload, pileup or processing, but are not assigned a physical label here.

The example figure includes event 60,542 with multiple large pulses, event 79,942 with a very low pre-baseline and an early large pulse, and event 77,849 with extensive floor clipping. These are zero-based indices. The file `flagged_events.csv` lists **1,283 records** flagged by a descriptive rule: multiple >10,000 excursions, any zero code, or pre-baseline below 3,200. This rule is for review, not an established pass/fail health test.

![Stability](05_stability.png)

![Correlations](06_correlations.png)

![Representative records](07_examples.png)

## 6. RNGLab assessment

| Candidate measurement | Promise | Main obstacles and next step |
|---|---|---|
| Raw analog noise or ADC innovations in genuinely quiet windows | Highest priority for a controlled noise-source study | Acquire free-running untriggered noise; establish bandwidth, coupling, ADC behavior and predictors. Quiet-window selection here is biased and sparse |
| Small-pulse arrival positions and spacings | Worth investigating as an event-based physical source | Identify the sensor/process; measure afterpulsing, crosstalk, pileup, threshold effects and electronic interference; retain timestamps and amplitudes |
| Main event gaps after documented clock quantization | Worth a separate study | Current distribution and lag-1 dependence require explanation; determine holdoff, lost events, clock scale and whether triggers identify the intended physical events |
| Unclipped pulse amplitude, charge or width | Secondary candidate after instrument correction | Reduce gain or increase dynamic range, calibrate baseline and charge, and quantify dependence on timing and other pulse features |
| Clipped maxima, fixed trigger index, constant first field, raw timestamp low bits | Poor direct choices in this export | Clipped maxima and constant fields are often fixed; trigger position is almost fixed; low timestamp bits are identically zero |

For example, bit 0 is one in **49.4258%** of the 11,796,200 pre-window ADC codes. A near-balanced bit histogram does not measure unpredictability, and the observed bias and sample dependence must be modeled. No Shannon-entropy histogram, compression result or randomness test has been substituted for a conditional min-entropy assessment.

An RNGLab pipeline should preserve the raw sequence and settings, define one candidate symbol stream in acquisition order, monitor it before conditioning, and keep a separate record of rejected events. Any adaptive window or rejection rule changes the source being assessed. Dividing timestamps by the observed factor 80 removes a storage factor; it does not make their remaining bits independent. Combining correlated pulse features or hashing the file does not create additional physical entropy.

NIST SP 800-90B frames evaluation around a defined noise source, entropy assessment, conditioning and health tests. Its collection requirements include at least one million raw sequential sample values and restart data from 1,000 restarts with 1,000 samples each. The appropriate source symbol must be defined first: 120.7 million waveform codes are not automatically that many independent source outputs, while a one-symbol-per-event stream here has only about 118,000 symbols. Use the non-IID route unless independence is justified, and assess worst-case operating conditions. See [SP 800-90B, especially Sections 3.1 and 4](https://nvlpubs.nist.gov/nistpubs/specialpublications/nist.sp.800-90b.pdf); the [official publication page](https://csrc.nist.gov/pubs/sp/800/90/b/final) also lists potential errata.

The [official NIST assessment software](https://github.com/usnistgov/SP800-90B_EntropyAssessment) provides IID, non-IID and restart tools. Its input symbols fit in one byte, so conversion from these 16-bit codes needs an explicit, documented representation. No NIST validation suite was run here, and no certification claim follows from this report.

For a one-time export, potential physical entropy belongs to the original acquisition. A fixed file is replayable; someone who has a copy can reproduce any deterministic output derived from it. It cannot serve as an ongoing supply of fresh entropy without fresh measurements.

### Monitoring signals to implement first

| Signal | What to monitor |
|---|---|
| Acquisition integrity | Record length, parse errors, timestamp reversals/repeats, counter wraps, dropped-event indicators and identical repeated waveforms |
| ADC and analog chain | Ceiling/floor fraction, longest rail run, stuck codes/bits, missing codes, baseline and noise-scale changes |
| Trigger and pulse shape | Crossing-index distribution, main-pulse amplitude/width, secondary/large-pulse count, undershoot and end-of-window recovery |
| Timing | Gap lattice/residue patterns, minimum and quantiles, rate in fixed time bins, serial correlation and bursts |
| Source health | Raw repetition-count and adaptive-proportion tests after source definition and entropy assessment; thresholds must follow that assessment |
| Dependence and environment | Conditional feature relationships, periodic pickup, temperature/bias changes, restart transients and correlations with external controls |

The numerical ranges measured in this run are a starting baseline for comparison, not validated alarm limits. In particular, high clipping is normal in this file but may signal unsuitable gain for future measurements.

## 7. Metadata and acquisitions needed

1. **Hardware identity and physical arrangement:** detector/sensor, scintillator if used, optical coupling, shielding, geometry, number of channels, coincidence logic and independent evidence that main triggers correspond to radiation or muons.
2. **Digitizer details:** model, ADC native bits, signed/unsigned encoding, voltage range, gain/offset, sample rate and clock, analog bandwidth, filtering/decimation/interpolation, any code scaling, compression, pedestal subtraction or clipping.
3. **Trigger settings:** source/channel, threshold and polarity, pretrigger length, holdoff, dead time, rearm behavior, pileup rejection, prescale, coincidence/veto windows, software cuts and any artificial pulser.
4. **Timestamp specification:** column meaning, unit, time origin, counter width/wrap, relation to waveform sampling, whether stored values are event times or host arrival times, and reason for the factor 80.
5. **Exporter specification:** application/version, exact column schema, why 1,023 values are stored, whether samples or records are omitted, acquisition start/stop and live time, buffer overflows, dropped events and ordering.
6. **Operating conditions:** bias/high voltage, temperature, power supply, light exposure, grounding, environmental logs and any changes during the run.
7. **Controlled follow-up data:** untriggered/free-running baseline, sensor disconnected or blocked where appropriate, light/source controls, lower-gain unclipped pulses, longer recovery windows, longer runs and independent restart datasets. Keep raw settings and event loss counters with each run.

The radiation claim remains open. Electronic and sensor noise may themselves be useful entropy mechanisms, but their presence does not establish that the export measures radiation, or specifically muons.

## 8. How this full export was inspected, and how to inspect another

The earlier four-row attachment was an excerpt; the desktop file used here was scanned in full. The workflow streams the text rather than opening a 606 MB document in an editor. It checks byte size, integer syntax, row widths, source hash, waveform range, timestamps and duplicate waveform digests. It then stores waveform codes in a temporary 16-bit array (about 241 MB) and processes features in batches. Plots use aggregates and declared subsets while integrity statistics cover the complete file.

For another nominal “500 MB” export, point the supplied `inspect_export.py` at the original file. It reports actual bytes, all row widths, leading-field distributions, waveform-code range, timestamp changes and SHA-256 without loading the whole text into RAM. It assumes the same two-leading-field layout; verify that assumption against the exporter. It does not run entropy validation. Compare reported record count and acquisition live time against the instrument log, and transfer the actual original file rather than a preliminary excerpt or pasted rows.

The `analysis_source.zip` archive contains the exact analysis workflow and a reproduction README. `event_features.csv` has one row per record and retains field 1 and the preceding raw gap. `file_manifest.json`, `analysis_metrics.json` and `supplement_metrics.json` retain machine-readable measurements, including full distributions and definitions encoded in the scripts. `flagged_events.csv` supports targeted manual review. All charts are also separate PNGs.

## Measurement definitions and limitations

- Pre-baseline: median of indices 0–99. Early and late post levels: medians of 160–219 and 900–1,022. Window differences are computed per record before taking their median.
- Raw pre standard deviation uses sample denominator n−1. MAD scale is median absolute deviation / 0.67448975. Adjacent-difference scale is MAD of first differences / (0.67448975 × √2); the white-Gaussian interpretation requires assumptions not established here.
- Main area: sum of codes minus the pre-baseline at indices 110–179; units are code×sample. It is clipped, window dependent and uncalibrated. Main maximum and clipping window: 110–159. Total `clip_samples` counts the entire waveform, including additional pulses.
- Threshold crossing: first code >10,000 in indices 100–159. Rise index: largest adjacent positive difference in that window, labeled by its left sample. Global peak index: first occurrence of the maximum across the entire record.
- Small-peak counts: SciPy local maxima, minimum prominence 150 and minimum separation three samples. Pre search 0–99, post search 220–1,022; endpoint peaks are not counted. Width is at half prominence, not necessarily physical FWHM. Threshold sensitivity was not exhaustively assessed.
- Small-peak widths/prominences use 6,000 evenly spaced records. Small-peak templates use up to 10,000 qualifying snippets from those records, prominence strictly between 180 and 300, with no other detected peak within 25 samples. They are selected examples, not a calibrated single-photoelectron template.
- Pre-window autocorrelation pools products after removing a separate fitted line from each record, normalized by total residual energy. It never joins waveform boundaries. Event autocorrelations use a global sample mean and energy denominator. Correlations are descriptive; broad testing and multiple comparisons were not used to declare validation.
- No detector manual, sampling-rate calibration, physical source control, event-loss log, environmental record, restart collection, formal min-entropy bound, or independent radiation tag was available.

Reproduction environment: Python 3.12.14, NumPy 2.5.4, SciPy 1.18.1, Matplotlib 3.11.2. Source remained read only; all generated material is separate from it.
