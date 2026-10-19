# GM dataset roadmap: week of 26 October 2026

1. Before acquisition, obtain instrument/tube identity, documented output type,
   isolated low-voltage interface, timestamp units and counter width, measured
   clock precision, trigger/threshold, holdoff, dead-time model, missed/overflow
   counters, acquisition firmware/export versions and environmental settings.
   SICHAS/iSES event output availability is unconfirmed. Count-only readouts stay
   count windows with start/end/live-time metadata; never fabricate individual
   arrival times. Technician-approved isolated interfaces only for HV equipment.
2. Save untouched raw transport/files alongside settings and provenance. Record
   startup/reset boundaries, sequence IDs, host time and synchronization uncertainty.
   Keep detector-off/electronic controls and background runs distinct. Radiation
   identity needs detector evidence; waveform morphology alone does not identify it.
3. Start with event timestamps and gap distributions. Compare clock lattices,
   repeats/reversals, bursts, rate drift and serial dependence. Fit documented
   nonparalyzable/paralyzable dead-time models with censored observations and
   independent controls. A minimum gap alone is not a dead-time calibration.
   The inherited v0.14 GM emulator treats its two modes similarly; do not use it
   as a validated dead-time reference. New defect work here tests supplied bytes,
   not a claimed detector physics model.
4. Register the source symbol before looking at confirmatory data. Collect at least
   1M sequential outputs of that defined source and 1000 x 1000 restart samples
   for the SP800-90B assessment plan; verify current standard/errata details with
   the assessor. Plan acquisition duration from observed event rate, rather than
   treating every timestamp bit or waveform code as independent. At 100 events/s,
   1M events alone needs about 2.78 hours; lower-rate GM backgrounds need longer.
5. Split by acquisition/restart/day/hardware conditions. Select windows, symbols,
   models, predictors, alarms and extraction budget on development runs. Freeze
   them and run official non-IID/restart tools on fresh held-out acquisitions.
   Assess worst-case conditions and known interference/predictors. IID needs
   specific supporting evidence. Do not sum correlated channel entropy budgets.
6. Calibrate continuous RCT/APT to the validated conservative symbol entropy and
   chosen false-alarm budget. Inject stuck codes, duplicate/lost pulses, clock
   quantization, holdoff, bursts, periodic interference and drift. Test latching,
   quarantine, restart behavior, and the absence of emitted credited output on
   failure. Tests on conditioned output cannot replace raw-source health tests.
7. Only after external entropy review, select an extractor with independent seed
   provenance and per-block entropy/security accounting. Keep direct physical
   extraction and an explicitly seeded DRBG as distinct experimental arms.
   Measure sustainable credited bit rate including startup, rejection and outages.
8. Use fresh independent runs for physical dirty/conditioned application confirmation.
   Increase finite input budgets and run horizons after pilot power/convergence
   studies; freeze pi/VMC spatial moments/Ising/PIMC endpoints and corrections.
   Publish only approved non-sensitive physical-source aggregates and source code.

Concrete current blockers: no GM files, no verified event protocol, no clock/trigger
calibration, no restart/condition coverage, no validated source entropy, no official
NIST executable installed, and no ESP32 RNG capture with entropy-source status.
