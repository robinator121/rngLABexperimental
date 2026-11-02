"""Paired CPU-only research campaign. Inputs and outcomes are synthetic."""
from __future__ import annotations

from collections import Counter
import hashlib
from pathlib import Path
import tempfile
import time

import numpy as np
from scipy.stats import binomtest
from cryptography.exceptions import InvalidTag
from cryptography.hazmat.primitives.asymmetric.ed25519 import Ed25519PrivateKey

from rnglab11.drbg import HMACDRBG
from rnglab_crypto.sources import ReplaySource, seed_for
from rnglab_crypto.faults import FaultAdapter, FaultSpec, EntropyUnavailable
from rnglab_crypto.monitoring import (calibrate, diagnostics, source_scores, health_screen,
                                      EntropyGate, SeedEvidence)
from rnglab_crypto.primitives.conditioning import condition, expand
from rnglab_crypto.primitives.aead import (ALGORITHMS, NonceLedger, GuardedAEAD, NonceReuse,
                                         cipher, _xor, _gcm_one_block_candidate, nonce_misuse_demo)
from rnglab_crypto.primitives.signatures import signature_demo
from rnglab_crypto.risk import event_rate, birthday_probability, survival_summary, paired_difference
from rnglab14.statistics.multitest import holm_bonferroni_correction


def _source(config, replicate, role):
    return ReplaySource(seed_for(config.namespace, replicate, role), max_bytes=1_048_576)


def weak_keys(config):
    rows = []
    for rep in range(config.repetitions):
        base = _source(config, rep, "weak-key-selection")
        healthy = _source(config, rep, "key-reference")
        control = [healthy.read_bytes(32) for _ in range(64)]
        for bits in config.weak_seed_bits:
            values = [int.from_bytes(base.read_bytes(4), "big") & ((1 << bits)-1) for _ in range(64)]
            keys = [condition(x.to_bytes(32, "big")) for x in values]
            chosen = values[0]
            nonce, plaintext = bytes(12), b"BLUE synthetic known plaintext"
            ct = cipher("AES-256-GCM", keys[0]).encrypt(nonce, plaintext, b"BLUE")
            recovered, guesses = False, 0
            # Deliberately bounded dictionary; no external ciphertext input.
            for candidate in range(1 << bits):
                guesses += 1
                key = condition(candidate.to_bytes(32, "big"))
                try:
                    recovered = cipher("AES-256-GCM", key).decrypt(nonce, ct, b"BLUE") == plaintext
                except InvalidTag:
                    continue
                if recovered:
                    break
            # The same 32-byte values are Ed25519 seeds. Equal public keys expose
            # seed duplication even though deterministic signing itself is sound.
            public = [Ed25519PrivateKey.from_private_bytes(k).public_key().public_bytes_raw() for k in keys]
            exposed = bits//2
            unknown = bits-exposed
            prefix = chosen >> unknown
            partial_guesses = 0
            for suffix in range(1 << unknown):
                partial_guesses += 1
                candidate = condition(((prefix << unknown)+suffix).to_bytes(32, "big"))
                if candidate == keys[0]:
                    break
            reused_key = condition(chosen.to_bytes(32, "big"))
            gate = EntropyGate()
            try:
                gate.check(SeedEvidence(bits, "synthetic_model"))
                blocked = False
            except RuntimeError:
                blocked = True
            rows.append({"replicate": rep, "seed_bits": bits, "keys_generated": 64,
                         "duplicates": 64-len(set(keys)), "ed25519_duplicate_public_keys": 64-len(set(public)),
                         "healthy_duplicates": 64-len(set(control)),
                         "recovered": recovered, "candidate_attempts": guesses,
                         "partial_seed_exposure_bits": exposed, "partial_exposure_attempts": partial_guesses,
                         "partial_exposure_recovered": candidate == keys[0],
                         "reinitialization_duplicates_key": reused_key == keys[0],
                         "startup_policy_blocked": blocked,
                         "entropy_semantics": "uniform bounded seed support, conditional on attacker not knowing replay seed"})
    return rows


def conditioning_blind_spot(config):
    rows = []
    for rep in range(config.repetitions):
        raw = _source(config, rep, "conditioning-input").read_bytes(32)
        weak = bytes(31)+raw[-1:]
        for conditioner in ("sha256", "hkdf"):
            for mechanism in ("hmac", "chacha20"):
                conditioned = condition(weak, conditioner)
                data = expand(conditioned, config.diagnostic_bytes, mechanism)
                ref = expand(condition(raw, conditioner), config.diagnostic_bytes, mechanism)
                prediction = expand(condition(weak, conditioner), config.diagnostic_bytes, mechanism)
                rows.append({"replicate": rep, "conditioner": conditioner, "expander": mechanism,
                             "weak_output": diagnostics(data), "reference_output": diagnostics(ref),
                             "raw_health_alarm": health_screen(weak)["latched_failure"],
                             "attacker_exact_prediction": prediction == data,
                             "output_sha256": hashlib.sha256(data).hexdigest(),
                             "seed_support_bits": 8, "conditional_bits_if_seed_known": 0,
                             "same_bytes_with_or_without_attacker_knowledge": True})
    return rows


def recovery(config):
    rows = []
    for rep in range(config.repetitions):
        seed = seed_for(config.namespace, rep, "recovery-initial")
        victim = HMACDRBG(seed, bytes(16))
        state = victim.export_sensitive_state(allow_sensitive_research_state=True)
        attacker = HMACDRBG.from_sensitive_state(state, allow_sensitive_research_state=True)
        old_key, guess = victim.read(32), attacker.read(32)
        old_cipher = cipher("AES-256-GCM", old_key)
        old_ct = old_cipher.encrypt(bytes(12), b"BLUE past message", None)
        victim.reseed(bytes(32))
        attacker.reseed(bytes(32))
        known_match = victim.read(64) == attacker.read(64)
        restored = seed_for(config.namespace, rep, "restored-model-secret")
        victim.reseed(restored)
        attacker.reseed(bytes(32))
        future, future_guess = victim.read(32), attacker.read(32)
        future_ct = cipher("AES-256-GCM", future).encrypt(bytes(12), b"BLUE rekeyed message", None)
        try:
            cipher("AES-256-GCM", future_guess).decrypt(bytes(12), future_ct, None)
            future_rejected = False
        except InvalidTag:
            future_rejected = True
        # Healthy retained state + known low-entropy reseed is NOT state replacement.
        unknown_state = HMACDRBG(seed_for(config.namespace, rep, "retained-secret"), bytes(16))
        guessed_state = HMACDRBG(bytes(32), bytes(16))
        unknown_state.reseed(bytes(32))
        guessed_state.reseed(bytes(32))
        retained_resists = unknown_state.read(32) != guessed_state.read(32)
        gate = EntropyGate()
        blocks = []
        for evidence in (SeedEvidence(None, "unvalidated"), SeedEvidence(0, "synthetic_model")):
            try:
                gate.check(evidence)
                blocks.append(False)
            except RuntimeError:
                blocks.append(True)
        gate.check(SeedEvidence(256, "synthetic_model"), recovery=True)
        starved = FaultAdapter(_source(config, rep, "starvation"), FaultSpec("starvation"), lab_only=True)
        try:
            starved.read_bytes(32)
            starvation_blocked = False
        except EntropyUnavailable:
            starvation_blocked = True
        rows.append({"replicate": rep, "known_state_prediction": guess == old_key,
                     "known_reseed_still_predictable": known_match,
                     "withheld_reseed_breaks_old_clone_prediction": future != future_guess,
                     "old_state_clone_rejected_after_rekey": future_rejected,
                     "past_ciphertext_still_decryptable": cipher("AES-256-GCM", guess).decrypt(bytes(12), old_ct, None) == b"BLUE past message",
                     "retained_unknown_state_resists_known_reseed_clone": retained_resists,
                     "unvalidated_startup_and_reseed_blocked": all(blocks),
                     "restored_model_budget_releases_quarantine": not gate.quarantined,
                     "starvation_blocks_output": starvation_blocked,
                     "interpretation": "clone mismatch is measured; entropy/secrecy of restored model input is assumed, not proven"})
    return rows


def _nonce_path(config, rep, name, directory):
    base = _source(config, rep, "nonce-source")
    kinds = {"serial": ("serial", .95),
             "drift": ("drift", 1.), "rollback": ("rollback", 1.)}
    if name.startswith("low_"):
        fault = FaultSpec("reduced_entropy", onset=config.fault_onset,
                          duration=config.operations-config.fault_onset,
                          entropy_bits=int(name[4:]), exposed_bits=0, target_layer="raw_source")
        source = FaultAdapter(base, fault, lab_only=True)
    elif name in kinds:
        kind, strength = kinds[name]
        source = FaultAdapter(base, FaultSpec(kind, onset=config.fault_onset,
                              duration=config.operations-config.fault_onset, intensity=strength), lab_only=True)
    else:
        source = base
    raw, nonces = [], []
    ledger = None
    if name == "counter_persistent":
        key = seed_for(config.namespace, rep, "ledger-key")
        path = Path(directory)/f"nonce-{rep}.sqlite"
        ledger = NonceLedger(path, key)
    try:
        for op in range(config.operations):
            block = source.read_bytes(32)
            raw.append(block)
            if ledger:
                if op == config.fault_onset:
                    ledger.close()
                    ledger = NonceLedger(path, key)
                nonce = ledger.allocate()
            elif name == "fixed" and op >= config.fault_onset:
                nonce = bytes(12)  # downstream fault; the raw source stays healthy
            elif name == "counter_reset":
                nonce = (op if op < config.fault_onset else op-config.fault_onset).to_bytes(12, "big")
            else:
                nonce = hashlib.sha256(b"BLUE nonce"+block).digest()[:12]
            nonces.append(nonce)
    finally:
        if ledger:
            ledger.close()
    return raw, nonces, source.stream_digest


def _aead_path(config, rep, name, algorithm, nonces, alarm, source_digest):
    key = seed_for(config.namespace, rep, "aead-key")
    if algorithm == "AES-128-GCM":
        key = key[:16]
    aead = cipher(algorithm, key)
    source_gated_aead = cipher(algorithm, key)
    guarded = GuardedAEAD(algorithm, key, limit=config.operations)
    emitted_nonces = []
    source_emitted_nonces = []
    seen, first, compromised, forged = {}, None, None, None
    repeats, blocked, guarded_emissions = 0, 0, 0
    for op, nonce in enumerate(nonces):
        plaintext = op.to_bytes(16, "big")
        encrypted = aead.encrypt(nonce, plaintext, None)
        if alarm is None or op < alarm:
            source_gated_aead.encrypt(nonce, plaintext, None)
            source_emitted_nonces.append(nonce)
        try:
            guarded.encrypt(nonce, plaintext)
            emitted_nonces.append(nonce)
            guarded_emissions += 1
        except NonceReuse:
            blocked += 1
        if nonce in seen:
            repeats += 1
            if first is None:
                first = op
                old_plain, old_cipher = seen[nonce]
                predicted_plain = _xor(_xor(old_cipher[:16], encrypted[:16]), old_plain)
                if predicted_plain == plaintext:
                    compromised = op
                if algorithm.startswith("AES"):
                    candidate = _gcm_one_block_candidate(old_cipher, encrypted)
                    try:
                        recovered = aead.decrypt(nonce, candidate, None)
                        if recovered != old_plain:
                            forged = op
                    except InvalidTag:
                        pass
        else:
            seen[nonce] = (plaintext, encrypted)
    counts = Counter(nonces)
    onset = None if name in ("healthy", "counter_persistent") else config.fault_onset
    raw_fault = name.startswith("low_") or name in ("serial", "drift", "rollback")
    raw_healthy_alarm = alarm is not None and (not raw_fault or alarm < config.fault_onset)
    source_prevented = first is not None and alarm is not None and alarm <= first
    return {"replicate": rep, "condition": name, "algorithm": algorithm,
            "fault_onset": onset,
            "fault_layer": ("raw_source" if name.startswith("low_") or name in ("serial", "drift") else
                            "drbg_state" if name == "rollback" else
                            "downstream_nonce_api" if name in ("fixed", "counter_reset") else "none"),
            "raw_source_fault_present": raw_fault,
            "nonce_risk_change_onset": onset,
            "nonce_uniqueness_breach_onset": first,
            "insecurity_onset_qualification": "predictable nonces alone are allowed; degradation increases risk, reuse breaches uniqueness",
            "first_source_alarm": alarm,
            "source_detection_delay_operations": None if not raw_fault or alarm is None or alarm < config.fault_onset else alarm-config.fault_onset,
            "prefault_or_healthy_alarm": raw_healthy_alarm,
            "first_policy_violation": first,
            "first_demonstrated_confidentiality_compromise": compromised,
            "first_demonstrated_integrity_forgery": forged,
            "integrity_demonstration_available": algorithm.startswith("AES"),
            "observation_horizon": config.operations, "repeated_operations": repeats,
            "collision_pairs": sum(n*(n-1)//2 for n in counts.values()),
            "affected_messages_including_earlier_use": sum(n for n in counts.values() if n > 1),
            "source_gate_mitigation_time": alarm,
            "source_gate_blocks_first_misuse": source_prevented,
            "source_gate_prevention_with_eligible_source_warning": source_prevented and raw_fault and not raw_healthy_alarm,
            "source_gate_emitted_operations": len(source_emitted_nonces),
            "source_gate_remaining_repeated_operations": len(source_emitted_nonces)-len(set(source_emitted_nonces)),
            "guard_first_block": first, "guard_blocked_operations": blocked,
            "guard_emitted_operations": guarded_emissions,
            "guard_repeated_emissions": len(emitted_nonces)-len(set(emitted_nonces)),
            "guard_capacity_entries": len(guarded.guard.seen),
            "alarm_to_compromise_operations": None if alarm is None or compromised is None else compromised-alarm,
            "source_output_sha256": source_digest}


def nonce_failures(config, progress=lambda _: None):
    paths = []
    for rep in range(config.calibration_repetitions):
        src = ReplaySource(seed_for(config.namespace+"/calibration", rep, "nonce-source"))
        paths.append([src.read_bytes(32) for _ in range(config.operations)])
    calibration = calibrate(paths)
    names = ["healthy"]+[f"low_{b}" for b in config.nonce_bits]+[
        "fixed", "serial", "drift", "rollback", "counter_reset", "counter_persistent"]
    rows, raw_metrics = [], []
    with tempfile.TemporaryDirectory(prefix="rnglab-blue-nonces-") as directory:
        for rep in range(config.repetitions):
            for name in names:
                raw, nonces, digest = _nonce_path(config, rep, name, directory)
                scores = source_scores(raw)
                hits = np.flatnonzero(scores > calibration["threshold"])
                alarm = int(hits[0]) if len(hits) else None
                raw_metrics.append({"replicate": rep, "condition": name,
                                    "max_score": float(scores.max()),
                                    "first_alarm": alarm,
                                    "first_postonset_alarm": next((int(i) for i in hits if i >= config.fault_onset), None),
                                    "scores": scores.tolist()})
                for algorithm in ALGORITHMS:
                    rows.append(_aead_path(config, rep, name, algorithm, nonces, alarm, digest))
            if rep % 16 == 0:
                progress(f"nonce trajectories: {rep+1}/{config.repetitions} paired replicates")
    return {"calibration": calibration, "rows": rows, "source_metrics": raw_metrics}


def _summarize_nonce(config, data):
    summary = {}
    # The three algorithms SHARE each nonce trajectory; do not triple the n.
    rows = [r for r in data["rows"] if r["algorithm"] == "AES-128-GCM"]
    healthy = [r for r in rows if r["condition"] == "healthy"]
    names = list(dict.fromkeys(r["condition"] for r in rows))
    for name in names:
        group = [r for r in rows if r["condition"] == name]
        events = [r["first_policy_violation"] for r in group]
        alarms = [r["first_source_alarm"] is not None for r in group]
        leads = [r["alarm_to_compromise_operations"] for r in group
                 if r["alarm_to_compromise_operations"] is not None and not r["prefault_or_healthy_alarm"]]
        delays = [r["source_detection_delay_operations"] for r in group if r["source_detection_delay_operations"] is not None]
        summary[name] = {
            "first_collision_survival": survival_summary(events, config.operations),
            "source_alarm_rate": event_rate(alarms),
            "prefault_or_healthy_alarm_rate": event_rate([r["prefault_or_healthy_alarm"] for r in group]),
            "source_gate_prevents_first_misuse_rate_all_paths": event_rate([r["source_gate_blocks_first_misuse"] for r in group]),
            "source_gate_eligible_warning_prevention_rate_all_paths": event_rate([r["source_gate_prevention_with_eligible_source_warning"] for r in group]),
            "missed_collision_by_source_rate_all_paths": event_rate([r["first_policy_violation"] is not None and not r["source_gate_blocks_first_misuse"] for r in group]),
            "guard_repeated_emissions": sum(r["guard_repeated_emissions"] for r in group),
            "source_gate_repeated_emissions": sum(r["source_gate_remaining_repeated_operations"] for r in group),
            "paired_source_alarm_change_from_healthy": paired_difference(alarms, [r["first_source_alarm"] is not None for r in healthy]),
            "conditional_lead_operations": {"n": len(leads), "median": float(np.median(leads)) if leads else None,
                                             "range": [min(leads), max(leads)] if leads else None,
                                             "selection": "both events observed, excluding prefault alarms; descriptive only"},
            "conditional_detection_delay_operations": {"n": len(delays), "median": float(np.median(delays)) if delays else None,
                                                         "range": [min(delays), max(delays)] if delays else None},
            "mean_affected_messages": float(np.mean([r["affected_messages_including_earlier_use"] for r in group])),
            "mean_source_gate_blocked_operations": float(np.mean([config.operations-r["source_gate_emitted_operations"] for r in group])),
        }
        if name.startswith("low_"):
            bits = int(name[4:])
            summary[name]["analytical_postfault_iid_birthday_probability"] = birthday_probability(config.operations-config.fault_onset, bits)
            summary[name]["analytical_scope"] = "uniform seed support; SHA256 truncation assumed injective on toy support; cross-prefix collisions neglected"
    return summary


def run_campaign(config, progress=lambda _: None):
    start = time.perf_counter()
    progress("bounded weak-key experiments")
    weak = weak_keys(config)
    progress("conditioning and known-seed prediction")
    blind = conditioning_blind_spot(config)
    progress("nonce and defensive-control trajectories")
    nonces = nonce_failures(config, progress)
    progress("state exposure, recovery, and signatures")
    reseed = recovery(config)
    signatures = [{"replicate": r, **signature_demo(seed_for(config.namespace, r, "signatures"))}
                  for r in range(config.repetitions)]
    controls = [{"replicate": r, "algorithm": a, "reuse": reuse,
                 **nonce_misuse_demo(a, seed_for(config.namespace, r, "aead-control"), reuse=reuse)}
                for r in range(config.repetitions) for a in ALGORITHMS for reuse in (False, True)]
    summary = {"nonce_conditions": _summarize_nonce(config, nonces), "weak_keys": {}, "conditioning": {}}
    for b in config.weak_seed_bits:
        selected = [r for r in weak if r["seed_bits"] == b]
        summary["weak_keys"][str(b)] = {"recovery_rate": event_rate([r["recovered"] for r in selected]),
            "mean_attempts": float(np.mean([r["candidate_attempts"] for r in selected])),
            "mean_duplicates_per_64": float(np.mean([r["duplicates"] for r in selected])),
            "mean_partial_exposure_attempts": float(np.mean([r["partial_exposure_attempts"] for r in selected]))}
    for c in ("sha256", "hkdf"):
        for e in ("hmac", "chacha20"):
            selected = [r for r in blind if r["conditioner"] == c and r["expander"] == e]
            summary["conditioning"][c+"/"+e] = {
                "weak_output_nonrejection": event_rate([not r["weak_output"]["rejected"] for r in selected]),
                "reference_output_nonrejection": event_rate([not r["reference_output"]["rejected"] for r in selected]),
                "exact_prediction": event_rate([r["attacker_exact_prediction"] for r in selected])}
    labels = list(summary["conditioning"])
    p_values = [float(binomtest(summary["conditioning"][label]["weak_output_nonrejection"]["events"],
                              config.repetitions, .9, alternative="greater").pvalue) for label in labels]
    corrected = holm_bonferroni_correction(p_values, .05)
    healthy_rate = summary["nonce_conditions"]["healthy"]["source_alarm_rate"]
    summary["decisions"] = {
        "conditioning_nonrejection_above_90percent": {
            label: {"one_sided_binomial_p": p_values[i], "holm_adjusted_p": corrected.adjusted_p_values[i],
                    "supported": corrected.reject_decisions[i]} for i, label in enumerate(labels)},
        "healthy_source_alarm_cp95_upper_at_most_10percent": healthy_rate["interval"][1] <= .10,
        "healthy_source_alarm_cp95_upper": healthy_rate["interval"][1],
        "calibration_interpretation": "scoped synthetic check; failure means insufficient evidence, never silently retune",
        "all_measured_gcm_collision_candidates_verified": all(
            r["first_demonstrated_integrity_forgery"] == r["first_policy_violation"]
            for r in nonces["rows"] if r["algorithm"].startswith("AES")),
        "no_fresh_nonce_control_recovery_or_forgery": all(
            not r["xor_plaintext_recovered"] and r["new_message_forgery_accepted"] is not True
            for r in controls if not r["reuse"])}
    results = {"schema_version": 1, "config": config.to_dict(), "summary": summary,
               "weak_keys": weak, "conditioning": blind, "nonces": nonces,
               "recovery": reseed, "signatures": signatures, "aead_controls": controls}
    return results, {"wall_seconds": time.perf_counter()-start,
                     "qualification": "single host campaign duration; not comparative performance validation"}
