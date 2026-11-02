"""Independent accounting/censoring audit plus labeled post-hoc finite-support checks."""
import argparse
import hashlib
import json
from pathlib import Path

from rnglab_crypto.monitoring import diagnostics
from rnglab_crypto.primitives.conditioning import condition, expand


def audit(run):
    run = Path(run)
    result = json.loads((run/"results.json").read_text())
    config = result["config"]
    checks = 0
    def require(condition, message):
        nonlocal checks
        if not condition:
            raise AssertionError(message)
        checks += 1
    manifest = json.loads((run/"manifest.json").read_text())
    for name, expected in manifest.items():
        require(hashlib.sha256((run/name).read_bytes()).hexdigest() == expected, f"hash {name}")
    rows = result["nonces"]["rows"]
    lookup = {(r["replicate"], r["condition"], r["algorithm"]): r for r in rows}
    names = result["summary"]["nonce_conditions"]
    require(len(lookup) == len(rows) == config["repetitions"]*len(names)*3, "unique experimental units")
    for r in rows:
        first = r["first_policy_violation"]
        require(first == r["first_demonstrated_confidentiality_compromise"], "actual plaintext recovery")
        require(r["guard_emitted_operations"]+r["guard_blocked_operations"] == config["operations"], "guard accounting")
        require(r["guard_repeated_emissions"] == 0, "guard blocked reuse")
        require(r["guard_blocked_operations"] == r["repeated_operations"], "duplicate accounting")
        if r["algorithm"].startswith("AES"):
            require(first == r["first_demonstrated_integrity_forgery"], "actual GCM forgery")
        else:
            require(r["first_demonstrated_integrity_forgery"] is None, "no unsupported ChaCha forgery")
        expected_raw_fault = r["condition"].startswith("low_") or r["condition"] in ("serial", "drift", "rollback")
        require(r["raw_source_fault_present"] == expected_raw_fault, "fault layer")
        if not expected_raw_fault:
            require(not r["source_gate_prevention_with_eligible_source_warning"], "no false credited detection")
        if r["condition"] in ("fixed", "counter_reset", "counter_persistent"):
            h = lookup[(r["replicate"], "healthy", r["algorithm"])]
            require(r["source_output_sha256"] == h["source_output_sha256"], "downstream fault preserves raw input")
            require(r["first_source_alarm"] == h["first_source_alarm"], "downstream source identifiability")
    threshold = result["nonces"]["calibration"]["threshold"]
    for m in result["nonces"]["source_metrics"]:
        first = next((j for j, x in enumerate(m["scores"]) if x > threshold), None)
        require(m["first_alarm"] == first, "causal alarm threshold")
    for name, summary in names.items():
        selected = [r for r in rows if r["condition"] == name and r["algorithm"] == "AES-128-GCM"]
        times = [r["first_policy_violation"] for r in selected]
        survival = summary["first_collision_survival"]
        require(survival["censored"] == times.count(None), "censoring count")
        for point in survival["curve"]:
            exact = sum(t is None or t+1 > point["operations"] for t in times)/len(times)
            require(abs(exact-point["survival"]) < 1e-12, "KM equals empirical common-horizon survival")
        require(summary["source_alarm_rate"]["events"] == sum(r["first_source_alarm"] is not None for r in selected), "alarm denominator")
    for r in result["aead_controls"]:
        if not r["reuse"]:
            require(not r["xor_plaintext_recovered"] and r["new_message_forgery_accepted"] is not True, "negative control")
    return {"checks_passed": checks, "results_sha256": manifest["results.json"],
            "scope": "independent accounting and censoring recomputation, not independent cryptographic implementation"}


def finite_support():
    # This is explicitly post-hoc verification, not a replacement confirmation sample.
    cells = {}
    for conditioner in ("sha256", "hkdf"):
        for expander in ("hmac", "chacha20"):
            rejects = []
            for value in range(256):
                raw = value.to_bytes(32, "big")
                data = expand(condition(raw, conditioner), 4096, expander)
                if diagnostics(data)["rejected"]:
                    rejects.append(value)
            cells[conditioner+"/"+expander] = {"nonrejected_seeds": 256-len(rejects),
                "seed_support": 256, "exact_nonrejection_fraction_for_this_finite_mapping": 1-len(rejects)/256,
                "rejected_public_seed_identifiers": rejects}
    injective = {}
    for bits in (8, 12, 16):
        nonces = {hashlib.sha256(b"BLUE nonce"+v.to_bytes(32, "big")).digest()[:12] for v in range(1 << bits)}
        injective[str(bits)] = {"input_count": 1 << bits, "unique_outputs": len(nonces), "injective_on_enumerated_support": len(nonces) == 1 << bits}
    return {"classification": "post-hoc exhaustive finite-domain diagnostic/mapping audit; not preregistered hypothesis evidence",
            "conditioning": cells, "nonce_seed_mapping": injective}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runs", nargs="+", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    report = {"runs": {p.name: audit(p) for p in args.runs}, "finite_support": finite_support()}
    with args.output.open("x") as handle:
        json.dump(report, handle, indent=2)
        handle.write("\n")
