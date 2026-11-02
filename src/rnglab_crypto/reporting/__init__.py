"""Machine-readable evidence, portable plots and deterministic report generation."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+"\n")


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def code_manifest():
    src = Path(__file__).resolve().parents[2]
    paths = [p for name in ("rnglab_crypto", "rnglab11", "rnglab14", "rnglab15", "rnglab17")
             for p in sorted((src/name).rglob("*.py"))]
    return {str(p.relative_to(src)): digest(p) for p in paths}


def environment():
    from cryptography.hazmat.backends.openssl.backend import backend
    return {"python": platform.python_version(), "platform": platform.platform(),
            "openssl": backend.openssl_version_text(),
            "dependencies": {name: importlib.metadata.version(name) for name in
                             ("numpy", "scipy", "cryptography", "matplotlib")}}


def freeze_run(out, config):
    """Freeze design and source hashes BEFORE generating any evaluation data."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=False)
    methodology = Path(__file__).resolve().parents[3]/"docs/crypto/EXPERIMENT_PROTOCOL.md"
    write_json(out/"protocol.json", {"status": "frozen_before_evaluation", "config": config.to_dict(),
               "source_sha256": code_manifest(), "methodology": "docs/crypto/EXPERIMENT_PROTOCOL.md",
               "methodology_sha256": digest(methodology) if methodology.exists() else None,
               "index_convention": "zero-based operations; null means not observed, never infinity",
               "seed_derivation": "SHA256 compact JSON [BLUE-v1, namespace, replicate, role]",
               "public_replay_seeds_are_not_real_secrets": True})
    return out


def _rate_text(rate):
    lo, hi = rate["interval"]
    return f'{rate["events"]}/{rate["n"]} ({100*lo:.1f}–{100*hi:.1f}% exact 95% CI)'


def report_markdown(results):
    config, summary = results["config"], results["summary"]
    lines = ["# RNGLab BLUE: measured synthetic campaign", "",
             f'Seed namespace: `{config["namespace"]}`. Independent paired replicates: {config["repetitions"]}. '
             f'Operations per trajectory: {config["operations"]}; injected onset: {config["fault_onset"]} (zero-based).', "",
             "These are local synthetic experiments. Statistical nonrejection does not establish security. "
             "All replay material is public; ‘healthy’ means the declared secret-unexposed simulation model.", "",
             "## Weak-key experiment", "",
             "| Uniform seed support | Successful bounded recovery | Mean guesses | Mean duplicate keys / 64 |",
             "|---|---|---|---|"]
    for b, s in summary["weak_keys"].items():
        lines.append(f'| {b} bits | {_rate_text(s["recovery_rate"])} | {s["mean_attempts"]:.1f} | {s["mean_duplicates_per_64"]:.2f} |')
    lines += ["", "Only the intentionally bounded seed dictionaries were searched. No production-strength AES key search was attempted.",
              "", "## Conditioning blind spot", "",
              "| Conditioner / expander | Known-seed output: three diagnostics do not reject | Reference output: do not reject | Exact attacker prediction |",
              "|---|---|---|---|"]
    for name, s in summary["conditioning"].items():
        lines.append(f'| {name} | {_rate_text(s["weak_output_nonrejection"])} | {_rate_text(s["reference_output_nonrejection"])} | {_rate_text(s["exact_prediction"])} |')
    lines += ["", "Monobit, block-frequency and runs diagnostics use Holm correction within each stream at 1%. "
              "This is a diagnostic subset, not the complete NIST STS battery. Intervals across rows are marginal, not simultaneous.",
              "", "## Nonce failure and detection", "",
              "| Condition | Observed collision trajectories | Source alarm trajectories | Repeated emissions with source gate | Repeated emissions with nonce guard |",
              "|---|---|---|---|---|"]
    for name, s in summary["nonce_conditions"].items():
        lines.append(f'| {name} | {_rate_text(s["first_collision_survival"]["event_rate"])} | '
                     f'{_rate_text(s["source_alarm_rate"])} | {s["source_gate_repeated_emissions"]} | {s["guard_repeated_emissions"]} |')
    lines += ["", "Source alarms during downstream fixed-nonce/counter-reset faults are alarms on an unchanged healthy raw source, "
              "not successful source-fault detections. Incidental/pre-fault quarantine can prevent misuse but is excluded from credited early warnings.", "",
              "The same nonce paths are paired across all three AEADs; they are not three independent samples. "
              "For each first collision, the pipeline tests exact known-plaintext recovery; GCM additionally tests a new one-block forgery. "
              "ChaCha20-Poly1305 forgery is not implemented or claimed as measured.", "",
              "| Algorithm | Collision paths | Demonstrated plaintext recovery | Accepted new-message forgery |",
              "|---|---|---|---|"]
    for a in ("AES-128-GCM", "AES-256-GCM", "ChaCha20-Poly1305"):
        rows = [r for r in results["nonces"]["rows"] if r["algorithm"] == a]
        collisions = sum(r["first_policy_violation"] is not None for r in rows)
        disclosure = sum(r["first_demonstrated_confidentiality_compromise"] is not None for r in rows)
        forgery = str(sum(r["first_demonstrated_integrity_forgery"] is not None for r in rows)) if a.startswith("AES") else "not attempted"
        lines.append(f'| {a} | {collisions} | {disclosure} | {forgery} |')
    decision = summary["decisions"]
    lines += ["", "Counts pooled over fault conditions are descriptive, not estimates of real-world prevalence.", "",
              f'Frozen healthy source-monitor calibration check: **{"PASS" if decision["healthy_source_alarm_cp95_upper_at_most_10percent"] else "INSUFFICIENT EVIDENCE"}** '
              f'(CP95 upper false-alarm bound {decision["healthy_source_alarm_cp95_upper"]:.4f}; criterion <=0.10).', "",
              'Four-cell conditioning criterion (one-sided exact binomial versus 90% nonrejection, Holm family alpha .05): '
              +", ".join(f'{k}: {"supported" if v["supported"] else "insufficient evidence"}' for k, v in decision["conditioning_nonrejection_above_90percent"].items())+".", "",
              "## Comparison matrix", "",
              "| RNG condition | Statistical tests | Cryptographic outcome | Security implication | Detection | Mitigation |",
              "|---|---|---|---|---|---|",
              "| Known 8-bit seed, SHA256/HKDF then HMAC/ChaCha expansion | Measured nonrejection rates above | Exact output prediction measured | Unpredictability absent when seed is known | Raw low-support health screen alarms; provenance identifies known seed | Require independently justified startup entropy; public salt adds no entropy |",
              "| Bounded weak key seed | Output appearance is not the recovery criterion | AES-256-GCM known-plaintext dictionary recovery and duplicate Ed25519 keys measured | Key length exceeds effective seed uncertainty | Model entropy-budget gate blocks | OS/library key generation with trustworthy entropy infrastructure |",
              "| Reduced nonce seed support | Raw monitor rates above | First collision, XOR disclosure and GCM forgery measured on colliding paths | Nonce uniqueness failure is distinct from key weakness | Source monitor plus per-key reservation guard | Durable counter or bounded guard |",
              "| Downstream fixed nonce / counter reset | Raw source remains the paired reference | Misuse and its consequences measured | Healthy raw randomness cannot prevent an API/state-management fault | Nonce policy monitor at attempted reuse | Guard before encryption; shared persistent ledger |",
              "| Unique persistent counter across simulated process restart | Counter output is intentionally structured; source reference unchanged | No observed reuse in this finite campaign | Unpredictable nonces are unnecessary for these AEADs if unique per key | Counter reservation invariants | Commit allocation before use; preserve state or rotate key |",
              "| Known DRBG state plus known reseed | Output tests cannot establish secrecy | Clone continues to predict, measured | Unvalidated reseed cannot restore unknown state | State provenance / policy; no statistical proof | Withheld model input plus rekey changes clone outcome; past ciphertext remains exposed |",
              "| ECDSA repeated signing nonce | No statistical battery needed | Synthetic private scalar recovered; new signature verifies | Message-specific nonce misuse compromises ECDSA key | Nonce policy; repeated r across distinct messages is evidence | RFC6979 deterministic signing; adequate key entropy |",
              "| RFC6979 / Ed25519 deterministic signing | Repeated signatures of same message expected | Same-message determinism and distinct-message behavior measured | No fresh signing entropy needed by these constructions | Library known-answer checks | Protect key generation and private state |", "",
              "## Signatures, recovery and controls", ""]
    for label, group in (("Signatures", results["signatures"]), ("Recovery", results["recovery"])):
        lines.append(f"**{label}:**")
        lines.append("")
        for key in group[0]:
            if isinstance(group[0][key], bool):
                lines.append(f'- `{key}`: {sum(r[key] for r in group)}/{len(group)}.')
        lines.append("")
    controls = [r for r in results["aead_controls"] if not r["reuse"]]
    lines += [f'Fresh-nonce negative controls: {sum(r["xor_plaintext_recovered"] for r in controls)}/{len(controls)} '
              'successful XOR recoveries; '
              f'{sum(r["new_message_forgery_accepted"] is True for r in controls)} accepted GCM forgeries. '
              'Bit-flip rejection is recorded separately from forged-tag acceptance.', "",
              "## Stability window and limits", "",
              "The source score cutoff is frozen from separate healthy paths using a maximum over all monitored times and score components. "
              "Whole-path false alarms are measured on the evaluation controls. A nonce guard acts at the attempted policy violation; "
              "its prevention is not forecasting. Null events are administratively right-censored. "
              "Kaplan–Meier curves and exact common-horizon binomial intervals are recorded without extrapolating censored tails.", "",
              "Alarm-to-compromise differences are conditional descriptive windows. A known seed is insecure immediately under the threat model, "
              "even before an observable collision. Identical byte histories can accompany different attacker knowledge or different future resets: "
              "no output-only universal warning time is identifiable. Birthday-risk curves are analytical results conditional on IID uniform finite support.", "",
              "The initial suite demonstrates infrastructure failures and scoped prevention. It provides no break of AES or ChaCha20, "
              "no certified entropy estimate, no validated physical TRNG, and no deployment-wide forecast guarantee. "
              "See the frozen protocol, per-replicate JSON, source hashes, audit, and security limitations.", ""]
    return "\n".join(lines)


def plots(results, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"svg.hashsalt": "rnglab-blue-v1", "font.size": 10})
    summary = results["summary"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), layout="constrained")
    for name, s in summary["nonce_conditions"].items():
        if name.startswith("low_") or name in ("healthy", "counter_persistent", "drift"):
            curve = s["first_collision_survival"]["curve"]
            axes[0].step([v["operations"] for v in curve], [v["survival"] for v in curve], where="post", label=name)
    axes[0].set(xlabel="Operations observed", ylabel="Fraction with no observed nonce collision", ylim=(-.02, 1.02))
    axes[0].axvline(results["config"]["fault_onset"], color="gray", ls=":", label="fault onset")
    axes[0].legend(fontsize=8)
    names = list(summary["nonce_conditions"])
    vals = [summary["nonce_conditions"][n]["source_alarm_rate"] for n in names]
    means = np_array([v["rate"] for v in vals])
    errors = np_array([[v["rate"]-v["interval"][0], v["interval"][1]-v["rate"]] for v in vals]).T
    axes[1].errorbar(means, range(len(names)), xerr=errors, fmt="o", color="#21618c", capsize=3)
    axes[1].set(yticks=range(len(names)), yticklabels=names, xlabel="Source alarm rate (exact 95% CI)", xlim=(-.05, 1.05))
    fig.suptitle("RNGLab BLUE · observed nonce collisions and source detection")
    for suffix in ("png", "svg"):
        fig.savefig(Path(out)/f"nonce_detection.{suffix}", dpi=160, metadata={"Date": None} if suffix == "svg" else None)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), layout="constrained")
    rows = [r for r in results["nonces"]["rows"] if r["algorithm"] == "AES-128-GCM" and r["replicate"] == 0]
    for j, row in enumerate(rows):
        for field, marker, color, label in (("fault_onset", "|", "gray", "injection"),
                 ("first_source_alarm", "o", "#e67e22", "source alarm"),
                 ("first_policy_violation", "x", "#c0392b", "nonce reuse / demonstrated compromise")):
            if row[field] is not None:
                axes[0].scatter(row[field], j, marker=marker, color=color, label=label if j == 1 else None)
    axes[0].set(yticks=range(len(rows)), yticklabels=[r["condition"] for r in rows], xlabel="Zero-based operation")
    axes[0].set_title("Prespecified replicate 0 (illustration only)")
    axes[0].legend(fontsize=7)
    names = list(summary["conditioning"])
    means = [summary["conditioning"][n]["weak_output_nonrejection"]["rate"] for n in names]
    axes[1].barh(names, means, color="#21618c", label="Statistical nonrejection")
    axes[1].scatter([1]*len(names), range(len(names)), color="#c0392b", marker="x", label="Exact known-seed prediction")
    axes[1].set(xlim=(0, 1.08), xlabel="Fraction of independent replicates")
    axes[1].legend(fontsize=8, loc="lower left")
    fig.suptitle("Detection times are distinct from attacker uncertainty")
    for suffix in ("png", "svg"):
        fig.savefig(Path(out)/f"timelines_conditioning.{suffix}", dpi=160, metadata={"Date": None} if suffix == "svg" else None)
    plt.close(fig)


def np_array(values):
    import numpy as np
    return np.asarray(values)


def finish_run(out, results, performance, *, make_plots=True):
    out = Path(out)
    protocol = json.loads((out/"protocol.json").read_text())
    if protocol["source_sha256"] != code_manifest():
        raise RuntimeError("source changed after protocol freeze")
    if protocol["config"] != json.loads(json.dumps(results["config"])):
        raise RuntimeError("results configuration differs from frozen protocol")
    write_json(out/"results.json", results)
    write_json(out/"environment.json", environment())
    write_json(out/"performance.json", performance)
    (out/"RESULTS.md").write_text(report_markdown(results))
    if make_plots:
        plots(results, out)
    write_json(out/"manifest.json", {p.name: digest(p) for p in sorted(out.iterdir()) if p.is_file()})


def verify(out):
    out = Path(out)
    manifest = json.loads((out/"manifest.json").read_text())
    if not manifest or any(Path(name).name != name for name in manifest):
        raise ValueError("invalid artifact manifest")
    for name, expected in manifest.items():
        if digest(out/name) != expected:
            raise ValueError(f"artifact checksum mismatch: {name}")
    actual = {p.name for p in out.iterdir() if p.is_file()} - {"manifest.json"}
    if actual != set(manifest):
        raise ValueError("artifact inventory mismatch")
    return {"verified_artifacts": len(manifest), "results_sha256": manifest["results.json"]}
