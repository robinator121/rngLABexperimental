# Opt-in source and fault interface

BLUE sources implement `read_bytes(n)` using RNGLab14's single accounting boundary. Every successful byte returned contributes once to `consumed_bytes` and `stream_digest`; a configured byte limit fails before over-consumption. Empty reads do not advance faults. Replay manifests disclose public seeds and zero validated entropy credit. OS output is available only through an explicitly instantiated, read-only reference source.

```python
from rnglab_crypto.sources import ReplaySource, seed_for
from rnglab_crypto.faults import FaultAdapter, FaultSpec

source = ReplaySource(seed_for("my-registered-study", 0, "seed-input"))
faulted = FaultAdapter(
    source,
    FaultSpec("bias", onset=64, duration=192, intensity=0.6),
    lab_only=True,
)
seed_input = faulted.read_bytes(32)
```

The adapter accepts only synthetic `ReplaySource`, with explicit `lab_only=True`. It cannot wrap the OS reference. One nonempty read attempt is one operation, including an entropy-starvation attempt. The caller's request schedule is part of replay. Each fault is active on `[onset, onset+duration)`; recovery resumes the simultaneously advancing healthy shadow. This is a counterfactual laboratory intervention and not an assurance of actual state recovery.

| Fault | Layer (default; explicit override recorded when used) | Mechanism and severity |
|---|---|---|
| `reduced_entropy` | DRBG initialization; raw source in nonce campaign | Mask input to a declared b-bit support, zero-pad to requested seed width |
| `seed_reuse` | DRBG initialization | Repeat first active seed block; fixed-width reads required |
| `fixed_nonce` | Downstream nonce API | Return fixed bytes; campaign injects downstream while preserving raw-source observations |
| `bias` | Raw source | Bernoulli bit probability .5+.5×intensity |
| `serial` | Raw source | Balanced stationary binary Markov model, P(next=1)=.5+rho×(previous-.5) |
| `periodic` | Raw source | Repeat a seeded byte pattern of `period_bytes`; preserve pattern offset |
| `rollback` | DRBG state | Restore initial synthetic HMAC-DRBG state at fault onset |
| `starvation` | Raw source | Raise `EntropyUnavailable`; no silent substitute bytes |
| `drift` | Raw source | Bias strength grows linearly across the active interval |
| `partial_exposure` | Bounded seed/state model | Emit b-bit supported seed, expose a known prefix through `.leak`; unknown suffix remains bounded |

Intensity zero disables any fault. For structural faults, positive intensity enables the specified mechanism; `entropy_bits`, `exposed_bits`, period and duration define severity. Bias/serial/drift intensity varies continuously in [0,1]. Onset/duration and support widths are strictly validated; unknown mechanisms/layers and nonfinite parameters are rejected.

Use `condition`/`expand` only in the research harness; no conditioner earns entropy credit. `EntropyGate` consumes explicit assumed/external evidence and quarantines on insufficient budget or a health alarm. `GuardedAEAD` reserves each nonce before invoking the primitive. `NonceLedger` is a separate transactional allocation model. These are local opt-in APIs, with no global monkeypatches, background hooks or changes to production cryptographic services.
