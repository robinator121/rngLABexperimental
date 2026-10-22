# TRIDENT frozen protocol

The authoritative executable protocol is `confirmation-run004/config.json`, hash 95b928faa4f04b16998c9965cf07c2e48efe1c0526d9c21da89f7bbb06cf9bc0. It was written before this run produced measurements. This is local timestamped registration, not external preregistration.

{
  "version": "0.16.0.dev3",
  "status": "bounded_confirmation_with_exploratory_extensions",
  "A": {
    "seed": 180001,
    "replicas": 256,
    "n": 2048,
    "parameters": [
      [
        0.1,
        0.2
      ],
      [
        0.01,
        0.02
      ],
      [
        0.7,
        0.8
      ]
    ],
    "epsilons": [
      0.0001,
      0.001,
      0.01,
      0.05
    ]
  },
  "B": {
    "seed": 181000,
    "replicas": 32,
    "steps": 1200,
    "burn": 300,
    "lattices": [
      2,
      3,
      4,
      8
    ],
    "temperatures": [
      1.5,
      2.269185314213022,
      3.5
    ],
    "faults": [
      "healthy",
      "bias",
      "coupled"
    ],
    "acceptance_bias": 0.15
  },
  "quantum": {
    "seed": 182000,
    "replicas": 32,
    "steps": 3000,
    "burn": 500,
    "beads": [
      4,
      16
    ],
    "faults": [
      "healthy",
      "bias",
      "coupled"
    ]
  },
  "C": {
    "seed": 183000,
    "L": 3,
    "T": 2.269185314213022,
    "steps": 500,
    "cutoff": 150,
    "horizon": 200,
    "window": 64,
    "stride": 20,
    "tolerance": 0.15,
    "persistent_checks": 2,
    "severity": 0.3,
    "development_per_family": 64,
    "validation_per_family": 64,
    "test_per_family": 128,
    "families": [
      "healthy",
      "abrupt",
      "drift"
    ],
    "OOD": [
      "periodic",
      "coupled",
      "correlation"
    ]
  },
  "limits": {
    "threads": 1,
    "wall_seconds": 900,
    "retries": 0,
    "no_cloud": true
  },
  "hypotheses": {
    "H1": "exact sensitivity remainder is quadratic away from singular a+b=0",
    "H2": "bias changes Ising energy; observable response differs by model",
    "H3": "theory features improve heldout Brier over source+application by at least .02"
  },
  "theory_domain": {
    "delta_max": 0.025,
    "C_amendment": "run002 had zero future events at tolerance .5; new tolerance .15 and cutoff150 frozen with fresh seeds before run004",
    "reason": "attempt001 lost stationary uniqueness at larger acceptance bias; subsequent replay is amended development, not independent confirmation"
  },
  "inference": {
    "unit": "independent SeedSequence child chain; paired arms within same base namespace",
    "alpha": 0.05,
    "multiple_testing": "Holm across Ising paired energy tests; Brier primary contrast paired trajectory bootstrap",
    "bootstrap_repetitions": 2000,
    "coverage": "pointwise exact binomial; no universal guarantee",
    "confirmation": "no threshold or features tuned on final test",
    "power": "32 paired B units detect about 0.51 paired SD at unadjusted 80% power; Holm reduces power. H3 n384 IID rows gives rough Brier difference SE<=0.051 worst-case; small gains may remain inconclusive."
  }
}

Run001 stopped on nonunique stationary theory. Run002 is development evidence; forecasting had zero future events with tolerance .5. Run003 was an accidental duplicate after a working-directory error and adds no independent evidence. Run004 uses fresh A/B/quantum/C seeds, with the stated C endpoint amendment. Do not pool attempts. Primary B tests use paired chain energy differences and Holm correction across 24 contrasts. C primary contrast is theory-minus-combined Brier on eligible test trajectories; OOD contrasts secondary. No guarantee follows from these finite simulations.

Runtime cost was a correctness pilot plus three approximately minute-long CPU executions, not a power-optimized large campaign. No running external campaign was stopped.
