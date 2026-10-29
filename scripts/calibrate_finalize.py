"""Link supplemental evidence into reports without changing frozen outcomes."""
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[2]


def append_once(name, heading, body):
    path = ROOT / name
    original = path.read_text().split('\n\n' + heading)[0]
    path.write_text(original.rstrip() + '\n\n' + heading + '\n\n' + body.strip() + '\n')


def main():
    probe = json.loads((ROOT / 'SUSCEPTIBILITY_PROBE.json').read_text())
    q = next(r for r in probe['rows'] if r['cell'] == 'ising_L4_T0')
    susceptibility = (
        f"At L4, T=1.5, the exact zero-field susceptibility is {q['reference']:.6f}. "
        f"Across 512 healthy chains the within-chain plug-in estimate averaged "
        f"{q['within_chain_plugin']['mean']:.6f}, a bias of "
        f"{q['within_chain_error']['mean']:.6f} (between-chain 95% interval "
        f"[{q['within_chain_error']['ci95'][0]:.6f}, {q['within_chain_error']['ci95'][1]:.6f}]). "
        "Subtracting the squared finite-chain magnetization is severely misleading when "
        "chains do not sample both magnetized phases. Accurate energy does not validate "
        "this response observable. This is healthy finite-run mixing/estimator error, "
        "not evidence of defective randomness.\n\n"
        f"An exploratory known-zero-mean symmetry estimator, N times mean(m²)/T, "
        f"averaged {q['known_zero_mean_symmetry_estimator']['mean']:.6f}, with bias interval "
        f"[{q['symmetry_error']['ci95'][0]:.6f}, {q['symmetry_error']['ci95'][1]:.6f}]. "
        "It uses exact finite-lattice zero-field spin-flip symmetry; this is not a "
        "generic estimator repair or a proof of second-moment equilibration. The "
        "post-confirmation calculation reuses saved traces, does not change any "
        "registered endpoint, and is not new confirmation. All nine cells are in "
        "SUSCEPTIBILITY_PROBE.json. Susceptibility and heat-capacity interval "
        "calibration still require separate validation."
    )
    for name in ['V017_ISING_VALIDATION.md', 'V017_INDEPENDENT_AUDIT.md',
                 'V017_NEGATIVE_FINDINGS.md', 'V017_SCIENTIFIC_REPORT.md',
                 'EXECUTIVE_SUMMARY.md']:
        append_once(name, '## Healthy susceptibility counterexample', susceptibility)

    tolerances = (
        "The frozen tolerances are transparent engineering choices for this benchmark, "
        "not externally established safety limits: energy error 0.1 J per spin; "
        "VMC x² error 0.1 against mean 0.5 (20%); PIMC energy error 0.1 in "
        "m=ħ=ω=1 units. Their practical importance depends on a downstream use case. "
        "They were fixed before confirmation, but broader healthy rank thresholds can "
        "make an inaccurate sampler appear normal relative to its own protocol. "
        "Absolute accuracy, valid uncertainty and calibrated anomaly probability must "
        "therefore be assessed separately. No user-level safety or application-loss "
        "threshold has been established."
    )
    append_once('V017_FAILURE_DEFINITIONS.md', '## Practical-tolerance limitation', tolerances)
    append_once('V017_HEALTHY_CALIBRATION.md', '## Practical-tolerance limitation', tolerances)

    warnings = json.loads((ROOT / 'WARNING_DIAGNOSTICS.json').read_text())
    lines = [
        '| Population | Warning | True warnings | False warnings | Missed events | Median declaration lead (sweeps) |',
        '|---|---|---:|---:|---:|---:|',
    ]
    for split, rows in warnings['splits'].items():
        for name, row in rows.items():
            lines.append(f"| {split} | {name} | {row['TP']} | {row['FP']} | {row['FN']} | "
                         f"{row['event_declaration_lead_sweeps']['median']} |")
    warning_text = '\n'.join(lines) + (
        '\n\nPost-confirmation descriptive reductions at previously fixed boundaries; '
        'no optimization or new confirmatory tests. Test contains 188 events and '
        '272 non-events; OOD contains 66 events and 33 non-events. Exact pointwise '
        'working binomial intervals and denominators are saved in WARNING_DIAGNOSTICS.json. '
        'These conditional mixture summaries do not provide distribution-free guarantees. '
        'The combined 0.1 boundary misses fewer test events but warns on 124/272 '
        'non-events; thresholds are not matched in false-positive probability. '
        'Do not claim a matched-FPR superiority from this table. All warnings are issued '
        'at sweep 1408. The reported 256-sweep median is time to a constructed '
        'two-window event declaration, with a 128–512 range among true warnings. '
        'Overlapping-window persistence creates mechanical lag: this does not '
        'establish 256 sweeps of usable warning before physical damage, successful '
        'mitigation, or an independent early-warning advantage. Support/abstention '
        'gates are applied in the separate four-action policy, not these raw alerts.'
    )
    append_once('V017_PREDICTIVE_RELIABILITY.md', '## Warning, miss and lead-time diagnostics', warning_text)

    (ROOT / 'ACADEMIC_SUMMARY.md').write_text('''# Supervisor discussion summary

RNGLab v0.17 is an executed finite-protocol calibration study, not a certified reliability system. It separates individual accuracy excursions from cohort excess loss, uncertainty undercoverage, structural transition failure and healthy-calibrated persistent anomalies. The original v0.16 healthy counts (49/128 early failures and 57/79 later failures) were reproduced from raw traces.

The new campaign used 12 primary Ising/VMC/PIMC cells, 399 healthy calibration and 512 independent healthy confirmation trajectories per cell, and 128 paired trajectories per fault arm. These are 18,996 simulation arms, not 18,996 independent replications: faulty arms share healthy seed groups. Seven mild-bias discrimination results pass both the multiplicity-adjusted effect gate and their healthy validation bound. Only 10/12 cells meet the frozen healthy validation ceiling, so global H1 is unsupported.

Independent reference arithmetic agrees for the finite physical models. This does not validate finite-run estimation: healthy PIMC CV intervals cover only 172/512 or 311/512 times with the tested short blocks. Healthy low-temperature L4 Ising susceptibility is estimated near 0.746 against exact 10.403 by the within-chain plug-in estimator, despite accurate energy. A model-specific symmetry diagnostic largely removes that bias, but is exploratory and not a general repair.

In a separately gated L3/Tc incident-forecast experiment, combined source/application Brier score improves from 0.24206 to 0.12078 on an eligible healthy/drift population. Theory features add no demonstrated benefit. Abrupt faults have already failed at issuance, and most unseen-fault trajectories are also excluded by early failure. The four-action policy continues 143/460 test paths with four later events, but costs 1.117 against always-stop cost 1 under the declared unmitigated loss. Unseen-family Brier advantage is inconclusive and continuation is zero. Full H4 and H5 remain unsupported.

The useful contribution is executable failure semantics, causal negative controls and preserved negative evidence. A stronger publication requires external independent review, defensible dependent-sample uncertainty, application-relevant tolerances and a held-out intervention study. Mathematical references and tools are established literature; the current evidence does not establish a new universal forecasting law.
''')
    (ROOT / 'V018_PRIORITIES.md').write_text('''# Scientifically justified v0.18 priorities

1. **Repair healthy estimator validity before richer forecasting** (3–5 engineering/research days). Add explicit symmetry-aware and generic susceptibility definitions, nonlinear time-series uncertainty and PIMC batch-size/convergence controls. Acceptance: derive definitions, reproduce original failures, freeze new controls and validate held-out coverage with stated uncertainty. Retain v0.17 negative results; do not patch its saved outcomes.
2. **Calibrate a narrower deployment domain** (2–4 days, dependent on item 1). Choose a downstream tolerance with a scientific user; independently recalibrate the failed L4/Tc and PIMC16 cells, plus cluster-specific endpoints. Acceptance: separate exact-accuracy, false-alarm and power gates all pass under a new seed namespace. Do not pool away failed cells.
3. **Measure whether intervention helps** (3–5 days, dependent on item 2). Randomize eligible paths to continued execution, source repair, reset and stop under a frozen cost model. Include reset overhead and loss of accumulated samples. Acceptance: a paired or randomized held-out decision-loss improvement over credible baselines; preserve a null result.
4. **Forecast useful incident failures rather than survivors of abrupt damage** (3–5 days). Freeze slower, nonstructural fault onset and multiple prospective issuance times with trajectory-grouped evaluation. Compare source/application/static/simple-hazard baselines at matched false-positive targets. Report whole-cohort exclusions, censoring, calibration and actual mitigation lead time. Require independent unseen-family confirmation; no post-test boundary tuning.
5. **External independent verification and packaging** (1–3 engineering days plus reviewer availability). A different reviewer reconstructs primary endpoints from delivered traces and inspects physics/statistics. Publish no universal guarantee, entropy certification or commercial operational claim from these results alone.

Effort estimates exclude waiting time and large simulation campaigns. Do not add physical hardware, a dashboard, another orchestration framework or extra predictor complexity before these gates. The present package is a reproducible development candidate with partial validation.
''')
    print('Supplemental findings and academic summary linked; frozen results unchanged.')


if __name__ == '__main__':
    main()
