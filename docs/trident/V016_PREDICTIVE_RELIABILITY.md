# Predictive reliability

Forecast at fixed sweep150 for failure within next200, not after a source alarm. Failure=two successive scheduled rolling64 energy errors above.15 per spin relative to exact L3 equilibrium energy; checks every20 sweeps beginning64. Initial failure before/equal cutoff makes a path ineligible and is counted separately. Complete follow-up is available to every requested horizon, so binary Brier needs no censor weighting; no future events means only censored event-time information, not known infinite survival.

Three training/validation/test families:healthy,abrupt,drift; three wholly excluded OOD families:periodic,coupled,correlation. Independent SeedSequence child replicas, distinct base seeds per split/family. No input normalization or fitting on test. L2-regularized logistic source-only, application-only, combined and theory models use fixed feature sets. Theory response is a deterministic exact-kernel lookup, not fitted to B/C outcomes; clamped beyond its domain for scoring, explicitly abstaining for continuation decisions. Other models share a conservative development-support range gate. Permit only when predicted risk<.1 and gate passes; this risk threshold is operational, not certified.

Single-class development data now rejects fitting. Run002 zero-event Brier scores are retained as a feasibility failure and establish no forecasting skill.

{
  "validation": {
    "n": 115,
    "events": 106,
    "scores": {
      "source": 0.07334723739874742,
      "application": 0.07603405833546628,
      "combined": 0.07620669480377015,
      "theory": 0.07633832106113791,
      "constant": 0.07327909321201466,
      "never_fail": 0.9217391304347826
    },
    "primary": {
      "difference": 0.0001316262573677688,
      "paired_bootstrap_ci95": [
        -5.659821584477547e-05,
        0.0003329560431553946
      ],
      "H3_supported": false,
      "unit": "one eligible trajectory; IID seed groups, conditional on fixed fitted model"
    }
  },
  "test": {
    "n": 226,
    "events": 204,
    "scores": {
      "source": 0.08754504868750128,
      "application": 0.08921045232050967,
      "combined": 0.08867247243982007,
      "theory": 0.0886943037353636,
      "constant": 0.08808584912609303,
      "never_fail": 0.9026548672566371
    },
    "primary": {
      "difference": 2.183129554352934e-05,
      "paired_bootstrap_ci95": [
        -0.0001077067391369006,
        0.00015344521899757652
      ],
      "H3_supported": false,
      "unit": "one eligible trajectory; IID seed groups, conditional on fixed fitted model"
    }
  },
  "OOD": {
    "n": 143,
    "events": 125,
    "scores": {
      "source": 0.1354145221580072,
      "application": 0.10740794731870285,
      "combined": 0.13126057716392628,
      "theory": 0.13912729428797752,
      "constant": 0.11022041276202985,
      "never_fail": 0.8741258741258742
    },
    "primary": {
      "difference": 0.007866717124051247,
      "paired_bootstrap_ci95": [
        0.005241278214905116,
        0.010416437039686751
      ],
      "H3_supported": false,
      "unit": "one eligible trajectory; IID seed groups, conditional on fixed fitted model"
    }
  }
}

![Scores](figures/forecast_brier.png)

![Calibration](figures/calibration.png)

Ablation primary: theory-minus-combined Brier, with paired trajectory bootstrap2000 and fresh bootstrap seed. H3 requires upper95%CI<-.02; no result is promoted automatically. A bin with few examples is not evidence of conditional calibration. Support and theory-domain abstention rates, permitted counts and failures permitted are in analysis.json. No external shift guarantee, no learned exact failure time. Abrupt onset after cutoff is intentionally not available in features; inability to anticipate an unannounced external shock is an identifiability limit.


## Final numerical assessment

- At L4 near Tc, healthy energy/spin was -1.56336, acceptance-bias energy -2, and coupled-draw energy -1.69727. Absorption explains the strong biased case; it is not a small stationary perturbation.
- 18/24 paired Ising energy contrasts pass the specified Holm family. This does not establish a universal cross-model law.
- VMC biased proposal still gives exact energy .5; x² bias and chain uncertainty are separately recorded.
- Final healthy M16 PIMC CV intervals cover 26/32. The exact 95% proportion interval is [0.63560769013591, 0.9279238345415024]. The nominal95% target is not established; the historical defect remains relevant.
- Test: theory Brier .08869430 versus combined .08867247; difference +.00002183, paired95% CI[-.00010771,.00015345]. The required improvement .02 was not reached.
- OOD: theory Brier .13912729 versus combined .13126058; difference+.0078667, secondary pointwise CI[.0052413,.0104164]. Theory was worse in this cohort.
- Every fitted/constant-risk conservative policy permitted0 continuations; this is not demonstrated useful reliability optimization.

The amended endpoint generated many future threshold crossings even in HEALTHY trajectories. `FORECAST_UNIT_TABLE.json` gives the split/family denominators and early failures. A rolling finite-sample fluctuation past .15 is an operational event, not proof of RNG-caused application damage. This endpoint confounds normal thermal variability and reliability loss and needs redesign before a practical forecasting claim. Source/application threshold scores and alarm timing were added only as clearly labeled post-hoc diagnostics. An application alarm defined by one threshold crossing will mechanically precede a two-crossing failure; such lead time is not genuine predictive discovery.
