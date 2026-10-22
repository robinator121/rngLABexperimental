# Supervisor summary

TRIDENT now has executable mathematical, physics and failure-risk experiments, with full traces and fresh held-out seeds. The most important mathematical boundary is that modest-looking acceptance corruption can create absorbing states, invalidating a unique stationary-response model. VMC energy can remain exact while position sampling fails. PIMC formula correctness still does not certify short-chain uncertainty.

The first forecasting endpoint generated no future failures and could not assess prediction. A separately frozen amendment with fresh seeds was executed. The final test and OOD Brier values and theory-feature ablation are reported in V016_PREDICTIVE_RELIABILITY.md; no universal early-warning or continuation guarantee follows.

This is a bounded v0.16 research candidate. Next: independently review the implementation and calibrate one observable/one reliability endpoint within a justified perturbation domain. Avoid scaling campaigns until the finite-run uncertainty and scientific acceptance criteria hold.


## Final numerical assessment

- At L4 near Tc, healthy energy/spin was -1.56336, acceptance-bias energy -2, and coupled-draw energy -1.69727. Absorption explains the strong biased case; it is not a small stationary perturbation.
- 18/24 paired Ising energy contrasts pass the specified Holm family. This does not establish a universal cross-model law.
- VMC biased proposal still gives exact energy .5; x² bias and chain uncertainty are separately recorded.
- Final healthy M16 PIMC CV intervals cover 26/32. The exact 95% proportion interval is [0.63560769013591, 0.9279238345415024]. The nominal95% target is not established; the historical defect remains relevant.
- Test: theory Brier .08869430 versus combined .08867247; difference +.00002183, paired95% CI[-.00010771,.00015345]. The required improvement .02 was not reached.
- OOD: theory Brier .13912729 versus combined .13126058; difference+.0078667, secondary pointwise CI[.0052413,.0104164]. Theory was worse in this cohort.
- Every fitted/constant-risk conservative policy permitted0 continuations; this is not demonstrated useful reliability optimization.

The amended endpoint generated many future threshold crossings even in HEALTHY trajectories. `FORECAST_UNIT_TABLE.json` gives the split/family denominators and early failures. A rolling finite-sample fluctuation past .15 is an operational event, not proof of RNG-caused application damage. This endpoint confounds normal thermal variability and reliability loss and needs redesign before a practical forecasting claim. Source/application threshold scores and alarm timing were added only as clearly labeled post-hoc diagnostics. An application alarm defined by one threshold crossing will mechanically precede a two-crossing failure; such lead time is not genuine predictive discovery.
