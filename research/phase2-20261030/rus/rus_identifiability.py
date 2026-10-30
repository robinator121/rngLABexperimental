#!/usr/bin/env python3
"""Synthetic, alarm-conditioned RUS identifiability counterexample.

This is a constructed Bernoulli Monte Carlo workload, not a PIMC result or a
successful prediction experiment. Both possible worlds have identical observed
prefixes and alarm times. An unobserved post-alarm law changes future reliability.
"""
from __future__ import annotations
import csv, hashlib, json, math, platform, random, statistics, sys
from datetime import datetime, timezone
from pathlib import Path

OUT=Path(__file__).resolve().parent
PREFIX_NAMESPACE="rnglab-phase2-rus-identical-prefix-20261030-v1"
FUTURE_NAMESPACE="rnglab-phase2-rus-independent-futures-20261030-v1"
T0=513                  # first degraded sample; 1-based
ALARM_WINDOW=64
ALARM_COUNT=48
BLOCK=64               # application Monte Carlo sample count per estimate
FAIL_LOW,FAIL_HIGH=16,48 # |block estimate - .5| >= .25
HORIZON_BLOCKS=32
REPLICAS=1024

def seed(label):
    digest=hashlib.sha256(label.encode()).hexdigest()
    return int(digest[:32],16),digest

def wilson(k,n):
    z=statistics.NormalDist().inv_cdf(.975); p=k/n; d=1+z*z/n
    c=(p+z*z/(2*n))/d; r=z/d*math.sqrt(p*(1-p)/n+z*z/(4*n*n))
    return [max(0.,c-r),min(1.,c+r)]

def failure_probability(p):
    return sum(math.comb(BLOCK,k)*p**k*(1-p)**(BLOCK-k)
               for k in range(BLOCK+1) if k<=FAIL_LOW or k>=FAIL_HIGH)

def first_failure(bits):
    for j in range(HORIZON_BLOCKS):
        count=sum(bits[j*BLOCK:(j+1)*BLOCK])
        if count<=FAIL_LOW or count>=FAIL_HIGH: return (j+1)*BLOCK
    return None

def make_svg(rows,path):
    w,h=820,460; left,top,pw,ph=82,75,680,300
    parts=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}">',
           '<rect width="100%" height="100%" fill="white"/>',
           '<style>text{font-family:Arial,sans-serif;fill:#202a35}.grid{stroke:#dde2e8}.axis{stroke:#59636e}</style>',
           '<text x="410" y="27" text-anchor="middle" font-size="18" font-weight="bold">Synthetic RUS counterexample: one prefix, two future laws</text>',
           '<text x="410" y="50" text-anchor="middle" font-size="13">Exact conditional failure CDF; no fitted forecast or PIMC measurement</text>']
    for tick in (0,.25,.5,.75,1):
        y=top+(1-tick)*ph
        parts += [f'<line class="grid" x1="{left}" y1="{y}" x2="{left+pw}" y2="{y}"/>',
                  f'<text x="{left-12}" y="{y+4}" text-anchor="end" font-size="12">{tick:.2f}</text>']
    for k in (0,8,16,24,32):
        x=left+pw*k/32
        parts.append(f'<text x="{x}" y="{top+ph+22}" text-anchor="middle" font-size="12">{k*BLOCK}</text>')
    colors={"restored":"#1f73b7","persistent":"#c54d24"}
    for row in rows:
        q=row["exact_per_block_failure_probability"]
        coordinates=[(left,top+ph)]
        for k in range(1,33):
            x=left+pw*k/32
            coordinates.extend([(x,top+ph*(1-q)**(k-1)),(x,top+ph*(1-q)**k)])
        points=" ".join(f'{x:.2f},{y:.2f}' for x,y in coordinates)
        parts.append(f'<polyline points="{points}" fill="none" stroke="{colors[row["scenario"]]}" stroke-width="3"/>')
    parts += [f'<text x="{left+pw/2}" y="{h-31}" text-anchor="middle" font-size="13">Additional sampling steps after the alarm (not ESS)</text>',
              '<text x="20" y="225" transform="rotate(-90 20 225)" text-anchor="middle" font-size="13">P(application failure by step r | observed prefix)</text>',
              '<text x="490" y="145" fill="#c54d24" font-size="14">Persistent p = 0.9</text>',
              '<text x="490" y="350" fill="#1f73b7" font-size="14">Restored p = 0.5</text>', '</svg>']
    path.write_text("\n".join(parts)+"\n")

def main():
    prefix_seed,prefix_digest=seed(PREFIX_NAMESPACE)
    rng=random.Random(prefix_seed); prefix=[]
    for t in range(1,10001):
        p=.5 if t<T0 else .9
        prefix.append(int(rng.random()<p))
        if len(prefix)>=ALARM_WINDOW and sum(prefix[-ALARM_WINDOW:])>=ALARM_COUNT:
            ta=t; break
    else: raise RuntimeError("alarm not reached")
    assert ta>=T0
    prefix_bytes=bytes(prefix)
    prefix_hash=hashlib.sha256(prefix_bytes).hexdigest()
    (OUT/"identical_observed_prefix.csv").write_text("sample_1based,bernoulli_observation\n"+"".join(f"{i+1},{v}\n" for i,v in enumerate(prefix)))
    future_seeds=[]; records=[]; scenarios={"restored":.5,"persistent":.9}
    for replica in range(REPLICAS):
        sd,digest=seed(f"{FUTURE_NAMESPACE}/{replica:04d}")
        future_seeds.append(dict(replica=replica,seed_sha256=digest,seed128_decimal=str(sd)))
        stream=random.Random(sd)
        uniforms=[stream.random() for _ in range(BLOCK*HORIZON_BLOCKS)]
        for scenario,p in scenarios.items():
            future=[int(u<p) for u in uniforms]
            remaining=first_failure(future)
            records.append(dict(replica=replica,scenario=scenario,alarm_time=ta,
                                remaining_steps=remaining if remaining is not None else "",
                                observed_time=min(remaining if remaining is not None else BLOCK*HORIZON_BLOCKS,BLOCK*HORIZON_BLOCKS),
                                failure_event=int(remaining is not None),
                                administrative_censoring=int(remaining is None)))
    rows=[]
    for scenario,p in scenarios.items():
        recs=[r for r in records if r["scenario"]==scenario]; events=sum(r["failure_event"] for r in recs)
        q=failure_probability(p)
        rows.append(dict(scenario=scenario,post_alarm_p=p,replicas=REPLICAS,events=events,
                         censored=REPLICAS-events,empirical_failure_by_horizon=events/REPLICAS,
                         wilson_95=wilson(events,REPLICAS),exact_per_block_failure_probability=q,
                         exact_failure_by_horizon=1-(1-q)**HORIZON_BLOCKS,
                         exact_mean_remaining_steps=BLOCK/q,
                         exact_median_remaining_steps=BLOCK*max(1,math.ceil(math.log(.5)/math.log1p(-q))),
                         observed_event_remaining_steps_median=statistics.median(r["remaining_steps"] for r in recs if r["failure_event"]) if events else None))
    with (OUT/"conditional_future_outcomes.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    (OUT/"future_seed_manifest.json").write_text(json.dumps(future_seeds,indent=2)+"\n")
    make_svg(rows,OUT/"conditional_failure_cdf.svg")
    report=dict(experiment_id="phase2-rus-identifiability-counterexample-v1",
                classification="Synthetic constructed Bernoulli workload; not a PIMC result; not successful forecasting.",
                timestamp_utc=datetime.now(timezone.utc).isoformat(),
                software=dict(python=sys.version,platform=platform.platform(),dependencies="stdlib only"),
                prefix=dict(degradation_onset_1based=T0,alarm_time_1based=ta,observed_samples=len(prefix),
                            observations_sha256=prefix_hash,seed_sha256=prefix_digest,
                            pre_onset_p=.5,onset_to_alarm_p=.9,
                            alarm_rule=f"First rolling {ALARM_WINDOW}-sample source count >= {ALARM_COUNT}"),
                target=dict(reference_mean=.5,application_observable="Mean of next disjoint 64 Bernoulli samples",
                            failure_rule="First post-alarm 64-sample block with absolute mean error >= .25",
                            failure_time="T = t_a + 64 K, K first failing block",
                            remaining_useful_samples="R = max(0,T-t_a) in additional sampling steps; not ESS",
                            forecast_issuance="Once, at first alarm, using the complete identical observed prefix",
                            horizon_steps=BLOCK*HORIZON_BLOCKS,censoring="Right censor at 2048 additional samples",
                            abstention="Required if a single informative forecast is requested without a post-alarm law or identified mixture probability"),
                scenarios=rows,
                experimental_units="1024 independently seeded future continuations conditional on one fixed observed prefix; two scenarios paired by uniforms; not 2048 independent prefixes",
                identifiability="Both worlds assign identical laws to all observations through the alarm, but different conditional laws to R. Prefix-only observations do not identify which future law applies or a mixing weight between them.",
                outcome_files=["identical_observed_prefix.csv","conditional_future_outcomes.csv","future_seed_manifest.json","conditional_failure_cdf.svg"])
    (OUT/"rus_identifiability_results.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(dict(prefix=report["prefix"],scenarios=rows),indent=2))

if __name__=="__main__": main()
