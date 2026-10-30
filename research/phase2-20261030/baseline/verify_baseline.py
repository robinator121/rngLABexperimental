#!/usr/bin/env python3
"""Independent, stdlib-only verification of the preserved Phase-I PIMC artifacts.

Does not import or execute the original simulator. Derives finite-M references
from the Gaussian precision eigenvalues; obtains Student-t quantiles by direct
density integration (adaptive Simpson), independent of the pilot beta-function
implementation. Reads the saved observations and recomputes every interval.
"""
from __future__ import annotations
import csv, gzip, hashlib, json, math, platform, statistics, sys, time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[1]
INPUTS = ["pimc_uncertainty_pilot.py", "pimc-uncertainty-pilot-results.json",
          "pimc-uncertainty-pilot-chain-summary.csv", "pimc-uncertainty-pilot-traces.csv.gz",
          "pimc-uncertainty-pilot-coverage.svg", "PIMC_uncertainty_pilot_methodology.md",
          "RNGLab_research_status_2026-10-30.md", "RNGLab_literature.bib",
          "physical-inventory-20261019-140519.json", "physical-inventory-20261019-143018.json"]

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""): h.update(block)
    return h.hexdigest()

def integrate(f, a, b, tol=2e-12):
    fa, fm, fb = f(a), f((a+b)/2), f(b)
    s = (b-a)*(fa+4*fm+fb)/6
    def rec(a,b,fa,fm,fb,s,tol,depth):
        m=(a+b)/2; lm=(a+m)/2; rm=(m+b)/2
        fl,fr=f(lm),f(rm)
        sl=(m-a)*(fa+4*fl+fm)/6; sr=(b-m)*(fm+4*fr+fb)/6
        err=sl+sr-s
        if depth<=0 or abs(err)<=15*tol: return sl+sr+err/15
        return rec(a,m,fa,fl,fm,sl,tol/2,depth-1)+rec(m,b,fm,fr,fb,sr,tol/2,depth-1)
    return rec(a,b,fa,fm,fb,s,tol,24)

def tcrit(df):
    c=math.exp(math.lgamma((df+1)/2)-math.lgamma(df/2))/math.sqrt(df*math.pi)
    f=lambda x:c*(1+x*x/df)**(-(df+1)/2)
    lo,hi=0.,5.
    for _ in range(48):
        mid=(lo+hi)/2
        if integrate(f,0,mid)<.475: lo=mid
        else: hi=mid
    return (lo+hi)/2

def wilson(k,n):
    z=statistics.NormalDist().inv_cdf(.975); p=k/n
    d=1+z*z/n; c=(p+z*z/(2*n))/d
    r=z/d*math.sqrt(p*(1-p)/n+z*z/(4*n*n))
    return [c-r,c+r]

def corr(x,y):
    a=statistics.fmean(x); b=statistics.fmean(y)
    xx=sum((v-a)**2 for v in x); yy=sum((v-b)**2 for v in y)
    return sum((u-a)*(v-b) for u,v in zip(x,y))/math.sqrt(xx*yy)

def acf(x,lag):
    a=statistics.fmean(x); d=sum((v-a)**2 for v in x)
    return sum((x[i]-a)*(x[i+lag]-a) for i in range(len(x)-lag))/d

def main():
    start=time.monotonic()
    original=json.loads((BASE/INPUTS[1]).read_text())
    conf=original["configuration"]; nrep=conf["replicas"]; nsweep=conf["sweeps_per_replica"]
    m,beta=conf["beads"],conf["beta"]
    eig=[beta/m+4*m/beta*math.sin(math.pi*k/m)**2 for k in range(m)]
    ex2=sum(1/v for v in eig)/m
    espring=sum(4*math.sin(math.pi*k/m)**2/eig[k] for k in range(m))
    td=m/(2*beta)-m*espring/(2*beta*beta)+.5*ex2
    cv=1/(2*beta)+ex2-.5/beta
    reference=ex2
    assert abs(td-cv)<1e-13
    assert abs(reference-original["reference"]["finite_bead_energy"])<1e-13
    traces={r:{"TD":[],"CV":[]} for r in range(nrep)}
    with gzip.open(BASE/INPUTS[3],"rt",newline="") as f:
        reader=csv.DictReader(f)
        nrows=0
        for row in reader:
            r=int(row["replica"]); i=int(row["sweep_1based"])
            assert 0<=r<nrep and i==len(traces[r]["TD"])+1
            for key,col in [("TD","thermodynamic_energy"),("CV","centroid_virial_energy")]:
                value=float(row[col]); assert math.isfinite(value)
                traces[r][key].append(value)
            nrows+=1
    assert nrows==nrep*nsweep
    assert all(len(v["TD"])==nsweep for v in traces.values())
    with (BASE/INPUTS[2]).open(newline="") as f: chainmeta=list(csv.DictReader(f))
    assert len(chainmeta)==nrep
    seeds=[]; full_seeds=[]
    for row in chainmeta:
        text=f"{conf['seed_namespace']}/{int(row['replica']):03d}"
        digest=hashlib.sha256(text.encode("ascii")).hexdigest()
        assert digest==row["seed_sha256"]
        seeds.append(int(digest[:32],16)); full_seeds.append(digest)
    assert len(set(seeds))==nrep
    trace_hashes=[hashlib.sha256(repr(traces[r]).encode()).hexdigest() for r in range(nrep)]
    assert len(set(trace_hashes))==nrep
    critical={df:tcrit(df) for df in (5,15,23,63,95,255)}
    qerrors={str(df):critical[df]-original["student_t_quantile_checks"][str(df)]["quantile"] for df in critical}
    assert max(map(abs,qerrors.values()))<1e-9
    coverage=[]; details=[]
    for saved in original["coverage"]:
        window,obs,b=saved["window"],saved["observable"],saved["batch_size"]
        lo,hi=conf["windows"][window]["sweeps_1based"]; estimates=[]; ses=[]; count=0
        for r in range(nrep):
            x=traces[r][obs][lo-1:hi]; n=len(x); a=n//b
            bm=[sum(x[j*b:(j+1)*b])/b for j in range(a)]
            mean=sum(bm)/a
            # Explicit centered sum of squares rather than statistics.stdev.
            varmean=sum((v-mean)**2 for v in bm)/(a-1)/a
            se=math.sqrt(varmean); q=critical[a-1]
            covered=mean-q*se<=reference<=mean+q*se
            count+=covered; estimates.append(mean); ses.append(se)
            details.append(dict(window=window,observable=obs,batch_size=b,replica=r,
                                estimate=mean,variance_estimate=varmean,lower=mean-q*se,
                                upper=mean+q*se,covered=int(covered)))
        emp=statistics.variance(estimates); est=statistics.fmean(s*s for s in ses)
        result=dict(window=window,observable=obs,batch_size=b,covered=count,replicas=nrep,
                    coverage=count/nrep,wilson_95=wilson(count,nrep),
                    mean_chain_estimate=statistics.fmean(estimates),
                    mean_bias=statistics.fmean(estimates)-reference,
                    empirical_variance_chain_means=emp,mean_estimated_variance=est,
                    estimated_to_empirical_variance_ratio=est/emp,
                    median_mcse=statistics.median(ses),
                    matches_saved_coverage=count==saved["covered"])
        assert result["matches_saved_coverage"]
        assert abs(result["mean_chain_estimate"]-saved["mean_chain_estimate"])<1e-13
        assert abs(result["median_mcse"]-saved["median_chain_mcse"])<1e-13
        assert max(abs(a-b) for a,b in zip(result["wilson_95"],saved["wilson_95"]))<1e-13
        coverage.append(result)
    lag_summary=[]
    for obs in ("TD","CV"):
        xs=[traces[r][obs][4096:] for r in range(nrep)]
        for lag in (1,8,32,128,512,1024):
            vals=[acf(x,lag) for x in xs]
            lag_summary.append(dict(observable=obs,lag=lag,mean_acf=statistics.fmean(vals),
                                    median_acf=statistics.median(vals),min_acf=min(vals),max_acf=max(vals)))
    # Descriptive synchronized correlations only. Serial correlation reduces the
    # effective sample count, so these are not IID-based independence tests.
    xs=[traces[r]["CV"][4096:] for r in range(nrep)]
    cross=[corr(xs[i],xs[j]) for i in range(nrep) for j in range(i+1,nrep)]
    cross_sorted=sorted(cross)
    cross_summary=dict(pair_count=len(cross),mean=statistics.fmean(cross),
                       rms=math.sqrt(statistics.fmean(v*v for v in cross)),
                       min=min(cross),max=max(cross),median=statistics.median(cross),
                       percentile_025=cross_sorted[int(.025*len(cross))],
                       percentile_975=cross_sorted[int(.975*len(cross))],
                       interpretation="Descriptive only; no calibrated independence test because each trace is strongly autocorrelated.")
    source_hashes={name:dict(bytes=(BASE/name).stat().st_size,sha256=sha(BASE/name)) for name in INPUTS}
    report=dict(experiment_id="phase2-baseline-independent-verification-v1",
                timestamp_utc=datetime.now(timezone.utc).isoformat(),
                software=dict(python=sys.version,platform=platform.platform(),dependencies="stdlib only"),
                source_revision=original["source_revision_reference"],
                validation_status="All saved coverage counts, means, median MCSEs, Wilson intervals and seed derivations verified.",
                source_files=source_hashes,trace_rows=nrows,unique_replica_traces=len(set(trace_hashes)),
                seed_count=nrep,unique_128bit_seed_count=len(set(seeds)),seed_sha256=full_seeds,
                configuration=conf,
                gaussian_reference=dict(precision_eigenvalues=eig,expected_x2=ex2,expected_TD=td,expected_CV=cv,
                                        pilot_reference_difference=reference-original["reference"]["finite_bead_energy"]),
                independent_student_t_quantile_differences=qerrors,coverage=coverage,
                long_window_acf=lag_summary,cross_replica_cv_correlations=cross_summary,
                acceptance=dict(mean=statistics.fmean(float(x["acceptance_rate"]) for x in chainmeta),
                                min=min(float(x["acceptance_rate"]) for x in chainmeta),
                                max=max(float(x["acceptance_rate"]) for x in chainmeta)),
                detected_issues=["Saved JSON describes 'one bead per sweep'; code actually proposes every bead once per sweep in fixed order.",
                                 "Short-window batch 512 uses only six batches; allowed by standalone four-batch threshold but below general RNGLab eight-batch minimum.",
                                 "Short and long windows share trajectories and change both warmup and retained length; no independent replication or isolated length effect.",
                                 "Distinct SHA-derived PRNG seeds and nonduplicate traces do not prove statistical independence."],
                elapsed_seconds=time.monotonic()-start)
    (HERE/"baseline_verification.json").write_text(json.dumps(report,indent=2)+"\n")
    with (HERE/"recomputed_intervals.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(details[0])); writer.writeheader(); writer.writerows(details)
    print(json.dumps(dict(coverage=coverage,reference=report["gaussian_reference"],acf=lag_summary,cross=cross_summary,elapsed_seconds=report["elapsed_seconds"]),indent=2))

if __name__=="__main__": main()
