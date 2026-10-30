#!/usr/bin/env python3
"""Measured fault-study plots, preserving first-alarm and censoring distinctions.

Reads the frozen fault-study-v1 schema. Does not alter experiment data. Delays
are explicitly conditional on observing both a post-onset threshold crossing
and a first post-onset application failure; they are not forecast lead times.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

ROOT=Path(__file__).resolve().parent
CONDITIONS=["healthy","acceptance_bias","channel_coupling"]
LABELS={"healthy":"Healthy","acceptance_bias":"Acceptance shift","channel_coupling":"Channel coupling"}
COLORS={"healthy":"#218C7C","acceptance_bias":"#C17F1C","channel_coupling":"#8A5C9C"}
MONITORS=["source_marginal","source_cross_channel","application_shape"]
MONITOR_LABELS={"source_marginal":"Source marginal means","source_cross_channel":"Source channel dependence","application_shape":"Application shape"}
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.titlesize":12,
                     "axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,
                     "grid.color":"#E0E4E8","grid.linewidth":.7,"axes.axisbelow":True,
                     "legend.frameon":False,"pdf.fonttype":42,"svg.fonttype":"none"})

def read(path):return json.loads(Path(path).read_text())

def save(fig,name,out,caption):
    fig.text(.5,.012,caption,ha="center",va="bottom",fontsize=8,color="#4E5761")
    for ext in ("png","pdf","svg"):
        fig.savefig(out/f"{name}.{ext}",bbox_inches="tight",**({"dpi":180} if ext=="png" else {}))
    plt.close(fig);print(f"wrote {name}.png/.pdf/.svg",flush=True)

def point_interval(ax,x,estimate,ci,color,vertical=True):
    """Draw intervals explicitly: percentile CIs need not contain the estimate."""
    low,high=ci
    if vertical:
        ax.vlines(x,low,high,color=color,lw=1.7);ax.hlines([low,high],x-.045,x+.045,color=color,lw=1.4)
        ax.plot(x,estimate,"o",color=color,ms=5)
    else:
        ax.hlines(x,low,high,color=color,lw=1.7);ax.vlines([low,high],x-.045,x+.045,color=color,lw=1.4)
        ax.plot(estimate,x,"o",color=color,ms=5)

def primary(rows,effects,out):
    chosen={r["condition"]:r for r in rows if r["observable"]=="CV" and r["method"]=="BM512"}
    fig,axs=plt.subplots(2,2,figsize=(10,7.5));x=np.arange(3)
    for i,c in enumerate(CONDITIONS):
        r=chosen[c];color=COLORS[c]
        point_interval(axs[0,0],i,r["coverage"],r["coverage_ci95"],color)
        axs[0,0].annotate(f"{r['covered']}/{r['replicas']}",(i,r["coverage"]),xytext=(0,10),textcoords="offset points",ha="center",fontsize=9)
        point_interval(axs[0,1],i,r["bias"],r["bias_ci95"],color)
    axs[0,0].axhline(.95,color="#9A423A",ls="--",lw=1.2);axs[0,0].set(title="CV intended-target interval coverage",ylabel="Empirical 95% interval coverage",ylim=(-.03,1.08))
    axs[0,1].axhline(0,color="#9A423A",ls="--",lw=1.2);axs[0,1].set(title="CV mean bias",ylabel="Mean energy − finite-M target")
    for ax in axs[0]:ax.set_xticks(x);ax.set_xticklabels([LABELS[c].replace(" ","\n") for c in CONDITIONS])
    faults=CONDITIONS[1:]
    for i,c in enumerate(faults):
        e=effects[c]
        point_interval(axs[1,0],i,e["paired_mean_energy_change"],e["paired_mean_change_CI95"],COLORS[c])
        point_interval(axs[1,1],i,e["paired_coverage_change"],e["paired_coverage_change_bootstrap_CI95"],COLORS[c])
    for ax in axs[1]:
        ax.axhline(0,color="#9A423A",ls="--",lw=1.2);ax.set_xticks([0,1]);ax.set_xticklabels([LABELS[c] for c in faults]);ax.set_xlim(-.5,1.5)
    axs[1,0].set(title="Paired energy change from healthy",ylabel="Fault − healthy mean energy")
    axs[1,1].set(title="Paired coverage change from healthy",ylabel="Fault − healthy coverage")
    fig.suptitle("Controlled offline faults after validated healthy sampling",fontsize=14)
    fig.tight_layout(rect=(0,.065,1,.95))
    save(fig,"fault_primary_endpoints",out,"512 paired replicas; 8,192 post-onset sweeps; BM512. Binomial coverage, t energy, and replica-bootstrap coverage-change marginal 95% intervals.")

def alert_rates(rows,out):
    fig,axs=plt.subplots(1,2,figsize=(11,4.9),gridspec_kw={"width_ratios":[1.6,1]})
    offsets=[-.17,0,.17]
    for offset,c in zip(offsets,CONDITIONS):
        for j,m in enumerate(MONITORS):
            r=next(r for r in rows if r["condition"]==c and r["monitor"]==m)
            point_interval(axs[0],j+offset,r["whole_run_rate"],r["whole_run_CI95"],COLORS[c])
        axs[0].plot([],[],"o",color=COLORS[c],label=LABELS[c])
    for j,m in enumerate(MONITORS):
        r=next(r for r in rows if r["condition"]=="healthy" and r["monitor"]==m)
        point_interval(axs[1],j,r["whole_run_rate"],r["whole_run_CI95"],COLORS["healthy"])
        axs[1].annotate(f"{r['events']}/{r['replicas']}",(j,r["whole_run_rate"]),xytext=(0,9),textcoords="offset points",ha="center",fontsize=8)
    labels=["Marginal\nmeans","Channel\ndependence","Application\nshape"]
    for ax in axs:
        ax.set_xticks(range(3));ax.set_xticklabels(labels);ax.axhline(.05,color="#9A423A",ls="--",lw=1.2)
    axs[0].set(ylim=(-.025,1.08),ylabel="Fraction of trajectories with ≥1 alert",title="Whole-run measured alert rates")
    axs[0].legend(ncol=1,loc="center right",fontsize=9)
    axs[1].axhline(.10,color="#757D85",ls=":",lw=1.2);axs[1].set(ylim=(0,.125),title="Healthy false alarms, expanded scale")
    fig.suptitle("Monitors share a nominal calibration target; actual healthy rates differ",fontsize=14)
    fig.tight_layout(rect=(0,.065,1,.94))
    save(fig,"fault_monitor_whole_run_rates",out,"Exact binomial 95% intervals; frozen thresholds from 399 independent healthy paths; 512 test paths per condition; no monitor-superiority claim.")

def as_int(row,key):
    value=row[key]
    return int(value) if value not in ("",None) else None

def alarm_partition(times,out):
    fig,axs=plt.subplots(1,3,figsize=(11.8,4.8),sharey=True)
    cats=[("Before onset", "#8C9BA7"),("After onset", "#355C7D"),("No alert by end", "#D9DFE4")]
    for ax,m in zip(axs,MONITORS):
        bottom=np.zeros(3)
        for idx,(label,color) in enumerate(cats):
            vals=[]
            for c in CONDITIONS:
                selected=[r for r in times if r["condition"]==c and r["monitor"]==m]
                a=np.array([as_int(r,"first_alert_sweep") for r in selected])
                count=int(((a>=0)&(a<=2048)).sum()) if idx==0 else int((a>2048).sum()) if idx==1 else int((a<0).sum())
                vals.append(count/len(selected))
            ax.bar(np.arange(3),vals,bottom=bottom,color=color,width=.62,label=label)
            bottom+=vals
        ax.set_title(MONITOR_LABELS[m]);ax.set_xticks(range(3));ax.set_xticklabels([LABELS[c].replace(" ","\n") for c in CONDITIONS]);ax.set_ylim(0,1)
    axs[0].set_ylabel("Fraction of 512 trajectories")
    handles,labels=axs[0].get_legend_handles_labels();fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,.94),ncol=3,fontsize=9)
    fig.suptitle("First whole-run alert: onset and administrative censoring",fontsize=14)
    fig.tight_layout(rect=(0,.07,1,.87))
    save(fig,"fault_first_alert_partition",out,"Fault onset: sweep 2,049; study end: 10,240. First-alert categories are mutually exclusive; absence of an alert is censored at study end.")

def failure_partition(effects,out):
    fig,ax=plt.subplots(figsize=(7.7,4.8));bottom=np.zeros(3)
    specs=[("before_onset","Before onset","#8C9BA7"),("new_after_onset","After onset","#C17F1C"),("right_censored","No failure by end","#D9DFE4")]
    for key,label,color in specs:
        vals=np.array([effects[c]["failure_summary"][key]/512 for c in CONDITIONS])
        ax.bar(np.arange(3),vals,bottom=bottom,color=color,width=.6,label=label)
        for j,value in enumerate(vals):
            if value>.06:ax.text(j,bottom[j]+value/2,f"{round(value*512)}/512",ha="center",va="center",fontsize=9)
        bottom+=vals
    for j,c in enumerate(CONDITIONS):
        summary=effects[c]["failure_summary"]
        point_interval(ax,j,summary["events"]/512,summary["CI95"],"#202833")
    ax.plot([],[],"o",color="#202833",label="Any failure, 95% CI")
    ax.set_xticks(range(3));ax.set_xticklabels([LABELS[c] for c in CONDITIONS]);ax.set_ylim(0,1)
    ax.set(ylabel="Fraction of trajectories",title="Prespecified numerical tolerance excursions and censoring")
    ax.legend(ncol=2,loc="upper center",bbox_to_anchor=(.5,-.1),fontsize=9)
    fig.tight_layout(rect=(0,.1,1,1))
    save(fig,"fault_application_failure_partition",out,"Failure: two successive 1,024-sweep CV errors >0.05 (stride 256); censored at 10,240. Points: exact binomial 95% CIs for any failure.")

def timing_summary(times):
    result=[]
    for c in CONDITIONS:
        for m in MONITORS:
            selected=[r for r in times if r["condition"]==c and r["monitor"]==m]
            first=np.array([as_int(r,"first_alert_sweep") for r in selected]);post=np.array([as_int(r,"first_postonset_alert_sweep") for r in selected]);failure=np.array([as_int(r,"failure_sweep") for r in selected])
            both=(post>=0)&(failure>2048);delay=post[both]-failure[both]
            firstboth=(first>2048)&(failure>2048);firstdelay=first[firstboth]-failure[firstboth]
            row={"condition":c,"monitor":m,"replicas":len(selected),
                 "preonset_first_alerts":int(((first>=0)&(first<=2048)).sum()),
                 "postonset_first_alerts":int((first>2048).sum()),"no_first_alert":int((first<0).sum()),
                 "postonset_crossings_including_realerts":int((post>=0).sum()),
                 "prealerted_then_crossed_again":int(((first>=0)&(first<=2048)&(post>=0)).sum()),
                 "preonset_failures":int(((failure>=0)&(failure<=2048)).sum()),
                 "postonset_failures":int((failure>2048).sum()),"failure_right_censored":int((failure<0).sum()),
                 "postonset_failure_and_crossing_both_observed":int(both.sum()),
                 "postonset_failure_without_postonset_crossing":int(((failure>2048)&(post<0)).sum()),
                 "conditional_crossing_before_failure":int((delay<0).sum()),
                 "conditional_simultaneous":int((delay==0).sum()),"conditional_crossing_after_failure":int((delay>0).sum()),
                 "conditional_crossing_minus_failure_median":float(np.median(delay)) if len(delay) else None,
                 "conditional_crossing_minus_failure_q25_q75":np.quantile(delay,[.25,.75]).tolist() if len(delay) else None,
                 "postonset_whole_run_first_alert_and_failure_both_observed":int(firstboth.sum()),
                 "whole_run_first_alert_before_failure_among_both":int((firstdelay<0).sum()),
                 "whole_run_first_alert_simultaneous_among_both":int((firstdelay==0).sum()),
                 "whole_run_first_alert_after_failure_among_both":int((firstdelay>0).sum()),
                 "whole_run_first_alert_minus_failure_median_among_both":float(np.median(firstdelay)) if len(firstdelay) else None}
            assert row["preonset_failures"]+row["postonset_failures"]+row["failure_right_censored"]==len(selected)
            assert row["postonset_failure_and_crossing_both_observed"]+row["postonset_failure_without_postonset_crossing"]==row["postonset_failures"]
            result.append(row)
    return result

def conditional_delays(times,out):
    fig,axs=plt.subplots(2,3,figsize=(12.5,7.9),gridspec_kw={"height_ratios":[1.3,1]})
    categories=[("Both observed","#355C7D"),("Failure; no later crossing","#C17F1C"),
                ("Failed before onset","#8C9BA7"),("Failure censored","#D9DFE4")]
    all_delays=[]
    for j,m in enumerate(MONITORS):
        parts=[]
        for c in CONDITIONS:
            selected=[r for r in times if r["condition"]==c and r["monitor"]==m]
            both=[r for r in selected if as_int(r,"failure_sweep")>2048 and as_int(r,"first_postonset_alert_sweep")>=0]
            delay=np.sort(np.array([as_int(r,"first_postonset_alert_sweep")-as_int(r,"failure_sweep") for r in both]))
            all_delays.extend(delay.tolist())
            if len(delay):axs[0,j].step(delay,np.arange(1,len(delay)+1)/len(delay),where="post",color=COLORS[c],lw=1.8,label=f"{LABELS[c]} (n={len(delay)})")
            else:axs[0,j].plot([],[],color=COLORS[c],label=f"{LABELS[c]} (n=0)")
            parts.append([len(both),sum(as_int(r,"failure_sweep")>2048 and as_int(r,"first_postonset_alert_sweep")<0 for r in selected),
                          sum(0<=as_int(r,"failure_sweep")<=2048 for r in selected),sum(as_int(r,"failure_sweep")<0 for r in selected)])
        bottom=np.zeros(3);p=np.array(parts)/512
        for idx,(label,color) in enumerate(categories):
            axs[1,j].bar(np.arange(3),p[:,idx],bottom=bottom,color=color,width=.65,label=label);bottom+=p[:,idx]
        axs[0,j].axvline(0,color="#9A423A",ls="--",lw=1.2);axs[0,j].set_title(MONITOR_LABELS[m]);axs[0,j].set_ylim(0,1.04)
        axs[0,j].legend(fontsize=8,loc="lower right");axs[0,j].xaxis.set_major_locator(MaxNLocator(5))
        axs[0,j].set_xlabel("Crossing time − failure time (sweeps)")
        axs[1,j].set_xticks(range(3));axs[1,j].set_xticklabels([LABELS[c].replace(" ","\n") for c in CONDITIONS]);axs[1,j].set_ylim(0,1)
    if all_delays:
        low=min(-256,min(all_delays));high=max(256,max(all_delays));pad=.04*(high-low)
        for ax in axs[0]:ax.set_xlim(low-pad,high+pad)
    axs[0,0].set_ylabel("Conditional empirical CDF");axs[1,0].set_ylabel("Fraction of all 512 trajectories")
    handles,labels=axs[1,0].get_legend_handles_labels();fig.legend(handles,labels,loc="upper center",bbox_to_anchor=(.5,.505),ncol=4,fontsize=8)
    fig.suptitle("Post-onset threshold-crossing delays, conditioned on observed failures",fontsize=14)
    fig.tight_layout(rect=(0,.075,1,.95));fig.subplots_adjust(hspace=.65)
    save(fig,"fault_conditional_crossing_failure_delays",out,"Descriptive ECDFs, no bands. Later crossings include previously alerted paths; these are not forecast lead times. Bottom bars expose exclusion and right censoring.")

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--root",type=Path,default=ROOT);args=parser.parse_args()
    root=args.root;source=root/"fault-study-v1";out=root/"figures";out.mkdir(exist_ok=True)
    names=["analysis.json","monitor_results.json","effects.json","alarm_failure_times.csv"]
    missing=[name for name in names if not (source/name).exists()]
    if missing:print("Fault results not complete; waiting for: "+", ".join(missing));return
    rows=read(source/names[0]);monitors=read(source/names[1]);effects=read(source/names[2])
    with (source/names[3]).open(newline="") as f:times=list(csv.DictReader(f))
    primary(rows,effects,out);alert_rates(monitors,out);alarm_partition(times,out);failure_partition(effects,out);conditional_delays(times,out)
    summary={"source_file":str(source/names[3]),"source_sha256":hashlib.sha256((source/names[3]).read_bytes()).hexdigest(),
             "interpretation":"Conditional timing only; later crossings can follow a preonset alert. Separate first-whole-run-alert fields exclude those realerts. No RUS forecast evaluation.","rows":timing_summary(times)}
    (out/"fault_timing_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    files={p.name:{"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(out.glob("fault_*")) if p.is_file() and p.name!="fault_figure_manifest.json"}
    inputs={name:hashlib.sha256((source/name).read_bytes()).hexdigest() for name in names}
    if (source/"protocol.json").exists():inputs["protocol.json"]=hashlib.sha256((source/"protocol.json").read_bytes()).hexdigest()
    (out/"fault_figure_manifest.json").write_text(json.dumps({"plot_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"matplotlib_version":matplotlib.__version__,"inputs_sha256":inputs,"files":files},indent=2)+"\n")

if __name__=="__main__":main()
