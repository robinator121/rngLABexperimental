#!/usr/bin/env python3
"""Plot completed RNGLab Phase-II study outputs without modifying study.py.

Run with the isolated environment's Python. Missing stages are skipped; rerun
after a stage finishes. All figures show measured results, with interval types
stated. Known-mean pooled ACFs use study.py's saved diagnostic, not a re-fit.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

ROOT=Path(__file__).resolve().parent
COLORS={"local_stationary_d04":"#355C7D","local_stationary":"#355C7D",
        "centroid_refresh":"#218C7C","iid_gaussian":"#C17F1C","reversible_scan":"#8A5C9C"}
LABELS={"local_stationary_d04":"Local Metropolis","local_stationary":"Local Metropolis",
        "centroid_refresh":"Centroid refresh","iid_gaussian":"Independent Gaussian",
        "reversible_scan":"Reversible scan"}
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.titlesize":12,
                     "axes.labelsize":10,"axes.spines.top":False,"axes.spines.right":False,
                     "axes.grid":True,"grid.color":"#E0E4E8","grid.linewidth":.7,
                     "axes.axisbelow":True,"legend.frameon":False,"pdf.fonttype":42,
                     "ps.fonttype":42,"svg.fonttype":"none","savefig.facecolor":"white"})

def read(path):
    return json.loads(Path(path).read_text())

def err(row,key,ci):
    center=row[key];low,high=row[ci]
    return np.array([[center-low],[high-center]])

def save(fig,name,out,caption):
    fig.text(.5,.012,caption,ha="center",va="bottom",fontsize=8,color="#4E5761")
    for ext in ("png","pdf","svg"):
        kwargs={"dpi":180} if ext=="png" else {}
        fig.savefig(out/f"{name}.{ext}",bbox_inches="tight",**kwargs)
    plt.close(fig)
    print(f"wrote {name}.png/.pdf/.svg",flush=True)

def bm_number(row):
    method=row["method"]
    return int(method[2:]) if method.startswith("BM") and method[2:].isdigit() else None

def exploratory_batch(rows,out):
    conditions=["local_stationary_d04","reversible_scan","centroid_refresh","iid_gaussian"]
    fig,ax=plt.subplots(figsize=(8,4.6));found=False
    for condition in conditions:
        data=sorted([r for r in rows if r["condition"]==condition and r["window"]=="n8192"
                     and r["observable"]=="CV" and bm_number(r) is not None],key=bm_number)
        if not data:continue
        found=True;x=np.array([bm_number(r) for r in data]);y=np.array([r["coverage"] for r in data])
        ci=np.array([r["coverage_ci95"] for r in data]).T
        ax.errorbar(x,y,yerr=np.vstack([y-ci[0],ci[1]-y]),fmt="o-",capsize=3,
                    color=COLORS[condition],label=LABELS[condition],lw=1.6,markersize=4.5)
    if not found:plt.close(fig);return
    ax.axhline(.95,color="#9A423A",ls="--",lw=1.3,label="Nominal 95%")
    ax.set(xscale="log",ylim=(0,1.02),xlabel="Non-overlapping batch size (sweeps)",
           ylabel="Empirical 95% interval coverage",title="Exploration: CV energy coverage versus batch size")
    ax.set_xticks([32,128,512,1024]);ax.set_xticklabels(["32","128","512","1024"])
    ax.legend(ncol=2,loc="lower right",fontsize=9)
    fig.tight_layout(rect=(0,.045,1,1))
    save(fig,"exploratory_coverage_batch",out,"128 paired replicas per condition; exact binomial 95% intervals; 8,192 retained sweeps, exact stationary starts.")

def exploratory_length(rows,out):
    fig,ax=plt.subplots(figsize=(7.7,4.6));palette={32:"#355C7D",128:"#BC7E25",512:"#218C7C",2048:"#8A5C9C"};found=False
    for b,color in palette.items():
        data=sorted([r for r in rows if r["condition"]=="local_stationary_d04" and r["observable"]=="CV"
                     and r["method"]==f"BM{b}"],key=lambda r:r["retained"])
        if not data:continue
        found=True;x=np.array([r["retained"] for r in data]);y=np.array([r["coverage"] for r in data]);ci=np.array([r["coverage_ci95"] for r in data]).T
        ax.errorbar(x,y,yerr=np.vstack([y-ci[0],ci[1]-y]),fmt="o-",capsize=3,color=color,label=f"Batch {b}")
    if not found:plt.close(fig);return
    ax.axhline(.95,color="#9A423A",ls="--",lw=1.3,label="Nominal 95%")
    ax.set(xscale="log",ylim=(0,1.02),xlabel="Retained sweeps per trajectory",ylabel="Empirical 95% interval coverage",
           title="Stationary local Metropolis: CV coverage versus length")
    ax.set_xticks([2048,8192,32768]);ax.set_xticklabels(["2,048","8,192","32,768"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.legend(ncol=2,loc="lower right",fontsize=9);fig.tight_layout(rect=(0,.045,1,1))
    save(fig,"exploratory_coverage_length",out,"Same 128 complete trajectories supply nested windows; exact binomial 95% intervals; methods require at least eight batches.")

def variance_ratio(rows,out):
    methods=["naive_IID","BM32","BM128","BM512","BM1024","Geyer_IMS_sensitivity","OBM512","Bartlett512"]
    labels=["IID","BM\n32","BM\n128","BM\n512","BM\n1024","Geyer\nIMS*","OBM\n512","Bartlett\n512"]
    conditions=["local_stationary_d04","reversible_scan","centroid_refresh","iid_gaussian"]
    fig,ax=plt.subplots(figsize=(10,5));found=False
    offsets=np.linspace(-.18,.18,len(conditions))
    for offset,condition in zip(offsets,conditions):
        selected={r["method"]:r for r in rows if r["condition"]==condition and r["window"]=="n8192" and r["observable"]=="CV"}
        data=[(i,selected[m]) for i,m in enumerate(methods) if m in selected]
        if not data:continue
        found=True;x=np.array([i+offset for i,r in data]);y=np.array([r["estimated_to_empirical_variance_ratio"] for i,r in data]);ci=np.array([r["variance_ratio_bootstrap_ci95"] for i,r in data]).T
        ax.errorbar(x,y,yerr=np.vstack([y-ci[0],ci[1]-y]),fmt="o",capsize=2.5,color=COLORS[condition],label=LABELS[condition],markersize=4)
    if not found:plt.close(fig);return
    ax.axhline(1,color="#9A423A",ls="--",lw=1.2);ax.set_yscale("log")
    ax.set_xticks(range(len(methods)));ax.set_xticklabels(labels)
    ax.set(ylabel="Mean estimated variance / empirical variance",title="Exploration: variance calibration of CV trajectory means")
    ax.legend(ncol=2,fontsize=9,loc="lower right");fig.tight_layout(rect=(0,.07,1,1))
    save(fig,"exploratory_variance_ratios",out,"128 paired replicas; bootstrap 95% intervals; *Geyer IMS is a sensitivity analysis for fixed-order nonreversible sweeps.")

def acf_panels(acfs,out):
    conditions=["local_stationary_d04","reversible_scan","centroid_refresh","iid_gaussian"]
    observables=[("CV","CV energy"),("centroid","Signed centroid"),("centroid_component","Centroid energy: ½ centroid²"),("shape","Shape: mean squared distance from centroid")]
    fig,axs=plt.subplots(2,2,figsize=(10.6,7.1),sharey=True);found=False
    for ax,(observable,title) in zip(axs.flat,observables):
        for condition in conditions:
            key=f"{condition}/n8192"
            if key not in acfs:continue
            found=True;data=acfs[key][observable]
            ax.plot(data["lag"],data["pooled_known_mean_acf"],color=COLORS[condition],label=LABELS[condition],lw=1.4)
        ax.axhline(0,color="#75808B",lw=.8);ax.set_title(title);ax.set(xlim=(0,2048),ylim=(-.16,1.03))
        if observable=="shape":ax.set_xlim(0,256)
    if not found:plt.close(fig);return
    axs[1,0].set_xlabel("Lag (sweeps)");axs[1,1].set_xlabel("Lag (sweeps; shorter scale)")
    axs[0,0].set_ylabel("Pooled ACF, known equilibrium mean");axs[1,0].set_ylabel("Pooled ACF, known equilibrium mean")
    handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc="upper center",ncol=4,bbox_to_anchor=(.5,.973),fontsize=9)
    fig.suptitle("Exploration: mode-specific dependence under stationary starts",y=1.005,fontsize=14)
    fig.tight_layout(rect=(0,.055,1,.945))
    save(fig,"exploratory_known_mean_acf",out,"Pooled known-mean estimates across 128 stationary trajectories; finite-lag diagnostics, not a convergence proof; no ACF confidence bands.")

def decomposition(rows,out):
    selected=[]
    for condition in ["local_stationary_d04","centroid_refresh","iid_gaussian"]:
        candidates=[r for r in rows if r["condition"]==condition and r["window"]=="n8192"]
        if candidates:selected.append(candidates[0])
    if not selected:return
    fig,axs=plt.subplots(1,2,figsize=(9.5,4.8),gridspec_kw={"width_ratios":[2.1,1]})
    x=np.arange(len(selected));w=.27
    axs[0].bar(x-w/2,[r["centroid_component_variance"] for r in selected],w,color="#355C7D",label="Centroid energy")
    axs[0].bar(x+w/2,[r["shape_variance"] for r in selected],w,color="#218C7C",label="Shape")
    axs[0].scatter(x,[r["CV_empirical_variance"] for r in selected],color="#202833",marker="D",s=28,label="Total CV",zorder=4)
    axs[0].set_yscale("log");axs[0].set_ylabel("Between-replica variance of time averages")
    axs[0].set_title("Component and total variances");axs[0].legend(fontsize=9)
    axs[1].bar(x,[r["twice_covariance"] for r in selected],color="#8A5C9C",width=.5)
    axs[1].axhline(0,color="#75808B",lw=.9);axs[1].set_title("Twice component covariance")
    axs[1].ticklabel_format(axis="y",style="sci",scilimits=(0,0))
    labels=[LABELS[r["condition"]].replace(" ","\n") for r in selected]
    for ax in axs:ax.set_xticks(x);ax.set_xticklabels(labels,fontsize=9)
    fig.suptitle("CV energy = constant + centroid energy + shape",fontsize=14)
    fig.tight_layout(rect=(0,.06,1,.94))
    save(fig,"exploratory_variance_decomposition",out,"128 replicas, 8,192 sweeps; descriptive component estimates. Total variance = centroid variance + shape variance + twice covariance.")

def confirmation(rows,decision,out):
    conditions=["local_stationary","centroid_refresh","iid_gaussian"]
    data=[]
    for condition in conditions:
        chosen=[r for r in rows if r["condition"]==condition and r["observable"]=="CV" and r["method"]=="BM512"]
        if chosen:data.append(chosen[0])
    if not data:return
    fig,axs=plt.subplots(1,3,figsize=(11.7,4.6),sharey=True)
    y=np.arange(len(data));labels=[LABELS[r["condition"]] for r in data]
    specifications=[("coverage","coverage_ci95","Coverage",.95,(.92,.98)),
                    ("estimated_to_empirical_variance_ratio","variance_ratio_bootstrap_ci95","Estimated / empirical variance",1,(.75,1.25)),
                    ("bias","bias_ci95","Mean energy bias",0,(-.01,.01))]
    for ax,(key,ci,title,target,gate) in zip(axs,specifications):
        ax.axvspan(*gate,color="#E6EFEA",zorder=0);ax.axvline(target,color="#9A423A",ls="--",lw=1.2)
        for i,row in enumerate(data):
            ax.errorbar(row[key],i,xerr=err(row,key,ci),fmt="o",capsize=4,color=COLORS[row["condition"]],markersize=5)
            if key=="coverage":ax.annotate(f"{row['covered']}/{row['replicas']}",(row[key],i),xytext=(0,10),textcoords="offset points",ha="center",fontsize=8)
        ax.set_title(title);ax.set_yticks(y);ax.grid(axis="y",visible=False)
    axs[0].set_xlim(0,1.015);axs[0].set_yticklabels(labels);axs[0].invert_yaxis()
    axs[1].set_xscale("log");axs[2].ticklabel_format(axis="x",style="sci",scilimits=(-2,2))
    n=data[0]["replicas"]
    fig.suptitle(f"Independent confirmation: CV energy, BM512, {n:,} paired stationary replicas",fontsize=14)
    fig.tight_layout(rect=(0,.07,1,.93))
    caption="8,192 sweeps; exact binomial coverage intervals, bootstrap variance-ratio intervals, and t bias intervals; shading: frozen acceptance ranges."
    save(fig,"confirmation_primary_endpoints",out,caption)

def initialization(rows,out):
    conditions=["local_uniform_d04","local_zero_d04","local_displaced_d04"]
    fig,axs=plt.subplots(1,2,figsize=(9.5,4.7));found=False
    colors=["#355C7D","#218C7C","#BC7E25"]
    for condition,color in zip(conditions,colors):
        data=sorted([r for r in rows if r["condition"]==condition and r["observable"]=="CV" and r["method"]=="BM512"],key=lambda r:r["burn"])
        if not data:continue
        found=True;label=condition.split("_")[1].capitalize();x=np.array([r["burn"] for r in data])
        for ax,key,ci in [(axs[0],"bias","bias_ci95"),(axs[1],"coverage","coverage_ci95")]:
            y=np.array([r[key] for r in data]);bounds=np.array([r[ci] for r in data]).T
            ax.errorbar(x,y,yerr=np.vstack([y-bounds[0],bounds[1]-y]),fmt="o-",capsize=3,color=color,label=label)
    if not found:plt.close(fig);return
    axs[0].axhline(0,color="#9A423A",ls="--",lw=1.2);axs[1].axhline(.95,color="#9A423A",ls="--",lw=1.2)
    axs[0].set(ylabel="CV mean bias against finite-M target",title="Initialization and warmup: bias")
    axs[1].set(ylabel="Empirical interval coverage",title="Initialization and warmup: BM512 coverage",ylim=(0,1.02))
    for ax in axs:ax.set_xlabel("Discarded warmup sweeps");ax.set_xticks([0,1024,4096]);ax.set_xticklabels(["0","1,024","4,096"])
    axs[0].legend(title="Initial path",fontsize=9);fig.tight_layout(rect=(0,.055,1,1))
    save(fig,"exploratory_initialization_warmup",out,"128 paired replicas per start; fixed 8,192 retained sweeps; t bias and exact binomial coverage intervals; overlapping windows are correlated.")

def main():
    parser=argparse.ArgumentParser();parser.add_argument("--root",type=Path,default=ROOT);args=parser.parse_args()
    root=args.root;out=root/"figures";out.mkdir(exist_ok=True)
    explored=root/"exploratory-v1"
    if (explored/"analysis.json").exists():
        rows=read(explored/"analysis.json")
        exploratory_batch(rows,out);exploratory_length(rows,out);variance_ratio(rows,out);initialization(rows,out)
        if (explored/"acf.json").exists():acf_panels(read(explored/"acf.json"),out)
        if (explored/"variance_decomposition.json").exists():decomposition(read(explored/"variance_decomposition.json"),out)
    else:print("Exploratory analysis not complete; skipped.")
    confirmed=root/"confirmation-v1"
    if (confirmed/"analysis.json").exists():
        decision=read(confirmed/"primary_decision.json") if (confirmed/"primary_decision.json").exists() else None
        confirmation(read(confirmed/"analysis.json"),decision,out)
    else:print("Confirmation analysis not complete; skipped.")
    manifest={p.name:{"bytes":p.stat().st_size,"sha256":hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(out.iterdir()) if p.is_file() and p.name!="figure_manifest.json"}
    inputs={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
            for stage in [explored,confirmed] for p in sorted(stage.glob("*.json"))
            if p.name in {"protocol.json","analysis.json","acf.json","variance_decomposition.json","primary_decision.json"}}
    (out/"figure_manifest.json").write_text(json.dumps({"plot_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"matplotlib_version":matplotlib.__version__,"inputs_sha256":inputs,"files":manifest},indent=2)+"\n")

if __name__=="__main__":main()
