"""Finite paired dirty/conditioned comparisons, separate from entropy claims."""
from __future__ import annotations
from dataclasses import asdict
import hashlib
from functools import lru_cache
import itertools
import math
from pathlib import Path
import random
import numpy as np
from scipy.stats import binomtest
from rnglab14.sources.base import RandomSource, SourceManifest, StochasticChannels
from rnglab14.models.pi_mc import run_pi_monte_carlo
from rnglab14.models.vmc_harmonic import HarmonicOscillatorVMC
from rnglab14.models.ising_2d import Ising2D
from rnglab14.models.pimc_harmonic import PathIntegralHarmonicOscillator
from .integrity import digest, write_json

SLOT_BYTES = 256
MODELS = ("pi", "vmc", "ising", "pimc")


class FiniteSource(RandomSource):
    """Compatibility adapter with an explicit 8-bit midpoint uniform grid.

    File bytes are never replayed, expanded or reseeded. This grid is for
    research comparisons only; it is not a full-precision uniform generator.
    """
    def __init__(self, payload):
        self.payload, self.position = payload, 0
        super().__init__(SourceManifest("finite-research", "measured-or-simulated-replay", True, False,
                                        metadata={"uniform_mapping":"(byte+0.5)/256"}), None)

    def _read_bytes(self, n: int) -> bytes:
        if self.position + n > len(self.payload):
            raise EOFError("finite physical stream exhausted; no replay or padding")
        data = self.payload[self.position:self.position + n]
        self.position += n
        return data

    def read_bits(self,n):
        if not isinstance(n,int) or n<0:
            raise ValueError("nonnegative integer bit count required")
        data=self.read_bytes((n+7)//8)
        return [(byte>>shift)&1 for byte in data for shift in range(7,-1,-1)][:n]

    def random_float(self):
        return (self.read_bytes(1)[0]+0.5)/256

    def read_int(self, low, high):
        if not 0 < high-low <= 256:
            raise ValueError("bounded integer range must be 1..256")
        size=high-low
        limit=256-(256%size)
        while True:
            x=self.read_bytes(1)[0]
            if x < limit:
                return low+x%size


def seed_for(namespace, cohort, identity, role):
    return int.from_bytes(hashlib.sha256(f"{namespace}/{cohort}/{identity}/{role}".encode()).digest()[:16],"big")


def model_run(model, proposal, seed):
    # Auxiliary seeds are exactly paired across treatments; role namespaces differ.
    acc=FiniteSource(random.Random(seed).randbytes(8192))
    init=FiniteSource(random.Random(seed^0x719AC5).randbytes(8192))
    src=FiniteSource(proposal)
    channels=StochasticChannels(src,acc,init)
    if model == "pi":
        result=asdict(run_pi_monte_carlo(channels,samples=128))
        outcome=result["estimated_pi"]
        failure=abs(outcome-math.pi)>0.30
    elif model == "vmc":
        result=asdict(HarmonicOscillatorVMC(alpha=0.5).run_simulation(channels,steps=256,burn_in=64))
        outcome=result["x2_mean"]
        failure=abs(outcome-0.5)>0.30
    elif model == "ising":
        result=asdict(Ising2D(L=4,temperature=3.0).run_simulation(channels,sweeps=8,burn_in=2))
        outcome=result["energy_mean"]
        failure=abs(outcome-ising_exact_energy(4,3.0))>0.50
    elif model == "pimc":
        result=asdict(PathIntegralHarmonicOscillator(beads=4,temperature=0.5).run_simulation(channels,sweeps=32,burn_in=0))
        outcome=result["virial_energy_mean"]
        failure=abs(outcome-finite_bead_energy(4,2.0))>0.30
    else:
        raise ValueError("unknown model")
    return {"outcome":outcome,"tolerance_failure":failure,"result":result,
            "proposal_bytes_consumed":src.position,"proposal_slot_sha256":hashlib.sha256(proposal).hexdigest(),
            "proposal_consumed_sha256":src.stream_digest,"auxiliary_seed":seed}


def finite_bead_energy(beads,beta):
    """Exact discrete harmonic ring covariance, matching centroid virial observable."""
    dt=beta/beads
    variance=sum(1/(dt+4/dt*math.sin(math.pi*k/beads)**2) for k in range(beads))/beads
    centroid_variance=1/beta
    return 1/(2*beta)+variance-0.5*centroid_variance


@lru_cache(maxsize=8)
def ising_exact_energy(L,T):
    if L != 4:
        raise ValueError("bounded exact reference supports L=4")
    # Enumerate 65,536 spin configurations once; stable Boltzmann weights.
    states=np.arange(1<<(L*L),dtype=np.uint32)
    spins=(2*((states[:,None] >> np.arange(L*L))&1).astype(np.int16)-1).reshape(-1,L,L)
    energies=-(spins*(np.roll(spins,-1,axis=1)+np.roll(spins,-1,axis=2))).sum(axis=(1,2))
    weights=np.exp(-(energies-energies.min())/T)
    return float(np.sum(weights*energies)/np.sum(weights)/(L*L))


def defects(data, kind, seed):
    """Targeted deterministic transforms of the identical reference bytes."""
    rng=random.Random(seed)
    out=bytearray(data)
    if kind == "stuck":
        return bytes(len(data))
    if kind == "bias":
        return bytes(x&0x7f for x in data)
    if kind == "quantized":
        return bytes(x&0xf0 for x in data)
    if kind == "drift":
        return bytes(min(255,int(x+i/max(1,len(data)-1)*96)) for i,x in enumerate(data))
    if kind == "periodic_dropout":
        return bytes(0 if i%64<16 else x for i,x in enumerate(data))
    if kind == "repeat":
        for i in range(1,len(out)):
            if rng.random()<0.8:
                out[i]=out[i-1]
        return bytes(out)
    raise ValueError("unknown defect")


def hash_condition(data):
    return b"".join(hashlib.sha256(b"RNGLab15 research\0"+data[i:i+64]).digest()
                    for i in range(0,len(data)-63,64))


def paired_summary(rows, baseline="os", namespace="confirmatory"):
    indexed={(r["identity"],r["model"],r["treatment"]):r for r in rows if r["cohort"]==namespace}
    treatments=sorted({k[2] for k in indexed}-{baseline})
    tests=[]
    for model, treatment in itertools.product(MODELS,treatments):
        pairs=[]
        for (identity,m,t),value in indexed.items():
            ref=indexed.get((identity,m,baseline))
            if m==model and t==treatment and ref:
                pairs.append((value["tolerance_failure"],ref["tolerance_failure"],value["outcome"]-ref["outcome"]))
        if not pairs:
            continue
        worse=sum(a and not b for a,b,_ in pairs)
        better=sum(b and not a for a,b,_ in pairs)
        discordant=worse+better
        p=float(binomtest(worse,discordant,0.5).pvalue) if discordant else 1.0
        tests.append({"model":model,"treatment":treatment,"pairs":len(pairs),"worse":worse,"better":better,
                      "mean_outcome_difference":sum(x[2] for x in pairs)/len(pairs),
                      "excess_failure_fraction":(worse-better)/len(pairs),"paired_sign_test_p":p})
    # Holm step-down adjusted p-values across every model/treatment comparison.
    previous=0.0
    for rank,index in enumerate(sorted(range(len(tests)),key=lambda i:tests[i]["paired_sign_test_p"])):
        previous=max(previous,min(1.,(len(tests)-rank)*tests[index]["paired_sign_test_p"]))
        tests[index]["holm_adjusted_p"]=previous
        tests[index]["reject_fwer_0_05"]=previous<=0.05
    return tests


def simulation_campaign(out, trials=32, namespace="physical-v015-frozen-01", reference_path=None):
    """OS bytes archived before outcomes. Separate cohort/identity/role seeds."""
    import os
    out=Path(out); out.mkdir(parents=True,exist_ok=False)
    if not 2<=trials<=1024:
        raise ValueError("bounded trial count 2..1024 required")
    protocol={"namespace":namespace,"trials_per_cohort":trials,"cohorts":["exploratory","confirmatory"],
              "slot_bytes":SLOT_BYTES,"models":list(MODELS),"paired_roles":["acceptance","initialization","defect"],
              "primary_endpoint":"paired excess tolerance failure vs archived OS CSPRNG",
              "tolerances":{"pi":0.30,"vmc_x2":0.30,"ising_energy":0.50,"pimc_energy":0.30},
              "holm_family":"all 4 models x all 7 treatments; confirmatory cohort only",
              "qualification":"short-horizon controlled sensitivity test; no convergence or general reliability claim",
              "pimc_reference":finite_bead_energy(8,2.),"ising_reference":ising_exact_energy(4,3.)}
    write_json(out/"protocol.json",protocol)
    # OS entropy is not reproducible from a seed: retain exact bytes for replay.
    required_bytes=2*trials*len(MODELS)*SLOT_BYTES*2
    if reference_path is None:
        fixture=os.urandom(required_bytes)
    else:
        if Path(reference_path).stat().st_size != required_bytes:
            raise ValueError("reference fixture length differs from protocol")
        with Path(reference_path).open("rb") as source:
            fixture=source.read(required_bytes)
    with (out/"os-reference.bin").open("xb") as f:
        f.write(fixture)
    rows=[]
    cursor=0
    for cohort,identity in itertools.product(protocol["cohorts"],range(trials)):
        for model in MODELS:
            base=fixture[cursor:cursor+SLOT_BYTES*2]; cursor+=SLOT_BYTES*2
            seed=seed_for(namespace,cohort,identity,model+"/aux")
            variants={"os":base[:SLOT_BYTES],"os_hash":hash_condition(base)}
            for kind in ("stuck","bias","repeat","quantized","drift","periodic_dropout"):
                variants[kind]=defects(base[:SLOT_BYTES],kind,seed_for(namespace,cohort,identity,model+"/"+kind))
            for name,data in variants.items():
                rows.append({"cohort":cohort,"identity":identity,"model":model,"treatment":name,
                             **model_run(model,data,seed)})
    result={"protocol":protocol,"os_fixture_sha256":digest(out/"os-reference.bin"),"rows":rows,
            "confirmatory":paired_summary(rows),"exploratory":paired_summary(rows,namespace="exploratory"),
            "limitations":["Finite horizon, 8-bit grid and small cohort; power is limited.",
                           "Within-chain standard errors in inherited models are descriptive, not inferential endpoints.",
                           "Source identities are paired auxiliary seeds; finite raw vs hash input budgets differ by 2x."]}
    write_json(out/"results.json",result)
    return result


def physical_campaign(streams, out, max_trials=32):
    """Reuse exact finite file segments, retain offsets and consumption. No padding."""
    if not isinstance(max_trials,int) or not 2<=max_trials<=1024:
        raise ValueError("bounded max_trials 2..1024 required")
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    protocol={"streams":{k:{"path":str(p),"sha256":digest(p),"bytes":Path(p).stat().st_size} for k,p in streams.items()},
              "slot_bytes":SLOT_BYTES,"max_trials":max_trials,"models":list(MODELS),
              "status":"exploratory only; all physical acquisitions are already observed",
              "qualification":"short-horizon finite replay, no fresh entropy or population inference"}
    write_json(out/"protocol.json",protocol)
    results,blockers=[],[]
    for name,path in streams.items():
        count=min(max_trials,Path(path).stat().st_size//(SLOT_BYTES*len(MODELS)))
        if count<2:
            blockers.append({"stream":name,"reason":"fewer than two disjoint four-model trials; no expansion/replay allowed"})
            continue
        with Path(path).open("rb") as f:
            for identity in range(count):
                for model in MODELS:
                    offset=f.tell();data=f.read(SLOT_BYTES)
                    results.append({"treatment":name,"identity":identity,"model":model,"offset":offset,
                                    **model_run(model,data,seed_for("physical-finite-v015","exploratory",identity,model))})
    result={"protocol":protocol,"rows":results,"blockers":blockers,
            "qualification":"unequal available trial counts; no physical confirmatory significance claims"}
    write_json(out/"results.json",result)
    return result
