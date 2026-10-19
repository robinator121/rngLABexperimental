from __future__ import annotations
import argparse
from pathlib import Path
from .integrity import inventory,write_json
from .pipeline import run
from .applications import simulation_campaign,physical_campaign


def main():
    p=argparse.ArgumentParser(description="RNGLab measured physical-source research")
    sub=p.add_subparsers(dest="command",required=True)
    inv=sub.add_parser("inventory"); inv.add_argument("root",type=Path);inv.add_argument("output",type=Path)
    r=sub.add_parser("run");r.add_argument("config",type=Path);r.add_argument("output",type=Path);r.add_argument("--vault",required=True,type=Path)
    s=sub.add_parser("simulate");s.add_argument("output",type=Path);s.add_argument("--trials",default=32,type=int);s.add_argument("--reference",type=Path)
    a=sub.add_parser("applications");a.add_argument("streams",type=Path);a.add_argument("output",type=Path);a.add_argument("--max-trials",default=16,type=int)
    args=p.parse_args()
    if args.command=="inventory":
        write_json(args.output,{"files":inventory(args.root)})
    elif args.command=="run":
        run(args.config,args.output,args.vault)
    elif args.command=="applications":
        import json
        streams=json.loads(args.streams.read_text())
        physical_campaign({k:Path(v) for k,v in streams.items()},args.output,args.max_trials)
    else:
        simulation_campaign(args.output,args.trials,reference_path=args.reference)


if __name__=="__main__":
    main()
