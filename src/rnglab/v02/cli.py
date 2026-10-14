import argparse,json
from pathlib import Path
from rnglab.v02.runner import run,summarize


def main():
    p=argparse.ArgumentParser(description='V0.2 physical Monte Carlo research runner')
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=1);p.add_argument('--resume',action='store_true')
    p.add_argument('--shard-index',type=int,default=0);p.add_argument('--shard-count',type=int,default=1)
    p.add_argument('--max-jobs',type=int);p.add_argument('--summarize-only',action='store_true')
    a=p.parse_args();config=json.loads(a.config.read_text())
    if a.summarize_only:
        saved=json.loads((a.output/'manifest.json').read_text())['identity']['config']
        if saved!=config:p.error('summary configuration differs from the saved run')
        summarize(a.output,config)
    else:run(config,a.output,a.workers,a.resume,a.shard_index,a.shard_count,a.max_jobs)

if __name__=='__main__':main()
