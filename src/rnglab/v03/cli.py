import argparse,json
from pathlib import Path
from rnglab.v03.runner import run


def main():
    p=argparse.ArgumentParser(description='RNGLAB V0.3 local campaign')
    p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--resume',action='store_true');p.add_argument('--max-jobs',type=int)
    p.add_argument('--shard-index',type=int,default=0);p.add_argument('--shard-count',type=int,default=1)
    a=p.parse_args();run(json.loads(a.config.read_text()),a.output,a.workers,a.resume,a.max_jobs,a.shard_index,a.shard_count)

if __name__=='__main__':main()
