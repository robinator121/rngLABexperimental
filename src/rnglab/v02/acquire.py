"""Record OS CSPRNG bytes once for public, reproducible experimental replay."""
import argparse,datetime,hashlib,json,os,platform
from pathlib import Path


def capture_os(path,byte_count):
    path=Path(path)
    if byte_count<32:raise ValueError('capture at least 32 bytes')
    path.parent.mkdir(parents=True,exist_ok=True)
    h=hashlib.sha256()
    with path.open('xb') as file:
        left=byte_count
        while left:
            block=os.urandom(min(left,1024*1024));file.write(block);h.update(block);left-=len(block)
    record={'source':'Python os.urandom / operating-system CSPRNG','bytes':byte_count,'sha256':h.hexdigest(),
            'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'platform':platform.platform(),
            'classification':'OS CSPRNG output; no raw physical-source claim','purpose':'public experimental corpus, not secret key material',
            'encoding':'packed-msb','exhaustion':'raise EOFError; no cycling'}
    path.with_suffix(path.suffix+'.provenance.json').write_text(json.dumps(record,indent=2))
    return record


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--bytes',type=int,default=16*1024*1024)
    a=p.parse_args();print(json.dumps(capture_os(a.output,a.bytes),indent=2))

if __name__=='__main__':main()
