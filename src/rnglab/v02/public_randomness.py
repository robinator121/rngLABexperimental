"""Download bounded public data only; authenticate historical drand rounds before use."""
import argparse,datetime,hashlib,json,urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

CHAIN='8990e7a9aaed2ffed73dbd7092123d6f289930540d7651336225dc172e51b2ce'
PUBLIC_KEY='868f005eb8e6e4ca0a47c8a77ceaa5309a47978a7c71bc5cce96366b5d7a569937c529eeda66c7293784a9402801af31'
BASE='https://api.drand.sh/'+CHAIN


def verify(beacon):
    from py_ecc.bls import G2Basic
    signature=bytes.fromhex(beacon['signature'])
    if len(signature)!=96 or len(bytes.fromhex(beacon['randomness']))!=32:raise ValueError('invalid beacon length')
    message=hashlib.sha256(bytes.fromhex(beacon['previous_signature'])+int(beacon['round']).to_bytes(8,'big')).digest()
    if not G2Basic.Verify(bytes.fromhex(PUBLIC_KEY),message,signature):raise ValueError('invalid BLS beacon signature')
    if hashlib.sha256(signature).hexdigest()!=beacon['randomness']:raise ValueError('incorrect beacon randomness digest')
    return bytes.fromhex(beacon['randomness'])


def fetch(round_number):
    url=BASE+'/public/'+str(round_number)
    with urllib.request.urlopen(url,timeout=30) as response:
        raw=response.read(20001)
    if len(raw)>20000:raise ValueError('oversized beacon response')
    record=json.loads(raw)
    if record['round']!=round_number:raise ValueError('wrong beacon round')
    return record


def capture(output,first_round,count):
    if count<1 or count>1000:raise ValueError('capture count outside 1..1000')
    output=Path(output)
    if (output/'drand-rounds.json').exists():raise ValueError('capture already exists; choose a new directory')
    output.mkdir(parents=True,exist_ok=True)
    with ThreadPoolExecutor(max_workers=4) as pool:
        records=list(pool.map(fetch,range(first_round,first_round+count)))
    data=bytearray()
    for index,record in enumerate(records):
        data.extend(verify(record))
        if index and record['previous_signature']!=records[index-1]['signature']:raise ValueError('broken chain continuity')
        if (index+1)%20==0:print(f'{index+1}/{count} beacon signatures verified',flush=True)
    (output/'drand-rounds.json').write_text(json.dumps(records,indent=2))
    (output/'drand-seeds.bin').write_bytes(data)
    provenance={'source':BASE,'first_round':first_round,'last_round':first_round+count-1,'rounds':count,
                'public_key':PUBLIC_KEY,'scheme':'pedersen-bls-chained','verification':'py_ecc G2Basic: every BLS signature, sha256(signature), round match and adjacent chain link',
                'trust_anchor':'public key obtained via HTTPS from drand /info and pinned in this capture tool; not an independent out-of-band audit',
                'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest(),
                'classification':'public distributed beacon; published values are known, not secret entropy','use':'one distinct historical round per trial as ChaCha20 initialization; no claim of independent physical measurements'}
    (output/'drand-provenance.json').write_text(json.dumps(provenance,indent=2))
    return provenance


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--first-round',type=int,required=True);p.add_argument('--count',type=int,default=100)
    a=p.parse_args();print(json.dumps(capture(a.output,a.first_round,a.count),indent=2))

if __name__=='__main__':main()
