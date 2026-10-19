"""Reproduce the explicit exploratory finite-source and dirty-treatment study."""
from pathlib import Path
import argparse
from rnglab15.applications import physical_campaign,defects,hash_condition
from rnglab15.integrity import write_json
from rnglab15.nist import assess


def generate(root):
    root=Path(root)
    streams={}
    for i in range(4):
        source=root/f'measured-run-01/source-{i:02d}'
        labels=['ch0-baseline_sample0_bit0','ch0-gap_divisor_bit0','ch0-timestamp_bit0'] if i<2 else ['ax','gx','ir_pulses']
        for label in labels:
            for treatment in ['raw','vn','sha256']:
                streams[f's{i}-{label}-{treatment}']=source/f'{label}-treatments/{treatment}.bin'
            assess(source/f'{label}.symbols',root/'nist-blockers'/f's{i}-{label}')
    dirty=root/'physical-defects';dirty.mkdir(exist_ok=False)
    for key,p in list(streams.items()):
        if key.endswith('-raw'):
            # Read at most 1MiB of a finite research stream; this source study is
            # explicitly bounded. Large raw waveform files never take this path.
            if p.stat().st_size>1_048_576:
                raise ValueError('physical defect-study budget exceeded')
            data=p.read_bytes()
            for kind in ['stuck','repeat','periodic_dropout','quantized']:
                out=dirty/f'{key}-{kind}.bin';out.write_bytes(defects(data,kind,1234))
                streams[f'{key}-{kind}']=out
                out2=dirty/f'{key}-{kind}-sha256.bin';out2.write_bytes(hash_condition(out.read_bytes()))
                streams[f'{key}-{kind}-sha256']=out2
    write_json(root/'physical-streams.json',{k:str(p.resolve()) for k,p in streams.items()})
    return physical_campaign(streams,root/'physical-applications-01',max_trials=16)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output_root',type=Path)
    generate(p.parse_args().output_root)
