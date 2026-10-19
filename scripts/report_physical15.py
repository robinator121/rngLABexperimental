"""Rebuild aggregate plots/report from completed local experiment artifacts."""
from __future__ import annotations
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import statistics
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def generate(output):
    output=Path(output)
    reports=[json.loads((output/f'measured-run-01/source-{i:02d}/report.json').read_text()) for i in range(4)]
    simulation=json.loads((output/'synthetic-campaign-01/results.json').read_text())
    physical=json.loads((output/'physical-applications-01/results.json').read_text())
    inventory=json.loads((output/'inventory.json').read_text())['files']
    resource=json.loads((output/'resource-usage.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for i in range(2):
        profile=list(csv.DictReader((output/f'measured-run-01/source-{i:02d}/ch0-sample_profile.csv').open()))
        x=np.array([int(r['sample_index']) for r in profile]);y=np.array([float(r['mean']) for r in profile])
        axes[0,i].plot(x,y,color='#006d77');axes[0,i].set(title=f'Export {i+1}: ensemble waveform mean',xlabel='Sample index (time scale unknown)',ylabel='Code')
        c=reports[i]['channels']['0'];blocks=c['baseline_512_event_blocks']
        axes[1,i].plot(np.arange(len(blocks))*512,[b['mean'] for b in blocks],color='#b5651d')
        axes[1,i].set(title='Exploratory baseline trend',xlabel='Event index',ylabel='Block mean of event baseline median')
    fig.suptitle('Measured waveforms: triggered/clipped acquisitions, unverified physical identity')
    fig.savefig(output/'waveform-diagnostics.png',dpi=160);plt.close(fig)
    labels=['baseline sample-0 LSB','raw timestamp LSB','raw gap LSB','gap / 80 LSB','peak LSB']
    keys=['baseline_sample0_bit0','timestamp_bit0','gap_bit0','gap_divisor_bit0','peak_bit0']
    fig,axes=plt.subplots(1,2,figsize=(12,4),layout='constrained')
    for i in range(2):
        values=[reports[i]['symbol_streams']['ch0-'+k]['screen']['ones_fraction'] for k in keys]
        axes[i].barh(labels,values,color=['#006d77','#b23a48','#b23a48','#83c5be','#e29578'])
        axes[i].axvline(0.5,color='black',linestyle='--',linewidth=1);axes[i].set_xlim(0,1)
        axes[i].set(title=f'Export {i+1}',xlabel='Fraction of ones (balance is not entropy)')
    fig.savefig(output/'source-bit-diagnostics.png',dpi=160);plt.close(fig)
    arms=['os','os_hash','stuck','bias','repeat','quantized','drift','periodic_dropout'];models=['pi','vmc','ising','pimc']
    matrix=[]
    for model in models:
        matrix.append([statistics.mean(r['tolerance_failure'] for r in simulation['rows'] if r['cohort']=='confirmatory' and r['model']==model and r['treatment']==arm) for arm in arms])
    fig,ax=plt.subplots(figsize=(12,4),layout='constrained');im=ax.imshow(matrix,vmin=0,vmax=1,cmap='YlOrRd',aspect='auto')
    ax.set_xticks(range(len(arms)),arms,rotation=20);ax.set_yticks(range(4),models)
    for i in range(4):
        for j in range(len(arms)):ax.text(j,i,f'{matrix[i][j]:.0%}',ha='center',va='center',color=('white' if matrix[i][j]>0.6 else 'black'))
    ax.set_title('Held-out synthetic cohort: finite-horizon tolerance failures (32 pairs per arm)')
    fig.colorbar(im,ax=ax,label='Failure fraction');fig.savefig(output/'application-sensitivity.png',dpi=160);plt.close(fig)
    lines=['# RNGLab v0.15 physical-source development results','',
           'Completed 19 October 2026. Local research results; **no validated secure TRNG output**.', '',
           '## Delivered implementation','',
           'Isolated branch `feat/v015-physical-sources`, based on the clean v0.14 commit `7f6f3b6`. '
           'The actual datasets are in `/home/robi/rnglab-data`; `~/rnglab` does not exist. '
           'Existing source modules and the 120-test baseline were preserved. This release adds raw snapshots, '
           'strict streaming ingest, measured diagnostics, separate symbol/treatment provenance, entropy screens, '
           'continuous health diagnostics, an official-tool adapter, and finite application comparisons.', '',
           'The physical acquisition pipeline integrates directly with the v0.14 platform architecture and measured data manifests.', '',
           '## Input inventory','',
           '| File | Bytes | Mode | Verified interpretation / limit |','|---|---:|---|---|']
    descriptions={'capture.py':'Host serial-to-sensor-CSV program; not a source stream',
                  'baseline_20261019_172913.csv':'Sensor registers, 2,633 rows; local firmware match is provisional',
                  'baseline_20261019T153146Z.csv':'Host-wrapped sensor registers, 3,179 rows; resets/missing record retained',
                  'meas_data_2026-10-19 14.05.19.txt':'117,962 records × (2 provisional metadata + 1,023 waveform codes)',
                  'meas_data_2026-10-19 14.30.18.txt':'103,336 records × (2 provisional metadata + 64 waveform codes)'}
    for r in inventory:lines.append(f"| {r['path']} | {r['size_bytes']:,} | {r['mode']} | {descriptions[r['path']]} |")
    lines+=['','All original raw-file SHA-256 values are in [inventory.json](inventory.json). '
            'Read-only content-addressed snapshots are in the private work evidence store, outside the source commit. '
            'Both original and snapshot hashes were checked after analysis. Execute permission on the text exports '
            'is only a filesystem attribute. No ESP32 RNG or ADC capture was found in this directory.', '',
            '## Measured limitations and useful channels','',
            '| Observation | Export 1 | Export 2 |','|---|---:|---:|']
    for label,key,fmt in [('Waveform codes per event','waveform_samples','d'),('Events reaching 65,535','clipped_event_fraction','.2%'),('Pooled pre-window lag-1 correlation','pooled_baseline_lag1_correlation','.4f')]:
        values=[format(r['channels']['0'][key],fmt) for r in reports[:2]]
        lines.append(f'| {label} | {values[0]} | {values[1]} |')
    for label,key,fmt in [('Gap GCD, unknown raw units','positive_gap_gcd_raw_units','d'),('Gap coefficient of variation','gap_cv','.4f'),('Consecutive positive-gap correlation','lag1_positive_gap_correlation','.4f')]:
        values=[format(r['channels']['0']['timing'][key],fmt) for r in reports[:2]]
        lines.append(f'| {label} | {values[0]} | {values[1]} |')
    lines+=['','The detector clock scale, waveform sampling rate, actual trigger/holdoff settings, sensor identity '
            'and radiation/muon classification remain unknown. No detector rate in seconds or dead-time estimate '
            'is reported. The second waveform is a different record length, not an independent validation of the first. '
            'All rows passed the explicitly configured integer/count/range schema; this does not verify physical semantics.', '',
            'Sensor timing is approximately 19.99 readings/s **conditional on the local firmware micros() interpretation**; '
            'the actual flashed-firmware identity is unverified. The host-wrapped file has two backward device-time steps. '
            'Accelerometer x/y/z codes are all even, so their LSB streams are constant. The gyroscope registers vary. '
            'Pressure, temperature, humidity and IR age are entirely NA; IR pulse counts are entirely zero. '
            'Those absent/constant channels cannot be relabeled entropy sources. The records describe motion-sensor '
            'registers, not muon events or an ESP32 HWRNG.', '',
            'Priority candidate research: free-running untriggered baseline noise under documented bandwidth/ADC '
            'conditions, and separately documented event gaps. Current sample-0 ADC LSB and gap/80 LSB are near balanced, '
            'but a finite predictor screen cannot validate unpredictability. Clipped peak LSB, raw timestamp LSB, '
            'raw gap LSB and accelerometer LSBs are poor direct choices. The timestamp/gap LSBs are identically zero. '
            'Separate channels and correlated pulse features receive no additive entropy credit.', '',
            '![Waveform and baseline diagnostics](waveform-diagnostics.png)', '',
            '![Bit diagnostics](source-bit-diagnostics.png)', '',
            '## Entropy budget and conditioning','',
            '**Validated source entropy: 0 credited bits. Validated output entropy and secure bit rate: unassigned.** '
            'IID-only marginal bounds are explicitly conditional on an unestablished IID assumption. '
            'Held-out majority/context predictors report empirical predictability, not non-IID lower bounds. '
            'The conservative credited non-IID lower bound stays zero. Every defined one-per-event stream has fewer '
            'than 1M sequential symbols; no restart data or official NIST executable is available. '
            'The machine-readable [NIST blockers](nist-blockers/) record why external assessment was not run.', '',
            'Continuous raw RCT/APT thresholds use an assumed H=0.5 bit/symbol and exact binomial APT cutoffs. '
            'These are exploratory alarms, not calibrated source acceptance criteria. Timestamp/gap/peak and '
            'constant sensor streams trigger failures. Raw, VN and fixed-block SHA-256 outputs retain separate hashes, '
            'input/output lengths, incomplete-tail counts and zero entropy credit. The budget-checked Toeplitz '
            'reference is tested but is not used to claim secure extraction without validated conditional entropy '
            'and independent seed evidence. Hashing stuck data can hide its appearance; it cannot create physical entropy.', '',
            'See [NIST SP 800-90B and errata](https://csrc.nist.gov/pubs/sp/800/90/b/final) and '
            '[official byte-symbol tools](https://github.com/usnistgov/SP800-90B_EntropyAssessment/blob/master/README.md). '
            'Future ESP32 comparison needs chip/firmware and entropy-source status: '
            '[Espressif documents the operating prerequisites](https://docs.espressif.com/projects/esp-idf/en/latest/esp32/api-reference/system/random.html).', '',
            '## Application outcomes','',
            f"The controlled campaign contains {len(simulation['rows']):,} histories: 32 exploratory and 32 confirmatory paired identities "
            'per model, eight arms, and pi/VMC/Ising/PIMC. Exact OS bytes are archived for replay. '
            'Primary tolerances, finite horizons and auxiliary seed namespaces were saved before outcomes. '
            'Confirmation uses paired discordant-failure tests and Holm correction over 28 model/treatment hypotheses. '
            'The 8-bit uniform grid and short horizons apply equally to every arm. Ising uses an exact finite L=4 reference; '
            'PIMC uses the exact finite-bead harmonic reference. VMC evaluates spatial x² as well as its identically '
            'constant optimum local energy. This is a source-defect sensitivity experiment, not a convergence study.', '',
            '| Defect / treatment | pi excess failures | VMC x² | Ising energy | PIMC energy |',
            '|---|---:|---:|---:|---:|']
    table={(r['model'],r['treatment']):r for r in simulation['confirmatory']}
    for arm in arms[1:]:
        vals=[]
        for m in models:
            r=table[(m,arm)];vals.append(f"{r['excess_failure_fraction']:+.1%}"+(' *' if r['reject_fwer_0_05'] else ''))
        lines.append('| '+arm+' | '+' | '.join(vals)+' |')
    lines+=['','An asterisk marks Holm-adjusted p ≤ 0.05 in the held-out synthetic cohort. '
            'Null findings are retained and do not establish equivalence. Broad confidence intervals and limited '
            'power follow from only 32 pairs; false-negative risk remains. Within-chain inherited standard errors '
            'are not used for these inferential comparisons.', '',
            '![Application sensitivity](application-sensitivity.png)', '',
            f"Finite physical replay produced {len(physical['rows']):,} short application histories. "
            f"{len(physical['blockers'])} requested stream arms had fewer than two disjoint four-model trials and were explicitly blocked. "
            'The arms include raw/VN/SHA, targeted stuck/repeat/dropout/quantization defects, and SHA after defects. '
            'Physical segments use paired auxiliary seeds and retained offsets; each stream is consumed without expansion or replay. '
            'Different output lengths produce unequal trial counts. These already-observed records are exploratory, '
            'not held-out physical confirmation; no significance or cryptographic quality is claimed.', '',
            'Full per-history outcomes, source byte ranges, hashes and limitations are in '
            '[synthetic results](synthetic-campaign-01/results.json) and [physical replay results](physical-applications-01/results.json).', '',
            '## Reproducibility and remaining work','',
            f"The complete 642.6 MB detector corpus plus both sensor files was scanned in {resource['elapsed_seconds']:.1f} seconds "
            f"with peak Linux resident memory {resource['max_rss_kib_linux']/1024:.1f} MiB in the recorded process. "
            'The raw store uses disk space approximately equal to the originals. Logs and code/input/output hashes '
            'are retained. The source package excludes raw datasets and raw-store objects. No sensitive biometric '
            'or genetic input was found in the explicitly selected acquisition files or published.', '',
            'Blockers are concrete: unknown detector clock/trigger provenance; correlated/clipped waveforms; '
            'sensor missing channels and resets; no GM/verified ESP32 captures; no restart data, validated entropy '
            'or official NIST executable; insufficient finite sensor bits for many application arms; and no future '
            'physical confirmation cohort. Follow the included GM roadmap for the week of 26 October 2026. '
            'The inherited RUS/forecasting and production DRBG stages need separate calibration/entropy evidence '
            'before physical reliability or reseed-security claims.', '',
            'The standalone prior acquisition package remains separate with its existing license. '
            'No finance work, services, original datasets or other working trees were altered. No remote publication was performed.']
    index=lines.index('## Reproducibility and remaining work')
    physical_table=['## Selected measured-stream replay outcomes','',
                    'Export 1 only. Percentages are descriptive tolerance-failure fractions. '
                    'The small, unequal trial counts do not support ranking source quality or claiming conditioning success.', '',
                    '| Stream / treatment | Trials per model | pi | VMC x² | Ising | PIMC |',
                    '|---|---:|---:|---:|---:|---:|']
    selected=[('baseline_sample0_bit0','raw'),('baseline_sample0_bit0','vn'),('baseline_sample0_bit0','sha256'),
              ('gap_divisor_bit0','raw'),('gap_divisor_bit0','vn'),('gap_divisor_bit0','sha256'),
              ('timestamp_bit0','raw'),('timestamp_bit0','sha256')]
    for label,treatment in selected:
        key=f's0-ch0-{label}-{treatment}'
        cells=[]
        n=0
        for model in models:
            rows=[r for r in physical['rows'] if r['treatment']==key and r['model']==model]
            n=len(rows)
            cells.append(f"{statistics.mean(r['tolerance_failure'] for r in rows):.1%}" if rows else 'blocked')
        physical_table.append('| '+label+' / '+treatment+' | '+str(n)+' | '+' | '.join(cells)+' |')
    physical_table+=['','The timestamp-bit input is all zero. Fixed-block hashing makes a repeated digest stream '
                     'that passes some short application tolerances; it still has no source entropy. '
                     'This is a direct example of why application outcomes cannot validate a TRNG.','']
    lines[index:index]=physical_table
    (output/'RNGLAB_V015_RESULTS.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('output',type=Path)
    generate(parser.parse_args().output)
