"""Post-confirmation descriptions at frozen warning boundaries; no fitting."""
from pathlib import Path
import json
import numpy as np
from rnglab17.endpoints import rate

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'confirmation-v2'


def main():
    analysis = json.loads((RUN / 'forecast_analysis.json').read_text())
    result = {
        'classification': 'post-confirmation descriptive diagnostics, not new hypothesis tests',
        'qualification': (
            'Warnings issued only at cutoff 1408, horizon 512. Incident survivors only. '
            'Lead time is event declaration minus issuance, not time to physical damage. '
            'Overlapping windows and two-crossing event construction create mechanical lag. '
            'No restart benefit is measured.'
        ),
        'splits': {},
    }
    for split in ['test', 'OOD']:
        times, features, expected_labels = [], [], []
        for row in analysis['registry']:
            if row['split'] != split:
                continue
            q = json.loads((RUN / f"forecast_{split}_{row['family']}.json").read_text())
            eligible = np.asarray(q['eligible'], bool)
            times.extend(np.asarray(q['events'])[eligible])
            features.extend(np.asarray(q['features'])[eligible])
            expected_labels.extend(np.asarray(q['labels'])[eligible])
        times, features = np.asarray(times), np.asarray(features)
        with np.load(RUN / f'forecast_predictions_{split}.npz') as z:
            y = z['labels'].astype(bool)
            assert np.array_equal(y, expected_labels)
            assert np.array_equal(y, (times > 1408) & (times <= 1920))
            alerts = {
                'static_source_warning': abs(features[:, 0]) > analysis['static_thresholds']['source'],
                'static_application_warning': abs(features[:, 2]) > analysis['static_thresholds']['application'],
                'combined_risk_above_0.1': z['combined'] > .1,
                'combined_risk_at_least_0.3': z['combined'] >= .3,
            }
        result['splits'][split] = {}
        for name, alert in alerts.items():
            tp = int(np.sum(alert & y))
            fp = int(np.sum(alert & ~y))
            fn = int(np.sum(~alert & y))
            tn = int(np.sum(~alert & ~y))
            lead = times[alert & y] - 1408
            assert tp + fp + fn + tn == len(y)
            result['splits'][split][name] = {
                'n': len(y), 'TP': tp, 'FP': fp, 'FN': fn, 'TN': tn,
                'sensitivity': rate(tp, tp + fn),
                'non_event_false_warning_fraction': rate(fp, fp + tn),
                'missed_event_fraction': rate(fn, tp + fn),
                'event_declaration_lead_sweeps': {
                    'n': len(lead),
                    'median': float(np.median(lead)) if len(lead) else None,
                    'minimum': int(lead.min()) if len(lead) else None,
                    'maximum': int(lead.max()) if len(lead) else None,
                },
            }
    (ROOT / 'WARNING_DIAGNOSTICS.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
