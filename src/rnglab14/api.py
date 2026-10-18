"""Bounded deterministic JSON interface, no server, provider or experiment daemon."""
import dataclasses
import json
import math
import sys
from .models.pimc_harmonic import exact_finite_bead_energy
from .forecast import FixedHorizonPredictor,evaluate_forecast_on_trajectory,summarize_forecasting_performance


def handle_request(request):
    if not isinstance(request,dict) or request.get('schema')!='rnglab.headless.v1':
        raise ValueError('schema rnglab.headless.v1 required')
    if set(request)-{'schema','operation','parameters'}:raise ValueError('unknown request fields')
    p=request.get('parameters',{});operation=request.get('operation')
    if not isinstance(p,dict):raise ValueError('parameters must be an object')
    if operation=='harmonic_reference':
        if set(p)!={'beads','beta'}:raise ValueError('beads and beta required')
        finite=exact_finite_bead_energy(p['beads'],p['beta'])
        result={'finite_bead_energy':finite,'continuum_energy':.5/math.tanh(p['beta']/2),
                'units':'m=hbar=omega=1','reference_scope':'harmonic primitive finite-bead action'}
    elif operation=='forecast_evaluate':
        if set(p)!={'trajectory','threshold','landmarks','horizon'}:raise ValueError('fixed horizon evaluation fields required')
        if not isinstance(p['trajectory'],list) or len(p['trajectory'])>10000 or not isinstance(p['landmarks'],list) or len(p['landmarks'])>100:
            raise ValueError('bounded trajectory/landmark lists required')
        fs=evaluate_forecast_on_trajectory(FixedHorizonPredictor(p['horizon'],p['threshold']),p['trajectory'],p['threshold'],p['landmarks'])
        result={'forecasts':[dataclasses.asdict(f) for f in fs],
                'summary':dataclasses.asdict(summarize_forecasting_performance(fs)) if fs else None,
                'scientific_confirmation':False}
    else:raise ValueError('unsupported deterministic operation')
    return {'schema':'rnglab.headless.result.v1','operation':operation,'result':result}


def main():
    try:
        raw=sys.stdin.buffer.read(1_000_001)
        if len(raw)>1_000_000:raise ValueError('request exceeds1MB')
        response=handle_request(json.loads(raw))
        print(json.dumps(response,allow_nan=False))
    except (ValueError,TypeError,KeyError) as e:
        print(json.dumps({'schema':'rnglab.headless.error.v1','error':str(e)}));return 2
    return 0

if __name__=='__main__':raise SystemExit(main())
