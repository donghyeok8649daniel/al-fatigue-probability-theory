"""Wait for the existing campaign, then analyze; never launch duplicate MD."""
import argparse
import json
from pathlib import Path
import time
from .run_weak_replica_v31 import analyze
from .run_collective_bridge_v31 import run as bridge
from .run_vector_registry_audit import save_json
from .run_low_frequency_forcing_v29 import sha


def finish(study,out,bridges,watch_dir,timeout_hours=24.):
    study,out,bridges,watch_dir=map(Path,(study,out,bridges,watch_dir))
    if out.exists() or watch_dir.exists():raise FileExistsError('fresh completion outputs required')
    if not 0<timeout_hours<=48:raise ValueError('bounded positive wait required')
    cases=json.loads((study/'cases.json').read_bytes())
    watch_dir.mkdir(parents=True);start=time.monotonic()
    try:
        while True:
            done=[c for c in cases if (study/c['name']/'summary.json').exists()]
            save_json(watch_dir/'status.json',dict(completed=False,stage='waiting_for_existing_MD',
                completed_records=len(done),required_records=len(cases),
                elapsed_seconds=time.monotonic()-start,production_clock_calibrated=False))
            if len(done)==len(cases):break
            if time.monotonic()-start>timeout_hours*3600:raise TimeoutError('MD incomplete; raw states preserved')
            for c in cases:
                log=study/(c['name']+'.log')
                if log.exists() and 'Traceback (most recent call last)' in log.read_text(encoding='utf8'):
                    raise RuntimeError('MD failure in '+c['name']+'; no partial calibration')
            print(f'Existing weak MD: {len(done)}/{len(cases)} complete',flush=True)
            time.sleep(30)
        analyze(study,out)  # checks completion, protocol and source binding
        for seed in sorted({c['seed'] for c in cases}):
            source=study/f'seed{seed}_null';dest=bridges/f'seed{seed}_new'
            if dest.exists():
                existing=json.loads((dest/'summary.json').read_bytes())
                if not existing['completed'] or existing['source_sha256']!=sha(source/'plane_coordinates.npz'):
                    raise ValueError('existing bridge report binding mismatch')
            else:bridge(source,dest)
        save_json(watch_dir/'status.json',dict(completed=True,stage='analysis_complete_not_calibration_approval',
            completed_records=len(cases),elapsed_seconds=time.monotonic()-start,
            production_clock_calibrated=False,material_accepted=False))
        print('MD response and generator analyses complete; physical calibration NOT automatically approved',flush=True)
    except Exception as error:
        save_json(watch_dir/'status.json',dict(completed=False,stage='failed',error=str(error),
            elapsed_seconds=time.monotonic()-start,production_clock_calibrated=False))
        raise


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__)
    for name in ('study','out','bridges','watch-dir'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();finish(a.study,a.out,a.bridges,a.watch_dir)
