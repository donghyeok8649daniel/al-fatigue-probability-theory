"""Read-only force audit of an incomplete research checkpoint, not a PASS."""
import argparse
import csv
import json
from pathlib import Path
import time

import numpy as np

from .current_material_rows import CurrentMaterialScrewCore
from .isolated_screw_core import ScrewFarField
from .run_current_material_core import ROOT,EV_J,load_current_material


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('case',type=Path)
    args=parser.parse_args();began=time.perf_counter()
    meta=json.loads((args.case/'metadata.json').read_bytes())
    model,tensor,binding=load_current_material(ROOT/meta['parameter_source'])
    if binding['parameter_sha256']!=meta['parameter_sha256']:
        raise ValueError('checkpoint material binding differs')
    far=ScrewFarField(tensor,model.geometry.b,tuple(meta['center_over_L0']))
    core=CurrentMaterialScrewCore(model,far,free_radius=meta['radius_over_L0'],ring=meta['ring'],
        tolerance=meta['reciprocal_tolerance'],burgers_sign=meta['burgers_sign'],
        shear_traction=meta['shear_traction_MPa']*1e6*binding['length_scale_m']**3/EV_J)
    with (args.case/'checkpoint.csv').open(encoding='utf8',newline='') as stream:
        rows=list(csv.DictReader(stream))
    by_index={(int(r['j']),int(r['l'])):np.array([float(r[k]) for k in ('ux','uy','uz')]) for r in rows}
    field=np.array([by_index[tuple(i)] for i in core.indices[core.free_ids]])
    result=core.evaluate(field)
    print(json.dumps(dict(energy=result['energy'],maximum_gradient=float(np.max(abs(result['gradient']))),
        declared_force_tolerance=meta['force_tolerance_eV_L0'],Morse_check_performed=False,
        completed_run_certified=False,inspection_seconds=time.perf_counter()-began)))


if __name__=='__main__':main()
