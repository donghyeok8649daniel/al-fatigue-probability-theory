"""Traceable static core scales, not an imposed yield or mobility law."""
import argparse
import json
from pathlib import Path

import numpy as np

from .aluminum_calibration import EV_J
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def mev_per_angstrom3_to_mpa(value):
    """1 meV / Angstrom^3 = 160.2176634 MPa; neither Hz nor per-cell force."""
    value=np.asarray(value,float)
    if np.any(~np.isfinite(value)):raise ValueError('finite energy-density stress required')
    return value*(1e-3*EV_J/1e-30/1e6)


def load_core_reference(path=None):
    path=Path(path) if path else Path(__file__).with_name('data')/'aluminum_core_validation.json'
    data=json.loads(path.read_bytes())
    if data['status']!='external_static_validation_only_not_fitted' or data['production_calibration']:
        raise ValueError('reference must not claim a production fit')
    return data


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('preserve previous reference audit')
    data=load_core_reference()
    rows=[dict(**row,peierls_MPa=float(mev_per_angstrom3_to_mpa(row['peierls'])),
        temperature_scope=data['condition'],method=data['method'],doi=data['doi'],
        source_location=data['source_location'],used_in_fit=False) for row in data['rows']]
    write_csv(args.out/'published_core_scales.csv',rows)
    save_json(args.out/'scope.json',data)
    print('published static core units audited; no yield or mobility fitted')


if __name__=='__main__':main()
