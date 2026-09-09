"""Actually audit a downloaded public trajectory, not a hypothetical clock."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .physical_time import BOLTZMANN_J_PER_K
from .public_aluminum_kinetics import correlation_audit
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def audit(directory, output):
    directory, output = Path(directory), Path(output)
    if output.exists():
        raise FileExistsError('fresh kinetic audit output required')
    metadata = json.loads((directory/'source_metadata.json').read_bytes())
    raw = (directory/'plane_coordinates.npz').read_bytes()
    with np.load(directory/'plane_coordinates.npz', allow_pickle=False) as data:
        times, q = data['time_seconds'], data['coordinates_m']
    if len(q) != metadata['frames'] or not np.allclose(
            np.diff(times), metadata['frame_seconds'], rtol=1e-10, atol=0):
        raise ValueError('source trajectory length/time mismatch')
    rows=[]; studies=[]; detailed={}
    # Fixed tests before looking at the measured covariance, not an optimized
    # choice of a convenient fit window. Prefix comparison is statistical/time
    # window sensitivity, not an atomistic box-size convergence certificate.
    plans=[('first_half',q[:len(q)//2],1,4),('full',q,1,4),
           ('eight_blocks',q,1,8),('stride2',q[::2],2,4)]
    for name, values, stride, blocks in plans:
        result=correlation_audit(values,frame_seconds=metadata['frame_seconds']*stride,
            max_lag=31//stride,blocks=blocks)
        detailed[name]=result
        for k,t in enumerate(result['lags_seconds']):
            rows.append(dict(study=name,lag_seconds=t,
                minimum_whitened_correlation_eigenvalue=result['whitened_eigenvalues'][k,0],
                maximum_whitened_correlation_eigenvalue=result['whitened_eigenvalues'][k,-1],
                block_sensitivity=result['block_sensitivity'][k],
                antisymmetric_residual=result['antisymmetric_residual'][k],
                empirical_error_floor=result['empirical_error_floor'][k],
                negative_beyond_floor=bool(result['negative_beyond_empirical_floor'][k])))
        k=int(np.argmin(result['whitened_eigenvalues'][:,0]))
        failures=np.flatnonzero(result['negative_beyond_empirical_floor'])
        studies.append(dict(study=name,frames=len(values),blocks=blocks,stride=stride,
            covariance_rms_coordinate_m=np.sqrt(np.diag(result['covariance_m2'][0])),
            first_resolved_negative_lag_seconds=None if not len(failures) else result['lags_seconds'][failures[0]],
            most_negative=result['whitened_eigenvalues'][k,0],
            lag_at_minimum_seconds=result['lags_seconds'][k],error_at_minimum=result['empirical_error_floor'][k],
            interpretation=result['interpretation']))
    write_csv(output/'correlation_modes.csv',rows)
    save_json(output/'correlation_details.json',detailed)
    save_json(output/'source_metadata.json',dict(metadata,projection_file_sha256=hashlib.sha256(raw).hexdigest()))
    C0=detailed['full']['covariance_m2'][0]
    save_json(output/'kinetic_decision.json',dict(completed=True,studies=studies,
        temperature_K=metadata['temperature_K'],source_doi=metadata['doi'],
        inverse_covariance_curvature_J_per_m2=BOLTZMANN_J_PER_K*metadata['temperature_K']*np.linalg.inv(C0),
        covariance_curvature_scope='local Gaussian PMF diagnostic of these periodic plane-mean coordinates, not rigid-interface W Hessian',
        coordinate_equivalent_to_production=False,
        reversible_harmonic_overdamped_consistent=not any(np.any(r['negative_beyond_empirical_floor']) for r in detailed.values()),
        M_a_phys=None,M_s_phys=None,t0_seconds=None,
        physical_seconds_enabled=False,physical_Hz_enabled=False,
        source_time_is_physical=True,source_time_is_solver_calibration=False,
        next_requirement='matched collective coordinates/PMF and a resolved Markov lag window; defect kinetics requires defect-containing data'))
    return studies


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(audit(args.source,args.out),default=lambda x:x.tolist(),indent=2))
