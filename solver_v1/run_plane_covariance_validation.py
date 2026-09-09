"""Validate known plane-count/energy units against actual 300K public MD."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from .reference_eam_targets import MishinRigidFCCReference
from .periodic_plane_covariance import source_predicted_covariance,empirical_plane_modes
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('fresh source validation required')
    metadata=json.loads((args.source/'source_metadata.json').read_bytes())
    with np.load(args.source/'plane_coordinates.npz',allow_pickle=False) as data:q=data['coordinates_m']
    reference=MishinRigidFCCReference(args.source/'source_Al99.eam.alloy',lattice_angstrom=metadata['lattice_angstrom'])
    predicted=source_predicted_covariance(reference,repeats=q.shape[1],temperature_K=metadata['temperature_K'])
    measured,local=empirical_plane_modes(q)
    block=[empirical_plane_modes(q[i*len(q)//4:(i+1)*len(q)//4]) for i in range(4)]
    error=np.max(abs(np.array([r[1] for r in block])-local),axis=0)
    mode_error=np.max(abs(np.array([r[0] for r in block])-measured),axis=0)
    rows=[]
    for k in range(1,q.shape[1]):
        for j,name in enumerate(('normal','slip_110','transverse')):
            source=predicted['mode_covariances_m2'][k-1,j,j]
            value=measured[k,j,j].real
            rows.append(dict(mode=k,coordinate=name,harmonic_source_variance_m2=source,
                actual_MD_variance_m2=value,relative_difference=value/source-1,
                empirical_four_block_range_m2=mode_error[k,j,j],
                harmonic_stiffness_J_m2=predicted['hessians_J_m2'][k-1,j,j]))
    write_csv(args.out/'mode_variance.csv',rows)
    write_csv(args.out/'local_variance.csv',[dict(coordinate=name,
        harmonic_source_rms_m=np.sqrt(predicted['local_covariance_m2'][j,j]),
        actual_MD_rms_m=np.sqrt(local[j,j]),variance_relative_difference=local[j,j]/predicted['local_covariance_m2'][j,j]-1,
        empirical_block_variance_range_m2=error[j,j]) for j,name in enumerate(('normal','slip_110','transverse'))])
    save_json(args.out/'normalization.json',dict(completed=True,source_doi=metadata['doi'],
        source_potential_sha256=reference.sha256,
        coordinate_sha256=hashlib.sha256((args.source/'plane_coordinates.npz').read_bytes()).hexdigest(),
        repeats=q.shape[1],atoms_per_plane=predicted['atoms_per_plane'],
        actual_lattice_angstrom=metadata['lattice_angstrom'],temperature_K=metadata['temperature_K'],
        local_MD_covariance_m2=local,harmonic_source_covariance_m2=predicted['local_covariance_m2'],
        empirical_block_variance_range_m2=error,
        finite_time_window_uncertainty='four blocks, not a statistical confidence interval or box-size study',
        independent_regions_assumed=False,empirical_scale_factor_fitted=False,
        static_test_only=True,anharmonic_thermal_renormalization_included=False,
        time_calibration_available=False,production_changed=False))
    print('source harmonic/actual MD normalization comparison complete',flush=True)


if __name__=='__main__':main()
