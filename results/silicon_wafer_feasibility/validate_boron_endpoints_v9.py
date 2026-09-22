"""Compare new stationary endpoints under actual periodic cubic symmetries.

Also compute a chemical-potential-free neutral interstitial association energy.
These fixed-volume ML energies are neither finite-T cluster fractions nor rates.
"""
import argparse
import csv
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment


def periodic_cubic_match(reference, candidate, *, proper_only=False):
    """Minimize atom-assignment distance over cubic rotations and B translations.

    Each candidate translation aligns one B with the first reference B. After
    assignment, subtract the small common least-squares translation. Both proper
    and improper cubic symmetries are allowed and explicitly recorded.
    """
    cell=reference.cell.array
    if not np.allclose(cell,np.eye(3)*cell[0,0],rtol=0,atol=1e-9):
        raise ValueError('cubic cell required')
    if not np.array_equal(np.sort(reference.numbers),np.sort(candidate.numbers)):
        raise ValueError('composition mismatch')
    if not np.allclose(cell,candidate.cell.array,rtol=0,atol=1e-9):
        raise ValueError('cell mismatch')
    length=cell[0,0]
    best=None
    for permutation in itertools.permutations(range(3)):
        for signs in itertools.product([-1,1],repeat=3):
            rotation=np.eye(3)[list(permutation)]*np.asarray(signs)[:,None]
            if proper_only and np.linalg.det(rotation)<0:continue
            rotated=candidate.positions@rotation.T
            for anchor in np.flatnonzero(candidate.numbers==5):
                translation=reference.positions[np.flatnonzero(reference.numbers==5)[0]]-rotated[anchor]
                shifted=rotated+translation
                differences=[];pairs=[]
                for z in np.unique(reference.numbers):
                    i=np.flatnonzero(reference.numbers==z);j=np.flatnonzero(candidate.numbers==z)
                    delta=shifted[j][None,:,:]-reference.positions[i][:,None,:]
                    delta-=np.rint(delta/length)*length
                    rows,cols=linear_sum_assignment(np.sum(delta**2,axis=2))
                    differences.extend(delta[rows,cols]);pairs.extend(zip(i[rows].tolist(),j[cols].tolist()))
                differences=np.asarray(differences)
                drift=differences.mean(axis=0);differences-=drift
                rms=float(np.sqrt(np.mean(np.sum(differences**2,axis=1))))
                if best is None or rms<best['RMS_displacement_A']:
                    best=dict(RMS_displacement_A=rms,maximum_displacement_A=float(np.max(np.linalg.norm(differences,axis=1))),
                              rotation=rotation.tolist(),determinant=int(round(np.linalg.det(rotation))),
                              translation_A=(translation-drift).tolist(),assignment_reference_to_candidate=pairs)
    return best


def main(root):
    from ase.io import read
    initial=root/'mace_boron';follow=root/'mace_boron_follow'
    reference=read(initial/'3boron_int_cluster_confi_1_final.extxyz')
    # Recognize independently constructed symmetry-equivalent coordinates,
    # including atom permutation and a periodic wrap, without mixing species.
    control=reference[::-1]
    rotation=np.array([[0.,0.,1.],[1.,0.,0.],[0.,1.,0.]])
    control.positions=control.positions@rotation.T+np.array([1.137,-.381,2.019])
    control.wrap()
    identity_control=periodic_cubic_match(reference,control)
    if identity_control['maximum_displacement_A']>1e-10:raise AssertionError('known symmetry match failed')
    matches=[]
    for path in sorted(follow.glob('*_final.extxyz')):
        name=path.name.removesuffix('_final.extxyz')
        match=periodic_cubic_match(reference,read(path))
        proper_match=periodic_cubic_match(reference,read(path),proper_only=True)
        h=json.loads((root/'mace_boron_follow_hessian'/(name+'_hessian.json')).read_text())
        r=json.loads((follow/(name+'_result.json')).read_text())
        matches.append(dict(structure=name,reference='source-seeded relaxed 3Bi configuration 1',
            symmetry_match=match,proper_rotation_only_match=proper_match,
            symmetry_scope='unloaded cubic cell; do not merge orientations or state multiplicities under directional loading',
            force_max_eV_A=h['max_force_eV_A'],
            min_curvature_eV_A2=h['Cartesian_projected_min_eigenvalue_eV_A2'],
            fixed_cell_minimum_verified=h['Cartesian_projected_min_eigenvalue_eV_A2']>0 and h['max_force_eV_A']<1e-4,
            energy_change_from_saddle_eV=r['energy_change_from_saddle_eV']))
    baseline=json.loads((initial/'run_summary.json').read_text())
    # Find the explicitly stored fixed source-lattice per-atom reference.
    size_records=[]
    for directory in ['mace_boron_216host','mace_boron_512host']:
        size_records.extend(json.loads((root/directory/'summary.json').read_text())['cases'])
    ebulk=size_records[0]['pure_host_energy_per_atom_eV']
    for name,nb in [('b_silicon_64_interstitial',1),('2boron_int_cluster',2)]:
        case=next(r for r in baseline['cases'] if r['structure']==name)
        size_records.append(dict(structure=name,host_atoms=64,B_atoms=nb,
            energy_eV=case['final']['energy_eV'],pure_host_energy_per_atom_eV=ebulk,
            chemical_B_cm3=nb/(10.862**3)*1e24,converged=case['converged']))
    association=[]
    for count in sorted(set(r['host_atoms'] for r in size_records)):
        one=next(r for r in size_records if r['host_atoms']==count and r['B_atoms']==1)
        two=next(r for r in size_records if r['host_atoms']==count and r['B_atoms']==2)
        if not one['converged'] or not two['converged']:raise AssertionError('unfinished force relaxation')
        if one['pure_host_energy_per_atom_eV']!=two['pure_host_energy_per_atom_eV']:
            raise AssertionError('bulk reference mismatch')
        pure=count*one['pure_host_energy_per_atom_eV']
        association.append(dict(host_atoms_per_cell=count,
            one_interstitial_chemical_B_cm3=one['chemical_B_cm3'],
            two_interstitial_chemical_B_cm3=two['chemical_B_cm3'],
            one_interstitial_excess_eV=one['energy_eV']-pure,
            two_interstitial_excess_eV=two['energy_eV']-pure,
            association_energy_eV=2*one['energy_eV']-two['energy_eV']-pure,
            convention='2 E(Si_N+B) - E(Si_N+2B) - E(Si_N); positive favors the paired endpoint',
            full_Hessian_checked=count==64,
            scope='neutral MACE fixed-volume periodic endpoint energy; not association barrier'))
    with (root/'boron_cell_size_association.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(association[0]));w.writeheader();w.writerows(association)
    output=dict(corrected_endpoints=matches,association=association,
        known_rotation_translation_permutation_control=identity_control,
        binding_difference_216_to_512_eV=association[-1]['association_energy_eV']-association[-2]['association_energy_eV'],
        large_cell_Hessians_checked=False,infinite_cell_convergence_certified=False,
        finite_T_cluster_population_calibrated=False,chemical_potential_of_B_required_for_this_difference=False,
        material_calibrated=False,kinetic_calibrated=False)
    (root/'boron_endpoint_validation.json').write_text(json.dumps(output,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(output,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True)
    main(p.parse_args().root)
