"""Independent real-space sinusoidal energy check at cubic q=(1/2,1/2,0).

Validation ONLY: the candidate generator remains infinite LJ/Bessel. Each
atom sums its own environment before F or tensor invariants. Four exact FCC
phase classes suffice here; no arbitrary simulation box or mass is introduced.
"""
import argparse
from pathlib import Path
import time
from types import SimpleNamespace

import numpy as np

from .angular_environment_reference import traceless_second, traceless_third
from .coordination_screening import CoordinationScreenedBulk, screening_factor
from .periodic_plane_covariance import SourceEAMBlochHessian
from .run_current_material_core import ROOT, load_current_material
from .run_source_core_reference import load_source_material
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json
from .static_bulk_stability import neighbors_in_plane_frame
from .tail_constrained_material import normalized_amplitude
from .vector_material_calibration import LENGTH_M


Q_CUBIC = np.array([.5, .5, 0.])


def phase_differences(geometry, positions):
    q = geometry.plane_basis_in_stacked_cubic_axes() @ (
        Q_CUBIC*2*np.pi/geometry.lattice_constant)
    phase = positions @ q
    if np.max(abs(phase/(np.pi/2)-np.rint(phase/(np.pi/2)))) > 2e-11:
        raise ArithmeticError('neighbors do not belong to the declared four FCC phase classes')
    theta = np.arange(4)*np.pi/2
    return np.cos(theta[:, None]+phase)-np.cos(theta[:, None])


def direct_candidate(model, bulk, positions, differences, amplitude, polarization):
    c = model.coefficients
    energies = []
    for difference in differences:
        R = positions+amplitude*difference[:, None]*polarization
        r = np.linalg.norm(R, axis=1)
        weights = [normalized_amplitude(k)*np.exp(-k*r) for k in bulk.decays]
        x = weights[0].sum()
        w1 = normalized_amplitude(bulk.rank1_decay)*np.exp(-bulk.rank1_decay*r)
        Q1 = w1 @ R
        Q3 = traceless_third(np.einsum('n,ni,nj,nk->ijk', weights[1], R, R, R))
        odd2 = traceless_second(np.einsum('n,ni,nj->ij', weights[1], R, R))
        Q2 = traceless_second(np.einsum('n,ni,nj->ij', weights[2], R, R))
        QE = (odd2-bulk.gauge['eta']*Q2)/bulk.gauge['normalization']
        I3, I2, IE = np.sum(Q3*Q3), np.sum(Q2*Q2), np.sum(QE*QE)
        g = float(screening_factor(x, model.screened_shape[5], law='power')[0])
        energy = (.5*np.sum(c[0]*r**-12-c[1]*r**-6)
            -c[2]*np.sqrt(x)+c[3]*(x-1)+c[4]*(x-1)**2
            +c[5]*I3+c[6]*g*(Q1@Q1)
            +c[7]*I2/(1+model.alpha*I2/model.even_reference_curvature)
            +c[8]*IE/(1+model.alpha*IE)+c[9]*I3**2)
        energies.append(float(energy))
    return float(np.mean(energies))


def direct_source(source, positions, differences, amplitude, polarization):
    energies = []
    for difference in differences:
        R = positions+amplitude*difference[:, None]*polarization
        r = np.linalg.norm(R, axis=1)
        inside = r < source.r[-1]
        rho = np.sum(source._rho(r[inside]))
        energies.append(.5*np.sum(source._rphi(r[inside])/r[inside])+source.F(rho))
    return float(np.mean(energies))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--amplitudes', type=float, nargs='+', default=[.002, .001, .0005])
    args = p.parse_args()
    if args.out.exists():
        raise FileExistsError('fresh independent energy validation required')
    began = time.perf_counter()
    path = ROOT/'results/core_interface_compatibility_v23/radial_full_validation/joint_LS/research_candidate_snapshot.json'
    model, _, binding = load_current_material(path)
    if model.quartic.saturation != 0:
        raise ValueError('this independent validation declares the existing unsaturated I3^2 law')
    context, _, source_binding = load_source_material()
    source = context.source; length_angstrom = LENGTH_M/1e-10
    operator = SourceEAMBlochHessian(source)
    target = operator.evaluate(Q_CUBIC*2*np.pi/source.geometry.lattice_constant)*length_angstrom**2
    bulk = CoordinationScreenedBulk(model.screened_shape, radius=16., law='power')
    columns, tails = bulk.evaluate(Q_CUBIC)
    analytic = np.einsum('c,cij->ij', model.coefficients, columns)
    bound = float(abs(model.coefficients) @ tails)
    eigenvalues, eigenvectors = np.linalg.eigh(analytic-target)
    polarizations = dict(screw=np.array([1., 0., 0.]),
        largest_matrix_discrepancy=eigenvectors[:, np.argmax(abs(eigenvalues))])
    amplitudes = tuple(args.amplitudes)
    if any(not np.isfinite(h) or h <= 0 for h in amplitudes):
        raise ValueError('strictly positive finite displacement amplitudes required')
    results = []
    for name, polarization in polarizations.items():
        for radius in (12., 16.):
            # Keep the same declared lattice shells throughout each amplitude
            # check. This direct truncation is assessed against the analytic
            # harmonic tail and a separate radius change, never made canonical.
            R = neighbors_in_plane_frame(SimpleNamespace(geometry=model.geometry, a0=model.h), radius)
            phases = phase_differences(model.geometry, R)
            zero = direct_candidate(model, bulk, R, phases, 0., polarization)
            for amplitude in amplitudes:
                hi = direct_candidate(model, bulk, R, phases, amplitude, polarization)
                lo = direct_candidate(model, bulk, R, phases, -amplitude, polarization)
                value = 2*(hi+lo-2*zero)/amplitude**2
                expected = float(polarization @ analytic @ polarization)
                results.append(dict(model='v23_radial_joint', polarization=name, radius_L0=radius,
                    amplitude_L0=amplitude, direct_curvature=value, analytic_curvature=expected,
                    absolute_error=abs(value-expected), relative_error=abs(value/expected-1),
                    analytic_operator_tail_bound=bound, energy_units='eV/atom', curvature_units='eV/L0^2'))
        R = neighbors_in_plane_frame(SimpleNamespace(geometry=source.geometry, a0=source.h),
            source.r[-1]+2*max(amplitudes)*length_angstrom+1e-5)
        phases = phase_differences(source.geometry, R)
        zero = direct_source(source, R, phases, 0., polarization)
        for amplitude in amplitudes:
            h = amplitude*length_angstrom
            hi = direct_source(source, R, phases, h, polarization)
            lo = direct_source(source, R, phases, -h, polarization)
            value = 2*(hi+lo-2*zero)/amplitude**2
            expected = float(polarization @ target @ polarization)
            results.append(dict(model='Al99_target_only', polarization=name, radius_L0=None,
                amplitude_L0=amplitude, direct_curvature=value, analytic_curvature=expected,
                absolute_error=abs(value-expected), relative_error=abs(value/expected-1),
                analytic_operator_tail_bound=0., energy_units='eV/atom', curvature_units='eV/L0^2'))
    write_csv(args.out/'energy_curvature.csv', results)
    save_json(args.out/'definition.json', dict(q_fractional_cubic=Q_CUBIC, exact_phase_classes=4,
        energy_expansion='mean delta E=amplitude^2*(e.H.e)/4+O(amplitude^4)',
        polarizations=polarizations, candidate=binding, source=source_binding,
        amplitudes_L0=amplitudes, analytic_operator_radius_L0=16.,
        analytic_tail_is_not_a_nonlinear_direct_energy_error_bound=True,
        cutoff_is_validation_only=True, candidate_production_energy='infinite LJ/Bessel',
        physical_Hz=False, fitting_performed=False, held_out_certification=False))
    save_json(args.out/'completion.json', dict(completed=True, actual_energy_curvatures=len(results),
        elapsed_seconds=time.perf_counter()-began, material_accepted=False,
        infinite_nonlinear_energy_validated_by_this_check=False,
        nonlinear_core_cause_proved=False))
    print('Independent phase-class energy checks:', len(results), flush=True)


if __name__ == '__main__':
    main()
