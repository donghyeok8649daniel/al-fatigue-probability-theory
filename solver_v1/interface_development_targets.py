"""Predeclared v15 interface-development data, separate from blind validation.

The v14 held-out shape failures are now DEVELOPMENT evidence. Reusing them
as a loss while still calling them held out would be leakage. The potential,
units, rounded 0 K bulk targets and reference potential remain unchanged.
No yield, fatigue, kinetic or source-length target is introduced here.
"""
from dataclasses import replace

import numpy as np

from .vector_material_calibration import MaterialObservation, IDEAL_H, UNITS


def development_observations(source, previous, states):
    """Return fixed source-state jets with an explicit development/test split.

    All previously inspected rows become development rows. New off-grid and
    off-path states are kept out of fitting. They are interpolation/extrapolation
    tests of this source surface, NOT independent experimental Al validation.
    The 10% scales are model-discrepancy tolerances, not measurement errors.
    Small/nonzero-force rows use an explicit 250 MPa absolute force scale.
    """
    rows = [replace(o, role='fit') if o.role == 'heldout' else o for o in previous]
    provenance = {o.name: ('previously inspected v14 data; development' if o.role == 'heldout'
                           else 'inherited v14 development') for o in previous}

    def append(name, q, component, role, *, weights=None):
        value = source.evaluate(q)
        jet = np.r_[value.energy, value.gradient, value.hessian[0],
                    value.hessian[1, 1:], value.hessian[2, 2]]
        w = np.eye(10)[component] if weights is None else np.asarray(weights, float)
        target = float(w @ jet)
        if component in (1, 2, 3):
            scale = float(UNITS.traction_mpa_to_force(250.))
            unit = 'eV/L0'
        elif component == 0:
            scale = max(.10*abs(target), .01/float(UNITS.energy_to_surface(1.)))
            unit = 'eV/cell'
        else:
            scale = max(.10*abs(target), .1)
            unit = 'eV/L0^2'
        rows.append(MaterialObservation(name, target, scale, unit, role,
                                        tuple(map(float, q)), tuple(w)))
        provenance[name] = ('new declared development jet' if role == 'fit'
                            else 'new off-grid/off-path validation; excluded from loss')

    # Resolve a/s coupling as well as diagonal stiffness at the measured saddle.
    direction = np.array([.5, np.sqrt(3)/6]); direction /= np.linalg.norm(direction)
    coupling = np.zeros(10); coupling[5:7] = direction
    append('saddle_normal_path_coupling', states['saddle'], 5, 'fit', weights=coupling)
    for ratio in (1.1, 1.3, 1.5, 2., 2.5, 3.):
        append(f'opening_{ratio:g}h_normal_force', [ratio*IDEAL_H, 0., 0.], 1, 'fit')
    for fraction in (.15, .35, .75):
        q = [IDEAL_H, fraction*.5, fraction*np.sqrt(3)/6]
        append(f'shockley_{fraction:g}_energy', q, 0, 'fit')

    # These locations are fixed before any v15 optimization. Do not select a
    # model on their residual and subsequently call the same points blind.
    held_states = [
        ('direct_017', (IDEAL_H, .17, 0.)),
        ('direct_039', (IDEAL_H, .39, 0.)),
        ('raised_direct', (1.065*IDEAL_H, .29, 0.)),
        ('off_path_positive', (1.045*IDEAL_H, .19, .085)),
        ('off_path_negative', (1.12*IDEAL_H, .37, -.095)),
        ('shockley_058', (IDEAL_H, .29, .58*np.sqrt(3)/6)),
    ]
    for label, q in held_states:
        for component, field in ((0, 'energy'), (1, 'normal_force'), (4, 'Haa')):
            append(f'new_{label}_{field}', q, component, 'heldout')
    for ratio in (1.17, 1.42, 1.78, 2.23, 3.4, 5.5):
        for component, field in ((0, 'energy'), (1, 'normal_force')):
            append(f'new_opening_{ratio:g}h_{field}', [ratio*IDEAL_H, 0., 0.],
                   component, 'heldout')
    names = [o.name for o in rows]
    if len(names) != len(set(names)):
        raise ValueError('observation keys must be unique')
    return rows, provenance
