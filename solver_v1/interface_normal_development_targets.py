"""v17 normal-response development and predeclared independent state jets."""
from dataclasses import replace

import numpy as np

from .interface_even_development_targets import even_development_observations
from .vector_material_calibration import IDEAL_H, UNITS, MaterialObservation


def normal_development_observations(source, previous, states):
    rows, provenance = even_development_observations(source, previous, states)
    for i, o in enumerate(rows):
        if o.role == 'heldout':
            rows[i] = replace(o, role='fit')
            provenance[o.name] = 'previously inspected v16 state; v17 development'

    def append(label, q, role, components):
        v = source.evaluate(q)
        jet = np.r_[v.energy, v.gradient, v.hessian[0], v.hessian[1, 1:], v.hessian[2, 2]]
        for index, field in components:
            target = float(jet[index])
            scale = (float(UNITS.traction_mpa_to_force(250.)) if index == 1 else
                     max(.1*abs(target), .01/float(UNITS.energy_to_surface(1.))) if index == 0 else
                     max(.1*abs(target), .1))
            units = {0: 'eV/cell', 1: 'eV/L0', 4: 'eV/L0^2', 7: 'eV/L0^2'}[index]
            name = f'v17_{label}_{field}'
            rows.append(MaterialObservation(name, target, scale, units, role,
                                            tuple(q), tuple(np.eye(10)[index])))
            provenance[name] = ('declared near-perfect normal development' if role == 'fit'
                                else 'v17 predeclared off-grid validation; excluded from loss')

    for ratio in (.985, 1.015, 1.065, 1.095):
        append(f'near_{ratio:g}', (ratio*IDEAL_H, 0., 0.), 'fit',
               ((1, 'normal_force'), (4, 'Haa')))
    states_to_hold = [('direct023', (IDEAL_H, .23, 0.)),
                      ('direct047', (IDEAL_H, .47, 0.)),
                      ('off_plus', (1.075*IDEAL_H, .27, .095)),
                      ('off_minus', (1.135*IDEAL_H, .41, -.07)),
                      ('shockley067', (1.035*IDEAL_H, .335, .67*np.sqrt(3)/6)),
                      ('compression', (.973*IDEAL_H, .08, -.025)),
                      *[(f'opening{r:g}', (r*IDEAL_H, 0., 0.))
                        for r in (1.045, 1.255, 1.555, 1.855, 2.455, 3.75)]]
    for label, state in states_to_hold:
        append(label, state, 'heldout', ((0, 'energy'), (1, 'normal_force'), (4, 'Haa'), (7, 'Hxx')))
    if len({o.name for o in rows}) != len(rows):
        raise ValueError('duplicate observation keys')
    return rows, provenance
