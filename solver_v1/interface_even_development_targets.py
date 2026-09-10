"""v16: inspected v15 failures are development, never relabeled blind tests."""
from dataclasses import replace

import numpy as np

from .interface_development_targets import development_observations
from .vector_material_calibration import IDEAL_H, UNITS, MaterialObservation


def even_development_observations(source, previous, states):
    rows, provenance = development_observations(source, previous, states)
    for i, o in enumerate(rows):
        if o.role == 'heldout':
            rows[i] = replace(o, role='fit')
            provenance[o.name] = 'v15 previously inspected validation; v16 development'
    # Declared before v16 fits; none selected from candidate errors.
    hold = [('direct011', (IDEAL_H, .11, 0.)),
            ('direct031', (IDEAL_H, .31, 0.)),
            ('raised043', (1.09*IDEAL_H, .43, 0.)),
            ('off_plus', (1.025*IDEAL_H, .23, .065)),
            ('off_minus', (1.16*IDEAL_H, .32, -.12)),
            ('shockley043', (IDEAL_H, .215, .43*np.sqrt(3)/6)),
            *[(f'opening{r:g}', (r*IDEAL_H, 0., 0.))
              for r in (1.21, 1.62, 1.93, 2.71, 4.3, 6.2)]]
    for label, state in hold:
        v = source.evaluate(state)
        jet = np.r_[v.energy, v.gradient, v.hessian[0], v.hessian[1, 1:], v.hessian[2, 2]]
        for index, field in ((0, 'energy'), (1, 'normal_force'), (4, 'Haa')):
            target = float(jet[index])
            scale = (float(UNITS.traction_mpa_to_force(250.)) if index == 1 else
                     max(.1*abs(target), .01/float(UNITS.energy_to_surface(1.))) if index == 0 else
                     max(.1*abs(target), .1))
            units = {0: 'eV/cell', 1: 'eV/L0', 4: 'eV/L0^2'}[index]
            name = f'v16_{label}_{field}'
            rows.append(MaterialObservation(name, target, scale, units, 'heldout',
                                            tuple(state), tuple(np.eye(10)[index])))
            provenance[name] = 'v16 predeclared off-grid source validation; excluded from loss'
    if len({o.name for o in rows}) != len(rows):
        raise ValueError('duplicate observation names')
    return rows, provenance
