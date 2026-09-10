"""Exact interior-force restriction of a declared frozen atomic configuration.

The smaller free set only removes unnecessary gradient components. Exterior
coordinates copy the ORIGINAL full core; they are not relaxed or replaced by
continuum positions inside the source free disk. All affected site energies
remain in the same per-atom assembler. This is source target evaluation only.
"""
import numpy as np

from .current_core_coefficient_basis import current_core_coefficients
from .current_material_rows import CurrentMaterialScrewCore


class FrozenCoreForceTarget:
    def __init__(self,core,field,*,radius):
        if not np.isfinite(radius) or radius<=0 or radius>core.free_radius:
            raise ValueError('interior target radius must fit the source free disk')
        core.full_field(field)
        self.source_core=core
        self.radius=float(radius)
        self.displacements={tuple(i):u.copy() for i,u in zip(core.indices[core.free_ids],field)}
        gradients=core.evaluate(field)['gradient']
        self.gradients={tuple(i):g.copy() for i,g in zip(core.indices[core.free_ids],gradients)}

    def coefficient_matrix(self,model,*,ring=7,tolerance=2e-12):
        reference=self.source_core
        if not (np.isclose(model.geometry.b,reference.rows.b,rtol=0,atol=1e-13)
                and np.isclose(model.h,reference.rows.h,rtol=0,atol=1e-13)):
            raise ValueError('source and candidate must use the same frozen atomic geometry')
        core=CurrentMaterialScrewCore(model,reference.far_field,free_radius=self.radius,
            ring=ring,tolerance=tolerance,shear_traction=reference.shear_traction,
            burgers_sign=reference.burgers_sign)
        for index,logical in enumerate(core.indices):
            if tuple(logical) in self.displacements:
                core.boundary[index]=self.displacements[tuple(logical)]
        field=np.array([self.displacements[tuple(i)] for i in core.indices[core.free_ids]])
        target=np.array([self.gradients[tuple(i)] for i in core.indices[core.free_ids]])
        columns=current_core_coefficients(core,field)
        return dict(design=columns['gradient'].reshape(-1,10),target=target.ravel(),
            logical_rows=core.indices[core.free_ids],free_radius=self.radius,
            source_free_radius=reference.free_radius,energy_sites=len(core.energy_ids),
            ring=ring,source_positions_frozen=True,states_relaxed=False)
