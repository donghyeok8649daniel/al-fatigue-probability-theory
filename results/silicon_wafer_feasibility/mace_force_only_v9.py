"""Optional version-bound force evaluation without unused stress derivatives.

Same loaded neutral MACE energy; independent ASE comparisons are required by
the caller. No model/weight modifications. This is not production code.
"""
import numpy as np


def force_only(calc,atoms):
    if (calc.num_models!=1 or calc.model_type!='MACE' or calc.use_compile
            or calc.energy_units_to_eV!=1 or calc.length_units_to_A!=1):
        raise ValueError('one uncompiled MACE model with eV/Angstrom units required')
    batch=calc._atoms_to_batch(atoms)
    out=calc.models[0](calc._clone_batch(batch).to_dict(),training=False,
        compute_force=True,compute_virials=False,compute_stress=False)
    energy=float(out['energy'].detach().cpu().reshape(-1)[0])
    force=out['forces'].detach().cpu().numpy().copy()
    if force.shape!=atoms.positions.shape or not np.isfinite(energy) or not np.all(np.isfinite(force)):
        raise ValueError('invalid force-only result')
    return energy,force


def force_calculator(model):
    """An optional ASE wrapper that exposes only the properties computed."""
    from ase.calculators.calculator import Calculator,all_changes
    class NeutralForceCalculator(Calculator):
        implemented_properties=['energy','free_energy','forces']
        def __init__(self):super().__init__();self.evaluations=0
        def calculate(self,atoms=None,properties=None,system_changes=all_changes):
            super().calculate(atoms,properties,system_changes)
            energy,force=force_only(model,atoms);self.evaluations+=1
            self.results=dict(energy=energy,free_energy=energy,forces=force)
    return NeutralForceCalculator()
