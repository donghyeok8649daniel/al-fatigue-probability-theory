"""Run matched-geometry Mishin TARGET cores; no production energy selection."""
from types import SimpleNamespace

from .fcc111_geometry import fcc111_geometry_from_b
from .nonlocal_interface_elasticity import rotate_elastic_tensor
from .reference_eam_targets import MishinRigidFCCReference, DEFAULT_CACHE, SOURCE_URL
from .run_current_material_core import argument_parser, run, ROOT, EV_J
from .source_core_reference import MishinScrewCoreReference, source_elastic_tensor
from .vector_material_calibration import LENGTH_M


def load_source_material(path=DEFAULT_CACHE, *, tolerance=2e-12):
    source = MishinRigidFCCReference(path)
    length_angstrom = LENGTH_M/1e-10
    geometry = fcc111_geometry_from_b(source.geometry.b/length_angstrom)
    context = SimpleNamespace(source=source, geometry=geometry, h=source.h/length_angstrom)
    tensor = source_elastic_tensor(source, LENGTH_M)
    cubic = rotate_elastic_tensor(tensor, geometry.plane_basis_in_stacked_cubic_axes().T)
    factor = EV_J/LENGTH_M**3/1e9
    metadata = dict(parameter_source='published NIST Al99.eam.alloy target-only',
        parameter_sha256=source.sha256, source_url=SOURCE_URL,
        material_accepted=False, energy_terms_omitted=False, length_scale_m=LENGTH_M,
        energy_unit='eV per infinite straight-row repeat b*L0', tensor_unit='eV/L0^3',
        elastic_tensor_same_candidate=True,
        C11_C12_C44_GPa=[float(cubic[0,0,0,0]*factor), float(cubic[0,0,1,1]*factor), float(cubic[0,1,0,1]*factor)],
        source_nominal_lattice_angstrom=source.geometry.lattice_constant,
        source_reference_only=True, reciprocal_tolerance_applicable=False,
        kernel='published finite-support source EAM, target only; not the LJ/Bessel model',
        physical_time=False, physical_Hz=False, actual_yield_calibrated=False, production_changed=False)
    return context, tensor, metadata


def build_source_core(context, far_field, **kwargs):
    return MishinScrewCoreReference(context.source, far_field, length_scale_m=LENGTH_M, **kwargs)


def main():
    parser = argument_parser()
    parser.description = __doc__
    parser.set_defaults(material=DEFAULT_CACHE, ring=4)
    run(parser.parse_args(), material_loader=load_source_material, core_builder=build_source_core)


if __name__ == '__main__':
    main()
