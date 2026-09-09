import numpy as np

from .report_scaling_validation import atomic_unit_ledger
from .vector_material_calibration import LENGTH_M
from .aluminum_calibration import EV_J


def test_static_unit_ledger_distinguishes_atomic_line_bulk_and_time_quantities():
    rows=atomic_unit_ledger();d={r['quantity']:r for r in rows}
    assert len(rows)==len(d)
    a_lat=np.sqrt(2)*LENGTH_M
    assert np.isclose(d['crystallographic_per_atom_volume']['value'],a_lat**3/4,rtol=1e-14,atol=0)
    force=d['one_MPa_interface_force']['value']
    assert np.isclose(force*EV_J,1e6*d['atomic_interface_cell_area']['value']*LENGTH_M,rtol=1e-14,atol=0)
    assert np.isclose(force*d['one_eV_per_reduced_coordinate_traction']['value'],.001,rtol=1e-14)
    assert d['mobility_time_scale_numerator']['units']=='m^2/J'
    assert all(not r['actual_strength_calibrated'] and not r['physical_PDE_time_calibrated'] for r in rows)
