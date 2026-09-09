import numpy as np
import pytest

from .public_core_reference import load_core_reference,mev_per_angstrom3_to_mpa


def test_core_pressure_conversion_keeps_mev_distinct_from_ev():
    np.testing.assert_allclose(mev_per_angstrom3_to_mpa([1,.02,1.6]),
        [160.2176634,3.204353268,256.34826144],rtol=2e-15)
    assert mev_per_angstrom3_to_mpa(-1)==-mev_per_angstrom3_to_mpa(1)
    with pytest.raises(ValueError):mev_per_angstrom3_to_mpa(np.nan)


def test_published_core_reference_is_not_yield_or_cell_mobility():
    data=load_core_reference()
    assert not data['production_calibration'] and not data['experimental_yield']
    assert data['collective_coordinate_mobility'] is None
    assert len(data['rows'])==8
    assert data['EAM_potential']=='Ercolessi-Adams (not Mishin Al99)'
    assert data['PN_burgers_angstrom']!=3.94/np.sqrt(2)
