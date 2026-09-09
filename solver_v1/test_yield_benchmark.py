import csv
import json
from pathlib import Path

import numpy as np
import pytest

from .run_yield_benchmark import marker_center, normalized_axes


def test_digitization_axis_unit_endpoints():
    np.testing.assert_array_equal(normalized_axes(2053,2140.5),[0,0])
    np.testing.assert_allclose(normalized_axes(2893,1300.5),[6e-5,6e-4],rtol=2e-16)


def test_marker_centroid_and_missing_marker_refusal():
    image=np.full((40,40,3),255,dtype=np.uint8)
    image[10:30,10:30]=[255,0,0]
    x,y,count=marker_center(image,(5,5,35,35))
    assert (x,y,count)==(19.5,19.5,400)
    image[:]=255
    with pytest.raises(ValueError): marker_center(image,(5,5,35,35))


def test_real_yield_criterion_not_old_large_strain_flow_or_guessed_MPa():
    folder=Path(__file__).resolve().parents[1]/'results/fcc111_active_interface/yield_bridge_v11/experimental_benchmark'
    with (folder/'yield_reference_points.csv').open(newline='') as f: rows=list(csv.DictReader(f))
    assert len(rows)==7
    assert all(float(r['plastic_shear_criterion'])==.002 for r in rows)
    assert all(r['CRSS_MPa']=='' and r['normalizing_G_Pa']=='' for r in rows)
    assert all(r['used_in_energy_fit']=='False' for r in rows)
    provenance=json.loads((folder/'provenance.json').read_text())
    assert provenance['supplement_md5']==provenance['supplement_repository_md5']
    assert not provenance['experimental_strength_used_in_fit']
