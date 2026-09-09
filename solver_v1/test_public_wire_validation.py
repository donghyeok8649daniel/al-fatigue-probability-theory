import io
import numpy as np
import pytest

from .public_wire_validation import (
    EXPECTED_COLUMNS, maximum_fcc_schmid, parse_relaxation, summarize_trace,
)
from .public_aluminum_kinetics import SequentialRangeReader


def test_fcc_projection_and_experimental_stress_units():
    assert maximum_fcc_schmid([0, 0, 1])['schmid_factor'] == pytest.approx(1/np.sqrt(6))
    assert maximum_fcc_schmid([1, 1, 1])['schmid_factor'] == pytest.approx(np.sqrt(6)/9)
    projection=maximum_fcc_schmid([7, -1, -2])['schmid_factor']
    assert 0 < projection <= .5
    diameter=14.7; force=np.array([1., .99, .98])
    stress=force*1000/(np.pi/4*diameter**2)
    values=np.column_stack([np.arange(3),np.ones(3),force,np.zeros(3),stress,np.zeros(3),stress*projection])
    result=summarize_trace(values,diameter_um=diameter,orientation=[7,-1,-2])
    assert result['max_stress_from_load_discrepancy_MPa'] == 0
    assert result['implied_diameter_median_um']==pytest.approx(diameter)
    assert result['measured_shear_to_axial_ratio'] == pytest.approx(projection)
    assert result['proof_yield_MPa'] is None
    assert result['collective_coordinate_mobility'] is None


def test_raw_duplicates_and_unlabelled_columns_are_exposed_not_repaired():
    lines=['\t'.join(EXPECTED_COLUMNS),'0\t0\t1\t0\t2\t0\t1\t7\t0',
           '0\t0\t1.1\t0\t2.2\t0\t1.1\t7\t0',
           '1\t0\t.9\t0\t1.8\t0\t.9\t7\t0']
    raw=('\r\n'.join(lines)).encode('cp1252')
    q=parse_relaxation(raw)
    assert q.shape==(3,9) and q[1,2]==1.1
    result=summarize_trace(q,diameter_um=14.7,orientation=[7,-1,-2])
    assert result['duplicate_timestamps']==1
    assert result['unnamed_extra_columns']==2
    assert result['unnamed_columns_max_absolute']==7
    with pytest.raises(ValueError,match='units'):
        parse_relaxation(raw.replace(b'EngStress[MPa]',b'EngStress[GPa]'))
    with pytest.raises(ValueError,match='decrease'):
        parse_relaxation(raw+b'\r\n-1\t0\t1\t0\t2\t0\t1\t7\t0')


def test_ordered_http_ranges_without_network(monkeypatch):
    raw=b'abcdefghijklmnopqrstuvwxyz'
    requests=[]
    def fetch(self,start,end):
        requests.append((start,end))
        return raw[start:end+1],len(raw)
    monkeypatch.setattr(SequentialRangeReader,'_fetch',fetch)
    reader=SequentialRangeReader('unused',20,chunk_bytes=5,workers=2)
    try:
        assert reader.read(3)==b'abc'
        assert reader.read(9)==b'defghijkl'
        assert reader.count==12
        assert reader.read(8)==b'mnopqrst'
        with pytest.raises(ValueError,match='budget'):
            reader.read(1)
    finally:
        reader.close()
    assert sorted(set(requests))==[(0,4),(5,9),(10,14),(15,19)]
