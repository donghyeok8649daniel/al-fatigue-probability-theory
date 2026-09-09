import io
import zipfile

import numpy as np
import pytest

from .public_activation_reference import (
    activation_unit_conversion, cached_workbook, parse_selected_activation,
)


def test_source_kinetic_derivative_units_not_specimen_area_or_clock():
    result=activation_unit_conversion(2988.898369,2988.898369,burgers_m=.286e-9,
        temperature_K=293.,boltzmann_J_K=1.3806488e-23)
    assert result['rate_log_slope_per_MPa']==pytest.approx(17.2845631366,rel=1e-9)
    assert result['apparent_activation_volume_m3']==pytest.approx(
        .286e-9*result['apparent_activation_area_m2'])
    assert result['specimen_correlation_area'] is None
    assert result['cell_mobility'] is None and result['t0_seconds'] is None
    with pytest.raises(ValueError):
        activation_unit_conversion(1,1,burgers_m=0,temperature_K=293)


def test_parse_preserves_source_signed_error_and_processed_stress():
    raw=('# Relax\tEff stress\tEff strain\ta [b²]\tV [b^3]\ta StError\tb²/a\tb²/a StError\n'
         'relax 1\t3.284\t-.00001\t2988.898369\t2988.898369\t295.6\t.000334571\t-.0000331\t\t\n').encode('cp1252')
    row=parse_selected_activation(raw)[0]
    assert row['source_signed_inverse_error']<0
    assert row['source_effective_strain']<0
    assert row['volume_area_identity_error']==0
    with pytest.raises(ValueError):parse_selected_activation(raw.replace(b'V [b^3]',b'volume'))


def test_cached_xlsx_does_not_execute_formulas():
    stream=io.BytesIO()
    ns='http://schemas.openxmlformats.org/spreadsheetml/2006/main'
    with zipfile.ZipFile(stream,'w') as z:
        z.writestr('xl/sharedStrings.xml',f'<sst xmlns="{ns}"><si><t>T</t></si></sst>')
        z.writestr('xl/workbook.xml',f'<workbook xmlns="{ns}" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"><sheets><sheet name="Numbers" r:id="x"/></sheets></workbook>')
        z.writestr('xl/_rels/workbook.xml.rels','<Relationships><Relationship Id="x" Target="worksheets/sheet1.xml"/></Relationships>')
        z.writestr('xl/worksheets/sheet1.xml',f'<worksheet xmlns="{ns}"><sheetData><row><c r="B1" t="s"><v>0</v></c><c r="C1"><f>UNTRUSTED_EXTERNAL()</f><v>293</v></c></row></sheetData></worksheet>')
    result=cached_workbook(stream.getvalue())['Numbers']
    assert result['B1']['value']=='T'
    assert result['C1']['value']==293
    assert result['C1']['cached_formula']=='UNTRUSTED_EXTERNAL()'
