"""Digitize a few declared experimental flow points, NOT a new material fit.

Uses the original embedded figure's pixels; .cache copy fetched separately.
The wide strain axis cannot resolve the paper's 0.2% CRSS intercept: that
quantity stays unavailable. Literature points are held out from energy fitting.
"""
from pathlib import Path
import hashlib
import numpy as np
from PIL import Image

from .fetch_strength_reference import ROOT, CACHE
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def trace_pixel_y(strip, *, offset=500):
    """Resolve one connected dark-blue trace; reject occlusion/ambiguity.

    JPEG2000 ringing can leave isolated blue pixels on the other curve. Keep
    connected vertical components spanning at least four pixels, consistent
    with the inspected source stroke width. Never interpolate a hidden trace.
    """
    strip = np.asarray(strip, dtype=int)
    selected=(strip[:,:,2]-strip[:,:,0]>45)&(strip[:,:,2]-strip[:,:,1]>45)
    selected &= (strip[:,:,2]<160)&(strip[:,:,0]<60)&(strip[:,:,1]<60)
    ys=np.where(selected)[0]+offset
    unique=np.unique(ys)
    groups=np.split(unique,np.where(np.diff(unique)>1)[0]+1)
    candidates=[g for g in groups if 4<=len(g)<=16]
    if len(candidates)!=1:
        return None, 0, 'trace occluded or ambiguous; no point inferred'
    kept=ys[np.isin(ys,candidates[0])]
    return float(np.median(kept)), len(kept), 'resolved digitized flow point; NOT yield'


def main():
    figure=CACHE/'article_page6_0.png'
    rgb=np.asarray(Image.open(figure).convert('RGB')).astype(int)
    if rgb.shape!=(2480,3508,3):
        raise ValueError('unexpected source figure dimensions; pixel calibration invalid')
    # Original embedded Figure2b, independently inspected axis pixels.
    x0,x1=2053.,2893.; y0,y20=1031.,191.
    rows=[]
    for g in (.1,.2,.5,.8):
        x=int(round(x0+g*(x1-x0)))
        strip=rgb[500:1010,x-2:x+3]
        # The dark-blue D=103um curve is distinct from bright-blue14.1um,
        # cyan18um, gray6.7um, black axes and the dashed historical curve.
        y, count, status = trace_pixel_y(strip)
        stress=None if y is None else (y0-y)*20/(y0-y20)
        rows.append(dict(source='Krebs et al. 2017 doi:10.1038/nmat4911 Fig2b',
            specimen='99.99% Al as-cast single-crystal wire, diameter103um',
            temperature='room temperature; numerical K not supplied',
            orientation='single-slip; exact loading vector not digitized from stereogram',
            metric='resolved_shear_stress_at_plastic_shear_strain',plastic_shear_strain=g,
            stress_mpa=stress,diameter_um=103.,
            strain_digitization_halfwidth=3/(x1-x0),stress_digitization_halfwidth_mpa=8*20/(y0-y20),
            pixel_x=x,pixel_y=y,selected_pixel_count=count,digitization_status=status,
            experimental_uncertainty='not supplied; pixel band is NOT specimen scatter',used_in_energy_fit=False))
    out=ROOT/'results/fcc111_active_interface/material_strength_v10/experimental_benchmark'
    write_csv(out/'digitized_flow_points.csv',rows)
    save_json(out/'provenance.json',dict(
        doi='10.1038/nmat4911',
        title='Cast aluminium single crystals cross the threshold from bulk to size-dependent stochastic plasticity',
        authors='J. Krebs; S. I. Rao; S. Verheyden; C. Miko; R. Goodall; W. A. Curtin; A. Mortensen',
        year=2017, source_url='https://eprints.whiterose.ac.uk/id/eprint/126662/1/Article%20Text.pdf',
        source_pdf_sha256=hashlib.sha256((CACHE/'article.pdf').read_bytes()).hexdigest(),
        figure_sha256=hashlib.sha256(figure.read_bytes()).hexdigest(),source_pdf_page=6,figure='2b',
        axis_pixels=dict(x_strain_zero=x0,x_strain_one=x1,y_stress_zero=y0,y_stress_20MPa=y20),
        extraction='dark-blue RGB inequalities; connected 4-to-16-row stroke in five-pixel x strip; median y',
        missing_points='gamma=.2 is occluded by the cyan trace and retained as unavailable, not interpolated',
        pixel_uncertainty_not_measurement_uncertainty=True,
        displacement_rate_m_s=300e-9,
        critical_resolved_shear_at_0_2_percent=None,
        unavailable_reason='0.002 strain is below the declared digitization halfwidth on this wide-axis figure',
        exact_specimen_loading_vector=None,source_geometry_population=None,
        numerical_strength_prediction_available=False,material_parameters_fitted_to_experiment=False,
        physical_mobility_calibrated=False))
    print(rows,flush=True)


if __name__=='__main__':
    main()
