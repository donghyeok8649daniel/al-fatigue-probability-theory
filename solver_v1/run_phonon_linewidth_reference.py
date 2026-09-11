"""Reproducible vector-marker extraction from a primary Al experiment figure.

Requires optional PyMuPDF for reading only. The source page and crop were
visually inspected. Plot resolution is not an experimental confidence interval.
No production mobility is inferred from these finite-wavevector linewidths.
"""
import argparse
import hashlib
from pathlib import Path

import numpy as np

from .run_vector_registry_audit import save_json


PDF_SHA256='07014983720a37e0ee4b908f36510fd988ab049ebdbb41c55de5a5e5021531cf'


def extract(pdf, output):
    import pymupdf
    pdf,output=Path(pdf),Path(output)
    if hashlib.sha256(pdf.read_bytes()).hexdigest()!=PDF_SHA256:
        raise ValueError('visually verified publisher PDF required')
    if output.exists():
        raise FileExistsError('fresh output required')
    document=pymupdf.open(pdf)
    page=document[2]
    drawings=page.get_drawings()
    # Exact vector axes of Fig2(b), bottom-right T(1/2,1/2,1/2).
    # Source page units are PDF points, not inferred data units.
    x0=486.95599365234375;x300=508.0570068359375;x900=550.2630004882812
    y0=623.9819946289062;y1=567.675048828125
    markers=[];dft_markers=[]
    for d in drawings:
        r=d['rect'];color=d['color']
        if (color is not None and np.allclose(color,[.4,.4,.4],atol=1e-6)
                and 488<r.x0<554 and 553<r.y0<625
                and abs(r.width-r.height)<.01 and 4<r.width<5):
            x=(r.x0+r.x1)/2;y=(r.y0+r.y1)/2
            temperature=300+(x-x300)*600/(x900-x300)
            linewidth=(y0-y)/(y0-y1)
            # Source defines tau=1/(pi Gamma), Eq2 uses cycle-frequency nu.
            g=np.pi*linewidth
            markers.append(dict(temperature_K_from_figure=float(temperature),
                linewidth_THz_from_figure=float(linewidth),
                damping_g_per_ps=float(g),envelope_time_ps=float(1/g),
                source_marker_rect_pdf_points=list(r),
                plotting_resolution_linewidth_THz=float(.5/(y0-y1)),
                plotting_resolution_temperature_K=float(.5*600/(x900-x300)),
                experimental_uncertainty_not_reported_by_this_extraction=True))
        fill=d['fill']
        if (fill is not None and np.allclose(fill,[1.,0.,0.],atol=1e-6)
                and 488<r.x0<554 and 553<r.y0<625
                and abs(r.width-r.height)<.01 and 5.5<r.width<5.8):
            x=(r.x0+r.x1)/2;y=(r.y0+r.y1)/2
            dft_markers.append(dict(temperature_K_from_figure=float(300+(x-x300)*600/(x900-x300)),
                linewidth_THz_from_figure=float((y0-y)/(y0-y1)),
                source_marker_rect_pdf_points=list(r),method='DFT MD GGA, not experiment',
                individual_supercell_assignment_not_resolved=True))
    if len(markers)!=4:
        raise ArithmeticError('expected four visually identified experimental markers')
    markers.sort(key=lambda r:r['temperature_K_from_figure'])
    if len(dft_markers)!=2:raise ArithmeticError('expected two visually identified DFT square markers')
    dft_markers.sort(key=lambda r:r['temperature_K_from_figure'])
    output.mkdir(parents=True)
    save_json(output/'reference.json',dict(
        title='Phonon Lifetimes throughout the Brillouin Zone at Elevated Temperatures from Experiment and Ab Initio',
        authors='Glensk, Grabowski, Hickel, Neugebauer, Neuhaus, Hradil, Petry, Leitner',
        year=2019,doi='10.1103/PhysRevLett.123.235501',
        source_url='https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevLett.123.235501/fulltext',
        pdf_sha256=PDF_SHA256,page=3,figure='2(b), bottom-right',
        method='inelastic neutron scattering; digitized publisher vector-marker centers',
        wavevector_cubic_reciprocal_units=[.5,.5,.5],polarization='transverse',
        linewidth_convention='Eq2 cycle-frequency Gamma; g=pi Gamma; envelope tau=1/g',
        axis_calibration_pdf_points=dict(x0=x0,x300=x300,x900=x900,y0=y0,y1=y1),
        records=markers,dft_comparison_records=dft_markers,
        fit_target=False,external_validation_only=True,
        source_material_is_experimental_Al_not_Al99=True,
        production_coordinate_mobility=None))
    print(markers,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('pdf',type=Path);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();extract(a.pdf,a.out)
