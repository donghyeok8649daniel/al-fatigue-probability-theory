"""A few isolated actual 0.2%-plastic-shear CRSS points, not large-strain flow.

Read Figure2d in the original lossless embedded Krebs2017 figure. No invented
shear modulus is used to convert its normalized CRSS/G axis into MPa. Marker
overlaps are not deconvolved. Scientific figures/PDFs stay in ignored cache.
"""
import hashlib
from pathlib import Path

import numpy as np
from PIL import Image

from .fetch_strength_reference import CACHE, ROOT, URLS
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json


def marker_center(image, box):
    """One pre-inspected unoccluded red square, centroid and pixel support."""
    left,top,right,bottom=box
    block=np.asarray(image)[top:bottom,left:right].astype(int)
    mask=(block[:,:,0]>180)&(block[:,:,1]<90)&(block[:,:,2]<90)
    ys,xs=np.where(mask)
    if len(xs)<300 or len(xs)>410 or max(np.ptp(xs),np.ptp(ys))>21:
        raise ValueError('isolated reference marker not recovered; no inferred point')
    return float(np.mean(xs)+left),float(np.mean(ys)+top),int(len(xs))


def normalized_axes(x,y):
    # Bottom/top and left/right spine CENTERS in3508x2480 embedded figure.
    return (x-2053)/840*6e-5, (2140.5-y)/840*6e-4


def main():
    file=CACHE/'article_page6_0.png'; image=np.asarray(Image.open(file).convert('RGB'))
    if image.shape!=(2480,3508,3): raise ValueError('source pixel geometry changed')
    boxes=[(2610,1394,2640,1425),(2637,1476,2666,1504),(2637,1537,2666,1566),
           (2610,1674,2640,1700),(2627,1698,2657,1727),(2309,1844,2339,1874),
           (2245,1916,2274,1945)]
    rows=[]
    for index,box in enumerate(boxes):
        x,y,count=marker_center(image,box); inverse_d,stress=normalized_axes(x,y)
        rows.append(dict(id=f'Krebs2017_Fig2d_isolated_{index+1}',
            observable='critical_resolved_shear_stress_at_0.002_plastic_shear',
            plastic_shear_criterion=.002,inverse_diameter_b_over_D=inverse_d,
            CRSS_over_G=stress,CRSS_MPa=None,normalizing_G_Pa=None,
            inverse_diameter_reading_halfwidth=12/840*6e-5,
            normalized_stress_reading_halfwidth=12/840*6e-4,
            pixel_x=x,pixel_y=y,pixels=count,box=str(box),
            specimen_state='99.99% as-cast Al single crystal wire',
            temperature='room temperature; numeric K not reported',
            source_geometry='unavailable; wire diameter is NOT source length',
            loading='tensile displacement300nm/s; resolved-shear projection per source',
            source='Krebs et al.2017 doi:10.1038/nmat4911 Fig2d',
            used_in_energy_fit=False,reading_band_not_experimental_scatter=True))
    out=ROOT/'results/fcc111_active_interface/yield_bridge_v11/experimental_benchmark'
    write_csv(out/'yield_reference_points.csv',rows)
    save_json(out/'provenance.json',dict(
        doi='10.1038/nmat4911',source_url=URLS['article'],source_pdf_page=6,figure='2d',
        pdf_sha256=hashlib.sha256((CACHE/'article.pdf').read_bytes()).hexdigest(),
        figure_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
        supplement_url=URLS['supplement'],supplement_sha256=hashlib.sha256((CACHE/'supplement.pdf').read_bytes()).hexdigest(),
        supplement_repository_md5='7787b4507b8f25b0da6f17a829b02abb',
        supplement_md5=hashlib.md5((CACHE/'supplement.pdf').read_bytes()).hexdigest(),
        axis_pixels=dict(x_zero=2053,x_6e_5=2893,y_zero=2140.5,y_6e_4=1300.5),
        marker_boxes=boxes,reading_halfwidth_pixels=12,
        selection='seven isolated red squares; legend, occluded/merged markers not inverted',
        marker_center_uncertainty='conservative full half-marker width plus spine reading; not specimen scatter',
        normalizing_shear_modulus_unavailable=True,physical_strength_MPa_not_invented=True,
        physical_time_calibrated=False,experimental_strength_used_in_fit=False,
        supplementary_model_assumptions_not_measurements=dict(
            source_span_L_over_D=1/3,alternative_L_over_D=1/2,
            single_arm_prefactor_range=[.068,.148],page_range=[23,25],
            used_as_our_material_law=False),
        followup_source=dict(doi='10.1016/j.scriptamat.2018.10.009',
            url=URLS['annealed_article'],pdf_sha256=hashlib.sha256((CACHE/'annealed_article.pdf').read_bytes()).hexdigest(),
            condition='99.99%Al annealed500C2h underAr, wire14/23um; not as-cast benchmark',
            adopted_as_our_fit=False)))
    print(rows,flush=True)


if __name__=='__main__': main()
