"""Reproduce a scoped physical line-velocity benchmark; no production clock."""
from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys

import numpy as np

from .dislocation_line_kinetics import fit_velocity_response
from .run_low_stress_cyclic_diagnostic import write_csv
from .run_vector_registry_audit import save_json

ROOT = Path(__file__).resolve().parents[1]
SOURCE_SHA = 'b4536b12a161babc565460facb83b4c64a685bced57cf0f48da9f7cb79df44d3'
REFERENCE = Path(__file__).with_name('data')/'gorman1969_velocity_benchmark.json'


def extract_original_page(source: Path, output: Path):
    """Lossless source-page image, not a regenerated/AI-created data figure."""
    if hashlib.sha256(source.read_bytes()).hexdigest() != SOURCE_SHA:
        raise ValueError('source PDF differs from the audited published copy')
    if output.exists():
        raise FileExistsError('do not overwrite an existing inspected image')
    import pymupdf  # optional research-only dependency, never imported by UI
    with pymupdf.open(source) as document:
        images = document[4].get_images()
        if len(images) != 1 or images[0][2:4] != (2490, 3305):
            raise ValueError('unexpected source page image geometry')
        output.parent.mkdir(parents=True, exist_ok=True)
        pymupdf.Pixmap(document, images[0][0]).save(str(output))
    return dict(image_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
                source_sha256=SOURCE_SHA, source_page=5, image_size=[2490, 3305])


def velocity_points(reference):
    """Invert the source's slightly skewed scanned axes, then convert CGS.

    A 1e6 dyn/cm2 tick equals 0.1 MPa, NOT 1 MPa. Source uncertainty is
    not reduced by treating repeated points/markers as independent specimens.
    """
    axes = reference['axes_pixels']
    origin = np.array(axes['origin'], float)
    inverse = np.linalg.inv(np.column_stack([
        np.array(axes[key])-origin for key in ('stress_max', 'velocity_max')]))
    pixel = np.array(reference['centers_pixels'], float)
    scaled = (pixel-origin) @ inverse.T
    factors = np.array([
        reference['stress_max_original']*reference['stress_unit_to_Pa'],
        reference['velocity_max_original']*reference['velocity_unit_to_m_s']])
    values = scaled*factors
    # Includes the full half-marker reading allowance, not a claimed
    # probabilistic error bar. Separate source systematics remain uncorrected.
    error = np.sum(abs(inverse), axis=1)*reference['reading_halfwidth_pixels']*factors
    return values, error


def validate_original_markers(image_file, reference):
    """Recover every isolated ring in the predeclared ROI; never infer overlap."""
    from PIL import Image
    from scipy.ndimage import find_objects, label
    if hashlib.sha256(image_file.read_bytes()).hexdigest() != reference['image_sha256']:
        raise ValueError('source image changed')
    image = np.array(Image.open(image_file).convert('L'))
    if list(image.shape[::-1]) != reference['image_size']:
        raise ValueError('unexpected source pixel geometry')
    components, _ = label(image[2300:2840, 365:720] < 128, np.ones((3, 3)))
    measured = []
    for index, box in enumerate(find_objects(components), 1):
        height, width = [s.stop-s.start for s in box]
        if 7 <= height <= 15 and 7 <= width <= 15:
            measured.append([365+(box[1].start+box[1].stop-1)/2,
                             2300+(box[0].start+box[0].stop-1)/2])
    measured.sort()
    if not np.array_equal(measured, reference['centers_pixels']):
        raise ValueError('source isolated markers do not match recorded digitization')
    return len(measured)


def benchmark(output, *, image_file=None):
    import json
    if output.exists():
        raise FileExistsError('fresh calibration result directory required')
    reference = json.loads(REFERENCE.read_bytes())
    if image_file is not None:
        validate_original_markers(image_file, reference)
    values, reading = velocity_points(reference)
    training = np.ones(len(values), bool)
    training[reference['heldout_zero_based_indices']] = False
    fit = fit_velocity_response(values[:, 0], values[:, 1], training)
    x, y = values[training].T
    dx, dy = reading
    if np.any(x <= dx) or np.any(y <= dy):
        raise ValueError('reading band overlaps zero in the chosen calibration')
    # A conservative ratio bound for a joint worst-case perturbation of
    # every fit coordinate. It does NOT include physical point scatter.
    reading_slope_bounds = [float(np.sum((x-dx)*(y-dy))/np.sum((x+dx)**2)),
                            float(np.sum((x+dx)*(y+dy))/np.sum((x-dx)**2))]
    rows = []
    for i, (stress, velocity) in enumerate(values):
        rows.append(dict(point_id=i, stress_Pa=stress, velocity_m_s=velocity,
            stress_reading_halfwidth_Pa=reading[0], velocity_reading_halfwidth_m_s=reading[1],
            pixel_x=reference['centers_pixels'][i][0], pixel_y=reference['centers_pixels'][i][1],
            role='fit' if training[i] else 'heldout_same_figure',
            prediction_m_s=fit['prediction_m_s'][i], residual_m_s=fit['residual_m_s'][i]))
    write_csv(output/'velocity_points.csv', rows)
    save_json(output/'calibration.json', dict(fit,
        temperature_K=reference['temperature_K'], temperature_accuracy_K=reference['temperature_accuracy_K'],
        source_doi=reference['doi'], source_pdf_sha256=SOURCE_SHA,
        data_sha256=hashlib.sha256(REFERENCE.read_bytes()).hexdigest(),
        source_pixels_replayed=image_file is not None,
        reading_only_slope_bounds=reading_slope_bounds,
        observed_stress_range_Pa=[float(values[:, 0].min()), float(values[:, 0].max())],
        validation_scope='selected isolated leading edge/mixed dislocation markers at 23 C; not all source data',
        full_population_calibrated=False, line_character_resolved=False,
        physical_a_s_mobility_available=False, t0_seconds=None))
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.2, 4), layout='constrained')
    for mask, label_text in ((training, 'fit: selected source markers'), (~training, 'held-out markers')):
        ax.errorbar(values[mask, 0]/1e6, values[mask, 1], xerr=reading[0]/1e6,
            yerr=reading[1], fmt='o', ms=5, capsize=2, label=label_text)
    stress = np.linspace(0, values[:, 0].max(), 100)
    ax.plot(stress/1e6, fit['velocity_per_stress_m_per_Pa_s']*stress, label='zero-intercept velocity fit')
    ax.set(xlabel='Resolved shear [MPa]', ylabel='Measured line velocity [m/s]',
           title='23 C, 99.999% Al: line translation, NOT a/s mobility')
    ax.legend(fontsize=7); ax.grid(alpha=.2)
    fig.savefig(output/'velocity_validation.png', dpi=150)
    fig.savefig(output/'velocity_validation.svg')
    plt.close(fig)
    return fit


if __name__ == '__main__':
    import json
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extract', type=Path, help='audited original PDF')
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--pdf-tools', type=Path)
    parser.add_argument('--image', type=Path, help='optional original-image replay for calibration')
    args = parser.parse_args()
    if args.pdf_tools:
        sys.path.insert(0, str(args.pdf_tools))
    result = (extract_original_page(args.extract, args.out) if args.extract else
              benchmark(args.out, image_file=args.image))
    print(json.dumps(result, default=lambda x: x.tolist()), flush=True)
