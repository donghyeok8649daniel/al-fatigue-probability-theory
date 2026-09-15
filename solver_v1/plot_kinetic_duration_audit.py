"""Plot saved duration/heating diagnostics; no fitting or simulation."""
import argparse
import csv
import io
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plot(report):
    report = Path(report)
    outputs = [report/('duration_heating.'+extension) for extension in ('png', 'svg')]
    if any(path.exists() for path in outputs):
        raise FileExistsError('fresh plot outputs required')
    with (report/'noise_scaling.csv').open(newline='', encoding='utf8') as stream:
        noise = list(csv.DictReader(stream))
    with (report/'precision.csv').open(newline='', encoding='utf8') as stream:
        precision = list(csv.DictReader(stream))
    with (report/'thermal_blocks.csv').open(newline='', encoding='utf8') as stream:
        thermal = list(csv.DictReader(stream))
    seeds = sorted({int(row['seed']) for row in noise})
    if len(seeds) != 2:
        raise ValueError('this campaign figure expects two initializations')
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.3), layout='constrained')
    fig.suptitle('Longer kinetic measurements: noise and energy costs', fontsize=17)
    for axis, title in enumerate(('Normal gap', 'Registry gap')):
        ax = axes[0, axis]
        for seed, color in zip(seeds, ('#1768ac', '#d66521')):
            rows = sorted([r for r in noise if int(r['axis']) == axis and int(r['seed']) == seed],
                          key=lambda r: float(r['duration_ps']))
            durations = np.array([float(r['duration_ps']) for r in rows])
            rms = np.array([float(r['rms_loss_A2_eV']) for r in rows])
            ax.loglog(durations, rms, 'o-', color=color, label=f'Seed {seed}: observed RMS', ms=4)
        candidates = [r for r in precision if int(r['axis']) == axis]
        # Every full-record FDT band/taper/seed choice is retained in the band.
        full_duration = max(float(r['duration_ps']) for r in noise)
        stds = np.array([float(r['predicted_loss_std_A2_eV']) for r in candidates])
        ax.fill_between(durations, stds.min()*np.sqrt(full_duration/durations),
                        stds.max()*np.sqrt(full_duration/durations), color='#777777', alpha=.2,
                        label='Conditional FDT prediction range')
        ax.set(title=title+' | unforced sine quadrature', xlabel='Window duration [ps]',
               ylabel=r'Loss RMS [$\AA^2$/eV]')
        ax.grid(alpha=.18, which='both')
        ax.legend(fontsize=8, loc='upper right')
    for column, seed in enumerate(seeds):
        ax = axes[1, column]
        rows = [r for r in thermal if int(r['seed']) == seed]
        names = list(dict.fromkeys(r['case'] for r in rows))
        for name in names:
            records = [r for r in rows if r['case'] == name]
            axis = records[0]['axis']
            sign = int(records[0]['sign'])
            color = '#666666' if axis == '' else ('#1768ac' if axis == '0' else '#d66521')
            label = 'Unforced' if axis == '' else ('Normal' if axis == '0' else 'Registry')+f' {sign:+d}'
            centers = [(float(r['start_ps'])+float(r['end_ps']))/2 for r in records]
            temperatures = np.array([float(r['mean_temperature_K']) for r in records])
            ax.plot(np.asarray(centers)/1000, temperatures-temperatures[0],
                    color=color, linestyle='--' if sign < 0 else '-', label=label, linewidth=1.6)
        ax.set(title=f'Seed {seed} | measured 100 ps temperature means',
               xlabel='Reference MD time [ns]', ylabel='Change from first analyzed block [K]')
        ax.axhline(0., color='#999999', linewidth=.7)
        ax.grid(alpha=.18)
        ax.legend(fontsize=8, ncol=2, loc='upper left')
    fig.supxlabel('Existing 20 ns / 10 records. Shared planes and windows are not independent samples.\n'
                  'Prediction bands are sensitivity ranges, not confidence intervals. Production Hz remains unavailable.',
                  fontsize=9)
    for path in outputs:
        if path.suffix == '.svg':
            stream = io.StringIO()
            with matplotlib.rc_context({'svg.hashsalt': 'kinetic_duration_audit'}):
                fig.savefig(stream, format='svg', metadata={'Date': None})
            # Matplotlib paths contain trailing spaces; normalize the text export.
            path.write_text('\n'.join(line.rstrip() for line in stream.getvalue().splitlines())+'\n',
                            encoding='utf8')
        else:
            fig.savefig(path, dpi=180)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--report', type=Path, required=True)
    plot(parser.parse_args().report)
