"""Plot saved training and excluded-state comparisons without any refitting."""
import argparse
import csv
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .run_low_frequency_forcing_v29 import sha


def plot(report, validation):
    report, validation = Path(report), Path(validation)
    outputs = [report/('calibration_comparison.'+ext) for ext in ('png', 'svg')]
    if any(path.exists() for path in outputs):
        raise FileExistsError('fresh plot outputs required')
    selected = json.loads((report/'training_selection.json').read_bytes())
    checked = json.loads((validation/'summary.json').read_bytes())
    if (not selected['completed'] or not checked['completed']
            or sha(selected['selected_report']) != checked['comparisons'][1]['source_report_sha256']):
        raise ValueError('completed validation must refer to the training-selected candidate')
    with (validation/'interface.csv').open(encoding='utf8', newline='') as stream:
        interface = list(csv.DictReader(stream))
    with (validation/'finite_q.csv').open(encoding='utf8', newline='') as stream:
        waves = list(csv.DictReader(stream))
    fig, axes = plt.subplots(2, 2, figsize=(11.8, 8.4), layout='constrained')
    fig.suptitle('Material calibration with exact anchors', fontsize=18)
    colors = ('#84909a', '#1768ac', '#d66521', '#24876d', '#8765b3')
    for row, color in zip(selected['comparisons'][1:], colors[1:]):
        saved = json.loads(Path(row['report']).read_bytes())
        times = np.array([r['elapsed_seconds'] for r in saved['records']])/60
        errors = np.minimum.accumulate([r['profile_eta'] if r['positive_LJ'] else np.inf
                                        for r in saved['records']])
        axes[0, 0].plot(times, errors, '-o', color=color, ms=3, label=row['model'])
    axes[0, 0].set(title='Training progress (vertical axis enlarged)',
                   xlabel='Elapsed time [min]', ylabel='Best eligible training eta')
    axes[0, 0].legend(fontsize=8)
    rows = selected['comparisons']
    eta = [row['training_eta'] for row in rows]
    axes[0, 1].bar(np.arange(len(rows)), eta, color=colors[:len(rows)], width=.6)
    axes[0, 1].set_xticks(np.arange(len(rows)), [r['model'].replace('_', '\n') for r in rows], fontsize=9)
    axes[0, 1].axhline(1., color='#333333', linestyle='--', linewidth=1, label='Required eta <= 1')
    axes[0, 1].set(title='Remaining discrepancy on full scale', ylabel='Maximum normalized residual', ylim=(0, max(eta)*1.17))
    for i, value in enumerate(eta):
        axes[0, 1].text(i, value+.13, f'{value:.4f}', ha='center', fontsize=9)
    axes[0, 1].legend(fontsize=8, loc='upper right')
    chosen_color = colors[[r['model'] for r in rows].index(selected['selected_model'])]
    for label, color, name in (('baseline', colors[0], 'v34 baseline'),
                                ('candidate', chosen_color, 'training-selected candidate')):
        ir = [r for r in interface if r['model'] == label]
        qr = [r for r in waves if r['model'] == label and r['role'] == 'new_excluded']
        axes[1, 0].plot([int(r['state_index']) for r in ir],
                         [100*float(r['relative_H_error']) for r in ir], 'o-', color=color, label=name)
        axes[1, 1].plot(np.arange(len(qr)), [100*float(r['relative_H_error']) for r in qr],
                         'o-', color=color, label=name)
    axes[1, 0].set(title='Ten new interface states', xlabel='State index (coordinates in interface.csv)',
                   ylabel='Full Hessian relative error [%]', ylim=(0, None))
    axes[1, 0].set_xticks(np.arange(10))
    axes[1, 1].set(title='Six new excluded wavevectors', xlabel='Wavevector index (coordinates in finite_q.csv)',
                   ylabel='Full Hessian relative error [%]', ylim=(0, None))
    axes[1, 1].set_xticks(np.arange(6))
    for ax in axes[1]:
        ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(alpha=.18, axis='y')
        ax.set_axisbelow(True)
    fig.supxlabel('Same energy family, targets, discrepancy scales and coefficient signs.\n'
                  'Excluded states do not select the fit. Static source comparison does not calibrate physical time.', fontsize=9)
    fig.savefig(outputs[0], dpi=170)
    stream = io.StringIO()
    with matplotlib.rc_context({'svg.hashsalt': 'anchored_calibration_v35'}):
        fig.savefig(stream, format='svg', metadata={'Date': None})
    outputs[1].write_text('\n'.join(line.rstrip() for line in stream.getvalue().splitlines())+'\n', encoding='utf8')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--report', type=Path, required=True)
    parser.add_argument('--validation', type=Path, required=True)
    args = parser.parse_args()
    plot(args.report, args.validation)
