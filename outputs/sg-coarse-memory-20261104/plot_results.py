from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root = Path(__file__).resolve().parent
moving = json.loads((root / 'MOVING_LIFTING_RESULT.json').read_text())
pair = json.loads((root / 'BOUNDARY_MIXTURE_RESULT.json').read_text())

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
truth = [case['truth'] for case in pair['levels'][0]['cases']]
markov = [level['cases'][1]['final_audit']['cdf_markov'] for level in moving['levels']]
axes[0].bar([0, 1], truth, color=['#4d6a85', '#c76c37'], width=.6)
axes[0].set_xticks([0, 1], ['Gibbs input', 'Boundary mixture'])
axes[0].set_ylim(0, .45)
axes[0].set_ylabel('Initiated raw mass at step 24')
axes[0].set_title('Same initial well masses, distinct joint laws')
for i, value in enumerate(truth):
    axes[0].text(i, value + .012, f'{value:.6g}', ha='center')

for i, level in enumerate(moving['levels']):
    result = level['cases'][1]
    axes[1].plot([i, i], [result['lower'], result['upper']], color='#4d6a85', lw=6, solid_capstyle='round')
axes[1].plot(range(3), markov, 'o--', color='#c76c37', label='Memoryless prediction')
axes[1].axhline(truth[1], color='#263238', lw=1.4, label='Full SG table product')
axes[1].set_xticks(range(3), [str(level['macro_states']) for level in moving['levels']])
axes[1].set_xlabel('Macro states including cemetery')
axes[1].set_ylabel('Initiated raw mass')
axes[1].set_ylim(-.02, 1.04)
axes[1].set_title('Event envelopes retain unresolved-state uncertainty')
axes[1].legend(loc='upper right', frameon=False)
for ax in axes:
    ax.spines[['top', 'right']].set_visible(False)
    ax.grid(axis='y', alpha=.18)
fig.suptitle('Finite-grid SG audit · explicit input laws · model time', fontsize=13)
fig.savefig(root / 'source_memory_cdf.png', dpi=180)
print(root / 'source_memory_cdf.png')
