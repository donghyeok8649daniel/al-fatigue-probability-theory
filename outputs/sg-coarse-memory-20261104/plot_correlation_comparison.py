from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=Path(__file__).parent
r=json.loads((root/'CORRELATION_MEMORY_RESULT.json').read_text())
w=json.loads((root/'CORRELATION_WITNESS_EXACT.json').read_text())
h=json.loads((root/'REACHABLE_LOAD_ORDER_RESULT.json').read_text())
fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
ax=axes[0]
vals=[float(w[k]['decimal_lower']) for k in ['corrected_exact_cdf_plus','corrected_exact_cdf_minus']]
ax.bar(['Preparation +','Preparation -'],vals,color=['#007d72','#d97732'],width=.5)
ax.set_ylabel('First-opening probability (one model cycle)')
ax.set_title('Same marginals and mean energy\nDifferent joint correlations')
for i,val in enumerate(vals):ax.text(i,val+3e-5,f'{val:.7f}',ha='center',fontsize=10)
ax.set_ylim(0,max(vals)*1.18)
ax=axes[1];cycles=[2,4,8]
for key,label,style in [('marginals_energy','Marginals + energy','o-'),('plus_simple_covariance','+ simple covariance','x--'),('plus_one_cycle_event_memory','+ one-cycle response','s-')]:
    vals=[row['summaries'][key]['width'] for row in r['comparisons'][1:]]
    ax.plot(cycles,vals,style,label=label)
noisy=[x['width'] for x in r['noisy_probe_comparisons'] if x['probe_absolute_uncertainty']==1e-4]
ax.plot(cycles,noisy,'^:',label='Response uncertainty +/-1e-4')
ax.set_yscale('log');ax.set_xticks(cycles);ax.set_xlabel('Held-out horizon (same model cycle)')
ax.set_ylabel('Numerical identification interval width');ax.legend(fontsize=8)
ax.set_title('Information retained for future prediction\nNot a physical accuracy or speed factor')
ax=axes[2]
for mode,label in [('original_load_masks','Original load masks'),('common_fixed_mask','Common fixed mask')]:
    rows=[x for x in h['cases'] if x['mode']==mode]
    ax.plot([x['each_window_model_time'] for x in rows],[x['probes'][-1]['total_gap'] for x in rows],'o-',label=label)
ax.axhline(0,color='#999999',lw=.7);ax.set_xlabel('Each A/B preparation window (model time)')
ax.set_ylabel('CDF(AB) - CDF(BA), common probe')
ax.set_title('Load order and boundary sensitivity');ax.legend(fontsize=8)
for ax in axes:ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2)
fig.suptitle('CRACK: original finite SG model comparisons | declared preparations | physical calibration open',fontsize=12)
fig.savefig(root/'CORRELATION_MEMORY_COMPARISON.png',dpi=170)
print(str(root/'CORRELATION_MEMORY_COMPARISON.png'))
