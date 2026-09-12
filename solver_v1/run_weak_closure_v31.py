"""Predict separate forced MD responses from unforced full C/K, without refit.

Finite-band K remains a proxy, not a zero-frequency calibration. The dynamic
records are not used to choose a band or fit any matrix in this comparison.
"""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .collective_generator_bridge import plane_zero_sum_basis, projected_susceptibility
from .run_low_frequency_forcing_v29 import sha
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def run(responses, bridges, out):
    responses, bridges, out = map(Path, (responses, bridges, out))
    if out.exists():
        raise FileExistsError('fresh closure comparison required')
    measured = json.loads((responses/'summary.json').read_bytes())
    if not measured['completed']:
        raise ValueError('completed weak-response campaign required')
    frequency = measured['protocol']['frequency_per_ps']
    omega = 2*np.pi*frequency*1e12
    rows, hashes = [], {}
    for seed in measured['protocol']['seeds']:
        report_path = bridges/f'seed{seed}_new'/'summary.json'
        report = json.loads(report_path.read_bytes())
        null = next(r for r in measured['checks'] if r['seed'] == seed and r['axis'] is None)
        if not report['completed'] or report['source_sha256'] != null['trajectory_sha256']:
            raise ValueError('unforced source binding mismatch')
        hashes[str(seed)] = sha(report_path)
        basis = plane_zero_sum_basis(report['plane_count'])
        thermal = 1.380649e-23*report['temperature_K']
        for item in report['matrices']:
            # Full 1900ps comparison; half-record variability remains in the
            # original campaign, not counted as new independent observations.
            if item['parts'] != 1 or item['inversion_failure'] is not None:
                continue
            H, M, C, K = (np.asarray(item[k]) for k in ('H_J_m2', 'M_m2_J_s', 'C_m2', 'K_m2_s'))
            for axis in (0, 1):
                v = basis[axis]
                full = projected_susceptibility(H, M, v, omega)*16.02176634
                c, k = float(v@C@v), float(v@K@v)
                static = c/thermal*16.02176634
                scalar_m = c*c/(thermal*k)
                scalar = 1/(thermal/c+1j*omega/scalar_m)*16.02176634
                obs = next(r for r in measured['pairs'] if r['seed'] == seed and r['axis'] == axis and r['parts'] == 1)
                observed = complex(obs['real_A2_eV'], obs['imag_A2_eV'])
                rows.append(dict(seed=seed, axis=axis, nw=item['nw'], lo_cycles_ps=item['band'][0],
                    hi_cycles_ps=item['band'][1], frequency_cycles_ps=frequency,
                    static_A2_eV=static, observed_real_A2_eV=observed.real, observed_imag_A2_eV=observed.imag,
                    full_real_A2_eV=full.real, full_imag_A2_eV=full.imag,
                    scalar_real_A2_eV=scalar.real, scalar_imag_A2_eV=scalar.imag,
                    full_phase_deg=float(np.degrees(np.angle(full))), scalar_phase_deg=float(np.degrees(np.angle(scalar))),
                    observed_phase_deg=obs['phase_deg'], observed_storage_over_static=observed.real/static,
                    full_relative_complex_error=abs(full-observed)/abs(observed),
                    scalar_relative_complex_error=abs(scalar-observed)/abs(observed),
                    full_loss_error_over_observed_null=abs(full.imag-observed.imag)/obs['imag_null'],
                    observed_loss_over_null=obs['loss_to_null'], no_response_refit=True,
                    zero_frequency_calibrated=False, production_clock_calibrated=False))
    if len(rows) != 24:
        raise ValueError('two seeds, two axes and six bands/windows required')
    out.mkdir(parents=True)
    write_csv(out/'response_prediction.csv', rows)
    # All bands are shown. Ranges are observed method sensitivity, NOT CIs.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with (responses/'fdt.csv').open(newline='') as stream:
        fdt=list(csv.DictReader(stream))
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    seeds=measured['protocol']['seeds']
    for axis in (0,1):
        for index,seed in enumerate(seeds):
            subset=[r for r in rows if r['seed']==seed and r['axis']==axis]
            obs=next(r for r in measured['pairs'] if r['seed']==seed and r['axis']==axis and r['parts']==1)
            for col,field in enumerate(('real_A2_eV','imag_A2_eV')):
                ax=axes[axis,col];sign=1 if col==0 else -1
                val=sign*obs[field]
                bound=obs['complex_null'] if col==0 else obs['imag_null']
                ax.errorbar(index-.15,val,yerr=bound,fmt='o',color='tab:blue',capsize=3,
                            label='Driven +/- pair and observed null' if index==0 else None)
                values=np.array([sign*r['full_'+field] for r in subset])
                mid=float(values.mean())
                ax.errorbar(index,mid,yerr=[[mid-values.min()],[values.max()-mid]],fmt='s',
                            color='tab:orange',capsize=3,label='Unforced full C/K: all bands' if index==0 else None)
                if col==0:
                    ax.plot(index+.15,subset[0]['static_A2_eV'],'x',color='black',
                            label='Unforced static covariance' if index==0 else None)
                else:
                    values=np.array([-float(r['imag_A2_eV']) for r in fdt if int(r['seed'])==seed
                        and int(r['axis'])==axis and int(r['parts'])==1])
                    mid=float(values.mean())
                    ax.errorbar(index+.15,mid,yerr=[[mid-values.min()],[values.max()-mid]],fmt='x',
                                color='black',capsize=3,label='Unforced FDT near drive frequency' if index==0 else None)
        for col in (0,1):
            ax=axes[axis,col]
            ax.set_xticks(range(len(seeds)),[f'seed {seed}' for seed in seeds])
            ax.set_ylabel(('Re susceptibility' if col==0 else '-Im susceptibility')+' [Angstrom²/eV]')
            ax.set_title(('Normal' if axis==0 else 'Registry')+(' storage' if col==0 else ' loss'))
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
            ax.grid(alpha=.2)
            ax.legend(fontsize=7,loc='best')
    fig.suptitle('Reference MD: 0.2 cycles/ps; production physical clock remains uncalibrated\n'
                 'Bars: observed null / band sensitivity, not confidence intervals',fontsize=11)
    fig.savefig(out/'unforced_vs_driven.png',dpi=160)
    plt.close(fig)
    save_json(out/'summary.json', dict(completed=True, response_sha256=sha(responses/'summary.json'),
        bridge_sha256=hashes, rows=rows, no_dynamic_refit=True, selected_best_band=False,
        condition='same reference potential/box and new initialization, NVE temperatures retained',
        finite_band_not_zero_frequency=True, material_accepted=False, production_clock_calibrated=False))
    print(json.dumps(dict(rows=len(rows), maximum_relative_complex_error=max(r['full_relative_complex_error'] for r in rows),
        storage_over_static_range=[min(r['observed_storage_over_static'] for r in rows),max(r['observed_storage_over_static'] for r in rows)])))


if __name__ == '__main__':
    p=argparse.ArgumentParser(__doc__)
    for field in ('responses','bridges','out'):
        p.add_argument('--'+field,type=Path,required=True)
    a=p.parse_args()
    run(a.responses,a.bridges,a.out)
