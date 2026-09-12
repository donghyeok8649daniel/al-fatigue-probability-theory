"""Compare completed v29 signed/amplitude/timestep controls, without a clock gate."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .impedance_mobility import inverse_disk_bounds
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def report(results):
    results=Path(results);out=results/'comparison.json'
    if out.exists():raise FileExistsError('existing comparison preserved')
    data=json.loads((results/'summary.json').read_bytes())
    if not data['completed']:raise ValueError('completed records required')
    with (results/'fdt_spectral.csv').open() as stream:fdt=list(csv.DictReader(stream))
    protocol=data['protocol'];rows=[]
    for axis in (0,1):
        rr=[r for r in data['paired'] if r['axis']==axis]
        base=next(r for r in rr if r['fraction']==4 and r['dt_ps']==.0025 and r['parts']==1)
        fine=next(r for r in rr if r['fraction']==4 and r['dt_ps']==.00125 and r['parts']==1)
        half=next(r for r in rr if r['fraction']==2 and r['parts']==1)
        z=lambda r:complex(r['real_A2_eV'],r['imag_A2_eV'])
        ddt=abs(z(base)-z(fine));damp=abs(z(base)-z(half))
        block=max(abs(z(r)-z(base if r['dt_ps']==.0025 else fine)) for r in rr
            if r['fraction']==4 and r['parts']==2)
        radius=fine['complex_null_radius_A2_eV']+ddt+damp+block
        C0=protocol['rms_m'][axis]**2
        static=C0/(1.380649e-23*protocol['temperature_K'])*16.02176634
        expected=[float(r['imag_A2_eV']) for r in fdt if int(r['axis'])==axis]
        rows.append(dict(axis=axis,frequency_per_ps=.2,
            base_real_A2_eV=base['real_A2_eV'],base_imag_A2_eV=base['imag_A2_eV'],
            fine_real_A2_eV=fine['real_A2_eV'],fine_imag_A2_eV=fine['imag_A2_eV'],
            base_phase_deg=base['phase_deg'],fine_phase_deg=fine['phase_deg'],
            half_force_phase_deg=half['phase_deg'],
            base_loss_to_null=base['loss_to_null_ratio'],fine_loss_to_null=fine['loss_to_null_ratio'],
            half_force_loss_to_null=half['loss_to_null_ratio'],
            dt_complex_difference_A2_eV=ddt,amplitude_complex_difference_A2_eV=damp,
            full_force_block_difference_A2_eV=block,
            dt_difference_to_sum_null=ddt/(base['complex_null_radius_A2_eV']+fine['complex_null_radius_A2_eV']),
            amplitude_difference_to_sum_null=damp/(base['complex_null_radius_A2_eV']+half['complex_null_radius_A2_eV']),
            base_storage_to_static=base['real_A2_eV']/static,fine_storage_to_static=fine['real_A2_eV']/static,
            fdt_imag_min_A2_eV=min(expected),fdt_imag_max_A2_eV=max(expected),
            fine_loss_interval_overlaps_fdt=(fine['imag_A2_eV']+fine['observed_null_envelope_A2_eV']>=min(expected)
                and fine['imag_A2_eV']-fine['observed_null_envelope_A2_eV']<=max(expected)),
            base_point_mobility_m2_J_s=base['mobility_m2_per_J_s'],fine_point_mobility_m2_J_s=fine['mobility_m2_per_J_s'],
            control_radius_A2_eV=radius,**inverse_disk_bounds(z(fine),radius,.2),
            max_even_harmonic2_to_null=max(r['parity_harmonic2_to_null'] for r in rr if r['parts']==1),
            max_odd_harmonic3_to_null=max(r['parity_harmonic3_to_null'] for r in rr if r['parts']==1),
            zero_frequency_limit_certified=False,local_coordinate_mapping_validated=False))
    write_csv(results/'control_comparison.csv',rows)
    save_json(out,dict(completed=True,rows=rows,production_clock_calibrated=False,
        no_independence_or_confidence_claim=True))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for axis,ax in enumerate(axes):
        for r in data['paired']:
            if r['axis']!=axis or r['parts']!=1:continue
            ax.errorbar(r['fraction']+(.12 if r['dt_ps']==.00125 else 0),r['imag_A2_eV'],
                yerr=r['observed_null_envelope_A2_eV'],fmt='s' if r['dt_ps']==.00125 else 'o',capsize=4,
                color='C2' if r['dt_ps']==.00125 else 'C0')
        row=rows[axis]
        ax.axhspan(row['fdt_imag_min_A2_eV'],row['fdt_imag_max_A2_eV'],color='C1',alpha=.3)
        ax.axhline(0,color='k',lw=.8)
        ax.set(xlabel='Force / thermal-RMS scale (fine dt offset +0.12)',ylabel='Im susceptibility [Å²/eV]',
            title=('Normal','Direct110 slip')[axis]+' — reference 0.2 cycles/ps')
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.grid(alpha=.2)
    fig.suptitle('Actual new MD; bars: observed null, shade: FDT sensitivity — NOT confidence intervals')
    fig.savefig(results/'loss_comparison.png',dpi=160);plt.close(fig)
    print(json.dumps(rows,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('--results',type=Path,required=True)
    report(p.parse_args().results)
