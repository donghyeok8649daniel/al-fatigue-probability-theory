"""Compute conditional harmonic quantum sensitivity, with no new material fit."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import eigvalsh

from solver_v1.silicon_conditional_research import KB_EV_K
from solver_v1.silicon_quantum_diagnostic import quantum_harmonic_difference
from .run_local_crack_audit import save_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--references',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    if any(args.output.iterdir()): parser.error('use an empty output')
    start=time.perf_counter()
    meta=json.loads((args.references/'summary.json').read_text(encoding='utf-8'))
    initial=np.load(args.references/'initial.npz',allow_pickle=False)
    spectrum0=eigvalsh(initial['bath_hessian']);logdet0=np.log(spectrum0).sum()
    rows=[];spectra={}
    for source in meta['rows']:
        data=np.load(args.references/(source['name']+'.npz'),allow_pickle=False)
        values=eigvalsh(data['bath_hessian']);spectra[source['name']]=values
        for temperature in [0.,100.,300.,600.]:
            quantum=quantum_harmonic_difference(values,spectrum0,temperature)
            classical=.5*KB_EV_K*temperature*(np.log(values).sum()-logdet0)
            rows.append(dict(reference=source['name'],gap_A=source['gap_A'],temperature_K=temperature,
                classical_harmonic_correction_eV=float(classical),quantum_harmonic_correction_eV=quantum,
                quantum_minus_classical_eV=float(quantum-classical),
                candidate_energy_difference_eV=float(source['energy_difference_eV']+quantum)))
        print(source['name'],rows[-2],flush=True)
    np.savez_compressed(args.output/'spectra.npz',**spectra)
    save_json(args.output/'summary.json',dict(rows=rows,elapsed_seconds=time.perf_counter()-start,
        interpretation='fixed classical q, local harmonic oscillator sensitivity on original empirical SW; not quantum Si validation',
        anharmonic_quantum_effects_included=False,tunneling_included=False,material_fit=False,
        code_sha256={name:hashlib.sha256(Path(name).read_bytes()).hexdigest() for name in (
            'solver_v1/silicon_quantum_diagnostic.py','results/silicon_wafer_feasibility/run_quantum_sensitivity.py')}))


if __name__=='__main__': main()
