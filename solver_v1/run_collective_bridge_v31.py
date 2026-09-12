"""Full zero-sum plane-generator audit from actual equilibrium trajectories."""
import argparse
from pathlib import Path
import json
import numpy as np
from scipy.signal.windows import dpss
from .collective_generator_bridge import plane_zero_sum_basis,covariance_mobility,energy_rescale,harmonic_generator
from .run_low_frequency_forcing_v29 import sha
from .run_vector_registry_audit import save_json
from .run_low_stress_cyclic_diagnostic import write_csv


def run(source,out):
    source,out=Path(source),Path(out)
    if out.exists():raise FileExistsError('existing bridge study preserved')
    m=json.loads((source/'summary.json').read_bytes())
    if not m['completed'] or m['ensemble']!='nve':raise ValueError('completed equilibrium NVE required')
    if m.get('conjugate_drive') is not None:raise ValueError('unforced covariance required')
    with np.load(source/'plane_coordinates.npz') as z:
        t=z['time_seconds'];mask=t>=100e-12;q=z['coordinates_m'][mask];T=float(z['thermo'][mask,0].mean())
    dt=m['frame_ps']*1e-12;basis=plane_zero_sum_basis(q.shape[1]);raw=q.reshape(len(q),-1)
    coords=raw@basis;closure=float(np.max(abs(raw-coords@basis.T)))
    kBT=1.380649e-23*T;rows=[];matrices=[]
    for parts in (1,2):
        width=len(coords)//parts
        for block in range(parts):
            X=coords[block*width:(block+1)*width];X=X-X.mean(axis=0);C=X.T@X/len(X)
            f=np.fft.rfftfreq(len(X),dt)
            for nw in (2.,3.):
                windows=dpss(len(X),nw,Kmax=2,sym=False)
                Fourier=[np.fft.rfft(X*w[:,None],axis=0) for w in windows]
                for lo,hi in ((.02,.04),(.04,.08),(.08,.16)):
                    keep=(f>=lo*1e12)&(f<hi*1e12)
                    if keep.sum()<3 or lo*1e12<=nw/(len(X)*dt):raise ValueError('unresolved band')
                    K=sum(dt*(Z[keep].T@Z[keep].conj()).real/(np.dot(w,w)*keep.sum())/2
                          for w,Z in zip(windows,Fourier))/len(windows)
                    failure=None;H=M=None;drift_error=diff_error=None
                    try:
                        H,M=covariance_mobility(C,K,kBT)
                        hbar,mbar,tbar=energy_rescale(H,M,kBT,m['atoms_per_plane'])
                        drift,diff=harmonic_generator(H,M,kBT)
                        drift2,diff2=harmonic_generator(hbar,mbar,tbar)
                        drift_error=float(np.linalg.norm(drift-drift2)/np.linalg.norm(drift))
                        diff_error=float(np.linalg.norm(diff-diff2)/np.linalg.norm(diff))
                    except ValueError as error:
                        failure=str(error)
                    matrices.append(dict(parts=parts,block=block,nw=nw,band=[lo,hi],C_m2=C,K_m2_s=K,
                        H_J_m2=H,M_m2_J_s=M,inversion_failure=failure))
                    for axis in (0,1):
                        v=basis[axis];Cq=float(v@C@v);Kq=float(v@K@v)
                        scalar=Cq*Cq/(kBT*Kq);full=None if M is None else float(v@M@v)
                        rows.append(dict(parts=parts,block=block,nw=nw,lo_cycles_ps=lo,hi_cycles_ps=hi,axis=axis,
                            scalar_eliminated_mobility_m2_J_s=scalar,full_projected_mobility_m2_J_s=full,
                            full_to_scalar=None if full is None else full/scalar,covariance_condition=float(np.linalg.cond(C)),
                            spectral_condition=float(np.linalg.cond(K)),
                            drift_rescale_relative=drift_error,diffusion_rescale_relative=diff_error,
                            spectral_rank=int(np.linalg.matrix_rank(K)),dimension=len(K),inversion_failure=failure,
                            wrong_unscaled_thermal_diffusion_factor=m['atoms_per_plane'],
                            production_clock_calibrated=False))
    out.mkdir(parents=True);write_csv(out/'bridge_metrics.csv',rows)
    save_json(out/'summary.json',dict(completed=True,source_sha256=sha(source/'plane_coordinates.npz'),
        plane_count=q.shape[1],atoms_per_plane=m['atoms_per_plane'],temperature_K=T,
        zero_sum_reconstruction_max_m=closure,rows=rows,matrices=matrices,
        local_coordinate_equivalence=False,zero_frequency_limit_certified=False,production_clock_calibrated=False))
    ratios=[r['full_to_scalar'] for r in rows if r['full_to_scalar'] is not None]
    print(json.dumps(dict(closure_m=closure,rows=len(rows),unresolved=sum(r['inversion_failure'] is not None for r in rows),
        ratio_range=[min(ratios),max(ratios)] if ratios else None)))


if __name__=='__main__':
    p=argparse.ArgumentParser(__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();run(a.source,a.out)
