"""Independently reconstruct nonlinear diagnostics from saved raw arrays.

No potential calls, path optimization, transition barriers or rates. The audit
checks the complete sampling schedule, source hashes, coordinate convention,
linear residual forces and even/odd energy decomposition. Partial runs are
explicitly partial and cannot satisfy the complete-schedule assertion.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.constants import Boltzmann, electron_volt


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(actual, expected, label):
    a=np.asarray(actual,dtype=float); b=np.asarray(expected,dtype=float)
    if a.shape!=b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError(label+': nonfinite or incompatible shape')
    if not np.allclose(a,b,rtol=1e-9,atol=1e-10):
        raise ValueError(label+': raw replay mismatch')
    return float(np.max(np.abs(a-b))) if a.size else 0.


def audit(results, probes):
    protocol=read(probes/'protocol.json'); summary=read(probes/'summary.json')
    points=read(probes/'points.json') if (probes/'points.json').exists() else []
    if protocol['modes']!=[0,1,2] or protocol['fixed_max_atom_displacements_A']!=[.01,.05,.1,.2]:
        raise ValueError('unknown probe schedule')
    if protocol['quadratic_coordinate_amplitudes_at_300K']!=[1.,2.] or protocol['maximum_atom_bound_A']!=.3:
        raise ValueError('unknown thermal amplitude schedule')
    if summary['first_crack_certified'] or summary['physical_clock'] is not None:
        raise ValueError('line probes cannot certify first passage or a physical clock')
    kb=Boltzmann/electron_volt; kt=300*kb
    replay=[]; manifest=[]; schedule=[]; skipped=[]; pairs={}; max_error=0.
    sources={s['state']:s for s in protocol['sources']}
    if set(sources)!={'loading8','return8'}:
        raise ValueError('both same-grip sources must be declared')
    replays=summary['completed_source_replays']
    for name in ('loading8','return8'):
        folder=results/('dense_'+name)
        hpath=folder/'raw_hessian.npz'; epath=folder/'spectrum.npz'
        for path,key in ((hpath,'raw_sha256'),(epath,'spectrum_sha256')):
            if sha(path)!=sources[name][key]:
                raise ValueError('source hash mismatch: '+name)
            manifest.append(dict(path=path.relative_to(results).as_posix(),sha256=sha(path)))
        with np.load(hpath) as d: source={k:d[k].copy() for k in d.files}
        with np.load(epath) as d: spectrum={k:d[k].copy() for k in d.files}
        h=(source['hessian']+source['hessian'].T)/2
        b=source['basis']; g=source['gradient']; r=source['positions']; energy=float(source['energy'])
        directions=[]; expected=[]
        for mode in range(3):
            v=spectrum['eigenvectors'][:,mode].copy()
            if v[np.argmax(abs(v))]<0:v=-v
            direction=(b@v).reshape(r.shape); directions.append(direction)
            lam=float(spectrum['eigenvalues'][mode])
            if lam<=0:raise ValueError('positive curvature required for thermal amplitudes')
            close(h@v,lam*v,'eigenvector')
            peak=float(np.linalg.norm(direction,axis=1).max())
            amplitudes=[('max_atom_'+str(a),a/peak) for a in (.01,.05,.1,.2)]
            amplitudes += [('quadratic_sigma_'+str(a),a*np.sqrt(kt/lam)) for a in (1.,2.)]
            for label,q in amplitudes:
                if q*peak>.3+1e-15:
                    skipped.append((name,mode,label))
                    continue
                for sign in (-1,1):expected.append((name,mode,label,sign,float(sign*q)))
        schedule.extend(expected)
        stored=[p for p in points if p['state']==name]
        rawfile=probes/(name+'_raw.npz')
        if not rawfile.exists():
            if stored:raise ValueError('reported points have no raw arrays')
            continue
        with np.load(rawfile) as d: raw={k:d[k].copy() for k in d.files}
        manifest.append(dict(path=rawfile.relative_to(results).as_posix(),sha256=sha(rawfile)))
        if len(stored)!=len(raw['energies']) or len(stored)!=len(raw['forces']) or len(stored)!=len(raw['evaluated']):
            raise ValueError('raw and reported point counts differ')
        for key in ('positions','numbers','cell','basis'):
            if not np.array_equal(raw[key],source[key]):raise ValueError('raw coordinates changed: '+key)
        close(raw['baseline_energy'],source['energy'],'baseline energy')
        close(raw['baseline_forces'],source['forces'],'baseline forces')
        close(raw['directions'],np.asarray(directions)[:len(raw['directions'])],'Cartesian direction')
        for index,p in enumerate(stored):
            mode=int(p['mode']); x=float(raw['evaluated'][index,1])
            close(raw['evaluated'][index,0],mode,'mode index')
            close(p['coordinate_A'],x,'coordinate')
            v=b.T@directions[mode].ravel(); lam=float(v@h@v); linear=float(g@v)
            f=raw['forces'][index]; delta=float(raw['energies'][index]-raw['baseline_energy'])
            actual_gradient=-b.T@f.ravel()
            residual_gradient=actual_gradient-g-x*(h@v)
            quadratic=linear*x+.5*lam*x*x
            values=dict(curvature_eV_A2=lam,baseline_directional_gradient_eV_A=linear,
                maximum_atom_displacement_A=float(abs(x)*np.linalg.norm(directions[mode],axis=1).max()),
                actual_delta_energy_eV=delta,linear_quadratic_delta_energy_eV=quadratic,
                nonlinear_energy_residual_eV=delta-quadratic,
                nonlinear_energy_residual_over_kBT300=(delta-quadratic)/kt,
                actual_directional_gradient_eV_A=float(actual_gradient@v),
                nonlinear_full_gradient_norm_eV_A=float(np.linalg.norm(residual_gradient)),
                nonlinear_directional_gradient_eV_A=float(residual_gradient@v))
            for key,value in values.items():max_error=max(max_error,close(p[key],value,key))
            key=(name,mode,p['label'])
            if p['sign'] in pairs.setdefault(key,{}):raise ValueError('duplicate signed point')
            pairs[key][p['sign']]=(delta,x,linear,lam,p)
            replay.append((name,mode,p['label'],p['sign'],x))
    if len(replay)!=len(points) or len(points)!=summary['actual_points']:
        raise ValueError('reported total or unknown state')
    if len(replay)>len(schedule):raise ValueError('more points than declared schedule')
    for observed,expected in zip(replay,schedule):
        if observed[:4]!=expected[:4]:raise ValueError('point schedule changed')
        close(observed[4],expected[4],'amplitude schedule')
    if summary['complete'] and len(replay)!=len(schedule):raise ValueError('incomplete schedule marked complete')
    expected_skipped=skipped if summary['complete'] else skipped[:len(summary['skipped'])]
    actual_skipped=[(p['state'],p['mode'],p['label']) for p in summary['skipped']]
    if actual_skipped!=expected_skipped:raise ValueError('unreported displacement-bound skip')
    for pair in pairs.values():
        if set(pair)!={-1,1}:continue
        negative,positive=pair[-1],pair[1]
        dm,xm,linear,lam,_=negative; dp,xp,_,_,_=positive
        close(xm,-xp,'symmetric amplitudes')
        even=(dp+dm)/2-.5*lam*xp*xp; odd=(dp-dm)/2-linear*xp
        for item in (negative,positive):
            max_error=max(max_error,close(item[4]['symmetric_nonlinear_energy_eV'],even,'even energy'),
                close(item[4]['antisymmetric_nonlinear_energy_eV'],odd,'odd energy'))
    if summary['new_model_calls']!=len(points)+len(replays):raise ValueError('new call accounting mismatch')
    return dict(complete=bool(summary['complete']),raw_replay_passed=True,
        verified_points=len(points),planned_allowed_points=len(schedule),skipped_pairs=len(actual_skipped),
        complete_signed_pairs=sum(set(p)=={-1,1} for p in pairs.values()),maximum_scalar_replay_error=max_error,
        numerical_comparison_rtol=1e-9,numerical_comparison_atol=1e-10,
        source_manifest=manifest,protocol_sha256=sha(probes/'protocol.json'),
        reported_summary_sha256=sha(probes/'summary.json'),new_model_calls=0,new_DFT=0,new_MD=0,
        first_crack_certified=False,physical_clock=None,
        scope='raw algebra and declared sampling schedule; no material, basin, path or free-energy certification')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for name in ('results','probes','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise ValueError('fresh audit output required')
    report=audit(args.results,args.probes)
    report['runner_sha256']=sha(Path(__file__))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(report,indent=2))
