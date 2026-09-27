"""Validate a saved direct-Hessian prefix before reusing it; no physics model change."""
from pathlib import Path
import hashlib,json
import numpy as np

IDENTITY=('state_sha256','geometry_sha256','model_sha256','ensemble','dimension',
          'grip_extension_basis_norm','external_nominal_stress_GPa','basis','units')

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def load_prefix(folder, current_protocol, basis, gradient, expected_sha256):
    folder=Path(folder); payload=folder/'partial.npz'
    if not expected_sha256 or sha(payload)!=expected_sha256:
        raise ValueError('partial data hash does not match pinned source')
    previous=json.loads((folder/'protocol.json').read_text(encoding='utf-8'))
    summary=json.loads((folder/'summary.json').read_text(encoding='utf-8'))
    for name in IDENTITY:
        if previous.get(name)!=current_protocol.get(name):
            raise ValueError('resume context changed: '+name)
    with np.load(payload,allow_pickle=False) as data:
        rows=data['hessian_rows'].copy();old_basis=data['basis'];old_gradient=data['gradient']
    dimension=current_protocol['dimension']
    if (rows.ndim!=2 or rows.shape[1]!=dimension or not 0<len(rows)<dimension
            or summary.get('complete') is not False
            or summary.get('completed_rows')!=len(rows) or summary.get('dimension')!=dimension):
        raise ValueError('partial row count or completion metadata invalid')
    if any(not np.isfinite(x).all() for x in (rows,old_basis,old_gradient,basis,gradient)):
        raise ValueError('nonfinite resume data')
    if old_basis.shape!=basis.shape or not np.array_equal(old_basis,basis):
        raise ValueError('Cartesian basis changed')
    if old_gradient.shape!=gradient.shape or not np.allclose(old_gradient,gradient,atol=1e-10,rtol=1e-10):
        raise ValueError('current gradient does not replay saved gradient')
    square=rows[:,:len(rows)]
    asym=float(np.max(abs(square-square.T)))
    if asym>1e-9*max(float(np.linalg.norm(square)),np.finfo(float).tiny):
        raise ValueError('saved principal block is asymmetric')
    return rows,dict(reused_rows=len(rows),partial_sha256=sha(payload),
        protocol_sha256=sha(folder/'protocol.json'),summary_sha256=sha(folder/'summary.json'),
        gradient_replay_max_error_eV_A=float(np.max(abs(old_gradient-gradient))),
        principal_raw_asymmetry_eV_A2=asym,source_elapsed_is_not_new_cpu_time=True)

def check_replayed_row(saved, recomputed):
    saved=np.asarray(saved);recomputed=np.asarray(recomputed)
    if saved.shape!=recomputed.shape or not np.isfinite(recomputed).all():
        raise ValueError('invalid independently recomputed row')
    error=float(np.linalg.norm(saved-recomputed))
    tolerance=1e-9*max(1.,float(np.linalg.norm(saved)))
    if error>tolerance:
        raise ValueError('saved Hessian row does not reproduce')
    return dict(error_norm_eV_A2=error,tolerance_eV_A2=tolerance)
