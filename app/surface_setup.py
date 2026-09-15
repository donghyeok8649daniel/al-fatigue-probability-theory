"""Versioned declarative setup exchange; never executes scripts or a solver."""
from dataclasses import asdict
import hashlib
import json
import numpy as np
from .tensor_load import compile_tensor_matrix


def mesh_id(mesh):
    return hashlib.sha256(mesh.vertices.astype('<f8').tobytes()+mesh.faces.astype('<i8').tobytes()).hexdigest()


def encode(mesh, loads, correction):
    return json.dumps(dict(schema='aft.surface-loads/1', mesh_sha256=mesh_id(mesh),
        units=dict(coordinates='mm', stress='MPa', frequency='cycles/model_time'),
        loads=[asdict(load) for load in loads],
        correction_faces=None if correction is None else correction[0].tolist(),
        scope='surface_setup_only_not_spatial_solution'), indent=2, allow_nan=False)


def decode(text, mesh):
    from .load_workflow import FaceLoad
    from .load_balance import correction_operator
    import ast
    if len(text) > 4_000_000: raise ValueError('setup too large')
    value = json.loads(text)
    if not isinstance(value, dict): raise ValueError('setup must be an object')
    if value.get('schema') != 'aft.surface-loads/1' or value.get('mesh_sha256') != mesh_id(mesh):
        raise ValueError('setup schema/mesh mismatch')
    if value.get('units') != dict(coordinates='mm', stress='MPa', frequency='cycles/model_time'):
        raise ValueError('setup units mismatch')
    loads = []
    for row in value['loads']:
        row = dict(row)
        ids = row['face_indices']
        if not ids or any(type(i) is not int or i < 0 or i >= len(mesh.faces) for i in ids) or len(set(ids)) != len(ids):
            raise ValueError('invalid faces')
        row['face_indices'] = tuple(ids)
        row['area_mm2'] = float(mesh.areas[ids].sum())
        load = FaceLoad(**row)
        if not np.isfinite([load.frequency, load.normal_mean_mpa, load.normal_amplitude_mpa,
                            load.shear_mean_mpa, load.shear_amplitude_mpa]).all() or load.frequency <= 0:
            raise ValueError('invalid load values')
        compile_tensor_matrix(load.tensor_expression)
        expressions = [r.split(',') for r in load.tensor_expression.split(';')]
        for i in range(3):
            for j in range(i):
                if ast.dump(ast.parse(expressions[i][j].strip(), mode='eval')) != ast.dump(ast.parse(expressions[j][i].strip(), mode='eval')):
                    raise ValueError('stress must be symmetric')
        loads.append(load)
    ids = value.get('correction_faces')
    # Return a proposal only. Import never approves a correction itself.
    if ids is not None:
        if any(type(i) is not int for i in ids): raise ValueError('invalid correction faces')
        correction_operator(mesh, ids)
    return loads, ids
