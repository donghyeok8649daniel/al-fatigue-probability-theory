"""Surface traction bookkeeping in MPa, mm, N, N mm; not a mechanics solver."""
import numpy as np


class LoadBalanceError(ValueError):
    def __init__(self, net, time=None):
        super().__init__('solid.unbalanced')
        values = {name: ', '.join(f'{v:.6g}' for v in vector) for name, vector in
                  (('force', net[:3]), ('torque', net[3:]),
                   ('opposite_force', -net[:3]), ('opposite_torque', -net[3:]))}
        values['time'] = '—' if time is None else f'{time:.8g}'
        self.ui_error_data = dict(key='solid.unbalanced_detail', values=values)


def balance_report(mesh, traction):
    traction = np.asarray(traction, dtype=float)
    if traction.shape != (len(mesh.faces), 3) or not np.isfinite(traction).all():
        raise ValueError('solid.traction')
    net = resultant(mesh, traction)
    length = np.linalg.norm(np.ptp(mesh.vertices, axis=0))
    scale = float(np.linalg.norm(traction*mesh.areas[:, None], axis=1).sum())
    relative = float(np.max(np.abs(np.r_[net[:3], net[3:]/length]))/scale) if scale else 0.
    return dict(net=net, relative=relative, balanced=relative <= 1e-9,
                unloaded=scale == 0., origin=mesh.vertices.mean(axis=0))


def applied_traction(mesh, loads, time, correction=None, *, compiled=None):
    """Only stored assignments apply forces. Draft stress fields apply nothing."""
    from .tensor_load import compile_tensor_matrix, evaluate_tensor
    if not np.isfinite(time): raise ValueError('load.invalid_time')
    result = np.zeros((len(mesh.faces), 3))
    if compiled is None:
        compiled = [compile_tensor_matrix(load.tensor_expression) for load in loads]
    _, normals = face_geometry(mesh)
    for load, expression in zip(loads, compiled):
        ids = np.asarray(load.face_indices, dtype=int)
        if not len(ids) or ids.min() < 0 or ids.max() >= len(mesh.faces):
            raise ValueError('solid.loads')
        stress = evaluate_tensor(expression, t=time, frequency=load.frequency,
            normal_mean=load.normal_mean_mpa, normal_amplitude=load.normal_amplitude_mpa,
            shear_mean=load.shear_mean_mpa, shear_amplitude=load.shear_amplitude_mpa)
        if not np.allclose(stress, stress.T, rtol=1e-12, atol=1e-12):
            raise ValueError('solid.traction')
        result[ids] += normals[ids] @ stress.T
    return balanced_traction(mesh, result, correction) if correction is not None else result


def preflight_balance(mesh, loads, correction, end):
    """Reject common imbalance before meshing; evaluated-time checks still apply."""
    from .tensor_load import compile_tensor_matrix
    if not loads: return
    compiled = [compile_tensor_matrix(load.tensor_expression) for load in loads]
    times = set(np.linspace(0., end, 17))
    for load in loads:
        times.update(t for t in np.arange(5)/(4*load.frequency) if t <= end)
    for time in sorted(times):
        report = balance_report(mesh, applied_traction(mesh, loads, time, correction, compiled=compiled))
        if not report['balanced']: raise LoadBalanceError(report['net'], time)


def face_geometry(mesh):
    if hasattr(mesh, '_face_geometry'):
        centers, normals = mesh._face_geometry
        return centers.copy(), normals.copy()
    # Reject open/inconsistently oriented surfaces; normals must point outward.
    mesh.enclosed_volume_mm3
    tri = mesh.vertices[mesh.faces]
    cross = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    center = mesh.vertices.mean(axis=0)
    sign = np.sign(np.einsum('ij,ij->', tri[:, 0] - center, cross))
    return tri.mean(axis=1), sign * cross / np.linalg.norm(cross, axis=1)[:, None]


def stress_traction(mesh, ids, stress):
    stress = np.asarray(stress, dtype=float)
    if stress.shape != (3, 3) or not np.isfinite(stress).all():
        raise ValueError('finite 3x3 stress required')
    if not np.allclose(stress, stress.T, rtol=1e-12, atol=1e-12):
        raise ValueError('Cauchy stress must be symmetric')
    _, normals = face_geometry(mesh)
    out = np.zeros((len(mesh.faces), 3))
    out[ids] = normals[ids] @ stress.T
    return out


def resultant(mesh, traction):
    centers, _ = face_geometry(mesh)
    forces = np.asarray(traction) * mesh.areas[:, None]  # MPa mm² = N
    origin = mesh.vertices.mean(axis=0)
    return np.r_[forces.sum(axis=0), np.cross(centers-origin, forces).sum(axis=0)]


def correction_operator(mesh, ids):
    """Minimum integral |traction|² correction, six global balance constraints.

    The returned fixed linear operator acts on the instantaneous resultant;
    therefore it balances arbitrary histories, not just a sampled phase.
    Supports/reactions are NOT inferred. User must choose correction faces.
    """
    ids = np.unique(np.asarray(ids, dtype=int))
    if not len(ids) or ids.min() < 0 or ids.max() >= len(mesh.faces):
        raise ValueError('load.correction_faces')
    centers, _ = face_geometry(mesh)
    origin = mesh.vertices.mean(axis=0)
    length = np.linalg.norm(np.ptp(mesh.vertices, axis=0))
    blocks = []
    for r, area in zip((centers[ids]-origin)/length, mesh.areas[ids]):
        x, y, z = r
        cross = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
        blocks.append(np.vstack((np.eye(3), cross))*np.sqrt(area))
    matrix = np.concatenate(blocks, axis=1)
    singular = np.linalg.svd(matrix, compute_uv=False)
    if singular.size < 6 or singular[-1] <= 1e-12*singular[0]:
        raise ValueError('load.correction_rank')
    inverse = np.linalg.pinv(matrix, rcond=1e-12)
    inverse[:, 3:] /= length
    inverse /= np.repeat(np.sqrt(mesh.areas[ids]), 3)[:, None]
    return ids, -inverse


def balanced_traction(mesh, traction, correction):
    ids, operator = correction
    result = np.array(traction, copy=True)
    result[ids] += (operator @ resultant(mesh, traction)).reshape(-1, 3)
    return result
