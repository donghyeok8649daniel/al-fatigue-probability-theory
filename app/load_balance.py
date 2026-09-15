"""Surface traction bookkeeping in MPa, mm, N, N mm; not a mechanics solver."""
import numpy as np


def face_geometry(mesh):
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
        raise ValueError('select correction faces')
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
        raise ValueError('selected faces cannot independently balance force and torque')
    inverse = np.linalg.pinv(matrix, rcond=1e-12)
    inverse[:, 3:] /= length
    inverse /= np.repeat(np.sqrt(mesh.areas[ids]), 3)[:, None]
    return ids, -inverse


def balanced_traction(mesh, traction, correction):
    ids, operator = correction
    result = np.array(traction, copy=True)
    result[ids] += (operator @ resultant(mesh, traction)).reshape(-1, 3)
    return result
