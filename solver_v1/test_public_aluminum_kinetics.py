"""Deterministic reader/correlation tests, not invented Al kinetic data."""
import io

import numpy as np
import pytest

from .public_aluminum_kinetics import (
    FCCPlaneProjection, HashedLimitedReader, correlation_audit,
    lag_covariances, lammps_frames,
    reversible_harmonic_covariance,
)


def crystal(n=3):
    ijk = np.stack(np.meshgrid(*[np.arange(2*n)]*3, indexing='ij'), axis=-1).reshape(-1, 3)
    ijk = ijk[ijk.sum(axis=1) % 2 == 0]
    return np.arange(1, len(ijk)+1), ijk/(2*n), np.array([[0., n*4.065]]*3)


def test_plane_geometry_translation_invariance_and_sign():
    ids, positions, bounds = crystal()
    projection = FCCPlaneProjection(ids, positions, bounds, repeats=3, lattice_angstrom=4.065)
    assert np.all(projection.counts == 36)
    np.testing.assert_allclose(projection.evaluate(ids, positions, bounds), 0., atol=1e-25)
    shifted = (positions+np.array([.001, -.002, .003])) % 1
    np.testing.assert_allclose(projection.evaluate(ids, shifted, bounds), 0., atol=1e-25)
    displacement = np.zeros_like(positions)
    displacement[projection.plane == 1] = .01*projection.basis[0]/np.diff(bounds).ravel()
    q = projection.evaluate(ids, positions+displacement, bounds)
    np.testing.assert_allclose(q[:, 0], [1e-12, -1e-12, 0.], atol=1e-25)
    np.testing.assert_allclose(q[:, 1:], 0., atol=1e-25)


def test_projection_rejects_npt_and_large_site_change():
    ids, p, bounds = crystal()
    obj = FCCPlaneProjection(ids, p, bounds, repeats=3, lattice_angstrom=4.065)
    with pytest.raises(ValueError, match='box'):
        obj.evaluate(ids, p, bounds*1.01)
    with pytest.raises(ValueError, match='ambiguous'):
        obj.evaluate(ids, p+np.array([.1, 0, 0]), bounds)
    with pytest.raises(ValueError, match='FCC sites'):
        FCCPlaneProjection(ids, p+np.array([1/6, 0, 0]), bounds,
                           repeats=3, lattice_angstrom=4.065)


def test_dump_ordering_and_truncation():
    text = ('ITEM: TIMESTEP\n5\nITEM: NUMBER OF ATOMS\n2\n'
            'ITEM: BOX BOUNDS pp pp pp\n0 4\n0 4\n0 4\n'
            'ITEM: ATOMS id type xs ys zs vx vy vz\n'
            '2 1 .5 0 .5 1 0 0\n1 1 0 0 0 0 0 0\n')
    step, ids, p, bounds = next(lammps_frames(io.StringIO(text), 1))
    assert step == 5
    assert ids.tolist() == [1, 2]
    assert p[0].tolist() == [0., 0., 0.]
    with pytest.raises(ValueError, match='truncated'):
        list(lammps_frames(io.StringIO(text.rsplit('\n', 2)[0]+'\n'), 1))
    with pytest.raises(ValueError, match='orthogonal'):
        list(lammps_frames(io.StringIO(text.replace('pp pp pp', 'xy xz yz pp pp pp')), 1))


def test_ordered_covariance_and_constant_offset():
    q = np.arange(60., dtype=float).reshape(10, 2, 3)**2
    q[:, :, 1] = q[::-1, :, 1]
    actual = lag_covariances(q, [0, 2])
    dq = q-q.mean(axis=0)
    expected = sum(np.outer(dq[t+2, p], dq[t, p]) for t in range(8) for p in range(2))/16
    np.testing.assert_allclose(actual[1], expected)
    assert not np.allclose(expected, expected.T)
    np.testing.assert_allclose(lag_covariances(q+17, [0, 2]), actual)
    with pytest.raises(ValueError):
        lag_covariances(q, [.5])


def test_oscillatory_measured_coordinate_does_not_get_clock():
    # Exact deterministic oscillator: no trajectory counting or random data.
    t = np.arange(2048)
    q = np.stack([np.cos(2*np.pi*t/32), np.sin(2*np.pi*t/32)], axis=-1)[:, None]*1e-12
    result = correlation_audit(q, frame_seconds=25e-15, max_lag=32, blocks=4)
    assert result['negative_beyond_empirical_floor'][16]
    assert not result['production_calibration_available']
    assert result['lags_seconds'][16] == 400e-15
    with pytest.raises(ValueError, match='four lag windows'):
        correlation_audit(q[:100], frame_seconds=25e-15, max_lag=32, blocks=4)


def test_prefix_reader_hard_budget_and_checksum():
    import hashlib
    reader = HashedLimitedReader(io.BytesIO(b'abcdef'), 4)
    assert reader.read(2) == b'ab'
    assert reader.read(5) == b'cd'
    assert reader.count == 4
    assert reader.digest.hexdigest() == hashlib.sha256(b'abcd').hexdigest()
    with pytest.raises(ValueError, match='budget'):
        reader.read(1)


def test_noncommuting_harmonic_covariance_uses_correct_whitening_basis():
    from scipy.linalg import expm
    H=np.array([[3.,.7],[.7,1.]])
    M=np.diag([.8,.06]);kT=.02;times=np.array([0.,.1,2.])
    C=reversible_harmonic_covariance(H,M,kT,times)
    for t,value in zip(times,C):
        np.testing.assert_allclose(value,expm(-M@H*t)@(kT*np.linalg.inv(H)),rtol=2e-14,atol=1e-17)
        np.testing.assert_allclose(value,value.T,atol=1e-17)
        assert np.linalg.eigvalsh(value)[0]>0
    h,V=np.linalg.eigh(H);white=(V*np.sqrt(h/kT))@V.T
    root=(V*np.sqrt(h))@V.T
    np.testing.assert_allclose(white@C[-1]@white,expm(-root@M@root*times[-1]),atol=1e-14)
    wrong=np.diag(np.sqrt(np.diag(M)))@H@np.diag(np.sqrt(np.diag(M)))
    assert not np.allclose(white@C[-1]@white,expm(-wrong*times[-1]))


def test_periodic_spatial_cancellation_does_not_hide_opposite_plane_drift():
    t=np.arange(1024.)
    values=np.stack([np.sin(t*.27)+t/150,np.cos(t*.37)-t/250],axis=-1)*1e-12
    q=np.stack([values,-values],axis=1)
    result=correlation_audit(q,frame_seconds=25e-15,max_lag=16,blocks=4)
    np.testing.assert_allclose(result['block_mean_shifts_in_rms_units'],0.,atol=1e-15)
    assert result['spatial_mean_shift_is_not_stationarity_evidence']
    assert result['maximum_plane_mean_shift_in_rms_units']>1
