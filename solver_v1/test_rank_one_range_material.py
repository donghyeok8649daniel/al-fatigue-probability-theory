import json

import numpy as np
import pytest

from .quadrupole_saturation import SaturatedQuadrupoleInterface, SaturatedQuadrupoleCache
from .rank_one_range_material import RankOneRangeInterface, RankOneRangeCache, RankOneRangeBulk, shape_for_bulk
from .isotropic_bulk_validation import IsotropicBulkBasis


C = np.array([.3, .5, 2., 1., .2, 3., .7, 2., .03, 1200.])
D = (3., 6., 3.5, 12.)


def test_tied_range_recovers_full_old_jet_and_only_rank1_changes():
    old = SaturatedQuadrupoleInterface(D, C)
    tied = RankOneRangeInterface((*D, D[1]), C)
    changed = RankOneRangeInterface((*D, 3.1), C)
    for q in ((old.h, 0., 0.), (1.08*old.h, .21, -.06)):
        a, b, z = [obj.evaluate(q) for obj in (old, tied, changed)]
        assert a.energy == b.energy
        np.testing.assert_array_equal(a.hessian, b.hessian)
        for name in a.components:
            if name != 'angular_1':
                np.testing.assert_array_equal(a.components[name], z.components[name])
        assert np.linalg.norm(a.components['angular_1']-z.components['angular_1']) > 1e-4
    np.testing.assert_allclose(changed.evaluate((old.h, 0., 0.)).components['angular_1'][1:4],
                               0., atol=1e-12, rtol=0.)
    # This synthetic C is not equilibrium calibrated; compare reference energy,
    # not zero total normal force for arbitrary positive test coefficients.
    assert abs(changed.evaluate((old.h, 0., 0.)).energy) < 1e-12


@pytest.mark.parametrize('q', [(1.03*np.sqrt(2/3), .23, .07),
                               (1.25*np.sqrt(2/3), .38, -.09)])
def test_full_analytic_derivative_refinement_and_registry_period(q):
    model = RankOneRangeInterface((*D, 3.1), C)
    q = np.asarray(q); v = model.evaluate(q); errors = []
    for step in (4e-5, 2e-5):
        grad, H = [], []
        for axis in range(3):
            d = np.eye(3)[axis]*step
            p, m = model.evaluate(q+d), model.evaluate(q-d)
            grad.append((p.energy-m.energy)/(2*step))
            H.append((p.gradient-m.gradient)/(2*step))
        errors.append(max(np.max(abs(np.array(grad)-v.gradient)), np.max(abs(np.array(H).T-v.hessian))))
    assert errors[1] < .35*errors[0] and errors[1] < 1e-5
    np.testing.assert_array_equal(v.hessian, v.hessian.T)
    repeat = model.evaluate(q+[0., 1., 0.])
    np.testing.assert_allclose(repeat.hessian, v.hessian, atol=3e-11, rtol=3e-11)


def test_per_site_direct_rank1_sum_converges_to_reciprocal():
    model = RankOneRangeInterface((*D, 3.1), C)
    q = (1.08*model.h, .21, -.06)
    moment = next(m for _, m in model.base.moments if m.invariant.rank == 1)
    exact = model.base._angular(q, moment)[0][0]
    coarse = model.direct_rank1_energy(q, radius=6, layers=6)
    fine = model.direct_rank1_energy(q, radius=16, layers=16)
    assert abs(fine-exact) < abs(coarse-exact)
    assert fine == pytest.approx(exact, rel=2e-10, abs=2e-12)


def test_new_dataset_and_exact_coefficient_column_mapping():
    from .run_vector_material_calibration import source_and_targets
    from .interface_normal_development_targets import normal_development_observations
    source, previous, states = source_and_targets()
    obs, provenance = normal_development_observations(source, previous, states)
    assert len(obs) == 157  # 5 bulk + 104 interface fit + 48 heldout
    assert len([o for o in obs if o.role == 'heldout']) == 48
    assert all(o.name.startswith('v17_') for o in obs if o.role == 'heldout')
    assert all(o.role == 'fit' for o in obs if o.name.startswith('v16_'))
    assert set(provenance) == {o.name for o in obs}
    shape = (*D, 3.1)
    old = SaturatedQuadrupoleCache(obs).matrix(D)
    new = RankOneRangeCache(obs).matrix(shape)
    np.testing.assert_array_equal(new[:5], old[:5])
    np.testing.assert_array_equal(new[:, [0, 1, 2, 3, 4, 5, 7, 8, 9]], old[:, [0, 1, 2, 3, 4, 5, 7, 8, 9]])
    model = RankOneRangeInterface(shape, C)
    for o, row in zip(obs[::13], new[::13]):
        if o.bulk_index is not None:
            continue
        v = model.evaluate(o.state)
        jet = np.r_[v.energy, v.gradient, v.hessian[0], v.hessian[1, 1:], v.hessian[2, 2]]
        assert row@C == pytest.approx(np.asarray(o.jet_weights)@jet, rel=2e-10, abs=2e-10)


def test_nonuniform_operator_uses_independent_range_and_direct_energy():
    shape = (*D, 3.1)
    new = RankOneRangeBulk(shape, stretch=1.007, radius=12.)
    old = IsotropicBulkBasis(D[:3], stretch=1.007, radius=12.)
    q = (1/3, 1/3, 1/3)
    columns, tails = new.evaluate(q); tied, _ = old.evaluate(q)
    assert np.linalg.norm(columns[6]-tied[6]) > 1e-4
    np.testing.assert_array_equal(columns[[0,1,2,3,4,5,7,8,9]], tied[[0,1,2,3,4,5,7,8,9]])
    H = np.einsum('c,cij->ij', C, columns)
    for v in np.eye(3):
        options = dict(planes=6, mode=2, polarization_plane=v, validation_radius=12., quadrupole_saturation=D[3])
        zero = new.direct_sinusoidal_energy(C, amplitude=0., **options)
        second = []
        for h in (1e-4, 5e-5, 2.5e-5):
            p = new.direct_sinusoidal_energy(C, amplitude=h, **options)
            m = new.direct_sinusoidal_energy(C, amplitude=-h, **options)
            second.append(2*(p+m-2*zero)/h**2)
        fourth = [(4*b-a)/3 for a, b in zip(second, second[1:])]
        sixth = (16*fourth[1]-fourth[0])/15
        assert abs(sixth-v@H@v) < 1e-5+tails@abs(C)
    coarse = RankOneRangeBulk(shape, stretch=1.007, radius=8.)
    hc, tc = coarse.evaluate(q)
    assert np.linalg.norm(hc[6]-columns[6], 2) < tc[6]+tails[6]


def test_rank1_is_higher_gradient_not_an_added_elastic_modulus():
    model = RankOneRangeBulk((*D, 3.1), radius=12.)
    values = [model.evaluate((q, 0., 0.))[0][6] for q in (.002, .001)]
    ratio = np.linalg.norm(values[0])/np.linalg.norm(values[1])
    assert ratio == pytest.approx(16., rel=1e-3)  # H1=O(q^4), not q^2


def test_loader_and_bulk_shape_refuse_dropping_rank1(tmp_path):
    from .validate_tail_calibration import load_material
    definition = dict(rank_one_range_extension=True)
    shape = [*D, 3.1]
    (tmp_path/'definition.json').write_text(json.dumps(definition))
    (tmp_path/'calibration.json').write_text(json.dumps(dict(completed=True,
        best=dict(decays=shape, coefficients=C.tolist(), strictly_positive_LJ=True))))
    _, de, model, cache_type, bulk_type = load_material(tmp_path)
    assert isinstance(model, RankOneRangeInterface) and cache_type is RankOneRangeCache
    assert shape_for_bulk(shape, de) == shape
    with pytest.raises(ValueError, match='five'):
        bulk_type(shape[:3])
    for invalid in (float('nan'), -1.):
        with pytest.raises(ValueError, match='five'):
            bulk_type([*D[:3], invalid, 3.1])
    assert bulk_type([*D[:3], 0., 3.1]).rank1_decay == 3.1
    assert bulk_type(shape).rank1_decay == shape[4]


def test_normal_report_keeps_all_independent_validation_quantities():
    from .report_normal_response import grouped_validation
    rows = [dict(observable=f'v17_state{i}_{field}', role='heldout',
                 normalized_residual=2. if i == 11 else 0.)
            for field in ('energy', 'normal_force', 'Haa', 'Hxx') for i in range(12)]
    groups = grouped_validation(rows)
    assert len(groups) == 4
    assert all(r['count_within_declared_scale'] == 11 for r in groups)
    assert all(r['normalized_rms'] == pytest.approx(2/np.sqrt(12)) for r in groups)
    with pytest.raises(ValueError, match='twelve'):
        grouped_validation(rows[:-1])


def test_source_cutoff_geometry_and_disputed_normal_curvature_are_independent():
    from .run_source_normal_audit import crossing_gaps
    from .run_vector_material_calibration import source_and_targets
    source, _, _ = source_and_targets()
    ref = source.reference
    for row in crossing_gaps(ref):
        counts = []
        for sign in (-1., 1.):
            a = (row['a_over_h']+sign*1e-8)*ref.h
            d = a+(row['layer']-1)*ref.h
            xy = ref.R+ref.geometry.abc_shift(row['layer'])
            counts.append(np.sum(d*d+np.sum(xy*xy, axis=1) < ref.r[-1]**2))
        assert counts[0]-counts[1] == row['inplane_degeneracy']
    q = np.array([1.62*source.h, 0., 0.]); step = 1e-5
    H = source.evaluate(q).hessian[0, 0]
    p, m = source.evaluate(q+[step, 0., 0.]), source.evaluate(q-[step, 0., 0.])
    assert (p.gradient[0]-m.gradient[0])/(2*step) == pytest.approx(H, abs=1e-6, rel=0.)
