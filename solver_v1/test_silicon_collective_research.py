import numpy as np

from solver_v1.silicon_crack_research import RelaxedCoordinates
from solver_v1.silicon_collective_research import RetainedBathDirection


def test_collective_coordinate_is_orthonormal_keeps_gap_and_has_exact_dual_forces():
    reference=np.array([[0.,-.4,0.],[0.,.4,0.],[1.,0.,0.]])
    base=RelaxedCoordinates(reference,[False,False,True],bond=[0,1])
    retained=RetainedBathDirection(base,[.2,.7,.3,-.6,.1])
    matrix=np.column_stack([retained.lift(v,0.) for v in np.eye(retained.dimension)])
    np.testing.assert_allclose(matrix.T@matrix,np.eye(retained.dimension),atol=7e-16)
    np.testing.assert_allclose(matrix.T@retained.direction,0.,atol=4e-16)
    y=np.array([.2,-.3,.1,.4]);q,u=1.1,.17
    positions=retained.positions(y,[q,u])
    np.testing.assert_allclose(positions[1,1]-positions[0,1],q,atol=2e-16)
    np.testing.assert_allclose(retained.direction@base.encode(positions),u,atol=2e-16)
    np.testing.assert_allclose(retained.encode(positions),y,atol=2e-16)
    gradient=2*positions
    g,reaction=retained.pullback(gradient)
    h=1e-6
    energy=lambda variables,values:np.sum(retained.positions(variables,values)**2)
    for index in range(len(y)):
        delta=np.eye(len(y))[index]*h
        numeric=(energy(y+delta,[q,u])-energy(y-delta,[q,u]))/(2*h)
        np.testing.assert_allclose(numeric,g[index],atol=3e-10)
    for index in range(2):
        delta=np.eye(2)[index]*h
        numeric=(energy(y,np.array([q,u])+delta)-energy(y,np.array([q,u])-delta))/(2*h)
        np.testing.assert_allclose(numeric,reaction[index],atol=3e-10)
    a=np.random.default_rng(41).standard_normal((5,5));a=a.T@a
    np.testing.assert_allclose(retained.restricted_hessian(a),matrix.T@a@matrix,atol=4e-15)
