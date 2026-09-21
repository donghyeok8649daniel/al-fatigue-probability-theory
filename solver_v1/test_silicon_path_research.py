import numpy as np
from solver_v1.silicon_path_research import relax_neb,reparameterize_path


def test_climbing_band_finds_known_curved_saddle_and_keeps_endpoints():
    # V=(x²-1)²+5[y-.3(1-x²)]² has minima (+/-1,0), saddle (0,.3), barrier1.
    def evaluate(index,point):
        x,y=point;offset=y-.3*(1-x*x)
        return (x*x-1)**2+5*offset*offset,np.array([4*x*(x*x-1)+6*x*offset,10*offset])
    initial=reparameterize_path([[-1.,0.],[0.,.05],[1.,0.]],13)
    result=relax_neb(evaluate,initial,max_steps=3000,climb_after=100,tolerance=2e-6)
    assert result['converged']
    np.testing.assert_allclose(result['images'][result['highest_image']],[0.,.3],atol=2e-6)
    np.testing.assert_allclose(result['energies'][result['highest_image']],1.,atol=2e-11)
    np.testing.assert_array_equal(result['images'][[0,-1]],initial[[0,-1]])
