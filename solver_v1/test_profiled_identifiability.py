import numpy as np
import pytest

from .profiled_identifiability import equality_tangent_sensitivity
from .vector_material_calibration import MaterialObservation


def test_exact_shape_correction_preserves_all_constraints_and_excludes_holdout():
    # y0=c0+k*c1 exact; y1=c2 exact. There is only ONE free coefficient.
    k=.7;c=np.array([1.,2.,3.])
    M=np.array([[1.,k,0.],[0.,0.,1.],[0.,1.,2.],[2.,-1.,1.],[9.,8.,7.]])
    dM=np.zeros((1,5,3));dM[0,0,1]=1.
    obs=[MaterialObservation(str(i),y,1.,'test','exact' if i<2 else 'heldout' if i==4 else 'fit')
         for i,y in enumerate(M@c)]
    value=equality_tangent_sensitivity(M,obs,c,dM,exact_rows=[0,1])
    assert value['exact_matrix_rank']==2 and value['coefficient_tangent_dimension']==1
    assert value['selected_rows'].tolist()==[2,3]
    assert value['exact_derivative_residual']<1e-13
    np.testing.assert_allclose(M[:2]@value['physical_coefficient_tangent'],0.,atol=1e-14)
    step=1e-5;dc=np.array(value['shape_coefficient_corrections'][0])
    predictions=[]
    for sign in (-1,1):
        moved=M+sign*step*dM[0]
        coefficients=c+sign*step*dc
        # Exact nonlinear correction, whose first derivative is the computed dc.
        correction=np.linalg.lstsq(moved[:2],M[:2]@c-moved[:2]@coefficients,rcond=None)[0]
        predictions.append(moved@(coefficients+correction))
    fd=(predictions[1]-predictions[0])/(2*step)
    np.testing.assert_allclose(value['jacobian'][:,-1],fd[[2,3]],rtol=1e-8,atol=1e-10)
    M[4]*=1e9
    second=equality_tangent_sensitivity(M,obs,c,dM,exact_rows=[0,1])
    np.testing.assert_array_equal(value['jacobian'],second['jacobian'])


def test_shape_cannot_silently_leave_overdetermined_exact_manifold():
    M=np.array([[1.,0.],[1.,0.],[0.,1.]])
    obs=[MaterialObservation(str(i),1.,1.,'test','exact' if i<2 else 'fit') for i in range(3)]
    dM=np.zeros((1,3,2));dM[0,0,0]=1.
    with pytest.raises(ValueError,match='incompatible'):
        equality_tangent_sensitivity(M,obs,np.ones(2),dM,exact_rows=[0,1])
