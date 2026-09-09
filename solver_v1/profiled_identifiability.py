"""Local sensitivity on the actual exact-calibration equality manifold.

Coefficient normalization is explicit, not a material parameter. Inequality
cones are NOT replaced by equality constraints here, so this is an equality
tangent diagnostic, not confidence intervals or a global identifiability proof.
"""
import numpy as np
from scipy.linalg import null_space


def equality_tangent_sensitivity(matrix, observations, coefficients,
                                 shape_matrix_derivatives, *, exact_rows):
    """Differentiate y=M(k)c with M_E(k)c=target_E kept exactly fixed.

    Put D=diag(max(1,|c|)), A=M_E D (rows scaled by target tolerances).
    Coefficient tangent directions are D Z, A Z=0. For each shape derivative
    M_,j, choose the minimum scaled-norm particular correction
        A z_,j = -M_E,j c / scale_E.
    Then y_,j=M_,j c + M D z_,j. Another particular solution differs only by
    the retained coefficient tangent columns. Held-out rows are never fitted.
    Shape derivatives use their caller-declared coordinates (e.g. log decay).
    """
    M=np.asarray(matrix,float);c=np.asarray(coefficients,float)
    dM=np.asarray(shape_matrix_derivatives,float);exact=np.asarray(exact_rows,int)
    if (M.ndim!=2 or c.shape!=(M.shape[1],) or dM.ndim!=3 or dM.shape[1:]!=M.shape
            or len(observations)!=len(M) or exact.ndim!=1 or len(exact)==0
            or len(np.unique(exact))!=len(exact) or np.any(exact<0) or np.any(exact>=len(M))
            or any(np.any(~np.isfinite(x)) for x in (M,c,dM))):
        raise ValueError('finite matching observation/coefficient/shape derivatives required')
    scale=np.asarray([o.scale for o in observations],float)
    if np.any(~np.isfinite(scale)) or np.any(scale<=0):raise ValueError('positive target scales required')
    selected=np.array([i for i,o in enumerate(observations) if o.role=='fit' and i not in exact],int)
    if not len(selected):raise ValueError('non-exact fit observations required')
    D=np.maximum(1.,abs(c));A=M[exact]*D/scale[exact,None]
    Z=null_space(A);rank=M.shape[1]-Z.shape[1]
    coefficient_columns=M@(D[:,None]*Z)
    shape_columns=[];corrections=[]
    for derivative in dM:
        rhs=-(derivative[exact]@c)/scale[exact]
        dz=np.linalg.lstsq(A,rhs,rcond=None)[0]
        if np.linalg.norm(A@dz-rhs,np.inf)>1e-9*max(1.,np.linalg.norm(rhs,np.inf)):
            raise ValueError('shape derivative incompatible with exact target manifold')
        corrections.append(D*dz)
        shape_columns.append(derivative@c+M@(D*dz))
    all_columns=np.column_stack([coefficient_columns,np.array(shape_columns).T])
    jacobian=all_columns[selected]/scale[selected,None]
    _,singular,right=np.linalg.svd(jacobian,full_matrices=False)
    norms=np.linalg.norm(jacobian,axis=0)
    cosines=np.divide(jacobian.T@jacobian,norms[:,None]*norms[None,:],
                      out=np.full((len(norms),len(norms)),np.nan),
                      where=(norms[:,None]*norms[None,:])>0)
    return dict(jacobian=jacobian,selected_rows=selected,exact_rows=exact,
        coefficient_normalization=D,exact_matrix_rank=rank,
        coefficient_tangent_dimension=Z.shape[1],shape_dimension=len(dM),
        physical_coefficient_tangent=D[:,None]*Z,shape_coefficient_corrections=corrections,
        exact_derivative_residual=float(np.max(abs(all_columns[exact]/scale[exact,None]))),
        singular_values=singular,right_singular_vectors=right,column_cosines=cosines,
        rank=int(np.linalg.matrix_rank(jacobian)),
        condition_number=float(singular[0]/singular[-1]) if singular[-1]>0 else float('inf'),
        active_inequality_cones_imposed=False,statistical_confidence_claimed=False,
        coordinate_system='all declared exact rows eliminated; normalized coefficient tangent plus caller-declared shape coordinates')
