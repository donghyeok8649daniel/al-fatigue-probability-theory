import numpy as np
from scipy.optimize import minimize
from .joint_shape_epigraph import constraints_and_jacobian


def test_joint_jacobian_includes_shape_coefficients_and_epigraph_signs():
    targets=np.array([1.,2.,3.]); scales=np.array([.2,.3,.4]); D=np.array([2.,.5])
    def evaluate(z):
        x=z[0]; M=np.array([[x*x,1],[np.exp(x),2],[1,x]])
        dM=np.array([[[2*x,0],[np.exp(x),0],[0,1]]])
        return constraints_and_jacobian(M,dM,targets,scales,[0],[1,2],z[1:3]*D,D,z[-1])
    z=np.array([.3,.5,.7,2.])
    a=evaluate(z); analytic=np.vstack([a[2],a[3]])
    numerical=[]
    for i in range(4):
        delta=np.eye(4)[i]*1e-6
        plus=evaluate(z+delta); minus=evaluate(z-delta)
        numerical.append((np.r_[plus[0],plus[1]]-np.r_[minus[0],minus[1]])/2e-6)
    np.testing.assert_allclose(analytic,np.array(numerical).T,rtol=2e-8,atol=2e-8)


def test_smooth_epigraph_solves_absolute_value_cusp_without_profile_derivative():
    # exact c=1; min |x*c| has a cusp at optimum x=0.
    def ev(z):
        return constraints_and_jacobian(np.array([[1.],[z[0]]]),np.array([[[0.],[1.]]]),
            np.array([1.,0.]),np.ones(2),[0],[1],np.array([z[1]]),np.ones(1),z[2])
    r=minimize(lambda z:z[2],[.8,1.,.8],jac=lambda z:np.array([0.,0.,1.]),method='SLSQP',
        bounds=[(-1,1),(0,None),(0,None)],constraints=[
            dict(type='eq',fun=lambda z:ev(z)[0],jac=lambda z:ev(z)[2]),
            dict(type='ineq',fun=lambda z:ev(z)[1],jac=lambda z:ev(z)[3])],
        options=dict(ftol=1e-12,maxiter=30))
    assert r.success
    np.testing.assert_allclose(r.x,[0,1,0],atol=1e-10)
