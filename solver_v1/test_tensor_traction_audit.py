import numpy as np

from .run_tensor_traction_audit import declared_tensor_scenarios
from .interface_static_scenarios import resolved_tensor_tractions


def test_declared_cubic_scenarios_keep_both_shears_and_work_signs():
    cases,n,m,m2=declared_tensor_scenarios()
    np.testing.assert_allclose(np.cross(n,m),m2,atol=1e-15)
    for name,target in [('pure_interface_normal50MPa',[50,0,0]),
                        ('pure_interface_shear1_4MPa',[0,4,0]),
                        ('pure_interface_shear2_4MPa',[0,0,4])]:
        np.testing.assert_allclose(resolved_tensor_tractions(cases[name],n,m),target,atol=1e-14)
    basis=np.array([n,m,m2]);dq=np.array([.002,-.003,.004])
    for tensor in cases.values():
        local=resolved_tensor_tractions(tensor,n,m)
        np.testing.assert_allclose(local@dq,(tensor@n)@(basis.T@dq),atol=1e-15)
        np.testing.assert_allclose(resolved_tensor_tractions(-tensor,n,m),-local,atol=1e-15)
