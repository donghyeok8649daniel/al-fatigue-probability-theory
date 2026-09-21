"""Independent analytic free energy checks on the thermal mean-force integral."""
import numpy as np
from results.silicon_wafer_feasibility.analyze_conditional_ensemble import profile_with_static_control


def test_known_static_energy_control_and_grid_refinement():
    # Exact conditional integral of U(q,z)=(q²-1)²+(2+q²) z²/2 at kBT=.1:
    # F(q)-F(0)=q^4-2q²+.05 log[(2+q²)/2].
    target=np.linspace(-1.5,1.5,1001)
    exact=target**4-2*target**2+.05*np.log((2+target*target)/2)
    errors=[]
    for n in [13,25,49]:
        q=np.linspace(-1.5,1.5,n)
        u=(q*q-1)**2;du=4*q*(q*q-1)
        force=du+.1*q/(2+q*q)
        report,free,interpolated_force=profile_with_static_control(q,force,u,du,0.)
        np.testing.assert_allclose(interpolated_force(q),force,atol=3e-15)
        assert len(report['stationary']) == 3
        assert abs(report['stationary'][1]['gap_A']) < 1e-12
        errors.append(float(np.max(abs(free(target)-exact))))
    assert errors[0]/errors[1]>10
    assert errors[1]/errors[2]>10
    assert errors[-1]<3e-6
