import numpy as np
from scipy.special import logsumexp
from solver_v1.silicon_conditional_research import KB_EV_K,AMU_EV_PS2_A2,SI_MASS_AMU
from solver_v1.silicon_quantum_diagnostic import HBAR_EV_PS,quantum_harmonic_difference


def test_oscillator_free_energy_matches_explicit_energy_level_partition():
    frequencies=np.array([20.,80.]);reference=np.array([35.,45.])
    mass=AMU_EV_PS2_A2*SI_MASS_AMU
    temperature=300.;thermal=KB_EV_K*temperature
    def logz(omega):
        levels=np.arange(5000)+.5
        return sum(logsumexp(-HBAR_EV_PS*w*levels/thermal) for w in omega)
    exact=-thermal*(logz(frequencies)-logz(reference))
    np.testing.assert_allclose(quantum_harmonic_difference(mass*frequencies**2,mass*reference**2,temperature),exact,atol=3e-16)
    np.testing.assert_allclose(quantum_harmonic_difference(mass*frequencies**2,mass*reference**2,0.),
        .5*HBAR_EV_PS*np.sum(frequencies-reference),atol=3e-17)


def test_high_temperature_limit_is_same_measure_classical_logdet():
    values=np.array([.2,1.5,8.]);reference=np.array([.3,1.,10.])
    temperature=1e7
    expected=.5*KB_EV_K*temperature*np.log(values/reference).sum()
    np.testing.assert_allclose(quantum_harmonic_difference(values,reference,temperature),expected,atol=8e-8)
