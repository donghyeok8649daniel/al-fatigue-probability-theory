# Normal-LJ Quasistatic Protocol Diagnostic

## Classification

**EXACT STATIC RESULT + CONTROLLED NUMERICAL DYNAMICAL DIAGNOSTIC**

For the homogeneous open chain under a constant tensile end force $f$,

$$
\Pi=\sum_{i=1}^{M}[\phi(\lambda_i)-f\lambda_i].
$$

On the stable branch $\phi''>0$, stationarity gives $\phi'(\lambda_i)=f$ for every spacing, hence all spacings equal the unique stable root $\lambda_s(f)$. Therefore the exact zero-temperature quasistatic empirical distribution is

$$
P_M^{\rm qs}(\lambda\mid f)=\delta[\lambda-\lambda_s(f)].
$$

The Milestone 13 snapshots are taken at integer cycle 2 for a zero-mean sine load. At that exact phase the applied force is zero, so the quasistatic reference is $\lambda_i=1$ and $C_0=0$. Any nonzero variance in those snapshots is therefore dynamical residual structure, not static material randomness.

## Numerical sweep

The sweep uses force amplitude `0.03`, sample cycle `2`, and decreases the protocol parameter $\alpha=\omega M$ while keeping the model strictly one-dimensional.

| M | alpha=omega M | C0 | sqrt(C0) | |mean-1| | rho1 | M_eff^(+) |
|---:|---:|---:|---:|---:|---:|---:|
| 31 | 0.62 | 2.564064e-05 | 5.063659e-03 | 1.140016e-02 | 0.933406 | 2.929680 |
| 31 | 0.31 | 1.467564e-06 | 1.211430e-03 | 2.701704e-03 | 0.934132 | 2.928631 |
| 31 | 0.155 | 7.643812e-09 | 8.742890e-05 | 1.980035e-04 | 0.930860 | 2.961337 |
| 31 | 0.0775 | 2.911425e-10 | 1.706290e-05 | 3.991731e-05 | 0.929173 | 2.978992 |
| 63 | 0.62 | 1.748917e-05 | 4.182005e-03 | 9.642913e-03 | 0.966267 | 2.986395 |
| 63 | 0.31 | 1.180593e-06 | 1.086551e-03 | 2.366335e-03 | 0.968346 | 2.926300 |
| 63 | 0.155 | 8.351185e-09 | 9.138482e-05 | 2.032952e-04 | 0.966781 | 2.956018 |
| 63 | 0.0775 | 2.964282e-10 | 1.721709e-05 | 3.915505e-05 | 0.965757 | 2.976604 |

The key diagnostic is that the fluctuation amplitude collapses toward the exact quasistatic state as $\omega M$ is reduced, while the normalized correlation shape and its positive-window effective count remain of order three. A normalized correlation length can therefore remain apparently system-scale even while the field it normalizes is disappearing.

This **supersedes the interpretation**, not the arithmetic, of the earlier $M_{\rm eff}^{(+)}\approx3$ result. That number remains a valid shape diagnostic for the selected deterministic residual snapshot, but it is not evidence for a finite material statistical-cell count or a material correlation length.

## Consequence for P(lambda,t)

A nontrivial fatigue probability distribution cannot be obtained from the adiabatic limit of one perfectly homogeneous deterministic zero-temperature chain alone. A physically broad $P$ requires a justified ensemble source, for example finite-temperature microstates, physically specified initial-condition uncertainty, or independently justified material heterogeneity. This does not authorize an arbitrary fitted distribution.

The next 1D target is therefore to separate

$$
P_{\rm spatial}^{\rm traj}(\lambda,t)=\frac1M\sum_i\delta(\lambda-\lambda_i(t))
$$

from an ensemble-averaged physical probability state and test the latter under cyclic driving.

---

# Normal-LJ 준정적 프로토콜 진단

## 분류

**정확한 정적 결과 + 통제된 수치 동역학 진단**

균질한 open chain에 일정한 인장 end force $f$를 가하면 spacing 좌표에서

$$
\Pi=\sum_{i=1}^{M}[\phi(\lambda_i)-f\lambda_i]
$$

이다. 안정 branch에서는 $\phi''>0$이므로 정지조건 $\phi'(\lambda_i)=f$의 안정해가 유일하고 모든 spacing은 같은 $\lambda_s(f)$가 된다. 따라서 zero-temperature 준정적 empirical distribution은

$$
P_M^{\rm qs}(\lambda\mid f)=\delta[\lambda-\lambda_s(f)]
$$

이다.

Milestone 13 snapshot은 zero-mean sine load의 정수주기 cycle 2에서 저장된다. 이 정확한 위상에서 applied force는 0이므로 준정적 기준은 $\lambda_i=1$, $C_0=0$이다. 따라서 기존 snapshot의 nonzero variance는 정적 물질 확률분포가 아니라 동적 잔류구조다.

## 수치 sweep

위 표와 동일한 수치에서 $\alpha=\omega M$을 낮추면 variance와 mean offset은 준정적 값 0으로 급격히 감소하지만, normalized correlation shape와 positive-window $M_{\rm eff}^{(+)}$는 약 3 수준을 유지한다.

따라서 기존 $M_{\rm eff}^{(+)}\approx3$의 **계산 자체가 틀린 것은 아니지만 해석은 수정해야 한다**. 그것은 선택한 deterministic residual snapshot의 normalized shape 진단값이지, 물질 고유 통계셀 개수나 물질 고유 correlation length의 증거가 아니다.

## P(lambda,t)에 대한 결과

완전히 균질한 deterministic zero-temperature chain의 adiabatic limit만으로는 폭을 가진 fatigue probability distribution이 생기지 않는다. 물리적인 broad $P$를 만들려면 finite-temperature microstate, 물리적으로 정의된 initial-condition ensemble, 또는 독립적으로 정당화된 material heterogeneity 같은 ensemble source가 필요하다. 그렇다고 임의의 분포 fitting을 허용하는 것은 아니다.

다음 1D 목표는 한 trajectory의 spatial empirical measure와 실제 ensemble-averaged probability state를 명확히 분리하고, 후자를 cyclic loading에서 검증하는 것이다.
