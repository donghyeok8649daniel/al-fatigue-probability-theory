# Candidate Local-Traction Symbol Index / 후보 Local-Traction 기호 Index

Status: **candidate-only symbol index** for `CANDIDATE_LOCAL_TRACTION_P0_PROPAGATOR.md`. These symbols are not promoted to the active exact finite-chain formulation.

### $\mu_a$

- English: effective local spacing inertia
- 한국어: 국소 층간거리 유효 관성
- Mathematical definition: inertia multiplying $\ddot a$ in $\mu_a\ddot a=A_0\sigma_n-dU/da$
- Physical meaning: reduced inertia of the chosen local relative-spacing coordinate
- Unit: kg for a particle-like representative unit; otherwise consistent generalized inertia for the selected area/unit
- Status: CANDIDATE PHYSICAL PARAMETER; must be derived from the kinetic unit, not fitted only to tune fatigue time

### $t_{0,a}$

$$
t_{0,a}=\sqrt{\frac{\mu_a a_0}{EA_0}}
$$

- English: local spacing time scale
- 한국어: 국소 층간거리 시간척도
- Unit: s
- Dependencies: $\mu_a,a_0,E,A_0$
- Status: CANDIDATE DEFINITION

### $F(\lambda,c,\tau\mid x)$

- English: local spacing/rate phase-space density
- 한국어: 국소 층간거리-속도 위상공간 밀도
- Mathematical definition: solution of the candidate local Liouville equation
- Physical meaning: distribution over local spacing and local spacing rate at material point $x$
- Scaling: $\lambda=a/a_0$, $c=d\lambda/d\tau$
- Status: CANDIDATE STATE

### $P_0(\lambda\mid x)$

$$
P_0(\lambda\mid x)=P(\lambda,0\mid x)
$$

- English: initial local spacing density
- 한국어: 초기 국소 층간거리 밀도
- Physical meaning: measured or independently derived initial population of local spacings at material point $x$
- Status: REQUIRED INPUT; physical construction remains open

### $\Lambda(\tau;\lambda_0,x)$

$$
\frac{d\Lambda}{d\tau}=C
$$

$$
\frac{dC}{d\tau}=q(x,\tau)-\phi'(\Lambda)
$$

- English: local spacing characteristic
- 한국어: 국소 층간거리 특성곡선
- Initial condition: $\Lambda(0;\lambda_0,x)=\lambda_0$
- Physical meaning: trajectory of one initial-spacing label under the prescribed local stress history
- Status: CANDIDATE EXACT CHARACTERISTIC within the reduced model

### $C(\tau;\lambda_0,x)$

- English: characteristic spacing rate
- 한국어: 특성곡선 층간거리 속도
- Mathematical definition: $C=d\Lambda/d\tau$
- Initial static-preparation condition: $C(0;\lambda_0,x)=0$
- Status: CANDIDATE STATE

### $\mathcal T_{\tau,0}^{q}$

$$
P(\tau)=\mathcal T_{\tau,0}^{q}[P_0]
$$

- English: local stress-history probability propagator
- 한국어: 국소 응력이력 확률 전파연산자
- Mathematical meaning: push-forward of $P_0$ through $\Lambda(\tau;\lambda_0,x)$
- Physical meaning: candidate map from initial spacing density and supplied stress history to later spacing density
- Status: CANDIDATE REDUCED PROPAGATOR

### $e^*$

$$
e^*=\frac12C^2+\phi(\Lambda)
$$

- English: characteristic intrinsic mechanical energy
- 한국어: 특성곡선 고유 기계에너지
- Mathematical identity: $de^*/d\tau=qC$
- Status: EXACT within candidate local model
