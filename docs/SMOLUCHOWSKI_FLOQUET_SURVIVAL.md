# Periodic first passage as a Floquet--Perron survival problem

## Status and scope

This result belongs to the active one-dimensional normal-tensile model. It
uses only the spacing coordinate $\lambda=a/a_0$, its intact density, the
normal-force waveform, and the already declared tangent-instability absorbing
boundary $\lambda_c$. It introduces no slip coordinate, plastic strain,
empirical damage variable, retained-energy fraction, fitted fatigue kernel, or
new physical coefficient.

The numerical values below are dimensionless consequences of an uncalibrated
demonstration. They are not a fatigue-life prediction for aluminum.

## 1. Periodic absorbing evolution

In reduced variables the intact density obeys

$$
\partial_t\rho=-\partial_\lambda J,\qquad
J=-\left[(\partial_\lambda\phi-f(t))\rho
+\beta^{-1}\partial_\lambda\rho\right],
$$

where $\beta=E_0/(k_BT)$, $f=Fa_0/E_0$, and time is scaled by the spacing
relaxation time

$$
t_r=\frac{a_0^2}{M_aE_0}=\frac{\gamma a_0^2}{E_0}.
$$

Thus the reduced period used below is $T_c/t_r$; its conversion cannot be made
until the spacing mobility is physically determined. The lower boundary
reflects and $\lambda_c$ absorbs. For a
periodic force $f(t+T)=f(t)$, let $\mathcal U(t_2,t_1)$ be the linear evolution
operator. The one-cycle operator is

$$
\boxed{\mathcal K=\mathcal U(T,0).}
$$

It is positive and sub-Markov: a nonnegative density remains nonnegative and
its total intact mass cannot increase. This operator is not an assumed
fatigue rule. It is the exact one-period propagation of the stated PDE (or of
its declared finite-volume discretization).

For one backward-Euler finite-volume substep,

$$
\rho_i^{j+1}-\rho_i^j
=-\frac{\Delta t}{\Delta\lambda}
\left(J_{i+1/2}^{j+1}-J_{i-1/2}^{j+1}\right).
$$

Summing over cells and using zero lower flux gives the discrete identity

$$
\boxed{S^{j+1}-S^j=-\Delta t\,J_{\rm out}^{j+1}.}
$$

The Chang--Cooper flux and positive diffusion make the implicit matrix an
irreducible M-matrix under the present finite-domain conditions. Its inverse
is nonnegative. A product of these substeps therefore yields a nonnegative,
mass-decreasing matrix $K$.

## 2. Principal survival multiplier

Perron--Frobenius theory applied to the finite-volume cycle matrix gives a
positive principal right eigenvector $q_0$ and multiplier $r$:

$$
\boxed{Kq_0=rq_0,\qquad q_0>0,\qquad
\int q_0\,d\lambda=1,\qquad 0<r<1.}
$$

$q_0$ is the conditional intact distribution at a chosen cycle phase. The
mass-normalized power iteration used by the code is simply

$$
\widetilde q_{j+1}=Kq_j,\qquad
r_j=\int\widetilde q_{j+1}\,d\lambda,\qquad
q_{j+1}=\widetilde q_{j+1}/r_j.
$$

No lifetime curve is fitted in this iteration.

Changing the chosen cycle origin cyclically permutes the invertible implicit
substep matrices. Cyclic products have the same nonzero eigenvalues, so $r$
must be phase-origin invariant although $q_0$ changes with phase. This identity
is now a regression test. It also exposed and corrected a timestep-counting
bug: floating-point values such as $2.000000000000001$ previously caused an
accidental extra substep in isolated output intervals.

Propagating the eigenstate inside a cycle defines

$$
g(\theta)=\int\mathcal U(\theta,0)q_0\,d\lambda,
\qquad
q(\lambda,\theta)=\frac{\mathcal U(\theta,0)q_0}{g(\theta)}.
$$

Then

$$
g(0)=1,\qquad g(T)=r,\qquad
\boxed{q(\lambda,T)=q(\lambda,0)}.
$$

Thus $q$ is a periodic conditional, or periodic quasi-stationary,
distribution. Because normalization removes the escaping mass, it obeys the
nonlinear conditioned equation

$$
\boxed{\partial_tq=-\partial_\lambda J[q]+h(t)q,\qquad
h(t)=\frac{J_{\rm out}[q](t)}{\int q\,d\lambda}=J_{\rm out}[q](t).}
$$

The $+hq$ term is not a source of intact material. It is the normalization
term required after conditioning on survival.

## 3. Exact cycle hazard and long-cycle survival

Since $\dot g=-J_{\rm out}$,

$$
h(t)=\frac{J_{\rm out}(t)}{g(t)}=-\frac{d}{dt}\log g(t).
$$

Integration through one period gives the central result

$$
\boxed{\mathcal H_c=\int_0^T h(t)\,dt=-\log r.}
$$

The three equivalent cycle-scale quantities are therefore

$$
\boxed{r=\text{survival per cycle},\qquad
1-r=\text{escape per cycle},\qquad
-\log r=\text{integrated cycle hazard}.}
$$

Starting exactly from $q_0$,

$$
\boxed{S_N=r^N.}
$$

If $K_{\rm init}\in\{1,2,\ldots\}$ is the cycle containing first passage,
then the principal-mode cycle distribution is not assumed Weibull; it follows
exactly from the Markov cycle map:

$$
\boxed{\Pr(K_{\rm init}=k)=r^{k-1}(1-r),\qquad
\mathbb E[K_{\rm init}]=\frac{1}{1-r},\qquad
\operatorname{Var}(K_{\rm init})=\frac{r}{(1-r)^2}.}
$$

This geometric law applies to cycle bins after the periodic conditional mode
has been reached. The phase-resolved outgoing flux retains the within-cycle
timing that the binned law discards.

For a generic nonnegative initial density, $S_{N+1}/S_N\to r$ after the
conditional shape approaches $q_0$. The prefactor before $r^N$ depends on the
initial density and is not silently set to one. Let $w_0$ be the positive left
Perron mode, normalized biorthogonally:

$$
w_0^TK=rw_0^T,\qquad
\int w_0(\lambda)q_0(\lambda)\,d\lambda=1.
$$

Then the full leading asymptotic law is

$$
\boxed{S_N\sim C(\rho_0)r^N,\qquad
C(\rho_0)=\int w_0(\lambda)\rho_0(\lambda)\,d\lambda.}
$$

$w_0$ is a survival-propensity weight of the starting coordinate, not a new
material field. If $r_2$ is the second eigenvalue in modulus, the leading
conditional transient contracts approximately as $(|r_2|/r)^N$. Thus both
the life multiplier and the number of cycles required before it is dominant
come from the same operator. In the principal mode only, a probability
quantile $p$ corresponds to

$$
N_p=\frac{\log(1-p)}{\log r}.
$$

This is a derived cycle count of the model, not a new damage law.

## 4. What accumulates, and what does not

The unconditioned intact mass decreases every cycle. That is the irreversible
quantity in the active model. In the principal regime, however, the
conditional density, its mean spacing, variance, and mean interaction energy
all repeat after each period:

$$
q(T)=q(0),\qquad
\langle U\rangle_{q(T)}=\langle U\rangle_{q(0)}.
$$

Therefore this Markov model does **not** predict indefinite cycle-to-cycle
storage of energy inside the survivors. Dissipated loop work is transferred
to the eliminated isothermal bath; it is not added to the pair potential.
The potential controls drift and barrier approach, and thereby the escape
flux, but dissipated work is not converted into an arbitrary tail source.

A continuing growth of the conditional tail would require additional
physical memory: for example a microscopically derived slow structural
coordinate, evolving bath state, or non-Markovian kernel. Adding such a state
without a derivation would contradict the present scope. Plastic slip is also
not produced by the normal-only coordinate.

## 5. Numerical verification

The committed implementation performs five independent checks:

1. nonnegativity and sub-Markov column sums of an explicitly assembled small
   cycle matrix;
2. agreement between its spectral radius and power iteration;
3. exact discrete equality between lost mass and integrated outgoing flux;
4. direct $N$-cycle propagation from $q_0$ versus $r^N$;
5. convergence of a generic initial density's cycle ratios to $r$, plus grid
   and timestep refinement; and
6. invariance of $r$ under quarter-cycle changes of the protocol origin.

The dense verification also computes $w_0$ and confirms that
$S_N/r^N\to C(\rho_0)$ for a generic initial density. A positive transpose
power iteration is used because direct left eigenvectors of the strongly
nonnormal absorbing matrix are poorly conditioned near the repulsive edge.

At the demonstration point

$$
m=12.19,\quad n=6,\quad \beta=2000,\quad
f(t)=0.008+0.007\sin(2\pi t/12),
$$

the refined run gives approximately

$$
r=0.9047876,\qquad 1-r=0.0952124,\qquad
\mathcal H_c=0.100055,\qquad N_{50}=6.928.
$$

The mass-balance residual is about $2\times10^{-15}$. Direct cycle ratios from
the original conditional-equilibrium initialization converge to $r$ by the
second to third cycle. On the 70-cell dense diagnostic, $|r_2|/r\simeq
9.12\times10^{-5}$ and the independently calculated initial-state coefficient
is $C\simeq0.98978$, quantitatively explaining the rapid transient collapse.
The low $N_{50}$ explains the earlier observation that
the demonstration fails at very few cycles: the chosen reduced inputs allow
roughly $9.5\%$ escape per cycle. It does not demonstrate that real
single-crystal aluminum has this lifetime.

The phase-resolved principal mode quantifies the hysteretic transport delay.
The tensile force peaks at phase $0.25$, the conditional mean spacing, mean
interaction energy and absorbing hazard peak at about $0.35$, and the variance
peaks at about $0.40$. These delays arise because finite mobility transports
the density over a nonzero relaxation time. Thermal diffusion fixes the width
and participates in the current, but diffusion by itself is not the cause of
the loading--unloading distinction.

Increasing frequency reduces the time available for first passage in each
cycle, so escape per cycle falls in this sweep. The mean hazard per unit
reduced time and the escape per cycle are distinct and are both reported.

This trend has two controlled operator limits. With phase
$\theta=t/T$, the PDE is $\partial_\theta\rho=T\mathcal L(\theta)\rho$.
For $T\to0$, the first Magnus term gives

$$
\mathcal K=I+T\overline{\mathcal L}+O(T^2),\qquad
1-r=O(T),
$$

so escape per cycle vanishes while $-\log(r)/T$ approaches the principal decay
rate of the phase-averaged generator. Numerically it changes only from about
$0.0071535$ to $0.0071561$ between $T=0.125$ and $0.25$. In the slow limit,
provided the instantaneous absorbing generator has a separated principal mode,
adiabatic following gives

$$
-\frac{\log r}{T}\longrightarrow
\int_0^1\kappa_0[f(\theta)]\,d\theta,
$$

where $\kappa_0$ is the instantaneous principal escape rate. This explains why
the sweep's mean hazard per reduced time approaches a finite plateau although
the probability lost during one increasingly long cycle approaches one.

## 6. Physical quantities still required for an aluminum prediction

Mapping this spectrum to a physical life requires independently justified
$E_0$, $A_0$, $T$, mobility or friction (hence the physical relaxation time),
the single-crystal loading-axis modulus, and validation of the operational
$\lambda_c$ initiation definition. The mechanical area $A_0$ is not a
correlation area or a FEM element area. No parameter in the present sweep was
fitted to an aluminum S--N curve.

The exact nondimensional groups show what can and cannot be identified. With
the loading-axis single-crystal modulus and cell calibration

$$E_0=E_{[hkl]}A_0a_0,$$

the reduced normal force simplifies to

$$
\boxed{f=\frac{Fa_0}{E_0}
=\frac{\sigma A_0a_0}{E_{[hkl]}A_0a_0}
=\frac{\sigma}{E_{[hkl]}}.}
$$

Thus $A_0$ cancels from the reduced stress input. It does not cancel from the
thermal and time groups:

$$
\boxed{\beta=\frac{E_{[hkl]}A_0a_0}{k_BT},\qquad
t_r=\frac{\gamma a_0}{E_{[hkl]}A_0},\qquad
T^*=\frac{T_{\rm physical}}{t_r}.}
$$

Consequently a measured loading-axis modulus and stress waveform are
insufficient to predict cycles to initiation. $A_0$ controls the thermal
coarse-graining and $\gamma$ controls the physical clock. Neither can be
replaced by FEM element area or tuned invisibly through the timestep. The
complete reduced multiplier has the identifiable structure

$$
r=\mathscr R\!\left(m,n,\frac{\sigma_m}{E_{[hkl]}},
\frac{\sigma_a}{E_{[hkl]}},\frac{1}{f_{\rm physical}t_r},
\frac{E_{[hkl]}A_0a_0}{k_BT};\lambda_c\right).
$$

## 한국어 요약

주기하중 아래 흡수형 Smoluchowski 방정식을 정확히 한 주기 적분하면
생존확률을 줄이는 선형 연산자 $K$가 생긴다. 여기에 별도 피로법칙을 붙인
것이 아니다. 최대 고유값 $r$은 장시간 한 cycle 생존비이고, $1-r$은 한
cycle 유출확률이며, $-\log r$은 한 cycle 누적 hazard다. 고유분포에서
$S_N=r^N$이 정확히 성립하고, 일반 초기분포도 transient 뒤에는 연속 cycle
생존비가 $r$로 수렴한다.

현재 무차원 예제에서는 $r\simeq0.9048$이므로 cycle당 약 $9.5\%$가
변곡점 흡수경계를 통과한다. 그래서 중앙 생존 cycle이 약 6.9회로 낮다.
이는 현재 선택한 무차원 온도, 하중, mobility 시간척도와 경계 정의의
결과이지 알루미늄 수명 예측이 아니다.

또 하나의 중요한 결론은 생존조건부 분포가 장기적으로 매 cycle 똑같이
돌아온다는 점이다. 따라서 현재 Markov 모델에서 생존 원자근방의 에너지가
cycle마다 계속 쌓여 tail을 영구히 키우지는 않는다. 누적되는 것은 경계를
통과해 빠져나간 확률이다. 계속되는 조건부 tail 성장을 표현하려면 실제
미시역학에서 유도한 느린 내부좌표나 memory kernel이 추가로 필요하며,
근거 없이 damage 또는 energy-storage 변수를 넣어서는 안 된다.

## References for the spectral framework

- N. Champagnat and D. Villemonais, *Quasi-limiting estimates for periodic
  absorbed Markov chains*, arXiv:2211.02706 (2022).
- D. Daners, *Existence and perturbation of principal eigenvalues for a
  periodic-parabolic problem*, Electronic Journal of Differential Equations,
  Conference 05, 51--67 (2000).
- N. Champagnat and D. Villemonais, *Uniform convergence of penalized
  time-inhomogeneous Markov processes*, ESAIM: Probability and Statistics 22,
  129--162 (2018), doi:10.1051/ps/2017022.
- E. Seneta, *Non-negative Matrices and Markov Chains*, 2nd ed., Springer
  (1981).
