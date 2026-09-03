# Milestone 25 — Exact Integral Representation / 정확한 적분 표현

이 문서는 active 1D $P$-$u$-$\Theta$ 이론의 exact integral representation을 정리한다.

핵심 결론은 다음과 같다.

- reduced PDE hierarchy가 autonomous하게 닫히지 않아도 exact integral representation은 존재한다.
- higher conditional moments는 임의의 phenomenological coefficient가 아니라 full deterministic flow의 projection으로 정의할 수 있다.
- specimen survival도 initial full-state measure가 주어지면 path integral로 쓸 수 있다.

## 1. Full deterministic state / 전체 결정론적 상태

finite chain의 full microscopic state를

$$
\Gamma=(x_1,\ldots,x_M,\dot x_1,\ldots,\dot x_M)
$$

로 둔다.

초기상태는

$$
\Gamma(\tau_0)=\Gamma_0
$$

이다.

주어진 forcing $q$ 아래 deterministic flow map을

$$
\Gamma(\tau)=\Phi^q_{\tau,\tau_0}(\Gamma_0)
$$

로 정의한다.

초기 full-state measure는

$$
\int\mu_0(d\Gamma_0)=1
$$

을 만족한다.

한 개의 deterministic realization만 사용할 경우 $\mu_0$는 한 점에 집중된 measure로 둘 수 있다. 이 정의에는 Boltzmann 또는 Gibbs 가정이 필요하지 않다.

## 2. Trajectory projections / 궤적 투영

spacing trajectory를

$$
\Lambda_i(\tau;\Gamma_0)=x_i(\tau;\Gamma_0)-x_{i-1}(\tau;\Gamma_0)
$$

로 정의한다.

spacing-rate trajectory는

$$
C_i(\tau;\Gamma_0)=\frac{d\Lambda_i}{d\tau}
$$

이고 spacing-acceleration trajectory는

$$
A_i(\tau;\Gamma_0)=\frac{d^2\Lambda_i}{d\tau^2}
$$

이다.

## 3. Exact push-forward phase-space density / 정확한 push-forward 위상공간 밀도

full deterministic flow를 spacing-rate phase space로 투영하면

$$
F(\lambda,c,\tau)=\frac{1}{M}\sum_{i=1}^{M}\int\delta(\lambda-\Lambda_i(\tau;\Gamma_0))\delta(c-C_i(\tau;\Gamma_0))\,\mu_0(d\Gamma_0)
$$

이다.

따라서 spacing marginal은

$$
P(\lambda,\tau)=\frac{1}{M}\sum_{i=1}^{M}\int\delta(\lambda-\Lambda_i(\tau;\Gamma_0))\,\mu_0(d\Gamma_0)
$$

이다.

## 4. Exact integral representation of $u$ and $\Theta$

첫 conditional rate moment는

$$
P(\lambda,\tau)u(\lambda,\tau)=\frac{1}{M}\sum_i\int C_i(\tau;\Gamma_0)\delta(\lambda-\Lambda_i(\tau;\Gamma_0))\,\mu_0(d\Gamma_0)
$$

이다.

second raw conditional moment는

$$
P(\lambda,\tau)[u(\lambda,\tau)^2+\Theta(\lambda,\tau)]=\frac{1}{M}\sum_i\int C_i(\tau;\Gamma_0)^2\delta(\lambda-\Lambda_i(\tau;\Gamma_0))\,\mu_0(d\Gamma_0)
$$

이다.

따라서 $P>0$인 곳에서

$$
u(\lambda,\tau)=\frac{\frac{1}{M}\sum_i\int C_i\delta(\lambda-\Lambda_i)\,\mu_0(d\Gamma_0)}{P(\lambda,\tau)}
$$

이다.

또

$$
\Theta(\lambda,\tau)=\frac{\frac{1}{M}\sum_i\int C_i^2\delta(\lambda-\Lambda_i)\,\mu_0(d\Gamma_0)}{P(\lambda,\tau)}-u(\lambda,\tau)^2
$$

이다.

## 5. Exact integral representation of higher terms / 상위항의 정확한 적분표현

one-point conditional acceleration은

$$
P\mathcal A=\frac{1}{M}\sum_i\int A_i\delta(\lambda-\Lambda_i)\,\mu_0(d\Gamma_0)
$$

로 쓸 수 있다.

third conditional central moment는

$$
PC_3=\frac{1}{M}\sum_i\int(C_i-u)^3\delta(\lambda-\Lambda_i)\,\mu_0(d\Gamma_0)
$$

이다.

rate-acceleration covariance source는

$$
P\Psi=\frac{1}{M}\sum_i\int(C_i-u)A_i\delta(\lambda-\Lambda_i)\,\mu_0(d\Gamma_0)
$$

이다.

따라서 $C_3$와 $\Psi$는 arbitrary fitting constants가 아니라 full mechanics에서 정의된 projection이다.

## 6. Volterra form / 볼테라 적분형

continuity equation을 시간에 대해 적분하면

$$
P(\lambda,\tau)=P_0(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(\lambda,s)u(\lambda,s)\,ds
$$

이다.

first moment balance는

$$
P(\lambda,\tau)u(\lambda,\tau)=P_0(\lambda)u_0(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(u^2+\Theta)(\lambda,s)\,ds+\int_{\tau_0}^{\tau}P(\lambda,s)\mathcal A(\lambda,s)\,ds
$$

이다.

second raw moment는

$$
P(u^2+\Theta)(\lambda,\tau)=P_0(u_0^2+\Theta_0)(\lambda)-\partial_\lambda\int_{\tau_0}^{\tau}P(u^3+3u\Theta+C_3)(\lambda,s)\,ds+2\int_{\tau_0}^{\tau}P(u\mathcal A+\Psi)(\lambda,s)\,ds
$$

이다.

이 식들은 exact hierarchical Volterra equations다.

## 7. Characteristic form for $P$ / $P$의 특성곡선형

characteristic을

$$
\frac{dX}{ds}=u(X(s),s)
$$

$$
X(\tau_0)=\alpha
$$

로 정의한다.

누적 mean-rate gradient는

$$
\mathcal I_u(s;\alpha)=\int_{\tau_0}^{s}\partial_\lambda u(X(r),r)\,dr
$$

이다.

continuity equation은 characteristic 위에서

$$
P(X(\tau),\tau)=P_0(\alpha)e^{-\mathcal I_u(\tau;\alpha)}
$$

를 준다.

## 8. Characteristic form for $u$ / $u$의 특성곡선형

mean-flow equation으로부터

$$
u(X(\tau),\tau)=u_0(\alpha)+\int_{\tau_0}^{\tau}\left[\mathcal A-\frac{1}{P}\partial_\lambda(P\Theta)\right]_{(X(s),s)}\,ds
$$

이다.

이 식은 closure가 아니라 exact integral representation이다. 오른쪽의 fields는 full mechanics 또는 hierarchy에서 공급되어야 한다.

## 9. Characteristic form for $\Theta$ / $\Theta$의 특성곡선형

variance equation은

$$
D_\tau\Theta+2\Theta\partial_\lambda u+\frac{1}{P}\partial_\lambda(PC_3)=2\Psi
$$

이다.

source를

$$
S_\Theta=2\Psi-\frac{1}{P}\partial_\lambda(PC_3)
$$

로 두면 characteristic 위에서

$$
\frac{d\Theta}{ds}+2\Theta\partial_\lambda u=S_\Theta
$$

이다.

integrating factor를 적용하면

$$
\Theta(X(\tau),\tau)=e^{-2\mathcal I_u(\tau;\alpha)}\left(\Theta_0(\alpha)+\int_{\tau_0}^{\tau}e^{2\mathcal I_u(s;\alpha)}S_\Theta(X(s),s)\,ds\right)
$$

를 얻는다.

## 10. Exact G1 and G2 push-forward forms / G1·G2 push-forward형

평균 spacing은

$$
\bar a(\tau)=\frac{a_0}{M}\sum_i\int\Lambda_i(\tau;\Gamma_0)\,\mu_0(d\Gamma_0)
$$

이다.

평균 intrinsic configurational energy는

$$
\bar U(\tau)=\frac{U_{\mathrm{ref}}}{M}\sum_i\int[\phi(\Lambda_i(\tau;\Gamma_0))-\phi(1)]\,\mu_0(d\Gamma_0)
$$

이다.

## 11. Exact local first-passage survival / 정확한 국소 최초통과 생존

각 trajectory의 first-passage event를 path functional로 직접 쓸 수 있다.

$$
S_{\mathrm{local}}(\tau)=\frac{1}{M}\sum_i\int I\left[\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c\right]\,\mu_0(d\Gamma_0)
$$

따라서

$$
F_{\mathrm{ci}}^{\mathrm{local}}(\tau)=1-S_{\mathrm{local}}(\tau)
$$

이다.

## 12. Exact specimen survival / 정확한 시편 생존

시편 전체가 아직 initiation되지 않았다는 사건은 모든 represented spacing이 threshold 아래에 머무르는 사건이다.

$$
S_{\mathrm{spec}}(\tau)=\int I\left[\max_i\sup_{s\in[\tau_0,\tau]}\Lambda_i(s;\Gamma_0)<\lambda_c\right]\,\mu_0(d\Gamma_0)
$$

따라서

$$
F_{\mathrm{ci}}^{\mathrm{spec}}(\tau)=1-S_{\mathrm{spec}}(\tau)
$$

이다.

이 식에는 independent-cell product가 필요하지 않다. spatial correlation은 $\mu_0$와 full deterministic trajectories에 포함될 수 있다.

## 13. Meaning of the remaining open problem / 남은 미해결의 의미

수학적으로 missing formula가 있는 것이 아니라 실제 재료의 initial full-state measure와 correlation scale을 어떻게 정의하고 calibration할지가 OPEN이다.

$$
\mathrm{exact\ integral\ representation}\ne\mathrm{autonomous\ low\ order\ closure}
$$

두 개념을 구분해야 한다.
