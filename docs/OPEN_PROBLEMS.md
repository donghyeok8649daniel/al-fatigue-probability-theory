# Open Problems — Active 1D Layer-LJ Mainline

## Scope

The active theory is strictly one-dimensional and normal-only. The represented microscopic coordinate is the normal separation between material layers. The effective layer interaction is the calibrated generalized Lennard-Jones model.

Archived three-dimensional FCC and shear work are not part of the active derivation.

## 1. Validate or falsify the derived distribution closure

The current large-$M$ candidate is

$$
\boxed{
p_\lambda(\lambda,t)
=
Z^{-1}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)].
}
$$

The highest-priority test is to compare it against deterministic 1D layer-LJ simulations at the same measured

$$
\mu(t)
$$

and

$$
\mathcal E(t).
$$

The equal-measure fixed-$(L,E)$ assumption must be rejected if the empirical distribution systematically disagrees.

## 2. Derive the time law for stored configurational energy

The closure converts

$$
\mu(t),\mathcal E(t)
$$

into a distribution, but it does not yet determine $\mathcal E(t)$.

The next mechanics problem is

$$
\boxed{
\sigma_n(t)
\rightarrow
W(t)
\rightarrow
\mathcal E(t)
}
$$

using an explicit energy balance, without a fitted retained-energy fraction.

Kinetic energy, reversible mean deformation energy, and energy stored in distributional broadening must be separated.

## 3. Determine whether external stress supplies an additional independent moment

A candidate normal force constraint is

$$
\frac{1}{A_0}
\int U'(a)P_a(a,t)\,da
=\sigma_n(t).
$$

Its exact validity depends on the precise reduced layer model and boundary conditions. It should first be derived from the 1D mechanics and then used as an independent test of the two-moment closure.

If it is genuinely independent, the distribution family may need an additional conjugate multiplier rather than silently forcing the two-moment form to fit it.

## 4. Finite-$M$ corrections

The exponential form is a large-$M$ saddle-point result. The exact finite-$M$ marginal under the stated ensemble is

$$
p_M(\lambda\mid L,E)
=
\frac{\Omega_{M-1}(L-\lambda,E-\psi(\lambda))}{\Omega_M(L,E)}.
$$

The magnitude of the saddle-point error as a function of represented layer count must be quantified.

## 5. Tail versus first passage

The instantaneous tail

$$
Q_c(t)
=
\int_{\lambda_c}^{\infty}p_\lambda(\lambda,t)\,d\lambda
$$

is not automatically the cumulative crack-initiation probability.

A later theory must connect this instantaneous distribution to a first-passage event

$$
\tau_c
=
\inf\{t:\text{mechanically defined normal-opening instability occurs}\}.
$$

## 6. Exact status of the compression side

The earlier exact feasibility theorem shows that normalization, mean, and energy alone cannot force a tensile tail because the LJ repulsive side can carry unbounded energy as $\lambda\to0^+$.

The new distribution closure selects a particular entropy-dominant state rather than solving that exact worst-case feasibility problem.

These are complementary statements and must not be confused:

- the exact feasibility bound asks what **any** admissible distribution could do;
- the saddle-point closure predicts which distribution is selected under the additional ensemble assumption.

## 7. Required falsification tests

Any next 1D model must satisfy:

1. normalization of $P$;
2. positivity of $P$;
3. exact recovery of the imposed mean and energy moments;
4. convergence with quadrature/grid refinement;
5. no fitted Gaussian, Weibull, or fatigue-damage law;
6. fixed LJ parameters through time;
7. explicit separation between exact results and closure assumptions;
8. direct comparison against microscopic 1D dynamics whenever a distribution form is proposed;
9. no claim that $\beta=1/(k_BT)$ without an equilibrium derivation;
10. no claim that instantaneous $Q_c$ equals cumulative crack-initiation probability without a first-passage derivation.

---

# 한국어 번역 — 활성 1D Layer-LJ Mainline 미해결 문제

## 범위

활성 이론은 엄격하게 1차원 수직변형만 다룬다. represented microscopic coordinate는 material layer 사이의 수직간격이며 layer 간 유효상호작용은 calibration된 generalized Lennard-Jones model이다.

archive된 3D FCC와 shear 연구는 active derivation에 포함하지 않는다.

## 1. 유도된 distribution closure 검증 또는 반증

현재 large-$M$ 후보는

$$
\boxed{
p_\lambda(\lambda,t)
=
Z^{-1}
\exp[-\alpha(t)\lambda-\beta(t)\psi(\lambda)]
}
$$

이다.

가장 우선순위가 높은 시험은 deterministic 1D layer-LJ simulation에서 측정한 동일한

$$
\mu(t)
$$

및

$$
\mathcal E(t)
$$

를 넣었을 때 empirical distribution과 이 closure를 비교하는 것이다.

empirical distribution이 체계적으로 다르면 equal-measure fixed-$(L,E)$ assumption은 기각해야 한다.

## 2. Stored configurational energy의 시간법칙 유도

현재 closure는

$$
\mu(t),\mathcal E(t)
$$

로부터 distribution을 계산하지만 $\mathcal E(t)$ 자체를 결정하지는 않는다.

다음 mechanics 문제는 fitted retained-energy fraction 없이 explicit energy balance로

$$
\boxed{
\sigma_n(t)
\rightarrow
W(t)
\rightarrow
\mathcal E(t)
}
$$

를 유도하는 것이다.

kinetic energy, reversible mean deformation energy, distributional broadening에 저장된 energy를 분리해야 한다.

## 3. External stress가 독립적인 추가 moment를 주는가

candidate normal force constraint는

$$
\frac{1}{A_0}
\int U'(a)P_a(a,t)\,da
=\sigma_n(t)
$$

이다.

정확한 유효성은 precise reduced layer model과 boundary condition에 의존한다. 먼저 1D mechanics에서 유도한 뒤 two-moment closure에 대한 independent test로 사용해야 한다.

실제로 독립적인 constraint라면 two-moment form을 억지로 맞추는 대신 additional conjugate multiplier가 필요한지 검토해야 한다.

## 4. Finite-$M$ correction

exponential form은 large-$M$ saddle-point result다. stated ensemble 아래 exact finite-$M$ marginal은

$$
p_M(\lambda\mid L,E)
=
\frac{\Omega_{M-1}(L-\lambda,E-\psi(\lambda))}{\Omega_M(L,E)}
$$

이다.

represented layer count에 따라 saddle-point error가 얼마나 되는지 정량화해야 한다.

## 5. Tail과 first passage 구분

instantaneous tail

$$
Q_c(t)
=
\int_{\lambda_c}^{\infty}p_\lambda(\lambda,t)\,d\lambda
$$

은 자동으로 cumulative crack-initiation probability가 아니다.

나중 이론에서는 이 instantaneous distribution을 first-passage event

$$
\tau_c
=
\inf\{t:\text{mechanically defined normal-opening instability occurs}\}
$$

와 연결해야 한다.

## 6. Compression side의 정확한 상태

기존 exact feasibility theorem은 normalization, mean, energy만으로 tensile tail을 강제할 수 없다는 것을 보였다. LJ repulsive side는 $\lambda\to0^+$에서 무한한 energy를 담을 수 있기 때문이다.

새 distribution closure는 그 exact worst-case feasibility problem을 푸는 대신 additional ensemble assumption 아래 entropy-dominant state 하나를 선택한다.

두 명제는 서로 보완적이며 혼동하면 안 된다.

- exact feasibility bound는 **어떤** admissible distribution이라도 할 수 있는지를 묻는다.
- saddle-point closure는 additional ensemble assumption 아래 어떤 distribution이 선택되는지를 예측한다.

## 7. 필수 반증시험

다음 1D model은 최소한 다음을 만족해야 한다.

1. $P$ normalization;
2. $P$ positivity;
3. imposed mean/energy moment의 정확한 복원;
4. quadrature/grid refinement convergence;
5. fitted Gaussian, Weibull, fatigue-damage law 금지;
6. 시간에 따른 LJ parameter 고정;
7. exact result와 closure assumption 명시적 분리;
8. distribution form 제안 시 microscopic 1D dynamics와 직접 비교;
9. equilibrium derivation 없이 $\beta=1/(k_BT)$라고 주장하지 않기;
10. first-passage derivation 없이 instantaneous $Q_c$를 cumulative crack-initiation probability라고 주장하지 않기.
