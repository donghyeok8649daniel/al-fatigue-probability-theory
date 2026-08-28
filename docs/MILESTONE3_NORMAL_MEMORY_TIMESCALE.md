# Milestone 3 — Normal-LJ time-scale falsification and projected-memory target

## Status

**Current result: a major candidate explanation is falsified.**

The active question is whether slow fatigue-scale memory under normal cyclic loading can emerge from a fixed generalized Lennard-Jones normal lattice without inserting an empirical relaxation time.

The present calculation shows that ordinary conservative acoustic modes and ordinary local LJ critical softening are far too fast to directly explain a 20 Hz fatigue-scale memory in the current baseline.

## 1. Atomic time scale

Using the existing calibration values

$$
a_0=2.8627442948\times10^{-10}\ {\rm m},
$$

$$
E=69\ {\rm GPa},
$$

$$
A_0=6.0338\times10^{-20}\ {\rm m^2},
$$

and the Al atomic mass, define

$$
\boxed{
t_0=\sqrt{\frac{M_{\rm Al}a_0}{EA_0}}.
}
$$

Numerically,

$$
\boxed{t_0=5.55046\times10^{-14}\ {\rm s}.}
$$

The associated lattice-speed scale is

$$
\boxed{c_0=\frac{a_0}{t_0}\approx5.16\times10^3\ {\rm m/s}.}
$$

## 2. Slowest fixed-free normal-chain mode

Linearize the generalized-LJ chain about equilibrium. Because the normalized potential satisfies

$$
\phi''(1)=1,
$$

the linearized nearest-neighbor stiffness is unity.

For $L$ moving atoms with the left end fixed and the right end free,

$$
q_j=\frac{(2j-1)\pi}{2L+1},
$$

and

$$
\boxed{
\omega_j^*=2\sin\left(\frac{q_j}{2}\right).
}
$$

Therefore

$$
\boxed{
\omega_{\min}^*=2\sin\left[\frac{\pi}{2(2L+1)}\right].
}
$$

This relation is **EXACT for the linearized fixed-free chain**.

For 20 Hz,

$$
\omega_{20}^*=2\pi(20)t_0=6.97492\times10^{-12}.
$$

Inverting the exact mode relation gives

$$
\boxed{L\approx2.2521\times10^{11}},
$$

corresponding to

$$
\boxed{\ell\approx64.5\ {\rm m}.}
$$

A portable or ordinary laboratory specimen therefore cannot obtain a 20 Hz microscopic memory merely from the slowest conservative acoustic mode of this baseline chain.

## 3. Can LJ critical softening create a 20 Hz local mode?

The local small-oscillation frequency is

$$
f_{\rm loc}(\lambda)=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}.
$$

At the normal LJ stability limit,

$$
\phi''(\lambda_c)=0,
$$

with

$$
\lambda_c=1.1077715386.
$$

Near $\lambda_c$,

$$
\phi''(\lambda)\approx|\phi'''(\lambda_c)|(\lambda_c-\lambda),
$$

and

$$
\phi'''(\lambda_c)\approx-2.78642.
$$

A 20 Hz local mode would require

$$
\boxed{
\phi''\sim(2\pi\,20\,t_0)^2=4.86495\times10^{-23}.
}
$$

Therefore

$$
\boxed{
\lambda_c-\lambda\sim1.75\times10^{-23}.
}
$$

This is not a credible ordinary-stress slow mechanism.

## 4. 100 MPa comparison

On the stable homogeneous branch,

$$
\phi'(\lambda_{100})=\frac{100\ {\rm MPa}}{69\ {\rm GPa}},
$$

which gives

$$
\boxed{\lambda_{100}=1.00147203.}
$$

The tangent stiffness is still

$$
\boxed{\phi''(\lambda_{100})=0.969214,}
$$

and the corresponding local harmonic frequency is

$$
\boxed{f_{\rm loc}\approx2.82\times10^{12}\ {\rm Hz}.}
$$

Thus 100 MPa is nowhere near the critical-softening regime required to create a 20 Hz local time scale.

## 5. Falsified explanation

The following explanation is rejected for the current normal-LJ baseline:

> ordinary conservative normal lattice modes or ordinary local LJ critical softening directly provide the slow 20 Hz fatigue memory.

This is useful because the roughly $10^{11}$–$10^{12}$ time-scale separation must not be hidden inside an arbitrary fitted relaxation parameter.

## 6. Next target

A slow reduced state may still emerge after projection if the reduced description eliminates genuinely slow normal structural variables. Such a state must be derived from the normal atomic mechanics.

The next candidates to test are:

1. pair-spacing correlations $P_2(a_1,a_2,t)$ and higher hierarchy terms;
2. localized normal energy density;
3. the spatial envelope of normal strain or energy localization;
4. free-surface normal-opening coordinates;
5. the distribution of local tangent stiffness $\phi''(a_i/a_0)$;
6. first-passage variables associated with normal stability loss.

The acceptance criterion is

$$
\boxed{
\text{a slow state is accepted only if its time scale follows from mechanics.}
}
$$

---

# 한국어 번역 — Normal-LJ 시간척도 반증과 projected-memory 목표

## 상태

**현재 결과는 중요한 후보 설명 하나를 반증했다.**

활성 연구문제는 경험적인 relaxation time을 넣지 않고 고정 generalized Lennard-Jones 수직 lattice에서 느린 피로 시간척도의 memory가 나올 수 있는가이다.

이번 계산에서는 일반적인 conservative acoustic mode와 일반적인 LJ local critical softening이 현재 baseline에서 20 Hz 피로 memory를 직접 설명하기에는 지나치게 빠르다는 것이 확인됐다.

## 1. 원자 시간척도

기존 calibration 값

$$
a_0=2.8627442948\times10^{-10}\ {\rm m},
$$

$$
E=69\ {\rm GPa},
$$

$$
A_0=6.0338\times10^{-20}\ {\rm m^2}
$$

와 Al 원자질량을 사용하여

$$
\boxed{
t_0=\sqrt{\frac{M_{\rm Al}a_0}{EA_0}}
}
$$

를 정의한다.

수치적으로

$$
\boxed{t_0=5.55046\times10^{-14}\ {\rm s}}
$$

이다.

이에 대응하는 lattice-speed scale은

$$
\boxed{c_0=\frac{a_0}{t_0}\approx5.16\times10^3\ {\rm m/s}}
$$

이다.

## 2. 가장 느린 fixed-free normal-chain mode

평형점 주변에서 generalized-LJ chain을 선형화한다. 현재 normalization에서는

$$
\phi''(1)=1
$$

이므로 nearest-neighbor stiffness가 1인 harmonic chain이 된다.

왼쪽이 고정되고 오른쪽이 자유로운 $L$개의 moving atom chain에서

$$
q_j=\frac{(2j-1)\pi}{2L+1}
$$

이고,

$$
\boxed{
\omega_j^*=2\sin\left(\frac{q_j}{2}\right)
}
$$

이다.

따라서 가장 느린 mode는

$$
\boxed{
\omega_{\min}^*=2\sin\left[\frac{\pi}{2(2L+1)}\right]
}
$$

이다.

이 관계는 **선형화된 fixed-free chain에서는 정확하다.**

20 Hz의 dimensionless angular frequency는

$$
\omega_{20}^*=2\pi(20)t_0=6.97492\times10^{-12}
$$

이다.

정확한 mode 식을 역으로 풀면

$$
\boxed{L\approx2.2521\times10^{11}}
$$

개의 moving atom이 필요하고, 이는

$$
\boxed{\ell\approx64.5\ {\rm m}}
$$

길이에 해당한다.

따라서 portable 또는 일반적인 laboratory specimen에서는 현재 baseline chain의 가장 느린 conservative acoustic mode만으로 20 Hz microscopic memory를 만들 수 없다.

## 3. LJ critical softening으로 20 Hz local mode를 만들 수 있는가?

stretch $\lambda$에서 local small-oscillation frequency는

$$
f_{\rm loc}(\lambda)=\frac{\sqrt{\phi''(\lambda)}}{2\pi t_0}
$$

이다.

LJ normal stability limit에서는

$$
\phi''(\lambda_c)=0
$$

이고

$$
\lambda_c=1.1077715386
$$

이다.

$\lambda_c$ 근처에서

$$
\phi''(\lambda)\approx|\phi'''(\lambda_c)|(\lambda_c-\lambda)
$$

이며

$$
\phi'''(\lambda_c)\approx-2.78642
$$

이다.

20 Hz local mode에 필요한 dimensionless tangent stiffness는

$$
\boxed{
\phi''\sim4.86495\times10^{-23}
}
$$

이다.

따라서

$$
\boxed{
\lambda_c-\lambda\sim1.75\times10^{-23}
}
$$

수준까지 inflection point에 가까워야 한다.

이는 ordinary-stress slow mechanism으로 보기 어렵다.

## 4. 100 MPa와 비교

stable homogeneous branch에서

$$
\phi'(\lambda_{100})=\frac{100\ {\rm MPa}}{69\ {\rm GPa}}
$$

를 풀면

$$
\boxed{\lambda_{100}=1.00147203}
$$

이다.

이때 tangent stiffness는

$$
\boxed{\phi''(\lambda_{100})=0.969214}
$$

이고 local harmonic frequency는

$$
\boxed{f_{\rm loc}\approx2.82\times10^{12}\ {\rm Hz}}
$$

이다.

즉 100 MPa 상태는 20 Hz 시간척도를 만드는 critical-softening 상태와 전혀 가깝지 않다.

## 5. 현재 반증된 설명

현재 normal-LJ baseline에서 다음 설명은 기각한다.

> 일반적인 conservative normal lattice mode 또는 일반적인 LJ local critical softening 자체가 느린 20 Hz 피로 memory를 직접 만든다.

이 결과가 중요한 이유는 약 $10^{11}$–$10^{12}$의 시간척도 차이를 임의의 relaxation parameter 안에 숨기면 안 되기 때문이다.

## 6. 다음 목표

축약상태가 실제로 천천히 진화하는 normal microscopic variable을 제거하고 있다면 projection 이후 slow state가 나타날 가능성은 남아 있다. 하지만 그 state는 반드시 normal atomic mechanics로부터 유도되어야 한다.

다음 후보들을 검증한다.

1. pair-spacing correlation $P_2(a_1,a_2,t)$ 및 higher hierarchy;
2. localized normal energy density;
3. normal strain/energy localization의 spatial envelope;
4. free-surface normal-opening coordinate;
5. local tangent stiffness $\phi''(a_i/a_0)$의 분포;
6. normal stability loss와 연결된 first-passage variable.

기준은 엄격하다.

$$
\boxed{
\text{slow state의 시간척도는 반드시 mechanics에서 나와야 한다.}
}
$$
