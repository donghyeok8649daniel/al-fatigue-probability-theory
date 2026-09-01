# Milestone 13 — One-Dimensional Statistical Cell from Correlation

## Scope

The active theory remains strictly one-dimensional. The goal of this milestone is to define, from the 1D layer-spacing statistics themselves, how many represented spacings should be treated as effectively independent, partially dependent, or completely identical when probabilities are later aggregated.

No FEM mesh, transverse correlation model, three-dimensional volume model, arbitrary block size, or empirical independence threshold is introduced.

## 1. Three logically different dependence classes

Let $X_i=\lambda_i$ and $X_j=\lambda_j$ denote two layer-spacing random variables.

### 1.1 Complete identical dependence

The exact condition that the two variables are the same almost surely is

$$
\boxed{
\mathbb E[(X_i-X_j)^2]=0
\iff
X_i=X_j\quad\text{almost surely}.
}
$$

This is stronger and cleaner than declaring two variables identical because a measured correlation coefficient is merely close to one.

For a stationary pair with equal finite variance $\sigma_\lambda^2$,

$$
\mathbb E[(X_i-X_j)^2]
=2\sigma_\lambda^2(1-\rho_{ij}),
$$

so $\rho_{ij}=1$ is equivalent to complete identical dependence in that restricted equal-variance setting.

Classification: **EXACT / IDENTITY**.

### 1.2 Full statistical independence

Independence requires factorization of the joint density,

$$
\boxed{
P_2(\lambda,\lambda';k)
=P_1(\lambda)P_1(\lambda').
}
$$

Zero covariance or $\rho_k=0$ is not sufficient in a non-Gaussian theory.

Classification: **DEFINITION**.

### 1.3 Partial dependence

Everything between the two limits above is partially dependent. A single correlation coefficient is useful for second-moment aggregation but does not contain the full joint dependence needed for arbitrary crack-initiation events.

## 2. Exact finite-$M$ variance identity

Assume a second-order stationary 1D spacing sequence with

$$
\operatorname{Var}(\lambda_i)=\sigma_\lambda^2,
$$

and lag correlation

$$
\rho_k
=\frac{\operatorname{Cov}(\lambda_i,\lambda_{i+k})}{\sigma_\lambda^2}.
$$

For

$$
\bar\lambda_M=\frac1M\sum_{i=1}^M\lambda_i,
$$

direct expansion of the double covariance sum gives

$$
\boxed{
\operatorname{Var}(\bar\lambda_M)
=
\frac{\sigma_\lambda^2}{M}
\left[
1+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)\rho_k
\right].
}
$$

Define the finite correlation factor

$$
\boxed{
\tau_M
=
1+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)\rho_k.
}
$$

Then

$$
\operatorname{Var}(\bar\lambda_M)
=
\frac{\sigma_\lambda^2}{M/\tau_M}.
$$

Classification: **EXACT / IDENTITY** when the true stationary correlations are used. If empirical $\rho_k$ from one finite snapshot are inserted, the result is an **ESTIMATOR / DIAGNOSTIC**.

## 3. Variance-equivalent independent count

This motivates the definition

$$
\boxed{
M_{\rm eff}
=\frac{M}{\tau_M}.
}
$$

$M_{\rm eff}$ is the number of independent equal-variance variables that would give the same variance of the sample mean.

It is therefore a **variance-equivalent independent count**, not a proof that the original variables are jointly independent.

Two limiting cases are exact.

### Independent limit

If

$$
\rho_k=0\qquad(k>0),
$$

then

$$
\tau_M=1,
\qquad
M_{\rm eff}=M.
$$

### Completely identical limit

If every represented spacing is the same random variable,

$$
\rho_k=1\qquad(0\le k<M),
$$

then

$$
\tau_M=M,
\qquad
M_{\rm eff}=1.
$$

Thus the definition exactly reproduces both extreme cases that the later probability aggregation must distinguish.

## 4. One-dimensional statistical characteristic length

Let $a_0$ be the equilibrium represented-layer spacing. The total represented axial length is approximately

$$
L_M=Ma_0.
$$

Define

$$
\boxed{
\ell_{\rm stat}^{(2)}(M)
=\frac{L_M}{M_{\rm eff}}
=a_0\tau_M.
}
$$

The superscript $(2)$ emphasizes that this is a **second-moment / variance-equivalent characteristic length**.

Its exact limits are

$$
\ell_{\rm stat}^{(2)}=a_0
$$

for independent spacings and

$$
\ell_{\rm stat}^{(2)}=Ma_0
$$

when all $M$ spacings are completely identical.

This is the current preferred 1D statistical-cell length because it is derived from the correlation structure instead of being chosen as a numerical mesh size.

## 5. Why this length does not prove event independence

Suppose one statistical cell has crack-initiation probability $q$.

For $N$ genuinely independent cells,

$$
\boxed{
P(\text{at least one initiation})
=1-(1-q)^N.
}
$$

But if $N$ nominal cells are completely identical copies of one random variable, then

$$
\boxed{
P(\text{at least one initiation})=q,
}
$$

not $1-(1-q)^N$.

For partial dependence, neither formula is generally exact. The required object is a joint survival/event distribution, ultimately connected to $P_2$, $P_3$, and higher dependence information.

Therefore $M_{\rm eff}$ is appropriate for variance-based coarse graining, but arbitrary event probabilities must not be aggregated as independent unless joint factorization has been justified.

## 6. Characteristic area and volume are deliberately not identified yet

Three different quantities must not be conflated:

1. $A_0$ — the representative layer-patch area entering the mechanical energy calibration

$$
E_0=EA_0a_0;
$$

2. $A_{\rm stat}$ — a transverse statistical correlation/independence area;
3. $V_{\rm stat}=A_{\rm stat}\ell_{\rm stat}$ — a future statistical correlation volume.

The present 1D theory can derive or test only the axial correlation length. It contains no transverse coordinates and therefore cannot derive $A_{\rm stat}$.

Hence

$$
\boxed{A_0\neq A_{\rm stat}\quad\text{unless an independent physical argument proves the identification}.}
$$

Likewise a three-dimensional characteristic volume is outside the present scope.

This distinction is essential because $A_0$ controls the physical energy scale and temperature mapping, while $A_{\rm stat}$ would control how many statistically independent transverse regions exist. They solve different physical problems.

## 7. Connection to the measured 1D spatial correlations

The existing deterministic chain already provides empirical

$$
C_k(t),\qquad \rho_k(t).
$$

The next calculation should therefore evaluate

$$
\tau_M(t),
\qquad
M_{\rm eff}(t),
\qquad
\ell_{\rm stat}^{(2)}(t)
$$

for the dynamically matched system-size sweep.

A material-like statistical length requires convergence with increasing represented system size. If

$$
\ell_{\rm stat}^{(2)}(M)
$$

keeps scaling with $M$, then the current 1D protocol still contains system-scale coherent motion and has not produced a local material correlation length.

That outcome would be a physical diagnostic, not a reason to impose an arbitrary cutoff.

## 8. Current conclusion

The active 1D probability theory now distinguishes:

$$
\boxed{
\text{complete identity}
\;\neq\;
\text{partial dependence}
\;\neq\;
\text{independence}.
}
$$

For second-moment aggregation, the physically derived finite-$M$ statistical-cell measure is

$$
\boxed{
M_{\rm eff}=\frac{M}{\tau_M},
\qquad
\ell_{\rm stat}^{(2)}=a_0\tau_M.
}
$$

Full event-probability aggregation remains a higher-order joint-distribution problem and must not be replaced by an independence product without justification.

---

# 한국어 번역 — 상관으로 정의하는 1차원 통계 셀

## 범위

활성 이론은 계속 엄격한 1차원으로 유지한다. 이번 milestone의 목표는 나중에 확률을 합칠 때 represented spacing들을 몇 개의 독립 확률변수로 볼지, 부분 종속으로 볼지, 또는 완전히 동일한 하나의 확률변수로 볼지를 1D layer-spacing 통계 자체에서 정의하는 것이다.

FEM mesh, 횡방향 상관모델, 3차원 특성부피, 임의 block size, 경험적 independence threshold는 도입하지 않는다.

## 1. 서로 다른 세 종류의 종속성

$X_i=\lambda_i$, $X_j=\lambda_j$를 두 layer-spacing 확률변수라고 하자.

### 1.1 완전히 동일한 종속

두 변수가 거의 확실하게 같은 변수라는 정확한 조건은

$$
\boxed{
\mathbb E[(X_i-X_j)^2]=0
\iff
X_i=X_j\quad\text{almost surely}
}
$$

이다.

측정된 상관계수가 단지 1에 가깝다는 이유로 두 변수를 동일하다고 선언하는 것보다 강하고 명확한 기준이다.

stationary pair이고 두 변수의 유한 분산이 같은 경우에는

$$
\mathbb E[(X_i-X_j)^2]
=2\sigma_\lambda^2(1-\rho_{ij})
$$

이므로 이 제한된 조건에서는 $\rho_{ij}=1$과 완전 동일 종속이 동치다.

분류: **EXACT / IDENTITY**.

### 1.2 완전한 통계적 독립

독립은 joint density의 factorization을 요구한다.

$$
\boxed{
P_2(\lambda,\lambda';k)
=P_1(\lambda)P_1(\lambda')
}
$$

비가우시안 이론에서는 covariance가 0이거나 $\rho_k=0$인 것만으로 독립을 보장할 수 없다.

분류: **DEFINITION**.

### 1.3 부분 종속

위 두 극한 사이의 상태가 부분 종속이다. 하나의 correlation coefficient는 second-moment를 합칠 때는 유용하지만 임의의 crack-initiation event에 필요한 full joint dependence를 담지는 못한다.

## 2. 정확한 finite-$M$ 평균분산 식

second-order stationary 1D spacing sequence에 대해

$$
\operatorname{Var}(\lambda_i)=\sigma_\lambda^2
$$

이고 lag correlation을

$$
\rho_k
=\frac{\operatorname{Cov}(\lambda_i,\lambda_{i+k})}{\sigma_\lambda^2}
$$

라고 하자.

$$
\bar\lambda_M=\frac1M\sum_{i=1}^M\lambda_i
$$

에 대해 covariance double sum을 직접 전개하면

$$
\boxed{
\operatorname{Var}(\bar\lambda_M)
=
\frac{\sigma_\lambda^2}{M}
\left[
1+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)\rho_k
\right]
}
$$

를 얻는다.

finite correlation factor를

$$
\boxed{
\tau_M
=
1+2\sum_{k=1}^{M-1}
\left(1-\frac{k}{M}\right)\rho_k
}
$$

라고 정의하면

$$
\operatorname{Var}(\bar\lambda_M)
=
\frac{\sigma_\lambda^2}{M/\tau_M}
$$

이다.

true stationary correlation을 사용할 때 **EXACT / IDENTITY**다. finite snapshot의 empirical $\rho_k$를 넣으면 **ESTIMATOR / DIAGNOSTIC**이다.

## 3. 분산 기준 유효 독립개수

따라서

$$
\boxed{
M_{\rm eff}=\frac{M}{\tau_M}
}
$$

를 정의한다.

$M_{\rm eff}$는 같은 평균분산을 만드는 독립 equal-variance 확률변수의 개수다. 따라서 **variance-equivalent independent count**이며, 원래 변수들의 joint distribution이 실제로 factorize된다는 증명은 아니다.

두 극한에서는 정확하다.

### 독립 극한

$$
\rho_k=0\qquad(k>0)
$$

이면

$$
\tau_M=1,
\qquad
M_{\rm eff}=M.
$$

### 완전 동일 종속 극한

모든 represented spacing이 하나의 동일한 확률변수라면

$$
\rho_k=1\qquad(0\le k<M)
$$

이고

$$
\tau_M=M,
\qquad
M_{\rm eff}=1.
$$

따라서 나중의 probability aggregation에서 반드시 구분해야 하는 두 극단을 정확히 재현한다.

## 4. 1차원 통계 특성길이

$a_0$를 equilibrium represented-layer spacing이라고 하자. 전체 represented axial length는 대략

$$
L_M=Ma_0
$$

이다.

다음을 정의한다.

$$
\boxed{
\ell_{\rm stat}^{(2)}(M)
=\frac{L_M}{M_{\rm eff}}
=a_0\tau_M
}
$$

위첨자 $(2)$는 이것이 **second-moment / variance-equivalent characteristic length**라는 점을 명시한다.

독립 spacing이면 정확히

$$
\ell_{\rm stat}^{(2)}=a_0
$$

이고, $M$개 spacing이 전부 완전히 동일하면

$$
\ell_{\rm stat}^{(2)}=Ma_0
$$

이다.

따라서 현재 1D 통계셀 길이는 numerical mesh size를 임의로 고르는 대신 correlation structure에서 직접 유도하는 이 정의를 우선 사용한다.

## 5. 이 길이가 event independence를 증명하지는 않는다

statistical cell 하나의 crack-initiation probability가 $q$라고 하자.

$N$개의 셀이 정말 독립이면

$$
\boxed{
P(\text{적어도 하나 개시})
=1-(1-q)^N
}
$$

이다.

하지만 nominal cell $N$개가 사실 하나의 확률변수를 완전히 동일하게 복제한 것이라면

$$
\boxed{
P(\text{적어도 하나 개시})=q
}
$$

이지 $1-(1-q)^N$이 아니다.

부분 종속이면 두 식 모두 일반적으로 정확하지 않다. 이때는 joint survival/event distribution이 필요하며 결국 $P_2$, $P_3$ 이상의 dependence 정보와 연결된다.

따라서 $M_{\rm eff}$는 variance-based coarse graining에는 적합하지만 임의 event probability를 독립식으로 합치려면 joint factorization을 별도로 정당화해야 한다.

## 6. 특성면적과 특성부피는 아직 동일시하지 않는다

다음 세 양은 서로 섞으면 안 된다.

1. $A_0$ — mechanical energy calibration에 들어가는 representative layer-patch area

$$
E_0=EA_0a_0;
$$

2. $A_{\rm stat}$ — 횡방향 statistical correlation/independence area;
3. $V_{\rm stat}=A_{\rm stat}\ell_{\rm stat}$ — 나중의 statistical correlation volume.

현재 1D 이론은 axial correlation length만 유도하거나 검증할 수 있다. 횡방향 coordinate가 없으므로 $A_{\rm stat}$를 유도할 수 없다.

따라서

$$
\boxed{
A_0\neq A_{\rm stat}
\quad\text{unless an independent physical argument proves the identification}
}
$$

으로 둔다.

3차원 특성부피도 현재 범위 밖이다.

이 구분이 중요한 이유는 $A_0$는 physical energy scale과 temperature mapping을 결정하지만, $A_{\rm stat}$는 횡방향으로 몇 개의 독립 확률영역이 존재하는지를 결정하기 때문이다. 서로 다른 물리문제다.

## 7. 기존 1D spatial correlation 결과와의 연결

현재 deterministic chain에서 이미 empirical

$$
C_k(t),\qquad\rho_k(t)
$$

를 계산하고 있다.

따라서 다음 계산은 dynamically matched system-size sweep에서

$$
\tau_M(t),
\qquad
M_{\rm eff}(t),
\qquad
\ell_{\rm stat}^{(2)}(t)
$$

를 계산하는 것이다.

material-like statistical length라고 부르려면 represented system size를 늘릴 때 수렴해야 한다.

만약

$$
\ell_{\rm stat}^{(2)}(M)
$$

가 계속 $M$과 같이 커진다면 현재 1D protocol에는 여전히 system-scale coherent motion이 남아 있고 local material correlation length가 아직 나타나지 않은 것이다.

이 경우 임의 cutoff를 넣는 것이 아니라 그 자체를 물리적 diagnostic으로 받아들인다.

## 8. 현재 결론

활성 1D probability theory는 이제

$$
\boxed{
\text{완전 동일}
\;\neq\;
\text{부분 종속}
\;\neq\;
\text{독립}
}
$$

을 명시적으로 구분한다.

second-moment aggregation에 대한 finite-$M$ statistical cell 척도는

$$
\boxed{
M_{\rm eff}=\frac{M}{\tau_M},
\qquad
\ell_{\rm stat}^{(2)}=a_0\tau_M
}
$$

이다.

반면 full event-probability aggregation은 higher-order joint-distribution 문제이며 정당화 없이 independence product로 바꾸지 않는다.
