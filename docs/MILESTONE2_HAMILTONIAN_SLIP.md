# Milestone 2 — Hamiltonian non-affine slip as a proof of secular structural evolution

## Status

**Proof of principle only. This is not yet a calibrated aluminum fatigue model.**

The purpose of this step is to test whether cycle-to-cycle structural evolution can arise from microscopic conservative dynamics without inserting a phenomenological damping law, damage variable, or empirical fatigue evolution equation.

---

## 1. Microscopic model

Introduce one resolved non-affine slip coordinate $s$ coupled to a long harmonic lattice bath with coordinates $u_j$. The total driven Hamiltonian is

$$
H(t)=\frac{P_s^2}{2M}+V_\gamma(s)
+\frac{k_c}{2}(s-u_1)^2
+\sum_{j=1}^{\infty}\frac{p_j^2}{2m}
+\frac{k}{2}\sum_{j=1}^{\infty}(u_{j+1}-u_j)^2
-F(t)s.
$$

The equations of motion are

$$
M\ddot s=-V_\gamma'(s)-k_c(s-u_1)+F(t),
$$

$$
m\ddot u_1=k_c(s-u_1)+k(u_2-u_1),
$$

and for interior bath sites

$$
m\ddot u_j=k(u_{j+1}-2u_j+u_{j-1}).
$$

There is **no viscous damping term** anywhere in the full system.

### Classification

- **EXACT under the stated model:** Newton/Hamilton equations above.
- **CONTROLLED APPROXIMATION:** the finite numerical chain used to approximate a semi-infinite bath before reflected waves return.
- **CONTROLLED APPROXIMATION:** the periodic slip landscape is represented by one Fourier harmonic,

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}\left[1-\cos\left(\frac{2\pi s}{b}\right)\right].
$$

For real fcc Al this function must eventually be replaced by a DFT/EAM generalized-stacking-fault energy surface $\gamma(\mathbf s)$ rather than treated as exact.

---

## 2. Why this coordinate is physically motivated

For fcc Al, relative displacement of two halves of a crystal across a $\{111\}$ plane defines the generalized stacking-fault energy surface. First-principles calculations show that this energy landscape controls dislocation-core structure and slip energetics. Therefore $s$ is not introduced as an empirical plastic-strain variable; it is a coarse-grained relative atomic displacement.

Relevant primary literature:

- G. Lu, N. Kioussis, V. V. Bulatov, E. Kaxiras, *Generalized-stacking-fault energy surface and dislocation properties of aluminum*, Phys. Rev. B **62**, 3099 (2000), DOI: 10.1103/PhysRevB.62.3099.
- C. Brandl, P. M. Derlet, H. Van Swygenhoven, *General-stacking-fault energies in highly strained metallic environments: Ab initio calculations*, Phys. Rev. B **76**, 054124 (2007), DOI: 10.1103/PhysRevB.76.054124.
- R. J. Rubin, *Momentum Autocorrelation Functions and Energy Transport in Harmonic Crystals Containing Isotopic Defects*, Phys. Rev. **131**, 964 (1963), DOI: 10.1103/PhysRev.131.964.
- R. Zwanzig, *Nonlinear generalized Langevin equations*, J. Stat. Phys. **9**, 215–220 (1973), DOI: 10.1007/BF01008729.

The Rubin/Zwanzig connection is important: eliminating the harmonic bath produces memory and radiation terms in the resolved coordinate even though the full system is conservative.

---

## 3. Exact energy balance of the driven conservative model

Define internal energy without the external potential $-F(t)s$:

$$
E_{\rm int}
=\frac{P_s^2}{2M}+V_\gamma(s)
+\frac{k_c}{2}(s-u_1)^2
+\sum_j\frac{p_j^2}{2m}
+\frac{k}{2}\sum_j(u_{j+1}-u_j)^2.
$$

Using the equations of motion,

$$
\boxed{\frac{dE_{\rm int}}{dt}=F(t)\dot s(t)}.
$$

Therefore the hysteresis work over one cycle is

$$
\boxed{A_H=\oint F\,ds}
$$

and this work is transferred into the unresolved lattice modes and/or retained in a changed structural state. It is not numerical or phenomenological dissipation.

---

## 4. Cycle map and structural accumulation

At cycle endpoints $t_N=NT$, define

$$
s_N=s(NT).
$$

If the trajectory is bounded inside a single basin, a periodic steady state gives

$$
s_{N+1}=s_N.
$$

If inter-basin transitions occur, the cycle map can instead satisfy

$$
\boxed{s_{N+1}\neq s_N}.
$$

For an ensemble of representative areas,

$$
P_s(s,t)=\lim_{N_{\rm RA}\to\infty}
\frac{1}{N_{\rm RA}}
\sum_{\alpha=1}^{N_{\rm RA}}
\delta\!\left(s-s_\alpha(t)\right),
$$

and the exact kinematic conservation law is

$$
\partial_tP_s+\partial_s(P_s v_s)=0.
$$

The natural extension of the original spacing theory is then

$$
\boxed{P(a,s,t)},
$$

with marginal

$$
P(a,t)=\int P(a,s,t)\,ds,
$$

and joint continuity equation

$$
\boxed{
\partial_tP+\partial_a(Pv_a)+\partial_s(Pv_s)=0.
}
$$

This preserves the original $P(a,t)$ framework while making the minimal non-affine state explicit.

---

## 5. Numerical reference experiment

The repository simulation uses nondimensional parameters

$$
M=m=k=k_c=b=1,\qquad \Delta_\gamma=0.1,
$$

with

$$
F(t)=F_a\sin(0.2t)
$$

after a smooth two-cycle ramp. A long finite chain is used only long enough that reflected waves cannot return to the resolved coordinate.

Three regimes were found:

| $F_a$ | Long-time cycle-end behavior | Interpretation |
|---:|---|---|
| 0.34 | $s_N\approx-0.024$ | bounded intrawell periodic state |
| 0.40 | $s_N\approx-1.965$ after an initial transition | finite transient structural relocation, then periodic |
| 0.50 | $s_N$ decreases by approximately one period per cycle | running/inter-basin secular state |

For $F_a=0.50$, representative cycle endpoints were

$$
-5.8529,\,-6.8542,\,-7.8523,\,-8.8538,\,-9.8519,\,-10.8534,
$$

so asymptotically

$$
\boxed{s_{N+1}-s_N\approx-1.00.}
$$

The relative global energy-balance error in the direct Newton integration was

$$
\boxed{1.8\times10^{-7}},
$$

which is far below the structural drift and rules out numerical energy loss as its origin in this reference run.

---

## 6. What has and has not been proved

### Established by this model

1. A conservative microscopic system can give a nonzero resolved hysteresis loop.
2. Adding a physically motivated periodic non-affine energy landscape permits deterministic inter-basin transitions.
3. A cycle map with $s_{N+1}\neq s_N$ can arise without a fitted fatigue evolution law.
4. The full energy balance remains conservative; energy is redistributed into propagating lattice modes and structural potential energy.

### Not established

1. The nondimensional running state at $F_a=0.50$ is **not** a prediction of fatigue life in Al.
2. The sinusoidal $V_\gamma$ is not the true Al $\gamma$-surface.
3. The running state should not yet be called a dislocation multiplication law or crack-initiation law.
4. The low-stress 20 Hz regime of real high-purity Al has not been reproduced.
5. A single homogeneous perfect-slip coordinate will generally have an ideal barrier far above ordinary macroscopic fatigue stresses; defects, surfaces, stress concentrations, thermal initial conditions, and multi-slip correlations must be derived rather than hidden in fitted parameters.

A particularly important numerical observation is that the transition from bounded to running motion is nonlinear and non-monotonic in forcing amplitude. Therefore fitting a single scalar threshold to these dynamics would throw away essential phase-space information.

---

# 한국어 번역

## 상태

**현재 단계는 원리 증명이다. 아직 보정된 알루미늄 피로수명 모델이 아니다.**

목적은 경험적인 감쇠계수, 손상변수, 피로 누적식을 넣지 않고도 보존적인 미시 역학에서 cycle마다 구조상태가 달라질 수 있는지 확인하는 것이다.

## 1. 미시 역학 모델

원자면 사이의 상대적인 비아핀 slip 좌표 $s$를 하나 두고, 이를 긴 조화격자 bath $u_j$와 결합한다.

$$
H(t)=\frac{P_s^2}{2M}+V_\gamma(s)
+\frac{k_c}{2}(s-u_1)^2
+\sum_{j=1}^{\infty}\frac{p_j^2}{2m}
+\frac{k}{2}\sum_{j=1}^{\infty}(u_{j+1}-u_j)^2
-F(t)s.
$$

전체계 어디에도 점성 감쇠항을 넣지 않는다. $s$의 주기 퍼텐셜은 현재

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}
\left[1-\cos\left(\frac{2\pi s}{b}\right)\right]
$$

로 두었지만, 이것은 **근사**다. 실제 Al 계산에서는 DFT/EAM으로 얻은 $\gamma(\mathbf s)$를 넣어야 한다.

## 2. 왜 $s$가 임의의 손상변수가 아닌가

FCC Al의 $\{111\}$ 면을 기준으로 결정의 위쪽과 아래쪽을 상대이동시키면 generalized stacking-fault energy surface가 정의된다. 따라서 $s$는 경험적 plastic strain을 발명한 것이 아니라 원자좌표에서 정의 가능한 상대변위다.

## 3. 에너지 수지

외력 퍼텐셜을 제외한 내부에너지는 정확히

$$
\boxed{\frac{dE_{\rm int}}{dt}=F(t)\dot s(t)}
$$

를 만족한다. 따라서

$$
\boxed{A_H=\oint F\,ds}
$$

는 사라진 에너지가 아니라 다른 격자모드로 전달되거나 변경된 구조상태에 남은 에너지다.

## 4. cycle 누적

cycle 끝의 상태를 $s_N=s(NT)$라 하면 한 well 안에서 주기상태가 만들어진 경우

$$
s_{N+1}=s_N
$$

이지만, basin을 넘는 전이가 발생하면

$$
\boxed{s_{N+1}\neq s_N}
$$

이 가능하다.

대표영역들의 집합에 대해 $P_s(s,t)$를 정의하면 원래 이론은 자연스럽게

$$
\boxed{P(a,s,t)}
$$

로 확장된다. 기존 $P(a,t)$는

$$
P(a,t)=\int P(a,s,t)\,ds
$$

라는 정확한 marginal로 그대로 남는다.

## 5. 수치결과

무차원 기준계에서 세 가지 응답이 확인되었다.

- $F_a=0.34$: 동일 well 내부의 주기응답. 히스테리시스는 있지만 secular accumulation은 없음.
- $F_a=0.40$: 초기에 약 두 개의 well을 이동한 뒤 새 주기상태에서 고정.
- $F_a=0.50$: 이후 거의 cycle당 한 period씩 계속 이동.

$F_a=0.50$에서 후반 cycle의 $s_N$은

$$
-5.8529,\,-6.8542,\,-7.8523,\,-8.8538,\,-9.8519,\,-10.8534
$$

였고,

$$
\boxed{s_{N+1}-s_N\approx-1}
$$

이었다. 전체 Newton 적분의 에너지 수지 상대오차는 약

$$
\boxed{1.8\times10^{-7}}
$$

이었다.

## 6. 현재 결론

이번 계산으로 **미시 보존역학 → 히스테리시스 → basin transition → cycle-to-cycle 구조상태 변화**라는 연결 자체는 가능하다는 것을 보였다.

하지만 이것을 곧바로 실제 Al의 저응력 피로라고 부르면 안 된다. 다음 단계는 실제 Al의 $\gamma$-surface, 표면/결함에 의한 국부응력 집중, 온도에 따른 미시 초기상태 분포, 여러 slip system과 $a$의 상관관계를 순서대로 넣고 어느 단계가 실제 저응력 secular evolution을 만드는지 반증 가능한 방식으로 확인하는 것이다.
