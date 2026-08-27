# Aluminum gamma-surface constraints for the next calibration step

## Purpose

The nonlinear Hamiltonian slip-bath model proves that a non-affine periodic energy landscape can generate inter-basin cycle-state changes without an empirical fatigue law. The next step is to replace the artificial nondimensional barrier with an Al-specific atomistic energy landscape.

## Literature constraints

First-principles and atomistic studies of fcc Al generalized stacking-fault energetics consistently place the unstable stacking-fault energy in the rough range of order $10^{-1}\,\mathrm{J/m^2}$, with reported values depending on relaxation and electronic-structure method. Representative literature values include approximately $175$–$224\,\mathrm{mJ/m^2}$ for unstable stacking-fault energies and roughly $120$–$166\,\mathrm{mJ/m^2}$ for intrinsic stacking-fault energies in common first-principles datasets.

Primary sources relevant to this project:

- G. Lu, N. Kioussis, V. V. Bulatov, E. Kaxiras, Phys. Rev. B **62**, 3099 (2000), DOI: 10.1103/PhysRevB.62.3099.
- C. Brandl, P. M. Derlet, H. Van Swygenhoven, Phys. Rev. B **76**, 054124 (2007), DOI: 10.1103/PhysRevB.76.054124.
- S. Ogata, J. Li, S. Yip, *Ideal Pure Shear Strength of Aluminum and Copper*, Science **298**, 807–811 (2002), DOI: 10.1126/science.1076652.

The numerical spread is not a nuisance to be fitted away. It reflects method, relaxation, strain state, and path dependence, so the future model should carry the chosen atomistic input explicitly.

## Consequence for ordinary fatigue stress

A homogeneous perfect-crystal slip coordinate driven directly over a DFT-scale gamma-surface barrier corresponds to ideal shear stresses of order gigapascals, consistent with first-principles ideal-strength calculations for Al.

Therefore this project should **not** attempt to explain a tens-of-MPa resolved cyclic stress by simply lowering the gamma-surface barrier until the simulation slips. That would amount to hidden fitting.

Instead, the low-stress fatigue problem must identify how microscopic mechanics changes the *local* generalized force or the accessible state space. Candidate mechanisms to derive and test are:

1. free-surface stress concentration and surface-modified gamma surfaces;
2. pre-existing non-affine defect fields;
3. spatially correlated multi-slip coordinates;
4. finite-temperature microscopic initial conditions and exact bath memory;
5. local geometric amplification near persistent slip structures;
6. interaction of slip with the spacing distribution $P(a,t)$.

Each mechanism must be introduced with its microscopic definition and independently calculable parameters.

## Recommended replacement of the one-harmonic potential

The current approximation

$$
V_\gamma(s)=\frac{\Delta_\gamma}{2}
\left[1-\cos\left(\frac{2\pi s}{b}\right)\right]
$$

should ultimately be replaced by

$$
V_\gamma(\mathbf s)=A_{\rm RA}\,\gamma_{\rm Al}(\mathbf s),
$$

where $\gamma_{\rm Al}(\mathbf s)$ is an atomistically calculated two-dimensional generalized-stacking-fault energy surface and $A_{\rm RA}$ is the mechanically defined representative slip area.

The conjugate resolved shear force is then

$$
\mathbf F_s=A_{\rm RA}\,\boldsymbol\tau,
$$

and local mechanical stability is determined by the Hessian of the driven landscape

$$
\Phi(\mathbf s,t)=A_{\rm RA}\gamma_{\rm Al}(\mathbf s)
-A_{\rm RA}\boldsymbol\tau(t)\cdot\mathbf s.
$$

Loss of a local minimum occurs when the smallest eigenvalue satisfies

$$
\boxed{
\lambda_{\min}\left[\nabla_{\mathbf s}^2\Phi\right]=0.
}
$$

This is a mechanically cleaner instability condition than prescribing a yield stress.

---

# 한국어 번역

## 목적

현재 비선형 Hamiltonian slip-bath 모델은 경험적 피로식을 넣지 않아도 주기적인 비아핀 에너지 지형에서 basin 간 이동과 cycle-to-cycle 구조변화가 가능하다는 것을 보여준다. 다음 단계는 임의의 무차원 barrier를 실제 Al의 원자수준 에너지 지형으로 교체하는 것이다.

## 문헌 제약조건

FCC Al의 generalized stacking-fault 계산에서는 계산방법과 relaxation 조건에 따라 차이가 있지만 unstable stacking-fault energy가 대략 $10^{-1}\,\mathrm{J/m^2}$ 규모이며, 대표적인 first-principles 자료에서 약 $175$–$224\,\mathrm{mJ/m^2}$ 수준의 값들이 보고된다. intrinsic stacking-fault energy 역시 대략 $120$–$166\,\mathrm{mJ/m^2}$ 범위의 계산값들이 흔하다.

이 차이는 fitting으로 없애야 할 오차가 아니라 계산방법, 변형상태, relaxation, slip path의 물리적 차이이므로 향후 모델에서 입력값의 출처를 명시해야 한다.

## 저응력 피로에 대한 중요한 결론

완전결정의 균일한 slip 좌표가 DFT 수준의 gamma-surface barrier를 직접 넘으려면 이상전단강도 수준, 즉 GPa 규모의 국부응력이 필요하다. 이는 Al의 first-principles ideal-strength 계산과도 일관된다.

따라서 수십 MPa 수준의 macroscopic cyclic stress에서 slip을 만들기 위해 gamma barrier를 임의로 낮추면 안 된다. 그것은 사실상 숨은 fitting이다.

대신 다음 항목들이 실제 국부 구동력 또는 접근 가능한 상태공간을 어떻게 바꾸는지 미시역학으로 유도해야 한다.

1. 자유표면의 응력집중과 surface-modified gamma surface;
2. 기존 비아핀 결함장;
3. 공간적으로 상관된 multi-slip 좌표;
4. 유한온도 미시 초기조건과 정확한 bath memory;
5. persistent slip 구조 주변의 국부 기하학적 증폭;
6. slip과 원자간격 분포 $P(a,t)$의 결합.

현재의 단일 cosine potential은 최종적으로 실제 Al의 2차원 $\gamma$-surface로 교체해야 한다.

$$
V_\gamma(\mathbf s)=A_{\rm RA}\gamma_{\rm Al}(\mathbf s).
$$

외부 resolved shear stress를 포함한 구동 에너지 지형은

$$
\Phi(\mathbf s,t)=A_{\rm RA}\gamma_{\rm Al}(\mathbf s)
-A_{\rm RA}\boldsymbol\tau(t)\cdot\mathbf s
$$

이고, 국부 최소점의 기계적 소실은

$$
\boxed{
\lambda_{\min}\left[\nabla_{\mathbf s}^2\Phi\right]=0
}
$$

으로 잡을 수 있다. 이는 경험적인 yield stress를 집어넣는 것보다 훨씬 역학적으로 명확한 기준이다.
