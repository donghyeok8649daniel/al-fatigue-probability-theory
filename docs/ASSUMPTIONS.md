# Assumptions and approximations

This file is normative for the active multiplicity-free multilayer theory.
Earlier normal-chain and two-row documents are retained as derivation history,
not as competing active total energies.

## Active physical assumptions

1. The target is a pure single crystal under repeated uniaxial tensile stress.
2. The microscopic state is `(a,s)`: one locally homogeneous normal layer
   spacing and one scalar collective/unwrapped registry for a declared slip
   system.  A 2D probability state does not imply 2D continuum mechanics.
3. The intrinsic potential is the generalized-LJ multilayer sum
   `U0(a,s)=sum_{k>=1} W(k*a,s)`, with no `k` multiplicity and with the same
   `s` for every normal layer.
4. Full absolute convergence requires `m>n>2`.  A single row requires only
   `q>1`; registry-excess reciprocal terms converge exponentially once the
   zero Fourier mode is removed.
5. Normal opening and slip are the exact decomposition of one `U0`.  The
   historical collinear `U_infinity` and row kernel `W` are not added.
6. `U0` contains atomic positional energy only.  The single applied stress
   enters as `Q_a=A0*sigma(t)` and `Q_s=A0*M*sigma(t)` in the probability
   current.  No independent shear-fatigue load is used.
7. Eliminated atomic/phonon coordinates form an isothermal bath, velocity
   relaxation is faster than resolved evolution, mobility is constant, and
   diffusion obeys `D_i=kBT*M_i`.  These are reduction assumptions requiring
   atomistic validation; no extra empirical diffusion is inserted.
8. No named family is imposed on `P(a,s,t)`.  A bonded-basin Gibbs density is
   only a conditional/metastable ensemble, never a global dead-load tensile
   equilibrium.
9. Finite-rate intrawell lag is not plasticity.  Plasticity requires a
   residual `Delta<z> != 0` after unloading/relaxation in
   `s=s0+z*b+s_tilde`.
10. Mean intrinsic energy and cumulative hysteresis are distinct.
    `E_hyst=integral dot(D)_irr dt` is nondecreasing dissipation but is not
    assumed to be entirely stored damage energy.
11. Crack initiation is first passage through the outer negative-curvature
    root of `partial_a U0=Q_a`, including relative-flux corrections when the
    boundary moves.  No arbitrary spacing or hysteresis-energy threshold is
    active.
12. `A0`, correlation area, slip homogenization thickness, and FEM element
    area remain distinct unless independently derived.
13. 2D/3D CAD/FEM remains geometry, mesh, scalar normal-stress transport, and
    visualization.  It does not activate shear or multiaxial constitutive
    physics.
14. EAM/DFT is future quantitative validation/extension only and does not
    replace the current generalized-LJ fundamental potential.

## Exact results versus numerical controls

The Mellin--Poisson identity, the Bessel--Lambert `H_q` representation, the
normal/slip identity decomposition, and the 12--6 polylog closure are exact.
`pmax`, `kmax`, reciprocal modes, layer modes, grid spacing, timestep, and
solver tolerance are numerical controls and must be refined.  They are not
material fitting parameters.

## Unresolved physical inputs

- `epsilon_LJ`, `sigma_LJ`, `b`, the reference `(a0,s0)`, and the chosen
  single-crystal slip system;
- representative mechanical area `A0`, mobilities/memory times, temperature,
  and `h_slip`;
- dislocation storage, hardening, backstress, and multiple-slip interactions;
- experimental/atomistic validation of the outer-barrier first-passage event
  as observed crack initiation.

---

# 가정과 근사의 한국어 정리

이 문서는 현재 활성화된 multiplicity-free 다층 이론의 기준 문서다. 과거의
normal-chain 및 two-row 문서는 유도 이력으로 보존하지만 현재 total energy와
경쟁하는 식으로 사용하지 않는다.

## 활성 물리 가정

1. 대상은 반복 단축 인장응력을 받는 순수 단결정이다.
2. 미시상태는 국소적으로 균일한 normal layer 간격 `a`와, 지정된 slip
   system의 하나의 scalar collective/unwrapped registry `s`다. `(a,s)`가
   2차원 확률공간이라는 사실은 2D continuum constitutive law를 뜻하지 않는다.
3. intrinsic potential은 `U0(a,s)=sum_{k>=1} W(k*a,s)`다. 앞에
   multiplicity `k`가 없고 모든 normal layer 항에서 같은 `s`를 사용한다.
4. 전체 absolute sum은 `m>n>2`를 요구한다. 단일 row는 `q>1`이면 되고,
   slip-excess는 zero Fourier mode가 제거되어 지수적으로 수렴한다.
5. normal opening과 slip은 하나의 `U0`에서 정확히 분해한다. 과거의
   collinear `U_infinity`와 row kernel `W`를 total energy로 더하지 않는다.
6. `U0`에는 외력 일을 넣지 않는다. 단 하나의 인장응력으로부터
   `Q_a=A0*sigma(t)`, `Q_s=A0*M*sigma(t)`를 만들어 확률 current에 넣는다.
   독립적인 shear fatigue 입력은 없다.
7. 생략한 원자/phonon 좌표가 등온 bath이고 속도완화가 빠르며 mobility는
   상수라고 가정한다. 확산은 fluctuation--dissipation으로 고정한다.
8. `P(a,s,t)`에 Gaussian 등 특정 분포족을 강제하지 않는다. bonded-basin
   Gibbs 분포는 조건부 metastable ensemble일 뿐 전역 인장평형이 아니다.
9. well 내부 phase lag는 소성이 아니다. `s=s0+z*b+s_tilde`에서 unloading과
   relaxation 뒤 `Delta<z> != 0`일 때만 잔류 소성으로 정의한다.
10. 평균 intrinsic energy와 누적 hysteresis를 구분한다. `E_hyst`는 비가역
    dissipation이지만 전부 저장 damage energy라고 단정하지 않는다.
11. 균열개시는 `partial_a U0=Q_a`의 음의 곡률 외측 장벽을 통과하는
    first passage다. 임의 거리나 energy threshold를 사용하지 않는다.
12. `A0`, 상관면적, slip 균질화 두께, FEM 요소면적은 서로 다르다.
13. 2D/3D CAD/FEM은 geometry, mesh, scalar normal stress 전달 및
    visualization용이다. shear/multiaxial 재료이론을 활성화하지 않는다.
14. EAM/DFT는 미래 정량 검증/확장용이며 현재 generalized-LJ governing
    potential을 대체하지 않는다.

수학적 항등식과 numerical cutoff를 혼동하지 않는다. `pmax`, `kmax`, Bessel
mode, grid, timestep은 모두 수치 수렴변수이며 재료 fitting parameter가 아니다.
