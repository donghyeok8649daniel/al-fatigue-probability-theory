# AGENTS.md — 이론, 구현, 검증 및 작업 원칙

## 0. 가장 먼저 읽을 것

이 저장소는 알루미늄 피로의 **결정론적 확률밀도 진화**와 원자적 에너지
지형을 연구한다. 보기 좋은 소성 변형이나 균열 확률을 만드는 것이 목표가 아니다.
수학적 일관성, 수치 검증, 물리적 정직성이 성공 기준이다.

- 사용자의 최신 지시가 현재 작업 범위를 결정한다. 오래된 미션을 임의로 다시 시작하지 않는다.
- 이 파일은 프로젝트 지침이다. 실행 중인 모델과 수치는 실제 코드/데이터로 재확인한다.
- 재개할 때 먼저 같은 폴더의 **CURRENT_WORK_HANDOFF.md**를 읽는다.
- 검증된 결과, 실험적 구현, 보정 후보, 미완료 계산을 구별한다.
- 이 파일만 읽고 모든 모델이 검증되었다고 가정하지 않는다.
- 한국어로 진행 상황을 간결하게 알린다. 긴 계산 중에도 상태를 알린다.
- 사용자가 요청하지 않으면 추가 에이전트를 생성하지 않는다.
- 추론 강도 변경은 사용자/실행 환경 설정이다. 에이전트가 임의로 변경했다고 주장하지 않는다.

## 1. Git과 사용자 작업 보호

대상 저장소: donghyeok8649daniel/al-fatigue-probability-theory
대상 작업 브랜치: probability-pde-solver-v1

작업 시작 전 실제 경로, HEAD, 원격 HEAD, 작업 트리를 확인한다:

    git fetch origin
    git rev-parse HEAD
    git rev-parse origin/probability-pde-solver-v1
    git branch --show-current
    git status --short
    git worktree list

- fetch가 실패했으면 캐시된 원격 참조를 새로 검증된 원격 HEAD라고 말하지 않는다.
- 예상 HEAD와 다르면 차이를 보고한다. 최신 원격을 보존하고 안전하게 이어간다.
- main 수정/병합 금지. force push, destructive reset, 사용자 파일 덮어쓰기 금지.
- detached worktree가 있으면 실제 작업 위치와 커밋 기반을 확인한다.
  원래 폴더와 실험 worktree의 코드가 같다고 가정하지 않는다.
- 기존 변경/미추적 파일은 사용자 작업으로 취급한다. 깨끗하게 만들려고 삭제하지 않는다.
- 변경에는 apply_patch를 사용한다. 검색은 rg / rg --files를 우선한다.
- 임시/개인 절대 경로, 개인 설정, 인증 정보를 커밋하지 않는다.
- 커밋 전 변경 파일을 읽고 git diff --check 및 관련 검증을 수행한다.
- push 직전 다시 fetch하여 분기를 확인하고 정상 fast-forward push만 한다.
- 미완료 과학 결과를 완성된 보정처럼 커밋/보고하지 않는다.
- 백업/인계 문서에는 완료, 미완료, 실제 테스트 결과와 미확인 결과를 구분한다.

## 2. 모델 계층 — 서로 바꾸어 부르지 말 것

| 모델 | 의미 / 상태 |
|---|---|
| TwoRowLJ | 검증된 reduced two-row / projected lattice reference |
| AnalyticLJEAM | LJ pair를 보존하는 reduced analytic LJ–EAM hybrid |
| Reduced Al-target best-feasible hybrid | 정적 Al 타깃을 사용했지만 식별성/표현력 한계가 있는 reduced 후보 |
| FCC111 full stack | 삼각 평면과 ABC 적층을 반영한 별도 homogeneous static reference |
| Full-FCC calibrated candidate | full geometry에서 다시 보정해야 하는 후보; reduced 파라미터 전이와 다름 |
| FCC111 active interface | 한 인터페이스만 열리고 미끄러지는 별도 static research reference |

기존 two-row W(a,s)는 정확한 full FCC 결정 에너지가 아니다.
더 완전한 기하학이라고 자동으로 실험적으로 정확해지지 않는다.
어느 연구 모델도 검증 없이 production default나 PDE 입력으로 승격하지 않는다.

UI 모델 추적의 기준은 solver_v1/energy_model_registry.py의 실제 생성 경로다:

    app.desktop_ui -> app.solver_adapter
      -> build_energy_model -> energy/gradient/Hessian -> probability PDE

현재 registry ID:
- two_row_lj_reference
- analytic_lj_eam_hypothetical
- al_target_best_feasible_hybrid

UI 레이블만 보고 사용 모델을 추론하지 않는다. 결과, 내보내기, 수렴 인증에
실제 class, model ID, parameter source, scale, calibration status를 기록한다.

## 3. 확률 이론의 기본 구조

N=1의 상태는 q=(a,s)이다. a는 수직 분리, s는 누적 registry/slip 좌표다.

    ∂t P = -∂a J_a - ∂s J_s
    J_i = -M_i [P ∂i G + kT ∂i P]
    M = diag(M_a, M_s)

기존 canonical N=1 기준값은 M_a=1, M_s=0.05, kT=0.02이다.
서로 다른 에너지 단위를 쓰는 모델은 온도/에너지 스케일을 별도로 확인한다.

- 실제 solver의 conservative Scharfetter–Gummel(SG) flux를 보존한다.
- 명시/암시 시간적분은 같은 generator를 근사해야 한다.
- 안정성, positivity, mass balance, Gibbs stationary state를 독립적으로 시험한다.
- 작은 음수 density 보정은 별도 수치 진단이며 물리적 흡수로 합산하지 않는다.
- Monte Carlo trajectory counting, empirical S–N, J2, Ramberg–Osgood,
  임의 hardening/yield cutoff를 이론의 대체물로 추가하지 않는다.
- 더 높은 차원/TT 구현을 N=1과 동일 수준으로 검증되었다고 가정하지 않는다.

## 4. Registry 분해와 변형률

    s = b n + xi
    n = floor(s/b + 1/2)
    -b/2 <= xi < b/2

n은 정수 registry-well index다. 음의 s와 정확한 ±b/2에서 규약을 시험한다.
n!=0 또는 epsilon_p!=0 자체가 실험적으로 분해 가능한 Al 소성의 증거는 아니다.

수치 intact density의 질량을 M_raw=∫P da ds라 하자.
물리적 생존 확률과의 차이는 수치 잔차로 분리한다.
잔존 ensemble의 구성응답은 실제 살아 있는 density에 대해 조건부 평균한다:

    <X>_surv = (∫ X P da ds) / M_raw
    epsilon_a  = <a-a0>_surv / a0
    epsilon_xi = chi <xi>_surv / a0
    epsilon_p  = chi b <n>_surv / a0
    epsilon_total = epsilon_a + epsilon_xi + epsilon_p

연속 방정식이 정확하면 M_raw=S_local이다. 수치적으로 다를 때 분모를
조용히 바꾸지 말고 기존 의미를 감사하여 보고한다.
생존 질량이 소진된 경우 조건부 평균은 정의되지 않으므로 0 strain으로 위장하지 않는다.
Unnormalized survivor moment ∫XP와 conditional moment를 구별한다.

매 저장 시점:

    R_epsilon = epsilon_total - epsilon_a - epsilon_xi - epsilon_p

의 최대 절댓값을 확인한다. xi와 n을 중복 합산하지 않는다.
흡수된 probability를 survivor strain에 포함하지 않는다.
소성 변형률을 clipping, artificial scaling 또는 specimen area로 확대하지 않는다.

## 5. Well population, flux, opening의 정확한 분리

    Omega_n = [(n-1/2)b, (n+1/2)b)
    P_n = ∫_{Omega_n} P da ds
    Σ_n P_n = M_raw
    J_{n+1/2} = ∫ J_s(a,(n+1/2)b) da

양의 J는 increasing s 방향이다. 모든 표현된 well을 포함한다.
Opening sink를 dot A_n이라 하면 reflecting outer-s domain에서:

    dot P_n = J_{n-1/2} - J_{n+1/2} - dot A_n
    M_n = Σ_n n P_n
    dot M_n = Σ_interfaces J_{n+1/2} - Σ_n n dot A_n

열린 outer-s 경계나 수치 repair가 있으면 별도 항을 명시한다.
정확한 연속 생존 질량 S에서:

    dot epsilon_p,flow = chi b/(a0 S) Σ_interfaces J
    dot epsilon_p,selective =
        chi b/(a0 S) [-Σ_n n dot A_n + <n>_surv Σ_n dot A_n]

두 번째는 **selective opening contribution — NOT plastic flow**다.
Conditional normalization의 시간 미분 항을 빠뜨리지 않는다.
실제 discrete balance와 snapshot 차분/시간평균의 차이를 구별한다.

SG 인터페이스 양방향 분해:

    j_net = j_plus - j_minus
    j_gross = j_plus + j_minus

j_gross를 abs(j_net)로 대체하지 않는다.
대칭 열활동은 gross>0, net≈0, epsilon_p≈0를 동시에 만들 수 있다.
다만 smooth continuum density에서 nearest-neighbor SG gross는 대략
2 D_s P_s(face)/ds로 발산하는 grid recrossing traffic이다. 이를 mesh-independent
committed hopping 횟수로 부르지 않는다. 물리 transition count에는 well core /
committor와 별도의 reactive-flux 유도가 필요하다. 기존 signed SG balance는 유지한다.

## 6. 소성의 증거 수준과 수렴

1. floating-point nonzero: 증거 아님.
2. 수렴한 gross interwell activity: configurational slip activity.
3. 수렴한 directional net transfer: model plastic flow.
4. unload + 충분한 zero-stress hold 후 남는 수렴한 registry shift: residual model plasticity.

수렴 확인: a-grid, s-grid, dt, 가능한 explicit/implicit 비교, well-face alignment,
well/domain 크기(예: 3/5/7 wells), flux/population balance, symmetry error.
Normal 및 intrawell relaxation 전에 residual plasticity를 선언하지 않는다.

이전 O(1e-11) apparent plastic strain이 aligned/refined grid에서 소멸한 결과를 보존한다.
해당 값은 검증된 소성으로 재해석하지 않는다.
epsilon_p, mean n, transferred mass, boundary flux의 수치 floor는 각각
실제 resolution changes/잔차/repair로 추정한다. 임의 threshold를 만들지 않는다.

반응은 분리한다:
- a intrabasin: normal reversible/dynamic response;
- xi intrawell: configurational/anelastic memory;
- n-changing net transfer: model plastic slip;
- a opening dividing-surface crossing: crack first passage.

## 7. Opening probability와 specimen aggregation

    A_abs = 누적 실제 opening boundary 흡수 질량
    P_init_local = A_abs
    S_local = 1 - A_abs
    R_mass = M_raw + A_abs - 1

A_abs는 nonnegative/nondecreasing, S_local는 nonincreasing이어야 한다.
흡수 increment와 시간적분 opening flux를 비교한다.
Raw 1-M_raw, mass residual, negative correction, flux discrepancy는 물리 확률이 아니다.
Mask 이동으로 제거되는 질량은 physical dividing-surface sweep인지
수치 truncation/repair인지 명시한다.

Rare-event floor는 grid/dt 변화, flux discrepancy, residual, correction으로 산정한다.
O(1e-16)을 단지 double precision의 nonzero라는 이유로 physical crack이라 하지 않는다.
하나의 정밀 실행만으로 convergence-certified라 하지 않는다.

    N_eff = A_stressed / A_c
    log S_spec = N_eff log1p(-P_init_local)
    P_spec = -expm1(log S_spec)

A_c는 외부 보정이 필요한 statistical correlation area다.
독립적·등가적 region 가정은 미검증이다. 임의 physical default를 정하지 않는다.
Local floor 이하 신호를 N_eff로 확대해 certified probability로 제시하지 않는다.
Mathematical extrapolation은 명확한 uncertified 레이블로 별도 표시할 수 있다.

    A_atomic_cell != A_c != mesh/element area

A_atomic_cell은 명시된 원자 에너지 per-area 변환에만 사용한다.
A_c는 stress mapping, energy, density, strain, mobility에 절대로 들어가지 않는다.

## 8. 검증된 reduced static stress calibration

TwoRowLJ reference(기본 chi=0.2, b=1)의 역사적 검증값:

    a0 = 0.7713438268704838
    W_aa = 122.78696473236548
    W_as ≈ 0
    W_ss = 50.347580146252554
    c = [1, chi]^T
    kappa_axial = a0 / (c^T H0^(-1) c)
                = 86.29296488740997
    f* = kappa_axial * sigma/E

응력/영률 단위를 일치시킨다(MPa와 GPa 혼동 금지).
이는 zero-frequency relaxed tangent다. 유한 주파수 strain이 sigma/E를
따르도록 kappa/M_a/M_s/chi를 조정하지 않는다.
Hybrid/full 모델은 자기 equilibrium와 Hessian으로 별도 kappa를 계산한다.
과거 a0/kappa/parameter 파일을 덮어쓰지 않는다.

주기 에너지 W(a,s+b)=W(a,s)여도 n→n+1은 누적 translation b를 뜻한다.
동일 local registry energy가 동일 macroscopic slip history를 뜻하지 않는다.
현재 minimal model에 물리 hardening이 존재한다고 주장하지 않는다.

## 9. Fast/slow dynamics와 finite-temperature reduction

기존 full overdamped 모델의 local linear diagnostic:

    R = sqrt(M) H0 sqrt(M)
    lambda_i = eigvalsh(R), tau_i = 1/lambda_i
    G(omega) = (kappa/a0) c^T (i omega I + M H0)^(-1) M c

Full Hessian을 사용한다. Eigenvector를 순수 a/s라고 가정하지 않는다.
Finite-a phase lag는 plastic hysteresis가 아니다.
M_a sweep은 singular-limit test이지 물리 mobility fitting이 아니다.

고정 reflecting fast-a domain:

    rho_eq(a|s,f) = exp(-G/kT) / Z_a
    Z_a = ∫ exp(-G/kT) da
    F_eff = -kT log Z_a
    ∂s F_eff = <G_s>_rho
    ∂t P_s = ∂s [M_s(P_s ∂s F_eff + kT ∂s P_s)]

F_eff는 free energy다. A_eff라는 기호를 사용하지 않는다.
이동 경계 L(s), U(s)에서는:

    ∂s F_eff = <G_s>_rho
       - kT [exp(-G(U)/kT) U_s - exp(-G(L)/kT) L_s] / Z_a

기계적 minimum a*는 finite-T conditional <a>의 대체물이 아니라 reference다.
Low-T Laplace:

    F_eff ≈ G(a*,s) + (kT/2) log(G_aa(a*,s)/(2π kT))

이는 코드의 dimensionless integration measure에 해당하며 dimensional measure이면
참조 척도와 additive constant를 명시한다.
Fixed-domain Gibbs, truncated-basin Gibbs, absorbing-process QSD는 구분한다.
Reduced crack sink는 fast absorbing operator/flux로 검증하지 않으면 production-valid가 아니다.

## 10. 물리 시간 — 현재 사용 불가

    a*=a_phys/L0, s*=s_phys/L0, G*=G_phys/E0, t*=t_phys/t0
    P*=L0^2 P_phys                         (N=1, 동일 길이 척도)
    M_i* = t0 E0/L0^2 M_i,phys
    [M_i,phys] = length^2/(energy·time)
    t0 = M_i* L0^2/(M_i,phys E0)
    t_phys = t0 t_model
    f_Hz = f_model/t0

한 t0에는 두 mobility가 일관된 ratio를 가져야 한다.
Energy/length만으로 overdamped time은 결정되지 않는다.
Atomic mass의 Newtonian time을 몰래 삽입하지 않는다.
M_phys는 명확히 대응되는 collective-coordinate MD relaxation,
coupled correlation modes 또는 diffusivity D=M kBT로만 보정한다.
Bulk self-diffusion을 근거 없이 a/s diffusivity로 대체하지 않는다.

현재 M_a,phys, M_s,phys, t0는 unavailable.
검증된 kinetic calibration이 로드되지 않으면 physical seconds/Hz는 disabled.
UI, plots, exports는 model time / cycles per model time을 사용한다.
Static Al calibration은 kinetic calibration이나 fatigue lifetime validation이 아니다.

## 11. LJ–EAM hybrid의 수학적 계약

LJ pair functional form은 유지한다:

    phi_LJ(r)=4 epsilon_LJ[(sigma_LJ/r)^12-(sigma_LJ/r)^6]
    rho_kernel(r)=rho0 exp[-beta(r/r_e-1)] = C exp(-kappa_rho r)
    C=rho0 exp(beta), kappa_rho=beta/r_e

rho는 environmental density이며 probability P와 다르다.
Reduced two-row site counting을 보존하고 full geometry로 무심코 전이하지 않는다.
General differentiable F(rho)에 대해:

    W_emb,i = F' rho_i
    W_emb,ij = F'' rho_i rho_j + F' rho_ij

site multiplicity를 별도로 곱한다.
Square-root F=-A sqrt(rho/rho_ref), analytic extension은 검증 후 최소한만.
EAM pair/embedding gauge와 density-scale redundancy를 분석한다.
F의 상수도 cohesion을 비교할 때 separated-atom reference와 일관돼야 한다.
예: F(x)=-A sqrt(x)+B(x-1)이면 F(0)=-B이다.

Canonical lattice는 infinite Poisson/Bessel 표현이다.
Finite real-space cutoff는 독립 검증용이며 production 이론으로 대체하지 않는다.
계열 truncation, reciprocal/layer tails와 derivative error를 노출한다.

## 12. Full FCC(111) reference

    e1=(1,-1,0)/sqrt(2)
    e2=(1,1,-2)/sqrt(6)
    e3=(1,1,1)/sqrt(3)
    b=a_lat/sqrt(2), h111=a_lat/sqrt(3)
    a1=b e1, a2=b(e1/2+sqrt(3)e2/2)
    A_atomic_cell=sqrt(3)b^2/2
    tau=(a1+a2)/3, 3 tau∈Lambda

ABC 위상을 보존하고 scalar path u=s e_s를 명시한다.
direct_110와 Shockley_112 경로를 혼동하지 않는다.
Partial translation과 full registry period는 다르다.

**+ABC 방향 감사(nonlocal_v5):** 기존 e1/e2/e3는 기하학적 구성 프레임이다.
`+tau,+h`로 생성된 FCC의 실제 cubic axes에 대한 plane basis는
`(-e1,-e2,e3)`이다. `plane_basis_in_stacked_cubic_axes()`를 사용해 탄성
텐서와 cubic Gamma-X/L/K 경로를 변환한다. 기존 원자 위치/에너지는 바꾸지 않는다.
잘못된 축 대응은 일부 acoustic matrix의 약37–40% 오차를 만들었다.
RegistryPath.direction_3d는 기존 기하학 프레임 벡터로 보존한다.
과거 finite-q CSV의 cubic path 이름은 이 감사 이전 기록이며, 수정된 경로 검사는
nonlocal_v5/finite_q_corrected_paths.csv에 별도 저장한다.

    R_mn=m a1+n a2
    |R_mn|^2=b^2(m^2+mn+n^2)
    r_mnl^2=(l a)^2+|R_mn+Delta_l|^2
    E_pair/atom=(1/2) Σ' phi(r_mnl)
    rho_full/atom=Σ' rho_kernel(r_mnl)
    E_emb/atom=F(rho_full)

자기항 제외. Density를 모두 합산한 다음 site마다 F를 적용한다.
Σ_l F(rho_l), Σ_l old W_hybrid(l a,Delta_l), old normal-chain+full-stack은 금지한다.

2D Fourier convention exp(-i G·r), reciprocal ai·bj=2π delta_ij:

    f_hat_p(G)=2π/Gamma(p) (G/(2d))^(p-1) K_(p-1)(dG)
    f_hat_p(0)=π/(p-1) d^(2-2p)
    S_p(d,delta)=(1/A_atomic_cell) Σ_G f_hat_p(G) exp(i G·delta)

Positive-layer mean sum은 zeta(2p-2):
LJ p=3 -> zeta(4), p=6 -> zeta(10).
Reduced row의 zeta(5),zeta(11)와 다르다.
1D Bessel order p-1/2를 2D의 p-1로 혼동하지 않는다.

Exponential density의 2D transform:

    q=sqrt(kappa_rho^2+|G|^2)
    f_hat_exp(G)=2π kappa_rho exp(-d q)(1+d q)/q^3

Full-FCC homogeneous reference의 Delta_l=l(tau+s e_s)는
local interface slip과 같지 않다. 기존 경로 선택/geometry를 검증 없이 바꾸지 않는다.

## 13. Active interface — 아직 static research

    lower l<=0: z_l=l h, u_l=0
    upper l>=1: z_l=l h+delta_a, u_l=s e_s
    a=h+delta_a

Lower-lower와 upper-upper는 정확히 상쇄된다.
Cross-interface layer separation k의 pair multiplicity는 k다:

    Delta E_pair = Σ_(k>=1) k[
        w_plane(kh+delta_a,k tau+s e_s)-w_plane(kh,k tau)]

깊이 r의 atom마다:

    Delta rho_r = Σ_(k>=r+1)[rho_plane(active)-rho_plane(bulk)]
    Delta E_emb = 2 Σ_(r>=0)[F(rho_bulk+Delta rho_r)-F(rho_bulk)]
    W_int = Delta E_pair + Delta E_emb

대칭 두 half-crystal인 경우에만 위 factor 2를 사용한다.
하나의 global crystal density에 F를 적용하지 않는다.
Infinite bulk energy를 수치적으로 크게 더했다 빼지 않는다.

검증:
W_int(h,0)=0, equilibrium gradients≈0, Hessian/stability,
registry periodicity, GSF path, separated half-crystal limit,
reciprocal/direct 및 neighborhood/tail convergence.

    gamma(s)=[W_int(h,s)-W_int(h,0)]/A_atomic_cell
    W_sep=lim_(a->infinity) W_int(a,0)
    W_sep/A_atomic_cell = 2 gamma_surface

마지막 관계는 같은 orientation, 같은 surface relaxation convention일 때 비교한다.
Unrelaxed interface와 relaxed DFT surface target의 차이를 숨기지 않는다.
Separation work는 finite-load saddle barrier와 동일하지 않다.

별도 physical interface work의 후보:

    G_int=W_int-A_atomic_cell[T_n delta_a_phys+tau_res s_phys]
    T_n=sigma(e·n_slip)^2
    tau_res=sigma(e·n_slip)(e·m)

이는 active-interface의 차원 있는 work다.
기존 reduced f*=kappa sigma/E를 몰래 이 식으로 교체하지 않는다.
A_c는 어느 식에도 들어가지 않는다.

Small-strain slip의 axial projection:

    epsilon_p_tensor=(gamma/2)(m⊗n_slip+n_slip⊗m)
    gamma=s/h
    epsilon_axial=(s/h)(e·m)(e·n_slip)
    chi_candidate=(a0/h)(e·m)(e·n_slip)

길이 단위를 일치시킨다. Orientation/좌표 mapping이 확정되기 전 default chi를 바꾸지 않는다.

## 14. Calibration과 identifiability

- 문헌 수치는 source title/authors/year/DOI 또는 stable ID, method,
  temperature, orientation, units와 함께 저장한다. 값을 만들어 넣지 않는다.
- 0 K DFT/EAM과 room-temperature experiment를 무표기 혼합하지 않는다.
- per atom, reduced 2-atom cell, per interface cell, J/m^2를 구분한다.
- 직접 대응, projection, coarse-grained, 현재 unrepresentable을 분류한다.
- Vacancy 없는 상태공간을 vacancy validation이라 부르지 않는다.
- Homogeneous shear를 local GSF, uniform expansion을 cleavage로 부르지 않는다.
- Parameter gauges를 먼저 제거하고 deterministic staged fit을 수행한다.
- 정규화 residual의 scale/uncertainty를 공개하고 held-out을 loss에서 제외한다.
- 모든 deterministic starts, bounds, losses, individual residuals, SVD를 저장한다.
- 실제 최적화를 재실행하는 코드와 저장된 fit vector로 결과를 재생하는 코드를 구분한다.
- 유한차분 Jacobian의 작은 singular value를 structural identifiability로 오해하지 않는다.
- 관측량의 analytic dependencies를 먼저 분석한다.
- 두 후보의 실패만으로 전체 analytic family 불가능을 증명했다고 주장하지 않는다.
- 추가 F 항은 정량적 baseline 실패와 gauge/rank 검토 이후 최소한만 허용한다.
- Stable bulk fit만으로 GSF/decohesion/kinetics가 검증되었다고 말하지 않는다.

특히 현재 [111] normal alpha / equal-biaxial beta 변형에서는:

    X=C11+2C12
    Haa=V(X+4C44)/3
    Hbb=4V(X+C44)/3
    Hab=2V(X-2C44)/3
    Hbb=2Haa+Hab

세 Hessian 성분은 **두 개 독립 cubic elastic 조합**만 제약한다.
C11-C12는 여기서 식별되지 않는다.
Perfect cubic configuration의 normal/lateral force도 독립 두 정보가 아니다.
Numerical full rank만으로 반대 결론을 내리지 않는다.
audited_v2는 engineering simple shear gamma=s/h를 추가해 독립 모드를 사용한다:

    H_hydro = 3V(C11+2C12)
    H_normal = V(C11+2C12+4C44)/3
    H_shear = V(C11-C12+C44)/3

이는 homogeneous elasticity이고 local GSF가 아니다.
불완전한 bulk fit은 b와 h를 함께 isotropic relaxation하여 실제 FCC equilibrium을
찾는다. Potential/L0/rho_ref를 strain 동안 바꾸지 않는다. 실제 평형의 원자
면적·체적으로 J/m^2와 GPa를 환산한다.
Bulk positive Hessian만으로 surface/cleavage 안정성을 인증하지 않는다.
Negative work of separation 후보는 물리적으로 채택하지 않는다.
좁은 첫 국소 장력을 놓칠 수 있는 전체 구간 unimodal 최적화를 branch tracking
대신 사용하지 않는다. Normal-only saddle은 coupled a-s saddle과 구별한다.
현재 결과에 충돌하는 주장이 있으면 인계 기록/원시 CSV를 보고 재감사한다.

## 15. UI 계약

- 구현은 한 벌, 중앙 app/i18n.py 또는 실제 localization layer의 stable keys 사용.
- ko/en 키의 완전성을 시험하고 기본 언어는 한국어.
- 언어 변경은 PDE 재실행, 입력/결과 변경, graph 선택/zoom/pan 초기화를 하지 않는다.
- 모델, 품질, 시간 기준, physical calibration status를 명시한다.
- Preview alone는 convergence certification이 아니다.
- Probability resolution과 plasticity resolution은 서로 다른 상태다.
- Tiny plastic strain은 별도 scientific-notation/detail plot으로 표시한다.
  원래 데이터를 증폭하지 않는다.
- Stress input을 임의 clamp하지 않는다. 큰 sigma/E는 해석 경고만 낸다.
- Ac/Astressed 변경은 post-processing이며 PDE를 다시 풀지 않는다.
- 긴 Pre/Solve 폼은 `app/scrollable_panel.py`의 세로 viewport로 접근한다.
  Solve/수렴 버튼은 스크롤 바깥 고정 action row에 둔다. 전역 wheel 바인딩으로
  Matplotlib zoom/다른 위젯을 가로채지 않는다. 스크롤/resize/Tab/언어 전환은
  입력값·결과·PDE 실행을 바꾸지 않아야 한다. 작은 창에서 버튼 가시성을 시험한다.

사용자가 요청한 후속 UI는 geometry import / 기본 편집 가능한 원통 시편 ->
meshing -> pre조건 -> solving/posting 순서다. 다만 각 솔버 단계의 실제 실행과
물리 검증을 먼저 하고, 솔버 준비 상태를 보고한 뒤 **사용자 확인을 받아야**
UI 재설계로 넘어간다. 현재 active interface는 그 gate를 통과하지 못했다.
자세한 기준은 solver_v1/SOLVER_VALIDATION_GATES.md를 따른다.
Mesh colormap은 실제 spatial mechanics/local probability 연결로 계산한 field만
사용한다. Global 확률 하나를 공간 예측처럼 칠하지 않는다. Mesh area와 A_c를
동일시하지 않는다. UI view/field/언어 변경이 계산을 재실행하지 않게 분리한다.

## 16. 문서 라우팅

추가 연구는 기존 canonical model과 구분한다:

- CURRENT_MATERIAL_CORE_BRIDGE_V22.md / current_material_core_v22:
  최신 per-atom law의22채널(pair,rho,Q1,Q2even,Q3,Q2odd)을 무한 Bessel row에
  모두 유지한다. Eg/rank3 quartic/even saturation/density-screening 항이 있는
  후보를 base.surface만으로 예전17채널 core에 넘기지 않는다. 이는 연구 연결부이며
  기존 canonical energy/PDE는 변경하지 않는다. 각 atom의 환경 합 후 비선형 에너지,
  전체 chain-rule Hessian 및 영향받는 exterior site energy를 유지한다.
  Target-only Al99 core는 같은 FCC/벡터/경계/단위 비교용이며 LJ/Bessel 대체물이 아니다.
  Force convergence만으로 stable이라고 하지 말고 Morse index/eigenpair/FD를 확인한다.
  Reconstruction saddle, full lattice glide, finite-loop activation, 실제 항복은 다르다.
  직선 row-repeat당 energy에 임의 선 길이/Ac를 곱해 activation energy로 쓰지 않는다.
  Static ±5/20/50MPa load-return은 kinetic hold/잔류시편변형 검증이 아니다.
  Projected x-phase triangle crossing 위치는 연속 core 위치/partial 분리와 다르다.
  Source frozen-core force의 fixed-shape exact7/null3 lower bound는 해당 shape만의
  호환성 진단이다. 전체 family 불가능성을 증명하지 않는다. Core-only 계수 fit이
  다른 계면오차를 키우면 채택하지 않는다. Local joint shape trial도 실제 optimizer
  endpoint/제외검증/재이완과 구별한다. Parameter/default/material/kinetic gate 불변.
  Source 좌표 multistart는 own potential/exterior를 바꾸지 않는 초기값 검사다.
  Budget/급수실패 checkpoint를 certified minimum으로 승격하지 않는다; 명시적
  recovery와 새 force/Morse 검증이 필요하다. Coefficient energy분해는 고정gauge에서만
  해석하고, 임의 보간경로의 비용을 MEP/activation barrier로 부르지 않는다.
  새 Hurwitz-zeta majorant는 infinite transverse **LJ pair-force tail만** 상계한다.
  비선형 environment tail/free-domain/Al material 인증과 구분한다. Absolute normal
  prestress tail과 shear-response/국소force convergence를 혼동하지 않는다.
  Full atomistic core에는 elastic far-field energy가 이미 있으므로 v21 outer-line
  energy에 그대로 더하지 않는다. Same-energy matching/subtraction이 먼저다.
  넓힌 shape 탐색의 force RMS 개선을 material 채택으로 혼동하지 않는다. 잘못된
  opening/registry curvature 부호, outer optimizer 예산종료, active bounds를 공개한다.
  기존 x^z screening의 고정-shape scan은 새 에너지법칙이나 전체family 적합이 아니다.
  N-atom superperiod B=Nb와 basis regrouping은 문서의 후속 유도이며 아직 finite-kink
  구현이 아니다. Coincident-row self항을 positive-radius Bessel에 넣지 않는다.

- SPECIMEN_YIELD_AND_SLIP_BUDGET_V21.md / specimen_yield_bridge_v21:
  실제 signed swept area에서 beta_slip=Σ A b⊗n/V와 symmetric strain을
  계산한다. Uniform stress의 V sigma:deps와 Σ(b·sigma·n)dA가 일치해야 한다.
  Initial defect slip은 새 transfer에서 빼고, 회복하는 휨의 gross 면적을
  영구소성으로 더하지 않는다. Saved area 사이 숨은 역전은 시간해상도 문제다.
  동일 LJ/Bessel bulk 탄성의 finite-source MPa fold는 시편 항복이 아니다.
  rho_line=NL/V와 eta=NL³/V의 strain budget에서 필요한 밀도는 수요 진단이지
  문헌 보정값이 아니다. 이를 실제 항복에 맞추는 adjustable multiplier로 쓰지 않는다.
  Pigato2026 XML의 reportedYS12개는 source/온도/초기조직을 분리하며,
  수치 offset이 명시되지 않았으므로 Rp0.2로 바꾸지 않는다. 전위 방출/상호작용,
  실측 source population/코어/kinetic/production gate는 미완료다.

- TANGENT_CONSTRAINED_CALIBRATION_V20.md / tangent_calibration_v20:
  새 energy항 없이 기존 shape/10개 coefficient에서 5개 exact bulk와
  source pristine Haa/Hxx를 동시에 등식으로 시험한다. MPa/Å/eV 변환과
  LJ/Bessel은 불변. Imposed tangent CONTROL과 original source target를
  분리하며 common102 loss와 original104 loss를 모두 보고한다.
  Fixed-shape LP equality/sign feasibility와 spectral QP 검증은 다른 주장이다.
  v=0 LJ-attraction closure를 finite positive LJ 물성으로 채택하지 않는다.
  Shape optimizer가 zero-pair에서 멈췄는데 과거 positive trial을 골랐다면
  그 trial을 optimizer-converged라고 하지 않는다. 48개 새 heldout은
  fit/selection에서 제외. 현재 결과·실행/검증 상태는 최신 인계를 읽는다.

- COORDINATION_SCREENING_V19.md / coordination_screening_v19:
  같은 무한 scalar 환경 x와 rank1 벡터 Q를 site마다 합한 뒤
  D1||Q||² g(x)를 적용하는 별도 연구 가설이다. g=1/(1-z+zx)와 g=x^p는
  서로 대체하는 한-parameter 비교이며 동시에 곱하지 않는다. z/p=0은
  기존 v17을 복원한다. Affine inversion의 Q=0 때문에 bulk C는 불변이나
  실제 cubic dilation의 finite-q는 g(x_bulk) H1이며 density tail을 전파한다.
  같은 reference density/gauge를 strain 중 재보정하지 않는다. Site 합을
  global density로 나누지 않는다. 독립 direct energy에도 같은 함수를 쓴다.
  초기 fixed-range profile은 추가항 개선이 작았고 재료 채택 안 됐다.
  Joint radial 재보정/독립 검증 완료:260실제 profiles,297정적 tensor상태.
  새 후보 Haa는41.3%높고50MPa개구변위는29.1%작다. 추가 shape는 채택하지
  않는다. Old control의 max-nfev종료, 새 p=-1/k_even=12 bound를 숨기지 않는다.
  GPa C11/C12/C44 exporter와 sensitivity에 변환 전 mode단위를 쓰지 않는다.
  수치 rank10/11과 condition229/888은 inequality cone/통계CI가 아닌
  exact-bulk equality-tangent 진단이다. 세부 수치/한계는 최신 인계를 본다.
  과거 검사한 heldout은 retrospective로 구분한다. 낮은 loss를 실제 항복,
  physical a/s mobility나 Hz, production gate 통과로 읽지 않는다.

- LINE_KINETICS_AND_FATIGUE_VALIDATION.md / strength_fatigue_kinetics_v18:
  실제23°C Gorman 전위속도–전단10개 isolated marker에서7fit/3heldout을
  분리했다. mu_tau=1.1627e-5 m/(Pa s), heldout RMSE3.50m/s(평균의23.6%).
  원본 이미지/축/CGS변환/선택편향을 저장한다. 전위선 좌표의 제한된 보정이지
  a/s 이동도나 production clock이 아니다. B_line=b/mu_tau는 b의 조건을
  명시하고, B_line L²/(pi² T_line)은 hypothetical pin과0K 탄성 비교값을 쓴다.
  높은 rolloff주파수는 상수drag 수치모델 검증이며 실험 broadband 보정 아님.
  Pinned-line 작은 휨은 같은 LJ/Bessel longwave Schur stiffness에서 유도;
  analyticseries/finitegrid/time/hold/dissipation을 실제 검사했다. 회복하는 휨의
  루프를 피로/영구소성으로 부르지 않는다. 미지 L/core로 실제yield에 맞추지 않는다.
  Deschanel2017의 Al fatigue Delta50/62MPa는 full range(진폭25/31MPa)다.
  PSB/균열관찰/AE multiplet/final fracture/local absorption을 구분하고 Nf를
  M이나S-Nloss로 쓰지 않는다. 실제yield/Al fatigue/a-s초·Hz는 계속미검증이다.

- NORMAL_ENVIRONMENT_RESPONSE_V17.md / normal_response_v17:
  v16의 normal force/Haa 오차를 실제 coefficient control로 분해한 후 rank1
  microscopic decay k1만 독립화하는 한-parameter nested 연구 가설이다.
  k1=k_odd는 이전 full jet를 정확히 복원한다. Rank3/radial Eg 안의 k_odd를
  k1로 바꾸지 않는다. Per-atom vector 환경 합 후 norm, LJ/무한 Poisson–Bessel
  구조와 normalization은 보존한다. Affine inversion에서 Q1=0이므로 bulk
  C는 불변이지만 finite-q O(q^4) 및 interface response는 바뀐다. Full shape
  5개를 bulk/tail/direct/MD 검증에 전달하며 [:3]으로 k1를 유실하지 않는다.
  원래 양의 pair를 고정한 실험은 물리 관측값이 아닌 inherited CONTROL이다.
  v16 validation은 이제 development다. 고정-shape control의 tradeoff는 전체
  nonlinear 가족 불가능성 증명이 아니므로 같은 자료의 shape 대조군을 실행한다.
  새 heldout/analytic derivative/direct/finite-q/dilation/normal-shear-mixed
  응답을 재검사하기 전 material이나 core/PDE에 채택하지 않는다. 실제 완료/
  예산종료/물리적 실패는 CURRENT_WORK_HANDOFF와 원시 결과를 확인한다.
  완료된 v17은365profile/4nonlinear budget-stops와3후보 독립검증에도
  perfectHaa오차~41%가 남아 **추가 k1를 물리보정으로 채택하지 않았다**.
  Source1.62h곡률차이는 sourcecutoff/미분오차로 설명되지 않는다. 기존
  48heldout도 이제 검사한 자료다. Global가족불가능성/실제yield/Hz를 주장하지 않는다.

- EVEN_ENVIRONMENT_CALIBRATION_V16.md / even_environment_v16:
  기존 rank2/Eg per-site I=Q:Q에 한 개 공통 무차원 alpha의
  `D I/(1+alpha I/N^2)`를 적용한 별도 정적 연구 후보. 모든 환경을 합한 뒤
  site별 비선형 함수를 적용한다. 고정 cubic shear gauge이며 alpha=0은 기존 모델.
  Cubic Q=0의 모든 harmonic finite-q 및 cubic dilation은 불변이나 noncubic은 아니다.
  큰 alpha의 비균일 극한을 near-equilibrium 로그+균일 bracket으로 검사한다.
  기존 holdout은 development로 재분류했고 새36개 heldout을 loss에 넣지 않았다.
  7runs/420profiles의 에너지 장벽 개선은 실제 결과지만 force/Haa와300K covariance
  오차가 남아 Al 계면/PDE/kinetics/UI 채택은 안 됐다. Fixed-pair는 새 측정값 아닌
  CONTROL; zero-LJ closure는 재료가 아니다. 실제 QP의 작은 음의 물리계수는 clipping
  아닌 정확 활성경계 재해+원래 KKT 검사로 수정했다. 정적 MPa 시나리오는 실제 항복,
  dynamic hold, 잔류소성 증거가 아니다. 다음은 잔여 normal-response 원인 검증이다.

- LOADING_CONNECTIONS_AND_CALIBRATION_V15.md / interface_development_v15 /
  public_aluminum_validation_v15: 실제 production은 여전히 scalar축응력 P(a,s).
  연구 tensorprojection/vectorstatic/reflectingSG와 혼동하지 않는다. 실제Al99
  MD의576atoms/plane,12periodic class 정규화를 확인하되, 개별mode오차와
  thermostat/통계오차를 무시하지 않는다. 공개 line drag, PN Peierls, 실험
  apparentactivation derivative는 서로 다른 척도이며 a/s mobility/A_c가 아니다.
  ZERO_FREQUENCY_MOBILITY_AUDIT.md의 Gamma0=kBT C0^-1 (int Cdt) C0^-1는
  동일coordinate/PMF와 수렴한 적분을 요구한다. 음의 짧은시간 상관만으로
  모든 저주파 reduction을 부정하지도, 양수인 적분창만 골라 M을 만들지도 않는다.
  QP near-active/Gram수치결함은 원래 certificate를 유지하여 수정했다.
  낮은fitloss라도 독립 finite-q/작은strain의 불안정·계면오차가 있으면 탈락이다.
  FINITE_SOURCE_BARRIER_DERIVATION.md의 유한장벽과 Morse index는 지정된
  anisotropic outerline 모델의 결과이지 검증된 atomisticcore/피로확률이 아니다.

- TAIL_CONTROLLED_AL_CALIBRATION.md / tail_calibration_v14: 실제 FCC Voronoi
  tail bound를 coefficient-linear finite-q PSD 제약에 포함한 정적 재보정이다.
  기존 두 radial quadrupole을 결합한 Eg 채널은 signed cubic response로
  유도하고 amplitude gauge를 고정한다. 같은 per-atom 환경을 모두 합친 후
  norm을 취한다. Quartic rank3는 per-site ||Q3||^4이며 global norm 제곱 아님.
  별도 cubic-density/rational-shape ablation을 자동으로 합치지 않는다.
  **C11/C12/C44와 cohesion의 벌크 보정은 완료**됐지만 full-vector 계면 강성,
  direct110/개구 곡선 검증은 실패했다. 이를 전부 미보정 또는 완전 Al 채택으로
  뭉뚱그리지 않는다. 155q/radius12,16/tail/연속검색 양성은 전 q 안정성 증명 아님.
  Radial bound 확장에서 training loss가 내려가도 held-out 오차가 커진다.
  실패한 checkpoint와 완료된 calibration.json을 구별한다. Static MPa 이완/
  unload는 kinetics/실측 항복을 보정하지 않는다. 새 후보를 이전 screw core,
  production PDE/UI에 몰래 대입하지 않는다. Physical seconds/Hz 상태는 불변.

- STABLE_CORE_AND_MATERIAL_DIAGNOSTICS.md / stable_core_v13: 안정 source의
  hash/L0/논리적 row/경계를 묶어 실제 domain continuation을 한다. 초기 변위
  전달은 수렴 인증이 아니다. 선택적 Newton-CG도 같은 analytic Hessian을
  사용하며 force와 Morse를 독립 확인한다. Optimizer warning과 실제 잔차를
  함께 보존한다. Smooth radial window는 고정된 수치 partition이며 정확한
  logarithmic subtraction 상수가 필요하다. 이를 material core radius로
  피팅하지 않는다. Domain/ring/window 차이는 무한영역 오차와 별도 보고한다.
  같은 에너지에서 coefficient-linear H(q)의 negative mode는 necessary
  stability halfspace를 준다. 이 제약은 물리적 허용성 검사이지 새 yield law가
  아니다. Fixed-range convex profile와 실제 radial 재최적화를 구분한다.
  한 q/polarization의 조건은 전 Brillouin-zone 안정성 인증이 아니며,
  cutoff 한 곳의 minH≈0도 양의 안정성 margin으로 인정하지 않는다.
  더 작은 fit loss가 finite-q 불안정성/held-out 악화에서 온 것인지 확인한다.
  정적 0→50→0 MPa screw 비교를 거시 0.2% 항복이나 kinetic hold로 부르지 않는다.

- RANGE_AND_CORE_REPAIR.md / range_core_v12: STF rank2 radial range와 ranks1/3
  range 분리는 별도 analytic 연구 가설이다. Common-range limit을 보존하며,
  이를 지원하지 않는 scalar FFT row로 몰래 전달하지 않는다. 새 loss 개선이
  Al 채택을 뜻하지 않는다. Exact C11/C12/C44 후보의 finite-q 불안정도 확인했다.
  Isolated straight screw는 고정 anisotropic far-field 경계와 모든 affected
  site의 per-atom F를 포함한다. Periodic 반대전위쌍 소멸과 다른 문제이며,
  fixed boundary의 image force와 free disk/ring 수렴을 별도 검사한다.
  Centered core는 force≈0이어도 negative-Hessian saddle일 수 있다. 실제
  에너지 방향차분과 양쪽 negative-mode 이완을 확인한다. 불안정 초기점에서
  하중 후 남은 변위는 zero-load control과 비교하기 전 소성 증거가 아니다.
  Partial raw site-energy 합의 Taylor-linear reference work는 탄성 quadratic
  annulus energy와 구별한다. 이 항의 전체 합/변분이 불변임을 확인한 뒤
  raw/linear/remainder를 함께 보고한다. Core finite part의 reference와
  outer logarithm을 일치시키고, straight core를 finite activation energy나
  모든 character의 curved source/yield로 승격하지 않는다.

- YIELD_STRENGTH_BRIDGE.md / yield_bridge_v11: 실제 항복의 검증 기준과
  finite-source leading-log reference를 분리한다. Krebs2017 Fig2d의 0.002
  plastic-shear CRSS/G 7개를 원본 픽셀/공식 SI와 연결했다. G의 수치 미확인으로
  MPa는 null; 이전 large-strain flow를 항복으로 바꾸지 않는다. Wire D와
  source pin spacing L은 다르다. 문헌 L=D/3, D/2는 가정이지 실측값이 아니다.
  C11/C12/C44와 hydro/normal/shear mode는 정확한 선형 변환이다. 기존 metric을
  좌표만 바꾸려면 공분산의 off-diagonal도 유지해야 한다. 독립 5%-C metric은
  새 discrepancy 가정이지 단위 버그 수정이 아니다. material_metric_refined가
  최종이며 초기 SLSQP material_metric은 보존된 superseded 기록이다.
  5개 bulk를 exact하게 맞춘 최종 연구 후보도 fault/saddle 및 held-out 실패로
  채택하지 않았다. Fixed-decay SVD를 radial 포함 식별성/신뢰구간으로 부르지 않는다.
  같은 전위의 K(q)=|q|K0에서 k_line=b^T ReK0 b/(2pi)를 유도한다.
  양단 고정선의 outer energy gamma=k_line ln(R/r_core), line stiffness는
  gamma+gamma''다. 이를 써서 실제 정적 graph를 풀고 독립 parametric branch와
  격자 수렴을 비교했다. L=.5–10um, r_core/b=.5–2는 명시적 가정이며 보정값 아님.
  MPa bow-out threshold는 유한 atomistic core/source 활성화나 0.2% 항복 예측이
  아니다. Single-arm 실험 모델과 double-ended reference도 구별한다.
  Core/finite nonlocal energy와 실제 source geometry가 빠져 있다. 임의 길이를
  Arrhenius energy로 곱하지 않는다. b*swept_area/specimen_volume은 기하학적
  변형률 관계이지 한 bow-out이 proof strain에 도달했다는 뜻이 아니다.
  기존 energy/PDE/UI/physical-time/A_c는 바꾸지 않았고 gate는 미통과다.

- MATERIAL_TO_SPECIMEN_STRENGTH.md / material_strength_v10:
  같은 기존 LJ/Bessel+STF 가족의 exact coefficient basis로 0K bulk와 vector
  interface를 공동 fit했다. Initial Powell은 예산소진이며 global optimum 아님.
  Sample 사이 음의 개구력은 W_aa 극값 기반 coefficient constraint exchange로
  수정했지만 탄성/heldout shape 불일치로 새 후보도 채택하지 않았다.
  SOURCE 자체도 -785.7MPa opening lobe가 energy FD와 일치하므로 양의 견인력
  조건은 이번 shape prior이지 보편 법칙이 아니다. Prior 제거 QP에서도 탄성
  오차가 남았다. 음의 견인력만으로 구현 bug/비물리라고 단정하지 않는다.
  Source 빈 neighbor sum은 정확히 zero jet로 처리한다. 이는 source-only
  빈 einsum NaN 재현에 대한 정확한 identity 처리이지 LJ cutoff/force clipping이 아니다.
  실온 Krebs2017 Al wire 흐름점은 validation-only이며
  v10에서 0.2%CRSS/finite-source geometry는 미지였다(후속 v11 참조).
  Flow stress, first slip, uniform
  ideal fold와 axial/resolved shear를 등치하지 않는다. 직경을 source length나
  A_c로 바꾸지 않는다. Static MPa probe는 실제 plasticity/physical-time hold가 아니다.

- VECTOR_REGISTRY_AND_STRENGTH_AUDIT.md / vector_registry_v9: 기존 후보 그대로
  q=(a,ux,uy)의 full registry 계면을 계산한다. 각 site 환경 합 뒤 F 및 STF
  norm을 적용하며, 원래 scalar PDE/production selector는 바꾸지 않는다.
  direct110 중간의 W_x=0만으로 saddle이라 하지 않는다. W_a/W_y가 남으면
  구속 반력이 필요한 경로다. 전체 Hessian index와 saddle 양방향 descent로
  연결된 minimum을 검증한다. 고정-x에서 a,y를 이완한 branch는 local constrained
  branch이지 global MEP 증명이 아니다. Schur curvature의 첫 root와 H_zz>0으로
  uniform-interface fold를 확인한다. 그 GPa 이상강도를 실험 항복으로 바꾸지 않는다.
  Al 기준 전위도 동일한 rigid-half geometry에서는 GPa 이상강도를 준다.
  후보의 forward saddle 일치만으로 material 통과 불가: reverse fault barrier와
  independent elastic tensor 오차가 남는다. Static unload는 dynamic hold가 아니다.
  SOURCE 전위는 target/검증 전용이며, cutoff 밖의 평탄한 zero-force 상태를
  intact minimum으로 승인하지 않는다. Nonfinite jet은 실패로 처리한다.
  Physical kinetics, finite source, production 및 UI gate는 여전히 미통과다.

- VECTOR_FCC_CORE_DERIVATION.md / vector_core_v8: 동일 LJ/Bessel 후보의 실제
  3성분 무한 원자열 코어 이완이다. 실제 y,z 반경이 변하므로 고정계수 FFT를
  재사용하지 않고 m=0 포함 원자열 급수/analytic vector Hessian을 계산한다.
  각 site의 모든 환경 변화를 합산한 뒤 F(rho_i), D_r||Q_ri||²를 적용한다.
  세 local displacement는 자유지만 transverse periodic cell shape는 고정이다.
  무한 직선 원자열이지 유한 전위 loop/source나 3D 동역학이 아니다.
  v7의 큰 omitted force를 줄였지만 가까운 반대 전위쌍은 소멸할 수 있고,
  작은 주기 cell에는 cell-spanning registry fault가 남을 수 있다.
  x winding=0만으로 full vector defect 소멸을 선언하지 않는다. 실제 벡터
  minimum이 x=b/2 scalar partition에 있을 수 있으므로 기존 scalar n을
  벡터 소성 지표로 쓰지 않는다. layer_registry의 두 slip 성분, 0/tau/2tau
  lattice-equivalent 거리, opening, row dispersion을 함께 본다.
  횡변위는 LJ zero-mode algebraic tail을 살린다. v7의 작은 tail error를
  복사하지 않는다. 에너지/force 및 전체 Hessian symbol의 ring 수렴을 본다.
  최소 고유값만 같아도 polarization이 교차해 전체 matrix 오차가 남을 수 있다.
  continuation은 source state/parameter binding을 기록하고 중단된 계산을
  완료라고 보고하지 않는다. static unload는 물리 시간 hold가 아니다.
  Al material fit, finite-source geometry, kinetics, production/UI gate는 미완료다.

- KINETICS_LOADING_AND_STRESS_AUDIT.md / kinetics_loading_audit:
  실제 covariance CSV에서 B=-log(Ct C0^-1)/t 및 M=B H^-1를 검증한다.
  file binding은 source의 물리적 진실성 인증이 아니다. Physical mode는
  보정의 T와 mobility ratio를 같은 PDE에 적용; model mode 기본값 불변.
  actual Al mobility/t0는 여전히 없다. Dislocation-line drag를 cell a/s에
  대입 금지. UI direction placeholder는 미연결로 명시한다.
  3x3 tensor traction projection/연구 혼합하중은 3D PDE나 spatial solver가 아니다.
  고정-a first traction peak와 normal-relaxed Schur spinodal, 나중 registry
  peak를 구별한다. 다중 curvature root를 포함하는 큰 bracket에서 brentq
  한 번으로 첫 물리 불안정점을 골랐다고 주장하지 않는다.
  Native Tcl 오류를 headless skip로 숨기지 말고 실제 GUI 결과를 확인한다.


- NONLINEAR_FCC_SCREW_DERIVATION.md / nonlinear_screw_v7: 같은 LJ/Bessel 후보의
  per-site nonlinear scalar anti-plane row energy. 전체 환경 합 후 F 적용,
  FFT는 retained row coefficient의 정확한 convolution이다. x-force 평형과
  full vector core 평형을 구별한다. 실제 frozen y/z force~2.28eV/L0가
  tail error~1.7e-10보다 커서 현재 core는 full mechanical equilibrium 아님.
  Transverse force는 x-slip에서 소거된 reciprocal m=0을 반드시 포함한다.
  고정 gamma=0은 전위가 있는 cell의 zero stress가 아니다. 실제 cell volume의
  -tau*V*gamma로 stress-control하고 initial defect content를 새 slip과 분리한다.
  24/32/48 원자열 domain, +/-50MPa 정적44개에서 새 registry 변화0.
  초기 winding pair가 있다는 이유로 생성/잔류소성/Al 강도 검증이라 하지 않는다.
  Static unload는 physical-time hold가 아니다. Energy/b-repeat는 J/m이며
  finite activation energy가 아니다. Full transverse/vector relaxation,
  partial core, finite-source geometry, material/kinetics 및 production gate 미완료.
  y,z를 실제 움직이면 radius가 site-dependent여서 고정계수 FFT를 그대로
  canonical vector energy로 재사용하지 않는다. Bessel/G=0/tails를 보존한다.

- DISCRETE_FCC_SCREW_DERIVATION.md / discrete_screw_v6: 같은 LJ/Bessel 전위의
  무한 원자열을 먼저 합산하고 transverse row/layer를 이산적으로 남기는
  anti-plane harmonic reference. q_x=0, line=e1에만 해당한다.
  Scalar F''가 사라지는 것은 각 완전한 row의 rho_x=0인 이 subspace에서만이다.
  Angular moment는 모든 row를 합친 후 제곱한다. per-atom EAM을 보존한다.
  논리적 인접 row jump와 spectral same-y jump는 finite-q에서 다른 제약이다.
  `K_jump=(B H^+ B^dagger)^-1`에는 local K0가 이미 들어 있으므로
  기존 gamma Hessian을 더해 이중 계산하지 않는다. K0를 뺀 후 rigid gamma를
  붙이는 것도 nonlinear relaxation convention 검증 없이 채택하지 않는다.
  같은 continuum reference의 전체 profile Euler 방정식은 별도로 풀었다.
  Mean-content 제약의 반력은 holding traction이고 내부 force residual과 다르다.
  고정 mean content와 실제 core separation은 다르다. Domain 비교는 실제
  crossing 간격을 맞추고 해야 한다. Mesh pinning을 atomistic Peierls로 부르지 않는다.
  Atomistic nonlinear/vector/normal-relaxed core, finite-loop activation 및
  독립 Al material fit/kinetics는 아직 미완료. Production/UI gate 미통과.

- NONLOCAL_INTERFACE_ELASTICITY.md / nonlocal_v5: 같은 LJ/Bessel 후보의 bulk
  C_ijkl에서 두 반무한 결정의 static Schur kernel을 유도한 long-wave reference.
  `K=|q|(Z_plus^-1+Z_minus^-1)^-1`, K(0)=0이며 local gamma Hessian을 중복하지 않는다.
  Straight screw 에너지는 J/m 또는 eV/L0(line)이며 유한 activation energy가 아니다.
  임의 선 길이/A_c를 곱해 Arrhenius 확률로 만들지 않는다. 현재 arctangent trial의
  폭~0.28b, full Euler residual~1.2GPa: 실제 atomistic core를 푼 것이 아니다.
  4–50MPa에서의 기존 전위쌍 구동력은 생성/속도/피로/잔류소성 검증과 구별한다.
  후보와 source 탄성값은 서로 다른 결과를 주므로 교체해 성공을 만들지 않는다.
  다음은 discrete half-crystal kernel, core 및 vector registry/finite-loop 검증이다.

- LOW_STRESS_CYCLIC_AUDIT.md / low_stress_v4: 10–50 MPa 축응력 및 4 MPa 전단
  실제 주기 SG 가설 검사를 무하중 대조군/격자/dt/hold와 비교한다. Production
  등록이 아닌 별도 reflecting-domain 수치 실험이며 opening_probability=null이다.
  Uniform infinite-interface per-cell energy를 곧바로 국소 activation energy로
  간주하는 정규화는 미검증이다. 작은 수렴한 순이동도 Al 피로 검증과 구별한다.
- MATCHED_INTERFACE_CALIBRATION.md / matched_v3: 같은 rigid FCC 조건의 NIST
  source를 target으로만 재계산하고 joint fit/응력 시나리오를 실제 실행한다.
- reference_eam_targets.py의 tabulated EAM은 target 생성/비교 전용이다.
  LJ/Bessel production potential이나 PDE energy selector로 사용하지 않는다.
- ANALYTIC_ANGULAR_ENVIRONMENT.md: per-atom STF 1/2/3차 moment invariant의
  별도 연구 가설. 홀수차는 affine centrosymmetric bulk에서 0이지만,
  짝수차 Q2는 cubic bulk에서만 0이다. noncubic Q_bulk 누락 금지.
  수치 검증/적합 개선만으로 calibrated Al/standard MEAM이라
  부르지 않으며 생산 모델에 자동 추가하지 않는다.
- STATISTICAL_CORRELATION_AREA.md: 공간 covariance에서 얻는 variance-matching
  area가 void/survival probability의 독립영역 area와 자동으로 같지 않다.
  A_c의 값/공간독립성/kinetics는 여전히 보정되지 않았다.

응력 시나리오의 root 실패는 확인된 spinodal과 구별한다. 정적 unload는
동적 zero-stress hold나 잔류 소성 검증을 대신하지 않는다. 같은 상태/이완/단위의
source 비교를 우선하고, fit에 쓴 점을 held-out RMSE에 섞지 않는다.
개구 energy endpoint만 맞추고 중간 traction overshoot를 숨기지 않는다.
coarse W_a>=0 grid는 연속 단조성 증명이 아니다. 실제 W_aa=0 극값과
독립 bracket/tolerance refinement를 검사한다. Force clipping 금지.
Finite-q K(q)도 static Hessian이다. 질량/kinetic calibration 없이 Hz로 바꾸지 않는다.

작업별 관련 문서를 실제 최신 코드와 함께 읽는다:

- solver_v1/RESULT_FIELDS.md
- solver_v1/CONFIGURATIONAL_PLASTICITY.md
- solver_v1/PLASTICITY_AND_STRAIN_AUDIT.md
- solver_v1/PHYSICAL_TIME_AND_MOBILITY.md
- solver_v1/ANALYTIC_LJ_EAM_HYBRID.md
- solver_v1/ALUMINUM_ANALYTIC_CALIBRATION.md
- solver_v1/FCC111_FULL_STACK_DERIVATION.md
- solver_v1/ALUMINUM_FULL_FCC_CALIBRATION.md (independent-mode 감사; 물리 검증 미통과)
- solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md (static research; production 아님)
- solver_v1/SOLVER_VALIDATION_GATES.md
- solver_v1/data/ 및 해당 results/의 machine-readable provenance/results

경로가 없거나 stale이면 새 경로를 검색하고 그 사실을 밝힌다.
문서의 주장과 코드/검증이 충돌하면 숨기지 않고 최소 재현 검사를 추가한다.
Full-FCC 독립 bulk 감사의 근거는 audited_v2이며, 후속 matched-condition
계면 보정/응력 검사는 matched_v3에 별도로 저장한다. 상위 legacy CSV/그림은
superseded 역사적 기록이다. CURRENT_WORK_HANDOFF.md의 최신 상태를 우선한다.

## 17. 검증과 최종 보고

먼저 targeted tests, 이어서:

    py -3 -m pytest solver_v1 -q
    py -3 -m pytest app -q
    py -3 -m app.desktop_ui --smoke
    git diff --check

py launcher가 없으면 실제 확인한 Python 3 interpreter로 동등 명령을 실행한다.
실행하지 않은 테스트를 PASS라 하지 않는다.
Lost session output / interrupted run은 결과 미확인이다. 필요하면 로그와 함께 재실행한다.
Targeted, full solver, app, smoke 결과와 timing을 각각 보고한다.
새 이론 구현은 analytic identity, independent direct sum, derivative checks,
grid/time/tail refinement 및 기존 reference 불변성을 검증한다.

최종 보고에는:
- 시작/원격/최종 HEAD, worktree, changed files, commit/push 여부;
- 실제 에너지 모델과 parameter status;
- derivation, normalization, residuals, truncation/convergence;
- resolved/unresolved/absent, 모델의 실패와 한계;
- physical mobility/seconds/Hz availability;
- 미완료 항목과 다음 검증
을 명시한다.

**소성을 만들지 말 것. Al 보정 성공을 만들지 말 것.
물리 시간/Hz를 만들지 말 것. 부정적인 과학 결과를 숨기지 말 것.**
