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
  0.2%CRSS/finite-source geometry는 미지다. Flow stress, first slip, uniform
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
