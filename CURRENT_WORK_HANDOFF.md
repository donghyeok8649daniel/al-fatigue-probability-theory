# CURRENT_WORK_HANDOFF.md — 단계별 검증 후 재개하기

기록일: 2026-09-09. 영구 이론/작업 원칙은 AGENTS.md를 따른다.
이 문서는 현재 계산의 근거, 부정적 결과, 미완료 과제와 재현 순서를 보존한다.

## 1. 최신 사용자 지시와 현재 결정

사용자는 AGENTS.md/인계 문서를 읽고 이어서, **각 단계 구현 후 실제 실행과
물리적 타당성 확인을 먼저 하고 다음 단계로 이동**하라고 지시했다.
캘리브레이션/단위를 정확히 하고 연구 이론과 Git을 보존해야 한다.

솔버가 정말 준비되면 먼저 사용자에게 보고한다. **사용자가 확인한 뒤에만**
모델 파일 가져오기 / 기본 원통 시편 -> meshing -> pre -> solving/posting UI로
넘어간다. 향후 mesh 위에서 사용자가 선택한 물리량을 컬러맵으로 표시하되,
근거 없는 공간 분포나 무거운 UI를 만들지 않는다.

현재 결정:

    full-FCC Al 재료/계면 물리 검증: 미통과
    기존 reduced 확률 PDE: 그대로 유지
    새 active interface -> production PDE: 연결하지 않음
    CAD / cylinder / mesh UI 구현: 시작하지 않음
    물리 mobility / seconds / Hz: unavailable

수치 테스트 PASS와 물리 검증 PASS를 혼동하지 않는다.
자세한 단계별 gate: solver_v1/SOLVER_VALIDATION_GATES.md.

## 2. Git / 실제 작업 위치 / 사용자 파일

이번 감사의 시작 HEAD와 실제 fetch로 확인한 원격 기준:

    b819abea7e8514d04698e6ab0e391b940262d16c

실제 과학 작업은 git worktree list에서 basename
aft-pde-bessel-38969ad인 별도 worktree에 있다.
이번에 해당 worktree의 detached 상태에서 같은 HEAD의 로컬
probability-pde-solver-v1 tracking branch를 안전하게 만들었다.
현재 작업/커밋 대상은 **그 worktree의 해당 브랜치**다.

원래 사용자 workspace는 별도 detached HEAD:

    c43d8e096f2a329c7cdfad31e546c70fb36fe592

원래 workspace의 docs/, examples/, fem1d/, libraries/, output/, paper/,
research/, results/, simulations/, tests/, theory/, tools/, requirements 계열
미추적 사용자 파일을 건드리지 않았다. 삭제/초기화/일괄 덮어쓰기 금지.
발견용 AGENTS.md와 이 인계 문서만 양쪽에 같은 내용으로 갱신한다.

main 기준 ref는 80cacb4180dbfbcd36a2964270703bc6cf1653ec이며 수정하지 않는다.
최종 commit/push 상태는 현재 branch의 git log/status와 새 fetch로 확인한다.
문서 안에 자기 자신을 포함하는 새 commit SHA를 예측해서 기록하지 않는다.
개인 절대 경로나 local test XML은 커밋하지 않는다.

## 3. 인계 당시 문제를 실제 수정한 내용

### 3.1 이전 six-observation rank 5 주장은 잘못됐다

[111] normal alpha, equal-biaxial beta 관측량:

    X = C11 + 2 C12
    Haa = V(X + 4 C44)/3
    Hbb = 4V(X + C44)/3
    Hab = 2V(X - 2 C44)/3
    Hbb = 2 Haa + Hab

Perfect cubic configuration에서 lateral force = 2 normal force.
따라서 기존 stage 1은 독립 rank <= 2, 기존 stage 2는 rank <= 4.
C11-C12가 제약되지 않았으므로 five-parameter 식별성을 주장할 수 없었다.
작은 numerical fifth singular value는 유한차분/중복 관측량의 오차였다.

legacy 후보를 독립 전단으로 다시 보면:

    C11 ~= 95.2225 GPa  (target 114)
    C12 ~= 71.3890 GPa  (target 62)
    C44 ~= 32.0002 GPa  (target 32)

기존 잘못된 rank assertion을 이 항등식/독립성 회귀 테스트로 대체했다.
원래 calibration 결과/fit vector는 삭제하지 않고 SUPERSEDED로 표시했다.

### 3.2 빠진 independent homogeneous simple shear 추가

    F = I + gamma m tensor n
    m = [1,-1,0]/sqrt(2)
    n = [1,1,1]/sqrt(3)
    gamma = s/h111

이는 기존 homogeneous stack이 표현하는 affine shear이고 local GSF가 아니다.

    H_etaeta = 3V(C11+2C12)
    H_alphaalpha = V(C11+2C12+4C44)/3
    H_gammagamma = V(C11-C12+C44)/3

세 독립 탄성 모드로 C11,C12,C44를 복원한다.
이것은 homogeneous local stability만 확인하며 phonon/global phase stability
검증을 대신하지 않는다.

### 3.3 실제 deterministic fit 재실행

새 실행:

    python -m solver_v1.run_full_fcc_calibration_audit

저장된 숫자를 재생하는 스크립트가 아니다.
고정 density decay d에서:

    u = 4 epsilon sigma^12
    v = 4 epsilon sigma^6
    y = Q(d) [u,v,A,B]^T
    sigma = (u/v)^(1/6)
    epsilon = v^2/(4u)

NNLS coefficient solve + d=0.35..8의 49개 deterministic log-grid 및
bracketed minima refinement. 이전 탐색 box에서는 고정된 세 starts의
bounded least_squares도 실제로 수행했다.
Search bounds는 물리적 Al parameter confidence interval이 아니다.

Stages: equilibrium/cohesion -> normal/hydrostatic -> independent shear.
GSF/opening은 held-out contextual comparison이며 loss에 넣지 않았다.
Scale: force .10 eV/strain, cohesion .0672 eV/atom, 각 독립 탄성 curvature 5%.
실험 불확도가 아니라 명시한 model-discrepancy normalization이다.
B=0 square와 기존 linear extension을 비교했으며 새 C 항은 넣지 않았다.

### 3.4 단위와 실제 FCC 평형

    L0 = b_target = 4.05/sqrt(2) angstrom
    E0 = 1 eV = 1.602176634e-19 J
    h_target/L0 = sqrt(2/3)
    A_atomic,target = 7.1024908428 angstrom^2
    V_target = 4.05^3/4 angstrom^3

불완전한 fit에서 b를 고정하고 a만 root 찾으면 무압력 FCC가 아니다.
현재는 b=lambda, a=sqrt(2/3)*lambda로 **isotropic lattice equilibrium**을 푼다.
물리 potential, L0, rho_ref는 strain 동안 고정한다.
실제 면적 lambda^2 A_target, 실제 체적 lambda^3 V_target을 각각 사용한다.

    GSF[J/m^2] = W[eV/cell] E0/A_atomic
    traction[Pa] = W_a E0/(L0 A_atomic)
    modulus[Pa] = strain Hessian[eV/atom] E0/V_atom

A_c/시편 면적/mesh 면적은 여기에 들어가지 않는다.
F(x)=-A sqrt(x)+B(x-1)이므로 F(0)=-B.
Cohesion은 atomized reference를 포함한 -B-W_bulk이다.

## 4. 현재 수치 결과 — audited_v2가 기준

독립 bulk targets는 0 K Mishin EAM reference:
a_lat=4.05 angstrom, cohesion=3.36 eV/atom,
(C11,C12,C44)=(114,62,32) GPa. 실험적 exact Al이라고 주장하지 않는다.

Parameter order:
(epsilon_LJ[eV], sigma_LJ/L0, density_decay_L0, A[eV], B[eV]).
기존 dataclass의 over_b/decay_b 이름은 고정 target L0 convention이다.

Square profile:

    (2.8462556966e-5, 1.737841252, 1.567711191, 3.862759460, 0)
    sum r^2 = 10.92396
    actual a_lat = 4.049133896 angstrom
    cohesion = 3.36139918 eV/atom
    C11=91.9543, C12=60.0708, C44=49.5093 GPa

이 탐색에서 선언한 탄성 scale 내 적합 실패.
전체 square-root family의 global 불가능성 증명으로 과장하지 않는다.

Linear unrestricted coefficient profile:

    (18.057395467, .69915514077, 5.1764034360, 3.7025097558, 55.334763179)
    sum r^2 ~= 6.74e-13
    actual a_lat ~= 4.050000002 angstrom
    C11=113.999992, C12=61.999998, C44=32.000003 GPa
    pair energy ~= -54.99225 eV/atom
    embedding relative atomized ~= +51.63225 eV/atom
    cancellation ratio ~= 31.7335

독립 log-J singular values:
~6797.18, 1334.55, 133.733, 3.97779, 3.41103; condition ~1992.70.
새 관측량에서는 local rank 5가 독립 모드로 확인되지만 global uniqueness/
물리 parameter 식별성/신뢰구간을 뜻하지 않는다.

Linear historical-box compromise:

    (.84661409719, .80696569538, 1.36764964209, 10.1664763961, 12)
    B upper bound active
    sum r^2 = .9293793254
    actual a_lat = 4.0469346273 angstrom
    C11=110.45729, C12=64.20875, C44=37.31294 GPa

## 5. Active-interface의 부정적 물리 결과

한 cut의 cross-half pair difference와 depth별 per-atom density/embedding
difference를 사용한다. Infinite bulk constants는 analytic cancellation.
모든 neighbor density 합산 후 원자마다 F를 적용한다.

현재 rigid-half 결과 (J/m^2):

| candidate | direct110 USF | Shockley USF | ISF | separation work |
|---|---:|---:|---:|---:|
| linear unrestricted | 1.35450 | .092901 | .004990 | 12.75507 |
| historical box | .812221 | .108680 | -.001748 | -12.28939 |
| legacy underidentified | .448626 | .076844 | -.000241 | 1.14555 |

문헌 contextual references: direct .250, Shockley .224, ISF .164 (Lu DFT);
separation 1.74 (Mishin EAM), 2.12 (Wei DFT).

**Mapping caveat**: Lu는 relaxed atoms/volume, DFT a_lat=3.94 angstrom이다.
여기서는 EAM-scale lattice의 rigid halves다. 같은 relaxation/geometry인
것처럼 residual을 해석하지 않는다. 정확한 mapping 대상 확보가 다음 단계다.

그럼에도 negative separation work는 명백한 candidate failure:
locally positive bulk/interface Hessian이어도 분리된 surfaces가 벌크보다 유리하다.
이 후보는 채택 불가. Bulk-exact 후보의 큰 decohesion/GSF mismatch도
검증 성공으로 볼 수 없다. 모든 family가 불가능하다고 증명한 것은 아니다.

## 6. 국소 opening barrier와 수치 감사

- reciprocal/layer/neighborhood tolerances: 2e-8,2e-10,2e-12 비교;
- reciprocal tolerance 최대2e-14;
- direct radius/layers: (24,24),(48,24),(24,48),(48,48),(80,80);
- analytic gradient/Hessian과 finite differences 비교;
- hydrostatic strain step 및 log sensitivity step 비교;
- zero opening reference와 large-a separation asymptote 검증.

Representative point a=1.05h,s=.17L0에서 coarse/fine reciprocal discrepancy:
최대 energy ~4.06e-9 eV/cell, gradient ~9.68e-8 eV/L0,
normal Hessian ~1.30e-7 eV/L0^2.
최대 direct(80,80) discrepancy는 energy~4.92e-6 eV, normal gradient~1.21e-4.
이는 finite direct LJ power-tail이고 refinement로 감소한다.
모든 direct 값이 machine precision이라고 주장하지 않는다.

estimated_tail_absolute라는 역사적 field는 density 등을 합친 internal
convergence proxy다. eV 단위 rigorous bound가 아니다.
실제 energy/derivative refinement differences를 authoritative로 사용한다.

Fixed-s=0 normal opening은 첫 W_aa=0까지 intact branch를 추적한다.
넓은 구간을 한 번의 unimodal minimization으로 찾으면 historical-box의
매우 좁은 초기 positive traction peak를 놓치는 버그를 발견해 고쳤다.
그 후보의 peak force~.0540409 eV/L0 (~.426324 GPa).
25%-95% peak의 local barriers~8.87e-4 -> 1.52e-5 eV/cell.
이 작은 metastable barrier는 negative separation work를 정당화하지 않는다.

Unrestricted 후보에서는 일부 constrained normal minima의 W_ss<0.
따라서 저장된 normal-only barrier를 안정한 2D crack barrier나 coupled
spinodal로 사용하면 안 된다. CSV에 registry curvature를 명시했다.
완전한 coupled saddle/spinodal 및 동적 event ordering은 미완료/미승격이다.

## 7. 결과/코드 위치

현재 자료:

    results/aluminum_full_fcc_calibration/audited_v2/
        independent_targets.csv, staged_profile.csv, optimization_log.json
        parameter_sets.csv, fit_residuals.csv, sensitivity_matrix.csv
        identifiability_svd.csv, derivative_refinement.csv, readiness.json

    results/fcc111_active_interface/audited_v2/
        heldout_validation.csv, gsf_curve.csv, opening_curve.csv
        hessian.csv, convergence_summary.csv, opening_barriers.csv
        active_interface_energy_grid.csv, physical_gate_comparison.png

상위 directories의 기존 CSV/그림은 historical underidentified evidence.
상위 summary.json에 SUPERSEDED와 current_results를 표시했다.
기존 run_aluminum_full_fcc_calibration.py는 ARCHIVAL saved-vector replay이며
별도 legacy_replay 하위 폴더로만 출력한다.

새 주요 API:

- full_fcc_calibration_audit.py: independent modes, coefficient profile,
  bounded fits, Jacobian, isotropic equilibrium, actual-volume elasticity.
- run_full_fcc_calibration_audit.py: 실제 fit + validation 결과 생성.
- fcc111_active_interface.py: static interface geometry and analytic derivatives.
- test_full_fcc_calibration_audit.py: 독립 elasticity/unit/gauge/fit/실패 검증.

fcc111_full_energy.py의 기존 동작은 유지. find_equilibrium=False일 때만
연구 strain evaluator가 root search를 생략한다. True가 기존 default.

## 8. 검증 기록

실제 재실행:
- full audited calibration/GSF/cleavage/refinement runner: 297.762 s, 완료;
- narrow-barrier correction 후 barrier table 12행 별도 재계산 완료;
- audit-only targeted tests: 17 passed, 13.93 s;
- 최종 combined targeted: 32 passed, 47.71 s;
- 최종 solver_v1 full regression: 177 passed, 759.98 s (12분 39초);
- app: 27 passed, 2 skipped, 64.24 s;
- desktop --smoke: PASS, a0=.7713438268704838, kappa=86.29296488740997;
- git diff --check: PASS, staged new files 포함.

동시에 실행한 CPU 작업이 있어 timing은 해당 실행의 실제 wall time이다.
모든 최종 테스트의 exit code=0을 직접 확인했다. 추정 수치가 아니다.
커밋은 이 부정적 물리 감사/백업의 완료이며 Al solver 완성 선언이 아니다.
원격은 최종 검사 전 재-fetch하여 b819abea 기준과 동일함을 확인했다.

Test XML은 local-only evidence로 ignore한다 (machine-specific path/hostname).
최종 검증 기록은 audited_v2/validation_record.json에도 저장한다.
숫자 추정으로 PASS라 하지 않는다.
이 환경에서 py -3 launcher가 없어 확인한 Python 3.13 interpreter로 동등 실행.

## 9. 남은 일 / 다음 재개 순서

1. 실제 worktree, 현재 branch/HEAD, fetch 원격 확인. 기존 사용자 작업 보존.
2. 이 문서와 두 derivation/calibration 문서 및 readiness.json 읽기.
3. 같은 lattice/relaxation/orientation의 atomistic interface targets를 마련하거나
   같은 relaxation convention을 수학적으로 구현해 GSF/cleavage를 비교.
4. 기존 family의 대상 관측량을 더 제약하고 degeneracy/양의 bulk·surface
   안정성을 함께 감사. 소수 local fit 실패를 global theorem으로 만들지 않는다.
5. 필요한 경우에만 최소 analytic embedding extension을 검토.
   arbitrary LJ/mobility/chi tuning으로 결과를 만들지 않는다.
6. 물리적으로 채택 가능한 interface 후보가 생긴 뒤 coupled saddles/GSF topology
   검증. 그 뒤에야 probability generator/drift/absorption 의미를 유도·검증.
7. 독립 MD collective-coordinate kinetic data 전에는 seconds/Hz 불가.
8. 공간 mechanical field와 local probability 연결을 유도·검증해야 실제 mesh
   parameter colormap이 의미 있다. 하나의 global P를 공간 해석인 양 칠하지 않는다.
9. 솔버 준비 증거를 사용자에게 보고하고 확인 받은 뒤에만 UI 작업.

현재 결론은 **수학/단위/수치 검증을 개선했지만 full-FCC Al 계면 솔버는
아직 물리 검증 미통과**다. 이 단계에서 UI로 넘어가지 않은 것이 의도된 gate다.
기존 TwoRowLJ/reduced hybrid/static Al files/PDE/strain bookkeeping/A_c/time
framework는 유지한다. 물리 mobility나 Hz를 만들지 않는다.
