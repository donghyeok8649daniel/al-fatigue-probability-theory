# CURRENT_WORK_HANDOFF.md — 단계별 검증 후 재개하기

## 최신 상태: matched_v3 연구 감사 (2026-09-09)

이 절이 아래 archived audited_v2보다 최신이다. 마지막 원격 checkpoint는
fee1da7a5e856bd9120a12f0c40e6afb79ec8c6e이며 시작 fetch로 일치를 확인했다.
이번 새 연구의 최종 검증은 완료했고 이 문서를 포함하는 후속 커밋에 보존한다.
commit/push의 최종 상태는 git log/status 및 새 fetch로 확인한다.
자기 자신의 commit SHA를 문서에 예측하지 않는다.

사용자의 최신 뜻:
- 물리적으로 말이 되도록 조건 안에서 실제 수정/계산할 것;
- 응력을 바꾸며 여러 시나리오를 실제로 실행하고 다음 단계 채택 여부를 판단할 것;
- LJ 급수합/Bessel 수식 틀을 보존하고 이론의 논문 가치도 중요하게 다룰 것;
- 가능한 경우 통계적 특성 상관 면적도 검토할 것;
- 검증된 solver 준비 상태를 먼저 보고하고, 사용자 확인 뒤 UI workflow로 이동할 것.

현재 결정:
- 기존 TwoRowLJ / reduced hybrid / production probability PDE / UI default: 변경 없음.
- 새 full-FCC angular/interface surface: STATIC RESEARCH ONLY.
- 정적 평형과 여러 안정성 검사는 개선/통과했지만, Al 독립 탄성/곡선 적합은 미통과.
- quantitative Al potential 또는 fatigue solver 완성 선언 없음.
- physical mobility/seconds/Hz: 모두 unavailable; A_c도 외부 미보정.
- UI CAD/cylinder/mesh redesign은 아직 시작하지 않음. Gate를 숨기지 않는다.

### A. 실제 작업 위치와 보호

과학 브랜치는 basename aft-pde-bessel-38969ad인 별도 worktree의
probability-pde-solver-v1이다. original workspace는 detached c43d8e0이며
많은 사용자 미추적 파일이 있다. 원본 코드를 reset/checkout/덮어쓰기하지 않는다.
발견용 AGENTS.md와 이 인계 문서만 소유/이전 버전을 확인한 뒤 동기화한다.
main local ref 80cacb4180dbfbcd36a2964270703bc6cf1653ec는 건드리지 않았다.
py launcher가 없으므로 실제 Python3.13 interpreter로 동등 명령을 사용했다.
개인 절대 경로나 test XML을 커밋하지 않는다. XML은 *.local.xml로 ignore한다.

### B. 같은 조건의 원자 기준을 먼저 마련했다

reference_eam_targets.py는 NIST Al99.eam.alloy를 target 생성에만 사용한다.
절대로 LJ pair/production potential/PDE selector로 등록하지 않는다.
캐시: ignored .cache/al-reference/Al99.eam.alloy
SHA256: 60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284
DOI:10.1103/PhysRevB.59.3393; NIST source URL은 코드/문서에 기록했다.
다운로드는 명시적 --download, checksum 확인, 다른 기존 파일 덮어쓰기 거부.
자료가 없으면 원자 수치를 만들어 넣지 말 것.

동일 rigid FCC half-crystal, 0 K, a_lat=4.05 angstrom:
cohesion=3.359999988239 eV/atom;
W_sep=1.741285308649 J/m2;
direct110 sampled maximum~.6030383534 J/m2;
Shockley maximum~.1897111890 J/m2;
ISF=.156630432472 J/m2=.06943467959308 eV/interface cell.
Source local normal tangent=126.44968004 GPa; registry tangent28.25750958 GPa.
원자모델 reference이지 실제 Al experimental truth가 아니다.
예전 relaxed Lu DFT .250/.224/.164, a_lat3.94와 조건을 혼합하지 않는다.
C11/C12/C44=114/62/32 GPa는 기존 rounded 0 K Mishin targets를 유지했다.

L0=4.05/sqrt(2) angstrom, E0=1eV, target atomic area7.1024908428 angstrom^2.
불완전 fit은 b와 h를 함께 isotropic relaxation한다. L0와 rho_ref는 고정.
실제 평형의 atomic area/volume으로 J/m2/GPa를 환산한다. A_c 사용 안 함.

### C. 정확한 first-plane cancellation

tau ABC에서 3tau는 lattice vector이고 radial plane energy는 inversion symmetric:
w(h,tau)=w(h,2tau). Shockley partial의 첫 상대 평면 ISF 기여는 정확히 0.
G=0도 fixed-opening registry difference에서 상쇄된다.
첫 비자명 plane k=2의 최소 reciprocal exponential attenuation:
exp(-2h Gmin)=exp(-8 pi sqrt(2)/3)=7.1550809226e-6.
이는 prefactor 없는 attenuation이며 완전한 에너지 bound가 아니다.

동일 target geometry에서 scalar joint hybrid:
pair ISF=-.0001497678823 eV, embedding=-.0000072921670 eV.
Mishin source:
pair=.06941007609 eV, embedding=.00002460350 eV.
즉 이 source는 farther-shell pair radial shape로 큰 부분을 만든다.
그것을 tabulated pair로 교체하지 않는다. 현 LJ+scalar 환경의 구체적 결핍 근거다.

### D. 실제 재보정/식별성 결과

고정 density decay에서 u=4 eps sigma^12, v=4 eps sigma^6.
F(x)=-A sqrt(x)+B(x-1)+C(x-1)^2; F(0)=-B+C 포함.
Bulk5 독립 관측량 + interface5 특정 상태로 시작했다.
규격화: force .10eV/unit strain, cohesion2%, independent curvature5%,
interface energies10%. 측정 불확도/통계적 표준편차가 아니라 model discrepancy scale.

실제 실행된 nested 결과:
- sqrt: loss300.38761;
- linear:121.38024, Wsep 개선~1.80754J/m2지만 ISF 음수~-.00036268;
- convex C 추가:C=0, 개선 없음;
- positive two exponential density: single exponential로 퇴화, 같은 loss;
- positive squared envelope q(1-dq)^2: loss109.52176, ISF 부호 실패;
- signed C probe도 실패하고 LJ attractive coefficient boundary 선호;
- 고정 decay LP upper-bound feasibility는 다른 targets를 3scale 내 허용해도
  decay2.8에서 ISF<=-8.40874e-5eV/cell. 연속 전 parameter space 불가능 정리 아님.

다음은 별도 analytic angular research hierarchy이며 canonical EAM 변경이 아니다.
per-atom environment moment를 모두 합산한 뒤 회전불변 norm을 취한다.
Fourier transform H=2pi k exp(-dQ)(1+dQ)/Q^3, Q=sqrt(k^2+G^2)를 G로
미분하여 moments를 얻는다. H는 Bessel K_3/2 형태와 동일하다. Canonical cutoff 없음.

I3=2 sum_depth||Delta STF rank3||^2, I1=2 sum_depth||Delta vector||^2.
홀수차는 ANY affine centrosymmetric FCC bulk에서 정확히 0이라 bulk 탄성 수정 불가.
I2=2 sum_depth||Delta STF rank2||^2는 cubic FCC에서 0이지만 anisotropic
affine strain에 반응한다. SAME angular decay, 추가 radial range 없음.
noncubic background Q_bulk!=0인 rank2는 구현 범위 밖이므로 명시적으로 거부한다.
D1/D2 signed study는 total energy 안정성을 따로 검사한다.

I3 tied/independent range losses44.37151/28.29377 (10 targets).
I1+I3는 actual normal interface curvature를 11번째 target으로 추가:
loss33.87411, epsLJ8.81eV, actual C44=54.37 vs32GPa. 채택 불가.
I1+I3+C exact force/cohesion fit33.10337, v~1e-16이 cancellation floor 이하:
positive LJ라고 주장하지 않고 거부.
I1+I2+I3 exact fit29.70446 역시 v=0: 거부.

I1+I2+I3+C signed B/D1/D2 sector의 finite-LJ 후보:
scalar decay1.9388928377620271, angular4.877232382551694,
eps=.411271272477eV, sigma=2.304473908549angstrom,
alat4.05, cohesion3.36 exactly;
C11/C12/C44=86.435833/64.684362/50.670331GPa.
ISF=.134108780770J/m2;
heldout direct/partial RMSE=.04689350/.01349077J/m2 (각46점, training fractions 제외).
normal interface tangent147.82995 vs126.44968GPa.
Loss29.19092; full signed-logabs SVD condition4685.63,
두 exact bulk constraints tangent condition1833.91. Confidence interval 아님.

그 후보는 intermediate opening overshoot도 발견:
a/h=3 W=1.88623J/m2, Tn=-.65538GPa, separated limit~1.7502J/m2.
같은 source curve와 맞지 않아 채택하지 않는다.

### E. 실제 opening-path 수정과 수렴

monotone_opening_calibration은 같은 targets/weights를 유지한 채,
analytic W_a>=0을 declared grid a/h=1.1,1.2,1.5,2,2.5,3,4,5,8에서 부과했다.
처음 grid는 a/h3.5에서 negative force -.002706 eV/coordinate를 놓쳤다.
run_monotone_opening_refinement는 analytic W_aa=0을 찾아 actual force minima를
추가하고 fixed radial parameters에서 coefficient QP를 다시 풀었다.
4회의 exchange, 독립80/160 bracket, reciprocal tolerance refinement를 수행했다.
마지막 min force~-1.03e-14, a/h3.57358204, observed arithmetic/refinement floor
수준이며 raw 값은 유지한다. Force clipping/strain clipping 없음.

최종 refined 후보 JSON:
monotone_opening_refined.json.
scalar decay1.287051801314886, angular4.898979485566356;
coeff order [u,v,A,B,C,D3,D1,D2], 전체 정밀값은 JSON을 읽는다.
Loss29.47407475, normal-force max~1.240518eV/coordinate at a/h1.1904995.
검사 구간1.001..12, 전체 continuum monotonicity interval-proof는 아니다.
거의 0 견인력 plateau가 reference와 동등하다고 주장하지 않는다.
Independent Al elasticity/heldout fit는 여전히 부족하다. PDE 승격 안 함.

최종 refined 후보의 실제 재실행:
alat4.05 angstrom, cohesion3.36eV/atom, epsLJ=.154457514727eV,
sigmaLJ=2.452464388775angstrom;
C11/C12/C44=87.8161695/64.8326341/49.2710486 GPa;
ISF=.131005541319 J/m2, separation1.763013746519 J/m2;
normal/registry local tangents147.2948633/31.8052565GPa.
Held-out direct/partial RMSE=.05270159259/.01601694045J/m2.
Full signed-logabs condition1093.33. Bulk-equality-only tangent condition580.08은
active opening inequalities까지 포함한 confidence 계산이 아니다.
Fixed-s normal traction peak9.7715342GPa, source12.9695677GPa.
이 수치로 coupled dynamic slip/opening ordering을 주장하지 않는다.

### F. 실제 물리/전산 시나리오

모든 stress 입력은 explicit normal traction/resolved shear GPa:
normal tension .05..10, compression -.1..-2, pure shear .05..4,
explicit45degree axis: normal=shear=sigma/2,
static load/unload [0,0]->[.5,.2]->[1,.4]->[.5,.2]->0.
direct110/partial112 두 경로, a/s Hessian, derivative checks, 2D surfaces를 실제 계산.
큰 GPa는 ideal mechanism probe이지 보정된 피로 실험 조건이 아니다.

Continuation step .2/.1/.05GPa 비교. 이것은 PDE dt refinement가 아니다.
실패한 root를 spinodal/registry jump로 위장하지 않는다.
Fixed s normal saddle는 coupled2D saddle 또는 dynamic event ordering이 아니다.
Quasistatic unload return은 동적 zero-stress hold/잔류 소성 증거가 아니다.

Finite-q STATIC second variation K(q):
LJ+scalar EAM pair redistribution, F'' term, odd/even moment gradients 포함.
odd derivative는 cos-1, even derivative는 sin; 그 차이를 시험했다.
Radius5/8/12 L0, GX/GL/GK 각20점; 지금까지 계산한 candidates 모두
sampled positive. 직접 sinusoidal displacement energy variation과 독립 검증.
전체 BZ/finite-amplitude stability 증명 아님, 원자질량/Hz 사용 없음.
마지막 후보는 direct radius20도 추가했다: 최소 eigenvalue .08074286975,
radius12의 .08074197056과 같은 양의 부호이며 차이~8.99e-7.
3개 high-symmetry directions의60개 sampled wavevectors만 검증한 것이다.

### G. 데이터·문서·재현 명령

기존 결과 audited_v2와 기존 Al calibration/static parameter 파일은 보존한다.
새 결과 ROOT=results/fcc111_active_interface/matched_v3.
하위 odd_moment/even_moment/monotone_opening은 각 후보의 새 실행이다.
fitted_candidates, nested profiles, residuals, SVD/correlations, stress_scenarios,
opening/registry curves, constrained barriers, finite-q, coupled grids를 저장한다.
원시 fit와 summary 재생을 구별한다.

문서:
- MATCHED_INTERFACE_CALIBRATION.md: source mapping, 실패/수정/fit/응력;
- ANALYTIC_ANGULAR_ENVIRONMENT.md: rank1/2/3 infinite transforms, derivatives, K(q);
- STATISTICAL_CORRELATION_AREA.md: covariance-derived area와 void probability는
  자동 동일하지 않음, local N=1은 A_c 식별 불가. K(q)^-1도 곧바로 crack covariance 아님.

재현 (실제 Python3 interpreter 사용):
python -m solver_v1.reference_eam_targets --download
python -m solver_v1.run_matched_interface_study --phase calibration
python -m solver_v1.odd_moment_calibration
python -m solver_v1.run_constrained_odd_calibration
python -m solver_v1.even_moment_calibration --convex --signed --free-linear
python -m solver_v1.monotone_opening_calibration
python -m solver_v1.run_monotone_opening_refinement
python -m solver_v1.run_matched_interface_study --phase monotone-scenarios
python -m solver_v1.run_extended_identifiability
python -m solver_v1.run_static_bulk_stability
python -m solver_v1.run_matched_research_summary

마지막 summary 명령은 저장된 fit 요약/실제 separation 계산이지 최적화 재실행 아님.
Source 없으면 다운로드 허가와 checksum 확인 후 실행하거나 unavailable 보고.
A_c, physical mobility, empirical fatigue/yield를 어떤 fit에도 넣지 않는다.

### H. 검증/다음 단계

이번 턴 actual records (최종 갱신은 matched_v3/validation_record.json 참조):
targeted35 passed4.90s;
최종 full solver212 passed491.96s;
최종 app27 passed2 skipped50.96s;
desktop smokePASS a0=.7713438268704838,kappa86.29296488740997.
이전 full199/202/209도 실행했으나 최신 tests의 대체가 아니다.
처음 rank2 perfect-force 테스트에서1.1e-20 roundoff vs1e-20 assertion 실패가
있었고 arithmetic tolerance로 수정했다. 물리 force=0 clamp는 없다.
최종 smoke와 변경 코어 불변 diff도 실제 통과했다.
stage된 새 파일까지 diff --check하고 정상 commit/push 한다.

물리적으로 수용할 수 있는 수준:
- 검증된 analytic framework + 개선된 STATIC mechanistic candidate.
- 아직 quantitative Al calibration/unique material parameters/kinetics 아님.
- Al C11/C44 및 held-out opening shape 충돌은 남아 있다.
- Future 방향은 독립 atomistic environment/relaxation targets와 최소 analytic
  environment의 표현력을 함께 검사하는 것. 단지 polynomial 차수를 늘려 fit 금지.
- 이 결과를 졸속 PDE/UI 승격하지 말고 사용자에게 현재 solver readiness를 보고.
- 진짜 spatial mechanics/correlation 없이 mesh에 global probability를 칠하지 않는다.

## 아래는 fee1da7 이전 audited_v2의 보존 기록

아래 '현재/다음'이라는 문구는 당시 기록이다. 후속 근거와 결정은 위 절을 우선한다.

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
