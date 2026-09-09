# CURRENT_WORK_HANDOFF.md — 단계별 검증 후 재개하기

## 최신 상태: kinetics_loading_audit — 시간·하중·GPa 원인 감사

이 절이 아래 discrete_screw_v6보다 최신이다. 요청은 “초/Hz 사용, 전단·차원
지원 여부, 반복 제외2개, 실제 소재 대비 GPa 응력의 원인 수정”이었다.
과학 worktree aft-pde-bessel-38969ad, probability-pde-solver-v1에서 수행.
시작 fresh fetch/HEAD/origin 모두09e0dc6ce85e0a81656a92995cb3c0816bb96e97,
작업 트리 깨끗함. 원래 detached 사용자 폴더의 다른 파일은 보존한다.

### 시간 연결에서 실제 수정한 결함

- 기존 변환식은 맞았지만 UI에 보정 파일 불러오기 기능이 없었다.
  이제 Load kinetic calibration / 동역학 보정 불러오기로 명시된 JSON을 읽는다.
- 기존 physical mode는 calibration의 온도/상대 이동도가 달라도
  generator를 M=(1,.05),kT=.02로 만들었다. 이제 선택한 물리 보정의
  kBT/E0 및 두 M*를 같은 energy/PDE에 적용한다. model mode는 기존값 그대로.
- 모델 ID, static parameter fingerprint, 동일 a/s cell 좌표, eV 및 Al-target
  L0를 검증한다. 다른 모델/온도누락/잘못된 ratio/t0를 거부한다.
  모델 변경시 incompatible physical mode 해제. 언어/보정 로드로 기존 결과의
  시간축/배열/zoom을 바꾸지 않는다. basis 변경은 동일 model frequency 유지.
- solver_v1/kinetic_calibration_workflow.py:
  C(t)=exp(-B t) C0, B=M H, C0=kBT H^-1,
  B=-log(C(t) C0^-1)/t, M=B H^-1의 행렬식으로 실제 CSV를 분석한다.
  여러 lag, 평형 covariance, diagonal-M, real log, 재구성 오차를 검사.
  t0는 M_a*=1로 정하고 M_s*에는 측정 mobility ratio를 쓴다.
- import_kinetic_calibration CLI: meta JSON + physical covariance CSV ->
  새 clock JSON + fit diagnostics. 덮어쓰기 금지. sample은 tests에서만
  exact synthetic C; Monte Carlo/가짜 Al 자료 없음.
- 문헌은 dislocation centre mobility/line drag를 다룬다. a/s cell mobility와
  에너지 정규화/단위부터 다르고 M_a도 주지 않으므로 대입하지 않았다.
  실제 repository Al M_a,phys/M_s,phys/t0는 여전히 unavailable,
  default kinetic JSON unchanged. 물리 시간이 보정됐다고 보고하지 말 것.

### 전단/차원 및 UI 실제 지원

- Desktop 여전히 scalar axial sigma/E + chi=.2, P(a,s,t).
  독립 normal/shear/phase는 low_stress_v4 연구용 reflecting SG에만 있다.
- UI tensile_direction은 검증만 하고 실제 config/PDE에 전달하지 않는
  placeholder였는데 적용된다고 설명했다. 설명 수정, 입력 disabled/not connected.
- 실제 capability flags를 result/conversion metadata에 저장한다.
  independent_shear_input, orientation_input_active, vector_registry_pde,
  spatial_specimen_solver는 현재 false.
- 연구용 resolved_tensor_tractions 추가: 대칭3x3 sigma와 명시된 동일frame
  n,m에서 [n.sigma.n, m.sigma.n, (n cross m).sigma.n] 반환.
  두 전단의 존재를 숨기지 않지만 이것을 3D PDE 완료라고 부르지 않는다.
- Full3D FCC geometry, transverse2D/one-component harmonic screw,
  2D probability state를 구별한다. CAD/mesh UI gate는 그대로 닫혀 있다.

### GPa 원인과 실제 새 계산

같은 연구 후보/파라미터와 같은 rigid-interface Mishin source를 비교했다.
원자적 work=T*Area_atomic*L0, MPa/Pa/eV 변환을4/15/25/50/100/2000MPa에
독립 재계산. factor1000 오류 없음; saved work error0, round-trip roundoff뿐.
Production kappa*sigma/E는 변경하지 않는다.

고정 normal gap만의 peak를 물리 임계로 쓰지 않고 W_a=f_n,W_aa>0의 branch,
Phi''=W_ss-W_as^2/W_aa로 첫 scalar-path shear spinodal을 계산한다.
Normal traction=0,33/65 bracketing 실제 실행:
- analytic direct110: fixed first7543.359226MPa -> normal-relaxed4768.897903MPa.
- source direct110: fixed6461.635050 -> relaxed4390.365017MPa.
- analytic Shockley112: fixed first3012.198345 -> relaxed2685.841510MPa.
- source Shockley112: fixed 약2405.245MPa -> relaxed2345.062485MPa.
뒤쪽 registry peak(Shockley에서~14/11GPa)를 첫 불안정점으로 사용하면 안 됨.
Source의 coarse fixed peak bracket에 curvature root가 여러 개 있어서 brentq
한 번으로 뒤쪽 근을 선택한 것도 발견. report script에서 구간을 세분화하고
첫 +->- curvature root를 검사한다. 최종 first_branch_comparison.csv 확인.

normal equilibrium residual<=1.7e-14 eV/reduced coordinate,
relaxed Schur residual<1e-12. 동일 제한조건의 source도GPa이다.
GPa 이상 강도와 defect-bearing 실험 시편의MPa 항복/피로를 동일시한 것이
핵심 물리 비교 오류. Normal accommodation/path 구속도 값을 높였다.
동시에 candidate C11/C12/C44=87.816/64.833/49.271GPa vs114/62/32GPa로
재료 적합은 불충분. 수직 이완으로 낮아졌다고 Al fatigue 채택하지 않는다.
기존 defect MPa 구동력은 static/line energy이지 kinetics나 nucleation 인증 아님.
에너지 파라미터/온도/M을 보기 좋은 소성을 위해 바꾼 적 없음.

### 재현 / 원시 기록 / 최종 실행된 검증

- solver_v1/KINETICS_LOADING_AND_STRESS_AUDIT.md 자세한 식/근거/현재 한계.
- python -m solver_v1.run_stress_scale_audit (실제8 curve runs,307.679s).
- python -m solver_v1.report_stress_scale_audit (첫 fixed peak 추가곡률 검증/그림).
- results/kinetics_loading_audit/ CSV/JSON/PNG. 출처·기존 파라미터 SHA 저장.
- targeted solver time+stress: **17passed /2.32s**.
- combined targeted35passed/61.78s checkpoint 뒤 native Tcl read 오류가 재현됨.
  mixed targeted35passed/1failed 및 app30passed/1failed(162.11s)는 실패 기록이며,
  그보다 이른 app31passed/152.50s만으로 최종 성공을 선언하지 않았음.
- 실제 Tk12회/App10회+12회 단독 생성/종료 정상이어도 mixed 실행에서
  두 번째/세 번째 Tcl 인터프리터 script read 접근 오류가 남았다.
  최종 수정은 한 module-scoped Tk + 테스트별 Toplevel/root injection.
  Production도 원래 한 Tk이다. 반복 retry나 Windows Tcl skip으로 숨기지 않는다.
  실제 platform 접근 오류의 모든 원인을 규명했다고 주장하지 않는다.
- 최종 lifecycle targeted: **10passed /2.37s**.
- 최종 full solver: **271passed /1236.63s**, 실패0, 제외0, exit0.
  .cache/time-shear-full.xml 및 실행 완료 출력을 직접 확인.
- 최종 full app: **31passed /161.26s**, 실패0, 제외0, exit0.
  .cache/time-shear-app-single-interpreter.xml 및 완료 출력 확인.
  실제 Tk3개 포함. 이전에 매번 제외되던2개도 실제 실행했다.
- 최신 desktop --smoke: **PASS /1.487s**. Smoke와 실제 Tk 검사는 별개.
- 위 결과는 results/kinetics_loading_audit/validation_status.json에도 기록한다.
  개인 machine 정보가 든 .cache JUnit XML은 커밋하지 않는다.
- static calibration/Bessel/production energy registry/kinetic default JSON unchanged.
- 이 절을 포함한 checkpoint의 commit/원격 상태는 git log/fetch/status로 재확인.
  최종 diff check 후 probability-pde-solver-v1에만 정상 commit/push한다.
  main 및 원래 detached 사용자 폴더의 다른 파일은 변경하지 않는다.

### 다음 단계 / 아직 해결되지 않은 것

- 사용자의 실제 a/s 집단좌표 covariance/relaxation 데이터가 필요하다.
  실제 Al M_a,phys/M_s,phys/t0를 제공한 적 없으며 seconds/Hz 기본은 disabled.
  소스의 물리적 대응까지 파일 validator가 자동으로 증명하는 것은 아니다.
- 독립 전단/normal 연구 입력은 있지만 production은 scalar axial load이다.
  독립 전단을 desktop에 켰다고, vector-registry/3D spatial PDE가 완성됐다고
  보고하지 않는다. UI의 가짜 orientation 적용 설명을 이번에 제거했다.
- 이번 완화된 ideal interface peak도 실험 시편 항복/피로 응력이 아니다.
  Matched bulk elasticity 오류, 실제 결함 생성/공간 연결 및 kinetics가 남았다.
  채택 gate는 미통과다. 단위를 바꾸거나 stress/M/온도/장벽을 조정해 통과시키지 않는다.


## 최신 상태: discrete_screw_v6 — 미해결 원인을 실제로 고치는 후속 감사

이 절이 아래 nonlocal_v5보다 최신이다. 사용자 최신 지시는
“왜 미해결인거 같냐 계속 원인 찾아서 고쳐봐 좀”이었다. 신규 agent 없음.
동일 a485bc4 기반 probability-pde-solver-v1 과학 worktree에서 이어서 수행한다.
원래 detached 폴더/사용자 파일은 보존. discovery 문서는 hash 일치 확인 후 동기화.
최종 코드 검증은 완료했다. 이 절을 포함하는 검증 checkpoint의 SHA와 원격 상태는
git log/status/fetch로 확인한다. 최종 응답에도 실제 commit/push 결과를 기록한다.

### 실제 고친 것

1. 폭 한 개의 arctangent trial -> 전체 profile의 constrained Euler 방정식을
   실제로 푼다. 기존 continuum 에너지는 그대로. Fourier preconditioning과
   L-BFGS 후 projected Newton-CG로 objective 종료와 실제 force 균형을 구별한다.
   `int s dx=Q` 제약의 반력 lambda는 holding traction; 내부 force error 아님.
   실제10개 profile의 최대 projected residual4.6914e-10MPa.
   기존 trial의 약1.2GPa 불균형을 인위적 mobility/energy tuning 없이 제거했다.
   하지만 이는 continuum PN의 stationary profile이지 atomistic core 인증 아님.
2. 새 discrete_fcc_screw.py는 line=e1을 따라 무한 LJ/Bessel sum을 정확히 하고,
   transverse row(j),layer(l)는 이산적으로 남긴다:
   R=(nb+(j+l)b/2, d(j+l/3),hl), d=sqrt(3)b/2.
   Scalar density와 rank1/2/3 STF derivative도 같은 Poisson 방식으로 전개한다.
   General scalar F''를 삭제한 것이 아니라 row rho_x=0인 anti-plane subspace만
   적용. full per-atom angular sum 후 norm 제곱; 기존 D1/D2/D3 모두 그대로.
3. 실제 bulk symbol K(q_y,theta), ABC phase, infinite-layer Green/Schur 및
   independent direct3D/finite-layer checks 구현.
   기존 W_int,ss=4.945211215219163과 row reconstruction=4.945211215219173:
   difference9.77e-15. R32 direct3D maxabs2.5865e-7,
   maxrelative5.2473e-6(작은 q 포함); R8/12/16/24/32 tail 감소 실제 비교.
   8rings/288rows/11modes, extra16rings 검증. Tail diagnostics는
   2 sum|Phi|+16 sum|D_r| S_r T_r의 같은 stiffness 단위로 산정한다.
   마지막ring1.4043e-14, omitted-validation bound3.3466e-17.
   이는 finite validation window bound + empirical infinite remainder 검증이지
   무증명 exact infinite error bound를 주장한 것이 아니다.
4. Relaxed K_jump(0)=4.795625552378889, rigid보다3.02486%낮다.
   논리적 인접 row와 spectral same-y interpolation은 finite-q에서 다른 제약이다.
   qL0=1:8.83393/8.08444; zone edge26.19934/22.34258.
   이전 continuum+rigid-local kernel을 원자 scale까지 쓰면 큰 오차다.
   Acoustic-pole slope도 따로 유도/검증:2.099137/1.089480(해당 reduced 단위).
   단순 mu|q|/2 coefficient 교체/fit은 하지 않았다.
5. K_jump는 local stiffness를 포함한다. 기존 gamma Hessian을 더하면 이중 계산.
   K0를 빼고 nonlinear rigid gamma를 넣는 것도 아직 정당화되지 않아 하지 않는다.

### 실제 낮은 응력 / refinement

- 4/15/25/50MPa, 4/8/16/32/128/512 atomic-row spatial periods,
  두 명시된 jump constraints의 실제 harmonic static response48개.
  Dynamic trajectory/Hz/피로 수명 계산이 아니다.
- 같은 continuum의 fixed-content d/b32/128/512/1024, L/d8, dx/b.0625:
  holding stress115.205590/28.231490/7.023611/3.508895MPa.
  외력0/4/5/10/15/25/50MPa에 대한 signed drive를 저장한다.
  이 제약에서 uniform 외력은 energy -tau*Q만 바꾸므로 여러 동적 실행인 척하지 않는다.
- d_content128,L/d8, dx/b.25/.125/.0625/.03125:
  9.842586/26.003904/28.231490/28.231613MPa.
  모든 force residual이 작아도 coarse mesh는 numerical pinning으로 틀렸다.
  마지막 refinement change.000123MPa. physical Peierls/소성 floor로 오해 금지.
- Mean content를 고정한 domain 증가가 실제 defect separation을 바꾼다는
  추가 protocol 오류를 찾아 수정했다. L/content4/8/16/32에서 실제 간격은
  127.6933/127.2580/126.4450/124.8193b였다. 서로 같은 결함이 아니었다.
  새 --matched-domains는 actual s=b/2 crossings를128b로 맞추는 outer bracket solve.
  현재 완료 L/d4/8/16/32:23.237335/28.049974/29.205530/29.491616MPa.
  L/d64도 완료:29.5629654769MPa 대 isolated29.5869013546MPa,
  남은 periodic effect.0809%. L/d8 실제 간격128b에서 dx/b.0625->.03125 변화
  1.95e-10MPa. 총38회 constrained profile solve,850.30s.
  실제 crossing 오차<=1.03e-9b. matched_separation_domain_refinement.csv 확인.

### 원시 파일과 재현

- solver_v1/DISCRETE_FCC_SCREW_DERIVATION.md에 자세한 식/정규화/오류 원인.
- solver_v1/discrete_fcc_screw.py: harmonic row/Bessel/Green reduction.
- nonlocal_registry_reference.py에 solve_fixed_registry_content 추가.
- run_discrete_screw_reference.py가 실제 연구 실행:
  `python -m solver_v1.run_discrete_screw_reference --kernel --profiles --matched-domains --report`
- results/fcc111_active_interface/discrete_screw_v6/에 작은 CSV/JSON/PNG/SVG.
- fitted candidate SHA256는 여전히
  9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655.
  production models/PDE/app/static fit/kinetic calibration 파일은 그대로다.

### 완료된 검증 체크포인트

- nonlocal_v5 전체 solver242passed(1117.29s)는 실제 완료했다.
- v6 targeted47passed(18.20s), tail-bound 별도15passed(11.66s).
- app27passed/2skipped(203.40s), desktop smoke PASS(2.626s).
- v6 첫 전체 실행은 tail-bound 단위 감사를 반영하기 위해 55%에서 중단;
  PASS 아님. 최신 코드를 fresh basetemp로 재실행:
  .cache/discrete-reviewed-full-regression.xml.
- 일부 첫 targeted 실패는 3D direct tail/작은q quadrature 해상도 부족이었다.
  tolerance를 풀지 않고 radius32 및 Ntheta8192 이상으로 검증했고 통과했다.
- matched domain 계산과 최종 전체 검사는 완료.
  최신 solver_v1 **260passed /903.34s**, app **27passed/2skipped /203.40s**,
  startup smoke PASS, a0=.7713438268704838, kappa=86.29296488740997.
  Tk의 init.tcl/display를 이 실행 환경에서 찾지 못해 rendering 관련2개만 제외됨.
  실제 GUI rendering을 검증했다고 하지는 않는다.
- git diff --check 및 --cached --check PASS. 새 Matplotlib SVG의 후행 공백을
  생성 코드에서 형식 정리하고 두 SVG XML parse도 확인했다. 수치/그림 내용 불변.
- 결과/명령/timing/status는 discrete_screw_v6/validation_status.json에 보존.
- 재확인 remote는 작업 전과 같은 a485bc4였다. main/origin-main 불변:
  80cacb4 /c43d8e0. 정상 fast-forward 외 push 방식은 사용하지 않는다.
- 이 checkpoint는 검증된 연구 단계의 완료일 뿐,
  nonlinear atomistic core / Al fatigue / physical kinetics 완료 선언이 아니다.

### 다음 진짜 물리 단계

Full nonlinear discrete row/cross-section 에너지 및 relaxed local GSF를
동일 constraint에서 유도해야 한다. 현재 exact harmonic + continuum profile을
그냥 붙여 atomistic core 완성이라 부르면 안 된다. Vector registry, normal
relaxation, finite-loop/patch의 유한 활성화 에너지, 독립 재료 적합/kinetics도 남아 있다.
현재 연구 후보 C11/C12/C44=87.816/64.833/49.271GPa 대 target114/62/32GPa:
재료 적합 자체도 여전히 불충분. 기존 파라미터 재보정/값 이동은 이번에 하지 않았다.
Straight line 에너지는 J/m이며 임의 선길이/A_c로 eV thermal barrier를 만들지 않는다.
M_a,phys/M_s,phys/t0 unavailable, seconds/Hz disabled, UI gate 미통과.
DISCRETE_FCC_SCREW_DERIVATION.md 마지막 절에 다음 nonlinear row functional을
명시했다. F'_bulk를 nonlinear 상태에 동결하면 안 된다. 각 site의 rho/Q를 먼저
합산해 F/norm을 적용하고, harmonic limit가 새 정확한 K와 일치해야 한다.
한 infinite row를 b만큼 옮기는 것은 atomic position relabeling과도 같으므로
unwrapped registry와 Burgers/winding/boundary 조건을 별도로 검증할 것.

## 최신 상태: nonlocal_v5 — 비국소 탄성 비용과 축 대응 감사

이 절이 아래 low_stress_v4 / matched_v3 기록보다 최신이다. 사용자 최신 요청은
남은 문제를 이어서 해결하라는 것이었다. 시작 fresh fetch에서 과학 worktree의
HEAD와 origin 모두 a485bc4ba977914a159acd2d59b256e5bfca399b, 깨끗한 상태였다.
실제 작업은 aft-pde-bessel-38969ad의 probability-pde-solver-v1에서 수행한다.
원래 사용자 폴더는 detached c43d8e0와 사용자 파일을 보존한다. 두 discovery 문서는
변경 전 hash 일치 확인 후 동기화한다. 최종 commit/push는 실제 git로 확인한다.

### 실제 구현 및 실행

- 같은 angular_monotone_opening 후보를 사용하며 재보정 없음. 기존 SHA256:
  9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655.
- 새 nonlocal_interface_elasticity.py: cubic C 회전, 두 half-space의 static
  ordered-Schur impedance, `K=|q|(Z_plus^-1+Z_minus^-1)^-1`.
  Bulk positivity, Hermitian/Riccati residual, translation K(0)=0 검사.
- 독립 isotropic 해석해, finite-depth FEM Schur, screw analytic reduction과 일치.
  LJ/Bessel의 continuum/long-wave 한계다. atomically exact finite-q 계면 kernel 아님.
- 새 nonlocal_registry_reference.py: 같은 analytic W_int의 검증된 Fourier gamma,
  에너지/gradient/Hessian, nonlinear intrawell Newton-CG static solve,
  명시된 기존 screw dipole의 제한 arctangent trial-width 최적화.
- 에너지 단위는 eV/L0 of dislocation LINE = J/m. 이를 eV activation barrier로
  사용하지 않는다. 임의 선 길이/집단 면적/A_c를 곱하지 않는다.
- 128 analytic samples, 64/128 비교, off-grid energy/force/Hessian 검사 실제 수행.
  최대 오차3.275e-15 /7.262e-13 /1.372e-10(해당 reduced 단위).
- bulk C=87.81616955/64.83263414/49.27104858 GPa. derived screw mu=23.79519794GPa.
  기존114/62/32GPa 타깃 불일치는 그대로이며 보정 성공으로 승격하지 않는다.

### 독립 검사로 찾아 고친 방향 버그

기존 +tau,+h FCC 원자 위치는 바꾸지 않았다. 다만 generated stack의 actual
cubic axes와 기존 geometric/lab e1/e2/e3를 동일시한 변환은 잘못이었다.
actual cubic에 대한 plane basis는 (-e1,-e2,e3)이며 proper rotation이다.
이 좌표에서 모든 sampled lattice vector는 (a_lat/2)*integer, even-sum FCC다.
`geometry.plane_basis_in_stacked_cubic_axes()` 추가,
`StaticBulkHessian.crystallographic_wavevector()` 수정.
옛 Gamma-X/L/K CSV는 그대로 두되 잘못 붙었던 labels를 문서에 명시했다.
RegistryPath.direction_3d는 원래 geometric frame; scalar path/에너지 불변.
잘못된 frame acoustic discrepancy37–40% -> 올바른 frame에서1.8525e-5 이하
(q L0=.005, independent radius32L0, 4directions).
새12/20/32L0의 60방향 점씩 corrected finite-q 검사는 모두 양수,
최소0.0669003eV/L0^2. full Brillouin zone/비선형 안정성 증명 아님.

### 낮은 응력과 기존 결함 — 완성된 피로가 아님

4/15/25MPa sinusoidal spatial traction, wavelength4/8/16/32/128/512b:
실제18개 static nonlinear solve. 최대force residual6.52e-8MPa,
static unload slip/b<=9.87e-16. short-wavelength atomistic 정확도는 미인증.
정적 unload는 동적 hold/잔류소성/주기 피로 검증이 아니다.

지정된 pre-existing opposite screw pair d/b32/128/512/1024:
d=9.164/36.656/146.626/293.251nm. L/d8, dx/b.0625에서 trial balance
112.2174/28.05040/7.01254/3.50627MPa.
각각에0/4/5/10/15/25/50MPa를 넣고 separation 에너지 미분을 실제 계산.
구동력이지 속도/생성/피로수명이 아니다. defect spacing은 주어진 시나리오다.
고정 normal-gap perfect registry의 ideal shear는7.543359GPa:
실험 yield나 mixed a-s spinodal 아님. 기존결함 이동과 pristine 생성은 다르다.

독립 수렴:
- d/b128,L/d8,dx/b=.25/.125/.0625/.03125:
  balance28.63133374/28.05266482/28.05040423/28.05040422MPa.
  마지막 차이1.143e-8MPa.
- d/b128,dx/b.0625,L/d4/8/16/32/64:
  23.23775308/28.05040423/29.20596108/29.49204783/29.56339690MPa.
  isolated far-field29.58690135MPa; 마지막periodic effect~.0794%.
  L/d8 한 번으로 domain-converged라고 하면 안 된다.
- analytic infinite Fourier elastic energy와도 비교; finest차이~1e-14.
- d/b1024, core widths.15/.28/.5/1b sensitivity spread~6e-5MPa.

하지만 trial width=.27789b(~.0796nm), full Euler residual~1.2GPa,
q90b~.67–1.23이다. 폭 하나를 최적화했을 뿐 실제 전위 core를 푼 것이 아니다.
숫자 mesh 수렴은 atomistic core 정확도/유한 nucleation barrier를 인증하지 못한다.
Candidate/0K source elastic mu=23.795/28.844GPa로 d293nm isolated attraction은
3.698/4.483MPa. 따라서4MPa에서 expansion 여부조차 material comparator에 따라
바뀐다. 낮은 응력의 Al 소성 검증 성공을 주장하지 않는다.

### 남은 다음 단계와 금지 사항

1. 같은 infinite LJ/Bessel + many-body에서 discrete half-crystal force constants/
   Green/Schur kernel을 유도. 실제 interface jump constraint B, ABC Bloch phases,
   zero-mode/gauge 및 relaxed local Hessian을 먼저 일치시킬 것.
   `(B H^+ B^T)^-1`은 명시된 constraint에 대한 선형식이지 현재 완성된 코드 아님.
2. 원자적 core, vector registry/partial splitting, normal relaxation 검증.
   현재 rigid direct110 scalar path를 최저에너지 2D GSF 경로라 하지 말 것.
3. finite loop/patch와 실제 spatial activation energy를 유도한 뒤 확률 이론 연결.
   현재 straight-line J/m에 임의길이를 곱해 thermal probability 생성 금지.
4. 독립 material fit/validation 및 collective kinetic data는 여전히 미완료.
   M_a,phys/M_s,phys/t0 unavailable, seconds/Hz disabled. A_c도 별도 미보정.
5. production PDE/energy registry/static parameter/UI는 그대로. UI gate 미통과.
   사용자 확인 없이 geometry-import/mesh UI로 넘어가지 않는다.

### 파일 / 재현 / 검증 상태

핵심 solver_v1/NONLOCAL_INTERFACE_ELASTICITY.md.
results/fcc111_active_interface/nonlocal_v5/에 작은 CSV/JSON/NPZ/PNG/SVG 저장.
`python -m solver_v1.run_nonlocal_interface_reference --prepare --scenarios`
가 실제 Bessel+탄성+정적 시나리오 실행이다.
`python -m solver_v1.report_nonlocal_interface_reference`는 저장 결과 summary/plot.
물리/수치 결과는 numerical_and_physical_status.json에 명시적으로 분리한다.

이번 targeted(새 이론+geometry+기존 finite-q)29passed(4.02s).
App27passed/2skipped(126.65s), desktop smoke PASS(1.7524s),
TwoRowLJ a0=.7713438268704838,kappa=86.29296488740997 불변.
첫 full-suite 실행은 새 안정성 보호 검사를 넣기 위해 중단했으므로 PASS 아님.
이 nonlocal_v5 full solver는 .cache/nonlocal-solver-final-regression.xml로
242passed(1117.29s) 실제 완료했다. 그 뒤 사용자 지시에 따라 위 discrete_screw_v6로
작업을 이어갔으므로 이 절만 보고 최신 검증/완료 상태를 추정하지 않는다.



## 이전 완료 상태: low_stress_v4 실제 저응력 주기 검사

이 절이 아래 matched_v3 / audited_v2 기록보다 최신이다. 이번 시작은 fresh
fetch로 local/origin 모두55a0d514eee99ec3f141a972b07d477135fe62a7임을 확인했고
과학 worktree는 깨끗했다. 실제 위치는 기존 aft-pde-bessel-38969ad worktree의
probability-pde-solver-v1. 원래 폴더는 detached c43d8e0 + 사용자 미추적 파일을
보존한다. 커밋/원격 최종 상태는 git log/status/fetch로 확인한다.

사용자의 최신 질문은 “혼합하중 결과를 일반적인 금속 피로라고 볼 수 있는가”,
“2 GPa는 너무 크지 않은가”였고, MPa 범위에서 실제 주기 검사를 하도록 승인했다.
과거 큰 static 미션을 다시 시작하거나 UI를 재설계하지 않는다.

### 이번에 실제 완료한 계산

- 기존 refined angular_monotone_opening 계수를 그대로 사용. 새 fit 없음.
- parameter file SHA256:
  9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655.
- LJ/Poisson/Bessel, 모든 production energy/PDE, UI, A_c, mobility, static fit 파일 불변.
- 새 low_stress_cyclic_diagnostic는 생산 energy selector 연결이 아닌 별도
  reflecting-box 가설/반증 검사. 검증된 SG energy-array generator만 재사용.
  Opening sink 없음: opening_probability=null, 균열 없음/확률0이라는 뜻 아님.
- 온도293.15 K, kT=.02526171246eV를 명시. 단위가 맞는 room-T 가설이지
  finite-T material calibration이 아니다. M*=(1,.05), 실제 초/Hz 없음.
- 0,4 MPa resolved shear,10/20/30/50 MPa axial의45도 명시 예제,
  15/15 MPa90도 비비례 혼합, normal-only/shear-only/reversed-shear controls.
- 실제69개 cyclic/hold 프로토콜. grid31x75,61x75,61x147,91x219;
  정상 domain2h->2.5h; well3/5/7; steps64/128/256/512; explicit/implicit 짧은 비교.
- static direct110/partial112 두 주기씩10–50 MPa: 모든 점 같은 stable well,
  repeat-coordinate error<=2.31e-15. 이것은 dynamic hold와 별개다.

### 해석의 핵심

period40,16cycles,61x147,128steps 조건에서 zero-load P_out=.02202894954.
30 MPa P_out=.02203167161, paired excess2.7220677e-6.
50 MPa excess6.6030413e-6. Raw outside population 대부분은 무하중 열적 퍼짐이다.
그러나 작은 directional transfer는 수렴해 남는다. 무조건 “전부 roundoff”도 틀리다.
91x219,256steps30 MPa paired excess2.6984433e-6; observed envelope1.6216931e-8.
hold n9.1177955e-7; observed envelope1.8402984e-8. 통계 confidence나 물리 인증 아님.
3/5/7well hold n=9.14542e-7/9.32945e-7/9.33098e-7.
7well zero hold320->1280model time 후9.33093e-7, xi/h~1.4e-13.
따라서 이 가설 안의 작은 persistent registry memory는 있으나 residual Al plasticity
검증이라고 할 수 없다. period4 짧은 hold endpoint는 xi 회복이 덜 되어 residual 선언 금지.
Shear를 뒤집으면 n부호가 뒤집힘; normal-only는 n~-9e-15.
Global Gibbs 무하중 control은 stationary. 중앙 well에 제한한 초기 Gibbs를 해제한
열확산을 fatigue accumulation으로 오인하지 않는다.

Gross SG plus+minus는 grid spacing ds에 대해~1/ds인 Brownian recrossing traffic이다.
Signed flux/occupation balance는 유효하지만 gross를 물리 hop 횟수라 하면 안 된다.
실제 모든 run 최대mass residual9.6541e-12, well step balance4.2262e-15,
registry balance9.4560e-16, decomposition2.0007e-16; repair0.

### 아직 미완료인 물리 문제 / 다음 선택

이 에너지는 rigid infinite-interface의 원자 cell당 에너지다. 이것을 국소
thermal activation barrier로 간주하는 집단좌표/공간 정규화는 아직 유도되지 않았다.
국소 slip nucleus/결함, 주변 탄성 cost, 독립 재료 적합, collective kinetic data가
필요하다. 이를 A_c나 임의 energy multiplier로 고치지 않는다. LJ/Bessel에서
비국소 interface Hessian/Schur kernel을 유도하는 경로는 문서의 제안이지 완성된 solver 아님.
새 모델을 “일반 Al 피로/균열 모델 완성”이라고 보고하지 않는다. UI gate 미통과.
생산 TwoRowLJ/reduced hybrid와 기존 opening bookkeeping은 그대로 유지한다.

### 파일 / 재현 / 백업

핵심 문서 solver_v1/LOW_STRESS_CYCLIC_AUDIT.md.
원시 출력 results/fcc111_active_interface/low_stress_v4/.
execution_manifest.json의 CLI로 실제 재실행; summary CLI는 계산 재실행이 아니다.
history.csv.gz는 무손실. 모든 cycle 통계와 선택된 full-precision density snapshots
저장. 기존 full snapshots57개는 .cache/low-stress-full-snapshots에 보존 후
명시된 cycle index로 repository output만 축소. 사용자/기존 결과 삭제 없음.

이번 새 targeted tests11passed(1.26s). 첫 tmp_path 실행은 OS Temp permission error
(10passed/1error)였고, worktree 내부의 별도 .cache basetemp에서11passed 재검증했다.
App27passed/2skipped(155.84s), desktop smoke a0/kappa historical값으로 통과.
Full solver223passed(1192.85s), app27passed/2skipped(155.84s), smoke1.1416s로
실제 완료했다. git diff --check 및 staged --check 통과. 모든57개 full-snapshot
백업의 SHA256도 독립 확인했다. 결과 validation_status.json에 수치/물리 상태를 분리했다.
이 문서의 후속 커밋에 보존하며 commit/push 최종 SHA는 실제 git로 확인한다.
물리적 피로 검증/solver 완성/production 승격은 여전히 NOT VALIDATED다.

## 이전 완료 상태: matched_v3 연구 감사 (2026-09-09)

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
