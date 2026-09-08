# CURRENT_WORK_HANDOFF.md — 현재 작업 상세 인계

기록일: 2026-09-09. 이 문서는 실행 상태/미완료 항목을 보존하는 인계 기록이다.
영구 이론/행동 원칙은 AGENTS.md를 따른다.

## 1. 이번 사용자 요청과 중단 지점

사용자는 진행 중인 full-FCC calibration/active-interface 작업 도중,
먼저 이론·수식·원칙을 담은 agent 지침 파일을 만들고 다음 실행을 더 높은
추론 강도로 진행할 수 있게 해 달라고 요청했다.

따라서 이번 문서 작성은:
- 이전 구현/결과를 삭제하거나 재시작하지 않는다.
- 미완료 보정 작업을 완료로 선언하지 않는다.
- PDE production 경로/기존 calibration/UI를 변경하지 않는다.
- 추론 설정을 변경했다고 주장하지 않는다.
- 과학 작업 전체가 끝나기 전 임의 commit/push하지 않는다.

## 2. Git / 작업 위치

기존 full-stack 완료 기준과 현재 calibration 기반:

    b819abea7e8514d04698e6ab0e391b940262d16c

문서 작성 중 실제 git fetch origin이 성공했고,
origin/probability-pde-solver-v1도 위 SHA임을 확인했다.

현재 작업은 별도의 detached worktree에 있다. git worktree list에서
basename이 aft-pde-bessel-38969ad인 worktree를 찾아 사용한다.
이 경로는 환경별 임시 디렉터리에 있으므로 이 문서에 개인 절대 경로를 고정하지 않는다.

원래 사용자 workspace는 별도 detached HEAD:

    c43d8e096f2a329c7cdfad31e546c70fb36fe592

이며 docs/, examples/, fem1d/, libraries/, output/, paper/, research/,
results/, simulations/, tests/, theory/, tools/ 및 requirements 계열의
미추적 파일/디렉터리가 있었다. **삭제/정리/덮어쓰기 금지.**

현재 과학 작업의 새 파일은 원래 workspace가 아니라 위 별도 worktree에 있다.
두 폴더를 혼동해 잘못된 baseline에서 작업하거나 main을 수정하지 않는다.
현재 AGENTS.md와 이 인계 문서는 발견 가능하도록 두 위치에 같은 내용으로 저장한다.

## 3. 완료된 이전 full FCC reference

기존 추적 파일:
- solver_v1/fcc111_geometry.py
- solver_v1/fcc111_lattice_sum.py
- solver_v1/fcc111_full_energy.py
- solver_v1/FCC111_FULL_STACK_DERIVATION.md
- results/fcc111_full_stack/

이전 full-stack 완료 커밋:
- 5ec204fc2eb68e30e9db785b1d914f9a8d6021d4
- 92310c97d6930a9be38bb9cf88248e4a951ead98
- b819abea7e8514d04698e6ab0e391b940262d16c

기존 reduced Al-target parameters를 full FCC에 대입하면 stable normal root가
없다는 이전 결과를 보존한다. 무리한 parameter transfer를 고치려고 reduced
보정 파일을 변경하지 않는다.

## 4. 현재 미커밋 변경

Modified:
- README.md
- solver_v1/README.md
- solver_v1/fcc111_full_energy.py

New implementation/tests:
- solver_v1/aluminum_full_fcc_calibration.py
- solver_v1/fcc111_active_interface.py
- solver_v1/run_aluminum_full_fcc_calibration.py
- solver_v1/test_aluminum_full_fcc_calibration.py
- solver_v1/test_fcc111_active_interface.py

New draft documents/results:
- solver_v1/ALUMINUM_FULL_FCC_CALIBRATION.md
- solver_v1/FCC111_ACTIVE_INTERFACE_DERIVATION.md
- results/aluminum_full_fcc_calibration/
- results/fcc111_active_interface/

이번 요청으로 AGENTS.md와 CURRENT_WORK_HANDOFF.md를 추가한다.
기존 파일의 체크섬/내용을 보존하고, 아래 문제는 다음 과학 감사에서 수정한다.

fcc111_full_energy.py의 변경은 find_equilibrium=True 옵션 추가다.
False이면 root search를 생략하고 a0는 NaN/명시적 상태가 된다.
기존 True default의 회귀를 유지해야 한다.

## 5. 현재 구현 개요

### Bulk calibration

aluminum_full_fcc_calibration.py:
- full FCC target/parameter/observable/fit dataclasses;
- density normalization을 ideal target FCC 환경으로 고정;
- log-parameter bounded deterministic least squares;
- normal alpha, equal-biaxial beta strain에 대한 energy finite differences;
- equilibrium check와 central log-parameter sensitivity Jacobian.

현재 일반적인 strain finite-difference step은 4e-4,
sensitivity log step은 2e-4다. 최적화/식별성 검증 시 step dependence를 확인한다.

run_aluminum_full_fcc_calibration.py는 **저장된 fit vectors를 사용해 결과표와
그림을 재생**한다. 전체 staged optimization을 처음부터 재실행하는 스크립트로
오해하면 안 된다. 실제 optimization starts/진행을 재현 가능한 별도 경로로
정리할 필요가 있다.

### Active interface

fcc111_active_interface.py:
- 한 interface의 opening/slip;
- cross-interface pair multiplicity k;
- reciprocal-zero power tail의 Hurwitz zeta 처리;
- depth별 per-atom density 변화 후 embedding 합산;
- analytic/semi-analytic gradient/Hessian;
- separated half-crystal limit;
- independent finite direct_reference evaluator.

현재 production PDE/desktop energy selector에 연결하지 않았다.

## 6. 중요: 다음 실행에서 먼저 바로잡아야 할 이론/검증 문제

### 6.1 Structural rank — 최고 우선순위

현재 draft와 test에 “rank 5이지만 practically non-identifiable”라는 주장이 있다.
특히 test_full_fcc_sensitivity_is_full_column_rank_but_ill_conditioned는
np.linalg.matrix_rank(jacobian)==5를 요구한다.

그러나 사용된 cubic strain mapping은:

    X=C11+2C12
    Haa=V(X+4C44)/3
    Hbb=4V(X+C44)/3
    Hab=2V(X-2C44)/3
    Hbb=2Haa+Hab

따라서 탄성 관측량 세 개는 독립 조합 두 개뿐이다.
C11-C12를 제약하지 못한다.

Perfect cubic configuration에서 equal-biaxial strain-force는
normal strain-force의 두 배다. 두 force residual도 독립 정보 두 개가 아니다.

현재 bulk residual의 독립 정보는 이상적인 symmetry limit에서 많아야:
1. equilibrium,
2. cohesive energy,
3. X,
4. C44

이다. Five-parameter family를 이것만으로 structurally identifiable이라고
부르면 안 된다. 작은 비영 singular value는 finite-difference/truncation에
의해 생성될 수 있다. Stage 1의 “rank 3” 주장도 같은 이유로 재검토한다.

해야 할 일:
- analytic observation dependencies를 regression test로 추가;
- derivative-step/tail refinement에서 smallest singular value 변화 확인;
- structural rank와 numerical threshold rank를 구분;
- 중복 residual의 weighting/correlation 처리;
- 필요하면 별도 symmetry-breaking strain mode를 도입해 C11-C12를 식별하되,
  geometry 확장의 정확성을 먼저 검증;
- 이론/테스트/CSV/summary의 불일치를 정정한다.

기존 잘못된 rank assert를 무조건 보존해서는 안 된다. 왜 잘못됐는지 입증하고
독립정보/수치 floor를 시험하는 더 정확한 regression으로 대체해야 한다.

### 6.2 Fit failure와 전체 family 불가능을 구분

현재 draft는 square-root family가 “structurally insufficient”라고 단정한다.
현재 기록된 local deterministic fit 실패만으로 global impossibility가
증명되는 것은 아니다. 추가 starts/continuation/analytic incompatibility가
없으면 “현재 시도한 설정에서 residual이 남았다”로 제한한다.

두 bulk-equivalent 후보의 GSF/surface held-out 실패 역시 전체 family의
불가능성 증명은 아니다. 과대주장하지 않는다.

### 6.3 더 필요한 검증

- Reciprocal shell/layer/tail tolerance 및 interface embedding neighborhood의
  독립 refinement 비교. 한 direct radius/layer 비교만으로 충분하지 않다.
- Real-space validation radius/layer 증가에 따른 수렴표.
- Static loaded slip/opening saddle 및 spinodal: 현재 zero-load GSF maxima와
  normal traction peak만으로 모든 loaded coupled barriers가 계산된 것은 아니다.
- Work of separation, finite-load opening saddle, maximum traction을 구분.
- Unrelaxed half-crystal surface와 relaxed atomistic target의 mapping 조건 확인.
- 소수 저장된 후보의 low residual이 unique Al potential을 뜻하지 않도록 수정.
- Staged fit의 실제 parameter changes/starts/convergence log 재현성 확보.
- Bulk elastic strain Hessian의 positivity가 모든 crystal modes의 안정성을
  증명하지 않는다는 점 명시.

## 7. 현재 수치 결과 — 모두 잠정 감사 자료

이 수치는 저장된 결과를 찾기 위한 참고다. 위 식별성/검증 문제가 해결되기 전
“최종 Al calibration” 또는 “물리적으로 검증된 interface”로 인용하지 않는다.

Target set:
- 0 K Mishin EAM bulk: a_lat=4.05 Å, cohesion=3.36 eV/atom;
- C11=114, C12=62, C44=32 GPa;
- L0=b=2.8637824638 Å, h111=2.3382685902 Å;
- A_atomic_cell=7.1024908428 Å²;
- Haa=12.6460389542, Hbb=37.3161805205, Hab=12.0241026121 eV.

문헌 수치는 기존 database/source documents를 다시 검증한다.
DFT/EAM surface/GSF 타깃은 method/relaxation/orientation를 혼합하지 않는다.

Parameter order:
(epsilon_LJ[eV], sigma_LJ/b, density_decay_b, A[eV], B[eV])

Square-root recorded candidate:
(0.0038583102832933933, 1.1681262990758785,
 0.5122315833354125, 3.68017690186849, 0)
- a0/b≈0.81673594;
- recorded objective≈29.8602;
- normalized elastic residuals≈(3.6100,0.2468,-6.8277).

Linear candidate start0:
(0.051463163185591194, 0.9278241602685452,
 4.6853176139362835, 5.679412614403284, 2.760290609784449)
- a0/b≈0.816496172881;
- recorded least-squares cost≈6.47e-10.

Linear candidate start1:
(0.013645377854871132, 1.02668545888573,
 4.711004774121304, 5.655895464976975, 2.30373517720161)
- a0/b≈0.816496182069;
- recorded cost≈6.16e-10.

Reported central-Jacobian condition numbers≈4.56e7 and 3.10e7.
**True structural rank requires the correction in section 6.1.**

Held-out interface results for those two candidates:
- direct_110 USF≈0.438–0.449 J/m² vs recorded DFT reference0.250;
- Shockley USF≈0.0764–0.0768 J/m² vs0.224;
- ISF≈-0.00027 to -0.00024 J/m² vs0.164;
- work of separation≈1.089–1.146 J/m² vs DFT2.12 or EAM1.74;
- local reduced-coordinate interface Hessian diagonals≈
  (16.3950,2.89469) and (16.3656,2.86957), mixed≈0.

조금 음수인 ISF를 rounding으로 “정확히 0” 또는 안정한 Al fault로 주장하지 않는다.
위 discrepancy는 bulk fit만으로 slip/opening이 검증되지 않음을 보여준다.

At a=1.05h, s=0.17b, direct radius50/layers40:
- candidate0 energy difference≈5.56e-7 eV/cell;
- normal derivative difference≈1.37e-5.
이것은 한 설정의 independent comparison이며 완료된 다중 refinement가 아니다.

## 8. 결과 파일

results/aluminum_full_fcc_calibration/:
- calibration_targets.csv
- parameter_sets.csv
- normalized_residuals.csv
- identifiability_svd.csv
- parameter_correlations.csv
- heldout_validation.csv
- model_comparison.csv
- summary.json

results/fcc111_active_interface/:
- active_interface_energy_grid.csv
- gsf_curve.csv
- opening_curve.csv
- hessian.csv
- barrier_summary.csv
- convergence_summary.csv
- active_interface_contour.png
- gsf_validation.png
- opening_validation.png
- summary.json

현재 summary.json의 “successful but practically non-identifiable”는
위 structural issue를 반영하기 전의 draft status다.
결과표/그림 존재 자체가 과학적 validation 완료를 뜻하지 않는다.

## 9. 테스트 기록과 현재 미확인 사항

이전 진행 기록에 있는 실행 결과:
- 당시 new calibration/interface tests:11 passed;
- old/new full-FCC targeted combined:38 passed,20.68 s;
- reproducibility test 추가 후 calibration subset:6 passed,3.64 s;
- existing hybrid/Al/fullFCC registry subset:40 passed,95.06 s;
- app:27 passed,2 skipped,81.08 s;
- desktop smoke:PASS, LJ a0=0.7713438268704838,
  kappa=86.29296488740997.

이는 과거 실행 기록이며 이번 문서 작성 중 재실행한 결과가 아니다.
특히 잘못된 theoretical rank assertion을 통과한 테스트는
그 이론의 타당성을 증명하지 않는다.

문서 요청 직전에 full solver suite를 실행했으나 session15773의 최종 출력을
이번 복구 과정에서 확보하지 못했다. Session은 현재 조회 불가다.
따라서 **현재 변경에 대한 full solver 최종 결과는 미확인**이다.
추정 개수로 PASS라 하지 않는다. 다음 실행에서 로그를 남겨 다시 수행한다.

이 환경에서는 py -3 launcher가 사용 불가했던 기록이 있다.
실제 Python3 interpreter를 확인하여 동등한 -m pytest 명령을 사용한다.

## 10. 재개 순서

1. AGENTS.md와 이 문서 읽기.
2. 실제 worktree/HEAD/remote/status 확인; 사용자 파일 보호.
3. 새 calibration source, target mapping, tests, draft docs 읽기.
4. Section6의 structural rank와 overclaim부터 수정.
5. 필요한 deterministic calibration/identifiability/refinement를 재실행.
6. Static interface limits/derivatives/GSF/cleavage와 mapping을 검증.
7. Loaded barrier/spinodal이 없으면 구현하거나 미완료를 명시.
8. 데이터/문서/summary를 동일한 evidence 수준으로 일치시키기.
9. Targeted -> solver_v1 -> app -> desktop smoke -> diff check.
10. 최종 변경 검토 후에만 target branch 정상 commit/push, remote 재확인.

물리 mobility/seconds/Hz는 계속 unavailable.
기존 TwoRowLJ/reduced hybrid/PDE/time framework/specimen semantics는 보존한다.
사용자가 다음 지시에서 범위를 바꾸면 최신 지시를 우선한다.
