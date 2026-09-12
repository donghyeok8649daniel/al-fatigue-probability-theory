# CURRENT_WORK_HANDOFF.md — 단계별 검증 후 재개하기

## v30 완료 — 내부 timestep work 검증

- 실제main12개162ps+legacy6개12ps완료. session10268/90005 exit0.
  results/work_resolution_v30는final_analysis재생본;초기endpoint-only분석은
  이번작업생성파일만갱신했다. raw cache및v29기록은보존.
- Native계측은모든내부step의power합산. dt2.5/1.25/.625fs에서RMS잔차:
  normal .002454/.0005808/.0001521eV,slip .002029/.0005094/.0001229eV.
  약2차수렴. coarse25fs power오차는이추세없음. 최종endpoint만고르지않음.
- Native/densework차max6.59e-14eV. 새observer ON/OFF6개legacy대조는
  savedcoordinate/thermo exact동일. 힘/동역학은변경하지않았다.
- target16PASS2.69s,fullsolver706PASS1824.40s,app34PASS250.16s,smokePASS.
  실행중MD/tests없음. 생산physics/kineticJSON불변,physicalsecondsHzdisabled.
- 해결한것:reference외력work의저장간격quadrature오차. 미완료:약한외력
  장기/독립restart,선형loss/lowfrequency/localPMF/materialmobility검증.
  이174ps검사를장기보정으로부르지않는다. git최종SHA는실제log/ref확인.

### 아래는 진행 당시 기록

- 사용자 ㄱㄱ,시작/remote795a6de6b77707b67b5068d5f9c34d2617420f41 clean/fetch확인.
- 약한외력장기화전에v29 power잔차원인검증. Native ave/time으로모든내부step
  conjugate power합산후trapezoid endpoint보정;힘/생산PDE변경없음.
- run_work_resolution_v30,INTERNAL_STEP_WORK_V30 참조. .cache/work_audit_v30/main
  실제12records162ps 실행중session10268. 3dt/normal-slip25ps+dense2ps대조.
  장기weak-probe/독립restart아직미실행. 이단계는mobility보정아님.
- targeted16PASS2.73s. fullsolver26020,app/smoke32474실행중. 완료후업데이트.

## v29 완료 체크포인트 — 아래 중간 로그보다 우선

- 실제12개10ns MD 모두완료(session2533 exit0). 실행중 MD없음.
  analysis/report 모두exit0, 결과results/low_frequency_forcing_v29와v29문서참조.
- normal lag dt2.5/1.25fs .090405/.164361deg,slip .237363/.180647deg.
  loss/null normal .975/1.771,slip1.870/1.422. FDT sensitivity와양립하지만
  정밀선형mobility인증아님. normal2차고조파6.414/6.260null,
  slip3차2.022/1.454null:강한probe의비선형오염증거. half-force loss미분해.
- 실제control disk에zero drag포함:normal/slip M상한null. 점추정만생산에
  넣지않는다. 25fs power잔차max.061785eV,parts잔차.009527eV;
  이를소산으로해석금지. 발열first/last100ps차 .068–5.024K.
- target27PASS3.79s,fullsolver700PASS1993.89s,app34PASS154.50s,smokePASS.
  모든결과/plot검사완료;git최종commit/push는실제log/ref확인.
  생산LJ/Bessel/static/PDE/MaMs/kineticJSON/UI/Ac불변. physicalHzdisabled.
- 다음:작은probe/긴record/독립restart와내부step work누적검증,
  lowfrequency/localPMF/material검증. 제안이지실행완료아님.

## v29 저주파 추가 MD — 과거 진행 로그

- 사용자 `해봐 임마`에따라실제새MD실행. 시작/원격5c519206ec6c32dd58014cf625cd96317b449102,
  fetch성공/clean/대상branch/worktrees확인. 기존이론/생산clock불변.
- `LOW_FREQUENCY_FORCING_V29.md`,run_low_frequency_forcing_v29.py.
  .cache/low_frequency_forcing_v29/main/protocol.json 먼저저장.
- .2cycles/ps,normal/slip±4RMS1000ps,±2RMS500ps,±4RMSdt/2 1000ps:
  총12개10ns. 첫100ps제외,25fs저장,기존restart/Al99reference/NVE유지.
  productionHz/이동도변경아님. half-force는같은loss정밀도가아닌진폭대조군.
- 실제MD session2533 실행중,첫normal±4RMS 로그에서25ps진행확인.
  완료전보정성공선언금지. 원본/partial덮어쓰기금지;summary검증완료case만재사용.
- 완료후 analyze --study .cache/low_frequency_forcing_v29/main --source
  .cache/modal_kinetic_v25/nve_N6_serial_5000ps_dt2p5 --out results/low_frequency_forcing_v29.
- 코드/문서미커밋. tests/전체회귀/actual분석은진행후갱신할것.
- 중간확인:normal4RMS±1000ps완료. pairedchi=.002367889772856313
  -3.7362257420419344e-6i A²/eV,phase-.0904053deg. 아직진폭/dt미완료이며
  최종보정인증아님. slip4RMS±1000ps가약500ps진행중.
  targeted27PASS1.92s,app34PASS154.50s,smokePASS. fullsolver53856실행중.
- 후속중간확인:fullsolver700PASS1993.89s(session53856 exit0).
  normalbase900ps관측null loss3.831129e-6 >pairedloss3.736226e-6:ratio.9752,
  현재normal은아직unresolved. slip4RMS±1000ps도완료,pairedchi
  .012192567670863129-5.051132031263749e-5i,phase-.237363deg.
  현재4/12완료,2RMS500ps진폭대조군실행중. dt/2미완료.
- 이후6/12완료(normal2RMS±500ps포함). slip2RMS진행후dt/2네개가남음.
  실행부/분석부는분리되어있고session2533유지. 현재report_low_frequency_forcing_v29.py
  추가됨:전체분석완료후 --results로진폭/dt/FDT/이동도오차범위와그림생성.
- 현재8/12완료,normaldt1.25fs±4RMS1000ps진행중;이후slipdt1.25fs두개.
  생산code/calibrationJSON은불변. 분석은whole/steady work와phasework의
  endpoint항등식을분리하며q사인적분일치를독립물리검증이라고하지않는다.

## v28 이동도 보정 — 연구 후보/검증 완료, 생산 보정 gate 미통과

- 사용자 `보정ㄱ`; 시작/원격8d222ea073c09461e72961e10be00813e1be33b9,
  fetch성공/clean/branch/worktree확인. main과모든다른worktree보존.
- `IMPEDANCE_MOBILITY_CALIBRATION_V28.md`,results/impedance_calibration_v28.
  기존5ns dt2.5/5fs spectra와완료v26/v27직접응답으로새WLS실행.
  새MDtrajectory는실행하지않음. 기존원본/결과불변.
- scalar저주파dragfit:normalM7.246372352e10,slip2.180655810e11 m²/(Js),
  ratio3.00930687. fitrange .02–.16cycles/ps;더낮은4band는loss제외.
  784관측/168controlfit;독립표본아님. controlfitrange normal6.523–8.846e10,
  slip1.859–2.682e11. 낮은band의일부오차큼:zero-frequency미인증.
- 복소disk를정확히1/chi로사상하여Gamma=Im(1/chi)/omega범위추정.
  null+dt/amp/block/sign범위에는fullforce4개모두zero drag포함:
  M상한null. 이는CI/엄밀확률bound아님. pointM만보고완료금지.
- 같은PMF상수overdamped저주파fit의고주파complex오차~8–10observednull.
  storage가static보다커서순수overdamped응답불가. 관성/메모리reference의
  상수drag까지기각하는것은아님. 생산관성도입/에너지변경금지.
- actual5ns두rawNPZ로별도scalarPSD재구성PASS,maxrelative2.22e-16/4.44e-16.
  arithmetic일치이지물리정확도아님. plot렌더/시각검사완료.
- 최종targeted20PASS .87s(새10tests포함),fullsolver693PASS769.26s
  (session72196 exit0),app34PASS143.03s,smokePASS. 실행중MD/tests없음.
  196개분해불가sourceband도원래사유와함께별도보존(172bandwidth/24bin).
  최종재생은.cache/impedance_v28_final_reproduction;runner의5개결과를
  이번turn생성results로반영. 원본/source결과는덮어쓰지않음.
- 저주파.2cycles/ps의SNR3 사전계획은1RMSforce에서normal13.235ns/
  slip7.310ns. force4배로시간1/16은선형가정의계획일뿐미실행/선형성미검증.
  필요한평균외력일은2kBT*SNR²라서단순force확대가발열문제를없애지않음.
- 생산LJ/Bessel/static/PDE/MaMs/kineticJSON/UI/A_c변경없음.초Hzdisabled.
- 최종diff/정상commit/push와remote SHA는실제git log/ref로확인.
  main참조80cacb4/originmain c43d8e0 보존. 아래v27는이전완료상태.

## v27 완료 체크포인트 — 아래 진행 로그보다 우선

- 시작/원격9afdc390239ae86b25ca5165d279ab5e723d569f, clean 확인.
  `PHASE_RESOLUTION_VALIDATION_V27.md`와 `results/phase_resolution_v27` 확인.
- 실제12×400ps reference MD 완료(session41642 exit0). 분석21355 exit0,
  그림시각검사 및 별도harmonic_controls 완료. 실행중MD없음.
- fullforce paired lag:normal2cycles/ps .63616→.88445deg(dt2.5→1.25fs),
  slip1cycles/ps1.36669→1.46510deg. 관측null envelope대비loss normal
  1.283/1.783배,slip1.748/1.872배. 유한주파수위상검출근거확보.
- **고정밀위상수렴/M보정완료아님**:normalIm dt변화39%,slip7.1%; thermal과
  순수dt오차분리못함. 다수180ps/half-force slip은null이하. v26저주파.05
  미분해유지. Observed null은CI아님. shared restart/planes독립으로세지않음.
- FDT spectral sensitivity와thermal envelope안에서양립. Storage/static
  ratio1.12079/1.15858:이band를같은PMF의상수과감쇠M로바로이전불가.
- covariance25/50fs trap/Simpson감사별도저장. coarse50fsnormalSimpson의
  허수부부호실패보존;정답으로선택하지않음. 소산예측은직접PSD로비교.
- raw2ndharmonic최대2.41%normal/9.98%slip. 무하중2f/3f와forceparity
  대조에서모든예상even2/odd3성분은관측null이하(maxratio.887).
  정확선형성증명/비선형성단정금지. 모든개별고조파도보존.
- clock최대5.68e-14ps,COM1.554e-14A/ps. 전체work잔차최대
  .009072eV(parts)/.008294eV(power). 평균T299.331–302.262K;
  first/last40ps T상승.469–4.542K. 위상적분/endpoint/내부에너지분리.
- 최종Target25 PASS1.65s,solver683 PASS1790.49s,app34 PASS70.04s(no skip),
  smokePASS. 신규syntheticengine2freq×3dt,legacyexact재현,protocol재현,
  incomplete결과생성거부 모두실행PASS. 최종git/diff는아래/실제log확인.
- 결과별pair재구성,work+endpoint항등식,독립고조파fit일치,production clock
  false재검사PASS. 최종diff검사PASS;commit/remote SHA는실제gitlog확인.
- 생산LJ/Bessel/static/PDE/Ma/Ms/kineticJSON/UI/A_c불변. 실제생산초Hz
  disabled. 원본캐시보존. 다음은저주파정밀도/동일local좌표PMF/Markov극한,
  materialkinetics 검증이며 이번유한주파수검출을전체mobility보정으로승격금지.

## v27 위상지연 검증 — 진행 중

- 사용자 요청: 위상지연 검증. 시작/원격 HEAD9afdc390, clean 확인.
- `PHASE_RESOLUTION_VALIDATION_V27.md`와 새 phase_resolution/study 모듈 참조.
- 무하중5ns PSD로 사전 선택: normal2/direct1101 cycles/ps. ±.5/1RMSforce,
  dt2.5fs 및fullforce± dt1.25fs, 각400ps 총12개. 생산주파수보정아님.
- 원본 `.cache/phase_resolution_v27/main`, protocol.json 먼저저장.
  MD session41642 실행중. 아직위상검증통과선언금지.
- fullsolver60933 완료:683 PASS1790.49s. 새엔진benchmark f1/f2 각3dt 완료,
  legacy기본benchmark exact재현PASS. 분석스크립트는실제MD완료후검사할것.
- 수직6개 완료. fullforce paired phase기본약-.636deg,dt/2약-.884deg.
  각각관측nullenvelope대비loss약1.28/1.78배. 이는고정밀M보정아님.
  slip절반진폭±완료,fullforce및dt/2진행중. 중간값으로전체PASS금지.
  전체진폭normal work약.75–1.14eV,parts에너지잔차.00014–.0029eV,
  first/last40ps평균T상승3–4.5K. 조건부해석필수.
- 시간간격변경후실제로그Time=0초기화확인;analyze는모든thermo행의
  engine time=step*dt 및끝400ps를검증. workendpoint와위상에서유도한
  periodicwork는분리한다. 동일q사인적분의일치는독립물리검증이아님.
- targeted25 PASS .75s. app34 PASS70.04s, smokePASS(session81942).
- source null360ps old.05 signal floor normal2.967e-5/slip1.660e-4 A²/eV;
  v26예상소산보다큼. observed envelope이지CI/독립replica아님.
- 완료후 analyze --study .cache/phase_resolution_v27/main --source
  .cache/modal_kinetic_v25/nve_N6_serial_5000ps_dt2p5 --out results/phase_resolution_v27.
  dt/amp/FDT/temperature/work잔차까지실제확인하고완료내용갱신.
- 코드/문서만현재미커밋. 생산LJ/Bessel/PDE/Ma/Ms/kineticJSON/UI불변.

## v26 완료 체크포인트 — 같은 좌표의 직접 강제응답 검증

### 최종 상태 (아래 실행중 로그보다 우선)

- 8개400ps 모두완료(session68476 exit0). 분석과시각검사완료:
  `results/collective_forcing_v26`. 전체원본 `.cache/collective_forcing_v26/main`.
- 수직 Rechi FDT오차.04–.25%, slip.13–1.19%; amplitude두배정규화응답
  변화.21%/1.06%. 동일좌표conjugate정규화/근선형응답은지지됨.
- 모든paired Imchi blockrange가0을가로지름. 음의마찰/상수M보정아님.
  최대work잔차.012083eV(power)/.005743eV(parts),864atoms전체.
  예상작은소산~.001–.002eV보다커서소산보정통과금지.
- 별도syntheticengineoscillator phase오차1.12e-10rad(dt.0025);
  3dt검증완료. Almass/drag보정값아님. grossclock/sign오류설명은지지되지않음.
- 계획SNR3기록길이normal156ns/slip97ns(.05cycles/ps,1RMSforce):
  stationarity/weakdrive/spectralnoise가정의계획추정이지CI/생산t0아님.
- 무외력새runner가이전5frames q/thermo정확재현. 생산LJ/Bessel/static/
  PDE/Ma/Ms/kineticJSON/UI/A_c불변. 초Hz여전히disabled.
- 실제테스트18targeted PASS.74s,full676 PASS968.59s(session79068),
  app34 PASS101.49s(no skip),smokePASS. 과학적Imchi분해능은미통과.
- 수정핵심은누락된강제응답검증경로추가. 전체mobility보정완료라고말하지않음.
- 시작fbaf4a2. 최종git상태는log/remote확인;아래는중간보존로그.

사용자 "고쳐보자 그럼"에 따라 v25의 누락된 independent response 검증을
추가했다. 임의 M비율수정/생산초Hz활성화가목표아님. 시작HEAD fbaf4a2,
origin동일/clean확인. `CONJUGATE_RESPONSE_VALIDATION_V26.md` 먼저읽기.

- `collective_forcing.py`: q=(u1mean-u0mean)·e, ±F/Np로정확히conjugate
  force/zero-net-force. complexlockin와canonicalFDT endpoint항유도.
- reference MD runner에optionaldrive만추가;기존none경로와생산PDE불변.
  실제LAMMPS .1ps smoke완료. 강제응답/정규화테스트3개포함focused18 PASS.74s.
- `run_collective_forcing_study.py`: 같은N6restart,±.5/1thermalRMSforce,
  normal/direct110총8개400ps,.05cycles/ps,dt2.5fs. 실제sourceMD단위;
  물리피로Hz아님. session68476에서최대2개동시실행중.
  ignoredcache `.cache/collective_forcing_v26/main`에protocol/원본보존.
- 완료후 analyze subcommand --study 위폴더 --source
  `.cache/modal_kinetic_v25/nve_N6_serial_5000ps_dt2p5` --out
  `results/collective_forcing_v26` 실행. 모든부호/amp/block/FDTCutoff,
  cross응답/COM/work-energy표검토. 아직결과완료/일치라고말하지말것.
- fullsolver session79068실행중(676예상;actual출력확인).
  app session86348완료34 PASS101.49s,smokePASS. no skip.
- M(omega)판정은smallImchi의실제오차중요. 실수응답일치만으로M통과금지.
  기존v25matrixM와single-coordinateinverseimpedance같지않음.
- 아직v26commit/push없음. 위세션끝난뒤실제결과/테스트/인계갱신,
  diffcheck/fetch후normalcommit/push. main/kineticJSON/production불변.

## v25 계산 완료 — 실제 모드 감쇠 / 이동도와 적용 한계

### 최종 과학 체크포인트 (아래 진행 로그보다 우선)

- 7 main records(5개1ns + 2개5ns), 별도thermal1ns/sampling250ps 모두완료.
  장기run50958와분석64535도exit0. 실행중MD없음.
- `all_control_comparison`: 7case,360 modal fit 모두성공/heldoutOU보다우수,
  activebound0. `control_comparison`은5개1ns비교로보존.
- `long_record_comparison`: 1/2/5ns, dt2.5/5fs, NW2/3, full/half/quarter
  전체표+시각검사된그림. 5ns초기3band M: normal6.70–8.05e10,
  direct1101.80–2.39e11 m²/(J s), **144atom plane-gap PMF**.
  최저band는상수plateau불충분: slip dt차이최대51.75%,최저band30.45%.
  이범위는CI아니며생산Ma/Ms아님. 기존576atomplane와좌표정규화다름.
- dt2.5/5fs5ns 에너지범위2.49533e-5/1.98929e-4eV/atom;
  T299.068/298.126K. 좁은phonon폭은~.120THz로잘일치하지만
  저주파M는그만큼일치하지않는다. 실험~.272/DFT~.186THz와별도비교.
- 1ns와5ns prefix,열관측추가trajectory 정확일치. 5fs/25fs sampling
  차이최대0.0241%. 같은기록이므로독립replica라고세지않음.
- 최종source-hash-boundleakage재실행이기존10records와정확히일치.
  입력은spectrum_lowlimit(5bands),초기3band폴더아님.
- source8192 DHO추정,직접spectralM,thermal19–41%결과/질량수정/
  에너지와열항동시정규화는MODAL_KINETIC_CALIBRATION_V25.md에유도.
- 최신테스트 targeted12 PASS1.42s,solver673 PASS906.20s,
  app34 PASS97.53s(skip0),desktop smoke PASS. 생산LJ/Bessel/static/PDE,
  A_c및kineticJSON불변. 물리Hz활성화는하지않음.
- 다음과학문제: 같은collective coordinate/PMF의독립소신호강제응답 또는
  올바르게projected constrained-force memory검증. 단순mass/phononlifetime,
  Np배율,실험linewidth배율을생산M로대입금지. 아직이실험안했음.
- Git v25시작4520d0a…;최종commit/push는git log/원격으로확인.
  아래는실행중상태의보존로그이며최종상태가아님.

### 2026-09-12 04:45 KST 확인된 최종 코드 테스트 / 남은 계산

- mass-aware 최종코드 full solver: **673 passed in 906.20s**,
  `.cache/modal_kinetic_v25/solver_regression_mass_checked.log`.
- app **34 passed in 97.53s**, skip0; desktop smoke PASS
  (`app_regression_mass_checked.log`, `smoke_mass_checked.log`).
- 최신 targeted12 PASS1.42s. 생산/static/physical kinetic JSON 불변.
- dt5 5ns 및 prefix1/2/5ns/mode 분석 완료(session24403 exit0).
- dt2.5 5ns session50958은 약4.4ns 진행. 분석pipeline64535 대기.
  완료 후 `run_long_record_comparison`과 7case control summary를 실행할 것.
- 250ps 5fs sampling control 완료: 기존25fs기록 prefix와 정확히 같음.
  동등225ps의5fs/25fs PSD 추정 M 대각성분 최대차이0.0241065%.
  `sampling_interval_check/summary.json`. sampling alias가DHO차이를설명하지않음.
- 새 `run_long_record_comparison.py`는 두dt별1/2/5ns prefix와전체band,
  half/quarter sensitivity를저장한다. freshoutput필수.
- full regression 끝났으므로 소스수정없이 결과/문서정리후diffcheck/최종fetch,
  normalcommit/push가능. **아직v25 commit/push안함**. main불변.

### 추가 실행/수학 감사 — 이 절이 아래 진행상태보다 최신

- 다섯1ns MD/분석 모두 완료. `control_comparison`과 `spectral_comparison`
  (그림시각검사완료)이최종1ns비교. N12 NVE dt5→2.5fs 에너지범위
  3.3871e-5→7.6155e-6eV/atom(약4.45배감소). MD온도는각각303.14/304.31K:
  동일300K라고속이지않는다. NVT~299.96K, N6~299.07K,
  lattice4.05Å N6~305.57K/압력9232bar(고정격자비교,무응력Al아님).
- source DHO보다직접저주파M는모든5cases에서작다. modal/half/cutoff
  fit288개모두성공·heldoutOU보다작음·activebound0. 이것은M보정성공아님.
- N6 5ns dt2.5 session50958/분석대기64535는계속진행.
  추가N6 5ns dt5 session60124: `nve_N6_serial_5000ps_dt5` 진행.
  dt5용후속analysis/prefix pipeline은아직설정안함. 두5ns를1/2/5ns
  prefix 및최저.001–.002cycles/ps까지비교해야함. 1ns동일restart부터
  다시시작한nested record라independent replica로세지않는다.
- thermal1ns N6 dt2.5 session83076 완료, analysis69792완료.
  `thermal_channel_mass_checked`가최종열적정규화;기존thermal_channel_audit보존.
  수직저주파power의19.46–40.60%를이웃2평면kinetic energy로temporal
  heldout예측. time-shift대조군은모두음수, slip은일관된효과없음.
  열확산인과증명/전체기억원인확정아님:PE/energy current/전자열수송미측정.
- 실제엔진으로확인한질량override:restart26.982, mass1 26.98후26.98,
  Al99 pair_coeff후26.982. 초기kinetic postprocess만26.98이었음.
  run_reference_thermostat_md가이제실제엔진mass조회/기록, helper에명시전달.
  기존원본thermal배열은보존하고analysis에26.982/26.98 보정명시.
  KE오차7.406e-5→약6.3e-8, M(C0,kBT,tau)와MD forces/timestamps영향없음.
- 정규화수식중요: Fbar=Fplane/Np, Mbar=Np Mplane은동일q에서
  **kBTbar=kBT/Np까지같이해야**같은확률generator. 물리T변경아님.
  M만Np배하고열항유지하면diffusion Np배오류. productionlocalcelltransfer
  가정과순수대수정규화를구별. 이관계회귀테스트추가; sourceJSONmetadata도명확화.
- Glensk동일그림DFT MDsquare도분리추출:300K .186069THz,
  900K .816578THz;experimentalgray와별도. `phonon_reference_with_dft`.
  물리초/Hz default는여전히disabled. LJ/Bessel/static/PDE/A_c불변.
- 첫full671 PASS1204.43s. 후속열적테스트추가후22focused PASS2.16s,
  massfix후newtests12 PASS1.58s. 두번째full session87774/로그
  solver_regression_final.log 진행;이실행도massfix전import일수있으므로
  최신source/helper에대해마지막full재실행권장. App34 PASS117.21s,
  smokePASS1.939s. 모든최종변경후targeted/full/app/smoke/diff확인후만commit.
- freshfetch재확인local/origin모두4520d0a…;main80cacb…/originmainc43d8e…불변.
  아직commit/push없음. 불필요한sourceMD새실행보다현재긴두기록완료/정리우선.

### 2026-09-12 약04:00 KST 추가 인계(계산 진행 중)

최신 사용자 추가 요청은 이동도 계산. 실제 SI 이동도까지 계산했다.
아래 초기3개 MD 안내보다 이 절/실제summary 상태를 우선한다.
- 완료1ns: N12 NVE dt5fs, N12 NVT dt5fs, N6 NVE dt2.5fs,
  N6 lattice4.05Å NVE dt2.5fs. 모두 `.cache/modal_kinetic_v25`에
  complete summary/plane_coordinates 보존. N12 NVE dt2.5fs는 약0.8ns 진행.
- 추가 N6 plain serial5ns dt2.5fs 실행 중:
  `nve_N6_serial_5000ps_dt2p5`, session50958. 같은 N6 equilibrium.restart를
  사용하므로 새독립replica가 아니라 nested reproducibility/record-length 연구다.
  완료 후1/2/5ns prefix 비교 및최저.001–.002cycles/ps band를 검증한다.
- 독립1ns 분석 pipeline session98363은 마지막 N12 dt2.5 완료 대기;
  lower-limit/mode-spectrum pipeline session22820도 같은 마지막case대기.
- full solver regression session94204 진행(64%까지출력), 로그
  `.cache/modal_kinetic_v25/solver_regression.log`. App34 PASS117.21s,
  smokePASS1.939s; focused20 PASS4.27s. 아직full통과라고말하지말것.
- 추가 low_frequency_mobility.py는 DPSS two-sided PSD, K_band=S_band/2,
  3x3 friction/mobility를 계산. Parseval/phase-leakage tests2개통과.
  기본3bands, 명시followuplower2bands, longrecord추가2bands를구별.
  band<taper bandwidth 또는3bin미만은invalid기록, 인증아님.
- run_mode_spectral_mobility.py는 공간mode별상수M가설 검사;
  run_spectral_leakage_control.py는순수진동선 leakage의worstphase2x2eig.
  첫1ns NVE leakage/실측power 최대6.898e-5: 이 leakage만으로 설명불가.
- source8192 DHO plane PMF M_a5.73–8.14e10,M_s1.17–1.68e11 SI.
  독립N12 NVE1ns 직접저주파는 M_a1.71–1.83e10,M_s4.42–5.53e10.
  DHO외삽 M_a7.06–7.27e10,M_s1.30–1.52e11와불일치. window선택금지.
  더낮은band포함normal1.47–1.87e10,slip4.42–7.47e10.
  Np144/576 extensivity환산normal약1e13,slip약3e13이나가설단계.
  공간mode별 M도수배차이: local평균값만으로constantM승인금지.
- Glensk2019 원논문Fig2b TL295K Gamma~.271619THz(g=pi Gamma),
  출처/벡터표식/hashes results experimental_linewidth에저장. N12 NVE
  ~303K TL Gamma .116–.120THz. 실험오차없는정밀보정이라고말하지않는다.
  Tang2010 PRB82,184301 PDF MD5abfe242aaeed0b72b086462c7dbef962 확인,
  5pages읽음/Fig2시각검사: 옛실험small-q선폭instrument floor경고.
  신규수치target로추출하지않았음. 설명은kinetic sources문서에추가.
- 마지막 완료되면 control_comparison (run_kinetic_control_summary,
  --sources-root .cache/modal_kinetic_v25로pressure도hashbound계산),
  spectral comparison/그림 생성,시각검사,문서최종수치/Git검증.
  `interim_control_comparison`은먼저끝난3cases중간결과이며최종아님.
  first_half_fit 및모든현재유용결과보존. 아직v25 commit/push없음.

사용자 요청: 보정값을 실제로 찾아낼 것, 아침까지 약6시간 작업.
시작/fresh remote4520d0a065815966461eefc7d117c217a13c5840,
probability-pde-solver-v1. v24는 완료/push됐으며 새 v25는 아직 미완료.
기존 source8192frame을 4096 training/4096 temporal heldout으로 나누고,
6개 독립k·3축·3lag cutoff의54개 DHO position-correlation fit을 실제 실행했다.
54개 모두 heldout RMSE가 단순 overdamped exponential보다 작다.
median .0752775 vs .3589432. Mode damping은 .020919..1.182374 /ps;
이것은 source modal parameter이지 production a/s mobility가 아니다.
Envelope time 1/g와 formal integral time 2g/(g²+wd²)를 혼동하지 않는다.
출력 results/modal_kinetic_calibration_v25/first_half_fit 및 mode_kinetic_calibration.py.
초기 targeted14 PASS1.21s; 이후 MD protocol test 추가, 전체회귀 아직 안함.

현재 독립 MD를 실제 실행 중(완료여부는 cache summary/checkpoint로 확인):
`.cache/modal_kinetic_v25/nve_1000ps_dt5`, `nvt_1000ps_dt5`,
`nve_1000ps_dt2p5`. 동일 nve_benchmark/equilibrated.restart 출발,
6912atoms,4.065Å,sourceAl99,nominal300K,1ns sampling,25fs frame.
원본 published trajectory 재현이 아니라 independently generated reference다.
Benchmark25ps equil +2ps NVE 실제완료17.45s,meanT303.1566K.
NVE/NVT/dt 차이와 실제온도를 비교해야 한다. 결과를 미리 PASS라 하지 않는다.
실행기 run_reference_thermostat_md.py. LAMMPS는 target-only, LJ/Bessel/PDE 불변.

도구는 시스템설치 없이 ignoredcache에 추출했다. PyPI wheel은 MS-MPI DLL없어
실행실패(남겨둠). 공식 serial22Jul2025update4 installer SHA256
47f1aeb0fcbeadcc211c9b1c258152c3545464e3381975a681bf558891a06986
배포 SHA256SUMS와 일치, 7zip MSI administrative extraction 후 NSIS파일만 추출.
Installer 실행/전역PATH/registry변경 없음. Python module은 캐시
lammps_serial/Python, LAMMPSDLLPATH는 lammps_serial/bin을 실행별로 지정.
실제 engine20250722/eam-alloy-omp 동작 확인. 개인절대경로는 commit하지 않는다.
아직 productionMa/Ms/t0없으며, 먼저 thermostat/time/coordinate projection 검증.
최신코드/결과 보존 후 targeted/full/app/smoke/diffcheck를 실제 실행하고만 commit.

## v24 kinetic validation — 2026-09-12, 실제 긴 궤적 검증 완료

최신 요청은 실제 항복/피로/production 초·Hz 검증을 계속하라는 것.
이번에는 kinetic 분기를 실제 실행했다. 아래 v23 재료/core 기록은 보존한다.
시작 local/fresh remote는 d222296c3fad369b74b6f4794761c427b16ee8d5,
probability-pde-solver-v1, clean. 실제 작업 위치는 기존 aft-pde-bessel-38969ad.
권한 변경으로 중간에 중단됐으나 완료 여부를 파일/도구로 확인해 이어갔다.
새 code는 검증 실행기와 protocol tests뿐; energy/PDE/UI/default/M/A_c 불변.

- 공개 Zenodo10014454 실제 Al99 300K NVT 자료8192frame/204.775ps 확보 완료.
  약1372.4MiB compressed prefix를 스트리밍했으며 전체3.5GB checksum 검증 아님.
  좌표만 ignored cache `.cache/kinetic_validation_v24/zenodo_8192`에 저장.
  Projection SHA256 0d4eead55eff67a4ec075152501db6c7cf39032cca04e34b32ce6c13ab7f3cb3.
  기존1024frame의 좌표와 시간은 새prefix와 exact equality 확인.
- 14개 predeclared prefix/cutoff/quarter 검사를 실제 실행. 동일3.15ps cutoff의
  floor는 .0382522ps→.00796563ps로 감소하지만 minimum integral eigenvalue는
  -.00719449ps. 8192frame에서 cutoff6.35ps는 모두양수지만 floor보다작고,
  25.55ps에는 다시음수. 일부양수cutoff를 골라 M을 만드는 것은 금지.
- 0.15ps lag correlation minimum -.569963, sensitivity .182309, 비3.126.
  짧은시간 direct harmonic overdamped fit은 이 좌표에서 맞지 않는다.
  이 비는 confidence/sigma가 아님. 충분히 느린 Markov limit의 부정도 아님.
- 같은 source/box/atomcount의 static harmonic covariance와 MD variance차이는
  normal -1.816%, slip +5.348%, transverse +9.568%; local four-block sensitivity
  이내. 개별mode 모두일치/비조화성없음/kinetics인증이라는 뜻은 아니다.
- 실제 Gorman1969 원본표식 재검사:7fit/3heldout, line velocity/shear
  1.1627233044e-5m/(Pa s), heldoutRMS3.501185m/s. a/s 이동도로 복사하지 않는다.
- 상세 수식/범위/수치는 solver_v1/KINETIC_RECORD_VALIDATION_V24.md 및
  results/kinetic_validation_v24. 이 단계에서 실제 항복/피로예측은 새로
  실행하지 않았다. 기존 finite-source strain budget와 피로관측량 mismatch를
  해결했다고 보고하지 않는다. Ma_phys/Ms_phys/t0/productionHz 불가 상태 유지.
- targeted42 PASS2.27s; full solver661 PASS1509.98s; app34 PASS128.91s,
  skip0; smoke PASS1.483s. 테스트6개 추가는 protocol 검증이고 실험인증이 아님.
  전체회귀는 코드변경 완료 후 실제 실행했다. 중단직전 record analysis는 완료,
  source-normalization 명령은 권한전환으로 실행 전 거부되어 새로 실행해 완료.
- Windows apply_patch.bat 인수처리 첫실패는 파일변경 없었고 native apply_patch로
  재시도했다. Git commit/push는 이 인계 뒤 fresh fetch/diffcheck 후 실시한다.

다음: 같은 적분을 편리한창으로 재피팅하지 말고 mode별 memory/thermostat와
production coordinate/PMF projection을 검증한다. 항복은 finite source 방출/
상호작용/실측 source population을 포함한 동일 proof-strain 관측량이 필요하다.
단순 GPa traction rescaling, 임의 line length/Ac/t0로 성공을 만들지 않는다.

## v23 최종 인계 — 2026-09-11, 실제 계산과 전체 회귀 완료

**이 절이 아래 모든 v23 진행/재시도중 기록을 대체한다.** 아래 기록은
실패 원인과 계산 경과를 보존한 이력이며 현재 실행 상태가 아니다.
사용자 최신 요청은 "계속해". 원래 LJ/Bessel 이론을 보존하며 v22의
source-core/계면 보정 충돌을 진단하고, 같은 함수족 후보를 실제로 풀었다.
Production/PDE/UI/default/응력 bridge/mobility/A_c는 바꾸지 않았다.

### Git / 재개 위치

- 시작 local과 fresh origin: `d1ac6695518d3738b8303687d8c6a7b58877203e`, clean.
- 실제 worktree: `aft-pde-bessel-38969ad`, branch `probability-pde-solver-v1`.
  OneDrive 원래 경로는 migration stub이다. 다시 preflight 후 재개한다.
- 이 절은 검증을 마친 커밋 직전에 작성한다. 이 인계를 포함하는 SHA는
  `git log -1 -- CURRENT_WORK_HANDOFF.md`로 확인한다. Push는 fresh fetch 후
  정상 fast-forward만 한다. 최종 도구 결과/사용자 보고가 실제 push 근거다.
- 기존 original8736cb9d…, v22wider2fd94b89…, targetAl9960c8a085…의
  보호 SHA256 세 개를 실제 재검증했다. Full binding은 validation_ledger에 있다.

### 이번 단계의 가장 중요한 결과

1. **Hxx>0은 slip saddle 부재의 증거가 아니다.** 실제 full3x3 Morse와
   연결성을 검사했다. v22wider 자기 saddle은 Hxx=+.068616인데도
   lambda_min=-1.387798eV/L0²이며 서로 다른 minimum으로 연결된다.
   Source/wider/newradial 9개 stationary state,3개 saddle connectivity,
   15개 signed static mixed-traction state를 실제로 확인했다.
2. 고정shape LP48개/독립primal-dual 검증: 두 offendingjet만은 동시에
   맞출 수 있다. 전체115jet minimax는 wider18.46445, sign-free18.33840.
   즉 가중치만 손보면 전체를 맞춘다는 근거는 없고, whole-family 불가능
   증명도 아니다. 같은 objective 두box QP는 otherloss41943.92로 악화된다.
3. 기존shape160 deterministic profiles(maxfev, optimum아님) 이후 새joint
   frozen-source force RMS는 .0995689eV/L0(R8), .0969927(R16excluded).
   Original .288834보다65.53%, wider .152417보다34.67% 작다. 그러나
   other-interface squared loss4520.30(wider2276.47)로 악화했다.
4. 새 radial 후보는 실제 stable vector core를 가진다:

   | R/ring |최대 force eV/L0|최소 Hessian eV/L0²|실제 시간 s|
   |---|---:|---:|---:|
   |8/5|4.343e-10|.172653874|487.823|
   |8/7|3.297e-8|.172236076|259.084|
   |10/7|2.771e-12|.130868746|461.814|

   모두 winding1과 FD곡률/replay 검사. R2 내부 변위차는 환경ring5→7에서
   2.1912e-5L0, 자유영역8→10에서1.9930e-3L0. Infinite-domain 인증 아님.
   Sampled quarter-span은 .885664/.885690/.880949L0, sourceR8은2.760327.
   이는 격자 profile 진단이며 물리 partial separation이나fit타깃이 아니다.
5. **새 단파장 원인 검사:** small-q matrix오차~.63%인데 q=(.5,.5,0)의
   한 vector curvature는 source60.0863에 비해375.8720eV/L0². 실제
   4개FCC phase-class 원자를 변위시키고 per-atom energy를 재계산했다.
   변위 .002→7.8125e-6L0 refinement에서 해석값으로 접근한다(최종오차
   4.86e-6relative). 단순 축/단위/행렬 prefactor 버그로 설명되지 않는다.
   Harmonic 양의고유값/장파장탄성 일치가 단파장 정확성을 보장하지 않는다.
   이것 하나로 nonlinearcore 오차의 유일한 원인이라고 단정하지 않는다.
   전체matrix470.8%오차와 screw projection -13.13%오차는 다른 지표다.
6. 기존 rational-quartic를 명시적5beta ablation으로 재검사했으나115jet
   bound개선없음, D3/K3상관.999929까지 악화. 채택하지 않는다.
   현재22채널 core가 미지원 nonzero quartic saturation을 조용히
   빠뜨리지 않도록 ValueError guard를 추가했다. 기본0은 완전 불변이다.

### 아직 채택하지 않는 정량적 이유

| J/m² |Al99 target-only|v22 wider|v23 radialjoint|
|---|---:|---:|---:|
|relaxed fault|.150479|.079457|.075309|
|adjacent saddle|.172002|.119750|.114757|
|finite W(40h)|1.741285|1.840583|2.524404|

안장점은 존재하지만 registry/개구/단파장/코어 형태를 함께 맞추지 못했다.
새minimum이나force개선을 Al 보정성공으로 승격하지 않는다. Al99는 같은
0K조건의 atomistic model target이지 모든 점의 실험적 정답은 아니다.
첫 fixed-registry normal peak12.970/8.924/6.537GPa는 idealcoherent
traction이지 실제항복이 아니다. Static ±50MPa normal, ±4MPa shear return은
동적hold나시편잔류소성검증이 아니다. Physicala/s M_a,M_s,t0/초/Hz,
유한loop activation, 실제yield/fatigue 및 A_c 보정은 여전히 미검증이다.

### 실패 / 수치 해상도 / 검증 기록

- 두box 후보 core는 root600s/descent900s 예산 종료; checkpoint 별도실제
  재평가force .242865/1.863155. 수렴minimum/unbounded energy증명 아님.
- 새radial 초기root는 force작지만 lambda_min=-.380259인 불안정정지점.
  L-BFGS 시도는 reciprocal series 실패로 중단했지만 같은energy Newton-CG는
  위minimum을 찾았다. 급수실패를 물리spinodal로 오해하지 않는다.
- Frozen-source ring7→9 force차 RMS5.864e-6, 잔여오차.099568eV/L0.
  Reciprocaltol2e-12→2e-14 변화0(표시정밀도). 모든비선형tail인증은 아님.
- `validation_ledger`: 7개core기록,9개 실패기록(launcher/preflight/report
  오류 포함), 미확인livecore0. Test/입력/result hashes 실제검사.
- Targeted48 PASS99.99s + topology2 PASS1.40s. 새 tests 총17개.
- **Full solver655 PASS1945.86s(32분25초), app34 PASS167.68s, skip0.**
- 실제 desktop smoke PASS1.523s: historicala0=.7713438268704838,
  kappa=86.29296488740997. 첫중단full은PASS가아니며별도기록했다.
- 최종 diffcheck/staging/freshfetch/commit/push는 이 기록 이후 도구로
  확인한다. 테스트를 다시했다고 문서만으로 주장하지 않는다.

### 다음 단계 (다시 동일 실패 fit을 반복하지 말 것)

1. `CORE_INTERFACE_COMPATIBILITY_V23.md`, 원시CSV, LP dual과 실제
   `bulk_dispersion_audited_partition`, `shortwave_direct_energy_refined`를 읽는다.
2. 같은LJ/Bessel per-atom law에서 단파장 vector 성분을 실제 source타깃으로
   명시하고, 일부분 fit / 나머지 excluded로 사전구분한다. 기존bulk타깃과
   겹칠 수 있는 이번 post-fit검사를 blindheldout이라 하지 않는다.
3. Exact7 뒤 nullspace에서 partial-coordination/finite-q derivatives의
   rank와 충돌을 먼저 분석한다. 단파장오차를 hiddenstiffness/큰saturation으로
   옮긴 채 core force만 줄이는 fit은 채택하지 않는다. 추가term은 수식/단위/
   per-atom합/독립미분/gauge/identifiability를 먼저 유도하고 최소한만 검토한다.
4. 검증후 actualcore/domain/loadreturn을 다시 수행한다. Straightrepeat당
   에너지를 임의line길이로 곱해 finiteactivation/초/항복을 만들지 않는다.
5. 아직 production/UI gate 미통과. 기존 UI 재설계와 physicalHz 활성화는 하지 않는다.

## v23 후속 백업 — 2026-09-11, 최종 회귀 실행 중

아래16:43 기록보다 이 절이 최신이다. 시작 fresh local/origin 모두
`d1ac6695518d3738b8303687d8c6a7b58877203e`, clean probability-pde-solver-v1.
실제 worktree는 `aft-pde-bessel-38969ad`; OneDrive 원래 폴더는 migration stub.
현재 변경은 이 turn의 v23 연구/테스트/인계뿐이며 아직 commit/push 전이다.

### 새로 확인한 핵심 해석

**Hxx 하나의 부호로 slip saddle 부재를 주장하면 안 된다.** Source saddle의
좌표와 후보의 실제 정지점은 다르며, full3x3 Hessian을 써야 한다.
실제 재계산에서 v22 wider의 자기 saddle은 Hxx=+.068616인데도
lambda_min=-1.387798eV/L0²이고 양쪽의 서로 다른 최소점까지 하강했다.
Source/wider/new-radial 후보 총9개 stationary state와3개 saddle connectivity
검사,15개 signed local-traction static 상태가 모두 검증됐다.
이것은 production PDE, 시간 hold, 잔류소성 또는 실제 항복 검증이 아니다.

그러나 에너지 값은 여전히 맞지 않는다(J/m²):

| Model |fault|saddle|W(40h)|
|---|---:|---:|---:|
| Al99 source |.150479|.172002|1.741285|
| v22 wider |.079457|.119750|1.840583|
| v23 radial joint |.075309|.114757|2.524404|

즉 saddle가 없어서가 아니라 **정량적 계면/코어 적합성이 부족해서** 채택 불가다.
Source의 모든 국소 curvature가 실험적 정확값이라는 뜻도 아니다.
Static source곡률 FD는1e-5에서2.475e-7/3.549e-9eV/L0²로 확인했다.

### 완료된 계산

- `fixed_shape_audit`:48 LP; exact7/signs 아래 두 offendingjet 동시맞춤 가능.
  그러나 wider115jet 고정shape minimax18.46445; 모든sign해제18.33840.
  Primal/dual/KKT 저장. 고정shape 필요조건일 뿐 전체함수족 불가능 증명 아님.
- `box_profiles`:4 동일v22 joint QP. 두원래scale box를 맞추면 otherloss
  2276.47→41943.92, source-force RMS.152417→.175221. Rejected후보별snapshot.
- `radial_minimax_corrected`:160 actualprofiles/908.01s, maxfev stop.
  best선택29jet eta9.98358, full115jet에그벡터를대면40.23339. 같은newshape에서
  full115jet다시minimax하면13.92314. 목표unitbox미달이며 globaloptimum아님.
- `radial_full_validation`:full285 관측량/분광tail/R8,R16 sourcecore force검사.
  full/subset replay0, 실제energyjet재현<8.9e-16. New joint RMS.0995689eV/L0
  (원래 .288834보다65.53% 작음), excludedR16 .0969927. 하지만 otherloss4520.30,
  source-coordinatecurvature부호불일치와 위fault/cleavage오차로채택하지않는다.
- `rational_quartic_ablation`:이전v14/v15 rational shape를 별도로명시해5beta검사.
  `K3 I3²/(1+alpha3 I3)`; default에는안넣음.115jetbound모두18.46445,K3=0.
  D3/K3상관 .9044→.999929, 식별성악화. 새로운core전달은미지원이므로
  CurrentMaterialSiteLaw가nonzero saturation을명시적으로거부하도록guard추가.
- `static_topology_corrected`:실제root/연결/개구33,65bracket/15static상태.
  첫 fixed-registry normal peak12.9696/8.92447/6.53654GPa는 idealcoherent
  traction이며 실제Al항복이아님. Source뒤쪽작은extrema는grid차이남음.

### 원자 코어 계산 상태 — 종료 확인 필요

- `box_core_R8r5`:600s root budget stop, 독립inspectionforce .242865eV/L0.
- `box_core_R8r5_descent`:900s descent budget stop, force1.863155,energy-48.332675.
  수렴minimum이나unbounded에너지의증명아님. Checkpoint검사hash저장.
- `radial_joint_core_R8r5`:7Newtoniteration후force1.7161e-8지만
  lambda_min=-.380258829, FD음의곡률확인. 완료한불안정정지점이지minimum아님.
- `radial_joint_escape_plus_R8r5`:L-BFGS trial에서reciprocal series소진,
  force보정없이중단. 물리spinodal로해석금지.
- `radial_joint_escape_plus_newton_R8r5`:같은energy/음의mode로Newton-CG 재시도중.
  종료summary/failure와실제tool출력을확인할것. max900s/120iteration.
- 최초dict버그/seed-domain/source-geometry helper실패도 별도failure.json보존.
  static_topology 최초helper는source Å→L0 geometry접근을고쳐새폴더재실행함.

### 테스트 / 남은 마무리

- targeted48 PASS99.99s, topology추가2 PASS1.40s 실제실행.
- app34 PASS167.68s, skip0. GUI테스트를headless skip으로숨기지않음.
- 첫full은collection후topology테스트가추가되어중단(exit1, PASS아님).
  `test_run_history.json` 보존. 최종full은
  `.cache/core_interface_compatibility_v23/solver_release_final.xml`로 실행중.
- 최종full/code확인→smoke실행→ledger생성→문서확정→diffcheck→freshfetch→
  정상commit/push. `report_v23_validation_ledger.py`는XML실제집계+smoke실제실행.
  report output은새폴더만허용. 연구artifact 변경후manifest재생성필수.
- 보호hash: original8736cb9d…, wider2fd94b89…, source60c8a085… 불변.
  canonicalenergy/PDE/temperature/mobilities/A_c/physicalHz/UIdefault 변경없음.

다음과학질문은source의partial-coordination response를exact7제약뒤남는
scalar/angular 자유도로어떻게표현할지다. 이전cubic/mixed/rational실패를무시하고
항을다붙이지말것. Actualcoreminimum/finite-loop/sourcepopulations/kinetics가
미완료이며, physicala/s M_a,M_s,t0/초Hz 및 실제yield/fatigue는여전히미검증.

## v23 진행 백업 — 2026-09-11 16:43KST (아직 최종 검증/commit 전)

사용자 "계속해"에 따라 재개. Fresh fetch/local/origin 모두
`d1ac6695518d3738b8303687d8c6a7b58877203e`, clean branch에서 시작했다.
실제작업은 동일 `aft-pde-bessel-38969ad` worktree. No main/default/PDE 변경.

`core_interface_compatibility.py`: 기존 exact coefficient matrix의 minimax LP,
독립 primal/dual certificate, source interval을 상수1좌표로 기존 spectral QP에
전달하는 helper. 새 material항 아님. 새 targeted10PASS1.58s까지 실제 실행.
아래 v22 전체638/app34 결과는 이전 turn; 이번 전체회귀는 아직 미실행이다.

완료: `fixed_shape_audit`48 LP/93.12s. 두 잘못된 곡률만은 exact7+부호조건에서
동시일치 가능. 하지만 wider shape의115 inspected curvatures에 대한 최소
worst normalized error=18.46445, sign-free18.33840, bulk5-only13.94048.
이는 해당shape 필요조건이며 전체radial family 불가능성 증명 아님.
Dual support는 saddle_Hxx,v19_new_4_Haa,v19_power_new_2_Haa와K3 bound.
`box_profiles`4actual profiles 완료: 두 source box(원래1scale)를 강제하면
Hxx=-.205067/Haa=.803279로 맞지만 다른interface loss2276.47→41943.92,
coreforce RMS.152417→.175221. PositiveLJ/sampled spectrum/KKT는 통과해도
material은 실패. Saddle-only는 zero-LJ closure라 snapshot채택 금지.

실행 중(완료로그확인 필수):

- `radial_minimax_corrected`: 기존5shape bounded Powell,27fit curvatures+
  위2opening witness, exact7/signs, max160fev/1200s. Full/subset matrix replay
  후 실제실행; 현재 eta18.464→약15.76, 미완료. Spectral/전체자료/force는 후속검증.
- `box_core_R8r5`: both-box 후보의 실제atomic stationary root/Morse,
  v22 wider R8 좌표만 seed, own material/exterior 유지, max600s.
- `radial_minimax/failure.json`: 최초runner duplicate completed dict 오류,
  물리실패 아님, 코드수정후 새dir retry. 덮어쓰기/삭제 안 함.
- `box_core_R4r5`, `box_core_seeded_R4r5`: 각각 과거metadata의 b누락,
  R8→R4 seed-domain guard의 preflight실패. 실제원자이완 전 중단, guard안바꿈.

상세수학/통계해석은 `CORE_INTERFACE_COMPATIBILITY_V23.md`.
남은일: radial/core 끝난실제결과 확인→필요한 전체자료/spectral/core검증→
수치후보와 물리채택 구별→targeted/fullsolver/app/smoke/diffcheck→상세인계/
freshfetch/정상commit/push. SourceAl99는target-only, LJ/Bessel은불변.
Physical a/s mobility/초/Hz 및 실제yield/피로는 이번static연구로 인증하지 않는다.

## v22 최종 연구 인계 — 2026-09-11 08:35KST

이번 약4시간 연구의 계산과 최종 검증은 끝났다. 아래08:10/07:31/06:10 기록은
실행 중 남긴 역사적 백업이며, **현재 상태는 이 단락과 validation_ledger가 우선**한다.
시작 fresh fetch/local/origin은 `a8e4130934ab8bb67f2a657a489d269d10b20b73`, clean,
branch `probability-pde-solver-v1`. 실제 작업은 이전에 이전된
`aft-pde-bessel-38969ad` worktree다. OneDrive의 migration stub에서 재개하지 않는다.
이 문서를 포함하는 커밋과 원격 상태는 재개 시 git으로 다시 확인한다.

### 완료된 이론·계산과 물리적 결론

- 최신 후보의 **22채널 무한 LJ/Bessel row + 원자별 전체 비선형 환경항**을
  vector screw core에 연결했다. 기존17채널 연구 인터페이스로는 전체 후보항을
  전달할 수 없던 연결 문제를 해결한 것이며, production에서 항을 누락해 실행했다는
  뜻이 아니다. 독립 직접 원자합: small-core energy6.66e-16eV,
  gradient2.73e-14eV/L0 오차. Gradient/Hessian 및 matrix-free/dense 교차검증 완료.
- 원래 candidate/source/첫 trial/넓힌 trial의 실제 static stress states **56개**,
  양방향 resolved shear0→5→20→50→0MPa를 계산했다. 각 물질의 own C/exterior 사용.
  Nominal50MPa의 직접 atomistic stress49.9891/49.9884/50.0034MPa:
  MPa/GPa1000배 환산오류로 강도 차이가 생긴 것은 아니다.
- 원래 source-core force RMS .2888336013→넓힌 후보 .1524166786eV/L0,
  **47.23% 감소**. 기존 exact7 bulk/tangent anchors 유지(<5.95e-14 scaled residual).
  다른 계면 normalized loss2215.203→2276.471, **2.77% 증가**.
  새 energy항/항복목표/임의 strength multiplier 없이 기존 family를 profile했다.
- 넓힌 후보 실제 R8/r5 index-one reconstruction saddle와 양쪽 descent 완료:
  barrier .009778537325eV/straight-row repeat, 원래 .014834588765보다 **34.08% 감소**.
  그러나 R10 sampled25–75% span .87783192L0 vs source2.79084554L0로 아직 너무 좁다.
  이 barrier는 core reconstruction이며 full translation Peierls/finite-loop activation이 아니다.
- 최종 후보는 **채택 불가**: registry-saddle Hxx +.182166 vs target-.186425,
  제외 opening Haa -1.234287 vs +.903279eV/L0²로 부호가 틀린다.
  Outer optimizer도 max_nfev15/bound-active 종료이며 수렴한 유일 Al parameter가 아니다.
  Local numerical sensitivity rank5/condition7.494는 구조적 식별성/신뢰구간이 아니다.
- Frozen force ring7→9 RMS1.216e-6/max3.240e-6eV/L0, tolerance refinement 표시차이0:
  현재 material force error .1524는 이 측정 수치 변화보다 훨씬 크다.
  LJ transverse tail의 analytic Hurwitz-zeta 상계도 유도/검증했다. 이 상계는
  **pair force만** 보증하며 nonlinear environment나 infinite free-domain 인증은 아니다.
- 넓은 source 초기값을 사용해도 넓힌 후보는 독립 회복계산에서 같은 좁은 core로
  돌아왔다(field8.13e-10L0). R8 ring5→7 field9.27e-7L0,
  R8→10 same-ring7 field.00185908L0. 단일 정밀 실행을 무한영역 수렴이라 하지 않는다.
- 넓힌 후보 ±50→0의 core return 오차7.98e-9/1.03e-9L0: 이 경계에서 지속하는
  registry shift는 확인되지 않았다. Source는 R16 j+1, R20 j+2 shifted zero-load
  endpoint를 독립 재현했지만 이동거리가 domain 의존한다. 실제 거시항복/잔류소성 인증 아님.
- Source R16의 두 shifted branch 사이 index-one saddle와 양쪽 descent도 완료:
  forward .00208837349/reverse .00033034417eV/repeat. Source-only 비교이며
  우리 potential 대체 아님. Static load-return은 시간적 hold/피로/Hz 결과가 아니다.
- 같은 법칙의 screening exponent z=0,±.5,±1 다섯 profile도 실제 실행했다.
  이 제한된 scan에서는 joint loss 개선 없음. 전체 analytic family 불가능성 증명 아님.
- Finite line-period B=Nb의 basis별 Bessel regrouping과 self-term 유도는 문서에
  추가했다. **유도만 완료**, kink/loop solver·finite activation·kinetics는 미구현이다.

### 최종 산출물과 검증 — 실제 실행 완료

상세 유도: `solver_v1/CURRENT_MATERIAL_CORE_BRIDGE_V22.md`.
최종 결과 안내: `results/current_material_core_v22/README.md`.
최종 비교는 `final_comparison`, 실제 에너지 재생은 `wider_replay_final`.
완전한 fit residual/제외 상태/parameter snapshot은 `wider_probe_validation`.
`validation_ledger`: 원자계산 항목107개, 실패 시도4개 보존, 미완료 atomic case0개.
실패한 Newton/L-BFGS budget-stop은 별도 성공 recovery로 덮어쓰지 않고 함께 남겼다.
Optimizer budget termination은 원자계산 완료 여부와 별개의 미해결 보정 상태다.

최종 실제 검증(모두 failures/errors/skips0):

- targeted77 passed,48.56s (`targeted_final.xml`).
- solver_v1 **638 passed,1195.90s** (`solver_final.xml`).
- app **34 passed,60.64s** (`app_release_final.xml`); 더 이상 GUI2개 skip 없음.
- ledger가 desktop smoke 재실행: exit0,0.90666s, 기존 a0/kappa 확인.
- 실행 XML 요약/hash와 smoke stdout은 `validation_ledger/executed_tests.json`.
- Artifact byte inventory는 `validation_ledger/artifact_hashes.csv`.
  이 manifest 생성 후 결과 파일 내용을 바꾸면 반드시 ledger를 다시 생성한다.
  Local .gitattributes가 scientific result byte/newline을 보존한다.

Production PDE/UI/default energy, 과거 static calibration, mobility와 A_c는 불변.
원래 candidate SHA8736cb9d28f1430e3991b1046ef42544cefaa3c9a4f5d743e0d02931329f5c2e
보존. 새 후보는 별도 research snapshot이며 SHA
2fd94b89168b173378f4c67d9d6faf52c0363292fc66d8d7e705f214a55de23a.
Source Al99는 기존 SHA60c8a085be79d273324ab421f5b1447578fef55c1acfc6492c0999f15ee8a284.
Physical a/s M_a,M_s,t0/seconds/Hz는 계속 unavailable; 생산/UI gate는 닫혀 있다.

### 다음 연구 순서

1. 같은 source 상태의 atomic force와 local opening/registry jets를 함께 제약하되,
   현재 sign conflict를 풀지 못한 원인을 analytic coefficient/shape sensitivity로
   더 조사한다. 제한된 local search 실패를 family 전체 불가능으로 확대하지 않는다.
2. 최소 추가 analytic term을 검토한다면 기존 family의 구체적 compatibility/rank
   실패와 독립 제외상태 개선을 먼저 요구한다. LJ/Bessel·per-atom law 유지.
3. 새 material이 나오면 반드시 실제 core 재이완→load-return→domain/ring/saddle
   검증을 반복한다. Frozen-source force error 개선만으로 core/강도 개선이라 하지 않는다.
4. Finite-period atom basis를 구현하기 전 이 문서의 self exclusion/row regrouping을
   시험한다. Straight-line energy에 임의 길이를 곱해 activation energy를 만들지 않는다.
5. 실제 source population/interactions, finite line/loop energy, 집단좌표 mobility가
   있어야 실제 항복·피로·시간으로 이어진다. 현재 static core 결과로 이를 인증하지 않는다.

## v22 추가 진행 백업 — 2026-09-11 08:10KST (아직 최종 commit 전)

아래07:31 이후 `core_informed_wider_search`가 실제70profiles/945.07s 완료.
선택profile=accepted endpoint이지만 max_nfev15로 outer stop, optimality5.98e-4,
even saturation log상한active. Core force RMS .15241668 (기존.28883360보다47.23%
감소), other-interface loss2276.471 (기준2215.203보다2.77% 증가). 모든exact7
scaled residual<5.95e-14. 물리채택 아님: registry saddle curvature 부호오류,
제외 opening curvature 부호오류, 중간opening force 약32%부족을 docs/CSV에 공개.
`wider_probe_validation`의 실제제외R10/R16 RMS .15133584/.1497450, perturbedR16
.1554813/.1548614. Local SVDcondition7.494, numerical rank5일 뿐 구조식별성 아님.

새wide candidate를 실제실행:

- `wider_core_near_R8r5`: stable force4.889e-10, minH.19770670, E.875860833256.
- Sourcewide 초기값 Newton은 reciprocal40 exhaustion. `wider_core_R8r5` 실패보존.
- `wider_core_lbfgs_R8r5`는480s budget stop, 미완료이며 checkpoint보존.
- 그checkpoint `wider_core_recovered_R8r5`는 독립stationary root로
  force3.752e-11, minH.19770670, E.875860833256 통과. 다른seed와 비교예정.
- `wider_core_R10r7`: 실제stable force3.680e-7(default5e-7이하), minH.13952580.
- `wider_replay`가 위3상태 에너지 독립재계산0차이 확인.
- `wider_frozen_refinement`: ring7→9 force최대3.24e-6/RMS1.22e-6eV/L0,
  재료오차.1524보다 작음. Tol2e-12→2e-13 표시정밀도에서0차이.
- `existing_screening_profile`: 기존 x^z의z0/±.5/±1 5개 실제QP. Zero기준replay,
  어느점도joint objective개선 못함. 새potential항/새보정기본값 없음.
- SourceR20의j+2 독립zero seed가 +50→0 결과와 E9.33e-15/inner2field1.63e-9
  일치. R16과 이동거리다름: 무한계Peierls/거시잔류변형률 인증 아님.

08:10 실행중(완료로그확인 필요): `wider_stress_R8r5` 양방향8상태,
`reconstruction_wider_R8r5` saddle양쪽실제descents, `wider_core_R8r7` 독립ring검사.
Targeted최종77PASS48.56s(`targeted_final.xml`), 전체최종재실행중(`solver_final.xml`).
이미 앞전체638PASS823.59s(`solver_release.xml`), app34PASS128.74s/skips0 실제완료.
Core/stress끝나면 `report_core_research_v22 --trial .../core_informed_validation/
research_candidate_snapshot.json --wider .../wider_probe_validation/research_candidate_snapshot.json`
로 새report 생성. source/초기trial/넓힌trial을 구분. 그 후 docs와README 최종수치,
ledger실제smoke/해시, stage검토/diffcheck/freshfetch/정상commit/push가 남음.
새스크립트 `run_core_screening_profile.py`는 기존shape별5점검사이며제로profile
replay를 강제한다. AGENTS/GATES/상세유도/README 갱신, default/production불변.
Finite-period Nbasis/B=Nb의 Bessel regrouping은 문서에 유도만 추가; finite-kink
구현/activation energy/새PDE로 잘못 말하지 말 것. 물리seconds/Hz 계속없음.

## v22 최신 진행 백업 — 2026-09-11 07:31KST (최종 commit 전)

사용자 최신 지시: 약4시간 동안 연구 계속. 시작04:49KST, 대략08:49까지.
실제 worktree는 migration 이후 `aft-pde-bessel-38969ad`, branch는
probability-pde-solver-v1. 이 turn fresh fetch/local/origin 모두
`a8e4130934ab8bb67f2a657a489d269d10b20b73`, 시작 clean. 이후 변경은 이번 연구;
아직 commit/push 안 함. Main/production/default/기존 material parameter는 불변.
이 단락 아래의06:10 백업은 중간 기록이며 여기와 충돌하면 여기가 최신이다.

핵심: 같은 LJ/Bessel의 **full current22채널 vector atomistic screw core**를
실제로 연결/검증했다. 이론/정량은 `solver_v1/CURRENT_MATERIAL_CORE_BRIDGE_V22.md`.
Al99는 target-only source이며 우리 potential/production 대체 아님.
Source R20까지 늘려 양방향0→5→20→50→0MPa8상태를 추가 실행했다.
Original R8, trial R8, source R8/R10/R16/R20 =48 static stress states.
Candidate/source 각자 자신의 C를 사용. Model-time PDE/kinetic hold는 실행 아님.

정량적으로 끝난 내용:

- Original core와 source core의 sampled25–75% span: R10에서.87408 vs2.79085L0,
  source R16/20는2.83015/2.83868. 단순 sampled span이며 continuum core radius 아님.
- Same-source-state original gradient RMS .2888336eV/L0. 고정shape exact7/null3
  부호까지 풀어도 하한.1132313이며 K3<0 필요. 선언부호/sampled spectrum QP는
  .1152431까지 되지만 다른계면 loss2215.20→17100.51(7.72배). 채택 금지.
- Joint local shape 첫38profiles/492s: core RMS .166348, 계면loss2473.438(+11.66%).
  best는 optimizer endpoint 아님. `core_informed_validation/research_candidate_snapshot.json`
  에 별도 보존했고 default를 바꾸지 않았다. Original parameter SHA8736…2e 불변.
- 그 trial을 실제 재이완: 첫R4 root는 saddle, mode descent 후R8/R10 안정.
  R10 Newton-CG와 full-matrix stationary root E차1.93e-14eV, field차3.19e-8L0.
  trial span .87750L0로 원자 force error 개선과 달리 core 폭은 거의 개선 안 됨.
  trial R8±load return field오차1.89e-9/1.35e-9L0. 잔류시편소성 아님.
- 넓은 source core를 좌표 초기값으로만 복사해도 두 후보 모두 좁은 원래 core로
  돌아왔다. Own material/exterior/C 유지. Original 첫 Newton은 reciprocal exhausted,
  L-BFGS는1000s budget stop(둘 다 failure.json 보존). 그 checkpoint를 명시적으로
  재개한 `source_seed_current_recovered_R8r5`는 force4.36e-11/minH.2162692 통과;
  original stable state와 field차1.48e-10L0. Trial independently차2.80e-9L0.
- Source R16의 +50→0 endpoint는 independent j+1 seed와 field2.98e-9L0/E1.78e-15
  일치. 실제 index1 연결saddle: forward .00208837349, reverse .00033034417
  eV/직선row repeat, 양쪽 descent endpoint 확인. Source R20에서는 제하 후 더 큰
  shift가 남아(domain 의존) 무한계 Peierls/실제 yield 인증 아님.
- 실제 farfield stress: nominal50MPa→original49.989066/trial49.988367/source50.003387.
  MPa/GPa 환산오류 아님. Normal prestress의 algebraic ring tail은 별도 미수렴:
  original ring7/11/15 -1.49/-.417/-.170MPa; 실제 load reaction과 구분한다.
- `core_energy_mechanism`: 좁은→source넓은 동일 내부경로(각자 own exterior)에서
  E변화 original+.210119, trial+.097749, source-.150291eV/row repeat.
  현재 gauge의D3/Eg/K3가 주된 양의 비용. 이 경로는 MEP/activation barrier 아님.
- `core_pair_tail_bound.py`: 무한 transverse LJ force remainder를 Hessian majorant,
  8k square shells, Hurwitz-zeta로 엄밀 상계. Force 수정/finite cutoff 도입 아님.
  Only pair tail: nonlinear environment/free-domain/material certificate 아님.
  `analytic_pair_tail_with_source`에 own/source frozen 상태 및ring3/5/7/9/13 검증.

계속된 동일목적함수 최적화:

`core_informed_shape_continuation` 실제27profiles 완료, xtol success지만
optimality2.26e-4>gtol이고 saturation 상한 active. Loss1.44825896,
core RMS .166719, static2470.133. 채택 아님. 제한된상자의 수치종료로 명시.
현재 `core_informed_wider_search`가 같은 source/목적함수/7anchors/부호/finite-q를
유지하고, 원래 baseline 중심 log범위만±.2→±.6으로 명시 확대해 실제 계산 중.
새 energy term이나yield 목표는 없다. max_nfev15, wall1000s; 끝나면 반드시
completion/실제외부검증/optimizer endpoint 여부 확인. 실패/미완료면 그대로 보고.
`continued_probe_validation`은 앞 continuation의 제외 force검증 계산 중.

실행·검증 상태(07:31기준):

- 최초 전체627 passed780.60s; 후속 전체632 passed1276.46s 실제 완료.
- 최신 targeted77 passed22.56s: `.cache/current_material_core_v22/targeted_release.xml`.
- 최신 전체(추가tail/ledger6개 포함) 실행 중: `solver_release.xml` 결과 확인 필요.
- app34 passed128.74s, skipped0: `app_full_final.xml`.
- 앞 smoke실제exit0; 최종ledger가 smoke를 다시 실제실행한다.
- diff check 현재통과; untracked를 stage한 뒤 cached diff check도 필수.
- `report_core_validation_ledger.py`는 실제XML+실제smoke+case상태/hash를 모은다.
  모든 계산 종료 후 새out예:`validation_ledger`에 실행할 것. hash manifest 이후
  같은 study결과파일을 수정하면 새ledger 필요. Tests PASS와material 채택을 구별.
- `study_summary` 현재 completed48stress states/독립seed확인/기존trial comparison.
  연구단계별 old `combined_analysis`, `final_analysis*`는 과거중간 replay이다.

남은 마무리: 넓힌 search 결과 확인→제외검증(좋아져도 재이완 없이는core개선 아님),
현재 실험/실패/수치결과를 docs/README에 요약, full regression 최종확인,
실제ledger 실행, 변경파일검토·stage·diff check, fresh fetch·정상commit/push·원격확인.
Main기준 local80cacb4180dbfbcd36a2964270703bc6cf1653ec,
remote c43d8e096f2a329c7cdfad31e546c70fb36fe592 유지 확인.
실제 Al 항복/피로/production a/s mobility·초Hz는 여전히 미검증. UI승격 금지.

## v22 진행 중 — 2026-09-11 약04:49KST 시작 (아직 최종 검증/commit 아님)

### 06:10KST 부근 추가 백업 — 최신 진행은 이 단락 우선

22채널 current-material/nonlinear vector row core 연결과 독립 실원자 검증 완료.
Source-only Al99도 같은 원자열/벡터/경계에서 실제 계산했다(우리 potential 대체 아님).
Current R4/6/8/10의 ring5/7 안정 branch를 실행; R8→10 inner2L0 최대변화
.00208094L0, ring5→7 .00000031636L0. 무한domain/Peierls 인증은 아직 아니다.
Morse index1 reconstruction saddle: current R8/r5 .01483458876eV/row repeat,
R10/r5 .014553901, R8/r7 .014834180; source R8/10/12/16 각각
.0003362987/.0006445512/.0008267211/.0010162721. Source domain효과 큼.
R8 current/source 양방향 descent는 실제 두 안정 endpoint를 확인했다.
이는 core symmetry reconstruction이며 full-lattice glide 또는 유한 loop activation과 다르다.
Current R8 및 source R8/R10의 0→±5→±20→±50→0MPa 모두8상태 안정.
정적 load-return만 실행; dynamic hold/잔류변형률/실험적 항복 인증 아님.
Source spline Hessian coarse energy FD 오차는 실제13-step refinement 수행:
gradient-directional curvature는 analytic eigenvalue에 약2.84e-11까지 일치.

`force_compatibility/`: source stable core에 동일 위치로 current energy를 적용.
원자 gradient RMS .2888336eV/L0(inner2). 기존 exact7조건을 유지한 고정-shape
계수nullspace는3차원; 부호제약마저 풀어도 최소RMS .1132313, K3<0요구.
`constrained_force_profile/`: 양의LJ/선언부호/기존 finite-q tail constraints
유지하면 RMS .1152431, K3=0, 다른계면 loss17100.51로 악화(기준2215.20).
따라서 어느 vector도 채택하지 않았다. 전체analytic family 불가능성 증명 아님.

`core_informed_shape_probe/` 실제 실행 중: 에너지항 추가 없이 기존5개log-shape만
±.2의 작은영역에서 source core+계면 공동profile. objective는 두 data block의
기준 squared loss로 나눈 합(weight1), 물리uncertainty/유일Al fit 아님.
모든trial/checkpoint보존, max_nfev8/max1800s. endpoint와 최적trial 구별필수.
신규+기존 core targeted57PASS10.80s. 전체 solver/app 회귀는 새로 실행중;
완료로그 나오기 전 PASS 금지. 현재 core/study files 외 생산/PDE/UI/물성/Hz 불변.

최신 사용자: 베셀 기본틀과 전위선의 연결을 설명한 뒤 아침까지 약4시간 연구 요청.
Fresh fetch 성공, local/origin a8e4130934ab8bb67f2a657a489d269d10b20b73, clean.
기존 문서/원자열 코어 코드를 확인했다. v21은 C와 b만 쓰는 outer line이고,
초기 isolated-core의17채널은 최신 v20의 Eg/Quartic/even saturation 전체를
지원하지 않는다. 최신 후보를 base.surface로만 코어에 전달하면 항 누락이다.
과거 production에서 이 누락을 실행했다는 뜻은 아니다; 아직 연결하지 않았다.

현재 작업: `current_material_rows.py`의22채널 무한 Bessel row와 원자별 전체
에너지/gradient/Hessian. `isolated_screw_core.py`에 kernel/site-response hook을
넣되 default energy는 유지한다. 새 `test_current_material_rows.py` 및 기존
core targeted 실행 중. 유도/계획 `CURRENT_MATERIAL_CORE_BRIDGE_V22.md`.
모든 결과는 연구 진단이며 v20 static-material gate 불통과는 그대로다.
초기 targeted29PASS15.53s, 확장 row/core25PASS16.56s,
protocol11PASS1.65s, source+legacy15PASS11.18s. 전체회귀는 아직이다.
`run_current_core_validation`은 실제 독립 실공간 원자/force 합산을 실행했다:
nonuniform small core fine energy error6.66e-16eV, force2.73e-14eV/L0.
이것은 에너지표현 검증이지 rare probability 또는 재료 정확도 인증이 아니다.

v22 현재계수 SHA8736cb9d28f1430e3991b1046ef42544cefaa3c9a4f5d743e0d02931329f5c2e
불변. 최초R2/r3 centered force1.19e-7이나 minH=-2.17874: saddle.
양쪽 mode descent에서 E=.68878052884, minH1.58268 stable fixed-boundary.
R3/r4, R4/r5, R6/r5, R8/r5 실제 안정해; R8 minH.216269.
R4 Newton/기존 L-BFGS E차4.6e-14eV. R4 j+1 결정학 translation 초기화는
다른 안정위치로 수렴했다(외부경계 고정; 무한계 Peierls/유한 activation 아님).
R4 0→50→0MPa 계산완료; 정적 제하이고 kinetic hold 아니다. 자세한 비교 진행중.
R4/R6/R8 ring5→7 독립refinement 실행 중. `run_current_material_core`의
mode/source stability gate를 준수한다. relative path `.resolve()` 오류수정;
최초 실패계산은 새결과를 만들기 전에 중단되었다.

추가 검증범위: 같은 FCC/벡터/경계 프로토콜의 **target-only Mishin Al99 core**.
`source_core_reference.py`, `run_source_core_reference.py`. published source cutoff는
source에만 사용, 우리 LJ/Bessel교체/생산selector등록 아님. 각자 자기C를 사용.
NIST 원문 페이지 재확인; 파일hash 기존60c8…e8a284 그대로. Source bulk/row
탄성항도 독립일치. Source R4/R6 안정이나 R8 centered minH=-.01736523였고,
그상태를 stable로 넘기려던 R10 실행은 guard가 거부했다(결과directory미생성).
원본R8의 path명 stable_R8r5는 요청label일 뿐 실제 summary는 unstable이다.
현재 양쪽 negative-mode descent 후 R10 확장검사 중. Source stable상태와
후보 core를 같은 물리상태/단위에서 비교한 뒤 다음 판단을 한다.
실제항복/피로/물성/초Hz 보정 성공 없음. 전체회귀/최종검증 전commit금지.

## 실제 항복 연결 v21 — 2026-09-10~11 (계산·전체 회귀 완료; 실제 항복 미보정)

최신 요청은 이상강도를 낮추지 말고 같은 LJ/Bessel에서 결함을 통한 실제
항복으로 나아가라는 지시 후 "그치 계속해라"이다. Fresh fetch 성공;
시작 local/origin 모두 `a0b611b1b3b29c6bef86dcb2e4b5b23c79228d04`, clean,
branch probability-pde-solver-v1, 실제 worktree aft-pde-bessel-38969ad.
OneDrive는 migration stub이다. Main/production/기존보정/이동도/Ac 불변.
새 agent 없음. 현재 사용자 우선순위는 **source motion→시편 변형률→항복**;
아래 v20의 scalar embedding refit만 무작정 다시 시작하지 않는다.

읽을 문서 `solver_v1/SPECIMEN_YIELD_AND_SLIP_BUDGET_V21.md`.
새 모듈 `specimen_slip_kinematics.py`는 signed swept area의 tensor 관측량이다:

    beta_slip=Σ A_k b_k⊗n_k / V
    eps_slip=sym(beta_slip)
    V sigma:deps=Σ(b·sigma·n)dA

실제 물리적 specimen volume/표면 기하학 정규화이며 activation volume이나
Ac가 아니다. 기존 positive scalar area helper를 바꾸지 않고 negative slip,
여러 slip system, frame rotation, initial-slip subtraction, reversible/gross
분리를 추가한다. 운동법칙/empirical plasticity/production PDE가 아니다.
확률밀도에서 같은 observable의 transport와 opening-selective 항도 유도했다.

실제 runner `run_specimen_yield_bridge.py`는 새 directory에만 실행한다.
v20 positive old-family의 이전 bulk C11/C12/C44=114/62/32GPa를 hash결합
재사용; 새 fit 아님. 기존 same-energy Schur line을 edge 방향으로 계산한다.
Hypothetical L=.2/1/5um, R=L, r_core=2b 그대로; 실측/항복 튜닝값 아님.
[100]/[110]/[111] crystal projection을 명시했다. 72 stress cases
(0,2,5,10,15,20,30,50MPa), 297 static cycle states, 36 independent grid checks,
9 dynamic line hold histories를 실제 실행했다. Source 각도32/64,
적분1e-9/1e-11 비교 포함. 수치 본체17.18s 완료.

새 정량적 진단: L=1um, [100]의 outer-only first-source 축응력24.03976MPa지만,
선언한 dilute eta=NL³/V=.01에서 fold 최대 축변형률6.89449e-7이다.
0.002를 순간 휨으로만 만들려면 eta=29.0087 (약2900.87배 strain demand)이어서
선언한 희박 독립-source 가정과 양립하지 않는다. 전체9경우 요구 eta=5.8017–217.5649.
이는 필요 density 진단이지 fit할 density가 아니다. 모든 전위 배치에 관한
불가능성 증명도 아니다. 임계 분해전단9.81419MPa와 실제 시편 YS를 등치하지 않는다.
Static 제하에서는 signed 변위0/gross>0. 별도 기존 linear-line 계산을 실제 실행:
6 cycles와 >24tau hold, finest peak3.87441e-10, final-1.25857e-20,
ratio3.2484e-11; decay이지 잔류소성이 아니다. Hold 끝시각은64 steps에서
24.05282tau, 128/256에서24.00373tau이므로 엄밀히 같은 끝시각 수렴표는 아니다.
V18 조건부 line drag를 재사용했으며 production a/s 초·Hz 보정은 아니다.

새 primary source Pigato et al.2026 doi10.3390/ma19061195 원본 XML을 다운로드했다.
`yield_reference_data.py`가 Table2의 YS12개를 직접 파싱한다. 실온6N16.4±1.3,
5N5 74.6±3.7, 5N34.1±1.1MPa는 서로 다른 초기조직이다. 원문은 YS로만 표기하고
정확한 수치 offset 명시가 없어 plastic_strain_criterion=null이다.
실측 source 간격/이동전위 밀도 없음; 이 값들은 loss에 넣지 않았으며,
측정된 재료 항복으로 candidate를 인증하지 않는다.
원본은 .cache/specimen_yield_bridge_v21/pigato2026.xml; XML SHA는 assessment.json.

결과 `results/specimen_yield_bridge_v21/`: strain_budget, stress_sweep,
static_cycle_transport, spatial_refinement, unload_hold, reported_yield_references
CSV와 assessment/plot. 도면을 실제로 보며 log-x label 중첩을 수정했다.
독립 재생을 .cache/specimen_yield_bridge_v21/reproduction에 완료(6.48s),
여섯 CSV가 원본 실행과 바이트 단위로 일치함을 확인했다. 도면만 최종 갱신했다.
Targeted62PASS4.01s, 최종 재실행62PASS1.75s; app34PASS65.89s, smoke exit0/1.28s.
첫 임시 test에서 kinetic JSON의 status 대신 실제 calibrated 필드를 써야 한 오류를
수정했다. 잘못된 기존 test 파일명으로1회 collection 실패 후 정확한 목록으로 실행했다.
전체 solver591PASS988.27s(16분28초), skip/failure 없음. 최종 JUnit은
.cache/specimen_yield_bridge_v21/, 검증 요약은 results/specimen_yield_bridge_v21/
verification.json. 여섯 결과 CSV hash와 독립 재생 일치를 보존했다.
Commit/push 직전 diff/원격 확인을 수행하며, 정확한 최종 SHA와 push 여부는
최종 응답 및 git log/fetch로 확인한다.

다음 필요한 기구는 minor branch 소실 뒤의 실제 전위 방출/유한선 전파/상호작용과
same-energy의 검증된 core 및 실측 source/obstacle population이다. 현 모델을
임의 밀도로 곱해 YS에 맞추지 않는다. Static material/interface gate 불통과도
독립 과제다. 실제 항복/피로 승인 미완료; Ma_phys/Ms_phys/t0 null,
production 초/Hz disabled.

## 보정 재개 v20 — 2026-09-10 (범위 내 계산·회귀 완료; 물성 채택하지 않음)

최신 요청 "계속해". Fresh fetch 성공; 시작 local/origin 모두
`56e1c3ab0f570bf7628a320e03b7cf4eb3f1cb22`, clean,
branch `probability-pde-solver-v1`. 실제 작업은 기존
`aft-pde-bessel-38969ad` worktree이며 OneDrive 폴더는 이동 안내 stub이다.
Main/production PDE/UI/kinetic/static parameter 불변. 새 agent 없음.
최종 commit/push SHA는 git log와 fresh fetch 및 최종 응답으로 확인한다.

### 이론과 실제 보정

읽을 파일: `solver_v1/TANGENT_CONSTRAINED_CALIBRATION_V20.md`,
`interface_tangent_calibration.py`. v19의 약41% 높은 initial Haa를
고치려면 기존 에너지 family 안에서 어떤 tradeoff가 생기는지 유도했다.
새 energy항 없이5개 exactbulk와 source pristine Haa/Hxx를 동시에 제약한다.
상태/단위는 E0=1eV, L0=4.05/sqrt(2)Å, atomic area=sqrt(3)L0²/2이다.
Source는0K Al99 rigid interface target-only, fatigue/yield/실험 kinetic 아님.
old_family라는 새 결과 레이블은 v19 p=0 연구 family이며 TwoRowLJ가 아니다.

Fixed-shape CONTROL16개 완료52.95s. Both-exact loss old2831.110021,
power2859.423588이고 u>0,v=0 closure이다. Phi=u/r^12-v/r^6에서
v=0은 finite positive LJ가 아니므로 물성으로 채택하지 않는다.
Source target와 imposed control은 CSV에서 별도로 유지한다.
Fixed-shape V(t)의 convexity를 유도하고 actual scan으로 확인했다.
이 결과만으로 nonconvex 전체 family 불가능성을 주장하지 않는다.

Actual shape refit82old+162power profiles=244개,3126.49s 완료.
Controls포함 총260 recorded profiles(추가 selected-profile replay 별도).
기존 bounds/동일102 loss/48새 사전선언 heldout분리, 임의 양수 floor없음.

- Old:17nfev/13njev,ftol종료,optimality1.408. Last accepted loss2012.434839;
  최저 difference trial2012.434779도 v=0이다. 선택된 positive trial
  loss2215.203159는 optimizer endpoint가 아니다. u=.00469129023710,
  v=.0107489011795,eps=.006157094889eV,sigma/L0=.8709396023.
- Power:30nfev/22njev **budget종료, optimizer수렴 아님**,optimality4.716.
  Last accepted loss2014.002599, 최저 difference trial2014.002476, v=0.
  선택된 positive trial loss2223.200136, u=.00673566805554,v=.0198388390094.
  이것 역시 optimizer endpoint가 아니다. 추가 shape가 이 실험을 개선하지 못함.
- 두 positive trial은 bulk5+Haa19.6609128633/Hxx4.39359303512의7등식,
  sign/샘플 spectral/KKT를 충족한다. Exact residual<4.45e-14,KKT<1.44e-11.
  Feasible coefficient profile과 채택 가능한 물성/전역 최적해는 다르다.

### 독립 검증과 실제 저응력 응답

`final_validation/`에서 source+v19parent+두v20후보를 동일하게 검사했다.
11tensor scenarios×9signedstatic states×4models=396root 모두 검증.
50MPa normal delta a[Å]: source.0009253659034,parent.0006562111178,
old.0009286879722,power.0009286672203. 오차 -29.09%→+.3590%/+.3568%.
4MPa shear1 slip source.0003309978847 vs old.0003309946228Å.
초기 tangent를 fit한 결과여서 independent yield검증으로 부르지 않는다.
Shear의 secondary normal변위 old8.67e-8Å vs source2.57e-8Å는 여전히 다름.
전체 성분/혼합하중/음의 하중 모두 CSV에 있고 staticreturn≠dynamic hold이다.

48새 excluded jet RMS(energy,force,Haa,Hxx):
parent=(1.0025,5.6153,6.0263,6.1650),
old=(.8786,7.8649,6.0596,5.1335),
power=(.8752,7.8780,6.0714,5.1376).
두 후보 모두 force는1/12만 선언 scale이내: **full material gate실패**.
Fault/saddle[J/m²] source.150479/.172002,old.122424/.162561,
power.123257/.162566. First fixed-registry ideal traction source12.96957GPa,
old10.33897,power10.33844GPa. 실제 항복/동적 eventordering이 아니다.
W(40h)[J/m²] source1.741285,old1.883279,power1.887106. Candidate는40h에서도
약한LJ attraction이 있어 정확한 무한분리라 하지 않는다. Source의 이후 cutoff
진동/negative traction도 그대로 기록했다. 두 bracket에서 source의 첫 peak는
일치하지만 이후 작은 minimum위치는 약.000588h차이가 있어 무조건 수렴이라 안함.

수치: finestFD Hessian error3.7525e-7eV/L0². Tolerance2e-11/2e-13,
per-site direct6/10/16, full nonlinear/BlochFD를 구별해 기록.
5dilation×243q×radius12/16 minimumtail-subtracted margin.00565074.
모든 q/strain 안정성 증명 아님. 최종validation walltime1221.61s.

### 원인 분해와 식별성

`final_tail_audit/`는3후보/sourceSHA/계수hash 모두 결합된 최종 tail감사다.
Al99 declared cutoff6.28721Å, effective support6.286581279Å.
LJ meanforce Hurwitz식/large-a~a^-3를 유도, exact mean오차5.42e-20eV/L0.
이는 확률정밀도 주장이 아니다. a/h2.5의 old total1242.61MPa 중
LJ29.86, 환경항1212.74이며 그중scalar1093.50MPa. Source4.06MPa.
따라서 큰 개구 force오차의 주원인을 unavoidable LJtail로 돌릴 수 없다.
Source cutoff를 candidate에 도입하거나 target를 버리지 않았다.

Exact7등식 tangent SVD: oldrank8/8,condition300.943;
powerrank9/9,condition4280.321. Saved-vector replay오차0,
2e-4→1e-4derivative step변화<1.925e-4. Activeinequalitycone/통계CI/
전체 family식별성증명은 포함하지 않는다. Powershape의 conditioning악화도 기록.

### 재개 경로 / 끝난 것과 남은 것

결과root `results/fcc111_active_interface/tangent_calibration_v20/`:
`controls/`,`shape_refinement/`는 실제 optimization;
`old_report/`,`power_report/`는 저장계수 replay/SVD;
`final_validation/`,`final_tail_audit/`,`response_report/`는 최종 비교다.
`old_independent_check/`,`tail_diagnosis/`,`tail_decomposition/`는 보존한 중간완료자료.
`checkpoint.json`,`progress.json`의 completed=false는 중간스냅샷이며
최종 completion/calibration와 optimizerstop이 우선한다. 원시CSV/JSON의
bytehash보존은 local .gitattributes로만 설정했다.

기존QP empty-inequality SciPy 오류를 최소재현 후 고쳤다. 모든 부등식이
등식manifold에서 상수일 때 해석적 least-squares해를 equality/KKT/원래상수
잔차/sign으로 검증한다. Varyingconstraint삭제/계수clipping은 없다.
최종 targeted47PASS59.99s, fullsolver565PASS923.85s,
app34PASS86.14s, desktop smokeexit0/.977s, 실패/오류/skip모두0.
JUnit은 `.cache/tangent_calibration_v20/*_final.xml`, 커밋요약verification.json.
Working/staged diffcheck PASS. 원시JSON/CSV/PNG102개 index/working bytes동일
검사 완료. Per-study 속성으로CRLF를 인식하되 다른 whitespace검사는 유지한다.
최종freshfetch에서도 origin은 시작56e1c3a와 동일, main/originmain불변이었다.
UI scroll은v18에서 이미 구현·검증되어 이번에 변경하지 않았다.

다음 연구는 **partial-coordination scalar embedding의 derivative 자유도**다.
먼저v14cubic-density와v15positive-mixture의 기존실패/계수예산을 읽고,
7등식하 남은 독립방향/rank/force·curvaturetradeoff를 유도한다. 무작정 polynomial
추가나 shapebounds확장 금지. v20에서 새term은 아직 도입하지 않았다.
부분 개선은 실제 수행했으나 실제Al항복/피로/검증된interface/물리시간은 미완료.
Ma_phys,Ms_phys,t0 null/초·Hzdisabled. Production/PDE/UI승격gate는 닫혀 있다.

## 보정 재개 v19 — 2026-09-10 (과학 계산 완료; 채택하지 않음)

최신 요청은 "그럼 보정 계속해". Fresh fetch 후 시작 actual HEAD와 origin은
모두 `53ced52f86e687bd3679bb85f09cd647cfc55913`, clean, branch
`probability-pde-solver-v1`이었다. 실제 worktree는 기존과 같은
`aft-pde-bessel-38969ad`; OneDrive는 이동 안내 stub이다. 이후 fresh fetch에서도
remote 변화가 없었다. Main/기존 PDE/UI/calibration/kinetic JSON/Ac 수정 없음.
추가 agent, reset, force push 없음. 최종 commit/origin은 git log/fetch로 확인한다.

### 이번 실제 작업과 이론

v17의 normal force/Haa 물성 오차가 단순 scalar density–rank1 coupling으로
줄어드는지 두 가지 **대체 가설**을 유도하고 실제 재보정했다:

    per atom E1 = D1 ||Q1||² g(x), x=rho_scalar/rho_ref
    rational g=1/(1-z+zx), z∈[0,1]
    power g=x^p, p∈[-1,1]

동일 LJ pair/Bessel 무한급수, per-site density 합산 후 비선형 함수 원칙을
유지한다. 두 g를 동시에 곱하거나 생산 모델에 등록하지 않았다. 원래 계수에서
pristine x=1,Q=0이므로 Hessian은 불변; 보정으로 Haa를 바꾸려면 다른 계수와
tradeoff가 생긴다. 이 점이 추가 항만으로 normal tangent를 고치기 어려운 이유다.

읽을 파일: `solver_v1/COORDINATION_SCREENING_V19.md`, 새 evaluator
`coordination_screening.py`, 실제 fit runners `run_coordination_screening.py`,
`run_coordination_shape_refinement.py`. 독립 validator/reporter는 저장계수 재생이며
optimization 재실행과 구분한다. 결과 root:
`results/fcc111_active_interface/coordination_screening_v19/`.

### 실제 완료 수치

- Rational fixed/free pair CONTROL55 profiles: best z=.093812/loss1255.7014,
  free z=0/loss1254.6575. Power27 profiles: fixed p=-.083511/loss1255.6106,
  free p=.252580/loss1253.7910. 모든 scalar optimizer 수렴, 개선 미미.
- Power fixed-shape 전체를 독립 경로로 실제 재실행: 모든 profile/best dict,
  summary/residual CSV bytes 동일. `reproducibility.json`에 hash와 범위 기록.
- Joint old five-shape90 + power six-shape88 =178 profiles 추가 실행.
  합계260 profiles. 동일104fit,5exactbulk,기존 spectral 제약 유지.
  Old loss1238.611009, power1227.634888(0.886%차이). Old는20nfev budget종료라
  **optimizer-converged 아님**. Power는16nfev ftol종료, gradient global증명 아님.
  Both k_even=12, power p=-1 bound; 임의 문헌 보정값으로 해석하지 않는다.
- Cohesion3.36eV/atom, C11/C12/C44=114/62/32GPa exact residual<8.5e-14.
  이들은 fit 제약이지 heldout검증이 아니다. Power Haa27.78327 vs19.66091
  (41.312%높음), Hxx4.44413 vs4.39359(1.150%높음).
- Joint 독립40jet RMS: energy1.0128, normalforce5.1668, Haa5.9073,
  Hxx8.4471. Force/Hxx는 old joint대조군보다 오히려 나쁨. **물성 채택 불가**.
- Fixed validation178~179s, joint validation178.98s 실제 완료. 각각
  source+두후보×11 tensor cases×9signedstates=297root 모두 검증.
  Joint force residual<1.762e-14eV/L0. 50MPa수직 개구변위 source.000925366Å,
  power.000656211Å(29.09%작음). 4MPa전단 slip source.000330998Å,
  power.000327231Å(1.14%작음). Static unload는 dynamic hold/잔류소성이 아니다.
- Power relaxed fault/saddle .173907/.197371J/m² vs source.150479/.172002.
  첫 fixed-registry ideal traction10.19465GPa vs source12.96957GPa.
  이것은 이상적인 coherent 계면 traction이지 실제 항복강도가 아니다.
- Joint final FD Hessian error3.26e-7eV/L0², per-site direct fine error3.33e-16.
  다섯 dilation×243q×radii12/16 sampled margin 최소.00560971eV/L0².
  전체BZ/모든strain 안정성 증명 아님. 새 unit-channel direct와 전체 harmonic
  displaced-site/analytic 검증을 구분해 rawCSV를 읽는다.
- 저장계수 exact replay 오차0. 전체 shape 포함 exact-bulk tangent SVD:
  old rank10/10 condition228.66, power rank11/11 condition887.53.
  추가항은 conditioning악화. Active inequality cones/통계CI는 포함하지 않는다.
- `final_report/normal_response_comparison.png`는 실제 저장curve에서 생성·확인.
  에너지뿐 아니라 traction/Haa불일치를 그대로 보여준다.

### 발견·수정한 보고 단위 오류

새 joint runner의 CSV3elastic행이 GPa prediction에 변환 전 mode 이름/target/
unit를 붙였다. **Optimizer는 처음부터 올바른 GPa target 사용**. Exporter와
sensitivity metadata를 고치고 두CSV 각3행 metadata만 수정했다. Raw fit JSON,
coefficients, losses 불변. 모든CSV에서 (prediction-target)/scale을 재계산하는
시험 추가. `final_report/*_audited_residuals.csv`는 계수 재생으로 독립 출력한다.

### 회귀와 인계

Targeted 최종31PASS48.19s, app최종34PASS65.27s, smokeexit0/2.257s.
기존 full solver553PASS614.20s 이후 exporter회귀2개를 추가하고 **전체 재실행:
555PASS975.44s(16분15초)**. 최종 모든 tests의 실패/오류/skip=0. 실제 session
종료와 JUnit XML 양쪽 확인. 긴 최종 실행은 sensitivity reporter와 일부 동시 실행.
회귀로그 `.cache/coordination_screening_v19/`, 커밋용 요약은 결과 root의
`verification.json`. Working/staged diffcheck PASS, scoped 파일/JSON finite수와
원본 hash/index bytes 검사 완료. 검증 후 대상 branch로 정상 commit/push한다.
현재 checkpoint의 completed=false는 과거 중간스냅샷이며
완료 판단은 completion/calibration.json이 우선한다. RawJSON/CSV hash보존용
local .gitattributes만 사용하고 전역git설정은 바꾸지 않는다.

### 다음 단계와 아직 불가

새 density-normalization/range변화만으로 물성 오차를 해결했다는 근거 없음.
다음은 normal/registry force와 curvature를 동시에 제약하는 최소 환경항의
독립성/표현력 감사이며, 반복적인 range bound 확대나 임의 mobility/yield튜닝을
하지 않는다. 이 bounded search만으로 전체 analytic family 불가능을 증명하지
않는다. 기존 v18의 선 kinetics/피로 source자료는 별개로 보존한다.

생산 Ma_phys/Ms_phys/t0/실제 초·Hz unavailable/disabled 유지. 실제 항복·피로·
잔류소성·spatial specimen model·A_c·active-interface PDE/UI gate는 열리지 않았다.
LJ/Bessel, 확률식, 기존 static parameter, v18 UI스크롤은 그대로다.

## 연구/사용성 인계 v18 — 2026-09-10

최신 연구 요청은 실제 항복·피로·물리 초/Hz의 검증/보정. 추가 요청은 UI를
켜고, 사용자가 확인한 뒤 해석 버튼을 가리는 긴 설정에 스크롤을 추가하는 것.
Fresh fetch 후 시작 local/origin은 모두
`35b2c5c94f9237bf871eeb6011d1fbdb28726aa9`, clean이었다.
실제 위치는 기존 연구 worktree `aft-pde-bessel-38969ad`, branch
`probability-pde-solver-v1`. OneDrive 경로는 이동 안내 stub이다.
추가 agent/리셋/main변경/force push 없음. 최종 commit과 origin은 git log와
fresh fetch로 확인한다. UI는 실제 실행하여 사용자가 확인했고 이후 수정했다.

### 완료된 범위와 채택하지 않은 것

- `solver_v1/LINE_KINETICS_AND_FATIGUE_VALIDATION.md`에 전위선 kinetics,
  pinned-line 모델의 유도/단위/수렴/hold/물리 한계를 정리했다.
- 실제 Gorman1969 원본 Fig5c(23°C)10개 isolated marker를 읽고7fit/3heldout
  보정을 실행했다. 새 사실은 **실험 전위선 v/전단 계수의 제한된 추정**이다.
  a/s 이동도나 전체 Al 물성 보정이 아니다. Heldout scatter도 크게 남는다.
- 같은 LJ/Bessel bulk Schur 탄성으로부터 작은 pinned-line bow 동역학을
  유도하고64/128/256/512 timestep,16/32/64/128 spatial grid,20tau hold를
  실제 풀었다. 이것은 가역 휨의 음성 대조군이며 fatigue/PDE 새 경로가 아니다.
- Deschanel2017의7개 실제 Al fatigue 실험군을 source/단위/endpoint별로
  정리했다. Delta stress50/62MPa는 full range, 진폭은25/31MPa다.
  Nf17200/5300은 final fracture이고 local opening initiation이 아니다.
- 실제 yield, fatigue, production M_a/M_s/t0는 **여전히 미검증**.
  v17 normal-force/Haa 물성 gate 실패를 뒤집거나 새 항을 채택하지 않았다.
  Production kinetic JSON, LJ/embedding/기존 static calibration/PDE/Ac 불변.

### 재현 경로와 정확한 결과

`results/strength_fatigue_kinetics_v18/`에 약0.75MB의 표/그림/범위를 저장했다.
두 runner는 기존 폴더를 덮어쓰지 않는다. 캐시에 별도 재실행하여6개 dynamics
CSV와 calibration JSON의 byte-identical 재현을 확인했다(그림은 별도 시각 확인).
`scope.json`에는 v17 원본 raw SHA와0K/23°C 불일치, 가상 pin 조건이 명시된다.

1. `run_line_kinetic_benchmark.py` / `data/gorman1969_velocity_benchmark.json`
   - Source DOI10.1063/1.1657472,99.999%Al leading edge/mixed 전위.
   - 원본 PDF SHA `b4536b12a161babc565460facb83b4c64a685bced57cf0f48da9f7cb79df44d3`.
   - Native2490×3305 page5 image를 lossless추출. 각도약간틀어진축을affine역변환.
   - 1e6dyn/cm²=1e5Pa,cm/s=.01m/s. 겹친marker는 추측하지 않고 제외.
   - mu_tau=1.1627233044e-5 m/(Pa s),trainRMSE2.38329m/s,
     heldoutRMSE3.50119m/s=heldout평균의23.62%. 정확한 universal mobility 아님.
   - Reading-only bound[1.040818,1.299157]e-5,LOO[1.104298,1.218419]e-5.
     둘 다 independent confidence interval이 아니다.
2. `dislocation_line_kinetics.py`, `run_line_timescale_validation.py`
   - B_line y_t=T_line y_xx+tau b;T_line=gamma+gamma'';fixed pinned endpoints.
   - G=(96/pi⁴)Σodd1/[n⁴(1+iomega tau1/n²)],tail<=32/(pi⁴N³).
   - v17의 변하지 않은bulk114/62/32GPa와 model0K b=2.8637824638Å를 사용한
     조건부 B=2.4629956698e-5Pa s. 실제296K b/line character 보정은 아니다.
   - 가상pin L=.2/1/5um,R=L,r_core=2b. yield/life에 맞춰 고르지 않았다.
   - tau1=1.913839e-10/3.753056e-9/7.718548e-8s. 이는 선 좌표/가상geometry,
     production time이 아니다. MHz/GHz rolloff는 상수drag 방정식 수치검증용.
   - 0.1/1/25/100Hz에서 거의평형휨. 가장큰100Hzlag도.00274244deg.
   - omega tau=1의공간error .00421475→.00105455→.000263690→.000065926.
     전체60case work/dissipation 상대잔차최대3.07e-15.
   - dtcomplexerror .0236542→.0119804→.00603896→.00304153(약1차수렴).
     마지막20tau hold meanbow4.108e-19m,static대비1.1446e-9; clipping안함.
   - .01/.1/1/4MPa stress sweep은 선형/비선형bow차이를 실제계산.5um,4MPa는
     outergraphfold밖이므로 해를만들지않음. fold38.49/9.81/2.39MPa는
     unknown source/core의 outer-only 연구값, 실제항복/핵생성장벽아님.
3. `fatigue_validation_reference.py`, `data/aluminum_fatigue_validation_v18.json`
   - DOI10.1038/s41598-017-13226-1,99.95%Al,roomT(numericK임의삽입없음),R=-1.
   - stress2군+totalstrain5군,printedrate와2fDelta확인. sourceN/f는 진짜실험초다.
   - PSB100cycle/균열징후600cycle/AE1200이후/파단을 동일 firstpassage로 취급금지.
   - Endpoint/clock/specimen/control/microstructure/convergence의6조건이
     충족되지 않아 현 PDE와수명오차/보정성공을 계산하지 않는다. S-N fit 없음.
4. 기존300K MD의1MiB range를 다시요청했으나20s에응답headers전timeout.
   추가frames없음. 기존1024frame25.575ps의low-frequency적분미해결판정보존.

### UI 수정 및 실제 확인

- `app/scrollable_panel.py`: 독립bindtag 세로스크롤; MouseWheel/PageUpDown,
  Tab으로필드보이기. Combobox휠이값을바꾸거나 Matplotlibzoom을가로채지않음.
- Pre/Solve설정과summary스크롤; 해석/수렴버튼과progress는고정아래row.
- 실제940×620및940×480Windows Tk창을검사하고nativewindowcapture로시각확인.
  최소창높이를480으로낮췄고 두높이×ko/en에서고정버튼가시성을시험했다.
  버튼가시성/입력·결과·zoom보존/언어전환/cleanup테스트추가. UI새physics없음.
- 테스트로그는 `.cache/research_validation_v18/`에보존. 실제python3.13으로실행.
- 최종 확인: **targeted22PASS8.57s / solver537PASS1708.52s /
  app34PASS172.38s / smokeexit0,2.49s / staged및working diffcheckPASS**.
  모든세트skips0. 최종XML은targeted_final.xml/solver.xml/app_final.xml이다.
  이전targeted21/app33는480높이추가전기록이고최종22/34검사로대체했다.
- `validation_manifest.json`에최종검증/범위/보존을집계했다. CSV6개및kinetic
  calibrationJSON 재현동일,rawJSON의index/worktree SHA결합동일도확인했다.
  두MatplotlibSVG의생성trailingwhitespace만정규화했고non-whitespace불변이다.
  연구용lineJSON은productionclock로더에서TypeError로거부됨을별도실행확인했다.

### 다음 실제 연구 과제

1. 전위선 이동계수가 있다고 a/s cell friction을 채울 수는 없다. 검증된
   collective-coordinate metric/slow reduction 또는 적합한장시간MD가 필요하다.
2. 실제yield는 core/source/domain/밀도/누적sweptarea와0.002 criterion이 필요하다.
   새외부mu나가상pin값으로 GPa이상강도를MPa실제강도로rescale하지 않는다.
3. v17계면force/Hessian오차와 defectcore gate는그대로다. 피로 데이터가 있다고
   localabsorption을specimenfracture로바꾸지않는다. 새명시모델/독립검증이 필요하다.
4. 이번 구현은 검증된수치reference와구체적실험비교조건을 추가한것이다.
   솔버완성/생산physicalHz/공간meshUI gate통과로보고하지않는다.

아래 v17 및 이전은 완료된 역사 기록이다. 최신 범위는 위 v18을 우선한다.

## 최종 연구 인계 v17 — 2026-09-10

최신 요청은 "남은거 해봐". Fresh fetch 후 시작 local/origin은 모두
`1818cf9cf630fd094597d919429a7c0499fcb0b5`, clean이었다. 실제 위치는 기존
연구 worktree `aft-pde-bessel-38969ad`, branch `probability-pde-solver-v1`.
OneDrive 경로는 이동 안내 stub이다. Main/reset/force/추가 agent 사용 없음.
최종 커밋과 원격은 이 문서를 포함하는 git log/fresh fetch로 확인한다.

### 완료 범위와 판정

v16 후속인 **normal force/Haa 불일치 감사 및 한 개 rank1 범위 가설**을
실제로 유도·구현·최적화·독립 검증했다. 이 단계는 완료했지만 **Al 계면
물성 채택은 실패**다. 추가 k1를 production/default에 채택하지 않았다.
연구 모듈/반례/수치 검증으로만 보존한다. "전체 analytic family 불가능",
"실제 항복 보정 완료", "솔버 완성"으로 읽으면 안 된다.

- 유도 및 재실행 명령: `solver_v1/NORMAL_ENVIRONMENT_RESPONSE_V17.md`.
- 전체 원시 결과: `results/fcc111_active_interface/normal_response_v17/`.
- 같은 자료로 재생한 비교/그림: `normal_comparison/`.
- 모든7run/종료이유/계수/SVD: `final_report/`.
- 테스트와 수치 검증 요약: `validation_manifest.json`.
- `nested_comparison/`은 먼저 완료한3모델의 중간 비교이며 삭제하지 않았다.
  최종 판정은 위4모델 `normal_comparison` 및7run `final_report`를 사용한다.

### 왜 안 맞는지 확인한 것

1. 두 v16 후보의8개 실제 fixed-shape coefficient control을 먼저 풀었다.
   마지막 perfect Haa27.9164중 rank1기여13.3095,source19.6609eV/L0².
   Haa를 정확히 맞추면 loss596.14->1373.60,Hxx4.7624->2.2565
   (source4.3936). D1=0이면5679.94. Spectral 제약 제거는 optimum 불변.
   이는 그 fixed shape의 tradeoff이며 전체 nonlinear 불가능성 증명 아님.
2. Source 자체의 cutoff5개 및31개 normal state 미분을 독립 검사했다.
   sourceHaa의 step1e-5 FD오차최대9.07e-5, 문제의1.62h는6.21e-9eV/L0².
   기존후보~1.4곡률 차이는 source cutoff/보간/미분 오차로 설명되지 않는다.
   source의 작은 total curvature는 pair/embedding cancellation이다.
   분해는 gauge-dependent이므로 tabulated pair를 LJ로 교체하자는 뜻이 아니다.
3. Rank1만 독립 exponential decay k1로 일반화했다. k1=k_odd면 기존 full
   jet가 정확히 복원된다. Per-atom 합 후 norm,LJ/무한 Poisson–Bessel 보존.
   Rank3와radial Eg의k_odd는 그대로다. Affine inversion에서Q1=0이므로
   bulkC는 불변이나 interface/finite-q O(q^4)는 달라져 모두 재검증했다.

### 실제 보정/대조군 — 종료 이유를 숨기지 말 것

- 5exactbulk+104interfacefit+미리 정한48heldout. 과거v16검증36개는
  development로 재분류했다. 새48개도 이제 검사했으므로 미래blind가 아니다.
  동일state/jet중복0개. 10%/0.01J/m²/0.1eV/L0²/250MPa척도는 측정오차가 아니다.
- v16양의pair(u=.01719039249,v=.10069343124)를 유지한 실험은 **CONTROL**,
  새 물성2개가 아니다. L0=2.8637824638Å,E0=1eV,A_atomic_cell=sqrt(3)L0²/2.
- 7run에 실제365coefficient profiles:4nonlinearstarts중360회+고정형상5회.
  Tied96회/nested72회/독립k1=3,9각96회. 모두 양의LJ이고 failedprofile0.
  **4nonlinearstarts 모두 평가예산 종료**. Profile-only 성공은 shape수렴이 아니다.
- 같은v17자료의 loss/heldoutRMS:
  oldv16재생1335.56/4.8765;
  tied1256.9468/4.9991;
  nested1256.8219/5.0089;
  독립2시작점1263.2421/4.9570.
  Old596.14는 이전자료의loss라 직접 비교하면 안 된다.
- Nested shape=[2.3172001610,6.2719556784,12,436.6576096,6.2801140650].
  독립시작점best=[2.4841497658,6.3099824750,11.9823531066,485.2569646,6.3310019895].
  k1가k_odd근처로 복귀. 추가자유도로 nested fit개선~0.01%,검증은 약간 악화.
  최종budget trial과optimizer의lastaccepted가 다를 수 있어 둘 다 저장했다.
- 같은shape에서 pairCONTROL해제: tied1255.4174/nested1254.6575/
  separate1258.2939. 정상곡률 불일치가 해소되지 않았다. 이들은 nonlinear
  joint재최적화나 독립fullvalidation을 한 후보가 아니라1회 coefficient control.
- k_even=12경계에 대한15/18고정shape검사 loss1293.02/1401.78로 악화.
  이것만으로 더 큰범위의 전체최적화도 불가능하다고 말하지 않는다.
- Exact equality tangent rank/cond: tied7/125.84,nested8/119.67,separate8/122.82.
  Bounds/cones/목표불일치를 고려한CI가 아니고 uniqueAl물성증명이 아니다.

### 물리 및 전산 검증 결과

- Source perfectHaa19.6609대비 tied27.7419/nested27.6537/separate27.8577,
  +40.7~41.7%오차. Nested fault/saddle energy=.175144/.198270J/m²,
  source=.150479/.172002. SaddleHaa24.3787vs12.7173: 힘/곡률 gate 미통과.
- Nested heldout RMS(energy/force/Haa/Hxx)=.8817/4.6549/5.8300/6.6274.
  Energy도8/12만척도이내이므로 RMS<1을 전부PASS로 읽지 않는다.
- 세후보 모두 독립243q×radii12/16+연속국소검색,6dilation×254q×2radii
  완료. 총36dilation기록은 tail이상양수, whole-zone안정성증명은 아니다.
- 세후보의 Hessian FD오차(step1e-5)최대8.63e-8,reciprocalH변화4.48e-13.
  Direct per-site unit-energy합은최대3.06e-16. 가장 작은FDstep에서 roundoff
  증가도 원시CSV에 그대로 저장했다. Canonical finite-neighbor cutoff 아님.
- 11개MPa tensor×9signedstate×(source+candidate)×3회=실제594states.
  Source99states가 후보마다 반복됐음을 명시. 모두root검증통과.
  Nested50MPa수직개구변위오차-28.76%,4MPa전단-2.05%.
  Tied는-28.98%/-3.11%,separate는-29.28%/-5.68%.
  Staticunload최대3.62e-15L0. PDE/zero-stressdynamic hold/소성검증이 아니다.
- 후보opening극값은[h,5h]의129/257혼합격자에서1개,source는5개.
  후보첫traction~10.27GPa/source12.97GPa는 이상 coherent 계면이지
  실제 시편 항복이 아니다. Source의 뒤쪽 작은 음의lobe도 숨기지 않는다.
- 기존source300K MD1024frames/576atoms/plane/12classes,실제boxstretch로
  고정물성비교. Nestedvariance오차(-12.41,-22.26,-14.84)%,
  tied(-12.62,-23.16,-15.83)%,separate(-13.59,-25.29,-18.16)%.
  새MD실행/운동보정은 아니며0Kloss에 넣지 않았다.

### 최종 코드 실제 회귀와 보존

**targeted29PASS129.63s / solver528PASS1102.12s /
app31PASS214.92s,skips0 / smokeexit0,최종0.98s / diffcheckPASS**.
첫526PASS2099.85s는 마지막2테스트추가전 실행이고 최종528완주로 대체했다.
XML은 `.cache/normal_response_v17/` 로컬로그, 집계는 trackedmanifest에 있다.
동시계산 중 측정한walltime이며CPU성능벤치마크가 아니다.
Productiona0=.7713438268704838,kappa=86.29296488740997불변.
App/PDE/과거staticparameter파일은 수정하지 않았다. Main도 건드리지 않았다.
v17 JSON은raw-byteSHA256으로 서로연결돼 있으므로 해당결과폴더에만
`.gitattributes`의 `*.json -text`를 설정했다. Git줄바꿈변환으로 원시검증
바이트가 바뀌는 것을 막는다. 과거폴더/계산값은 바꾸지 않았고 raw/index
동일성을 별도확인한다. SemanticJSONhash나 타플랫폼재실행동일성 인증은 아니다.

### 다음 단계 — 이 실험을 다시 무작정 반복하지 말 것

1. 추가range만으로 해결되지 않은 **total nonlinear normal environment**
   force/Hessian형상을 항별budget/equalitytangent와 함께 분석한다.
   새항은 명시적유도/최소자유도/검증을 먼저 하고LJ/Bessel틀은 보존한다.
   같은가족continuation을 더 하면 예산/대조군을 기록하고v17heldout은재분류한다.
2. Static재료gate를 통과한 뒤에만 nonlinear/vector/normal-relaxed discrete
   core와finite-source/loop 및collective PMF로 실제강도 문제를 진행한다.
   Source자체가13GPa인uniforminterface를 맞춰도 실제항복이 되는 것은 아니다.
   A_c/임의line길이를 곱해activationenergy를 만들지 않는다.
3. Ma_phys/Ms_phys/t0와실제초/Hz는 여전히 unavailable/disabled다.
   Tensor는연구3×3traction/staticstate경로이고 production다축PDE/UI연결아님.
   실제yield/fatigue/A_c/생산solver완성/UI모델링-meshing gate는 아직미통과.
4. 완료 계산/부정적 결과는 위원시파일과문서에 보존했다. 새 kinetic자료나
   새로운수학적환경가설 없이 같은몇개숫자를고쳐 성공으로보고하지 않는다.

아래 v16 및 이전은 완료된 역사 기록이다. 최신 판정은 위 v17을 우선한다.

## 최종 연구 인계 v16 — 2026-09-10

최신 요청은 "남은 작업 ㄱㄱ". Fresh fetch 후 시작 local/origin은 둘 다
`6afa55ebe60e2808523219e03da4f4eb5ea276c4`, clean이었다. 실제 worktree는
이전과 같은 `aft-pde-bessel-38969ad`, branch는 `probability-pde-solver-v1`.
OneDrive 경로는 저장소 이동 안내용 stub이며 실제 연구 파일이 아니다.
추가 agent는 생성하지 않았다. Main/reset/force 사용 없음. 최종 커밋/원격은
이 문서를 포함하는 git log와 fresh fetch로 확인한다.

### 실제 완료 및 보존 경계

- v15의 계면 Haa 오차를 항별로 분해한 뒤 rank2와 결합 Eg의 비선형 응답에
  **한 개 무차원 shape**를 추가한 별도 연구 가설을 유도/실행했다.
  `D I/(1+alpha I/N^2)`이며 각 원자의 무한 환경을 먼저 합한다. Cubic Q=0의
  조화 탄성 및 모든 finite-q Hessian은 보존되지만 비cubic 상태까지 불변은 아니다.
  LJ, Poisson–Bessel, 기존 보정 파일/생산 기본값/확률 이론은 바꾸지 않았다.
- 유도/검증/재실행: `solver_v1/EVEN_ENVIRONMENT_CALIBRATION_V16.md`.
  원시 결과: `results/fcc111_active_interface/even_environment_v16/`.
  합산 표/그림: 그 아래 `final_report/`; 자동 요약: `validation_manifest.json`.
- v15에서 이미 본 holdout은 v16 development로 재분류했다. 5개 exact bulk,
  60개 fitted interface, 36개 새 off-grid/off-path heldout(12상태×3관측량).
  후자는 loss에서 제외했지만 이제 이미 검사한 검증 자료다. 계속 blind라 하지 않는다.
- **7개 run / 420개 실제 성공 coefficient profiles**를 저장했다. 여기에는
  한 개 fixed-shape 재프로파일도 있으므로 7개 독립 재료/완료 최적화라 하지 않는다.
  Baseline loss31834.68/heldoutRMS19.725 → free saturated616.16/4.0069 →
  positive-pair refinement596.14/4.1397. 마지막은 train은 개선, heldout은 악화.
  모든 종료 이유/failed trial/예산소진도 보존했고 global optimum은 주장하지 않는다.
- 초기 진행 기록과 달리 free saturated에서 **20개 양의 LJ profile**을 찾았다.
  나머지172개의 zero-pair closure와 더 낮은 loss553.91은 비허용 진단이다.
  전부 비허용일 때도 runner는 진단 결과를 저장하고 material loader는 거부한다.
  Old-pair 고정 대조군은 v15 값, 마지막 positive-pair section은 free saturated의
  양의 pair를 고정했다. 이 두 정확 제약은 Al 측정값이 아닌 CONTROL이다.

### 물리적으로 개선된 것과 아직 실패한 것

- 정확한 bulk cohesion3.36eV, C11/C12/C44=114/62/32GPa, L0=2.86378246Å 유지.
  원자 면적은 sqrt(3)L0²/2, 에너지는 eV/interface cell. A_c와 무관하다.
- 정상 이완 fault energy[J/m²]: source0.150479, baseline0.373858,
  saturated0.167822, 마지막0.163016. Full-vector saddle: source0.172002,
  baseline0.421187, saturated0.190537, 마지막0.186420.
  Free saturated의 비적합 energy12개는 모두 선언한10%/0.01J/m² 척도 안이다.
  따라서 에너지 수준의 개선은 확인됐다. 이전 "아무 보정도 안 됨"과는 다르다.
- 그러나 perfect Haa는 source19.6609 대비 마지막27.9164eV/L0²(+42%),
  saddle Haa는12.7173 대비24.1127(+90%)다. Free saturated의 heldout force/Haa
  정규화 RMS5.0716/4.7053. **힘/곡률까지 맞춘 Al 계면으로 채택하지 않는다.**
- 개구 traction 극값을 near-h 로그격자+전체 균일격자129/257로 재검사했다.
  Baseline의 -4.284GPa 포켓은 두 양의 포화 후보에서 검출되지 않았다([h,5h]).
  최고 traction10.57/10.49GPa, source12.97GPa는 **이상 균일계면 견인력**이지
  실험 항복이 아니다. Source 자체의 뒤쪽 -0.786GPa lobe도 숨기지 않았다.
  큰 alpha의 비균일 극한/좁은 tangent 가능성을 수식과 synthetic test로 검사했다.
- 독립243q×radius12/16 및 연속 국소 검색, 고정재료6dilations×254q×2radii를
  baseline/free saturated/마지막 후보에서 완료했다. 모든36개 dilation 기록은
  analytic tail보다 양수지만 whole-zone stability 증명은 아니다.
- Direct/reciprocal, analytic jets와3단계FD, reciprocal tolerance 비교 완료.
  첫 saturated의 최종 H 차분오차8.00e-8, 직접합 unit-energy 오차1.04e-17.
  새 alpha를 적용한 per-site direct sinusoidal energy로 finite-q도 독립 검증했다.
- 실제297개(source/baseline/saturated)+198개(source/마지막) 정적 tensor 시나리오:
  축별 수직/전단/혼합응력과 정적 unload 실행. 첫 saturated는4MPa 전단 변위
  source오차-0.47%, 50MPa normal변위-29.7%. PDE/동적 hold/실제 소성 결과 아님.
- 기존 실제 source300K MD1024frames도 고정재료로 비교. 첫 saturated의
  normal/slip/transverse variance오차(-14.4%,-20.9%,-13.3%), 마지막은
  (-14.3%,-26.9%,-19.9%). 0K loss에 넣지 않았고 kinetic fit도 아니다.

### 독립 재현한 수치 결함과 최종 테스트

- 실제 fixed_saturated QP에서 K3=-1.714886e-9인데 whitened KKT3.55e-15인
  변환 오차를 `boundary_reproduction/actual_qp.npz`로 bitwise 재현했다.
  Constructor의 거부는 옳았다. 계수를 clip하지 않고 정확한 활성 경계에서
  다시 풀어 원래 full-space KKT/feasibility를 재검증하도록 수정했다.
  Fixed-shape 재프로파일 K3=0, prediction변화5.90e-13, exact잔차2.04e-13.
  수치 수정이 그 old-pair 후보의 재료적 오차를 고치지는 않았다.
- 최종 코드로 **target40PASS39.02s / solver518PASS585.55s /
  app31PASS63.25s, skips0 / smokeexit0,1.01s** 실제 완료.
  `.cache/even_environment_v16/*_final.xml`은 로컬 검증 로그이며 개인 절대경로는
  커밋하지 않는다. 앞서 중단한 full run은 PASS가 아니며 최종 완주만 센다.
  Smoke의 기존 a0=.7713438268704838, kappa=86.29296488740997 불변.

### 다음 작업 — 처음부터 다시 하지 말 것

1. 우선 잔여 normal force/curvature 오차의 항별 원인을 두 양의 후보로 비교한다.
   작은 bulk/energy residual만으로 surface Hessian 채택을 선언하지 않는다.
   새 물리 항은 기존 가족 실패/식별성 검토 뒤 최소한만. 이미 본 heldout 재분류 필수.
2. Shape optimization은 예산종료/경계 및 pair tradeoff가 남는다. fixed control을
   새 물리 정보라 하지 않는다. Exact-tangent cond146.4/122.5도 유일성/CI가 아니다.
3. 재료 검증 뒤에만 finite core/source와 올바른 collective-coordinate PMF/kinetics로
   이어간다. 이번 새 후보로 옛 core 결과를 재생 없이 재해석하지 않는다.
4. Ma_phys/Ms_phys/t0, 실제 초/Hz, A_c, 실제 항복/피로 및 생산/PDE/UI gate는
   여전히 미완료다. Mesh UI 재설계는 solver 검증과 사용자 확인 전 진행하지 않는다.

아래는 완료된 v15와 그 이전 역사 기록이다. 최신 판단은 위 v16을 우선한다.

## 최종 연구 인계 v15 — UTC 2026-09-09 22:50

**이 블록은 v15 종료 당시 인계다. 아래22:43/21:28/20:38도 역사 기록이다.**
약4시간의 외부자료·단위·보정 감사에서 실제 실행과 검증을 마쳤다.
연구/수치 구현의 검증 완료와 Al 계면/실제강도/kinetics 채택은 다르다.
후자는 아래 이유로 아직 통과하지 못했다. 기존 이론/production을 바꾸지 않았다.

- 최종 **21 runs /1770 coefficient profiles**; 마지막112profiles는 frozen
  code의 독립 재최적화다.532.31s 실제 실행에서 원래 strain_stable_cross와
  decays/coefficients/predictions/residuals/Cubic/loss 차이 모두0.
  독립243q×radius12/16,stationary/curve/SVD를 다시 실행(46.22s)했다.
  최신 합산 결과는 `interface_development_v15/release_report/`.
  이전 `final_report/`는20run snapshot이며 curvature term budget도 그곳에 있다.
- 코드 고정 후 **target73PASS18.73s / solver504PASS530.11s /
  app31PASS66.83s,skips0 / smokeexit0,1.21s**. 기록은
  `results/public_aluminum_validation_v15/validation_manifest.json`.
  보고 이후 Python 수정 없음. 마지막Gitdiffcheck는 커밋 시 다시 수행한다.
- Source thermo log의 6prefix/101whole coarse samples를 추가 확인.
  Source300K는 setpoint, prefix snapshotmean299.108K,압력+8.464MPa.
  500ps scalarlog를500ps좌표분석으로 부르면 안 됨. 실제좌표분석1024frames만.
- 원자면적·peratomvolume·lineenergy·cellenergy·traction/mobility units를
  모두 명시했다. 특히 coherentNp 모드에서 percellG*를 쓰면
  M*=t0*Np*E0*M_N/L0²와 kT*=kBT/(Np*E0)가 **동시에** 필요하다.
  소스576atoms/plane을 A_c나 현재PDE 독립영역 개수로 대체하지 않는다.
- Production thermal audit:TwoRowLJ kT=.02는 무차원이라Kelvin미정.
  Reducedhybrid의 kT_ev=.02는232.09036K에 상당하며300K가 아니다.
  이를 조용히실온으로고치거나 M과동시에재조절하지 않았다.

### 최종 물리적 결론

1. 정확 bulkC114/62/32GPa,cohesion3.36eV,scale은 보존/검증됐다.
   그러나 최신 안정 후보도 source계면 heldoutRMS71–79:재료 채택 불가.
   DirectregistryHaa source36.32 vs369–424eV/L0². 부정확한 수치가 아니라
   같은단위 에너지항들의 잘못된 지형이며, 보정 후 물리검증 실패다.
2. 더 낮은loss 후보는 finite-q 또는 작은dilation에서 불안정으로 탈락.
   QP 결함은 수정했지만 그 수정으로 물리 실패를 숨기지 않았다.
3. 실제MD,실험wire,문헌line/core의 값/단위/좌표를 구분했다.
   Shortlag와저주파적분 둘다 현재 a/s시계 보정으로 채택하지 못한다.
4. 유한전위원의 외부 탄성기하학으로 MPa 규모가 유도되지만 현재 pin/core
   가설값이다. 실험yield/작동segment/속도/잔류소성/피로 검증으로 채택 금지.
5. 전단/방향 tensor는 연구static/vector/reflectingSG까지만 연결돼있다.
   Production/UI는 여전히scalaraxial P(a,s). DefaultTwoRowLJ 불변.
6. Ma_phys/Ms_phys/t0는없고 실제PDE 초/Hz 사용불가. A_c도미보정.
   Mesh/UI gate는 닫힌 상태이며 사용자승인없는 UI재설계를 하지 않았다.

### 재개할 때

`V15_REPRODUCTION_AND_STATUS.md`에 실제재최적화/오프라인MD재생 명령이 있다.
다음 연구는 정확한 안정 계면/core 에너지와 유한 활성화좌표, matchedPMF/
kinetics를 우선한다. 같은holdout을 새loss에쓰면 development로 재분류할 것.
이미 실패한 family를 설명없이 다시 돌리거나 알고리즘실패와 물리실패를 섞지 않는다.
이번 동안 유지된 브랜치는 `probability-pde-solver-v1`, 시작 HEAD d7ef95cf...
이다. 이 문서가 담긴 커밋과 최종origin상태는 git log/fetch로 확인할 것.
main/reset/forcepush는 사용하지 않았다. 상세 scientific기록은 아래와6개v15문서.

## 최신 진행 인계: UTC 2026-09-09 22:43 — 검증 완료, 독립 재최적화 진행 중

**아래21:28 및 이전 진행 기록은 역사적 checkpoint다. 처음부터 재시작하지 말 것.**
사용자 최신 우선순위는 전단/방향 연결 상태를 먼저 브리핑하고 약4시간 동안 실제
외부 자료로 보정/단위/척도를 검증하는 것. 브리핑은 이미 전달했고 연구는19:15경
시작,23:15경까지 마무리 목표다. 새 에이전트는 생성하지 않았다.

Git: 시작 local/fresh-origin 둘 다
`d7ef95cf7d8c814053d0caa741b2f2b1e3f4d338`, clean,
`probability-pde-solver-v1`. 실제 worktree는 OneDrive 밖
`aft-pde-bessel-38969ad`다. main/local=`80cacb4180dbfbcd36a2964270703bc6cf1653ec`,
origin/main=`c43d8e096f2a329c7cdfad31e546c70fb36fe592`. Main 수정/reset/force 없음.
모든 v15 변경은 이 턴에서 만든 것. **현재 아직 commit/push 전**.

### 실제 완료

- 최종 코드 고정 후 targeted **73PASS18.73s**, full solver **504PASS530.11s**,
  app **31PASS66.83s**, skips0, desktop smoke **exit0/1.21s**,
  a0=.7713438268704838,kappa=86.29296488740997. `git diff --check` 통과.
  `.cache/public_aluminum_v15/`의 `targeted_final_confirmed.log`,
  `solver_frozen.log`, `app_frozen.log`, `smoke_frozen.log`에 실제 결과.
  이 결과 이후 연구 Python 코드는 바꾸지 않았다(문서/데이터 보고만 추가).
- 현재20개 실제 최적화 run/1658 성공 coefficient profiles를
  `interface_development_v15/final_report/`에 합쳤다. 이 수는 독립 재료 수가
  아니다. 실패 profile/optimizer 종료/모든 정확 residual과 SVD 보존.
- QP zero-dual near-active 및 Gram condition 제곱 문제를 같은 검증 tolerance로
  수정. 단일 probe:32cuts,-5.93e-9 실패→1cut,-1.59e-12,KKT1.87e-12.
- 낮은loss1084.85 basin은 독립 finite-q에서 Hmin=-5.08215,tail=.00205967로
  실제 불안정. 그 q를 추가한 spectral 후보도 실제MD dilation4.065/4.05에서
  불안정. 수치오차/whole-zone proof/재료 실패를 서로 구별하고 모두 기록.
- `.99/1/MD/1.01` strain을 포함해 다시 보정한 결과:
  strain_stable_quartic62profiles396.01s, loss1251.068874,heldoutRMS78.29824;
  strain_stable_cross112profiles720.95s,loss1128.234352,RMS71.20608;
  strain_cross_bloch95profiles503.31s,loss1957.053393,interfaceRMS78.90034,
  BlochRMS8.70063. 세 후보 모두 독립6stretches×254q×radius12/16 통과.
  Reference continuous searches도 양성. **계면 오차로 모두 미채택**.
  exact-tangent conditions311.58/920.60/207.26, 전체zone증명/CI 아님.
- 실제 same-unit heldout direct(.39L0)의 Haa budget:source36.3168 vs
  candidates403.575/369.201/424.413eV/L0². Pair/scalar/angular 합의 큰오차로,
  단위변환/strain plot문제가 아님. `heldout_curvature_term_budget.csv`.
- 실제 source1024MD frame의 short-lag + zero-frequency integral을 모두 검사.
  후자3.15ps eigen=-.01078,-.00676,-.00419ps vs floor.03825ps: **unresolved**,
  음의물리마찰/양수window선택으로 M 만들지 않음. Offline projected NPZ로
  실제 재실행했으며 `decision.json` 문자열까지 동일했다.
- 고정재료의 실제 MD box covariance: v14합산normal이잘맞아도 individualmode
  -25%~+37%오차.3후보MD normalvariance -8.89/+7.94/+30.94% 등 공개.
  Fixed density gauge/L0/coeff 유지, positive Hessian만 covariance역산.
- Source thermo log도 실제 확인:projection내6coarse samples 평균T299.108K,
  pressure+8.464MPa;전체500ps log101samples T299.930K,+9.832MPa.
  NVT는zero-pressure아님. 완전log가500ps좌표분석을 뜻하지 않음.
- 외부wire692raw traces/29specimens,643selected apparentactivation rows/30labels,
  물리line-drag6행,Lu2000 PN core/Peierls8행/units/DOI/SHA 보존.
  686/692/23excluded,114/112um불일치와raw/processedoffset을 숨기지 않음.
- Actual mixed-load saddle84cases, 추가cubic tensor **297staticstates70.86s**:
  x/y/z normal,xy/xz/yz shear,localtwo shear 및mixedtensor를 실제로 실행.
  LocalH positive,maximumstaticreturnerror3.72e-15L0. PDE/dynamichold 아님.
- SameLJ/Bessel elasticSchur의 finite pinned-line min/saddle/barrier를 유도.
  54HYPOTHETICAL source geometry/loadcases, angular/quadrature/FDrefinement.
  L.2/1/5um,rcore2b의critical59.139/15.079/3.666MPa. 실제yield/rate 미보정.
- 정확unitledger:1eV/cell=2.255795J/m²,tractionunit7.876979GPa(**yield아님**),
  1MPa=.0001269522eV/reducedcoordinate. CoherentNp mode는 drift와diffusion을
  함께 rescale해야한다는 유도 추가. Np=576은 source원자수,A_c가 아니다.

### 지금 실행 중 / 다음 처리

- session **18427**: 최종동일코드에서 `frozen_reoptimization_cross`를 실제
  2start/24nfev로 새로 최적화. 단순 best-vector replay 아님.22:42경47profiles.
  `.cache/public_aluminum_v15/frozen_reoptimization_cross.log` 확인.
- 완료 후 이전 strain_stable_cross와loss/coefficients/predictions 비교,
  새 `validate_tail_calibration ... --grid-step .08` 독립검증 실행.
  보고 directory는 새 이름으로 만들고20/1658 count를 actual최종수로 갱신.
- Python코드가 바뀌면 해당tests/full다시. 현재 data/docs추가만 있으므로위
  frozenfulltests가 해당코드기준이다. 마지막diffcheck/문서review는 다시 할 것.
- Beforecommit 관련파일검사. Beforepush freshfetch→topic/main확인→normalpush만.
  아직완료하지않은 Git결과를 완료라쓰지 말고 실제SHA/remoteHEAD를 보고.

### 절대 바뀌지 않은 adoption 상태

전단/다른방향은 연구 tensorprojection/vectorstatic/reflectingSG에만 연결.
Production `desktop_ui -> solver_adapter -> build_energy_model -> P(a,s)`는
scalaraxial이고 독립shear/direction입력은 미연결. DefaultTwoRowLJ 불변.
과거 matched_v3 연구SG를 최신v15 후보PDE처럼 말하지 않는다.
**Ma_phys/Ms_phys/t0 unavailable, physicalseconds/Hz disabled**.
실제강도/유한 core·source state/계면 정확성/공간A_c 미해결. 원자 LJ/Bessel,
확률/흡수/registry이론과 기존보정파일,productionsettings는 보존했다.

문서: `V15_REPRODUCTION_AND_STATUS.md`(실제실행재현),
`INTERFACE_CALIBRATION_DEVELOPMENT_V15.md`, `PUBLIC_ALUMINUM_CALIBRATION_EVIDENCE.md`,
`LOADING_CONNECTIONS_AND_CALIBRATION_V15.md`, `ZERO_FREQUENCY_MOBILITY_AUDIT.md`,
`FINITE_SOURCE_BARRIER_DERIVATION.md`. 위의미완료/검증실패를 지우지 말 것.

## 추가 진행 백업: UTC 2026-09-09 21:28 — 아직 최종 완료 아님

사용자 요청 약4시간 중19:15경 시작,23:15경까지 작업 목표. **재시작하지 말고
이 진행을 이어갈 것.** 아래20:38 기록보다 최신. 시작 HEAD/fetched origin
`d7ef95cf7d8c814053d0caa741b2f2b1e3f4d338`, 브랜치 동일, main/reset/force/agent 없음.
모든 v15 변경은 아직 uncommitted이며 실제 작업 worktree는 OneDrive 밖이다.

새 완료 사항:

- 기존 후보의 물리 척도를 유지한 채4.065/4.05배 cubic dilation에서 MD 비교.
  `isotropic_bulk_validation.py`: density gauge/L0/coefficients를 재설정하지 않음.
  H_A=-R/sqrt(x)+gg/(4x^1.5), H_B=2R, H_C=4(x-1)R+2gg;
  cross항 harmonic=(x-1)H_D3를 포함. 직접 site-energy/진폭/tail검증 통과.
  `fixed_material_modes/`: v14합산normal RMS4.00332e-13m vs MD4.00226e-13,
  shear9.84199e-13 vs1.00902e-12/9.64049e-13. 그러나 normal개별mode오차는
  -25%~+37%; block편차도커서 시간/완전한재료검증으로 채택 안 함.
- periodic plane sum=0은 stationarity증거가 아님을 명시하고 plane별 block평균
  진단 추가. `kinetic_plane_stationarity/`는 실제1024frame재실행 완료.
  Covariance/kinetic decision불변, M_a/M_s/t0 여전히null.
- `finite_source_barrier.py`: 같은 Schur/탄성 line모형의 minorminimum/overhanging
  saddle 및 양의 유한DeltaG, -dDeltaG/dtau=bDeltaArea를 유도. Second variation
  p*integral[(eta_theta)^2-eta^2]로 index0/1 검증. 원자core검증 아님.
  `finite_source_refined/`:54hypothetical spans/core/loadcases, 각도32/64,
  quadrature/stressFD refinements 완료. L.2/1/5um,r_core2b의 critical59.139/
  15.079/3.666MPa. **실측source길이/항복/확률로 채택금지**. 별도유도문서 작성.
- Lu2000 PRB62,3099 원문TableII 실제시각확인/PDFSHA보존.
  `data/aluminum_core_validation.json`:8DFT/EAM입력PN core/Peierls행.
  1meV/Å³=160.2176634MPa. DFT입력screw256.35/edge3.204MPa는 PN저항,
  실험yield/normalcellmobility가 아니다. EAM은Ercolessi–Adams, Mishin아님.
- Exact-stage sensitivity는5bulk등식 전체를소거. `validation_exact_tangent/`.
  quarticrank8/cond311.58,mixture9/6.06e6,cross9/921.72,
  rational9/3006.47, jointnonexactC11/3497.86. Inequalitycone/confidence미포함.
- `joint_bulk_interface`114profiles387.81s, loss311.22124/holdRMS14.82148지만
  C83.425/70.943/34.248GPa로부적합. `rational_development`96profiles258.22s,
  loss1146.17112/hold79.79577. 둘다독립검증완료, 채택안함.
- `completed_report/`는 지금까지11runs의실제fit/heldout/parameters/SVD/figure.
  `report_interface_development.py`는재최적화아닌보고. 추후완료run추가해서
  **새report디렉터리**로갱신할것. `scaling_report/`actualmode/barrier그림완료.

추가 실제 QP 수치결함 감사:

- 기존Gram(A A^T) 제약면투영은조건수를제곱하여 거의평행한spectralcut정밀도상실.
  특정실패probe k=(3.08175732286109,3.73956565371123,6.84896287148422)를
  재현한 `precision_probe_before/`:32cuts실패,margin-5.93447e-9.
  directaffineSVD+nullspace로같은probe는1cut,margin-1.59333e-12,KKT1.87e-12.
  원래sameprimal/KKT/spectralroundoff/tail기준유지. 독립near-parallel test추가.
- 단, rawSVDface는over-screenedface일수있으므로 primal뿐아니라dualstationarity도
  검사하여받아들여야함. 이추가guard를현재코드에반영. 첫 `grid_direct_svd`는
  이guard전의실패기록이며성공으로보면안됨.

현재 돌아가는 실제 세션(로그 먼저 확인; 이미끝났으면 결과수집):

- `77672`: `grid_svd_certified` (guard포함최신QP) profiles53근처/loss1084.90,
  아직completed아님. 이전막힌grid출발점이더낮은lossbasin으로진행중.
- `57259`: `rational_direct_svd` (직접SVD이지만guard전import),profiles136근처,
  best1146.17. 완료후 실제실패/종료조건을보존하고필요시최신guard로재실행.
- `49805`: app regression, `.cache/public_aluminum_v15/app_release.log`.

실제시험결과:

- targeted63PASS18.51s (최종SVD수정/새coretest이전).
- core-unit2PASS.65s, SVD/cross target14PASS9.70s.
- full494PASS520.77s **최종수정이검사실행중생겼으므로최종회귀로간주하지말것**.
- 최종target/fullsolver/app/smoke/diffcheck/문서/commit/push는계속수행해야함.

다음: 위보정종료후새best의독립validation,summary갱신. Fullmaterial가진짜검증되지
않으면PDE/UI에추가하지말것. 사용자처음요청한연결경로브리핑은이미전달:
전단/타축normal은연구tensorprojection/vectorstatic/reflecting연구SG까지만,
productionUI/PDE는scalaraxial P(a,s),독립전단/방향미연결. PhysicalHz불가.

## 진행 중 추가 인계: v15 실제 자료·보정 감사 (UTC 2026-09-09 20:38)

아래 초기 진행 기록은 **역사 기록**이며 다음이 더 최신이다. 사용자가 약4시간
지정한 작업은 UTC19:15경 시작, 대략23:15까지를 목표로 계속 중이다.
현재 모든 변경은 아직 uncommitted. 시작/fresh origin d7ef95cf... 보존.
실제 scientific worktree/branch 동일, main/reset/force/agent 사용 없음.

### 완료된 실제 결과

- 1024frame 공개300K Al MD 감사 완료. Source와 동일 EAM의 Bloch Hessian을
  독립 sinusoidal-displacement energy로 검증하고 정확한576atoms/plane,
  12periodic plane mode normalization을 유도했다. 정적 예측/실제MD RMS:
  normal4.18743e-13/4.00226e-13m, shear8.90817e-13/1.00902e-12m,
  transverse8.90817e-13/9.64049e-13m. Scale fitting 없음. Sample/anharmonic
  차이를 남기며 이로 kinetic t0를 얻었다고 하지 않는다. `plane_normalization/`.
- MD8192frame 확대는 두 번 네트워크 실패. `zenodo_range_8192`에는128frame
  checkpoint만 있으며 **완료 아님**.1024frame 결과/원본은 그대로 보존.
- Wire 최신 결과는 `wire_relaxation_audited/`:692parsed,0parsefailures,
  23source `Not considered/`. Published686과 같다고 하지 않는다. D114/112um
  불일치 및 중복 timestamps를 유지. 원자료 코드/notebook 실행 없음.
- 실제 final-selected643 activation rows/30labels 파싱. Raw29와 다른 이유,
  source analysis T293K,b.286nm, signed inverse-error fields를 그대로 기록.
  대표 Al20r10의 apparent stress derivative1021.36–2988.90b³를 source slope에서
  재계산(relative4e-10). `apparent_activation/`. 이것은 A_c/임의체적이 아니다.
- Source+v14+rejectedv15 총63 actual mixed-load saddle/minimum 검사를 수행.
  Tn=-25/0/+25MPa, shear0–10MPa. Envelope derivative
  -dDeltaG/dT=A_atomic_cell L0(q_saddle-q_min)를 독립FD로 검증.
  Source uniform-cell shear derivative.321647b³, barrier.0762491eV.
  실험 apparent derivative보다 수천 배 작다. 임의cell곱셈 금지; finite collective
  event/core/source state가 없는 것이 실제강도 연결의 핵심 한계.
- 실제 dislocation line drag 원문을 검토하고6source records를 저장.
  Olmsted2005MD/Ercolessi–Adams, Gorman1969experiment를 source별로 분리.
  B_line[Pa s]는 a/s mobility[m²/(Js)]가 아님. Conditional core-friction
  M_s=integral(s')²dx/(B_line A_atomic_cell) 유도/시험, 가정미검증으로 명시.
  M_a,t0 생성/생산JSON 변경 없음.

### v15 보정 실행 결과

- Baseline continuation154profiles 완료, loss1251.068873, holdoutRMS78.30.
- Declared100shape grid85verified/15numericfailed. Sign constraint 해제도 best
  loss1251.068874: 그 shape에서 residual 원인이 sign 제한은 아님.
- Positive full-density mixture1weight,274profiles/761s, w→.0001234,
  loss1251.096986, independent validation 완료. 개선없어 미채택.
- Source staticBloch 7fit+7heldout 추가보정41profiles,
  loss2300.28957(다른loss), independent validation36.58s. Interface 미채택.
- One analytic site term E_xI(x-1)||Q3||²를 별도 최소연구 확장으로 유도/시험.
  [[C,E/2],[E/2,K3]]PSD를 실제profile에 제약. Pair/원자당count불변.
  초기2runs는 active-face numerical failure, 잘못된 성공으로 부르지 않는다.
- 실제 QP 결함 수정: zero-dual near-active constraint를 등식으로 강제하던
  physical-face polish가 verified KKT~1e-14해를 KKT.003/33.9로 훼손.
  positive-dual face/열스케일링/factoredinverse 사용, 반환해는 기존 SAME
  primal/KKT tolerance를 통과해야 한다. 물리제약/오차허용을 완화하지 않음.
  독립 단순QP 재현 포함 관련8tests PASS4.01s.
- `cross_verified_polish`:112profiles actual completed, loss1128.23432566,
  heldoutRMS71.20138. Exact C114/62/32 유지, source-shape는 여전히 부적합.
  독립155q+3continuous checks positive above tail; validation49.26s 완료.
  decays(2.61928664,4.70406697,9.98631862), E_xI=-983.74934, PSD경계 근처.
  **full Al/actual yield/kinetics 채택 아님**.
- FixedQP `grid_verified_polish`도 실행종료.2starts는 이번에는 spectral cut
  residual(-5.9e-9/-1.6e-8)에서 거부됨. best2978.85는 미수렴 연구checkpoint.
  초기 numericalfail들을 최적점/physicalfailure로 재해석하지 않는다.

### 코드/테스트/다음

새핵심문서 `PUBLIC_ALUMINUM_CALIBRATION_EVIDENCE.md`,
`INTERFACE_CALIBRATION_DEVELOPMENT_V15.md`. 아직 finaltable/최종상태 추가 필요.
현재초기full449PASS532.65s/app31PASS63.77s 이후 많은 새코드가 생겼다.
최신full regression,smoke,diffcheck/최종commit/push는 **아직 실행해야 함**.
새수학targeted들은 실제 통과(공개reader,mixture,line-drag,activation,
Bloch/covariance,crossderivatives,QP 등); 최종수는 finalrun으로 갱신할 것.
다음은 actual exact-bulk tangent SVD, source/fit/heldout summary, 검증/문서와
물리적 한계 정리. Production PDE/UI/time/calibration defaults는 바꾸지 않는다.

## 초기 진행 기록: public_aluminum_v15 / interface_development_v15 (2026-09-10)

**이 구간은 아직 미완료 인계다.** 사용자 최신 우선순위는 모든 보정/단위/척도의
근거를 외부 원자료까지 찾아 검증하는 것. 약4시간 작업 요청(UTC2026-09-09
19:15 부근 시작)을 받았다. UI 재설계가 우선이 아니다. 전단/방향 상태는
이미 브리핑했다: tensor traction projection 및 W_int(a,s1,s2) STATIC은 있으나
production UI/PDE는 scalar axial, P(a,s), 독립전단/방향 입력 미연결이다.

시작 local/fresh-origin=`d7ef95cf7d8c814053d0caa741b2f2b1e3f4d338`, clean.
실제 작업은 OneDrive 밖 기존 scientific worktree. Main/reset/force/agent 없음.
현재 아래 작업은 아직 커밋/푸시하지 않았다. 기존 bulk v14/production은 보존.

- 실제 공개 crystalline Al MD: Fransson/Erhart Zenodo10014454,300K,Al99,
  6912atoms,4.065Angstrom,12 cubic repeats,NVT damping1ps,5fs integration,
  25fs saved sampling.1024frames(25.575ps)를 bounded HTTP range로 실제 읽고
  adjacent(111) periodic-plane mean3-vector를 추출했다. 전체3.5GB 파일은
  저장하지 않았고 full gzip checksum도 검증했다고 주장하지 않는다.
  source main MD5는 확인했고 potential SHA60c8a085...는 기존 static source와 같다.
  full whitened covariance min eigen=-.594207 at.150ps, empirical block+
  antisym floor.172295. First-half/8blocks/stride2에서도 negative mode를 확인.
  직접 overdamped a/s clock으로 채택 불가. 실제MD의fs/ps는 우리 solver t0 아님.
  `results/public_aluminum_validation_v15/crystalline_kinetics/`에 실제 값.

- 실제 Al wire 원자료: Verheyden/Deillon/Mortensen2018,
  DOI10.1016/j.dib.2018.11.047,PMC6265499. 공개83MB supplementary zip을 받아
  29specimen692rawrelaxation traces를 읽었다. Macro room-T monotonic tension+
  60s holds이며 fatigue/cell mobility가 아님. 대표Al_20_r_10(4N,D14.7um,
  axis[7,-1,-2]) 5.67–18.93MPa, Schmid.4838498257. Downloaded notebooks/code는
  실행하지 않았고 데이터만 파싱했다. 최종parse는692/0failures; 초기strict
  `wire_relaxation/`보고는 superseded, `wire_relaxation_verified/`가 최신.
  **추가 감사 필요:** 일부 원자료는 `Not considered/`하위 source제외기록이다.
 692를 published accepted count686과 같은 것으로 부르면 안 된다.
  5N_100_r_4의 rawforce/stress impliedD112um지만 table/notebook은114um;
  최대1.006MPa 차이를 숨기거나 임의로 고치지 말 것. No derivatives at duplicate
  timestamps;60000여 repeatedtimestamps와 unnamedextra2columns는 보존했다.

- v14 inspected heldout9rows는 이제 development로 재분류.30newoffgrid/offpath
  validation points를 loss에서 제외. Same quartic family actual joint v15fit
  71profiles loss1251.06887, heldoutRMS78.30, interface shape FAILED.
  Two starts numerical rejected;1xtolstop optimality.4778. Bulk114/62/32 유지.
  Ownfault/saddle.218260/.218464J/m2로 reversebarrier도 부적합. 범위 전체의
  불가능성 증명으로 과장하지 않는다. Independent validation47.82s완료.
- 이를 단순 최적화 초기점 문제와 구별하려고 v14best/v15best+fixed2starts에서
  max_nfev100 continuation을 실행 중. `continuation_starts.json`과
  `continuation_quartic/`, ignoredcache continuation_fit.log 참조.
  에너지 추가항/생산모델 변경 없이 기존family에 대한 추가 실제보정이다.

현재 실제 테스트: targeted42PASS24.26s; app31PASS63.77s(skip0);
full solver449PASS532.65s. 이 이후 code edits는 다시 검증해야 한다.
Smoke/diff/문서완료/최종Git은 아직. 원자료는 ignored `.cache/public_aluminum_v15/`;
대형rawzip/gzip을 커밋하지 말 것. CURRENT에 적힌 완료와 실행 중을 구분한다.

다음: continuation확인 → constrained/unconstrained residual 및 descriptor
표현력 진단 → 정량 실패가 확인되면 최소 analytic exponential-density extension만
별도 연구 가설로 검토. LJ/Bessel/원자당환경합 우선순위 불변. 새로운 density
형상은 gauge/미분/finite-qtail/heldout를 통과해야 하며 자동채택 금지.
MD/실험 source provenance, 제외기록/지름불일치, 실제 보정 실패를 문서화.
물리 M_a/M_s/t0, 실제항복완료, A_c는 여전히 얻었다고 주장할 근거 없음.

## 최신 작업: tail_calibration_v14 — 실제 벌크 보정 완료, 계면 형상은 미채택

최신 사용자 요청: “몇번째 미보정이야 보정좀 해라”. 진단만 반복하지 않고
실제 정적 Al 재보정을 수행했다. 시작 HEAD는
`a47869cdfbd9d425fa66f89380553f36e6f0e118`, fresh fetch 후 local/origin 동일,
작업 트리 clean이었다. 이전 이관된 OneDrive 밖 과학 worktree의
probability-pde-solver-v1에서 계속했다. Main/reset/force-push/agent 생성은 없다.

### 이번에 실제로 달라진 것

1. **FCC 무한 tail을 보정 제약 안에 넣었다.** `FCCTailEnvelope`는 FCC
   Voronoi cover로 power/exponential tail 적분을 상계한다. 각 coefficient의
   finite-q Hessian에 대해 H_R - sum|c|error I >=0를 polarization halfspace로
   풀었다. 유한 radius 하나에서0인 eigenvalue를 안정이라고 하지 않는다.
   코드는 무한 Poisson/Bessel 에너지를 보존하고 direct operator는 독립
   허용성 검사다. 같은8-coefficient family actual fit loss31.6274는 여전히
   C98.3205/68.0351/30.6616GPa라서 독립 탄성 오차를 해결하지 못했다.

2. **최소 대칭 채널을 유도하여 벌크를 실제 보정했다.** 기존 두 radial Q2의
   signed T2g strain response를 제거하는 Q_E=(Qodd-eta Qeven)/N을 구성.
   ONE amplitude D_E>=0가 Eg(C')를 독립 제약한다. Site 환경을 모두 더한 후
   norm을 취하고, 새 fitted orientation/length/yield cutoff는 없다.
   C11/C12/C44=114/62/32GPa, cohesion3.36eV/atom, force≈0인 후보를 실제로
   얻었다. 이 부분을 다시 “보정 안 됨”이라고 뭉뚱그리지 말 것.

3. **계면 형상 실패 후 별도 최소 ablation을 실제로 실행했다.**
   Eg만 exact-bulk loss52.3792. Density cubic는52.3156로 개선 미미하여 미선택.
   Per-atom K3||Q3||^4는 pristine harmonic stiffness를 바꾸지 않는 최소
   nonlinear angular test로 loss14.2312까지 개선. Rational I²/(1+alpha I)는
   loss13.4472지만 held-out 개선이 거의 없어 미선택. 이들을 모두 합친
   potential을 자동으로 채택하지 않았다. 새3 radial bounds/starts도 기록했다.
   Even range bound12→24 continuation은 loss9.85820이나 heldout18.39→34.02로
   악화. 범위에 걸린 fitted value/작은 training loss를 물리 최적값으로 위장하지 않는다.

4. **선택은 scoped bulk/static reference일 뿐 full Al 채택은 아니다.**
   Quartic range12 decays=(2.7829540002,5.6468406898,11.9999999974),
   coefficients(u,v,A,B,C,D3,D1,D2,D_E,K3)=
   (.3265081091,1.5095995104,10.7472347065,16.5382804733,.2218083326,
    4.3273088458,8.4877976488,.7980971836,.2647734651,139933.994322).
   L0=4.05/sqrt(2)Angstrom, energy=eV, per-atom density normalization 고정.
   K3는 작은 fourth-power invariant의 계수이며 그 자체가 bond energy 아님.
   기존 소재 파일 hash9d00fbf5...와 source60c8a085...는 바꾸지 않았다.

5. **독립 정적 검증 결과.** 실제 perfect/fault/saddle roots의 Morse0/0/1,
   force norm<2e-15. Own relaxed ISF .149229 vs source.150479J/m²(-.831%),
   saddle .179127 vs.172002(+4.142%). PerfectHaa24.0919 vs19.6607(+22.54%),
   faultHaa30.1750 vs14.1406(+113%), saddleHaa29.2931 vs12.7173(+130%).
   Direct110 curve도 크게 틀림. Work40h1.73238 vs1.74129J/m²는 약-.51%지만
   중간 opening shape는 실패. Candidate negative-traction lobe -1278.06MPa,
   source에도 -785.675MPa lobe가 있다. 129/257 curvature-root bracket으로
   같은 극값 확인; clipping이나 source-negative-lobe 은폐는 없다.

6. **과거 finite-q 불안정은 이번 후보에서 재현되지 않았다.** 독립155q,
   radius12/16의 min robust margins+.00620815/+.00642097eV/L0². 세 연속
   lambda_min/q² 검색이 L 부근의 positiveH14.4943, tail.000739로 수렴.
   이것은 전 Brillouin-zone/finite-strain/defect stability 증명은 아님.
   Local11-direction SVD condition1955.64, log step반감시Jac change5.12e-5.
   Rank11을 material confidence/unique parameter set으로 부르지 않는다.

7. **실제 낮은 MPa 시나리오도 수행.** Source/candidate 합32 static states:
   pure shear, mixed normal+2shear, compression/shear, unload, -150~150MPa.
   Candidate force residual<=1.64e-14eV/L0, unload registry<=1.61e-16L0.
   국소 탄성 상태로 회복한다. 실측 항복/피로, dynamic hold, kinetics가 아님.
   새 후보로 기존 screw core나 production PDE 결과를 재명명하지 않았다.

### 파일과 재개 시 먼저 볼 것

- `solver_v1/TAIL_CONTROLLED_AL_CALIBRATION.md`: 모든 유도/값/분류/한계.
- `tail_constrained_material`, `symmetry_resolved_material`, `quartic_angular_material`:
  기존 LJ/PDE와 분리된 검증된 수학/연구 energy helpers.
- `run_tail_constrained_calibration`: 실제 최적화. 기존 out 덮어쓰기 거부.
- `validate_tail_calibration`: completed fit replay+독립 root/q/sensitivity 검사.
- `report_tail_calibration`: 실제 MPa 정적 이완 및 saved fits 비교. refit 아님.
- `results/fcc111_active_interface/tail_calibration_v14/README.md`: 결과 라우팅.
  `calibration.json`은 실제 종료, checkpoint는 미완료. 초기3 numerical-failure
  디렉터리는 `attempt_status.json`으로 구별한다. 수정된 QP는 반환된 최종
  계수에 대해 KKT/primal/complementarity를 다시 검사한다.
  `quartic_exact_bulk/validation_refined`가 최신 독립 검증.
  `scoped_report`/존재하면 `final_report`에 portable scoped parameter JSON,
  단위/소재범위/적용불가 상태, 실제 CSV/검사한 curve 그림이 있다.

### 검증/Git — 완료 계산, 커밋 직전 인계

새 targeted22PASS67.94s. App31PASS108.85s(skip0), desktop smokePASS2.17s.
전체 solver437PASS792.45s(13분12초), 제외0. 모든 실제 계산이 종료됐다.
`.cache/tail_calibration_v14/tests_run1/`에 full solver/app 원본 로그가 있다.
`final_report`가 최종 완료-run inventory/portable parameter/정적 MPa 결과다.
Source/candidate curve를 실제로 시각 검사했고, 좋은 relaxed fault 숫자와
큰 fixed-a direct110 오류를 구별했다. Commit/push SHA는 이 문서 자기참조로
꾸며 쓰지 말고 실제 git log와 fresh origin을 확인한다. 이 문서는 커밋 직전
인계이며 후속 최종 응답이 push 검증 결과를 보고한다. Main에는 손대지 않았다.
최종 artifact JSON41/CSV24/SVG2 파싱, 이미지 시각 검사, working/staged
diff --check를 통과했다. Matplotlib SVG의 무의미한 trailing whitespace만
포맷 정리했고, 수정된 보고기 재실행도 완료했다. 수치 CSV는 동일했다.
완료 run의 중복 checkpoint6개는 final JSON의 profiles/optimizer가 정확히
동일함을 확인한 뒤 ignored `.cache/tail_calibration_v14/optimizer_checkpoints/`
로 복구 가능하게 이동했다. 원자료 삭제는 없다. 초기 실패 checkpoint는 보존.
커밋 전 fresh origin은 시작 a47869c와 동일했고 main local80cacb4,
origin/main c43d8e0도 변경하지 않았다.

### 미완료/다음 판단

벌크 보정은 완료, full-vector 계면 강성/direct110/개구 shape가 미채택 사유다.
이 검증 자료를 새 loss에 쓰면 development targets로 재분류하고 독립 검증점을
따로 남겨야 한다. 단순히 최고차항을 계속 늘리거나 mobility/source length를
맞춰 결과를 만들지 않는다. 새 소재를 core/finite-source 경로에 전달하려면
per-site angular 구현과 material gate를 먼저 통과해야 한다. 현재 production,
UI, 물리 M_a/M_s/t0/초/Hz, A_c는 변경하지 않았다. 실측 항복 완료도 아니다.

## 최신 작업: stable_core_v13 — 안정 branch 수렴과 소재 최적화/불안정성 분리

최신 “ㄱㄱㄱ” 요청으로 fab84ddc6ad922923b5f83af076eda92f56cc318에서 시작.
fresh fetch 후 local/origin 동일, 작업 트리 clean이었다. OneDrive 밖의 기존
과학 worktree, probability-pde-solver-v1에서 작업했다. 추가 agent, main 수정,
reset, 기존 파일 삭제, production/PDE/UI/kinetics/A_c 변경은 없다.

### 완료된 실제 연구 계산 — Al 채택/실제 항복 완료는 아님

1. **안정 코어에서만 domain continuation.** `core_continuation.py`는 소재
   hash/L0/center/논리적 row를 검증하고 기존 내부 변위를 새 문제의 초기값으로
   전달한다. `--extend-domain`, `--require-stable-source`가 명시적 gate다.
   새 영역/경계의 모든 affected-site F와 실제 analytic force/Hessian을 재계산.
   기존 L-BFGS default는 그대로이고 선택적인 analytic Newton-CG만 추가했다.
   R6/r6의 두 알고리즘 에너지는1.8e-13 eV/row 차이, 내부 변위는1.48e-7 L0 차이.

2. **새 8개 실제 코어 최소화 + 기존 4개를 포함한 12개 재평가 완료.**
   Ring8의 R6/R8/R10 최소 고유값은+.405864/+.261534/+.172356 eV/L0².
   무하중 R10/r8 force2.80e-10. 모든 새 상태의 force/Morse/에너지 방향차분과
   winding1을 확인했다. 원 데이터를 다시 평가한 energy 차이는 모두0이었다.
   이는 검증한 fixed boundary의 안정성이며 infinite-domain 인증은 아니다.
   Ring8에서 R6→8→10 내부 18행의 최대 변화는.004714→.002672 L0.
   Ring6→8 영향은 R6/R8/R10에서.000198/.000269/.000303 L0다.

3. **수치 partition과 물리 core radius를 분리.** `core_matching.py`의 smooth
   window에는 직접 유도한 logarithmic 상수 c_w를 뺀다. 원 에너지/힘은 불변.
   Ring8, 각 free radius의 finite part(alpha=.75)는 R6/R8/R10에서
   .377033/.372631/.369518 eV/row다. 하지만 같은 적분 반경에서 domain을
   바꾼 차이도 남으므로 완전한 outer-core matching으로 채택하지 않았다.

4. **큰 영역에서도 실제 0→50→0 MPa 전단 검사.** R8/r8의 내부 변위 증분은
   .0080933 L0, unloading 후7.10e-10 L0. Final force4.20e-10, minH+.261534.
   안정 초기점으로 회복했고 잔류 변화는 분해되지 않았다. 이를 거시 epsilon_p,
   실측0.2% 항복, kinetics/zero-stress time hold라 하지 않는다. 50MPa 단계는
   Newton optimizer precision-loss warning이 있지만 실제 force4.87e-10과
   독립 Morse/energy curvature는 통과했다. Warning도 원본 summary에 보존한다.
   작은 R4 비교는 symmetry-related 반대 stable branch이므로 동일 branch의
   하중-domain 인증으로 오해하지 않는다. Dirichlet boundary는 실측 pinning 아님.

5. **같은 소재식의 최적화 재실행.** 새 F 항/target/scale을 추가하지 않았다.
   고정 range의 coefficient problem은 exact force/cohesion + 기존 부호 bounds
   하에서 active-face convex profile로 풀었다. 실제 log-range least squares
   3 starts,588 profiles. Loss38.9096→38.9065 및 다른 basin23.9148을 찾았다.
   세 번째는65nfev budget stop(loss49.1963); global optimum이라 하지 않는다.
   최저 loss 후보 C99.649/66.016/31.187GPa, heldout RMS19.56, finite-q min
   약-76.96 eV/L0²: **불안정성을 이용한 낮은 손실이어서 채택하지 않았다.**
   Local9-direction condition3676, A=0 active face. Rank9가 물리 식별성/CI는 아니다.

6. **수정: spectral necessary constraint + 실제 range 재최적화.**
   H(q;c)=sum c_j H_j(q)에서 negative polarization v는 vᵀH v>=0라는 선형
   halfspace를 준다. Energy/force를 clipping한 것이 아니다. 고정 두 range의
   실제 constraint profile 후, 원인 모드를 유지하고 radial 재최적화도 수행했다
   (2 starts,99 profiles). 고정-range loss103.622→재최적화31.6247.
   최종 decays(5.3707449113,5.6320174464,3.2556291352), 계수 순서
   (u,v,A,B,C,D3,D1,D2): (.4130652816,3.4917299050,4.0910969697,
   23.4600618913,0,65.6450471980,-2.1169476521,-58.3030883841).
   C98.3205/68.0336/30.6618GPa, C'=15.1434 vs source26, heldout17.3383.
   Own fault/saddle .141521/.170555 J/m²는 그럴듯해도 Al 채택 근거로 부족하다.

7. **중요: cutoff 하나에서의 constrained zero는 안정성이 아니다.**
   위 후보 minH는 radius8≈0, radius10=-.0023415. 추가 radius12/16/20 검사에서
   -0.00354150, 마지막 matrix 변화8.00e-5로 음수가 약44배 크다. 고정-range
   constrained 후보도 -0.00281457 vs6.36e-5. 둘 다 음의 모드가 수렴하므로 거부.
   `material_comparison.csv`는 radius10 상태와 final refined 상태를 구별한다.
   실제 다음 해결점은 tail-controlled stability margin/adaptive q 검사를
   radial fit 내부에 포함하는 것이다. 검사 한 점/유한 radius의 조건으로
   whole-zone 안정성을 인증하지 않는다. 관측한 몇 basin만으로 전체 analytic
   family가 불가능하다고 단정하지 않는다. 임의 positive margin도 새 경험항도 금지.

### 파일 / 재개 경로

- `solver_v1/STABLE_CORE_AND_MATERIAL_DIAGNOSTICS.md`: 수식, 수치, 단위, 제한.
- `results/fcc111_active_interface/stable_core_v13/`: 실제 자료, old v12는 보존.
  `isolated_core/*/{metadata,summary}.json,state.csv`가 actual minimized state.
  `core_*` CSV와 SVG는 동일 state 재평가/수렴 비교. `core_decision.json` 참조.
  `material_reprofile`, `material_spectral`, `material_stable_ranges`의 completed
  JSON은 실제 최적화 기록. `material_*` 비교/잔차/고파수/반경12–20 검사가 있다.
- `run_profiled_range_calibration`, `run_spectral_material_profiles`,
  `run_spectral_range_continuation`은 실제 fit을 수행한다. `report_stable_*`는
  저장된 결과를 검증/재평가할 뿐 fit 재실행이 아니다. 실행 기본 경로는 기존
  결과 덮어쓰기를 거부한다. 중복 중간 optimizer snapshot은 ignored .cache의
  stable_core_v13 아래로 보관했고, complete profiles는 final JSON에 모두 남았다.
- Historical 소재 hash9d00fbf5...과 source Al99 hash60c8a085...는 불변이며
  문서/데이터에 full SHA256을 기록한다. 새 fitted 후보로 코어를 바꾸지 않았다.

최종 실제 검증: 새/기존 코어 및 profile targeted34PASS42.53s,
전체 solver415PASS1810.82s(30분10초), app31PASS200.25s(skip0),
desktop smokePASS2.51s. 실제 연구 계산/전체 테스트가 모두 종료됐다.
git diff --check와 staged diff check도 통과했다. 새 자료 JSON36개/CSV45개/
SVG1개를 파싱/검사했고, curve를 실제로 확인했다. 모든 기존 테스트는 유지했다.
이번 source/code/data 변경은 과학 검증을 마친 연구 결과로 커밋하는 것이며,
소재 채택 성공으로 커밋하는 것이 아니다. 이 문서는 해당 커밋에 포함되는
커밋 직전 인계다. 자기 자신의 SHA를 이 문서에 쓰지 말고 실제 git log 및
fresh origin/probability-pde-solver-v1로 최종 commit/push를 확인한다.
커밋 전 다시 fetch한 origin은 fab84ddc6ad922923b5f83af076eda92f56cc318로
시작과 같았다. main local80cacb4/remote c43d8e0는 건드리지 않았다.
일부 초기 계산의 긴 wall-clock 공백 원인은 미확인이다. CPU benchmark로
해석하지 말고 실제 iteration/evaluation 수를 함께 본다.

### 다음에 할 일 / 아직 넘지 못한 gate

- 먼저 아래 최종 검증/Git 기록 또는 실제 Git log를 확인하고 중복 재실행하지 않는다.
- 소재: finite-q tail/전 q 안정성을 포함한 견고한 constrained optimization,
  독립 탄성/계면 곡률 검증. 기존 조건 실패만으로 새 F 항을 무작정 늘리지 않는다.
- 코어: 더 큰 domain 또는 일관된 outer response, 모든 line character/partial
  core와 finite-source matching. 현재 straight screw만으로 finite-loop/source
  energy를 완성하지 않는다. 새 core-radius나 가상 pin length로 실측 항복을 맞추지 않는다.
- 실제 source/결함 geometry,0.2% strain까지의 거시 연결, physical mobility는 여전히
  필요하다. Static source/core와 probability/PDE kinetics를 혼동하지 않는다.
- 실측 Al 항복/피로, physical seconds/Hz, A_c 보정은 미완료. Production/UI gate는
  닫혀 있으며 사용자 확인 없는 meshing UI 재설계는 하지 않았다.

## 최신 작업: range_core_v12 — 소재 range 분리와 단일 전위 코어 경계

요청: “이제 좀고쳐라 걔네도”. b84eaf4d9a449a16476ab69d8d3de30496c6c00c에서
fresh fetch, local/origin 동일, clean으로 시작했다. 같은 OneDrive 밖 과학
worktree/branch에서 작업하며 추가 agent, main 변경, reset은 없다.

### 완료된 계산과 물리 판정

이번 단계의 실제 보정, 21개 코어 계산, 원 상태 재평가, 전체 회귀 검사를
완료했다. **실제 Al 항복강도/유한 전위원/피로 솔버 검증 완료가 아니다.**
새 소재 후보는 채택하지 않았고 production/PDE/UI는 바꾸지 않았다.

**핵심 원인과 수정:** centered Volterra seed는 force가 작아도 대칭 안장점에
멈춘다. Historical R4/r4, R6/r4, R6/r8, R8/r6의 최소 고유값은 각각
-0.83311, -0.74230, -0.75406, -0.71322 eV/L0^2였다. 새 소재 R4/r4도
-1.02819였다. `core_stability.py`에서 analytic Hessian-vector 최소 모드와
독립 실제 에너지 방향차분을 검사한다. 새 core 실행은 이 Morse 검사가 필수다.
음의 모드 양쪽을 명시적 deterministic seed로 검사하되 전위/경계/힘은 불변이다.

Historical R6/r4의 양쪽 이완 결과는 energy 1.032538443572 eV/row,
minH +0.3977618, force≤2.38e-7 eV/L0로 낮은 안정점에 도달했다.
원 안장점보다 0.00922177 eV 낮으며 양쪽 에너지는 6.4e-14 eV 이내로 같다.
R6/r6 같은 안정 branch는 minH +0.4045673, force 1.94e-7이다.
단, R4→R6 stable-domain 변화는 내부 동일 18개 row에서 0.013723 L0,
R6의 ring4→6 변화는 0.0012741 L0이다. 완전한 core/domain 수렴이 아니다.

**가짜 잔류 원인 분리:** 원 R4/r4 안장점에서 25→50→0 MPa 후 변위는
초기점과 0.125795 L0 달랐다. 그러나 무하중 negative-mode control과는
5.32e-8 L0 이내로 같았다. 이것은 stress-induced residual plasticity의
증거가 아니라 잘못된 초기 안장점 이완이다. 안정 무하중점에서 다시
25→50→0 MPa를 실행하니 내부 최대 변화는 각각 0.0037046, 0.0073759,
1.0143e-7 L0였다. 세 단계 모두 minH>0.64, force≤2.05e-7이다.
이 static branch에서는 분해 가능한 잔류 변형이 확인되지 않았다.
이는 core 변위이지 거시 epsilon_p가 아니며, 물리 시간 hold도 아니다.

- `range_resolved_material.py`: STF ranks1/3와 rank2의 analytic exponential
  range를 분리했다. 추가 microscopic shape 변수는 한 개이며 LJ/Bessel,
  scalar density/per-atom embedding과 기존 파라미터는 불변. Common-range
  limit의 전체 energy/gradient/Hessian이 기존 코드와 일치한다.
- 실제 고정39 starts와 두 metric의 Powell 프로파일292개 실행(438.55s).
  C별5% loss46.7495→38.9096, heldout normalizedRMS12.2948→9.4344.
  새 C는87.5166/69.9978/32.8059GPa: target114/62/32를 아직 만족하지 못한다.
  relaxed fault=.151279, index-one stationary energy=.168804J/m2;
  source=.150479/.172002. Opening40h=1.827637 vs source1.741285J/m2.
  Source saddle Hxx heldout 부호/크기도 부정확하다. **채택하지 않았다.**
- Full local sensitivity rank9, condition196.70(명시한 coefficient/log-range
  좌표), radial step refinement 차이1.85e-4. 이를 confidence나 global unique
  calibration이라 하지 않는다. C별 Powell180평가 예산소진; bulk_exact
  Powell34평가 종료점은 incumbent보다 나빠 전체 best를 보존했다.
- 중요한 새 원인: 기존 exact-C 후보는 positive bulk/perfect-interface H에도
  finite-q Gamma-K halfway에서 eigenvalue=-25.829eV/L0^2. 반경6/8/10으로
  재확인한 불안정이며 exactC가 채택 근거가 될 수 없다. 새 C별 후보는
  검사한 finite-q에서 양수지만 global stability/material validation은 아니다.
- `isolated_screw_core.py`: 같은 원자 전위의 anisotropic single-screw
  far-field Dirichlet 경계, 내부3성분 이완. 영향 받는 **고정 원자의 F까지**
  에너지 closure에 포함한다. x는 무한 Poisson/Bessel, transverse environment
  변화의 ring은 독립 수렴 대상. 주기 반대전위쌍 소멸과 다른 실험이다.
- 새 row kernel은 rank별 range를 정확히 사용한다. 기존 common-range scalar
  FFT에 새 range surface를 넣으면 명시적으로 거부한다. 기존 같은-range
  연구/production 계산은 바꾸지 않는다.
- R2/3/4/6/8, ring3/4/6/8, 무하중 ±mode control 및 안정/불안정 초기점의
  25→50→0 MPa를 포함한 21개 실제 코어 계산을 완료했다. force는 대략
  1e-7–1e-6 eV/L0이고 winding 1을 유지한다. 9개 case는 검사한 고정 경계에서
  안정 후보다. 미검사된 예전 소규모 case는 stability_unchecked로 남겼다.
- 코어 raw partial-site energy와 탄성 annulus 불일치의 선형 경계 에너지항을
  유도했다. `e_remainder=e_raw-sum g_ref·Du`의 전체 합은 원래 에너지와 같고
  free gradient/Hessian은 불변이다. Partial radial sum만 다른 partition이다.
  Raw/linear/remainder를 모두 보존한다. 임의 energy 보정/force clipping이 아니다.
  무한 row affine Hessian은 같은 full-bulk tensor와 독립적으로 일치한다.

### 재개 위치 / 파일

`solver_v1/RANGE_AND_CORE_REPAIR.md` 수식, 새 runner/report/test를 읽는다.
결과 `results/fcc111_active_interface/range_core_v12/`:
material/calibration.json은 실제 최적화 결과(모든 평가 포함),
material/validated_summary.csv,finite_q_validation.csv,static_load_unload.csv,
isolated_core/각case/{metadata,summary,progress}.json,state.csv.
`report_range_core_repair.py`는 저장 state의 analytic audit이며 최적화 재실행이
아니다. 모든 계산 종료 후 실제 실행하여 core_* CSV, decision.json과 두 SVG를
생성했다(최종 재평가 68.31s). core_morse_audit.csv에는 독립 에너지 차분을
저장한다. 원 안장점 기록/raw partition도 삭제하지 않았다.
한 코어 force 수렴이 all-domain/core/finite-source/실제항복 수렴을 뜻하지 않는다.
Straight x-line core만으로 curved source 전체 character energy를 만들지 않는다.
실측 source geometry, 실제0.2%plastic strain, kinetics, material gate 미완료.
물리 seconds/Hz/A_c 및 production/UI default는 변경하지 않는다.

### 실제 검증 결과 / 재현

- 새 stability/core/range targeted: **23 passed, 16.52s**.
- 새+기존 vector row targeted: **36 passed, 52.08s** (최종 두 Morse test 추가 전).
- 최종 `pytest solver_v1 -q`: **394 passed, 631.95s**.
- `pytest app -q`: **31 passed, 148.56s**, skip 0.
- `app.desktop_ui --smoke`: **PASS, 4.66s**. LJ a0/kappa 역사값 그대로.
- `git diff --check` 및 새 파일을 포함한 `git diff --cached --check`: **PASS**.
- py launcher 대신 설치된 Python 3.13을 사용했다. 최종 solver 실행의
  OPENBLAS_NUM_THREADS=1, OMP_NUM_THREADS=1은 해당 프로세스 성능 설정뿐이다.
  로그는 ignored .cache/range_core_v12_*tests.log 및 smoke.log에 있다.
- 기존 파라미터 SHA256는
  9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655로 불변.
- 이 인계 문서 자체가 v12 변경에 포함된다. 최종 commit/push 여부와 SHA는
  Git 실제 log/remote로 확인한다. 테스트 통과를 소재 채택으로 해석하지 않는다.

### 다음 단계의 정확한 미완료 문제

1. 소재: 독립 탄성 세 값, vector saddle curvature, finite-q 안정성을 동시에
   만족하는 작은 analytic family/후보가 아직 없다. 이번 range 한 개 추가만으로
   해결되지 않았으며, 전체 family 불가능을 증명한 것도 아니다.
2. Core: 낮은 안정 branch를 R8 이상/ring8 등으로 이어서 domain/tail/Morse
   수렴을 확인한다. 안장점의 그럴듯한 annular plateau를 채택하지 않는다.
3. 실제 항복: edge/mixed character, 유한 전위원/실측 source geometry,
   안정 초기상태에서의 방출과 실제 plastic-strain 기준 연결이 아직 필요하다.
4. 실제 Al kinetics/초/Hz, A_c 보정은 여전히 없다. UI 재설계 gate 미통과.

## 이전 완료 상태: yield_bridge_v11 — 실제 항복 기준, 독립 탄성, 유한 전위원 reference

최신 요청: “실제 항복강도에 점점 맞춰가야지”. 실제 과학 worktree는
OneDrive 밖 al-fatigue-probability-worktrees/aft-pde-bessel-38969ad,
branch probability-pde-solver-v1이다. 시작 fresh fetch에서 local/origin 모두
**062c83d14f33f8ac0c1c4d6bdff25804394ff2bd**, clean. 추가 agent 없음.
다른 worktree, migration 백업, main은 변경하지 않았다.

### 이번에 무엇이 진전되었고, 무엇이 아직 안 되었는가

1. 실제 실험의 **0.002 plastic-shear CRSS/G** 점을 확보했다. 과거 .1/.5/.8
   large-strain flow와 달리 항복 판정 기준이 명시되어 있다. 다만 공개된
   정규화 계수 G의 수치를 확인하지 못해 MPa는 null이다.
2. 같은 LJ/Bessel 에너지의 독립 C11/C12/C44 metric을 새로 검사하고 실제
   최적화를 재실행했다. C44만 개선하거나 bulk 전체를 정확히 맞춰도 계면과
   함께 Al을 만족하지 못한다. 새 파라미터는 **채택하지 않았다**.
3. 같은 전위의 탄성 tensor에서 **유한 양단 고정 전위원의 leading-log
   외부 탄성 reference**를 유도하고 실제 shape 평형을 풀었다. 가정한
   micrometre 크기 전위원이면 MPa bow-out 척도가 나온다. 이것은 새 경험적
   yield law, atomistic core 검증, 실제 시편 항복 또는 fatigue 예측이 아니다.

자세한 수식은 solver_v1/YIELD_STRENGTH_BRIDGE.md,
결과는 results/fcc111_active_interface/yield_bridge_v11/SCHEMA.md를 읽는다.
Production energy/PDE/UI, 기존 파라미터 파일, kinetic calibration, A_c 변경 없음.

### 실제 항복 기준과 출처

Krebs2017 doi10.1038/nmat4911 Fig2d 원본 3508x2480 그림의 서로 겹치지 않는
빨간 marker 7개를 centroid/pixel box/축 좌표와 저장했다. 조건은99.99% as-cast
Al single-crystal wire, 실온, tensile displacement300nm/s, gamma_p=.002.
CRSS/G는1.5085e-4–5.2216e-4, 판독 halfwidth8.57e-6. 이는 전체 표본 범위나
실험 scatter가 아니다. MPa/G 숫자, axial Schmid factor, source geometry를
추정해 넣지 않았다. 과거 Fig2b flow 점의 의미는 그대로 보존한다.

공식 EPFL ORIGINAL bundle의 SI를 실제 획득했다. Legacy URL 실패를 우회한
정상 공개 API이며 fetch_strength_reference.py의 URL만 갱신했다.
SI SHA256=16cec3051cd27303f3a77c20779c5aaf19c2714d1c2664c72e6667b663ebdb74,
공식 MD5=7787b4507b8f25b0da6f17a829b02abb와 일치한다. SI23–25쪽의L=D/3와
대안D/2는 single-arm 모델 가정이지 실측 pin 길이가 아니다. 논문의 DD 이동도를
우리 집단좌표 M으로 옮기지 않는다. 2019 annealed-wire 논문도 확보했으나
열처리/산화막/하중 조건이 달라 별도 provenance로 남겼고 fit에 사용하지 않았다.
PDF/원본 이미지는 ignored cache, 7개 읽은 점과 출처만 커밋한다.

### 탄성 metric 감사와 실제 재보정

H=Q C의 정확한 rank3 변환을 사용한다. 기존 diagonal H metric을 순수하게
좌표 변경하려면 Sigma_C=Q^-1 Sigma_H Q^-T의 상관을 유지해야 한다.
이번 diagonal 5%-C metric은 **다른 discrepancy 가정**이지 단위 버그 수정이
아니다. 새로운 target/weight를 실제 실험 불확실성이라고 주장하지 않는다.

26개 고정 radial starts와 Powell을 실행했다. cubic_5pct는70평가 후 local
수렴, bulk_exact는100평가 예산소진. 후자는 optimizer 반환223.0933보다
마지막 평가점223.0087이 조금 더 좋아 전체 evaluated best를 저장했다.
5bulk exact일 때 남는 3개 null 방향은 active-face 선형대수 QP로 해결한다.
A=0 같은 활성 경계는 제약식이지 energy/force clipping이 아니다.
최초 SLSQP run의 경계 roundoff 문제 때문에 active-face로 재실행했으며,
초기 material_metric은 superseded, 최종은 **material_metric_refined**다.

| 항목 | 기존 후보 | C별5% metric | 5bulk exact | 0K source |
|---|---:|---:|---:|---:|
| C11/C12/C44 GPa |87.82/64.83/49.27|85.73/70.68/32.46|114/62/32|114/62/32|
| relaxed fault J/m2 |.126800|.145773|.027060|.150479|
| index-one stationary energy J/m2 |.173238|.178987|.042449|.172002|
| opening at40h J/m2 |1.762768|1.798845|1.741053|1.741285|
| heldout normalizedRMS |15.6863|12.2948|25.4113|—|

모두 equilibrium/cohesion을 맞추고 perfect-interface H는 positive지만,
C별 fit은 C11 약24.8% 오류, bulk_exact는 fault/안장 energy 약82.0%/75.3%
과소평가다. Opening endpoint만 잘 맞는 것으로 채택하지 않는다. Finite-q
전체 안정성/새 global MEP는 이번에 인증하지 않았다. 후보 두 개의 실패가
전체 analytic family 불가능의 증명도 아니다. 두 radial과 활성경계를 제외한
conditional coefficient SVD를 물리 parameter confidence로 부르지 않는다.
기존 SHA9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655 불변.

### 유한 source의 수학과 실제 응력 실행

기존 Schur K(q)=|q|K0와 step disregistry b/(iq)에서 +/-q를 함께 적분:

    k_line(theta)=b^T Re[K0(n_perp)] b/(2pi) [J/m]
    gamma(theta)=k_line(theta) ln(R/r_core)
    E[y]=integral gamma(theta) dl - tau*b*integral y dx
    T=gamma+gamma'', Q=gamma sin(theta)+gamma' cos(theta)
    p=tau*b=2Q(theta_e)/L
    tau_c,outer=2gamma(pi/2)/(b L)

양의 line stiffness와 mirror symmetry를 확인한 local-line leading-log 문제다.
R/r_core는 변분 동안 고정한다. Finite core와 nonlocal finite part는 미포함이며,
물리적으로0이라고 보정한 것이 아니다. 같은 원자 전위 tensor를 유지하며
source 탄성값은 별도 comparator로만 둔다. 경험적 line prefactor를 fit하지 않는다.

기존 후보, R=L, r_core=b라는 **명시적 가정**의 outer-only 값:

| 가정한 pin 간격L um | critical resolved shear MPa |
|---:|---:|
|.5|29.0441|
|1|15.8705|
|2|8.60943|
|5|3.80027|
|10|2.03498|

L=1um에서 r_core/b=.5/1/2이면17.2189/15.8705/14.5221MPa. 이것은 실측
source나 Al 항복값이 아니다. 5모델 x5길이 x2/4/10/25/50MPa의125사례 중
53개 subcritical graph를 실제 풀고,72개는 이 branch의 fold 위라 제외했다.
Fold 위라고 atomistic source operation/multiplication을 확인했다고 하지 않는다.
32/64/128/256segment 독립 Newton과 analytic parametric shape를 비교했다.
기존 후보 theta=.8의 max-bow/L 오류2.033e-4→5.134e-5→1.287e-5→3.220e-6.
theta1.4 근처256segment는 model에 따라6.33e-4–9.92e-4L 오차가 남는다.
80개 refinement force residual 최대1.99e-10(normalized). 각도64→128 k오류
≤1.16e-23J/m, k+k''오류≤4.82e-20J/m. 실제 작은 수치오차와 아직 큰
물리 core/source 불확실성을 구별한다. Subcritical bow는 정적 제하 시 가역이다.

b*swept_area/V_specimen은 slip의 부피평균 운동학일 뿐, source 첫 작동이
0.2%plastic shear에 도달한다는 뜻이 아니다. 실제 pin/표면/전위밀도와
loop traffic/상호작용이 추가로 필요하다. V는 실제 specimen geometry일 때만
사용하며 임의 activation volume/A_c로 대체하지 않는다.

### 실제 검증과 다음 단계

최초 targeted19개 중 새 spectral second-derivative test1개가 Schur/FFT의
N^2 증폭 roundoff 때문에 실패했다. 에너지/힘을 바꾸지 않고 epsilon*N^2*k
오차 유도를 시험에 반영했다. 기존 테스트를 완화한 것은 아니다.
추가 benchmark를 포함한 최종 targeted22PASS3.99s, fullsolver371PASS767.52s,
app31PASS91.74s, smokePASS2.14s. 전부0skip. Solver/app은 연구 실행과 일부
겹쳐 고립 performance benchmark가 아니다. 실제 JUnit output 확인 완료.
초기 material run517.66s, 최종402.17s, source28.71s, stationary/FD/report18.43s.
FD step4e-5→2e-5에서 오차약4배감소; fine max gradient2.07e-8eV/L0,
Hessian5.91e-7eV/L0². 새 두SVG는 실제 결과로 생성/시각 확인했다.
JSON/CSV parse 및 private-path/nonfinite 검사, working/staged git diff check
통과. Commit 직전 fresh origin도 시작062c83d와 같았고 main은 불변이다.

현재 결론: **MPa source 메커니즘 reference 진전, 실제 항복 재현은 미완료**.
다음은 joint material compatibility, 독립적으로 검증한 vector/partial core,
유한 선방향 source 및 실제 미세조직/시편 조건이다. 원자 이상강도를 낮추거나
실험에 맞는 L을 역산해 끝내지 않는다. 물리 M_a/M_s/t0와 초/Hz는 여전히
unavailable. UI redesign gate와 production 승격은 열리지 않았다.
새 source/재료 모듈과 실험 데이터는 별도 연구이며 완료된 calibratedAl로
보고하지 않는다. Git 최종 해시는 포함 커밋/실제 remote 확인을 따른다.

## 이전 상태: material_strength_v10 — 실제 강도와 소재 보정의 검증 조건 분리

최신 요청은 “실제강도에 가까워야” 및 “해봐”였다. 실제 과학 worktree는
OneDrive 밖 al-fatigue-probability-worktrees/aft-pde-bessel-38969ad,
branch probability-pde-solver-v1이다. 시작 fresh fetch에서 local/origin 모두
0f0a3b43b6dd3b39d3f898671702d309d8e85c94, clean. 추가 agent 없음.
이전 migration과 다른 worktree/main은 변경하지 않았다.

### 실제 실행한 것과 채택하지 않은 것

같은 LJ/Bessel + 기존 per-atom scalar embedding/STF1/2/3 가족을 사용했다.
vector_material_calibration.py에서 고정 두 radial decay에 대해 전체 energy,
force, Hessian을 8개 coefficient의 **정확한 선형 basis**로 계산한다.
Finite table surrogate가 아니다. Bulk force/cohesion 두 제약을 정확히 제거,
14 fit/9 held-out +2 exact 관측량으로 실제 constrained fit을 수행했다.
Saddle/fault는 source full-vector stationary state, reverse barrier는 종속량으로
중복 fit하지 않는다. 0K source와 실온 실험 강도는 분리했다.

25-point grid+old decay와 Powell160회를 실제 실행했다. Powell은 평가 예산을
소진하여 global/local radial optimum 수렴을 주장하지 않는다. Inner coefficient
QP는 수렴했다. Runner는 fit/validate를 분리한다. staged coefficient ablation은
최종 decay에서 실행했으며 각 stage의 independent radial fit이 아니다.

첫 joint fit은 sample 사이 a/h=2.79634에서 -46.4765MPa 개구력을 냈다.
이를 숨기지 않고 W_aa=0 실제 극값을 찾아 세 차례 coefficient constraint
exchange를 수행했다. Raw 힘을 clipping하지 않았다. 최종 zero-registry
음의 값은 -5.72e-15eV/L0, 측정/roundoff floor3.77e-14보다 작다.
121/241 bracket와 tolerance2e-11/2e-13로 확인했다. 전체 registry/무한 구간의
전역 단조성 증명으로 과장하지 않는다.

**추가 반례:** 최종 그림 검사에서 SOURCE 자체의 a/h2.386631 개구력이
-785.6754MPa였다. 181/361 bracket, energy FD, 별도 scalar-sourceenergy가
일치한다. 따라서 nonnegative-traction/단조 energy 제약은 이번 fit의 명시적
shape prior이지 보편 물리 법칙이나 source의 검증 성질이 아니다.
그 제약을 빼는 추가 QP도 실제 실행: loss28.4043, C11/C12/C44
85.16/64.07/52.70GPa로 탄성 문제는 그대로다. Source cutoff/interpolation
때문인지 본래 potential인지 더 검증 없이 단정하지 않는다. 음의 힘만으로
bug/불가능한 물리라고 주장하지 않는다. Raw source curve를 보존했다.

최종 연구 후보 계수와 모든 결과는 material_strength_v10/opening_exchange에 있다.
과거 parameter SHA256
9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655는 **그대로**다.
Production energy/PDE/UI 변경 없음.

| 항목 | source | 이전 후보 | 새 연구 후보 |
|---|---:|---:|---:|
| C11/C12/C44 GPa | 114/62/32 | 87.82/64.83/49.27 | 85.26/64.24/52.54 |
| forward saddle J/m2 | .172002 | .173238 | .166823 |
| intrinsic fault J/m2 | .150479 | .126800 | .151991 |
| reverse barrier J/m2 | .021523 | .046438 | .014831 |

Loss87.4386->28.6052지만 heldout normalizedRMS9.7749이며 source-state saddle Hxx도
부호가 다르다. 자기 saddle의 Morseindex1과는 구별한다. **채택하지 않는다.**
전체 analytic family 불가능을 증명한 것은 아니다. Initial joint Jacobian8개
singular values최소1.592/최대328.949,condition206.60; A/C sensitivity cosine.9913.
활성 부등식/model discrepancy 때문에 이를 parameter confidence로 해석하지 않는다.

### 실제 강도에 관한 진전/한계

Krebs2017 doi10.1038/nmat4911 원문 Fig2b의 99.99%Al,103um 단결정 와이어,
실온300nm/s 조건에서 plastic shear strain .1/.5/.8의 resolved shear flow
4.595/6.929/9.881MPa를 픽셀/출처/hash와 저장했다. 판독오차약.1905MPa.
이 값은 항복강도가 아니다. .2점은 cyan curve에 가려 null,0.2%CRSS는
이미지 strain해상도 부족으로 null이다. Orientation/source geometry도 없다.
와이어 직경을 source length로 바꾸지 않는다. 실험 강도를 energy loss에 넣지 않았다.

같은 MPa 응력을 old/new/source 세 rigid-interface 모델에 실제 적용했다.
최종 후보 ux/L0는 약.000121/.000183/.000261이고 static unload 차이8.72e-16L0.
33개 저응력 상태와 32개25–150MPa mixed/static 상태를 저장했다.
단일 homogeneous 계면은 finite source가 아니므로 실제 와이어 소성 재현이 아니다.
새 ideal fold2.58174GPa를 실제 항복 개선이라고 말하지 않는다.
StrengthConditions 비교기는 observable/응력성분/straincriterion/microstructure/
온도/protocol이 달라지거나 미지이면 prediction error를 null로 둔다.
유한 source/core, 소재 곡면, kinetics가 여전히 필요하다. 초/Hz는 unavailable.

### 실패 기록과 수치 검증

첫 root폴더 fit-report는 bulk-only zeroD1/D3 sensitivity column 때문에 실패했다.
Zero column을 stage에서 제외하여 고쳤고 미완료 파일은 SCHEMA에 명시해 보존했다.
최종 계산 경로는 joint_fit와opening_exchange다. 첫 exchange launch는 SOURCE의
완전히 빈 이웃 배열에서 einsum gradient가 NaN인 사례를 포착했다(NumPy2.5).
빈 합의 정확한0을 명시 반환하도록 **source-only evaluator**를 수정했다.
라이브러리 내부 allocation 버그의 독립 재현/근본 해결을 주장하지 않는다.
Nonempty nonfinite guard유지. 빈 einsum을 금지하는 test와30회 separated반복검증.
수정 후 source fit target 재계산차이0이다.

독립 direct24/48/72 maxenergyerror4.01e-5->5.24e-6->1.58e-6eV/cell,
reciprocal refinementHerror6.39e-14eV/L0². FD fine gradient7.24e-9/H1.18e-7.
Finite-q8/12컷오프 sampledpositive이나 전BZ 증명 아님. Independent fold31/61.
최종 targeted27PASS18.92s, fullsolver349PASS1429.02s(23분49초),
app31PASS161.37s, 최종smokePASS2.73s. 모두0skip. Full/app은 다른 연구 검증과
동시 실행되어 timing은 고립된 성능 benchmark가 아니다. git diff --check 및
staged diff check PASS. CSV/JSON syntax, source target변화0, historical parameter
hash보존 확인. 실제 verification.json과 포함 커밋의 Git이 최종 근거다.

자세한 수식/결과: MATERIAL_TO_SPECIMEN_STRENGTH.md,
results/fcc111_active_interface/material_strength_v10/SCHEMA.md 및verification.json.
다음은 탄성+held-out vector shape의 공통 material compatibility/정규화 민감도,
그 뒤 finite-source/core다. 소재 gate를 건너뛰어 무리하게 UI로 넘어가지 않는다.

## 이전 상태: vector_registry_v9 — 지정 직선의 장벽과 실제 벡터 경로를 분리

최신 사용자 요청은 “이제 다시 연구해.”였다. 연구 전체를 OneDrive 밖으로
옮긴 뒤, 실제 probability-pde-solver-v1 과학 worktree에서 재개했다.
Fresh fetch에서 local/origin 모두 f603fb81674333356d84a0bf959c4a9c95100ae2였고
작업 트리는 깨끗했다. 이전 migration 백업/다른 worktree/사용자 파일은 건드리지
않았다. 아래 OneDrive/Temp 경로 언급은 이전 단계의 역사적 기록이다.
추가 agent 없음. 현재 연구 위치는 로컬 GitHub 트리 안의
al-fatigue-probability-worktrees/aft-pde-bessel-38969ad이다.

### 이번에 실제로 한 단계

v8의 작은 셀 fault를 해석하기 전에 같은 후보의 full registry 계면을 검사했다.
vector_interface_reference.py: q=(a,ux,uy), 전체 analytic gradient/Hessian.
기존 plane LJ/Bessel, exact Hurwitz-zeta G=0, density 합 뒤 per-atom F,
STF1/2/3 norm을 그대로 사용한다. 새 경험식/파라미터 fit/kinetics 없음.
MishinVectorInterfaceReference는 checksum 검증 NIST Al99 source의 동일 rigid
half-crystal 비교 전용이다. Production energy/PDE/UI에는 연결하지 않았다.

vector_registry_audit.py: full Hessian index를 확인한 minimum/saddle,
saddle 양쪽 downhill 연결, 고정 x에서 a,y 이완 후 Schur curvature를 구한다.
이 local constrained branch는 global MEP 증명이 아니다.
run_vector_registry_audit.py와 run_vector_registry_rechecks.py는 실제 계산을
다시 실행한다. 결과 재생을 최적화/실행이라고 보고하지 않는다.

### 실제 결과 / 해석

기존 direct110 중간(h,.5,0)은 W_x≈0여도 candidate W_a=-3.4654,
W_y=-1.82985eV/L0가 남아 full stationary saddle이 아니다. 정상/횡방향
구속 반력이 필요한 경로였다. 새 결과:

| [J/m2] | 후보 | 동일 조건 Al source |
|---|---:|---:|
| a 이완 후 연결된 forward saddle | .173237970 | .172002365 |
| intrinsic fault | .126800169 | .150479477 |
| fault에서 돌아가는 barrier | .046437801 | .021522888 |

Forward saddle은 약0.72% 차이지만 reverse barrier는 약2.16배다. 이것만으로
Al material gate를 통과시키지 않는다. C11/C12/C44≈87.82/64.83/49.27GPa의
기존 탄성 오차도 그대로다. 모든 fit 파일 및 candidate SHA
9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655를 보존했다.

구속 조건만 바꾼 첫 shear-traction 최대값:

| [GPa] | 후보 | source |
|---|---:|---:|
| a,y 고정 | 7.54336 | 6.46164 |
| a 자유, y 고정 | 4.76890 | 4.39037 |
| a,y 자유 | 3.12402 | 2.76750 |

마지막 행만 H_zz>0 및 full Hessian zero-mode를 확인한 uniform-interface fold.
41/81 및31/61 독립 bracket을 검사했다. Source도 GPa가 나온다는 결과는
단위가 틀렸다는 뜻이 아니라 flawless rigid-half ideal loading의 결과다.
이를 실험 항복으로 낮추는 임의 factor/mobility/barrier fitting은 금지한다.
Finite source, 실제 core/material, kinetic mapping이 다른 문제로 남아 있다.

25–150MPa pure shear / 두 전단성분 mixed / 압축+shear의32개 정적 상태 실행.
모든 경우 intact positive-Hessian branch; 제하 후 registry 차이≤1.76e-15L0.
소성/피로가 생겼다고 하지 않는다. 정적 제하는 physical-time hold가 아니다.

### 수치 검증과 저장

vector_registry_v9/SCHEMA.md와 VECTOR_REGISTRY_AND_STRENGTH_AUDIT.md를 읽는다.
주 계산122.69s, 제약/stencil 재검사29.46s. 882 grid points,22 special states,
32 load/unload states. 직접합24/48/72의 일반 off-path energy 오차는
8.85e-6 -> 1.16e-6 -> 3.48e-7eV/cell: algebraic tail을 machine epsilon으로
위장하지 않았다. reciprocal/depth tol2e-11->2e-13의 Hessian 차이≤4.50e-12.
Source finite-difference의 spline-knot crossing 오차를 숨기지 않고 별도
stencil refinement에서 2.48e-10eV/L0²까지 감소한 결과를 저장했다.

새 targeted17 tests PASS(23.63s), 재실행17 PASS(23.92s), 기존 vector와 합쳐
34 PASS(31.04s). App31 PASS(58.43s), smoke PASS, 기존 LJ a0/kappa 불변.
첫 full suite는321 passed / 새 source FD test1 failed(497.22s): FD 배열 NaN.
단독 및9600 derivative /3000 allocation 반복은 재현하지 못했다.
원인을 지레 단정하지 않고 source plane/최종 jet의 nonfinite fail-fast 진단을
추가했다. 전체 재실행의 최종 결과와 이 anomaly의 상태는 verification.json을
확인한다. 허용오차/skip으로 실패를 없애거나 원인이 밝혀졌다고 꾸미지 않는다.

최종 full 재실행322 PASS,0 skipped(508.27s), 이후 targeted17 PASS(19.82s).
최종 smoke1.26s PASS, git diff --check PASS. 새 fail-fast 코드로 특별 상태22,
하중32, fold4, energy grid882를 재계산했을 때 저장값과 차이0이었다.
postguard_data_revalidation.json에 현재 kernel hash와 검사 내용을 기록한다.
이 반복 성공을 단발 NaN의 원인 규명/완전 제거 증명으로 과장하지 않는다.
기존 app/data/PDE/LJ 및 v7/v8 코드와 후보 parameter는 시작 HEAD 대비 변경0.

### 다음 단계 / gate

먼저 vector_registry_v9의 최종 verification.json과 git log를 확인한다.
독립 탄성, intrinsic-fault 및 reverse barrier/위치/곡률을 함께 제약하는
material compatibility/identifiability를 감사한 뒤 필요한 최소 analytic
extension만 검토한다. 다음 finite source/core 연구는 그러한 surface 검증과
별도로 요구된다. 현재 source는 검증용이지 production LJ 대체물이 아니다.
Mobility/seconds/Hz/A_c 및 UI redesign gate는 미완료 상태 그대로다.

## 이전 상태: vector_core_v8 — 수직·횡방향 구속을 실제로 해제

사용자의 최신 ㄱㄱ는 v7에서 발견한 omitted transverse force를 해결하는
다음 단계 승인으로 해석했다. 추가 agent 없음. 아래 v7 기록은 역사적 근거다.
이번 시작 fresh fetch에서 branch/origin은 모두
e964d824b6ffbff330df5d3c2410cd131b3aa4eb, 과학 worktree는 깨끗했다.
작업 위치는 기존 aft-pde-bessel-38969ad / probability-pde-solver-v1이다.
원래 OneDrive 폴더는 detached c43d8e0 및 사용자 파일을 그대로 보존한다.
AGENTS/인계 discovery 사본의 변경 전 SHA가 두 폴더에서 같음을 확인했다.

### 구현과 이론

- vector_fcc_rows.py: 원자열마다 실제 U=(u,v,w), affine gamma를 둔다.
  실제 원자열 반경이 바뀌므로 Bessel m=0과 모든 필요한 양의 mode를
  재평가한다. LJ pair/밀도/embedding/angular 후보를 refit하지 않았다.
  같은 행의 거리 및 완전 FCC bulk reference는 상쇄/정확한 배경으로 처리한다.
  각 site 환경을 합산한 뒤 nonlinear F 및 angular norm을 적용한다.
- 반대 row의 parity를 이용해 half-neighbor만 계산하되 per-site counting과
  odd/even moment 부호를 보존한다. analytic gradient/Hessian-vector 및
  3 translation gauge를 제거한 안정성 검사가 있다.
- 세 local displacement는 자유지만 transverse periodic cell vector는 고정이다.
  무한 직선 line 연구이지 finite source/3D specimen/물리 시간 dynamics가 아니다.
- vector_fcc_validation.py는 독립 direct atom 값/기울기/Hessian 확인 전용.
  canonical x방향 원자 합은 계속 infinite Poisson/Bessel이다.
- runner는 실제 최적화, reporter는 저장 상태 재검증이다. 둘을 혼동하지 않는다.
  root parameter SHA256는 계속
  9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655이다.

### 이미 실제 실행된 핵심 결과

24x24 ring6에서 기존 scalar 평형을 풀자 E=2.22178949에서 약1.4e-13eV/
line repeat로 내려가고 x winding 쌍과 모든 nonuniform vector displacement가
사라졌다. full force residual6.09e-8eV/L0, minimum curvature0.25513782.
0,4,25,50,0,-50,0MPa 정적 continuation을 실행했다. +50MPa shear 증가
약0.002076, 새로운 registry winding 없음. 이것은 주어진 가까운 반대 전위쌍의
이완/소멸이지 Al yield=0, 생성속도, 피로 검증이 아니다.

24x24 ring10 및32x32 ring8도 독립 실행했다. L-BFGS가 거의 무결정 상태의
작은 탄성 잔차를 느리게 줄이던 중 저장한 iteration100/80에서 동일 함수의
analytic Newton 단계로 이어갔다. 원래 partial 폴더/기록을 삭제하지 않았고
summary.completed=false와 continuation_case를 기록했다. *_newton 폴더의
parameter/source binding 및 actual force/stability 결과를 확인한다.

중요한 반례: 8x8, seed3b에서는 x winding=0이어도 vector fault가 남았다.
ring4 E=.575047975, force7.22e-8, lambda_min1.33224;
ring10 E=.585289091, force8.17e-8, lambda_min1.39727.
한 layer에 추가 shift≈(.5,.25579)L0, opening .01553L0, 다른 layer의
횡변위가 fixed cell shape를 보상한다. j 방향으로 거의 일정한 cell-spanning
registry fault이지 localized screw pair가 아니다. extra tau=(.5,.288675)
근처지만 독립 Al partial/core 검증은 아니다. 16x16에서 같은 seed protocol은
초기 scalar 단계부터 pair가 소멸했다. 따라서 동일 final core의 domain 수렴이라고
하지 않는다. 작은 영역의 fault 결과를 material residual plasticity로 승격하지 않는다.

이 반례를 보고 새 테스트의 잘못된 기대(“x winding=0이면 모든 성분=0”)를
수정했다. 힘 residual/positive Hessian은 통과하되 실제 남은 vector 상태를
확인하는 회귀 테스트로 만들었다. 기존 production/기존 테스트는 약화하지 않았다.
특히 fault의 x slip≈b/2라 scalar n=floor(x/b+.5)가 수치적으로 모호하다.
`layer_registry`는 slip 두 성분/격자 동치 0,tau,2tau 거리/normal change/
row dispersion/scalar partition margin을 저장한다. Raw case의 registry_shear는
오직 x projection이며 full vector plastic strain이 아니다. Consolidated CSV는
이를 x_projected_ 접두어로 표시한다. Production s=bn+xi 규약은 그대로다.

### 수치 검증과 한계

- 12개 고정점 direct ±384 atoms vs reciprocal: max abs value1.78e-14,
  gradient1.42e-13, Hessian5.12e-13. 서로 다른 channel units를 한 물리 floor로
  합치지 않는다. 최대14 reciprocal modes, 마지막 envelope8.28e-15.
- 실제 nonzero transverse checkpoint에서 ring8->16 energy difference
  .000196158eV, force1.10746e-4eV/L0; ring12->16 force1.11542e-6.
  v7의 transverse tail≈machine precision 결과는 여기로 전이되지 않는다.
- LJ m=0의 force/energy tail bound를 유도했고 actual shell difference로 시험했다.
  density/angular/nonzero modes 전체에 대한 rigorous bound라고 하지 않는다.
- 독립 3x3 bulk Fourier Hessian assembly와 actual vector Hessian product의
  차이2.84e-14. 24x24 최소값은 ring6 .255138 -> ring10 .257829 ->
  ring16/20/28 .257935. 그러나 최소값이 같아도 polarization이 y/z에서 x로
  교차해 전체 Hessian 오차가 남는다. full symbol ring20->28 norm차3.81e-6도 저장했다.
- 8x8 vector fault를 ring16에서 실제 재이완: E=.585337628233541,
  force1.36e-9, lambda_min1.39755637. 0->50->0MPa 제하 후 실제 x phase 변화
  1.20e-10L0, transverse 변화1.51e-10L0지만 scalar index 변화는 -.07654655로
  나올 수 있다(ring10에서는 +.07654655). 이것은 b/2 partition의 가짜 변화다.
  최종19 static states/7 completed cases; precursor2개는 partial로 보존했다.
  최종 숫자는 vector_core_v8/scientific_status.json 및 실행 summary가 기준이다.

### 실제 검증 결과

새 vector tests17개 + 기존 nonlinear17개 =34 passed/100.54s를 실제 실행했다.
app31 passed/268.00s,0 skipped. desktop smoke PASS(2.80s).
Full solver305 passed/1154.60s,0 skipped. 연구 계산과 일부 동시 실행한 실제 wall time이다.
py launcher가 없어 설치된 Python3.13 interpreter로 동등 명령을 실행했다.
모든 새 파일까지 포함한 git diff --check PASS. 새 JSON40개 및 압축 상태54개의
parse/unique site index 검사도 통과했다. 그림은 실제 저장 상태를 그려 육안 확인했다.
연구 결과 폴더는 약1.1MB이며 중단 precursor와 잘못된 x-index 판정의 반례도 남긴다.
변경은 새 vector 연구/검증/runner/reporter/results와 지침/인계/gate 문서뿐이다.
기존 core/PDE/Al static parameter/Ac/time/UI selector는 변경하지 않았다.
최종 commit SHA는 이 변경을 담은 git log와 사용자 최종 보고를 확인한다.
push 직전 remote 재확인과 fast-forward만 허용한다. main은 수정/병합하지 않는다.

### 다음에 이어갈 때

VECTOR_FCC_CORE_DERIVATION.md, vector_core_v8 결과, 실행/partial summary를 먼저 읽는다.
누락된 local normal/transverse force 문제는 실제 vector 이완으로 해결되는 방향을
확인했지만, 이를 actual Al strength/fatigue solver 완성으로 부르지 않는다.
독립 Al material/core/GSF 보정, finite line/source geometry, 실제 kinetics가 필요하다.
Mobility, physical seconds/Hz, A_c, production energy selector, UI gate는 그대로다.
임의 pin/holding force, stress threshold, 선 길이, mobility fitting으로 다음 결과를
만들지 않는다. 셀을 통과하는 fault와 고립된 partial/core 및 기존 결함의 소멸을 구분한다.

## 이전 상태: nonlinear_screw_v7 — 이상강도에서 결함 지배 강도로 가는 첫 비선형 단계

이 절이 아래 kinetics_loading_audit보다 최신이다. 최신 요청은
“이상강도를 임의 보정할 것이 아니라 실제 결함이 있는 금속의 강도로 나아갈
방법을 구현하자”였고 사용자가 ㄱㄱ로 진행을 승인했다. 추가 agent 없음.

시작 fresh fetch, branch HEAD와 origin은 모두
b8e91777a4e32d695ecca04f93eb4be338f12add, 작업 트리는 깨끗했다.
과학 worktree aft-pde-bessel-38969ad / probability-pde-solver-v1에서만 구현했다.
원래 폴더는 detached c43d8e0과 사용자 파일들을 그대로 보존한다.
원래 폴더와 과학 worktree의 AGENTS.md/이 문서는 변경 전 SHA가 일치함을
검사했다. 마지막 검증 후 discovery 문서만 같은 내용으로 동기화한다.

### 구현한 것 — 단순 계획이나 harmonic 재생이 아님

- nonlinear_fcc_screw.py: 무한 e1 원자열은 기존 Poisson/Bessel로 합산하고,
  각 단면 원자열에 독립적인 비선형 x slip을 허용한다.
  전체 site 환경을 합친 후 F(rho_i), D_r||Q_ri||²를 적용한다.
  scalar F'_bulk를 고정하지 않고 F''를 실제 nonlinear Hessian에 포함한다.
  scalar/각도 density kernel, LJ 계수 및 이전 fit은 전혀 바꾸지 않았다.
- FFT는 보존된 reciprocal row coefficient의 정확한 convolution이다.
  새로운 empirical PN/pinning law나 continuum stiffness splice가 아니다.
  rigid cut 비선형 에너지/힘은 독립 full-plane W_int와 일치하고,
  harmonic 한계는 기존 FCCScrewRowHessian의 symbol과 일치한다.
- per-cell u_i와 affine engineering shear gamma를 함께 풀 수 있다.
  응력 제어에서 G=E-tau*V_cell*gamma이며
  V_cell=N*b*d*h*L0³는 실제 원자 cell의 기하학적 부피다.
  임의 activation volume/길이를 도입한 것이 아니다.
  1MPa -> .00014659180163330851 eV/L0³; L0=2.863782463805517e-10m.
- gamma=0은 strain-control이다. 전위쌍이 있는 상태에서 zero stress와 다르다.
  initial zero-stress equilibrium를 먼저 구해 이후의 새 registry 변화와
  원래 있던 결함의 slip content를 분리했다.
- integer row relabeling u_i->u_i+b*k_i는 같은 원자 위치다.
  b-step만 있는 field는 core로 세지 않는다. plaquette winding ±1,
  실제 core 위치/간격, half-period margin을 계산한다.
- layer-bond b*k+xi 분해는 gamma=gamma_registry+gamma_intrabond를 만족한다.
  기존 production chi axial-strain bridge/P_n/SG flux와 같은 양이라 부르지 않는다.
- optimizer 종료와 실제 force convergence를 구별한다.
  필요시 analytic Newton-CG polish. material mode를 고치는 damping/clip 없음.
  scalar+affine Hessian의 최소 nongauge eigenvalue를 별도 검사한다.
  gauge shift=1을 잘못 최소 물리 고유값으로 보고하지 않는 test도 추가했다.

### 실제 실행한 ±50MPa / domain / unload

runner: python -m solver_v1.run_nonlinear_screw_study --sizes 24 32 48
656.0834 wall seconds, 실제 정적 평형44개.
하중 MPa: 0,4,15,25,50,25,0,-25,-50,-25,0.
무결함24x24와 같은 초기 screw pair의24/32/48 cross-section을 비교했다.
모든 경우 scalar/affine force balance 성공:
max field force residual2.401589e-11 eV/L0,
max shear decomposition error2.507218e-18.

중요: atomic sites를 더 촘촘하게 만든 것이 아니다. 영역의 원자열 수/주기
image 간격을 늘렸다. 실제 pair separation은 모두1.98408669165nm다.
seed의7.3b 값을 실제 defect separation으로 보고하지 않는다.

| domain | 초기 E [eV/b-repeat] | 초기 registry shear | +50MPa affine 증가 |
|---|---:|---:|---:|
|24x24|2.22178949153|.017010345436|.002083009966|
|32x32|2.26973529530|.009568319308|.002080243615|
|48x48|2.30313892346|.004252586359|.002078031671|

모든 단계에서 새 net registry 변화=0, 각 ±core 위치 동일.
제하 후 affine 변화 최대2.12538e-13, 새 registry 변화0.
이전부터 있던 registry content를 residual plasticity 생성으로 부르지 않는다.
이것은 athermal static unload이지 finite-T/time zero-stress hold가 아니다.

32->48 energy 차이 .03340363eV (약1.45%)가 남는다.
따라서 isolated dipole energy가 domain-converged라고 하지 않는다.
+50MPa incremental response 차이2.21194e-6 (약.106%).
±50MPa 구간에서 이동 없음은 세 domain에서 유지되지만 Peierls/yield threshold
전체를 구한 것은 아니다. 실제 Al 항복강도/피로 검증 아님.

24x24 scalar+affine 최소고유값:
0,+50,-50MPa에서 .252714014,.252835457,.252579727.
eigen residual<=3.20e-8. 물리 Hz가 아닌 static curvature다.

### 핵심 원인 발견: scalar 힘 평형은 full core 평형이 아니다

nonlinear_screw_transverse_audit.py로 고정한 실제 y,z 방향 힘을 독립 계산했다.
LJ radial derivative/analytic exponential Bessel derivative/STF polynomial
derivative를 사용한다. x-slip에서 소거되던 m=0 term을 반드시 포함했다.
m=0 scalar density force가 site별 다른 F'(rho_i)와 결합한다.
잘못 생략하면 vector extension이 틀린다.

실제 zero-stress 최대 omitted gradient:
24:2.2298290, 32:2.2610920, 48:2.28338234 eV/L0 (최대약1.28nN).
+50MPa에서도2.2146541,2.2448068,2.2661844가 남는다.
14->20 transverse ring change<=1.69867e-10 eV/L0.
실제 y,z를 미소 변경해 real-index nonlinear energy를 재평가한 독립
directional FD check error9.03391e-11.
rigid-plane normal-force, zero-force perfect bulk, total force balance도 시험했다.

즉 residual이 작고 Hessian이 양수인 것은 scalar 구속 부분공간에 한정된다.
큰 빠진 힘이 core 근처에 집중되며 잡음/급수 truncation이 아니다.
현재 해를 완전한 Al screw core나 “전위를 넣어서 실제강도 해결”로 채택하면 안 된다.
이것은 명확한 다음 수정 원인을 찾은 결과다. Mobility/온도/장벽/재료 계수를
낮춰서 해결한 것이 아니다.

### 수렴과 실제 데이터

- infinite-row default6rings/168rows/11modes, extra14rings 확인.
- 실제 core6->8rings energy<=1.60e-14eV,
  x-force<=2.71e-15eV/L0, stress<=4.22e-13MPa.
- nonlinear FFT vs real-index channel discrepancy2.22e-16.
- directional energy gradient FD1.24e-9, Hessian-vector7.66e-10.
- bounds는 per-channel 단위로 노출; finite extra-ring bounds를
  무한합 rigorous theorem으로 부풀리지 않는다.
- transverse audit는 별도로6/10/14/20rings를 실제 실행했다.

원시 자료:
results/fcc111_active_interface/nonlinear_screw_v7/
metadata.json, stress_history.csv, static_states.csv.gz,
stability.csv, tail_refinement.csv, relaxed_core_tail_refinement.csv,
domain_and_unload_summary.csv, transverse_force_audit.csv,
transverse_energy_derivative_validation.csv, execution_summary.json,
physical_and_numerical_status.json, static_defect_audit.png.

report_nonlinear_screw_study는 saved state를 읽고 omitted force를 실제 재평가한다.
plot을 만든 것을 새 time dynamics 실행으로 부르지 않는다.
저장된 모델 SHA는
9d00fbf54831c134fe9961e17f0603defdc7958bc3590743b2094264a31b6655.
기존 candidate C11/12/44~87.816/64.833/49.271GPa의 재료 fit 한계도 그대로다.

### 검증 기록

- 처음 신규 targeted: 1failed/10passed; rank3 STF trace가2.9976e-15인
  norm-scaled arithmetic error를 fixed2e-15로 검사한 새 test만 실패했다.
  8 eps ||tensor|| 기준으로 새 assertion을 바로잡고 재실행했다.
  기존 검증을 약화시킨 것이 아니다.
- 중간 scalar-only11passed/5.00s; 후속15passed/17.32s.
- 최종 신규 nonlinear+transverse targeted:17passed/9.99s.
- 최종 combined nonlinear/discrete/nonlocal:42passed/16.59s.
- app:31passed/89.23s, 제외0. 실제 Tk tests 포함.
- desktop --smoke: PASS, wall1.464s; LJ a0/kappa 역사값 그대로.
- full solver regression: **288passed /716.12s**, 실패0/제외0/exit0.
  실행 완료 출력과 .cache/nonlinear-screw-full.xml을 직접 확인했다.
- git diff --check 및 git diff --cached --check: 신규23개 변경 파일까지
  모두 통과했다. 최종 검사 결과는 validation_status.json에 기록한다.
- 최종 수치검증은 완료됐으며 본 문서가 포함된 checkpoint의 commit/push
  상태는 실제 git log/fetch/status로 확인한다. 물리 채택 gate는 미통과다.

### 다음 단계 — 무엇을 해결해야 하는지 이제 구체적임

1. 같은 비선형 per-site energy에서 y,z/vector core를 실제 이완한다.
   A_(iR)가 위치마다 달라지므로 fixed-coefficient FFT를 그대로 쓰지 않는다.
   Infinite row Bessel sum은 보존하되 state-dependent radius/STF/G=0 terms와
   tails, analytic vector forces/Hessian을 검증한다.
   단순 frozen harmonic correction만으로 큰 core reconstruction을 인증하지 않는다.
2. vector/normal-relaxed core와 partial splitting, slip path/GSF/source material
   검증 후 실제 이동 saddle/임계하중을 구한다. 현재 ±50MPa 비이동만으로
   실제 Al 강도/Peierls를 단정하지 않는다.
3. finite curved line/source/loop 및 결함의 실제 밀도·길이·경계조건이 필요하다.
   무한 직선 energy는 J/m이고 임의 line length/A_c로 activation eV를 만들지 않는다.
4. 그 후 coarse probability variables/finite-event energy/kinetic mobility를
   정당화한다. 현재 물리 M_a/M_s/t0/초/Hz는 여전히 없다.
5. Production PDE/energy registry/UI/calibration/Ac는 그대로이며
   solver/UI gate 미통과다. UI 재설계는 여전히 사용자 확인 전 진행하지 않는다.

상세 유도: solver_v1/NONLINEAR_FCC_SCREW_DERIVATION.md.
단계적 연구진전은 검증됐지만 실험적 Al 강도/피로의 완성 선언은 아니다.

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
