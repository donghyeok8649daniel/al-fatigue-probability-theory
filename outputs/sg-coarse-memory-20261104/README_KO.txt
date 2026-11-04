CRACK 상관기억 보존형 피로개시 이론 - 2026-11-04 연구 공개 묶음

현재본 읽기 순서
1. CORRELATION_MEMORY_FINAL_REPORT_KO.txt: 실제 이번 결과와 비교 범위.
2. CORRELATION_MEMORY_THEORY_KO.txt: 정리·증명·정확준비쌍·보류horizon·noise.
3. CORRELATION_MEMORY_PRIOR_ART_DELTA_KO.txt: 기존연구와의구체적인차이/한계.
4. CORRELATION_MEMORY_COMPARISON.png: 결과그림.

핵심 수치
두좌표전체marginal과meanenergy를정확히맞춘75commoninterior상태의준비쌍에서
한modelcyclefirst-opening확률차이0.00038402609384002008...를확인했다.
한주기response정보추가후같은cycle의보류2/4/8horizon수치식별폭이약605/249/132배
줄었다. 속도·정확도배수가아니며probe오차를넣으면이득이줄어든다.

이 묶음의 상태
실제Al/물리시간/시편균열event보정은미완료다. 세계최초수학·모든기존모형우위·
산업피로수명정확도·시험비절감완료를주장하지않는다. exactcertificate는captured
binary64finite표와그명시보정에관한것이며continuum/전체source오차인증은아니다.
이전후보/기각검토/초기이론파일은연구과정의이력이다. 현재결과는위현재본을우선한다.
원형선행prototype의비공개local계약파일은이공개묶음에포함되지않는다.

재현(저장소 root, Python+numpy+scipy+matplotlib)
python outputs/sg-coarse-memory-20261104/compare_correlation_memory.py
python outputs/sg-coarse-memory-20261104/compare_reachable_load_order.py
python outputs/sg-coarse-memory-20261104/verify_correlation_artifacts.py
python outputs/sg-coarse-memory-20261104/plot_correlation_comparison.py

입력과 재현 경계
SOURCE_OPERATOR_PAIR.npz는이전원N1SG에서캡처한136상태24단계전이와격자/초기법칙이다.
compare_correlation_memory는이캡처표를재사용하고원model의energy를평가한다.
compare_reachable_load_order는같은격자와capturedmask에서원SGgenerator를조립한다.
공개branch는origin/probability-pde-solver-v1의2026-11-04 fetch확인본을base로했다.
이연구의원작업HEAD는5f878370c1c74be99e1e9b727a8ab9a11872b174이고원격base는
c8ccff7553cbb2e66b0c37a76ee8a24808478824다. 해당두commit사이solver_v1차이는없다.
원작업probability_pde_2d.py에는다른용도의initialbasin예외와stresscallback변경이있었지만,
이번2개비교가사용한SGgenerator/energy경로는그변경을사용하지않는다.
공개전원격base의solver_v1사본에서새비교/산출물검사를실행해확인한다.

주의
실제초기집단측정자료가아닌선언수치준비집단이다. heldout은동일연산자의다른horizon이다.
여러columnassertion을독립실험수로집계하지않는다. 원source/GUI/실물/계정정보/팀개인정보는
이연구묶음에새로추가하지않는다. 원격공개는연구산출물보존이며물리검증승인이아니다.
