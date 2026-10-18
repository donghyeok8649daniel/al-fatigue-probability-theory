# 기호·단위 정리

기준 원격: `c8ccff7553cbb2e66b0c37a76ee8a24808478824`의 `probability-pde-solver-v1`. N=1 UI 모형은 `two_row_lj_reference`입니다. 실제 소재 결합 모형은 아직 선택되지 않았습니다.

- [기호 시트 Excel](symbols.xlsx): 기호·단위·코드 대응과 입력 근거 2개 탭.
- [기호 시트 CSV](symbols.csv): GitHub에서 조회·비교할 수 있는 동일 기호표.
- [지배방정식과 적용 예시 1개](governing_equations_and_example.md): 에너지→Smoluchowski→SG face→원래 binary absorber.
- [확인 입력 JSON](input_evidence.json): 프로파일별 확인값과 미정 실제 입력을 구별.
- [출처와 필요한 원자료](sources_and_missing_inputs.md).

## 표기 규칙

`*`는 native 정규화 변수입니다. 물리량은 SI로 기록합니다. 응력·탄성계수는 Pa, 길이는 m, 에너지는 J, 시간은 s로 통일합니다. native 숫자를 SI 값으로 읽으려면 해당 모형의 L0/E0/t0 근거가 필요합니다.

`P_N/rho`는 density, `p_k`는 적분된 raw bin 질량입니다. `nu_B`는 구성공간 bin 적분측도, `nu_F`는 full reference 질량, `V_geom`은 실제 재료 부피입니다. raw p에 측도·부피를 다시 곱하지 않습니다.

`M_i`는 mobility, `K_k`는 binary keep mask입니다. `B_SG`는 Bernoulli 함수입니다. `sigma`와 `sigma_LJ`, `E_x`와 `E0`, diffusion `D_i`와 pre-reference `D_KL`을 구분합니다.

## 확인한 값의 범위

입력 JSON에는 20개 source reference 기록을 채웠습니다. N=1 UI 기준과 일반 ModelParams 기본값은 별도 프로파일입니다. 선언된 kinetic 파일의 L0/E0도 별도 범위입니다. 실제 결합 입력 15개는 근거가 부족해 null로 유지했습니다. kinetic 파일은 calibrated=false이며 물리 mobility·온도·t0는 null입니다.

적용 예시는 기존 N=1 native 모형의 한 내부 SG face를 기호로 계산합니다. actual load/P0/grid/time/geometry 값을 새로 지정하지 않습니다.

이 공개 시트의 기호 48개와 기준 입력은 아래 출처 목록의 기존 공개 저장소 코드·문서에서만 정리했습니다. 로컬 보관 자료 및 미공개 분석 기호는 공개본에 포함하지 않습니다.
