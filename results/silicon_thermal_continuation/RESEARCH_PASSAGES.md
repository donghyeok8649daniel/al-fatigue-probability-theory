# 기존 논문에서 우리 연구에 사용할 부분

이미 확보한 여덟 편에서 필요한 부분만 추렸다. 아래 내용은 원문을 길게 옮긴 것이 아니라
해당 쪽·표·그림의 요지와 사용 목적을 정리한 것이다. 직접 인용 문구와 검수 기록은
[기존 인용 목록](../silicon_force_thermal_comparison/literature_and_harmonic/citation_addendum.json),
[PDF 검수](../silicon_force_thermal_comparison/literature_and_harmonic/citation_proof.json),
[출판사 본문 검수](../silicon_force_thermal_comparison/literature_and_harmonic/publisher_quote_checks.json)에 있다.
이 메모에는 새로 검색한 후보 논문을 넣지 않았다.

## 먼저 꺼내 쓸 근거

| 연구에서 확인할 것 | 먼저 볼 자료 | 뽑아 사용할 부분 |
|---|---|---|
| 베셀 기반 Si 에너지 후보가 물성을 함께 맞추는가 | Bartók 등 | 에너지·힘·탄성·진동·표면과 파괴 경로를 함께 검증하는 부분 |
| 도핑에 따라 탄성이 얼마나 바뀌는가 | Jaakkola 등 | Table III의 농도별 탄성계수와 온도 적합식 |
| 전하 효과와 실제 도펀트 원자 효과를 구분하는가 | Noda 등 | Methods의 excess-charge 설정과 [111] 이상강도 결과 |
| 산화층 응력 분담과 공정 손상을 구분하는가 | Tsuchiya 등 | 산화·제거 시편의 강도와 FEM의 상별 국소 응력 |
| 웨이퍼 최종 파단을 어떻게 대조할 것인가 | Chen 등 | 시험 형상·도핑군 조건·Weibull 특성강도·군간 통계 |
| 결정 방향과 습도를 시험 조건에 어떻게 넣는가 | Ikehara·Ando | 방향별 응력 환산, 온습도 조건과 환경 영향 |
| 파단 전에 어떤 신호를 기록할 것인가 | Muhlstein 등 | 공진 주파수 변화, 측정 간격과 최종 파단의 구분 |

## 1. 원자 에너지 후보를 평가하는 기준

**Bartók, Kermode, Bernstein, Csányi (2018)**,
*Machine Learning a General-Purpose Interatomic Potential for Silicon*.
[원문](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.041048), DOI 10.1103/PhysRevX.8.041048.

- **볼 곳:** §I.A, 저널 041048-2; §III/IV의 검증 결과.
- **사용할 요지:** 소수 물성값을 맞춘 것만으로 원자 모델이 옳다고 판단할 수 없다. 에너지와 힘,
  진동·탄성·표면 및 파괴에 관련된 구성을 함께 대조해야 한다.
- **우리 연구에서 할 일:** 베셀 후보의 같은 파라미터로 이 관측량을 계산하고, 적합에 쓴 자료와
  제외한 자료의 오차를 함께 공개한다. GAP의 성공을 우리 베셀 후보의 검증으로 옮기지 않는다.
- **자료:** 논문의 원자 자료 DOI [10.17863/CAM.65004](https://doi.org/10.17863/CAM.65004).
  실제로 사용할 구조·에너지·힘과 계산 조건을 확인한 후 대응시킨다.

## 2. 도핑에 따른 탄성 변화

**Jaakkola, Prunnila, Pensala, Dekker, Pekko (2014)**,
*Determination of Doping and Temperature Dependent Elastic Constants of Degenerately Doped Silicon from MEMS Resonators*.
[원문](https://arxiv.org/pdf/1401.1363), arXiv:1401.1363.

- **볼 곳:** 원고 PDF 7쪽의 도핑 논의, Table III 8쪽, Eq. (8).
- **사용할 수치:** 25°C에서 P 농도 4.1→7.5×10¹⁹ cm⁻³의 C₁₁ 163.0→161.4,
  C₁₂ 65.4→66.1, C₄₄ 79.2→78.5 GPa.
- **조건:** B·As·P 시편을 공진으로 대조한 자료이며 농도는 저항에서 환산한 캐리어 농도다.
  온도 범위는 −40~85°C다. Sb 값은 이 자료에서 확보하지 않았다.
- **우리 연구에서 할 일:** 농도별 평형 곡률과 탄성 모드를 대조한다. 몇 퍼센트의 탄성 변화를
  같은 크기의 파단강도 변화나 피로 수명 변화로 환산하지 않는다.

## 3. 전하를 넣은 이상강도와 실제 도핑의 구분

**Noda 등 (2023)**, *Electronic strengthening mechanism of covalent Si via excess electron/hole doping*.
[원문](https://www.nature.com/articles/s41598-023-42676-z), DOI 10.1038/s41598-023-42676-z.

- **볼 곳:** Methods, 출판 PDF 2쪽; Fig. 7/8, 5–6쪽.
- **사용할 수치:** [111] LDA 이상강도는 무전하 20.98 GPa, 전자 5×10²¹ cm⁻³에서 11.19 GPa,
  정공 같은 농도에서 27.23 GPa.
- **조건:** excess electron/hole와 보상 배경을 넣은 결함 없는 결정의 정적 DFT다.
  B/P/As 등의 실제 불순물 원자를 넣은 계산이 아니다.
- **우리 연구에서 할 일:** 전하만의 효과와 원자 종·자리·국소 환경의 효과를 분리한다.
  이 큰 전하 밀도의 경향을 10¹⁹ cm⁻³ 웨이퍼나 MPa 개시값으로 외삽하지 않는다.

## 4. 산화층과 Si의 하중 분담

**Tsuchiya, Miyamoto, Sugano, Tabata (2016)**,
*Fracture behavior of single crystal silicon with thermal oxide layer*.
[원문 정보](https://repository.kulib.kyoto-u.ac.jp/handle/2433/258961), DOI 10.1016/j.engfracmech.2015.08.029.
사용자가 제공한 PDF와 같은 논문이다.

- **볼 곳:** 제공 원고 5쪽 Fig. 5 논의, Tables IV/V 14쪽과 FEM 설명.
- **사용할 수치:** 산화층 0/50/100/200 nm의 평균 공칭 파단강도 4.09/3.50/3.77/3.27 GPa.
  200 nm 산화 후 제거한 시편은 2.55 GPa.
- **조건:** (100)〈110〉, 약 4×5 µm 단면, 1100°C 산화 공정, 시험 26°C·50% RH.
  Si와 산화층 응력값은 실험 하중을 넣은 FEM의 특정 위치 결과다. 산화층 압축응력
  460 MPa는 FEM 입력이며 독립 잔류응력 측정값으로 쓰지 않는다.
- **우리 연구에서 할 일:** Si/SiO₂ 응력 분담, 표면 공정과 산화 후 남는 결함을 함께 검사한다.
  50→100 nm에서 강도가 올라가는 결과도 보존한다. 두께에 따라 위험도가 단조 증가하도록
  임의 함수를 만들지 않는다. 이 표는 첫 균열 시간이 아니라 단조 최종 파단이다.

## 5. 실제 웨이퍼에서의 최종 파단 통계

**Chen 등 (2025)**,
*Fracture strength of Czochralski silicon wafers: Statistically rigorous measurement and dopant effects*.
[원문](https://pubs.aip.org/aip/jap/article/138/5/055703/3357407/Fracture-strength-of-Czochralski-silicon-wafers),
DOI 10.1063/5.0270955.

- **볼 곳:** Tables I/III와 §III.C.
- **사용할 수치:** 여섯 군의 Weibull 특성 파단강도 2.07–2.12 GPa, 각 군 50개.
  특성강도의 표준오차 0.06–0.08 GPa. 군간 차이는 이 시험에서 p>0.05다.
- **조건:** 300 mm (100) Cz 웨이퍼, 두께 775 µm에서 17×17 mm 시편,
  ball-on-ring 0.2 mm/min. 산소 농도는 약 7.9×10¹⁷ cm⁻³다.
- **우리 연구에서 할 일:** 시험 형상과 표면·도핑 조건을 일치시켜 최종 파단 분포를 대조한다.
  특성강도를 산술평균으로 읽거나, 통계적 비유의성을 모든 도핑 조건의 효과 부재로 해석하지 않는다.
  시편별 첫 개시·회복 원자료는 아직 없다.

## 6. 결정 방향별 피로와 응력 환산

**Ikehara, Tsuchiya (2016)**,
*Crystal orientation-dependent fatigue characteristics in micrometer-sized single-crystal silicon*.
[원문](https://www.nature.com/articles/micronano201627), DOI 10.1038/micronano.2016.27.

- **볼 곳:** Fig. 4, Tables 2/3, Methods; Discussion의 실험 피로강도 비교 및 Conclusion.
- **사용할 부분:** 현재 시험 116개, 23°C·50% RH·R=−1에서 다섯 방향/형상군의
  최종 파단 반복 수를 비교했다. 최대주응력이 저자들의 비교에서 좋은 기준이었다.
- **우리 연구에서 할 일:** 결정 방향, 시편 형상, 실제 국소 응력 환산을 유지한다.
  처짐각을 MPa로 직접 읽지 않는다. 반복 표시된 이전 자료를 새 독립 시편으로 세거나,
  ramp의 등가 반복 수를 일정 진폭 원자료로 바꾸지 않는다.
- **경계:** 논문의 균열 성장 지수 n=27을 첫 개시 법칙에 넣지 않는다.
  최대주응력이 모든 환경과 형상에 보편적인 개시 기준이라는 주장도 하지 않는다.

## 7. 온도와 습도를 함께 기록해야 하는 근거

**Ando, Shikida, Sato (2016)**,
*Effect of Temperature and Humidity on Degradation of Single-Crystal Silicon Microbeam in MEMS Resonator*.
[출판사 초록](https://sensors.myu-group.co.jp/article.php?ss=1257), DOI 10.18494/SAM.2016.1257.

- **볼 곳:** 현재 확보한 출판사 초록.
- **사용할 요지:** 높은 습도에서 피로강도가 감소하고, 높은 온도에서는 열화가 억제되는
  결과를 제시한다. 온도와 습도의 상호작용을 시험에서 분리해 확인할 이유다.
- **우리 연구에서 할 일:** T와 RH를 함께 기록·대조한다. 현재 초록만 확인했으므로
  시험의 모든 숫자나 습도별 수명식을 만들어 채우지 않는다.

## 8. 파단 전 관측과 첫 균열 시점의 구분

**Muhlstein, Brown, Ritchie (2001)**, *High-Cycle Fatigue of Single-Crystal Silicon Thin Films*.
[원문](https://www2.lbl.gov/ritchie/Library/PDF/Muhlstein12172001.pdf), DOI 10.1109/84.967383.

- **볼 곳:** §III.A, 저널 595–596쪽/PDF 3–4쪽, Figs. 3–5.
- **사용할 부분:** 12개, 약 20 µm 두께 (110) 단결정 Si 미세 시편,
  30±0.1°C·50±2% RH·40/50 kHz·R=−1에서 4–10 GPa와 10⁶–10¹¹회의 최종 파단 기록.
  파단 전에 공진 주파수가 감소했다. 확인 간격은 1–5분이다.
- **우리 연구에서 할 일:** 하중·반복 수·영상과 동기화된 공진/강성 변화를 보조 관측으로 기록한다.
  공진 변화만으로 첫 균열의 위치·시간·비회복 여부를 확정하지 않는다.

## 이 자료로 아직 채울 수 없는 것

여덟 편은 원자 모델·탄성·이상강도·산화층·최종 파단과 환경 조건을 검증하는 근거다.
우리가 필요한 **개별 첫 균열 형성 시점과 같은 조건에서의 비회복 확인**을 함께 제공하는 자료는
아직 확보하지 못했다. 따라서 이 자료를 개시 확률이나 비회복 전이율의 정답으로 적합하지 않는다.
수치의 전사·조건·원래 표는 [52행 대조표](../silicon_force_thermal_comparison/literature_and_harmonic/numeric_comparison.csv)와 함께 확인한다.
