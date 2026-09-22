# v9 출처와 원본 파일

## Durham B 자료

Dhamotharan, Wang, Hayes, Ramasse, Phillips, Major, Zoppi, Clark, Mendis (2025),
*Electronic and vibrational properties of interstitial clusters in degenerately
boron-doped silicon*, J. Phys.: Condens. Matter 37,465901.
[논문 DOI](https://doi.org/10.1088/1361-648X/ae191e),
[데이터 DOI](https://doi.org/10.15128/r1mw22v5500), **CC BY 4.0**.

- `durham_boron_2025.zip`: 받은 원본 archive. 내부 member별 hash는 `source_audit.json`.
- `durham_original/`: archive에서 그대로 꺼낸 `.cell`7개, `.xlsx`2개, README.
  워크북의 수식/캐시·측정값을 수정하지 않았다.
- `boron_formation_energies_2025.csv`, `boron_local_modes_2025.csv`:
  논문 Table1/2의 전사. PDF 화면으로 대조했다. per-B 값은 인쇄된 반올림 값을 보존한다.
  source table label의 charge와 neutral MACE 계산을 구별한다.
  DFPT dispersion은 Gamma–X 경로이며 본문은 해당 결함 modes를 nondispersive라고
  설명한다. MACE에서는 Gamma만 계산했다. 높은3NB modes만으로 전체 안정성을 인증하지 않는다.

Archive SHA256:
`5ed01b79a17bda9f59cd680bc965d8679bf25151a01d1e8c0ea81da4767fa9dd`.
저장소의 MD5 `8d0c05222168c004374692040e91bc5f`도 일치했다.

확인한 공개 PDF는
[White Rose 원문](https://eprints.whiterose.ac.uk/id/eprint/234205/7/Dhamotharan_2025_J._Phys.__Condens._Matter_37_465901%20%281%29.pdf),
1,531,855bytes, SHA256
`fc799544cdf08b355fbc0cedd3ba14f05a19cc3796e1059776ccb3d3e264b3ca`다.
표는 PDF 0-based page6/7(논문 인쇄 page6/7)에 있다. PDF는 local cache에 보존했다.

### 워크북 읽기

- Hall `BDS raw data!A1:H9`: 298K LDS/HDS 각각4회 측정. 측정 전류도 다르다.
- Hall `BDS values!B3,E3,H3`: 원래 외부 참조가 깨진 #REF! 셀.
  별도 `hall_recomputed_summaries.csv`는 raw sheet에서 재계산한 결과다.
- Hall `HDS Temperature data!A1:G33`:10–320K,32행.
  원본의 n 열은 Hall carrier number다. 이 p-type 자료를 전자 농도 n으로 입력하지 않는다.
  수치로 역산한 사용 전하상수는1.602e−19C다. 정확한 SI 상수와의 약0.011% 차이를
  재료 측정 오류로 해석하지 않는다. Hall factor/compensation은 별도 미확인이다.
- SIMS `Sheet1!B:N`: 시간(msec),B11 counts,Si29 counts. 감도 계수나 depth 환산이 없다.
  count 비를 원자분율/농도/이온화율로 바꾸지 않는다. 모든 측정행을 보존했다.

## 독립 MACE reference

[공식 배포](https://github.com/ACEsuit/mace-foundations/releases/tag/mace_mp_0b3),
`mace-mp-0b3-medium.model`, MIT,
SHA256 `2f2be696351ac9e94fbe01cdfb6f017679acdbd2db7645209ef55fec9826b012`.
모델 원본은 ignored cache에 보존하며 Git에 큰 weight 파일을 넣지 않는다.
실제 모델 버전·패키지·정밀도는 각 run_summary/summary에 기록했다.
학습 수준은 MPTrj의 PBE+U 계열이다. B-2025의 LDA와 같은 전자구조 방법이 아니다.
화학 조성은 입력되지만 별도의 초과 전자/전하 제약을 조절하는 모델로 사용하지 않았다.
공식 모델표와 b3 release를 재확인했다. b3의 phonon 수정 설명을 Si/B 정확도 승인으로
읽지 않는다. 전하/분극을 다루는 OMOL/Polar 후보의 실제 학습 범위는 분자 자료이며,
이번에 periodic charged Si/B/P 검증 결과를 확보한 것은 아니다. 실행하지 않은 후보와
확인한 범위는 `charge_model_screening.json`에 따로 기록했다.

## Cambridge Si DFT

[Bartók et al., PRX2018](https://doi.org/10.1103/PhysRevX.8.041048),
[Cambridge 데이터](https://doi.org/10.17863/CAM.65004).
기존 v5와 같은 archive/member hash와 같은2,475원자배열을 사용했다. 원본은 GAP
학습 자료이며 MACE 학습 중복은 독립 감사하지 않았다. PW91/PBE/미표기를 구분한다.
새 DFT 계산, fit, MD는 아니다.

추가 B/P 파괴·피로 1차 문헌의 확보 수준과 미확보 수치는
`dopant_fracture_source_screening.json`에 기록했다. 접근하지 못한 숫자는 null이다.
