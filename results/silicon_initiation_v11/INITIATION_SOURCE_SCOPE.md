# 개시·표면·경계 조건의 1차 문헌 범위 감사

확인일 2026-09-24. 아래 문헌은 물리 가설과 비교 조건을 고르는 근거이며,
현재 MACE 시편의 개시강도 또는 rate를 보정한 타깃 목록이 아니다.
전체 원문을 읽지 않은 자료는 초록 확인이라고 표시한다.

## 1. 실제 균열이 생기는 경로와 표면 상태

Godet, El Nabi, Brochard, Pizzagalli (2015),
*Surface effects on the mechanical behavior of silicon nanowires: Consequence on the brittle to ductile transition at low scale and low temperature*,
[DOI10.1002/pssa.201500001](https://doi.org/10.1002/pssa.201500001).

- 확인 범위: 출판사 초록과 서지정보.
- 저자들은 새로 절단한 표면, 300K annealing으로 재구성한 표면, 비정질화한
  표면을 MD에서 비교했다. 표면 상태에 따라 균열 형성과 소성 경로가 달라졌다.
  초록의 wedge 모양 개방은 해당 계산에서 사용한 구조 판별이다.
- 우리 연구에 대한 추론: 거리 기준 이웃 손실 하나만으로 모든 Si의 첫 균열을
  정의하기보다 새 분리면, 전단/재구성, 비정질화 후보를 구분해야 한다.
  이 논문의 MD 기구를 현재 MACE 또는 실제 웨이퍼의 검증된 기구로 전이하지 않는다.

Kang and Cai (2010), *Size and temperature effects on the fracture mechanisms of silicon nanowires: Molecular dynamics simulations*,
[DOI10.1016/j.ijplas.2010.02.001](https://doi.org/10.1016/j.ijplas.2010.02.001).

- 확인 범위: 출판사 초록. [110] Si nanowire의 MEAM/일정 변형률속도 MD.
- 해당 모델에서 크기와 온도에 따라 표면 균열 핵생성과 표면 전위 핵생성이
  다른 파손 경로를 만들었다. 작은 wire의 전단 파손도 보고했다.
- 우리 연구에 대한 추론: 현재 폭1.16×1.64nm bare 시편은 제조된 micrometre
  시편과 같은 기구/강도를 갖는다고 가정할 수 없다. 이 문헌의 크기 전이값을
  우리 모델의 물리 cutoff로 입력하지 않는다.

## 2. 성장 경로를 개시 장벽으로 잘못 가져오지 않는다

Huang, Zhang, Belytschko, Terdalkar, Zhu (2009),
*Mechanics of nanocrack: Fracture, dislocation emission, and amorphization*,
[DOI10.1016/j.jmps.2009.01.006](https://doi.org/10.1016/j.jmps.2009.01.006),
[저자 소속 기관 원문](https://sites.esm.psu.edu/wiki/_media/research:suz10:2009_jmps.pdf).

- 확인 범위: 원문 초록, 1-3절과 Fig1/2의 설명.
- 제목/본문에 formation이라는 단어가 있지만, 주요 반응 경로 비교는
  central crack의 준안정 상태와 crack-tip bond의 절단/회복, 선단의 전위 방출과
  비정질화를 다룬다. 그 barrier를 무균열 시편의 최초 개시 장벽으로 쓰지 않는다.
- 응력제어/변형률제어와 유한 크기에 따라 에너지 경로가 바뀐다는 결과는
  비교 조건을 엄밀히 맞춰야 한다는 근거다. 논문 SW 재료 상수나 장벽을
  현재 MACE에 복사하지 않는다.

## 3. 경계 조건의 영향은 비교한 조건 안에서 읽는다

Xu and Kim (2020), *Role of boundary conditions and thermostats in the uniaxial tensile loading of silicon nanowires*,
[DOI10.1016/j.commatsci.2020.109636](https://doi.org/10.1016/j.commatsci.2020.109636).

- 확인 범위: 출판사 초록과 highlights. 세 MEAM을 써서 periodic/fixed ends와
  NVE/Langevin 조건을 비교했다. 그 계산에서 영률·강도의 의존성은 작았고
  failure strain 및 기구에는 조건별 차이가 있었다.
- 모든 경계조건 효과가 크거나, 모두 무시해도 된다는 보편 결론이 아니다.
  현재 3.184GPa의 초기 길이 반력은 실제 저장 원자력에서 얻은 값이다.
  이 문헌으로 그 숫자를 검증했다고 주장하지 않는다.

## 4. GPa 단위 자체가 수치 오류라는 뜻은 아니다

Tang 외 (2012), *Mechanical Properties of Si Nanowires as Revealed by in Situ Transmission Electron Microscopy and Molecular Dynamics Simulations*,
[DOI10.1021/nl204282y](https://doi.org/10.1021/nl204282y).

- 확인 범위: 출판사 초록과 supporting-information 목록.
- in situ TEM 인장과 굽힘, MD를 비교했고 인장에서 GPa 수준의 크기 의존 강도를
  보고했다. 이는 nanospecimen 강도를 실제 wafer의 MPa 범위에 임의 배율로
  낮춰서는 안 된다는 비교 근거다.
- 현재 시편의 재료 모델·표면·치수·하중과 일치하는 검증 데이터는 아니다.
  숫자를 맞추는 목적의 barrier 축소, correlation-area 조정, 이미 존재하는
  균열 길이의 역맞춤을 사용하지 않는다.

## 5. 논문에서 probability라는 이름을 써도 첫 통과 확률은 아닐 수 있다

Zare Pakzad, Nasr Esfahani, Alaca (2023),
*The Role of Native Oxide on the Mechanical Behavior of Silicon Nanowires*,
[DOI10.1016/j.mtcomm.2022.105002](https://doi.org/10.1016/j.mtcomm.2022.105002),
[대학 저장소 accepted manuscript](https://eprints.whiterose.ac.uk/id/eprint/195147/).

- 확인 범위: 원문 방법, Table1, 실패 지표 정의. PDF30쪽 중 파일17쪽/본문16쪽의
  식은 이미지로도 확인했다. PDF SHA256
  `248668b4b8a6650ad7ae99efd82abc2c585852ed41451eb8441ead9e2233aa65`.
- 두 고전 포텐셜과 native oxide를 비교한 MD다. 설정은10K와5e9s^-1이므로
  실온 저주파 피로나 현재 모델의 물리 clock 타깃이 아니다.
- 본문에서 ductile failure probability라고 부르는 p는
  **p=(W0-Wfail)/W0**, 즉 파손 시 폭 감소의 정규화 지표로 정의된다.
  이 양을 시편 ensemble의 균열 개시확률 또는 우리 PDE의 흡수 질량으로 쓰지 않는다.
- 표면층/단면 규약과 재료 포텐셜에 대한 비교 근거는 얻을 수 있지만,
  문헌의 necking 판별값을 새 개시 absorbing boundary로 채택하지 않는다.

## 6. 현재 실제 채택 상태

위 자료 중 새 물리 보정값으로 채택한 수치: **없음**.
문헌 범위 확인은 새로운 MD/DFT 실행이 아니다. 현재 정량 감사에 실제 사용한
원자료는 제공된 Tsuchiya 열산화층 논문 표, 공개 QE/PBE Si/O1159상태,
공개 CP2K Si/O/H 자료의 사전 선택 구조다. 이들 역시 실제 개시 장벽/동역학을
직접 보정한 것은 아니며 각 출처와 원자료 SHA는 별도 기록으로 유지한다.
