# Open Results Sharing Platform

## Status

**PLATFORM DESIGN — DOES NOT MODIFY THE ACTIVE 1D FATIGUE THEORY**

The intended website is not merely a private cloud storage service for individual simulation projects. Its primary scientific purpose is to let users **publish simulation and fatigue-test results as open, reusable research data**.

A useful analogy is:

- GitHub shares source code and its history;
- this platform shares tensile/fatigue simulation and experimental results together with enough setup metadata to reproduce and compare them.

Strictly speaking, published numerical and experimental results are **open data**, not source code. The software itself may be open source, while contributed datasets need an explicit data license.

## Core public object: a result contribution

The fundamental public unit should be a versioned **Result Contribution**, not an arbitrary uploaded ZIP file.

A contribution should contain at minimum:

```text
contribution/
├── manifest.json
├── setup.json
├── provenance.json
├── results/
│   ├── summary.json
│   ├── nodes.csv
│   ├── elements.csv
│   └── optional probability outputs
├── experiment/
│   ├── raw_measurements.csv        # optional
│   ├── calibration.json            # optional
│   └── machine_metadata.json       # optional
├── figures/
├── README.md
└── LICENSE-DATA
```

The contribution may contain simulation-only data, experiment-only data, or a matched simulation-plus-experiment pair.

## Reproducibility metadata

A public result is scientifically useful only if its origin can be reconstructed. Therefore `provenance.json` should record, when applicable:

- repository / software version;
- exact solver Git commit;
- theory Git commit;
- project schema version;
- operating-system / architecture information relevant to reproducibility;
- C compiler and relevant Python package versions;
- whether the result came from the official solver or an external solver;
- checksum of every primary data file;
- date/time of generation;
- contributor identifier or chosen public attribution;
- for experiments: tester firmware version and calibration identity.

A result produced by an unknown code version may still be uploaded, but the website should visibly label its reproducibility status.

## Public result states

Recommended publication states are:

1. `draft` — visible only to the contributor;
2. `public-unverified` — public, schema-valid, but not independently reproduced;
3. `reproduced` — the official platform worker regenerated the simulation from the submitted setup and matched the submitted numerical result within stated tolerances;
4. `experiment-attached` — real fatigue-test data is attached;
5. `curated` — optional later status after additional review/quality checks.

`reproduced` must not mean that the physical theory is experimentally validated. It means only that the submitted computational artifact was reproducible from its declared inputs and code version.

## Website model

The main website should behave more like an open research repository than a personal dashboard.

Users should be able to:

- browse public contributions without an account;
- search by material, specimen geometry, loading condition, frequency, stress amplitude, mean stress, solver/theory version, and experiment/simulation type;
- open a contribution page and inspect its setup and results;
- download the complete contribution bundle;
- download individual CSV/JSON files;
- reproduce a compatible simulation using the official solver;
- compare several public results on common axes;
- publish their own simulation result from the desktop application;
- publish fatigue-tester results from the local tester agent;
- attach experimental results to an existing simulation setup;
- create a derived/forked contribution while retaining provenance to the original result;
- cite a stable public version of a contribution.

## Contribution page

A contribution page should clearly separate **inputs**, **outputs**, and **interpretation**.

Example:

```text
Al 1D Tensile Fatigue — 100 MPa amplitude, 20 Hz
Contributor: ...
License: CC BY 4.0
Status: reproduced / experiment attached

SETUP
- material
- geometry
- sigma_mean
- sigma_amplitude
- frequency
- waveform
- solver commit

SIMULATION
- displacement
- strain
- axial stress
- later: P(lambda,t), Qc, first-passage outputs

EXPERIMENT
- measured force/stress
- displacement
- DCPD
- crack-initiation observation
- temperature / machine state when available

FILES
- setup.json
- raw CSV
- processed CSV
- figures
- provenance

RELATED
- parent/fork contribution
- same setup at other frequencies
- same material at other amplitudes
```

## Open-data licensing

Software licensing and result-data licensing should be separated.

Examples of possible dataset choices are:

- CC BY 4.0 — reuse is allowed with attribution;
- CC0 — contributor dedicates the dataset as broadly as possible to the public domain;
- another explicitly supported open-data license if needed later.

The website should not silently publish data without an explicit contributor license selection and confirmation.

## Immutable published versions

A scientific result must not silently change after others have cited it.

Recommended rule:

- drafts are editable;
- once a result version is published, its primary files and checksums are immutable;
- corrections create a new version, e.g. `v1 -> v2`;
- old versions remain accessible;
- the contribution page points to the latest version while preserving the full version history.

This is analogous to tagged releases rather than mutable cloud documents.

## Fork / derived-result model

A user should be able to take a public setup or result and create a derived contribution.

Examples:

- same geometry, different stress amplitude;
- same loading, finer FEM discretization;
- same simulation setup, real fatigue experiment attached;
- same experimental dataset, processed with a newer probability-theory version.

The platform should store explicit provenance links:

```text
parent_contribution_id
parent_version
change_summary
```

This makes result families traceable rather than producing disconnected duplicate files.

## Automatic validation on upload

Before public publication, the server should at minimum verify:

- required schema fields;
- units and physical dimensions;
- finite numeric values;
- CSV column definitions;
- declared file checksums;
- solver/theory version metadata;
- absence of executable content in dataset-only uploads unless explicitly supported;
- consistency between setup metadata and result metadata.

For official simulations, the strongest useful check is a **reproduction job**:

$$
\text{submitted setup}
\rightarrow
\text{official versioned solver}
\rightarrow
\text{recomputed result}
\rightarrow
\text{numerical comparison}.
$$

If the comparison passes the defined tolerances, the site can assign a `reproduced` badge.

## Experimental contribution quality

Experimental uploads need additional metadata because two CSV files are not comparable merely because both contain stress and cycles.

Recommended required or strongly encouraged metadata includes:

- specimen material and preparation;
- specimen dimensions;
- loading waveform;
- mean/amplitude stress or force;
- frequency;
- load-cell calibration;
- displacement-sensor calibration when used;
- crack detection method;
- fatigue tester model / firmware version;
- sampling rate;
- test termination reason;
- raw data availability.

The platform should distinguish raw measurements from post-processed results.

## Public aggregation without destroying raw provenance

The website can later generate community-level plots such as:

- crack-initiation probability versus cycles;
- frequency comparisons;
- stress-amplitude sweeps;
- simulation-versus-experiment error distributions;
- comparisons between theory versions.

However, every aggregate point must remain traceable to its source contribution IDs and versions.

Therefore:

$$
\boxed{\text{aggregate database} \text{ is derived from } \text{immutable public contributions}}
$$

and not the other way around.

## Relationship to the desktop application

The desktop application should eventually expose two distinct actions:

```text
Save locally
Publish result
```

`Publish result` should:

1. validate the local project/result schema;
2. create the contribution manifest and provenance records;
3. show exactly which files and metadata will become public;
4. require the user to choose/confirm a data license;
5. upload to the website;
6. receive a stable contribution ID/version;
7. optionally request server-side reproduction.

Local use must remain possible without publication.

## Relationship to the fatigue tester

The local tester agent should be able to package a completed fatigue test as an experimental contribution, but publication must still be an explicit user action.

The real-time control and safety boundary remains local to the tester MCU/hardware. Public-data sharing is a post-test data/provenance function, not part of the safety-critical feedback loop.

## Recommended first implementation

Before building the full website, implement the smallest open-result unit:

1. define `public_result_manifest.json` schema v1;
2. define simulation-result and experiment-result required fields;
3. add `Export contribution` to the desktop application;
4. implement a local validator that checks a contribution bundle;
5. create a minimal website/API that accepts a validated contribution and exposes a public read-only page;
6. add downloading and searching;
7. add server-side reproduction badges;
8. add fatigue-tester experiment publication;
9. add result comparison and public aggregate plots.

The central design rule is:

$$
\boxed{
\text{publish the inputs + raw/derived outputs + exact provenance, not just screenshots or summary numbers}
}
$$

---

# 공개 결과 공유 플랫폼

## 상태

**플랫폼 설계 — 활성 1D 피로 이론 자체를 변경하지 않음**

목표 웹사이트는 사용자의 시뮬레이션 파일을 개인 클라우드에 저장하는 서비스가 아니다. 핵심 과학적 목적은 사용자가 **시뮬레이션 결과와 실제 피로시험 결과를 공개하고, 다른 사람이 재사용·재현·비교할 수 있는 오픈 연구 데이터로 기여하도록 하는 것**이다.

비유하면 다음과 같다.

- GitHub는 소스코드와 그 변경이력을 공유한다.
- 이 플랫폼은 인장/피로 시뮬레이션 및 실험 결과를 재현에 필요한 설정정보와 함께 공유한다.

엄밀히 말하면 공개되는 수치/실험 결과는 `오픈소스`라기보다 **오픈 데이터**다. 프로그램 코드는 오픈소스 라이선스를 쓰고, 사용자가 기여한 데이터에는 별도의 데이터 라이선스를 명시해야 한다.

## 공개의 기본 단위: Result Contribution

사용자가 임의 ZIP을 올리는 구조보다, 버전이 지정된 **Result Contribution**을 공개 기본 단위로 삼는 것이 좋다.

최소 구조 예시는 다음과 같다.

```text
contribution/
├── manifest.json
├── setup.json
├── provenance.json
├── results/
│   ├── summary.json
│   ├── nodes.csv
│   ├── elements.csv
│   └── optional probability outputs
├── experiment/
│   ├── raw_measurements.csv        # 선택
│   ├── calibration.json            # 선택
│   └── machine_metadata.json       # 선택
├── figures/
├── README.md
└── LICENSE-DATA
```

하나의 contribution은 시뮬레이션만 포함할 수도 있고, 실험만 포함할 수도 있고, 동일 조건의 시뮬레이션+실험 쌍을 포함할 수도 있다.

## 재현성 metadata

공개 결과는 생성 경로를 복원할 수 있어야 과학적으로 의미가 있다. 따라서 `provenance.json`에는 가능한 경우 다음을 기록해야 한다.

- 저장소/소프트웨어 버전;
- 정확한 solver Git commit;
- theory Git commit;
- 프로젝트 schema version;
- 재현성에 필요한 OS/architecture 정보;
- C compiler 및 관련 Python package version;
- 공식 solver 결과인지 외부 solver 결과인지;
- 주요 데이터 파일 checksum;
- 결과 생성 시각;
- contributor 식별자 또는 공개 attribution 이름;
- 실험이면 시험기 firmware version과 calibration 식별정보.

코드 버전을 알 수 없는 외부 결과도 업로드할 수는 있지만, 웹사이트에서 재현성 상태가 낮다는 것을 분명히 표시해야 한다.

## 공개 상태

권장 공개 상태는 다음과 같다.

1. `draft` — 기여자만 볼 수 있음;
2. `public-unverified` — schema는 통과했지만 독립 재현되지 않은 공개 결과;
3. `reproduced` — 공식 platform worker가 동일 setup과 선언된 software version으로 다시 계산하여 정해진 tolerance 안에서 결과를 재현함;
4. `experiment-attached` — 실제 피로시험 데이터가 연결됨;
5. `curated` — 추후 추가적인 검토/품질관리 후 선택적으로 부여.

여기서 `reproduced`는 이론이 실제 물리와 맞다는 뜻이 아니다. **제출된 계산 결과가 선언한 입력과 코드로 다시 생성됨**을 의미할 뿐이다.

## 웹사이트 성격

웹사이트는 개인 대시보드보다 공개 연구 저장소에 가깝게 동작해야 한다.

사용자는 다음을 할 수 있어야 한다.

- 로그인 없이 공개 contribution 탐색;
- 재료, 시편형상, 하중조건, 주파수, 응력진폭, 평균응력, solver/theory version, 실험/시뮬레이션 종류로 검색;
- contribution 페이지에서 설정과 결과 확인;
- 전체 contribution bundle 다운로드;
- CSV/JSON 개별 다운로드;
- 호환되는 시뮬레이션을 공식 solver로 재실행;
- 여러 공개 결과를 동일 축에서 비교;
- 데스크톱 앱에서 자기 시뮬레이션 결과 공개;
- 실제 피로시험기 agent에서 실험결과 공개;
- 기존 시뮬레이션 조건에 실제 실험결과 연결;
- 공개 결과를 기반으로 파생/fork contribution 생성하면서 원본 provenance 유지;
- 안정된 특정 version을 인용.

## Contribution 페이지

한 페이지에서는 **입력 / 출력 / 해석**을 명확히 분리해야 한다.

예시:

```text
Al 1D Tensile Fatigue — 100 MPa amplitude, 20 Hz
Contributor: ...
License: CC BY 4.0
Status: reproduced / experiment attached

SETUP
- material
- geometry
- sigma_mean
- sigma_amplitude
- frequency
- waveform
- solver commit

SIMULATION
- displacement
- strain
- axial stress
- 추후 P(lambda,t), Qc, first-passage 결과

EXPERIMENT
- measured force/stress
- displacement
- DCPD
- crack initiation 관측값
- 가능하면 temperature / machine state

FILES
- setup.json
- raw CSV
- processed CSV
- figures
- provenance

RELATED
- parent/fork contribution
- 같은 setup의 다른 frequency
- 같은 material의 다른 amplitude
```

## 오픈 데이터 라이선스

프로그램 라이선스와 결과 데이터 라이선스를 분리해야 한다.

예시 선택지는 다음과 같다.

- CC BY 4.0 — attribution을 유지하면 재사용 가능;
- CC0 — 가능한 범위에서 데이터를 공공영역에 가장 넓게 개방;
- 추후 필요하면 다른 명시적 오픈 데이터 라이선스 지원.

사용자가 어떤 라이선스로 공개되는지 명시적으로 선택하고 확인하지 않은 채 결과를 자동 공개하면 안 된다.

## 공개된 version의 불변성

다른 사람이 이미 인용한 과학 결과가 조용히 바뀌면 안 된다.

권장 규칙은 다음과 같다.

- draft는 수정 가능;
- 한 번 공개한 result version의 주요 파일과 checksum은 immutable;
- 수정하면 `v1 -> v2`처럼 새 version 생성;
- 과거 version도 계속 접근 가능;
- contribution 페이지는 최신 version을 보여주되 전체 version history를 보존.

즉 일반 클라우드 문서보다 tagged release와 비슷한 성격이다.

## Fork / 파생 결과

공개된 setup이나 결과를 기반으로 다른 사용자가 파생 contribution을 만들 수 있어야 한다.

예를 들면:

- 같은 시편에서 stress amplitude만 변경;
- 같은 하중에서 FEM mesh만 세분화;
- 같은 simulation setup에 실제 피로시험 추가;
- 같은 실험 raw data를 새로운 확률이론 version으로 후처리.

플랫폼에는 다음 provenance를 명시적으로 저장한다.

```text
parent_contribution_id
parent_version
change_summary
```

그러면 중복 파일이 흩어지는 대신 결과 계보를 추적할 수 있다.

## 업로드 시 자동검증

공개 전 서버가 최소한 다음을 검사해야 한다.

- 필수 schema field;
- 단위와 물리차원;
- NaN/inf가 아닌 수치;
- CSV column 정의;
- 파일 checksum;
- solver/theory version metadata;
- dataset-only upload 안에 허용되지 않은 executable이 포함되지 않았는지;
- setup metadata와 result metadata의 일치.

공식 시뮬레이션의 경우 가장 강한 검사는 **서버 재현 job**이다.

$$
\text{submitted setup}
\rightarrow
\text{official versioned solver}
\rightarrow
\text{recomputed result}
\rightarrow
\text{numerical comparison}.
$$

정해진 tolerance를 만족하면 `reproduced` badge를 부여할 수 있다.

## 실제 실험 contribution의 품질정보

실험 CSV 두 개가 둘 다 stress와 cycle을 갖는다고 바로 비교 가능한 것은 아니다. 따라서 실험은 추가 metadata가 중요하다.

필수 또는 강력 권장 항목은 다음과 같다.

- 시편 material 및 preparation;
- 시편 치수;
- 하중 waveform;
- 평균/진폭 stress 또는 force;
- frequency;
- load-cell calibration;
- 사용한 displacement sensor calibration;
- crack detection method;
- fatigue tester model / firmware version;
- sampling rate;
- 시험 종료 이유;
- raw data 제공 여부.

또한 raw measurement와 post-processed result를 명확히 구분한다.

## 공개 데이터 집계

나중에는 커뮤니티 전체 결과를 이용해서 다음과 같은 plot을 만들 수 있다.

- cycle에 따른 crack-initiation probability;
- frequency 비교;
- stress amplitude sweep;
- simulation-vs-experiment error distribution;
- theory version별 결과 비교.

단, 모든 aggregate point는 원본 contribution ID와 version까지 역추적 가능해야 한다.

따라서

$$
\boxed{\text{공개 aggregate database는 immutable contribution들에서 파생된다}}
$$

가 원칙이다.

## 데스크톱 앱과의 연결

데스크톱 앱에는 나중에 명확히 두 버튼을 둔다.

```text
Save locally
Publish result
```

`Publish result`는 다음을 수행한다.

1. 로컬 프로젝트/결과 schema 검증;
2. contribution manifest 및 provenance 생성;
3. 어떤 파일과 metadata가 공개되는지 사용자에게 정확히 표시;
4. 데이터 라이선스 선택/확인;
5. 웹사이트 업로드;
6. stable contribution ID/version 수신;
7. 원하면 서버 재현검증 요청.

공개하지 않고 로컬에서만 사용하는 것도 항상 가능해야 한다.

## 실제 피로시험기와의 연결

local tester agent는 완료된 피로시험을 experimental contribution으로 포장할 수 있게 한다. 그러나 실제 공개는 사용자의 명시적 동작이어야 한다.

실시간 제어와 안전은 계속 시험기 MCU/hardware가 담당한다. 공개 데이터 공유는 시험 후의 데이터/provenance 기능이지 safety-critical feedback loop의 일부가 아니다.

## 첫 구현 순서

전체 웹사이트보다 먼저 가장 작은 공개 결과 단위를 구현한다.

1. `public_result_manifest.json` schema v1 정의;
2. simulation-result / experiment-result 필수 field 정의;
3. 현재 데스크톱 앱에 `Export contribution` 추가;
4. contribution bundle local validator 구현;
5. 검증된 contribution을 받아 공개 read-only page를 보여주는 최소 API/site 구현;
6. 다운로드와 검색 추가;
7. server-side reproduction badge 추가;
8. 실제 피로시험 결과 publication 연결;
9. 결과 comparison 및 공개 aggregate plot 추가.

중앙 설계 원칙은 다음과 같다.

$$
\boxed{
\text{스크린샷이나 요약 숫자만이 아니라 입력 + raw/derived output + 정확한 provenance를 공개한다}
}
$$
