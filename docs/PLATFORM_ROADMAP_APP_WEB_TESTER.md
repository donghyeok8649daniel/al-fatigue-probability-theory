# Platform Roadmap — Desktop App, Web Simulation Portal, and Fatigue-Tester Integration

## Status

**SYSTEM ARCHITECTURE ROADMAP — DOES NOT MODIFY THE ACTIVE FATIGUE THEORY**

The current repository already separates the standalone C 1D tensile FEM solver, Python orchestration, CSV outputs, and 2D/3D tensile-only presentation. This makes it possible to extend the project into a desktop product, a web simulation portal, and a real fatigue-test platform without rewriting the verified solver core.

The platform must preserve the present scope rule: visual 2D/3D geometry does not imply 2D/3D mechanics, and the active theory remains strictly one-dimensional normal tension until the probability law is physically closed and validated.

## Target architecture

The intended long-term flow is

$$
\text{simulation/test setup}
\rightarrow
\text{versioned project bundle}
\rightarrow
\begin{cases}
\text{C FEM / theory worker},\\
\text{physical fatigue tester}
\end{cases}
\rightarrow
\text{versioned result bundle}
\rightarrow
\text{desktop/web visualization and comparison}.
$$

The same setup schema should be consumable by the desktop app, the web backend, and the physical tester gateway. The same result schema should be viewable locally or uploaded to the web portal.

## Phase 1 — Package the current desktop application

### Goal

Turn the current Python application into a distributable desktop executable without changing the mechanical model.

### Recommended first implementation

- Keep the existing C solver as a separate native executable bundled with the application.
- Keep Python as the orchestration and visualization layer initially.
- Build a Windows distributable using a packaging tool such as PyInstaller.
- Later, if a more polished native UI is required, migrate only the presentation layer to PySide6/Qt while keeping the same solver and file interfaces.

### Required release metadata

Every application build should expose:

- application version;
- C solver version or Git commit;
- theory/repository Git commit;
- project-file schema version.

This prevents an old simulation result from becoming detached from the exact code that generated it.

## Phase 2 — Define a versioned simulation/test project bundle

A single portable file should contain both the setup and the reproducibility metadata. A practical implementation is a ZIP-based container with a project-specific extension such as `.aftproj`.

Example logical structure:

```text
project.aftproj
├── manifest.json
├── setup.json
├── geometry.json
├── solver/
│   ├── metadata.json
│   └── command.json
├── results/
│   ├── nodes.csv
│   ├── elements.csv
│   ├── metadata.csv
│   └── summary.json
├── figures/
│   ├── tension_2d_peak.png
│   └── tension_3d_peak.png
├── tester/
│   ├── machine_config.json
│   ├── calibration.json
│   └── raw_measurements.csv
└── logs/
    └── run.log
```

The exact contents will evolve, but the bundle must be schema-versioned from the start.

### `setup.json` should eventually describe

- specimen geometry;
- material inputs;
- loading waveform;
- mean stress;
- stress amplitude;
- frequency;
- number of cycles or stopping rule;
- discretization settings;
- selected solver/theory module;
- display settings that do not affect mechanics.

### `manifest.json` should include

- schema version;
- project UUID;
- creation timestamp;
- software version;
- solver version;
- theory/repository commit;
- result state (`setup-only`, `simulated`, `tested`, `simulation+test`);
- checksums for important generated files.

**DESIGN RULE:** visualization settings must be stored separately from physical inputs whenever possible so changing a camera angle or deformation display scale cannot silently alter the simulation definition.

## Phase 3 — Dedicated simulation website

The website should not reimplement the mechanics in JavaScript. It should submit versioned jobs to the same validated solver interface.

Recommended logical components:

```text
Web browser
   |
   | HTTPS / WebSocket
   v
API service
   |
   +---- PostgreSQL: users, projects, job metadata
   |
   +---- Object storage: .aftproj bundles, CSV, figures
   |
   +---- Job queue
              |
              v
        isolated solver worker
        (C FEM + Python theory/post-processing)
```

A practical software stack is:

- frontend: React/Next.js or equivalent;
- API: FastAPI or equivalent;
- database: PostgreSQL;
- object storage: S3-compatible storage;
- background jobs: a queue/worker system;
- solver execution: isolated container or worker process containing the native C solver and the required Python modules.

These are implementation choices, not physical-model assumptions.

### Website capabilities

The portal can provide:

- create a simulation setup in the browser;
- save it as a project bundle;
- upload an existing project/result bundle;
- launch a server-side simulation;
- show job progress;
- visualize stress/strain/time histories;
- compare multiple simulation runs;
- attach experimental results to the same project;
- share a read-only result page;
- preserve the exact software/theory version that generated each result.

## Phase 4 — Real fatigue-tester integration

The physical tester must not use the cloud/web application as the real-time control loop.

The recommended control hierarchy is

```text
Web / desktop application
        |
        | high-level test plan, start/stop request,
        | monitoring, file synchronization
        v
Local tester agent on lab PC/Raspberry Pi
        |
        | authenticated serial/CAN/Ethernet protocol
        v
MCU real-time controller
        |
        +---- actuator drive
        +---- load cell acquisition
        +---- DCPD / crack measurement
        +---- local control loop
        +---- local safety state machine
        +---- hardware interlocks / E-stop
```

### Safety boundary

**HARD REQUIREMENT:** force-control timing, actuator limits, emergency stop, travel limits, overload shutdown, and other machine-protection logic must remain local to the MCU/hardware. A network outage, browser crash, cloud outage, or website bug must not remove the tester's safety protections or destabilize the control loop.

The website may request a high-level test such as

$$
\{\sigma_m,\sigma_a,f,\text{waveform},N_{\rm cycles}\},
$$

but the MCU must validate local machine limits before arming the run.

### Local tester agent responsibilities

- discover/connect to the tester;
- transfer a validated test plan;
- read firmware and calibration versions;
- stream telemetry to the desktop/web UI;
- save raw measurements locally first;
- reconnect/resume cloud synchronization after network loss;
- package completed measurements into the same project/result bundle;
- never bypass MCU safety state or hardware interlocks.

## Simulation-to-experiment workflow

A useful final user workflow is:

1. create a specimen and loading setup once;
2. save `project.aftproj`;
3. run the numerical simulation;
4. inspect predicted tensile fields and later probability-derived outputs;
5. send the same physical loading definition to the local tester agent;
6. run the real fatigue test;
7. store load-cell, displacement, DCPD/crack, temperature, and machine-state histories in the project bundle;
8. upload/synchronize the completed bundle;
9. compare simulation and experiment using the same project identifier and loading definition.

This is preferable to maintaining unrelated simulation and experiment file formats.

## Important distinction: reproducibility versus database indexing

The portable project bundle should remain independently readable and reproducible. The web database should index projects for fast search and permissions, but it should not become the only copy of the physical setup or raw experiment data.

In other words:

$$
\text{database record} \neq \text{scientific source file}.
$$

The project bundle is the portable scientific artifact; the database is an index and service layer around it.

## Recommended implementation order

1. Freeze `setup.json` and `manifest.json` schema version 1.
2. Add project-bundle save/load to the current desktop application.
3. Package the desktop application as a Windows executable.
4. Build a small FastAPI service that accepts a project bundle and executes the same C solver.
5. Build a minimal browser frontend for setup/upload/result viewing.
6. Add user/project persistence and object storage.
7. Define the tester serial/CAN/Ethernet protocol and local agent.
8. Connect the local agent to the existing firmware HAL/controller architecture.
9. Add simulation-versus-experiment comparison.
10. Add probability/fatigue outputs only after the active 1D theory is physically validated.

## Architectural invariant

The platform should maintain the following separation:

$$
\boxed{
\text{UI} \neq \text{solver} \neq \text{theory} \neq \text{hardware control}
}
$$

They communicate through explicit, versioned interfaces and files. This allows the mechanics, probability theory, web interface, and test hardware to evolve independently without silently changing one another.

---

# 플랫폼 로드맵 — 데스크톱 앱, 웹 시뮬레이션 포털, 피로시험기 연동

## 상태

**시스템 아키텍처 로드맵 — 활성 피로 이론 자체를 변경하지 않음**

현재 저장소는 이미 독립 C 1D 인장 FEM solver, Python 실행 계층, CSV 출력, 2D/3D tensile-only 표현 계층이 분리되어 있다. 따라서 검증된 solver core를 다시 작성하지 않고도 데스크톱 제품, 웹 시뮬레이션 포털, 실제 피로시험 플랫폼으로 확장할 수 있다.

플랫폼을 확장해도 현재 범위 규칙은 유지해야 한다. 화면의 2D/3D 형상은 2D/3D 역학을 의미하지 않으며, 확률 법칙이 물리적으로 닫히고 검증되기 전까지 활성 이론은 엄격한 1D normal tension으로 유지한다.

## 목표 아키텍처

장기적으로 다음 흐름을 목표로 한다.

$$
\text{시뮬레이션/시험 설정}
\rightarrow
\text{버전이 지정된 프로젝트 파일}
\rightarrow
\begin{cases}
\text{C FEM / 이론 worker},\\
\text{실제 피로시험기}
\end{cases}
\rightarrow
\text{버전이 지정된 결과 파일}
\rightarrow
\text{데스크톱/웹 시각화 및 비교}.
$$

동일한 설정 schema를 데스크톱 앱, 웹 backend, 실제 시험기 gateway가 모두 읽을 수 있어야 한다. 동일한 결과 schema를 로컬에서 열거나 웹사이트에 업로드할 수 있어야 한다.

## 1단계 — 현재 데스크톱 앱 패키징

### 목표

역학 모델을 변경하지 않고 현재 Python 앱을 배포 가능한 데스크톱 실행파일로 만든다.

### 첫 구현 권장안

- 기존 C solver를 별도 native executable로 유지한 채 앱에 함께 포함한다.
- 초기에는 Python을 실행/시각화 계층으로 유지한다.
- PyInstaller 같은 패키징 도구를 사용해 Windows 배포본을 만든다.
- 더 완성도 높은 native UI가 필요해지면 presentation layer만 PySide6/Qt로 옮기고 solver/file interface는 그대로 유지한다.

### 필수 release metadata

모든 앱 빌드는 다음을 표시해야 한다.

- 앱 버전;
- C solver 버전 또는 Git commit;
- 이론/저장소 Git commit;
- 프로젝트 파일 schema 버전.

이를 통해 과거 시뮬레이션 결과가 어떤 코드에서 만들어졌는지 잃어버리는 문제를 막는다.

## 2단계 — 버전이 지정된 시뮬레이션/시험 프로젝트 파일

설정과 재현성 metadata를 하나의 이동 가능한 파일에 넣는 것이 좋다. 실제 구현은 `.aftproj` 같은 프로젝트 전용 확장자를 가진 ZIP 기반 container가 적합하다.

예시 논리 구조:

```text
project.aftproj
├── manifest.json
├── setup.json
├── geometry.json
├── solver/
│   ├── metadata.json
│   └── command.json
├── results/
│   ├── nodes.csv
│   ├── elements.csv
│   ├── metadata.csv
│   └── summary.json
├── figures/
│   ├── tension_2d_peak.png
│   └── tension_3d_peak.png
├── tester/
│   ├── machine_config.json
│   ├── calibration.json
│   └── raw_measurements.csv
└── logs/
    └── run.log
```

정확한 내용은 나중에 바뀔 수 있지만 schema version은 처음부터 넣어야 한다.

### `setup.json`이 장기적으로 포함할 항목

- 시편 형상;
- 재료 입력;
- 하중 파형;
- 평균 응력;
- 응력 진폭;
- 주파수;
- cycle 수 또는 종료 조건;
- 이산화 설정;
- 사용할 solver/theory module;
- 역학에 영향을 주지 않는 화면 설정.

### `manifest.json`이 포함할 항목

- schema version;
- project UUID;
- 생성 시각;
- software version;
- solver version;
- theory/repository commit;
- 결과 상태 (`setup-only`, `simulated`, `tested`, `simulation+test`);
- 주요 생성파일 checksum.

**설계 규칙:** camera angle이나 deformation display scale을 바꿨다고 시뮬레이션 정의가 몰래 바뀌지 않도록, 가능한 한 시각화 설정은 물리 입력과 별도로 저장한다.

## 3단계 — 전용 시뮬레이션 웹사이트

웹사이트가 역학식을 JavaScript로 다시 구현하면 안 된다. 이미 검증된 동일 solver interface에 버전이 지정된 job을 제출해야 한다.

권장 논리 구조:

```text
웹 브라우저
   |
   | HTTPS / WebSocket
   v
API service
   |
   +---- PostgreSQL: 사용자, 프로젝트, job metadata
   |
   +---- Object storage: .aftproj, CSV, figure
   |
   +---- Job queue
              |
              v
        격리된 solver worker
        (C FEM + Python 이론/post-processing)
```

실제 구현 스택 예시는 다음과 같다.

- frontend: React/Next.js 또는 동급 기술;
- API: FastAPI 또는 동급 기술;
- database: PostgreSQL;
- object storage: S3-compatible storage;
- background job: queue/worker 구조;
- solver 실행: native C solver와 필요한 Python module이 들어있는 격리 container 또는 worker process.

이 선택들은 구현 선택일 뿐 물리 모델 가정이 아니다.

### 웹사이트 기능

- 브라우저에서 simulation setup 생성;
- 프로젝트 파일 저장;
- 기존 프로젝트/결과 파일 업로드;
- 서버측 simulation 실행;
- job progress 표시;
- stress/strain/time history 시각화;
- 여러 simulation run 비교;
- 같은 프로젝트에 실제 실험 결과 부착;
- read-only 결과 공유 페이지;
- 각 결과를 만든 정확한 software/theory version 보존.

## 4단계 — 실제 피로시험기 연동

실제 시험기의 실시간 제어 loop를 cloud/web application에 두면 안 된다.

권장 제어 계층은 다음과 같다.

```text
웹 / 데스크톱 앱
        |
        | high-level 시험 계획, 시작/정지 요청,
        | 모니터링, 파일 동기화
        v
실험실 PC/Raspberry Pi의 local tester agent
        |
        | 인증된 serial/CAN/Ethernet protocol
        v
MCU real-time controller
        |
        +---- actuator drive
        +---- load cell acquisition
        +---- DCPD / crack measurement
        +---- local control loop
        +---- local safety state machine
        +---- hardware interlock / E-stop
```

### 안전 경계

**필수 요구사항:** force-control timing, actuator limit, emergency stop, travel limit, overload shutdown 등 machine-protection logic은 MCU/hardware에 로컬로 남아 있어야 한다. 네트워크 단절, browser crash, cloud 장애, website bug가 발생해도 시험기의 안전보호가 사라지거나 control loop가 불안정해지면 안 된다.

웹사이트는 예를 들어

$$
\{\sigma_m,\sigma_a,f,\text{waveform},N_{\rm cycles}\}
$$

형태의 high-level 시험을 요청할 수 있지만, MCU가 자체 machine limit를 검사한 뒤에만 arm해야 한다.

### Local tester agent의 역할

- 시험기 탐색/연결;
- 검증된 test plan 전송;
- firmware/calibration version 읽기;
- telemetry를 desktop/web UI로 stream;
- raw measurement를 우선 로컬에 저장;
- 네트워크가 끊겨도 시험을 유지하고 나중에 cloud sync 재개;
- 완료된 측정결과를 동일한 project/result bundle에 포함;
- MCU safety state나 hardware interlock 우회 금지.

## 시뮬레이션-실험 통합 workflow

최종적으로 유용한 사용자 흐름은 다음과 같다.

1. 시편/하중 설정을 한 번 만든다.
2. `project.aftproj`로 저장한다.
3. numerical simulation을 실행한다.
4. tensile field와 추후 probability-derived 결과를 확인한다.
5. 동일한 물리 하중 정의를 local tester agent로 보낸다.
6. 실제 피로시험을 실행한다.
7. load cell, displacement, DCPD/crack, temperature, machine-state history를 같은 project bundle에 저장한다.
8. 완성된 bundle을 업로드/동기화한다.
9. 동일한 project identifier와 loading definition 기준으로 simulation과 experiment를 비교한다.

simulation용 파일과 experiment용 파일을 따로 관리하는 것보다 이 구조가 훨씬 낫다.

## 재현성과 database indexing의 구분

이동 가능한 project bundle은 웹서비스 없이도 독립적으로 읽고 재현할 수 있어야 한다. 웹 database는 빠른 검색과 권한 관리를 위한 index 역할을 하되 물리 설정이나 raw experimental data의 유일한 원본이 되어서는 안 된다.

즉,

$$
\text{database record} \neq \text{scientific source file}.
$$

project bundle이 이동 가능한 과학적 artifact이고 database는 그 주변의 index/service layer다.

## 권장 구현 순서

1. `setup.json`, `manifest.json` schema version 1을 고정한다.
2. 현재 desktop app에 project bundle save/load를 추가한다.
3. desktop app을 Windows executable로 패키징한다.
4. project bundle을 받아 동일한 C solver를 실행하는 작은 FastAPI service를 만든다.
5. setup/upload/result view만 있는 최소 browser frontend를 만든다.
6. user/project persistence와 object storage를 추가한다.
7. tester serial/CAN/Ethernet protocol과 local agent를 정의한다.
8. local agent를 기존 firmware HAL/controller architecture에 연결한다.
9. simulation-vs-experiment 비교 기능을 추가한다.
10. probability/fatigue output은 활성 1D 이론이 물리적으로 검증된 뒤에만 추가한다.

## 아키텍처 불변조건

플랫폼 전체에서 다음 분리를 유지해야 한다.

$$
\boxed{
\text{UI} \neq \text{solver} \neq \text{theory} \neq \text{hardware control}
}
$$

각 계층은 명시적이고 versioned된 interface/file로 통신한다. 이렇게 해야 mechanics, probability theory, web interface, test hardware가 서로를 몰래 바꾸지 않고 독립적으로 발전할 수 있다.
