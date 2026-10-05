# Fatigue Tester

This directory contains the physical fatigue-tester implementation and is maintained on the dedicated `fatigue-tester` branch.

## Structure

- `docs/` — tester hardware sizing, firmware architecture, safety and integration documents.
- `firmware/` — hard real-time loading, force control, sensing, safety and telemetry core.
- `hardware/` — procurement/BOM, budget source data, mechanical/electrical hardware notes.
- `pc/` — tester-side PC telemetry/logging adapters and the boundary to the probability solver.

The intended boundary is

\[
\boxed{\text{tester + MCU} \rightarrow \text{timestamped telemetry} \rightarrow \text{PC probability solver}}
\]

The MCU does not solve `P(a,s,t)` or the four governing equations in the hard real-time loop. The probability model is imported/reused on the PC side as the theory evolves.

## Core documents

- `docs/FATIGUE_TESTER_HARDWARE.md` — hardware sizing equations, sensor chains, DCPD, safety, telemetry and procurement order.
- `docs/FIRMWARE_ARCHITECTURE.md` — real-time control architecture and HAL boundary.
- `hardware/bom/fatigue_tester_bom.csv` — editable BOM source with quantities, spare quantity, club-available quantity, purchase quantity, price and links.

## Procurement

The budget spreadsheet is generated from the BOM source. The intended purchase calculation is

\[
N_{\rm buy}=\max(N_{\rm required}+N_{\rm spare}-N_{\rm club},0).
\]

Changing the club-available quantity therefore updates the required purchase count and total budget.

---

# 한국어

이 디렉토리는 실제 피로시험기 구현물을 이론 개발과 분리하여 관리하며, 전용 `fatigue-tester` 브랜치에서 유지한다.

- `docs/`: 하드웨어 요구조건, 펌웨어 구조, 안전 및 연동 문서.
- `firmware/`: 반복 인장하중, force feedback, 센싱, 안전, telemetry.
- `hardware/`: 부품/BOM, 예산, 기계·전기 하드웨어 자료.
- `pc/`: MCU telemetry/logging 및 확률 solver 연결부.

확률 시뮬레이션 `P(a,s,t)`와 4방정식 계산은 MCU에 넣지 않고 PC에서 수행한다.


## CAD and published research handoff

- [CAD handoff guide (Korean)](docs/CAD_RESEARCH_HANDOFF_KO.md): physical
  measurement units, independent gauge channels and the MCU/PC/CAD boundary.
- [Research reference index](research/README.md): ten exact published documents
  from the probability and silicon branches, with local links and attribution.
- [Reference manifest](research/reference_manifest.json): full source commits,
  original paths, local paths, SHA-256, byte sizes and pinned GitHub URLs.
- [Collected source documents](research/published/): immutable reference
  snapshots; this does not merge the research solver or desktop applications.
- [MCU-to-CAD CSV adapter](pc/measurement_csv_adapter.py): preserve the exact
  original CSV and all sensor/fault channels while creating a CAD measurement
  CSV and provenance sidecar.

Convert a completed physical telemetry log on a filesystem that supports hard
links, such as NTFS. Choose a new output name; the adapter does not overwrite logs.

```text
python fatigue_tester/pc/measurement_csv_adapter.py raw_telemetry.csv measurements.csv --time-basis physical_seconds
```

The declaration states that `time_s` is laboratory seconds. It does not calibrate
sensors or map the probability solver's model time. Confirm that
`displacement_m` is crosshead/actuator travel. Only add
`--strain-channel-is-gauge` when `strain` is an independently measured specimen
gauge channel. The three resulting files (`measurements.csv`,
`measurements.source.csv`, `measurements.metadata.json`) belong together.

No actual CAD assembly is included in this reference collection. Linking the
repository to CAD supplies design reference material; it does not run a solver,
drive hardware, upload firmware or provide a board-complete HAL. The intended
MCU/PC boundary above is an architecture goal, not a completed solver connection.

Maximum machine force, frequency, stroke, stiffness, final actuator and sensor
calibration remain unresolved. Historical BOM capacities and prices are planning
candidates, not current ratings or current quotations. Al round grips and Si
wafer fixtures/loading methods require separate designs and validation.

## CAD 연결 상태와 제작 범위

[연결 가이드](docs/CAD_RESEARCH_HANDOFF_KO.md)와
[참고자료 목록](research/README.md)은 게시된 커밋의 자료만 모은다.
원문 링크와 해시는 [manifest](research/reference_manifest.json)에 기록한다.
개인 CAD 파일, 실측 로그, 진행 중인 연구 코드 변경은 포함하지 않는다.

실측 CSV 변환은 원본과 모든 센서·고장 채널을 보존하며, 실험의 초/Hz를
확률 PDE 시간으로 자동 환산하지 않는다. 현재 자료 연결에는 실제 CAD 조립,
솔버 실행, 보드별 완성 펌웨어나 실제 힘 제어 검증이 포함되지 않는다.
최대 하중·주파수·왕복거리와 장비 정격은 미정이고, BOM 수치는 과거 후보이다.
Al 원통 시편과 Si 웨이퍼의 고정구·하중 방식을 별도로 확인해야 한다.
