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
