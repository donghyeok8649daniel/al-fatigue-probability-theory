# Fatigue Tester

This directory groups the physical fatigue-test machine implementation separately from the probability/fatigue theory.

## Structure

- `firmware/` — hard real-time loading, force control, sensing, safety and telemetry core.
- `hardware/` — hardware architecture, procurement/BOM and budget source data.
- `pc/` — reserved for tester-side PC telemetry/logging adapters and later integration with the probability solver.

The intended boundary is

\[
\boxed{\text{tester + MCU} \rightarrow \text{timestamped telemetry} \rightarrow \text{PC probability solver}}
\]

The MCU does not need to solve `P(a,s,t)` or the four governing equations in the hard real-time loop. Those calculations remain PC-side and can be imported/reused as the theory evolves.

## Procurement

The editable BOM source is:

- `hardware/bom/fatigue_tester_bom.csv`

The spreadsheet version is generated from the same source and contains formulas for club inventory, required purchase quantity and total budget.

---

# 한국어

이 디렉토리는 실제 피로시험기 구현을 이론 코드와 분리해서 묶는다.

- `firmware/`: 반복 인장하중, force feedback, 센싱, 안전, telemetry.
- `hardware/`: 시험기 하드웨어 구조, 부품/BOM, 예산 자료.
- `pc/`: MCU telemetry/logging 및 향후 확률 solver 연결부.

확률 시뮬레이션 `P(a,s,t)`와 4방정식 계산은 MCU에 넣지 않고 PC에서 수행한다.
