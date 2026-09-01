# Inactive archived shear and ideal-slip research

> [!CAUTION]
> **INACTIVE — NO SHEAR MODEL IS USED BY THE ACTIVE SOLVER.**

This directory is deliberately isolated from the active one-dimensional
normal-tensile theory. Nothing here is imported by the default simulations,
FEM coupling, or UI. The material in this directory is documentation and
audited research only; it is not an executable shear solver.

The active model accepts one declared tensile-normal scalar,
$\sigma_{nn}=\mathbf n^T\boldsymbol\sigma\mathbf n$, and evolves a
one-dimensional normal-spacing probability state. A 2D/3D geometry mesh does
not reactivate resolved shear, slip-system kinetics, crystal plasticity,
von-Mises measures, or multiaxial crack criteria.

The two-row registry landscape is a mathematically useful ideal-slip test
geometry. It is not yet a quantitative plasticity model for single-crystal
aluminum. Activation requires, at minimum, a crystallographic FCC interface,
a validated EAM or first-principles generalized-stacking-fault surface,
calibrated kinetic coefficients, and a demonstrated source of residual slip
and hardening.

See `docs/SLIP_LATTICE_ENERGY_REVIEW.md` for the audit of the imported 23-page
source and `docs/slip_lattice_energy_corrected.tex` for the conservative
corrected note.

---

# 보관된 shear 및 이상적 slip 연구

> [!CAUTION]
> **비활성 — 현재 active solver는 shear model을 사용하지 않음.**

이 디렉터리는 활성 1차원 normal tensile 이론과 의도적으로 분리되어 있다.
여기의 자료는 문서와 검토용 연구이며 실행 가능한 shear solver가 아니다.
기본 simulation, FEM coupling 또는 UI에서 import하지 않는다. 활성 모델은
선언한 tensile-normal scalar 하나와 1D normal-spacing 확률상태만 사용한다.
2D/3D geometry mesh가 있어도 resolved shear, slip-system kinetics,
crystal plasticity 또는 multiaxial crack criterion은 활성화되지 않는다.

두 원자열 registry 에너지면은 이상적 slip을 검사하는 데 유용한 수학적
기하다. 그러나 아직 단결정 알루미늄의 정량적 소성모델은 아니다. 이를
활성화하려면 최소한 결정학적으로 정의된 FCC 계면, 검증된 EAM 또는
first-principles generalized-stacking-fault surface, 보정된 동역학 계수,
잔류 slip과 경화가 생기는 원인의 검증이 필요하다.

반입된 23쪽 원본의 검토는 `docs/SLIP_LATTICE_ENERGY_REVIEW.md`, 보수적으로
교정한 연구 노트는 `docs/slip_lattice_energy_corrected.tex`에 있다.
