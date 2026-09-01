# Ideal-registry source audit and reference implementation (historical kernel)

> The active model now treats this two-row result as the kernel $W(d,s)$ and
> sums it without multiplicity as $U_0(a,s)=\sum_{k\ge1}W(ka,s)$. This folder
> remains an independent single-row audit; separate-branch language below is
> pre-multilayer history.

This directory preserves the source audit and independent direct-sum reference
checks. The optional active implementation lives under
`theory/registry_lattice.py` and `theory/registry_plasticity.py`; the primary
normal-tensile workflow, FEM coupling, and UI do not import this audit library.
The optional solver uses one declared scalar resolved registry force; it does
not activate von-Mises measures, multiaxial crack criteria, or conventional
crystal plasticity. A 2D/3D geometry mesh does not activate it automatically.

The two-row registry landscape is an active ideal single-registry mechanism,
but it is not a quantitative plasticity model for single-crystal aluminum.
Quantitative use requires a crystallographic FCC interface, a validated EAM or
first-principles generalized-stacking-fault surface, calibrated kinetic
coefficients, a homogenization thickness, and dislocation hardening physics.

See `docs/SLIP_LATTICE_ENERGY_REVIEW.md` for the audit of the imported 23-page
source and `docs/slip_lattice_energy_corrected.tex` for the conservative
corrected note.

---

# 이상적 registry 원본 검토와 기준 구현

이 디렉터리는 원본 검토와 독립 direct-sum 기준 계산을 보존한다. 선택적 활성
구현은 `theory/registry_lattice.py`와 `theory/registry_plasticity.py`에 있으며,
주 normal-tensile workflow, FEM coupling 및 UI는 이 검토용 library를
import하지 않는다. 선택 solver는 하나의 scalar resolved registry force만
사용하며 von-Mises, multiaxial crack criterion 또는 통상적인 crystal
plasticity를 활성화하지 않는다. 2D/3D mesh가 이를 자동 활성화하지 않는다.

두 원자열 registry 에너지면은 활성 ideal single-registry mechanism이다.
그러나 아직 단결정 알루미늄의 정량적 소성모델은 아니다. 정량화하려면
결정학적으로 정의된 FCC 계면, 검증된 EAM 또는 first-principles GSF surface,
보정된 동역학 계수, homogenization thickness 및 전위 경화 물리가 필요하다.

반입된 23쪽 원본의 검토는 `docs/SLIP_LATTICE_ENERGY_REVIEW.md`, 보수적으로
교정한 연구 노트는 `docs/slip_lattice_energy_corrected.tex`에 있다.
