# Single-crystal loading direction in the active 1D model

## Scope

The material is pure single-crystal aluminum. Every simulation must declare a
nonzero cubic crystallographic loading direction `[h k l]`. The active fatigue
input remains one scalar normal stress along that direction.

This does not activate the archived FCC pair-sum model or a multiaxial failure
criterion. The optional active ideal-registry branch separately requires one
slip-plane normal and one in-plane slip direction and uses the signed Schmid
projection of this uniaxial load. It is not yet connected to the FEM/UI path.

## Directional scalar modulus

For direction cosines $l_1,l_2,l_3$ parallel to `[h k l]`,

$$
\frac{1}{E_{[hkl]}}=S_{11}
-2\left(S_{11}-S_{12}-\frac{S_{44}}{2}\right)
(l_1^2l_2^2+l_2^2l_3^2+l_3^2l_1^2).
$$

The cubic compliances are

$$
S_{11}=\frac{C_{11}+C_{12}}
{(C_{11}-C_{12})(C_{11}+2C_{12})},\quad
S_{12}=\frac{-C_{12}}
{(C_{11}-C_{12})(C_{11}+2C_{12})},\quad
S_{44}=\frac1{C_{44}}.
$$

$C_{44}$ appears because cubic anisotropy affects projected normal compliance.
No shear stress is supplied to the fatigue solver.

Two input modes are supported. `user_supplied_axis_modulus` uses a justified
$E_{[hkl]}$ directly. `cubic_direction_projection` requires all three of
$C_{11},C_{12},C_{44}$ and calculates $E_{[hkl]}$. The code intentionally has
no hard-coded aluminum elastic constants. Providing only some constants is
rejected.

The representative area $A_0$ remains unresolved. Direction alone does not
justify assigning $A_0$ from an FCC plane without a separate coarse-graining
definition.

Legacy examples use `[100]` and `E_axis=69 GPa` only as numerical test inputs,
not calibrated `[100]` single-crystal aluminum properties.

## 한국어 요약

순수 단결정 알루미늄이므로 모든 계산은 `[h k l]` 인장방향을 기록한다.
해당 방향의 Young 계수만 scalar normal constitutive input으로 사용한다.
cubic elastic constant가 주어지면 방향별 Young 계수를 계산한다. 선택적 활성
registry branch에서는 별도로 slip plane normal과 slip direction을 지정해
signed Schmid factor를 계산한다. 이는 full crystal plasticity나 multiaxial
fatigue criterion이 아니며 아직 FEM/UI에는 연결하지 않았다.
