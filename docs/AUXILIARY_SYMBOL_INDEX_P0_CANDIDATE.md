# Auxiliary Symbol Index — candidate physical P0 construction

Status: **CANDIDATE auxiliary index.** These symbols belong to `CANDIDATE_PHYSICAL_P0_CONSTRUCTION.md` and are not promoted to the active exact finite-chain theory.

## P0 structural / thermal distinction

### $P_0^{\mathrm{str}}(a)$

- English: structural/coarse-grained initial spacing density
- 한국어: 구조적/조대화 초기 층간거리 밀도
- Mathematical definition:

$$
P_0^{\mathrm{str}}(a)=\int_{\Omega}w(x)\,\delta[a-a_0^{\mathrm{str}}(x)]\,dx.
$$

- Physical meaning: probability measure obtained by spatially pushing forward a prepared slow structural spacing field after fast phonon motion is averaged out
- Unit: 1/m in physical spacing coordinates
- Status: CANDIDATE DEFINITION
- Depends on: $a_0^{\mathrm{str}}(x)$, $w(x)$, $\Omega$

### $P_0^{\mathrm{th}}$

- English: instantaneous thermal spacing marginal
- 한국어: 순간 열적 층간거리 주변밀도
- Mathematical meaning: positional marginal of a thermal phase-space preparation
- Physical meaning: includes fast thermal/phonon displacement fluctuations
- Status: CANDIDATE DIAGNOSTIC; not used for strict P0-only initialization
- Important: generally requires a nonzero initial rate distribution and therefore does not imply $F_0=P_0\delta(c)$

### $F_0^{\mathrm{str}}(a,c)$

$$
F_0^{\mathrm{str}}(a,c)=P_0^{\mathrm{str}}(a)\delta(c).
$$

- English: structural initial phase-space density
- 한국어: 구조적 초기 위상공간 밀도
- Physical meaning: static coarse-grained preparation used by the P0-only local-traction candidate
- Status: CANDIDATE DEFINITION

### $F_0^{\mathrm{th}}(a,c)$

- English: thermal initial phase-space density
- 한국어: 열적 초기 위상공간 밀도
- Physical meaning: joint thermal distribution of spacing and rate
- Status: CANDIDATE DIAGNOSTIC
- Important: generally $F_0^{\mathrm{th}}\ne P_0^{\mathrm{th}}\delta(c)$

## Spatial structural fields

### $\Omega$

- English: representative specimen region
- 한국어: 대표 시편 영역
- Mathematical meaning: spatial domain whose prepared structural spacings are represented by one P0
- Unit: depends on spatial dimension
- Status: CANDIDATE DEFINITION

### $w(x)$

$$
\int_\Omega w(x)\,dx=1.
$$

- English: normalized spatial sampling weight
- 한국어: 정규화 공간 샘플링 가중치
- Physical meaning: weighting assigned to each location when forming the spatial push-forward P0
- Unit: reciprocal measure of $\Omega$
- Status: CANDIDATE DEFINITION

### $a_{\mathrm{ref}}$

- English: reference structural spacing
- 한국어: 기준 구조적 층간거리
- Physical meaning: reference value relative to which structural microstrain is defined
- Unit: m
- Status: CANDIDATE DEFINITION

### $a_0^{\mathrm{str}}(x)$

- English: prepared local structural spacing field
- 한국어: 준비된 국소 구조적 층간거리장
- Mathematical definition:

$$
a_0^{\mathrm{str}}(x)=a_{\mathrm{ref}}[1+\epsilon_0^{\mathrm{str}}(x)].
$$

- Unit: m
- Status: CANDIDATE DEFINITION

### $\epsilon_0^{\mathrm{str}}(x)$

$$
\epsilon_0^{\mathrm{str}}(x)
=\frac{a_0^{\mathrm{str}}(x)-a_{\mathrm{ref}}}{a_{\mathrm{ref}}}.
$$

- English: prepared structural residual microstrain
- 한국어: 준비상태 구조적 잔류 미소변형률
- Unit: dimensionless
- Status: CANDIDATE DEFINITION

### $\lambda_0(x)$

$$
\lambda_0(x)=\frac{a_0^{\mathrm{str}}(x)}{a_{\mathrm{ref}}}=1+\epsilon_0^{\mathrm{str}}(x).
$$

- English: normalized prepared structural spacing
- 한국어: 정규화 준비상태 구조적 층간거리
- Unit: dimensionless
- Status: CANDIDATE DEFINITION

### $\sigma_0^{\mathrm{res}}(x)$

- English: prepared residual normal-stress field
- 한국어: 준비상태 잔류 수직응력장
- Unit: Pa
- Status: CANDIDATE INPUT

### $q_0^{\mathrm{res}}(x)$

$$
q_0^{\mathrm{res}}(x)=\frac{\sigma_0^{\mathrm{res}}(x)}{E}.
$$

- English: normalized residual normal traction
- 한국어: 정규화 잔류 수직 트랙션
- Unit: dimensionless
- Status: CANDIDATE DEFINITION

## Thermal diagnostic symbols

### $\xi$

$$
\xi=a-a_{\mathrm{ref}}.
$$

- English: local spacing deviation
- 한국어: 국소 층간거리 편차
- Unit: m
- Status: DIAGNOSTIC DEFINITION

### $K_a$

$$
K_a=\left.\frac{d^2U}{da^2}\right|_{a_{\mathrm{ref}}}=\frac{EA_0}{a_{\mathrm{ref}}}
$$

for the retained normalization $\phi''(1)=1$.

- English: harmonic normal spacing stiffness
- 한국어: 조화근사 수직 층간거리 강성
- Unit: N/m
- Status: DERIVED DIAGNOSTIC

### $k_B$

- English: Boltzmann constant
- 한국어: 볼츠만 상수
- Unit: J/K
- Status: STANDARD CONSTANT

### $T$

- English: thermodynamic temperature
- 한국어: 열역학 온도
- Unit: K
- Status: THERMAL PREPARATION INPUT

### $\sigma_a$

$$
\sigma_a=\sqrt{\frac{k_BT}{K_a}}
$$

under the explicitly stated classical single-coordinate harmonic-canonical assumptions.

- English: harmonic thermal spacing standard deviation
- 한국어: 조화 열적 층간거리 표준편차
- Unit: m
- Status: DERIVED DIAGNOSTIC, not adopted as structural P0 width

### $\sigma_\lambda$

$$
\sigma_\lambda=\frac{\sigma_a}{a_{\mathrm{ref}}}.
$$

- English: normalized harmonic thermal spacing standard deviation
- 한국어: 정규화 조화 열적 층간거리 표준편차
- Unit: dimensionless
- Status: DERIVED DIAGNOSTIC

## Diffraction inversion symbols

### $\lambda_X$

- English: X-ray wavelength
- 한국어: X선 파장
- Unit: m
- Status: EXPERIMENTAL INPUT

### $n_B$

- English: Bragg diffraction order
- 한국어: 브래그 회절 차수
- Status: integer experimental index

### $\theta$

- English: Bragg angle
- 한국어: 브래그 각
- Unit: rad
- Status: EXPERIMENTAL COORDINATE

### $\psi$

$$
\psi=2\theta.
$$

- English: diffraction angle coordinate
- 한국어: 회절각 좌표
- Unit: rad
- Status: EXPERIMENTAL COORDINATE

### $p_\theta(\theta)$

- English: corrected strain-only Bragg-angle profile
- 한국어: 보정된 변형률 기여 브래그각 프로파일
- Mathematical meaning: angular density only after instrumental, size, overlap, mosaic/orientation and other non-strain contributions have been treated sufficiently for the adopted inversion
- Unit: 1/rad after normalization
- Status: CANDIDATE EXPERIMENTAL INPUT
