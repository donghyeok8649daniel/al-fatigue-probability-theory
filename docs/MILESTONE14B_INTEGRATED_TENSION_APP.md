# Milestone 14B — Integrated Tensile-Only FEM Application

## Status

**PREPARED SCAFFOLD / UI TOOL — NOT AN ACTIVE FATIGUE LAW**

This milestone adds a desktop application that joins the already validated C 1D bar FEM solver with the 2D/3D tensile-only visualization layer. It does **not** change the active mechanics or probability theory.

## Scope rule

The displayed specimen may look two- or three-dimensional, but the mechanical state remains strictly one-dimensional:

$$
u=u(x,t),
\qquad
\epsilon_x=\frac{\partial u}{\partial x},
\qquad
\sigma_x=E\epsilon_x.
$$

No transverse degree of freedom, Poisson contraction, shear stress, von-Mises stress, multiaxial constitutive law, or multiaxial failure criterion is introduced.

The GUI asks for specimen width $w$ and thickness $h$ because they determine the axial cross-sectional area

$$
A=wh,
$$

and because the same dimensions are useful for display extrusion. Their presence in the GUI must not be interpreted as 2D/3D mechanics.

## User workflow

The integrated app accepts:

- gauge length;
- width;
- thickness;
- Young's modulus;
- number of 1D bar elements;
- mean axial stress;
- axial stress amplitude;
- loading frequency;
- cycle count;
- temporal samples per cycle;
- display deformation scale.

Pressing `Run FEM` performs the following sequence:

1. validate the GUI inputs;
2. locate or build the C FEM executable;
3. convert mm to m and GPa to Pa explicitly;
4. calculate $A=wh$;
5. call the C solver through its command-line interface;
6. verify `nodes.csv`, `elements.csv`, and `metadata.csv` exist;
7. load the result into the same window;
8. enable time-step, 2D/3D, and stress/strain result controls.

`Save views` writes peak-tension 2D and 3D preview figures using the same strictly axial scalar field.

## Separation from the probability theory

The app currently consumes only

$$
\sigma_x(x,t),\quad \epsilon_x(x,t),\quad u(x,t).
$$

It deliberately contains no assumed $P(\lambda,t)$, no damage variable, and no crack-initiation law. Once the active 1D probability theory is physically closed, an additional post-processing channel can be attached without changing the verified C continuum solver.

The intended future interface is therefore

$$
\{\sigma_x,\epsilon_x\}
\longrightarrow
\text{1D probability module}
\longrightarrow
\text{probability-derived scalar field},
$$

but the second arrow is not implemented in this milestone.

## Validation

Unit tests cover:

- geometric unit conversion;
- $A=wh$;
- GPa-to-Pa conversion;
- rejection of invalid geometry and discretization inputs;
- exact construction of the C solver command.

CI additionally builds the C solver, runs its analytical uniform-bar self-test, invokes the integrated application in headless smoke mode, reloads the generated CSV files, and renders both 2D and 3D tensile-only previews.

---

# 마일스톤 14B — 통합 인장 전용 FEM 애플리케이션

## 상태

**준비된 스캐폴드 / UI 도구 — 활성 피로 법칙이 아님**

이 마일스톤에서는 이미 검증된 C 1D bar FEM solver와 2D/3D tensile-only 시각화 계층을 하나의 데스크톱 애플리케이션으로 연결한다. 활성 역학이나 확률이론 자체는 변경하지 않는다.

## 범위 규칙

화면의 시편은 2차원 또는 3차원처럼 보일 수 있지만 역학 상태는 끝까지 엄격한 1차원이다.

$$
u=u(x,t),
\qquad
\epsilon_x=\frac{\partial u}{\partial x},
\qquad
\sigma_x=E\epsilon_x.
$$

횡방향 자유도, 포아송 수축, 전단응력, von-Mises 응력, 다축 구성방정식, 다축 파손기준은 넣지 않는다.

GUI에서 시편 폭 $w$와 두께 $h$를 받는 이유는 축방향 단면적

$$
A=wh
$$

를 계산하기 위해서이고, 같은 치수를 화면 extrusion에도 사용하기 때문이다. GUI에 폭과 두께가 있다는 사실을 2D/3D 역학으로 해석하면 안 된다.

## 사용자 흐름

통합 앱 입력은 다음과 같다.

- 게이지 길이;
- 폭;
- 두께;
- Young's modulus;
- 1D bar element 수;
- 평균 축응력;
- 축응력 진폭;
- 하중 주파수;
- cycle 수;
- cycle당 시간 샘플 수;
- 화면상의 변형 확대배율.

`Run FEM`을 누르면 다음 순서로 실행된다.

1. GUI 입력 검증;
2. C FEM 실행파일 탐색 또는 빌드;
3. mm→m, GPa→Pa 단위변환;
4. $A=wh$ 계산;
5. C solver CLI 실행;
6. `nodes.csv`, `elements.csv`, `metadata.csv` 생성 확인;
7. 같은 창으로 결과 로드;
8. 시간 step, 2D/3D, stress/strain 결과 전환 활성화.

`Save views`는 동일한 순수 축방향 scalar field를 사용해 peak-tension 2D/3D preview를 저장한다.

## 확률이론과의 분리

현재 앱이 사용하는 역학량은

$$
\sigma_x(x,t),\quad \epsilon_x(x,t),\quad u(x,t)
$$

뿐이다. 임의의 $P(\lambda,t)$, damage variable, crack-initiation law는 넣지 않았다. 활성 1D 확률이론이 물리적으로 닫힌 뒤에는 검증된 C continuum solver를 바꾸지 않고 post-processing channel만 추가할 수 있다.

따라서 미래 인터페이스는

$$
\{\sigma_x,\epsilon_x\}
\longrightarrow
\text{1D 확률모듈}
\longrightarrow
\text{확률 유도 scalar field}
$$

형태를 목표로 하지만 두 번째 화살표는 이 마일스톤에서 구현하지 않는다.

## 검증

unit test에서는 다음을 검사한다.

- 기하 단위변환;
- $A=wh$;
- GPa→Pa 변환;
- 잘못된 기하/이산화 입력 거부;
- C solver 명령행 인자의 정확한 생성.

CI에서는 추가로 C solver를 빌드하고 uniform-bar 해석해 self-test를 실행한 뒤, 통합 앱의 headless smoke mode로 C 해석을 실제 수행하고 CSV를 다시 읽어 2D/3D tensile-only preview까지 렌더링한다.
