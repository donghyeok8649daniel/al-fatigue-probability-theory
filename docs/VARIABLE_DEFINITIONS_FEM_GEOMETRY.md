# Variable Definitions — FEM Geometry, CAD Mesh, and Normal Projection

## Classification

- **EXACT / IDENTITY**: exact mathematical operation under the declared inputs.
- **DEFINITION**: chosen representation or coordinate.
- **CONTROLLED APPROXIMATION**: a mapping that must not be interpreted as a mechanics solution.
- **NUMERICAL INPUT**: discretization or file-import parameter.

| Symbol / code field | Definition | Meaning | Unit | Classification |
|---|---|---|---|---|
| $d_{\rm mesh}$ | topological cell dimension | 2 for triangle/quad, 3 for volume cells | dimensionless | DEFINITION |
| $d_{\rm probability}$ | dimension of spacing state | current value 1 | dimensionless | DEFINITION |
| $\mathbf x_I$ | mesh node coordinate | position of node $I$ | m | DEFINITION |
| $\mathbf x_c$ | mean of cell-node coordinates | cell centroid used for field mapping | m | DEFINITION |
| $\mathbf n$ | $\|\mathbf n\|=1$ | declared tensile-axis direction | dimensionless | DEFINITION |
| $\boldsymbol\sigma$ | continuum Cauchy stress tensor if supplied | possible future continuum output | Pa | external mechanics result |
| $\sigma_{nn}$ | $\mathbf n^{\mathsf T}\boldsymbol\sigma\mathbf n$ | only stress scalar admitted to the current probability model | Pa | EXACT / IDENTITY |
| $s_c$ | $\mathbf n\cdot\mathbf x_c$ | cell coordinate along the tensile axis | m | DEFINITION |
| $\xi_c$ | $(s_c-s_{\min})/(s_{\max}-s_{\min})$ | normalized mesh axial coordinate | dimensionless | DEFINITION |
| $q_e(t)$ | scalar on 1D element $e$ | axial stress/strain or probability-derived scalar | field-dependent | DEFINITION |
| $q_c(t)$ | $q_{e(\xi_c)}(t)$ | scalar copied to multidimensional mesh cell $c$ | field-dependent | CONTROLLED APPROXIMATION |
| `coordinate_scale_to_m` | multiplier applied to imported coordinates | explicit CAD/mesh unit conversion | m per source coordinate unit | NUMERICAL INPUT / unit metadata |
| `characteristic_length_m` | requested Gmsh target size | CAD meshing resolution; not correlation length | m | NUMERICAL INPUT |
| `nx, ny, nz` | element counts along generated axes | structured mesh resolution | counts | NUMERICAL INPUT |
| $N_{\rm cell}$ | total top-dimensional cell count | numerical storage/discretization count | count | DEFINITION |

The meshing `characteristic_length_m` is unrelated to the microscopic correlation length $\ell_c$ or statistical volume $V_c$. Likewise, $N_{\rm cell}$ is not an effective independent microscopic-unit count.

---

# 한국어 번역 — FEM geometry, CAD mesh 및 normal projection 변수

## 분류

- **EXACT / IDENTITY**: 선언한 입력 아래 정확한 수학연산.
- **DEFINITION**: 선택한 표현 또는 좌표.
- **CONTROLLED APPROXIMATION**: 역학해로 해석하면 안 되는 mapping.
- **NUMERICAL INPUT**: 이산화 또는 파일 import parameter.

| 기호 / code field | 정의 | 의미 | 단위 | 분류 |
|---|---|---|---|---|
| $d_{\rm mesh}$ | topological cell dimension | triangle/quad는 2, volume cell은 3 | 무차원 | DEFINITION |
| $d_{\rm probability}$ | spacing state의 차원 | 현재 값 1 | 무차원 | DEFINITION |
| $\mathbf x_I$ | mesh node coordinate | node $I$의 위치 | m | DEFINITION |
| $\mathbf x_c$ | cell-node coordinate 평균 | field mapping에 쓰는 cell centroid | m | DEFINITION |
| $\mathbf n$ | $\|\mathbf n\|=1$ | 선언한 인장축 방향 | 무차원 | DEFINITION |
| $\boldsymbol\sigma$ | 공급될 경우 continuum Cauchy stress tensor | 향후 가능한 continuum output | Pa | 외부 mechanics 결과 |
| $\sigma_{nn}$ | $\mathbf n^{\mathsf T}\boldsymbol\sigma\mathbf n$ | 현재 probability model에 허용되는 유일한 응력 scalar | Pa | EXACT / IDENTITY |
| $s_c$ | $\mathbf n\cdot\mathbf x_c$ | 인장축을 따른 cell 좌표 | m | DEFINITION |
| $\xi_c$ | $(s_c-s_{\min})/(s_{\max}-s_{\min})$ | normalized mesh axial coordinate | 무차원 | DEFINITION |
| $q_e(t)$ | 1D element $e$의 scalar | 축응력/축변형률 또는 확률 유도 scalar | field에 따름 | DEFINITION |
| $q_c(t)$ | $q_{e(\xi_c)}(t)$ | multidimensional mesh cell $c$에 복사한 scalar | field에 따름 | CONTROLLED APPROXIMATION |
| `coordinate_scale_to_m` | imported coordinate에 곱하는 계수 | 명시적 CAD/mesh 단위변환 | source coordinate unit당 m | NUMERICAL INPUT / unit metadata |
| `characteristic_length_m` | Gmsh target size | CAD meshing 해상도이며 correlation length가 아님 | m | NUMERICAL INPUT |
| `nx, ny, nz` | generated axis별 element 수 | structured mesh 해상도 | count | NUMERICAL INPUT |
| $N_{\rm cell}$ | top-dimensional cell 총수 | numerical storage/discretization count | count | DEFINITION |

meshing의 `characteristic_length_m`는 microscopic correlation length $\ell_c$ 또는 statistical volume $V_c$와 무관하다. 마찬가지로 $N_{\rm cell}$은 유효 독립 microscopic-unit 수가 아니다.
