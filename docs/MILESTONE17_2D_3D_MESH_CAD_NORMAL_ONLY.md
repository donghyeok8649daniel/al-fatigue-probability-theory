# Milestone 17 — 1D/2D/3D Geometry Mesh and CAD Input with a 1D Normal-Only Theory

## Status

**ACTIVE GEOMETRY/MESH INTERFACE — NOT A 2D/3D ELASTICITY OR MULTIAXIAL FATIGUE SOLVER**

This milestone permits one-, two-, and three-dimensional specimen meshes without changing the active microscopic theory. The dimensionalities are deliberately separated:

$$
d_{\rm mesh}\in\{1,2,3\},
\qquad
d_{\rm probability}=1.
$$

The mesh describes specimen geometry, cell adjacency, locations for result storage, and visualization. The probability state at every retained cell remains the one-dimensional normal-spacing density $P(a,t)$.

## 1. Only one stress scalar enters the probability equation

Let $\mathbf n$ be the declared tensile-axis unit vector. If a future continuum solver supplies a stress tensor $\boldsymbol\sigma$, the exact normal projection is

$$
\boxed{
\sigma_{nn}=\mathbf n^{\mathsf T}\boldsymbol\sigma\mathbf n.
}
$$

This is an **EXACT / IDENTITY** tensor projection. The current probability solver receives only $\sigma_{nn}(\mathbf x,t)$:

$$
\frac{\partial p_{c}}{\partial\tau}
=
\partial_\lambda
\left[
(\phi'(\lambda)-\sigma_{nn,c}/E)p_c
+\chi^{-1}\partial_\lambda p_c
\right].
$$

Shear components, von-Mises stress, principal-stress combinations, and multiaxial fatigue criteria are not independent inputs to $p_c$. A multidimensional mesh therefore does not make the active probability theory multidimensional.

## 2. Current 1D-to-mesh projection

The current C solver still produces a one-dimensional axial element field $q_e(t)$, where $q_e$ can be axial stress, axial strain, or a scalar derived from $P_e$. For mesh cell $c$ with centroid $\mathbf x_c$, define

$$
s_c=\mathbf n\cdot\mathbf x_c,
\qquad
\xi_c=\frac{s_c-s_{\min}}{s_{\max}-s_{\min}}.
$$

The same normalized coordinate is constructed for the 1D bar, and the piecewise-constant value of the containing 1D element is copied to cell $c$:

$$
\boxed{q_c(t)=q_{e(\xi_c)}(t).}
$$

This is a **CONTROLLED APPROXIMATION FOR VISUALIZATION/POST-PROCESSING**. It preserves axial ordering and cross-sectional constancy. It does not solve force equilibrium, Poisson contraction, stress concentration, or shear on the 2D/3D mesh.

## 3. Actual mesh connectivity

The dependency-free generators create:

- a 1D line-element mesh for the tensile axis;
- a 2D quadrilateral mesh for a rectangular tensile region;
- a 3D hexahedral mesh for a rectangular tensile volume.

Unlike the earlier graphical extrusion, each generated cell has explicit node connectivity. Boundary-face extraction removes shared interior faces before 3D rendering, while all 1D/2D cells remain visible. Mesh cell count is a numerical discretization count and is not identified with the number of statistically independent microscopic regions.

## 4. Lightweight visibility UI

`simulations/fem_mesh_app.py` uses only NumPy and Matplotlib in the default path. It supports:

- 1D/2D/3D dimension switching;
- generated or imported CAD/mesh geometry;
- node and cell-edge visibility toggles;
- surface opacity control;
- an axial clipping slider that removes cells beyond a selected centroid position and compacts hidden nodes;
- NPZ mesh and PNG preview export.

The clipping operation changes only the displayed subset. It does not delete FEM state or modify the physical specimen.

## 5. CAD and mesh input

The geometry loader distinguishes surface geometry from a volume mesh.

| Input | Backend | Result | Limitation |
|---|---|---|---|
| STL, OBJ | built-in reader | triangle surface mesh | not a volume FEM mesh |
| STEP/STP, IGES/IGS, BREP | optional Gmsh/OpenCASCADE | generated mesh of the requested 1D, 2D, or 3D topological dimension | requires `requirements-cad.txt` |
| MSH, VTK, VTU, XDMF and other supported mesh files | optional meshio | imported top-dimensional cells | requires `requirements-cad.txt` |
| NPZ mesh bundle | built-in reader/writer | exact saved `GeometryMesh` | project interchange format |

STL and OBJ contain tessellated surfaces. Reading them does not automatically create valid tetrahedral volume cells. Requesting a 3D volume mesh from a surface-only import is rejected instead of silently treating triangles as solid elements.

CAD coordinates have an explicit conversion factor `coordinate_scale_to_m`. For example, a millimetre-based STL uses

$$
\texttt{coordinate\_scale\_to\_m}=10^{-3}.
$$

The factor must be supplied from the CAD unit convention; it is not inferred from specimen size.

## 6. Commands

Launch the lightweight mesh UI:

```bash
python -m simulations.fem_mesh_app
```

Generate a real 3D hex mesh and map the peak critical-tail diagnostic:

```bash
python -m simulations.run_tensile_mesh_projection \
  --history-csv results/data/fem_probability_demo/probability_elements.csv \
  --field critical_tail_probability \
  --step peak-tension \
  --dimension 3 \
  --output-dir results/data/fem_probability_demo/mesh_projection_3d
```

Read a millimetre-based STL as a surface mesh:

```bash
python -m simulations.run_tensile_mesh_projection \
  --history-csv results/data/fem_probability_demo/probability_elements.csv \
  --field critical_tail_probability \
  --dimension 2 \
  --geometry specimen.stl \
  --coordinate-scale-to-m 0.001 \
  --axis 1,0,0 \
  --output-dir results/data/cad_surface_projection
```

For STEP/IGES volume meshing, install the optional dependencies:

```bash
python -m pip install -r requirements-cad.txt
```

## 7. Remaining boundary

A future true 2D/3D continuum solve may supply $\boldsymbol\sigma(\mathbf x,t)$, but the current research hypothesis still consumes only $\sigma_{nn}$. Extending the microscopic state to shear spacing, slip, crystal plasticity, or a multiaxial crack criterion requires a separate physical derivation and is outside this milestone.

---

# 한국어 번역 — 1D normal-only 이론을 유지한 1D/2D/3D 메시 및 CAD 입력

## 상태

**활성 geometry/mesh 인터페이스 — 2D/3D 탄성 또는 다축 피로 solver가 아님**

이 마일스톤은 활성 microscopic theory를 바꾸지 않고 1차원·2차원·3차원 시편 메시를 허용한다. 두 차원 개념은 의도적으로 분리한다.

$$
d_{\rm mesh}\in\{1,2,3\},
\qquad
d_{\rm probability}=1.
$$

메시는 시편 형상, cell adjacency, 결과 저장 위치 및 시각화를 표현한다. 각 retained cell의 확률상태는 계속 1차원 normal-spacing density $P(a,t)$이다.

## 1. 확률방정식에는 하나의 응력 scalar만 입력된다

$\mathbf n$을 선언한 인장축 unit vector라고 하자. 향후 continuum solver가 응력 tensor $\boldsymbol\sigma$를 공급하면 정확한 normal projection은

$$
\boxed{
\sigma_{nn}=\mathbf n^{\mathsf T}\boldsymbol\sigma\mathbf n
}
$$

이다. 이는 **EXACT / IDENTITY** tensor projection이다. 현재 확률 solver에는 $\sigma_{nn}(\mathbf x,t)$만 입력한다.

$$
\frac{\partial p_{c}}{\partial\tau}
=
\partial_\lambda
\left[
(\phi'(\lambda)-\sigma_{nn,c}/E)p_c
+\chi^{-1}\partial_\lambda p_c
\right].
$$

전단성분, von-Mises 응력, 주응력 조합 및 다축 피로기준은 $p_c$의 독립 입력이 아니다. 따라서 multidimensional mesh를 사용해도 활성 확률이론이 multidimensional해지는 것은 아니다.

## 2. 현재 1D-to-mesh projection

현재 C solver는 여전히 1차원 축방향 element field $q_e(t)$를 출력한다. $q_e$는 축응력, 축변형률 또는 $P_e$에서 유도한 scalar일 수 있다. centroid가 $\mathbf x_c$인 mesh cell $c$에 대해

$$
s_c=\mathbf n\cdot\mathbf x_c,
\qquad
\xi_c=\frac{s_c-s_{\min}}{s_{\max}-s_{\min}}
$$

를 정의한다. 1D bar에도 동일한 normalized coordinate를 만들고, cell이 속하는 1D element의 piecewise-constant 값을 cell $c$로 복사한다.

$$
\boxed{q_c(t)=q_{e(\xi_c)}(t)}.
$$

이는 **시각화/post-processing을 위한 CONTROLLED APPROXIMATION**이다. 축방향 순서와 단면 내 일정성을 유지하지만 2D/3D mesh에서 force equilibrium, Poisson contraction, stress concentration 또는 shear를 풀지는 않는다.

## 3. 실제 mesh connectivity

외부 의존성이 없는 generator는 다음을 만든다.

- 인장축을 위한 1D line-element mesh;
- 직사각형 인장영역의 2D quadrilateral mesh;
- 직육면체 인장체적의 3D hexahedral mesh.

기존의 단순 graphical extrusion과 달리 각 cell에는 명시적인 node connectivity가 있다. 3D rendering 전에는 공유 interior face를 제거하고 boundary face만 추출하며, 1D/2D cell은 모두 표시한다. mesh cell 수는 numerical discretization count이며 통계적으로 독립인 microscopic region 수와 동일하지 않다.

## 4. 경량 visibility UI

`simulations/fem_mesh_app.py`의 기본 실행경로는 NumPy와 Matplotlib만 사용한다. 지원 기능은 다음과 같다.

- 1D/2D/3D 차원 전환;
- generated 또는 imported CAD/mesh geometry;
- node 및 cell-edge visibility toggle;
- surface opacity 조절;
- 선택한 centroid 위치 이후 cell을 숨기고 hidden node를 compact하는 axial clipping slider;
- NPZ mesh 및 PNG preview export.

clipping은 화면에 보이는 subset만 바꾸며 FEM state를 삭제하거나 physical specimen을 변경하지 않는다.

## 5. CAD 및 mesh 입력

geometry loader는 surface geometry와 volume mesh를 구분한다.

| 입력 | backend | 결과 | 제한사항 |
|---|---|---|---|
| STL, OBJ | 내장 reader | triangle surface mesh | volume FEM mesh가 아님 |
| STEP/STP, IGES/IGS, BREP | 선택 Gmsh/OpenCASCADE | 요청한 1D, 2D 또는 3D 위상 차원의 생성 mesh | `requirements-cad.txt` 필요 |
| MSH, VTK, VTU, XDMF 및 기타 지원 mesh | 선택 meshio | top-dimensional cell import | `requirements-cad.txt` 필요 |
| NPZ mesh bundle | 내장 reader/writer | 저장된 `GeometryMesh` 그대로 | 프로젝트 interchange format |

STL과 OBJ는 tessellated surface를 담는다. 이를 읽었다고 해서 유효한 tetrahedral volume cell이 자동 생성되는 것은 아니다. surface-only import에 3D volume mesh를 요청하면 triangle을 solid element로 몰래 취급하지 않고 오류를 낸다.

CAD 좌표에는 명시적인 `coordinate_scale_to_m` 변환계수가 있다. 예를 들어 millimetre 기반 STL에는

$$
\texttt{coordinate\_scale\_to\_m}=10^{-3}
$$

를 사용한다. 이 계수는 CAD unit convention으로부터 사용자가 공급해야 하며 시편크기에서 추측하지 않는다.

## 6. 실행 명령

경량 mesh UI 실행 명령은 다음과 같다.

```bash
python -m simulations.fem_mesh_app
```

실제 3D hex mesh를 만들고 peak critical-tail 진단값을 매핑하는 명령은 다음과 같다.

```bash
python -m simulations.run_tensile_mesh_projection \
  --history-csv results/data/fem_probability_demo/probability_elements.csv \
  --field critical_tail_probability \
  --step peak-tension \
  --dimension 3 \
  --output-dir results/data/fem_probability_demo/mesh_projection_3d
```

millimetre 기반 STL을 surface mesh로 읽는 명령은 다음과 같다.

```bash
python -m simulations.run_tensile_mesh_projection \
  --history-csv results/data/fem_probability_demo/probability_elements.csv \
  --field critical_tail_probability \
  --dimension 2 \
  --geometry specimen.stl \
  --coordinate-scale-to-m 0.001 \
  --axis 1,0,0 \
  --output-dir results/data/cad_surface_projection
```

STEP/IGES volume meshing에는 선택 의존성을 설치한다.

```bash
python -m pip install -r requirements-cad.txt
```

## 7. 남아 있는 경계

향후 실제 2D/3D continuum solve가 $\boldsymbol\sigma(\mathbf x,t)$를 공급할 수는 있지만 현재 연구가설은 계속 $\sigma_{nn}$만 사용한다. microscopic state를 shear spacing, slip, crystal plasticity 또는 multiaxial crack criterion으로 확장하려면 별도의 물리 유도가 필요하며 이 마일스톤 범위 밖이다.
