# Desktop Application Layer

This directory owns the end-user desktop application and distribution layer built on top of the validated `numerical-fem` implementation.

The current entry point is:

```powershell
py -3 -m app.desktop_ui
```

`Theory Core v1` is always active. The Solve workspace exposes only the spatial discretization choice, `FVM` or `FEM`; FVM is the default. The application is organized into responsive Pre/Mesh, Solve, and Post workspaces, and Post provides normal stress, axial strain, diameter change, first passage, survival, hazard, life-distribution, and S--N probability-quantile result modes.

The Pre workspace exposes `Theory stress scale` in MPa per model-force unit. Mean stress and stress amplitude are converted through this explicit scale before the theory solve; compressive portions are clipped because the current probability solver is opening-only. The default scale is a numerical demonstration value, not an experimentally calibrated aluminum fatigue time-scale bridge.

For a declared crystal loading direction `[h k l]`, the axial normal stress is also projected exactly onto all twelve FCC $\{111\}\langle110\rangle$ slip systems. The generated `slip_systems.csv` records signed and absolute Schmid factors and the mean, amplitude, minimum, and maximum resolved shear stress. This projection introduces no applied macroscopic shear.

The material bridge evaluates the stress-controlled Tanaka--Mura--Wu dislocation-dipole initiation relation on those exact slip systems. Its high-purity-Al baseline uses the declared lattice parameter, surface energy, elastic constants, and zero ideal lattice-friction stress; it does not fit a Basquin curve. The TMW value is used only as the physical `N50` scale. The probability shape comes directly from the existing Theory Core v1 trajectory ensemble: every trajectory's opening first-passage cycle is retained, non-events are right-censored, and no Gaussian, Weibull, or smoothing kernel is imposed. `life_distribution.csv` stores the discrete probability mass, CDF, and survival at the selected load; `sn_curve.csv` stores the `N10`, `N50`, and `N80` S--N quantile curves.

This dimensional bridge is a **conditional model assumption** supplying the irreversible dislocation mechanism and cycle scale that the dimensionless LJ solver does not derive. Mean-stress effects, a measured physical initial measure, specimen-scale aggregation, and experimental calibration remain unresolved and are not silently added. In particular, the censored tail is not extrapolated, so an expectation is not reported when the ensemble does not identify it.

The default specimen is a 50 mm long, 6 mm diameter cylinder loaded along a user-entered Cartesian tensile axis. The included `examples/meshes/default_tensile_cylinder.stl` is three-dimensional presentation geometry. Current mechanics remain uniaxial: only the tensile-axis normal stress is applied. Transverse diameter change is the declared Poisson kinematic relation $\varepsilon_\perp=-\nu\varepsilon_{nn}$, with applied transverse stress fixed at zero. Registry coordinate $s$ is not converted into macroscopic diameter without an additional scale-bridging law.

## Branch responsibility

The `desktop-app` branch owns:

- Windows executable packaging (`.exe`);
- application entry point and runtime bootstrap;
- `.ftgsim` file opening and OS file association;
- application icon, version metadata, product name, and release layout;
- installer/uninstaller design;
- portable build layout when useful;
- user-facing crash/error reporting and log location;
- migration of saved project files between supported format versions;
- release smoke tests on a clean Windows environment.

It does **not** own the governing physics or numerical implementation. Those remain in `theory-core` and `numerical-fem`.

## Dependency direction

```text
theory-core
    ↓
numerical-fem
    ↓
desktop-app
    ↓
integration
    ↓
main (reviewed release checkpoint)
```

The desktop app may consume stable APIs from `numerical-fem`. Numerical and theory branches must not import packaging code back from this branch.

## Planned application modes

The same desktop application should support three data sources behind a common interface:

1. **Simulation** — run FEM/probability calculations locally.
2. **Replay** — open a saved `.ftgsim` project or logged experimental data.
3. **Live tester** — connect to the fatigue tester through the versioned telemetry/command interface.

The live tester path must not bypass MCU-side safety. E-stop, force limits, travel limits, watchdog, and drive disable remain hardware/firmware responsibilities.

## Packaging policy

Initial Windows packaging should be local/manual. Do not add or modify GitHub Actions, repository write permissions, auto-push workflows, secrets, PATs, deploy keys, or background release automation unless explicitly requested.

## Target layout

```text
app/
├─ README.md
├─ desktop_ui.py
├─ packaging/
│  └─ windows/
│     └─ README.md
└─ file_association/
   └─ README.md
```

Windows build and installer definitions live under `app/packaging/windows/`.

---

# 한국어 요약

`desktop-app` 브랜치는 사용자가 실제로 실행하는 프로그램의 배포 계층을 담당한다.

- Windows `.exe`
- `.ftgsim` 더블클릭 실행/파일 연결
- 설치/제거
- 앱 아이콘/버전/제품명
- portable 배포본
- 사용자 로그/오류 처리
- 저장 프로젝트 버전 호환

수치해석과 물리이론 자체는 여기서 소유하지 않는다. `numerical-fem`의 검증된 기능을 가져와 앱으로 포장한다.
