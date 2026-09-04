# Desktop Application Layer

This directory owns the end-user desktop application and distribution layer built on top of the validated `numerical-fem` implementation.

The current entry point is:

```powershell
py -3 -m app.desktop_ui
```

`Theory Core v1` is always active. The Solve workspace exposes only the spatial discretization choice, `FVM` or `FEM`; FVM is the default. The application is organized into responsive Pre/Mesh, Solve, and Post workspaces, and Post provides normal stress, axial strain, first passage, survival, and hazard result modes.

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
