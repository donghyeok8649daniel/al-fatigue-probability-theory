# Windows Packaging Plan

## Objective

Produce a reproducible Windows desktop build of the Pre/Solve/Post application without changing the physics/numerics. Theory Core v1 is always active; the user selects FVM or FEM only for the spatial reference field.

## Initial deliverables

1. `AlFatigue.exe` (working product name; can be renamed later)
2. portable ZIP layout
3. installer package
4. uninstall entry
5. application icon and version metadata
6. `.ftgsim` file association
7. clean-machine launch test

## Candidate build path

For the current Python/Matplotlib desktop stack, the first practical packaging route is PyInstaller or an equivalent frozen-Python packager. The exact tool is not frozen until the current imports, native dependencies, data files, and startup behavior have been audited.

Do not treat a successful local `.exe` build as release-ready. Verify at minimum:

- launch on a clean Windows machine without a development Python installation;
- opening a `.ftgsim` file by command-line path;
- opening a `.ftgsim` file by Explorer double-click after association;
- FEM solver subprocess/resource lookup;
- writable user-data/log directory;
- no assumptions about repository-relative paths;
- graceful error if required optional CAD dependencies are missing;
- application version displayed in UI/logs;
- compatibility with project bundle format versioning.

## Local build and installer

Run the packaging script from the repository root:

```powershell
.\app\packaging\windows\build.ps1
```

If `fem1d\bin\fem1d_solver.exe` exists, it is included as the native FEM backend. Otherwise the packaged application uses the equivalent bundled Python two-node linear-bar FEM implementation. FVM remains the default spatial backend in both cases.

The build also bundles `examples\meshes\default_tensile_cylinder.stl`, so the default 3D tensile specimen is available without a repository checkout.

With Inno Setup installed, compile `app/packaging/windows/AlFatigue.iss` to
produce `dist/installer/AlFatigue-Setup-0.1.0.exe`. The installer registers
the per-user `.ftgsim` association and creates Start Menu/Desktop shortcuts.

## Installer responsibilities

The installer should own only machine/user integration tasks:

- install files under a stable application directory;
- create Start Menu shortcut if desired;
- register `.ftgsim` association;
- register icon and application display name;
- add uninstall metadata;
- preserve user project files on uninstall unless explicitly selected otherwise.

## Security / repository policy

No automatic signing, secret storage, PAT, deployment key, GitHub Actions write permission, or automatic release upload is added by this branch without explicit authorization.

Code signing can be added later as a separate release-hardening task.

---

# 한국어

이 디렉토리는 Windows 실행파일과 설치 패키지를 담당한다.

초기 목표는 `AlFatigue.exe`, portable ZIP, installer, uninstall, 아이콘/버전정보, `.ftgsim` 파일 연결이다.

현재 Python/Matplotlib UI 구조에서는 우선 PyInstaller 계열이 현실적이지만, 실제 dependency audit 후 고정한다. 개발 PC에서 exe가 켜지는 것만으로 완료 처리하지 않고, Python이 없는 깨끗한 Windows 환경과 `.ftgsim` 더블클릭 실행까지 검증한다.
