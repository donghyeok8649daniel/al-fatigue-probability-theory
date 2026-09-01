# `.ftgsim` File Association

## Purpose

Allow users to open a saved fatigue/FEM project by double-clicking a `.ftgsim` file in Windows Explorer.

## Canonical launch contract

The packaged application must continue to support the existing command-line form:

```text
AlFatigue.exe path\to\project.ftgsim
```

Windows file association should invoke this same entry point rather than creating a second project-loading path.

## Registry / installer behavior

The installer should register a per-user association by default where practical, with a stable ProgID such as:

```text
AlFatigue.Project
```

The association should define:

- extension: `.ftgsim`
- human-readable type: `Al Fatigue Simulation Project`
- application icon
- open command pointing to the installed executable with the selected file path

Do not modify Windows registry automatically during ordinary source-code execution or unit tests. Registration belongs to installer/setup logic or an explicit user action.

## File-format compatibility

`.ftgsim` remains a versioned project container. The desktop app must:

1. inspect the bundle format version before loading;
2. reject unsupported future major versions clearly;
3. migrate supported older versions explicitly rather than silently guessing;
4. never execute code embedded in a project bundle;
5. preserve checksum/integrity checks already defined by the project format.

## Future extension policy

If additional file extensions are introduced, distinguish their roles. For example:

- `.ftgsim` — editable/openable project bundle;
- future result/export extension — optional immutable/shareable result bundle;
- `.csv` — telemetry/data interchange, not an application-owned project format.

Avoid creating multiple proprietary extensions until their semantics are genuinely different.

---

# 한국어

`.ftgsim`은 Windows에서 더블클릭하면 앱으로 열리는 **프로젝트 확장자**로 유지한다.

핵심은 별도의 로딩 코드를 만들지 않고

```text
AlFatigue.exe project.ftgsim
```

이라는 동일한 실행 경로를 Windows 파일 연결에서도 사용하게 하는 것이다.

소스코드를 실행했다고 레지스트리를 자동 수정하지 않는다. 파일 연결은 installer 또는 사용자가 명시적으로 요청한 setup 과정에서만 등록한다.
