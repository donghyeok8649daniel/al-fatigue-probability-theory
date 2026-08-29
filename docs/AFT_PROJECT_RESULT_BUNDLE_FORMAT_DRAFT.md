# AFT Project / Result Bundle Format — Draft

## Status

**FORMAT DESIGN DRAFT — NOT YET FROZEN OR IMPLEMENTED**

This document defines a proposed portable file format for the current tensile/FEM application, the future probability module, the fatigue tester, and the open-results platform. It does not modify the active 1D fatigue theory.

The main design objective is to provide dedicated application file extensions without creating a proprietary opaque binary format.

## Proposed extensions

Two user-facing extensions are recommended:

- `.aftproj` — mutable local project/workspace;
- `.aftres` — result contribution intended for sharing, publication, citation, or immutable archival.

Both extensions use the same open container principles and share common schemas. Their semantic roles differ.

### `.aftproj`

A project may contain setup only, simulation results, tester configuration, experimental data, display preferences, and unpublished working notes. It is expected to change over time.

### `.aftres`

A result bundle is a versioned scientific artifact. Once a specific public version is published, its primary files and checksums should be immutable. Corrections create a new version rather than silently replacing the old one.

The website may index an `.aftres` bundle in a database, but the database record is not the scientific source artifact.

## Why two extensions instead of one

A single extension would be simpler, but it would blur the distinction between an editable workspace and a citable result. The two-extension design gives the operating system and user interface an immediate semantic distinction while still allowing one shared parser/validator implementation.

The application should therefore implement one common bundle library with a `bundle_kind` discriminator rather than maintaining two unrelated file formats.

## Container format

**DESIGN PROPOSAL:** use a standard ZIP/ZIP64 container with UTF-8 path names and ordinary JSON/CSV/PNG files.

Advantages:

- no custom binary parser is required;
- files remain inspectable with standard archive tools;
- JSON and CSV remain usable outside the application;
- compression is effective for text-heavy time histories;
- future desktop, web, and tester software can share the same implementation;
- long-term scientific accessibility is better than an opaque serialized Python object or database dump.

The application must never depend on Python `pickle` for a portable/public bundle.

## Root manifest

Every bundle should contain the following root file:

```text
aft-manifest.json
```

Minimal proposed fields:

```json
{
  "format": "aft-bundle",
  "schema_version": "0.1.0",
  "bundle_kind": "project",
  "bundle_id": "UUID",
  "created_utc": "ISO-8601 timestamp",
  "generator": {
    "application_version": "...",
    "repository_commit": "...",
    "solver_commit": "...",
    "theory_commit": "..."
  },
  "paths": {
    "setup": "setup.json",
    "geometry": "geometry.json",
    "provenance": "provenance.json"
  },
  "checksums": {}
}
```

For a public result bundle:

```json
"bundle_kind": "result"
```

The manifest is a machine-readable routing/index file. Scientific inputs must remain in explicit domain files rather than being hidden in application state.

## Proposed `.aftproj` structure

```text
example.aftproj
├── aft-manifest.json
├── setup.json
├── geometry.json
├── provenance.json
├── display.json
├── solver/
│   ├── metadata.json
│   └── command.json
├── results/
│   ├── nodes.csv
│   ├── elements.csv
│   ├── metadata.csv
│   └── summary.json
├── probability/
│   └── ...                    # only after the 1D P theory is validated
├── experiment/
│   ├── raw_measurements.csv
│   ├── calibration.json
│   └── machine_metadata.json
├── figures/
└── logs/
```

Not every directory is required. A setup-only project may contain only the manifest, setup, geometry, and provenance files.

## Proposed `.aftres` structure

```text
example.aftres
├── aft-manifest.json
├── setup.json
├── geometry.json
├── provenance.json
├── results/
│   ├── summary.json
│   ├── nodes.csv
│   ├── elements.csv
│   └── optional probability outputs
├── experiment/
│   ├── raw_measurements.csv
│   ├── calibration.json
│   └── machine_metadata.json
├── figures/
├── README.md
└── LICENSE-DATA
```

A result can be simulation-only, experiment-only, or simulation-plus-experiment.

## Separation of physical inputs and visualization state

Physical setup and display state must be stored separately.

For example:

```text
setup.json      -> physical loading/material/discretization inputs
display.json    -> camera, 2D/3D selection, deformation scale, selected field
```

Changing `display.json` must never change the scientific simulation definition.

## Versioning rule

Use semantic schema versioning:

```text
MAJOR.MINOR.PATCH
```

Proposed interpretation:

- MAJOR — incompatible schema change;
- MINOR — backward-compatible field/file additions;
- PATCH — clarifications or compatible corrections.

A reader should reject an unsupported newer MAJOR version. Unknown optional fields from a newer MINOR version may be preserved/ignored only when the schema explicitly permits that behavior.

## Reproducibility and provenance

The bundle should record, where applicable:

- exact application version;
- C solver commit;
- probability-theory commit;
- repository commit;
- compiler and relevant Python package versions;
- operating system / CPU architecture where relevant;
- simulation command/configuration;
- tester firmware version;
- calibration identifiers;
- timestamps;
- SHA-256 checksum for primary scientific files.

The file format must distinguish generated results from user-entered interpretation/notes.

## Public-result immutability

For a published `.aftres` version:

1. calculate checksums for the primary files;
2. store them in the manifest/provenance record;
3. assign a stable contribution ID and version;
4. never overwrite that published version;
5. corrections produce a new version with a link to the preceding version.

A future optional digital-signature layer may sign the manifest/checksum set, but signatures are not necessary for the first implementation.

## Security rules

Especially for publicly uploaded `.aftres` files, the validator should:

- reject path traversal such as `../`;
- reject absolute archive paths;
- reject symbolic links;
- enforce decompressed-size and file-count limits to mitigate ZIP bombs;
- reject executable content by default;
- validate JSON against schemas;
- validate declared units/columns;
- verify checksums;
- never execute scripts embedded in a contribution bundle.

External solver results may be represented through provenance metadata without embedding executable solver binaries.

## Large experimental data

The first format should favor simplicity and openness. CSV inside ZIP is acceptable for v1 and compresses well for many repetitive numeric histories.

If large community datasets later make CSV inefficient, a backward-compatible MINOR revision may allow an additional columnar format such as Parquet while retaining a documented open schema. The initial application should not block on that optimization.

## Desktop file association

When the desktop application is packaged, the installer can register:

- `.aftproj` -> AFT application project;
- `.aftres` -> AFT result bundle.

Double-click should open both in the same application but in different modes:

- `.aftproj`: editable project mode;
- `.aftres`: result/read-only mode by default, with an explicit `Create derived project` action.

`Create derived project` should create a new `.aftproj` while preserving parent result ID/version provenance.

## Website behavior

The open-results website should accept `.aftres` directly. It may accept `.aftproj` only through an explicit export/publish conversion step.

Recommended publication path:

```text
.aftproj
   -> local validation
   -> Export / Publish Result
   -> .aftres
   -> server validation
   -> optional reproduction job
   -> immutable public contribution version
```

This prevents unpublished workspace state, private notes, temporary logs, or unrelated local files from being accidentally published.

## Fatigue-tester behavior

The tester/local agent may open the relevant physical subset from `.aftproj` as a high-level test plan, but the MCU remains responsible for local limit validation and safety.

After a test, the local agent can append experimental files to the project and export a separate `.aftres` contribution. Publication remains an explicit user action.

## Recommended implementation sequence

1. Keep `.aftproj` and `.aftres` names as a draft until the first working loader exists.
2. Freeze `aft-manifest.json` schema v1 only after one real save/load round-trip is tested.
3. Implement Python `aft_bundle` save/load/validate functions with no GUI dependency.
4. Add round-trip tests and malicious-ZIP/path tests.
5. Add `Save Project` / `Open Project` to the current tensile app.
6. Add `Export Result` to create `.aftres`.
7. Register file associations in the packaged Windows application.
8. Reuse the same validator on the future website.
9. Reuse the same schema in the tester local agent.

## Decision summary

The recommended architecture is:

```text
                  shared open bundle library
                         /        \
                        /          \
                  .aftproj        .aftres
                 mutable work    immutable result
```

The dedicated extension is therefore a presentation/association layer over a documented open scientific container, not a proprietary data lock-in mechanism.

---

# AFT 프로젝트 / 결과 번들 포맷 — 초안

## 상태

**포맷 설계 초안 — 아직 동결되거나 구현된 규격이 아님**

이 문서는 현재 인장/FEM 애플리케이션, 향후 확률모듈, 피로시험기, 공개 결과 플랫폼에서 함께 사용할 이동 가능한 파일 포맷을 제안한다. 활성 1D 피로이론 자체는 변경하지 않는다.

핵심 목표는 전용 확장자를 가지면서도 독점적이고 불투명한 바이너리 포맷을 만들지 않는 것이다.

## 제안 확장자

사용자에게 보이는 확장자는 두 개를 권장한다.

- `.aftproj` — 계속 수정하는 로컬 프로젝트/작업공간;
- `.aftres` — 공유·공개·인용·불변 보관을 위한 결과 contribution.

둘은 동일한 공개 container 원칙과 공통 schema를 사용하지만 의미가 다르다.

### `.aftproj`

설정만 있을 수도 있고, 시뮬레이션 결과, 시험기 설정, 실험 데이터, 화면 설정, 공개하지 않을 작업 메모 등이 함께 들어갈 수 있다. 계속 바뀌는 파일이다.

### `.aftres`

버전이 있는 과학적 결과물이다. 특정 공개 version이 게시된 뒤에는 주요 파일과 checksum을 수정하지 않는다. 수정이 필요하면 이전 결과를 덮어쓰지 않고 새 version을 만든다.

웹사이트 database는 `.aftres`를 검색/표시하기 위해 index할 수 있지만 database record 자체를 과학적 원본으로 보지 않는다.

## 확장자를 두 개로 나누는 이유

하나만 쓰면 구현은 단순하지만 수정 중인 작업공간과 인용 가능한 결과의 의미가 섞인다. 확장자를 두 개 사용하면 OS와 UI에서도 둘을 즉시 구별할 수 있다. 그러나 내부 parser/validator는 공통으로 사용한다.

따라서 서로 다른 두 파일포맷을 따로 구현하지 말고 공통 bundle library에서 `bundle_kind`만 구분하는 구조가 적합하다.

## 컨테이너 형식

**설계 제안:** UTF-8 path와 일반 JSON/CSV/PNG 파일을 담는 표준 ZIP/ZIP64 container를 사용한다.

장점:

- 별도 custom binary parser가 필요 없음;
- 일반 압축도구로 내용 확인 가능;
- JSON/CSV를 앱 밖에서도 사용할 수 있음;
- text 중심 time history 압축효율이 좋음;
- 데스크톱/웹/시험기 프로그램이 같은 구현을 공유 가능;
- opaque serialized object/database dump보다 장기 연구자료 보존에 유리함.

이동/공개 가능한 bundle에는 Python `pickle`을 사용하지 않는다.

## 루트 manifest

모든 bundle의 root에는 다음 파일을 둔다.

```text
aft-manifest.json
```

최소 필드 초안:

```json
{
  "format": "aft-bundle",
  "schema_version": "0.1.0",
  "bundle_kind": "project",
  "bundle_id": "UUID",
  "created_utc": "ISO-8601 timestamp",
  "generator": {
    "application_version": "...",
    "repository_commit": "...",
    "solver_commit": "...",
    "theory_commit": "..."
  },
  "paths": {
    "setup": "setup.json",
    "geometry": "geometry.json",
    "provenance": "provenance.json"
  },
  "checksums": {}
}
```

공개 결과에서는 다음을 사용한다.

```json
"bundle_kind": "result"
```

manifest는 파일 위치와 형식을 알려주는 machine-readable index이고, 실제 과학 입력값은 명시적인 domain file에 둔다.

## `.aftproj` 구조 초안

```text
example.aftproj
├── aft-manifest.json
├── setup.json
├── geometry.json
├── provenance.json
├── display.json
├── solver/
│   ├── metadata.json
│   └── command.json
├── results/
│   ├── nodes.csv
│   ├── elements.csv
│   ├── metadata.csv
│   └── summary.json
├── probability/
│   └── ...                    # 1D P 이론 검증 후에만
├── experiment/
│   ├── raw_measurements.csv
│   ├── calibration.json
│   └── machine_metadata.json
├── figures/
└── logs/
```

모든 폴더가 필수인 것은 아니다. 설정만 있는 프로젝트라면 manifest/setup/geometry/provenance만 있어도 된다.

## `.aftres` 구조 초안

```text
example.aftres
├── aft-manifest.json
├── setup.json
├── geometry.json
├── provenance.json
├── results/
│   ├── summary.json
│   ├── nodes.csv
│   ├── elements.csv
│   └── optional probability outputs
├── experiment/
│   ├── raw_measurements.csv
│   ├── calibration.json
│   └── machine_metadata.json
├── figures/
├── README.md
└── LICENSE-DATA
```

simulation-only, experiment-only, simulation+experiment 모두 가능하다.

## 물리 입력과 화면 설정의 분리

물리 setup과 display state는 분리해 저장해야 한다.

```text
setup.json      -> 하중/재료/이산화 등 물리 입력
display.json    -> camera, 2D/3D 선택, deformation scale, 표시 field
```

`display.json`을 바꿔도 과학적 simulation 정의는 바뀌면 안 된다.

## 버전 규칙

schema는 semantic versioning을 사용한다.

```text
MAJOR.MINOR.PATCH
```

- MAJOR — 호환되지 않는 schema 변경;
- MINOR — 하위호환되는 field/file 추가;
- PATCH — 호환되는 수정/명확화.

reader는 지원하지 않는 미래 MAJOR version을 거부한다. 새로운 MINOR의 unknown optional field를 무시/보존하는 동작은 schema가 허용할 때만 사용한다.

## 재현성과 provenance

가능한 경우 다음을 기록한다.

- 앱 버전;
- C solver commit;
- 확률이론 commit;
- 저장소 commit;
- compiler 및 관련 Python package version;
- 필요한 경우 OS/CPU architecture;
- 실제 solver command/configuration;
- 시험기 firmware version;
- calibration identifier;
- timestamp;
- 주요 과학 데이터 파일의 SHA-256 checksum.

생성된 결과와 사용자가 작성한 해석/메모도 구분한다.

## 공개 결과의 불변성

게시되는 `.aftres` version에 대해서는:

1. primary file checksum 계산;
2. manifest/provenance에 checksum 저장;
3. 안정된 contribution ID/version 부여;
4. 이미 공개된 version 덮어쓰기 금지;
5. 수정본은 이전 version 링크를 가진 새 version으로 생성.

나중에 manifest/checksum set에 digital signature를 추가할 수 있지만 초기 구현에 필수는 아니다.

## 보안 규칙

특히 공개 업로드 `.aftres` validator는 다음을 검사한다.

- `../` 같은 path traversal 거부;
- absolute archive path 거부;
- symbolic link 거부;
- ZIP bomb 방지를 위한 압축해제 크기/file count 제한;
- 기본적으로 실행파일 거부;
- JSON schema 검증;
- 단위/column 정의 검증;
- checksum 검증;
- bundle 내부 script를 절대로 실행하지 않음.

외부 solver 결과는 executable을 묶지 않고 provenance metadata로 표현할 수 있다.

## 대용량 실험 데이터

첫 version은 단순성과 개방성을 우선한다. ZIP 내부 CSV는 반복 numeric history에서 압축도 잘 되고 v1에 적합하다.

커뮤니티 데이터가 커져 CSV 효율이 문제가 되면 하위호환 MINOR version에서 Parquet 같은 columnar format을 추가로 허용할 수 있다. 초기 앱은 이 최적화 때문에 지연시키지 않는다.

## 데스크톱 파일 연결

앱 패키징 시 installer가 다음을 등록할 수 있다.

- `.aftproj` -> AFT application project;
- `.aftres` -> AFT result bundle.

둘 다 같은 앱으로 열되 모드를 다르게 한다.

- `.aftproj`: 편집 가능한 project mode;
- `.aftres`: 기본 read-only result mode + 명시적인 `Create derived project` 기능.

`Create derived project`는 parent result ID/version provenance를 유지한 새 `.aftproj`를 만든다.

## 웹사이트 동작

공개 결과 사이트는 `.aftres`를 직접 받는다. `.aftproj`는 명시적인 export/publish 변환을 거쳐야 공개할 수 있게 하는 것이 좋다.

```text
.aftproj
   -> local validation
   -> Export / Publish Result
   -> .aftres
   -> server validation
   -> optional reproduction job
   -> immutable public contribution version
```

이렇게 하면 작업 중 개인 메모/임시 log/관련 없는 local file이 실수로 공개되는 위험을 줄일 수 있다.

## 피로시험기 동작

tester/local agent는 `.aftproj`의 필요한 물리 setup subset을 high-level test plan으로 읽을 수 있다. 하지만 실제 MCU는 local machine limit과 safety를 독립적으로 검사해야 한다.

시험 종료 후 local agent가 experiment file을 project에 붙이고 별도 `.aftres` contribution을 export할 수 있다. 공개 여부는 항상 사용자가 명시적으로 선택한다.

## 구현 권장 순서

1. 실제 loader가 하나 나오기 전까지 `.aftproj`/`.aftres` 이름은 draft로 유지.
2. 실제 save/load round-trip을 한 번 검증한 뒤에만 `aft-manifest.json` schema v1 동결.
3. GUI와 독립적인 Python `aft_bundle` save/load/validate 함수 구현.
4. round-trip test와 악성 ZIP/path test 추가.
5. 현재 tensile app에 `Save Project` / `Open Project` 추가.
6. `.aftres`를 만드는 `Export Result` 추가.
7. Windows 패키지에서 file association 등록.
8. 미래 웹사이트에서도 같은 validator 재사용.
9. 시험기 local agent에서도 같은 schema 재사용.

## 결정 요약

권장 구조는 다음과 같다.

```text
                  shared open bundle library
                         /        \
                        /          \
                  .aftproj        .aftres
                 mutable work    immutable result
```

전용 확장자는 독점 데이터 lock-in이 아니라, 문서화된 공개 scientific container 위의 사용자 친화적인 file association layer로 사용한다.
