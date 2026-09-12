# LJ/Bessel 무한 격자와 확률 동역학

팀 세미나 및 연구 리마인드용 한국어 이론 자료, 40장.

- [PowerPoint](theory_core_kr_v31.pptx): 편집 가능한 텍스트·수식, 발표자 노트와 출처.
- [PDF](theory_core_kr_v31.pdf): 발표 화면용 고정 레이아웃.
- [상세 설명과 출처](theory_notes.md): 슬라이드 순서대로 정리한 수식과 유도 조건.

## 구성

| 슬라이드 | 내용 |
|---|---|
| 1–10 | 집단좌표, LJ 무한합, Gamma–Poisson–Bessel 유도, 환경 밀도와 EAM |
| 11–19 | FCC(111), ABC 적층, 활성 계면, 원자별 환경 확장, 외부 일 |
| 20–27 | 정적 접선, Smoluchowski, SG flux, 변형률, 소성·선택적 개구, 최초통과 |
| 28–33 | 유한온도 fast-a 축약, QSD, 완화·응답, 차원화, 결합 이동도 |
| 34–40 | 비균일 slip·전위, 비국소 kernel, 유한 source, 시편 확률, 식별성 |

슬라이드는 연구의 **수학적 구조**를 설명한다. 생산 TwoRowLJ 기준,
연구용 full-FCC/활성 계면/공간 결함 확장, 재료 보정 후보의 지위를 구별한다.
실제 Al 항복·피로 또는 생산 PDE의 물리 Hz가 검증되었다는 자료가 아니다.
현재 수치 연구 결과는 [v31 연구 기록](../../solver_v1/WEAK_KINETICS_AND_MATERIAL_V31.md)을 참조한다.

## 원본 및 재생성

`slides.json`이 텍스트·수식·노트의 원본이다. 수식의 `_{...}`, `^{...}`는
PowerPoint의 편집 가능한 subscript/superscript 텍스트로 변환한다.
그림으로 붙인 수식이나 별도의 Korean/English UI 구현을 만들지 않는다.

이 호스트에는 표준 artifact-tool 런타임이 없어 설치된 PowerPoint의 문서 API로
생성했다. Windows와 PowerPoint, Python, PDF 메타데이터용 PyMuPDF가 필요하다.
원본 참고 PPT나 기존 빌드 결과를 덮어쓰지 않으며 PowerShell 보안 정책을 바꾸지 않는다.

```powershell
py -3 presentations/theory_core_v31/build_native.py --out .cache/theory_deck_new
py -3 presentations/theory_core_v31/verify_native.py --build .cache/theory_deck_new
```

개별 의존성을 private 디렉터리에 설치했다면 두 명령에 `--pdf-tools`를 지정할 수 있다.
생성 후 모든 페이지의 표시와 수식·출처를 확인하고 출판 파일을 별도로 복사한다.
기존 다운로드 발표자료는 내용 참고에만 사용했으며 변경하지 않았다.

## 출판 검증

최종 40장 모두 native PNG로 시각 검토했고, PDF 40페이지도 별도 렌더러로
재생했다. 원본 JSON과 PPTX의 제목·본문·수식·노트·출처가 일치하고 텍스트
넘침은 없었다. PDF/native 화면의 차이는 글꼴 안티앨리어싱 범위로 확인했다.
개인 Office 작성자 메타데이터는 제거했다. PPTX의 수식은 이미지가 아니라
편집 가능한 텍스트이며, 자세한 유도 조건은 발표자 노트와 theory_notes.md에 있다.
