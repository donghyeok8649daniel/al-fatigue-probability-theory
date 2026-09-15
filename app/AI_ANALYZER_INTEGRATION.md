# AI analyzer / chat controller 통합 인계

## Main integration note

The independent module below has now been integrated as app/aft_ai_controller.py.
app/ai_chat_view.py wires the Tk window to app/setup_language_view.py. Local
mock integration tests exercise the actual Tk event loop, session isolation,
opt-in result summaries, stale-preview refusal and language switching after
window closure. The standalone delivery's "Tk not tested" statement below
describes its original delivery scope, not these subsequent integration tests.
Live API authentication/model access remains untested. No network inference
was performed by these tests. This work does not change physical validation.

## 완료 범위

`aft_ai_controller.py`는 Python 3.10+ 표준 라이브러리만 사용하는 독립 모듈이다.
Tk, NumPy, OpenAI SDK, 저장소의 solver/mesh/setup 모듈을 import하지 않는다.
준비 → 전송 내용 확인 → 사용자의 Send → worker HTTPS → queue → UI poll 구조다.
세션 history는 메모리에서 분리되고, 다음 요청에 포함할 과거 메시지도 명시적으로
선택해야 한다. 새 세션/새 preview/자동 poll은 API 요청을 만들지 않는다.

현재 `aft.surface-loads/1`의 고정 상태를 매 요청 instructions/metadata 및
완료 event의 로컬 안내문에 담는다:

- JSON + 실제 mesh SHA256에 묶인 **surface setup only**.
- 공간 FVM/PDE 다중응력 미연결.
- material/physical seconds/Hz 미검증. 주파수는 cycles/model_time.
- 응답은 검증 대기 텍스트/초안. correction 승인, import, 코드/도구/PDE 자동실행 없음.
- atomic physics, LJ/Bessel, mobility, A_c 변경 없음.
- 단결정 범위 고정. 무결함 단결정과 전위를 포함한 실제 단결정을 구분하며,
  다결정/입계 모델을 자동 제안하거나 추가하지 않는다. 생산 TwoRowLJ reduced
  reference와 full-FCC/dislocation 연구를 구분하고 자동 생산 승격하지 않는다.
  재료 검증 실패를 단순히 입계 누락 때문이라고 설명하지 않는다.

## UI 연결 예시

모듈을 통합 브랜치의 적절한 위치에 복사하고 import 경로만 맞춘다.
아래 코드는 연결 형태를 보이는 예시이며 실제 요청을 실행하지 않았다.

```python
from aft_ai_controller import ChatController, OpenAIResponsesTransport, Summary

# UI 메인 스레드에서 구성. key/model을 history, JSON, 로그에 넣지 않는다.
# credential_accessor는 worker-safe 메모리/credential store 접근 함수다.
# Tk StringVar.get 같은 widget 접근 함수를 넘기지 않는다.
controller = ChatController(
    OpenAIResponsesTransport(credential_accessor),
    model=explicitly_selected_model,  # 기본 모델 없음; 빈 값은 prepare에서 거부
)
sid = controller.new_session()

# 선택 summary 체크박스 기본 OFF. 사용자 검토 후에만 준비한다.
# 원본 mesh/파일/환경변수/전체 연구결과를 자동 수집하지 않는다.
preview = controller.prepare(
    sid,
    message=user_question,
    summaries=(Summary('load_summary', user_selected_non_secret_load_summary),),
    history_ids=tuple(user_selected_history_message_ids),  # 기본 ()
    confirmed_non_secret=True,  # 실제 UI 검토 결과; 자동으로 항상 True 금지
)
show_readonly_preview(preview.payload_json)  # 정확한 HTTP body, 키 없음

def on_send_button():
    # 이 함수만 Send 버튼에 연결. 분석/편집/선택/탭 이동 이벤트에 연결 금지.
    request_id = controller.send(preview, user_clicked_send=True)
    mark_busy(sid, request_id)

def poll_from_tk():
    for event in controller.poll():
        # 현재 탭에 무조건 삽입하지 말고 event.session_id로 대상 대화를 찾는다.
        # widget.insert 등은 여기, 즉 UI 메인 스레드에서만 한다.
        render_plain_text(event.session_id, event.status, event.text, event.metadata)
    root.after(80, poll_from_tk)

root.after(80, poll_from_tk)
```

예시의 UI 함수/변수는 호스트 앱이 제공한다. 모든 public controller 메서드는
생성한 UI 스레드에서만 호출해야 한다. 창 종료 시 after 예약을 취소하고
`controller.close()`를 호출한다. 외부 dependency나 widget callback을 worker에
넘기지 않는다. worker는 transport 호출과 queue 삽입만 수행한다.

### API 요약

| 메서드 | 의미 |
|---|---|
| `new_session()` | 네트워크 없는 새 대화 ID |
| `history(sid)` | 해당 대화의 불변 Message tuple, 최근 40개까지 |
| `prepare(...)` | 사용자가 선택한 내용의 로컬 snapshot, 이전 preview 무효화 |
| `send(preview, user_clicked_send=True)` | 변경되지 않은 preview를 한 번 전송 |
| `poll()` | 완료/error/configuration_error event를 UI 스레드에서 수거 |
| `cancel(sid)` | pending 결과 및 preview 무효화, history 유지 |
| `clear(sid)` | cancel + history 삭제 |
| `delete_session(sid)` | cancel + 세션 삭제 |
| `close()` | 전체 세션/preview 해제, 이후 호출 거부 |

세션당 요청 한 개, 전역 동시 worker 기본 4개다. 취소해도 실제 worker가 끝날 때까지
슬롯을 점유한다. 취소/삭제 뒤 늦은 응답은 history/UI에 반영하지 않는다.
HTTP가 이미 시작되면 cancel은 전송을 되돌리지 못하며 비용 취소도 보장하지 않는다.
30초 socket timeout을 사용하고 재시도/redirect/서버 대화 저장 기능을 사용하지 않는다.
socket timeout은 전체 작업의 절대 wall-clock deadline은 아니다.

## 셋업 언어 연결

메인 작업의 `app/setup_language.py` AFT 1 컴파일러는 이 모듈 밖에 둔다.
모델은 문자열만 반환하며 AFT 1 문법을 자동으로 안다고 가정하지 않는다.
문법 설명이 필요하면 사용자가 검토한 비밀 아닌 `setup_summary`에 선택해 넣는다.
원하는 문법 정보/mesh fingerprint/face IDs를 보내지 않았다면 모델이 발명해서는 안 된다.

사용자가 초안을 별도 편집기로 복사한 후, 호스트 앱이 명시적 검증/검토 동작에서
기존 mesh-bound 컴파일러/decoder로 검사한다. `app/surface_setup.decode`의
mesh hash, 단위, face 범위, tensor 제한 문법 및 대칭 검사를 유지한다.
문법 통과는 물리 검증이 아니다. import와 correction 재승인은 기존 사용자 흐름으로
진행한다. AI 응답을 `eval`, `exec`, shell 또는 Solve 입력으로 바로 넘기지 않는다.
이 산출물에는 import/apply/solve 메서드나 draft JSON 파서가 의도적으로 없다.

## 데이터 및 credential 경계

허용 summary 종류: `mesh_summary`, `load_summary`, `setup_summary`,
`result_summary`, `research_summary`. 값은 선택한 짧은 문자열이다.
선택하지 않은 summary나 history는 전송되지 않는다. preview 생성 이후 UI 내용을
수정하면 새 preview를 생성해야 한다. payload는 UTF-8 64 KB 이하이며 메시지/개별
summary/모델 원문 응답은 12,000자 제한이다. 원문 또는 metadata를 디스크에 저장하거나
로그로 출력하는 기능은 없다. exception 본문/HTTP 오류 body도 UI/history에 넣지 않는다.

키 전용 도구는 이 세션의 callable tool 목록에서 발견되지 않았다. 실제 key/model은
설정하지 않았고 환경의 secret을 검색하지 않았다. credential getter는 사용자가
Send한 worker 안에서만 호출된다. 키가 비어 있으면 HTTPS 연결 전에
`configuration_error`를 반환한다. 모델 미설정은 prepare 단계에서 거부한다.

일반적인 credential 패턴은 입력/출력에서 거부하고, 실제 provider key가 payload나
응답에 포함돼도 거부한다. **패턴 검사가 임의의 비밀을 완벽히 판별하지는 못한다.**
`confirmed_non_secret`와 `user_clicked_send`는 호스트 UI가 지켜야 하는 계약이며,
적대적인 Python 호출자를 격리하는 보안 sandbox가 아니다. key는 요청 처리 중 메모리와
HTTPS Authorization에만 사용되며, 호스트 앱의 debugger/외부 로깅은 별도로 관리한다.

## 공식 API 확인

2026-09-15 OpenAI Docs 스킬로 [공식 Responses API create 문서](https://developers.openai.com/api/reference/python/resources/responses/methods/create)를 검색 후 열어 확인했다.
구현은 `POST /v1/responses`, `model`, `instructions`, `input`, `store:false`,
`tools:[]`, `tool_choice:none`, `stream:false`를 사용한다. 응답의 assistant
`output_text`/`refusal`만 문자열로 읽고 tool call 등은 거부한다.
`store:false`는 응답의 후속 API 조회용 저장을 끄는 설정이며 모든 사업자 측 보존이
없다는 보장은 아니다. 실제 계정의 모델 접근권한/지원 여부는 확인하지 않았다.

## 검증 및 미연결 항목

18개 표준 unittest mock 테스트 통과. 단결정 지침과 미설정 key의 명시적 거부도 검사했다.
최종 테스트 결과는 `TEST_RESULTS.txt`를 참고한다. HTTP 테스트도 HTTPSConnection을
mock했고 **실제 네트워크 추론 호출은 0회**다. 문서 조회는 별도로 수행했다.

```powershell
Set-Location '<이 산출물 폴더>'
python -B -m unittest -v test_aft_ai_controller
```

검증 항목: 자동전송 금지, immutable/single-use preview, 선택한 history만 전송,
세션 혼합 거부, worker 분리, 취소/clear/delete 후 늦은 결과 폐기, 동시요청 제한,
secret/error 차단, 재시도/redirect 금지, 불완전 응답/tool call 거부, draft 비실행.

미검증: 실제 API 연결/계정 key/model, 실제 Tk 화면, AFT 1 컴파일러와의 end-to-end
통합, setup draft의 유효성, 공간 FVM/PDE, 물성/physical seconds/Hz.
기존 저장소의 solver/app 전체 회귀와 GUI smoke는 실행하지 않았다.
대상 저장소는 읽기 전용으로 취급했고 checkout/commit/push/병합을 수행하지 않았다.
읽기 시 HEAD는 `c8b88d98acc6b21d39e69cb45e16ed5afd1bbd67`이었다.
메인 작업이 계속 수정하므로 이 hash를 통합 시 최종 HEAD로 가정하지 않는다.
