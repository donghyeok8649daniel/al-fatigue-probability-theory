import json
import threading
import time
import unittest
from unittest.mock import patch

from app.aft_ai_controller import (BOUNDARY, BoundaryError, ChatController,
                               OpenAIResponsesTransport, Preview, Summary)


class MockTransport:
    def __init__(self, answer="검증 대기 초안", gate=None):
        self.calls = []
        self.answer = answer
        self.gate = gate
        self.entered = threading.Event()

    def complete(self, payload):
        self.calls.append((threading.get_ident(), json.loads(payload)))
        self.entered.set()
        if self.gate is not None:
            if not self.gate.wait(2):
                raise RuntimeError("test gate timeout")
        if isinstance(self.answer, Exception):
            raise self.answer
        return self.answer


def finish(c):
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline:
        events = c.poll()
        if events:
            return events
        time.sleep(.002)
    raise AssertionError("completion timeout")


class ControllerTests(unittest.TestCase):
    def setUp(self):
        self.t = MockTransport()
        self.c = ChatController(self.t, "user-selected-model")
        self.sid = self.c.new_session()

    def preview(self, sid=None, **kwargs):
        return self.c.prepare(sid or self.sid, "선택 내용을 설명", confirmed_non_secret=True, **kwargs)

    def test_no_request_on_start_prepare_poll_or_cancel(self):
        p = self.preview()
        self.c.poll()
        with self.assertRaises(BoundaryError):
            self.c.send(p)
        self.c.cancel(self.sid)
        self.assertEqual(self.t.calls, [])

    def test_non_secret_review_required_and_missing_model(self):
        with self.assertRaises(BoundaryError):
            self.c.prepare(self.sid, "hi")
        c = ChatController(self.t)
        with self.assertRaises(BoundaryError):
            c.prepare(c.new_session(), "hi", confirmed_non_secret=True)
        self.assertEqual(self.t.calls, [])

    def test_exact_snapshot_single_use_worker_and_metadata(self):
        p = self.preview(summaries=(Summary("load_summary", "2 loads, model time"),))
        self.c.send(p, user_clicked_send=True)
        event, = finish(self.c)
        self.assertEqual(event.status, "completed")
        self.assertEqual(event.metadata, BOUNDARY)
        self.assertIn("미검증", event.text)
        tid, payload = self.t.calls[0]
        self.assertNotEqual(tid, threading.get_ident())
        self.assertEqual(payload, json.loads(p.payload_json))
        self.assertFalse(payload["store"])
        self.assertEqual(payload["tools"], [])
        self.assertEqual(payload["tool_choice"], "none")
        self.assertNotIn("previous_response_id", payload)
        self.assertNotIn("conversation", payload)
        for keyword in ("aft.surface-loads/1", "SHA256", "FVM", "UNVALIDATED", "A_c",
                        "SINGLE CRYSTAL", "defect-free", "dislocation-bearing", "TwoRowLJ",
                        "missing grain boundaries"):
            self.assertIn(keyword, payload["instructions"])
        with self.assertRaises(BoundaryError):
            self.c.send(p, user_clicked_send=True)

    def test_sessions_and_explicit_history_selection(self):
        self.c.send(self.preview(), user_clicked_send=True)
        finish(self.c)
        history = self.c.history(self.sid)
        self.assertEqual(len(history), 2)
        sid2 = self.c.new_session()
        self.assertEqual(self.c.history(sid2), ())
        with self.assertRaises(BoundaryError):
            self.preview(sid2, history_ids=(history[0].id,))
        p = self.preview()
        self.assertEqual(len(json.loads(p.payload_json)["input"]), 1)
        p = self.preview(history_ids=tuple(m.id for m in history))
        self.assertEqual(len(json.loads(p.payload_json)["input"]), 3)

    def test_mutated_superseded_and_cleared_previews_rejected(self):
        p = self.preview()
        fake = Preview(p.token, p.session_id, p.payload_json + " ")
        with self.assertRaises(BoundaryError):
            self.c.send(fake, user_clicked_send=True)
        self.preview()
        with self.assertRaises(BoundaryError):
            self.c.send(p, user_clicked_send=True)
        p = self.preview()
        self.c.clear(self.sid)
        with self.assertRaises(BoundaryError):
            self.c.send(p, user_clicked_send=True)

    def test_secret_and_unknown_summary_rejected_without_network(self):
        for value in ("sk-abcdefghijk", "Authorization: Bearer abc", "api_key=abc", "password: x"):
            with self.assertRaises(BoundaryError):
                self.c.prepare(self.sid, value, confirmed_non_secret=True)
        with self.assertRaises(BoundaryError):
            self.preview(summaries=(Summary("credentials", "hidden"),))
        with self.assertRaises(BoundaryError):
            self.preview(summaries=(Summary("result_summary", "secret: abc"),))
        self.assertEqual(self.t.calls, [])

    def test_failure_not_logged_or_recorded_and_no_retry(self):
        self.t.answer = RuntimeError("sk-sensitivecredential private server body")
        self.c.send(self.preview(), user_clicked_send=True)
        event, = finish(self.c)
        self.assertEqual(event.status, "error")
        self.assertNotIn("sensitive", event.text + repr(event))
        self.assertEqual(self.c.history(self.sid), ())
        self.assertEqual(len(self.t.calls), 1)

    def test_output_secret_blocked(self):
        self.t.answer = "sk-abcdefghijk"
        self.c.send(self.preview(), user_clicked_send=True)
        self.assertEqual(finish(self.c)[0].status, "error")

    @patch("app.aft_ai_controller.http.client.HTTPSConnection")
    def test_missing_key_configuration_event_without_connection(self, connection):
        c = ChatController(OpenAIResponsesTransport(lambda: ""), "mock")
        sid = c.new_session()
        c.send(c.prepare(sid, "hi", confirmed_non_secret=True), user_clicked_send=True)
        event, = finish(c)
        self.assertEqual(event.status, "configuration_error")
        self.assertIn("API key", event.text)
        self.assertEqual(c.history(sid), ())
        connection.assert_not_called()

    def test_large_valid_answer_can_be_selected_as_history(self):
        self.t.answer = "a" * 12000
        self.c.send(self.preview(), user_clicked_send=True)
        finish(self.c)
        ids = tuple(m.id for m in self.c.history(self.sid))
        self.assertEqual(len(json.loads(self.preview(history_ids=ids).payload_json)["input"]), 3)

    def test_code_and_draft_remain_inert_text(self):
        self.t.answer = "```python\nraise RuntimeError('must not execute')\n```\nAFT 1"
        self.c.send(self.preview(), user_clicked_send=True)
        event, = finish(self.c)
        self.assertEqual(event.status, "completed")
        self.assertIn("raise RuntimeError", event.text)

    def test_late_response_cancel_clear_delete(self):
        for operation in ("cancel", "clear", "delete_session"):
            with self.subTest(operation=operation):
                gate = threading.Event()
                transport = MockTransport(gate=gate)
                c = ChatController(transport, "mock")
                sid = c.new_session()
                p = c.prepare(sid, "hi", confirmed_non_secret=True)
                c.send(p, user_clicked_send=True)
                self.assertTrue(transport.entered.wait(1))
                getattr(c, operation)(sid)
                gate.set()
                deadline = time.monotonic() + 2
                while c._completed.empty() and time.monotonic() < deadline:
                    time.sleep(.002)
                self.assertFalse(c._completed.empty())
                self.assertEqual(c.poll(), ())
                if operation != "delete_session":
                    self.assertEqual(c.history(sid), ())

    def test_busy_and_capacity_do_not_launch_extra_request(self):
        gate = threading.Event()
        t = MockTransport(gate=gate)
        c = ChatController(t, "mock", max_workers=1)
        sid = c.new_session()
        c.send(c.prepare(sid, "hi", confirmed_non_secret=True), user_clicked_send=True)
        self.assertTrue(t.entered.wait(1))
        try:
            with self.assertRaises(BoundaryError):
                c.prepare(sid, "hi", confirmed_non_secret=True)
            c.cancel(sid)
            p = c.prepare(sid, "hi", confirmed_non_secret=True)
            with self.assertRaises(BoundaryError):
                c.send(p, user_clicked_send=True)
            self.assertEqual(len(t.calls), 1)
        finally:
            gate.set()

    def test_wrong_thread_and_closed_controller(self):
        result = []
        def call():
            try:
                self.c.poll()
            except BoundaryError:
                result.append("blocked")
        t = threading.Thread(target=call)
        t.start()
        t.join()
        self.assertEqual(result, ["blocked"])
        self.c.close()
        with self.assertRaises(BoundaryError):
            self.c.new_session()


class TransportTests(unittest.TestCase):
    def transport(self):
        return OpenAIResponsesTransport(lambda: "testcredential")

    def fake(self, connection, data, status=200):
        response = connection.return_value.getresponse.return_value
        response.status = status
        response.read.return_value = json.dumps(data).encode()

    @patch("app.aft_ai_controller.http.client.HTTPSConnection")
    def test_http_contract_and_output(self, connection):
        self.fake(connection, {"status": "completed", "output": [{"type": "message", "role": "assistant",
                  "content": [{"type": "output_text", "text": "draft"}]}]})
        self.assertEqual(self.transport().complete('{"input":"hello"}'), "draft")
        connection.assert_called_once_with("api.openai.com", timeout=30)
        args, kwargs = connection.return_value.request.call_args
        self.assertEqual(args, ("POST", "/v1/responses"))
        self.assertEqual(kwargs["body"], b'{"input":"hello"}')
        self.assertEqual(kwargs["headers"]["Authorization"], "Bearer testcredential")
        connection.return_value.close.assert_called_once()

    @patch("app.aft_ai_controller.http.client.HTTPSConnection")
    def test_key_absent_or_in_payload_never_connects(self, connection):
        for transport, payload in ((OpenAIResponsesTransport(lambda: ""), "{}"),
                                   (self.transport(), '{"input":"testcredential"}')):
            with self.assertRaises(BoundaryError):
                transport.complete(payload)
        connection.assert_not_called()

    @patch("app.aft_ai_controller.http.client.HTTPSConnection")
    def test_tool_incomplete_malformed_and_echoed_key_rejected(self, connection):
        values = [
            {"status": "completed", "output": [{"type": "function_call"}]},
            {"status": "incomplete", "output": []},
            {"status": "completed", "output": []},
            {"status": "completed", "output": [{"type": "message", "role": "assistant",
             "content": [{"type": "output_text", "text": "testcredential"}]}]},
        ]
        for data in values:
            self.fake(connection, data)
            with self.assertRaisesRegex(BoundaryError, "^Provider request failed$"):
                self.transport().complete("{}")

    @patch("app.aft_ai_controller.http.client.HTTPSConnection")
    def test_redirect_http_error_timeout_no_retry(self, connection):
        for status in (302, 401, 429, 500):
            connection.reset_mock()
            self.fake(connection, {"private": "testcredential"}, status=status)
            with self.assertRaisesRegex(BoundaryError, "^Provider request failed$"):
                self.transport().complete("{}")
            self.assertEqual(connection.return_value.request.call_count, 1)
        connection.return_value.request.side_effect = TimeoutError("testcredential")
        with self.assertRaisesRegex(BoundaryError, "^Provider request failed$"):
            self.transport().complete("{}")


if __name__ == "__main__":
    unittest.main()
