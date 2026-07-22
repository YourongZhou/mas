from __future__ import annotations

import io
import json
import time
from pathlib import Path

from fastapi.testclient import TestClient
from langchain_core.messages import AIMessage
import pytest

from bioagent.config import AgentConfig, load_model_settings, model_settings_path
import bioagent.webapp.app as webapp_module
import bioagent.webapp.state as web_state_module
from bioagent.webapp.app import TEXT_PREVIEW_LIMIT_BYTES, TEXT_PREVIEW_TRUNCATED_MESSAGE, create_app, _read_text_preview
from bioagent.webapp.state import TaskRecord, TaskStore, apply_event_to_record, build_file_tree, record_from_payload, stop_active_docker_containers


def _config(tmp_path: Path) -> AgentConfig:
    repo = tmp_path / "repo"
    mas2 = repo / "mas_2"
    bio = repo / "BioAgent"
    (mas2 / "workflows").mkdir(parents=True)
    (mas2 / "docker").mkdir(parents=True)
    (bio / "logs").mkdir(parents=True)
    (bio / "runs").mkdir(parents=True)
    return AgentConfig(
        repo_root=repo,
        project_root=bio,
        mas2_root=mas2,
        workflows_root=mas2 / "workflows",
        docker_root=mas2 / "docker",
        logs_dir=bio / "logs",
        runs_dir=bio / "runs",
        api_key="test-key",
        base_url="http://example.test/v1",
        model_name="test-model",
        temperature=0.0,
        request_timeout=10.0,
        mimo_thinking_type="",
        chat_template_enable_thinking=None,
    )


def _fake_stream(task: str, *, data_path: str | None = None, result_dir: str | None = None, max_turns: int = 20):
    run_dir = Path(result_dir or data_path or ".").parent / "run_web_test"
    output_dir = run_dir / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis_summary.json").write_text('{"n_cells": 10, "n_clusters": 2}', encoding="utf-8")
    (output_dir / "umap_leiden.png").write_bytes(b"png")

    yield {"type": "run_start", "run_id": "run_web_test", "run_dir": str(run_dir), "log_path": str(run_dir / "run.log")}
    yield {"type": "assistant_message", "content": "I will inspect the skill."}
    yield {
        "type": "tool_call",
        "tool_name": "inspect_workflow_skill",
        "call_id": "call_skill",
        "args": {"skill_id": "scrnaseq-scanpy-core-analysis"},
    }
    yield {
        "type": "tool_result",
        "tool_name": "inspect_workflow_skill",
        "call_id": "call_skill",
        "ok": True,
        "result": {"scripts": ["scripts/qc_metrics.py"]},
    }
    yield {"type": "final", "content": "Analysis completed.", "status": "completed"}
    yield {
        "type": "run_end",
        "run_id": "run_web_test",
        "run_dir": str(run_dir),
        "log_path": str(run_dir / "run.log"),
        "status": "completed",
        "result": {"final_answer": "Analysis completed.", "status": "completed"},
    }


def _fake_stream_with_web_result(task: str, *, data_path: str | None = None, result_dir: str | None = None, max_turns: int = 20):
    web_name = Path(result_dir or "web_unknown").name
    run_dir = Path(result_dir or ".").parent / "run_web_test"
    output_dir = run_dir / "outputs"
    web_dir = run_dir / web_name
    output_dir.mkdir(parents=True, exist_ok=True)
    web_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "analysis_summary.json").write_text("{}", encoding="utf-8")
    (web_dir / "analysis_summary.json").write_text('{"preferred": true}', encoding="utf-8")
    (run_dir / "scripts").mkdir()
    (run_dir / "scripts" / "code.py").write_text("print('internal')", encoding="utf-8")

    yield {"type": "run_start", "run_id": "run_web_test", "run_dir": str(run_dir), "log_path": str(run_dir / "run.log")}
    yield {"type": "final", "content": "Analysis completed.", "status": "completed"}
    yield {
        "type": "run_end",
        "run_id": "run_web_test",
        "run_dir": str(run_dir),
        "log_path": str(run_dir / "run.log"),
        "status": "completed",
        "result": {"final_answer": "Analysis completed.", "status": "completed"},
    }


def _failing_stream(task: str, *, data_path: str | None = None, result_dir: str | None = None, max_turns: int = 20):
    raise RuntimeError("model service is loading")
    yield


def _needs_input_stream(task: str, *, data_path: str | None = None, result_dir: str | None = None, max_turns: int = 20):
    run_dir = Path(result_dir or ".").parent / "run_needs_input"
    run_dir.mkdir(parents=True, exist_ok=True)
    yield {"type": "run_start", "run_id": "run_needs_input", "run_dir": str(run_dir), "log_path": str(run_dir / "run.log")}
    yield {"type": "assistant_message", "content": "Which species should I use?"}
    yield {"type": "final", "content": "Which species should I use?", "status": "needs_user_input"}
    yield {
        "type": "run_end",
        "run_id": "run_needs_input",
        "run_dir": str(run_dir),
        "log_path": str(run_dir / "run.log"),
        "status": "needs_user_input",
        "result": {"final_answer": "Which species should I use?", "status": "needs_user_input"},
    }


def _pausable_stream(
    task: str,
    *,
    data_path: str | None = None,
    result_dir: str | None = None,
    max_turns: int = 20,
    pause_requested=None,
):
    run_dir = Path(result_dir or ".").parent / "run_paused"
    run_dir.mkdir(parents=True, exist_ok=True)
    yield {"type": "run_start", "run_id": "run_paused", "run_dir": str(run_dir), "log_path": str(run_dir / "run.log")}
    assert pause_requested is not None
    deadline = time.time() + 3
    while not pause_requested() and time.time() < deadline:
        time.sleep(0.01)
    assert pause_requested()
    yield {"type": "final", "content": "Task paused. Add instructions to continue.", "status": "paused"}
    yield {
        "type": "run_end",
        "run_id": "run_paused",
        "run_dir": str(run_dir),
        "log_path": str(run_dir / "run.log"),
        "status": "paused",
        "result": {"final_answer": "Task paused. Add instructions to continue.", "status": "paused"},
    }


def _fake_resume(resume_id: str, user_answer: str, *, max_turns: int = 20) -> dict:
    return {
        "final_answer": f"Resumed with {user_answer}",
        "run_dir": "/tmp/resumed_run",
        "log_path": "/tmp/resumed_run/run.log",
        "turns": 2,
        "status": "completed",
    }


def _fake_chat(context: dict, user_message: str) -> str:
    assert context["sessionId"]
    assert context["runStatus"] == "paused"
    assert "availableSkills" in context
    return f"Available skills for this session: {user_message}"


def test_model_settings_api_persists_redacts_and_applies_configuration(tmp_path: Path) -> None:
    tested: list[AgentConfig] = []

    def model_test_runner(config: AgentConfig) -> dict:
        tested.append(config)
        return {"ok": True, "latencyMs": 12, "preview": "OK"}

    config = _config(tmp_path)
    app = create_app(config=config, stream_runner=_fake_stream, model_test_runner=model_test_runner)
    client = TestClient(app)

    initial = client.get("/api/settings/model").json()
    assert initial["provider"] == "openai_compatible"
    assert initial["baseUrl"] == "http://example.test/v1"
    assert initial["apiKeyMasked"] == "te***"
    assert "apiKey" not in initial

    payload = {
        "provider": "anthropic",
        "base_url": "http://10.119.1.246:9010/v1",
        "api_key": "new-secret-key",
        "model_name": "new-qwen-model",
        "temperature": 0.15,
        "request_timeout": 900,
        "mimo_thinking_type": "",
        "chat_template_enable_thinking": False,
    }
    saved = client.put("/api/settings/model", json=payload)
    connection = client.post("/api/settings/model/test", json={**payload, "api_key": ""})
    health = client.get("/api/health").json()

    assert saved.status_code == 200
    assert saved.json()["provider"] == "anthropic"
    assert saved.json()["modelName"] == "new-qwen-model"
    assert saved.json()["apiKeyMasked"] == "new-se...-key"
    assert "new-secret-key" not in saved.text
    assert connection.json()["ok"] is True
    assert tested[-1].api_key == "new-secret-key"
    assert tested[-1].provider == "anthropic"
    assert tested[-1].model_name == "new-qwen-model"
    assert health["model"] == "new-qwen-model"
    assert health["baseUrl"] == "http://10.119.1.246:9010/v1"
    assert load_model_settings(config.project_root)["api_key"] == "new-secret-key"
    assert load_model_settings(config.project_root)["provider"] == "anthropic"
    assert model_settings_path(config.project_root).stat().st_mode & 0o777 == 0o600


def test_model_connection_rejects_invalid_tool_call_payload(tmp_path: Path, monkeypatch) -> None:
    class InvalidToolModel:
        def bind_tools(self, tools, **kwargs):
            return self

        def invoke(self, messages):
            return AIMessage(
                content="",
                invalid_tool_calls=[
                    {
                        "name": "bioagent_connection_probe",
                        "args": '{}{"value":"OK"}',
                        "id": "toolu_invalid",
                        "error": "JSONDecodeError: Extra data",
                    }
                ],
                response_metadata={"finish_reason": "tool_calls"},
            )

    monkeypatch.setattr(webapp_module, "build_llm", lambda config: InvalidToolModel())

    with pytest.raises(RuntimeError, match="invalid tool call"):
        webapp_module._test_model_connection(_config(tmp_path))


def test_model_settings_reset_restores_startup_configuration(tmp_path: Path) -> None:
    config = _config(tmp_path)
    app = create_app(config=config, stream_runner=_fake_stream, model_test_runner=lambda _: {"ok": True})
    client = TestClient(app)
    client.put(
        "/api/settings/model",
        json={
            "provider": "openai_compatible",
            "base_url": "http://new.test/v1",
            "api_key": "new-key",
            "model_name": "new-model",
            "temperature": 0.2,
            "request_timeout": 30,
            "mimo_thinking_type": "disabled",
            "chat_template_enable_thinking": None,
        },
    )

    response = client.delete("/api/settings/model")

    assert response.status_code == 200
    assert response.json()["source"] == "environment"
    assert response.json()["modelName"] == "test-model"
    assert client.get("/api/health").json()["model"] == "test-model"
    assert not model_settings_path(config.project_root).exists()


def test_workbench_exposes_a_dedicated_model_settings_page() -> None:
    root = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static"
    markup = (root / "index.html").read_text(encoding="utf-8")
    script = (root / "assets" / "app.js").read_text(encoding="utf-8")

    assert 'id="settingsNavButton"' in markup
    assert 'id="modelSettingsPage"' in markup
    assert 'id="modelSettingsForm"' in markup
    assert 'id="modelProviderInput"' in markup
    assert 'id="modelBaseUrlInput"' in markup
    assert 'id="modelApiKeyInput" type="password"' in markup
    assert 'id="testModelButton"' in markup
    assert 'id="resetModelSettingsButton"' in markup
    assert 'location.pathname === "/settings"' in script
    assert 'api("/api/settings/model")' in script
    assert 'api("/api/settings/model/test"' in script


def test_session_api_creates_an_open_conversation_with_an_initial_run(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)

    created = client.post("/api/sessions", json={"task": "Run GWAS"}).json()
    snapshot = _wait_for_status(client, created["sessionId"], "completed")

    assert created["sessionId"] == created["taskId"]
    assert snapshot["sessionStatus"] == "open"
    assert snapshot["runStatus"] == "completed"
    assert snapshot["runs"][0]["kind"] == "agent"
    assert snapshot["messages"][0]["role"] == "user"
    assert snapshot["messages"][0]["content"] == "Run GWAS"


def test_session_title_is_generated_asynchronously_and_persisted(tmp_path: Path) -> None:
    def title_runner(config: AgentConfig, task: str) -> str:
        assert config.model_name == "test-model"
        assert task == "Analyze PBMC3K cells"
        return "标题：PBMC3K 细胞分析"

    config = _config(tmp_path)
    app = create_app(config=config, stream_runner=_fake_stream, title_runner=title_runner)
    client = TestClient(app)
    session_id = client.post("/api/sessions", json={"task": "Analyze PBMC3K cells"}).json()["sessionId"]

    deadline = time.time() + 3
    snapshot = client.get(f"/api/sessions/{session_id}").json()
    while time.time() < deadline and snapshot["title"] != "PBMC3K 细胞分析":
        time.sleep(0.02)
        snapshot = client.get(f"/api/sessions/{session_id}").json()

    assert snapshot["title"] == "PBMC3K 细胞分析"
    assert any(event["type"] == "session_title" for event in snapshot["events"])
    stored = json.loads((config.runs_dir / "web_tasks" / f"{session_id}.json").read_text(encoding="utf-8"))
    assert stored["title"] == "PBMC3K 细胞分析"


def test_session_title_keeps_local_fallback_when_generation_fails(tmp_path: Path) -> None:
    def title_runner(config: AgentConfig, task: str) -> str:
        raise RuntimeError("model unavailable")

    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream, title_runner=title_runner)
    client = TestClient(app)
    session_id = client.post("/api/sessions", json={"task": "Analyze PBMC3K cells"}).json()["sessionId"]

    snapshot = _wait_for_status(client, session_id, "completed")

    assert snapshot["title"] == "Analyze PBMC3K cells"
    assert not any(event["type"] == "session_title" for event in snapshot["events"])


def test_session_events_have_stable_protocol_ids_and_tool_parent_links(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    session_id = client.post("/api/sessions", json={"task": "Run GWAS"}).json()["sessionId"]
    snapshot = _wait_for_status(client, session_id, "completed")

    assert snapshot["events"]
    assert len({event["eventId"] for event in snapshot["events"]}) == len(snapshot["events"])
    assert all(event["sessionId"] == session_id for event in snapshot["events"])
    assert all(event["createdAt"] for event in snapshot["events"])
    tool_call = next(event for event in snapshot["events"] if event["type"] == "tool_call")
    tool_result = next(event for event in snapshot["events"] if event["type"] == "tool_result")
    assert tool_result["parentEventId"] == tool_call["eventId"]
    assert tool_result["toolCallId"] == tool_call["toolCallId"] == "call_skill"


def test_explicit_plan_and_structured_question_are_exposed_in_snapshot() -> None:
    record = TaskRecord(id="session-plan", task="Analyze cells", title="Analyze cells", data_path=None, max_turns=20)
    plan = {
        "plan_id": "plan_1",
        "goal": "Analyze cells",
        "status": "awaiting_approval",
        "assumptions": [],
        "steps": [{"id": "step_1", "title": "Inspect data", "status": "pending", "success_criteria": "Input verified"}],
    }
    apply_event_to_record(
        record,
        {
            "type": "tool_result",
            "tool_name": "propose_plan",
            "call_id": "call_plan",
            "ok": True,
            "result": {
                "status": "needs_user_input",
                "interaction_type": "plan_approval",
                "question_id": "question_1",
                "question": "Approve this plan?",
                "options": ["Approve", "Request changes"],
                "allow_free_text": True,
                "plan": plan,
            },
        },
    )

    snapshot = record.snapshot()
    assert snapshot["planState"] == plan
    assert snapshot["plan"][0]["title"] == "Inspect data"
    assert snapshot["activeQuestion"]["questionId"] == "question_1"
    assert snapshot["activeQuestion"]["options"] == ["Approve", "Request changes"]


def test_snapshot_separates_semantic_plan_from_deduplicated_tool_activity() -> None:
    record = TaskRecord(id="session-activity", task="Analyze cells", title="Analyze cells", data_path=None, max_turns=20)
    record.plan = [{"id": "legacy", "title": "Call read_file", "status": "completed"}]
    record.traces = [
        {
            "id": "call_1",
            "eventId": "event_1",
            "toolCallId": "tool_1",
            "toolName": "read_file",
            "status": "running",
            "runId": "run_1",
            "turnId": "run_1:turn:1",
        },
        {
            "id": "call_1",
            "eventId": "event_1",
            "toolCallId": "tool_1",
            "toolName": "read_file",
            "status": "completed",
            "runId": "run_1",
            "turnId": "run_1:turn:1",
        },
        {
            "id": "call_2",
            "eventId": "event_2",
            "toolCallId": "tool_2",
            "toolName": "start_job",
            "status": "running",
            "runId": "run_1",
            "turnId": "run_1:turn:2",
        },
    ]

    snapshot = record.snapshot()

    assert snapshot["plan"] == []
    assert [item["toolName"] for item in snapshot["activity"]] == ["read_file", "start_job"]
    assert snapshot["activity"][0]["status"] == "completed"
    assert snapshot["activity"][0]["traceId"] == "call_1"


def test_paused_session_message_is_answered_without_resuming_the_run(tmp_path: Path) -> None:
    resume_calls: list[tuple[str, str]] = []

    def resume_runner(resume_id: str, user_answer: str, *, max_turns: int = 20) -> dict:
        resume_calls.append((resume_id, user_answer))
        return _fake_resume(resume_id, user_answer, max_turns=max_turns)

    app = create_app(
        config=_config(tmp_path),
        stream_runner=_pausable_stream,
        resume_runner=resume_runner,
        chat_runner=_fake_chat,
    )
    client = TestClient(app)
    session_id = client.post("/api/sessions", json={"task": "Run GWAS"}).json()["sessionId"]
    _wait_for_status(client, session_id, "running")
    client.post(f"/api/sessions/{session_id}/pause")
    _wait_for_status(client, session_id, "paused")

    response = client.post(
        f"/api/sessions/{session_id}/messages",
        json={"content": "what skills do you have?"},
    )
    deadline = time.time() + 3
    snapshot = client.get(f"/api/sessions/{session_id}").json()
    while snapshot["interactionStatus"] and time.time() < deadline:
        time.sleep(0.01)
        snapshot = client.get(f"/api/sessions/{session_id}").json()

    assert response.status_code == 202
    assert snapshot["runStatus"] == "paused"
    assert snapshot["interactionStatus"] == ""
    assert snapshot["messages"][-2]["content"] == "what skills do you have?"
    assert "Available skills" in snapshot["messages"][-1]["content"]
    assert resume_calls == []


def test_running_session_message_is_queued_for_the_active_agent(tmp_path: Path) -> None:
    received: list[str] = []

    def steerable_stream(
        task: str,
        *,
        data_path: str | None = None,
        result_dir: str | None = None,
        max_turns: int = 20,
        pause_requested=None,
        take_pending_messages=None,
    ):
        yield {"type": "run_start", "run_id": "run_steer", "run_dir": str(tmp_path), "log_path": str(tmp_path / "run.log")}
        deadline = time.time() + 3
        while time.time() < deadline and not received:
            received.extend(take_pending_messages())
            time.sleep(0.01)
        yield {"type": "final", "content": "Steering applied.", "status": "completed"}
        yield {"type": "run_end", "run_id": "run_steer", "run_dir": str(tmp_path), "status": "completed"}

    app = create_app(config=_config(tmp_path), stream_runner=steerable_stream)
    client = TestClient(app)
    session_id = client.post("/api/sessions", json={"task": "Run GWAS"}).json()["sessionId"]
    _wait_for_status(client, session_id, "running")

    response = client.post(
        f"/api/sessions/{session_id}/messages",
        json={"content": "Use only chromosome 22."},
    )
    snapshot = _wait_for_status(client, session_id, "completed")

    assert response.status_code == 202
    queued_message = next(message for message in snapshot["messages"] if message["content"] == "Use only chromosome 22.")
    assert received == [
        {
            "type": "user_message",
            "messageId": queued_message["id"],
            "content": "Use only chromosome 22.",
        }
    ]
    assert queued_message["delivery"] == "consumed"


def test_completed_session_message_starts_a_follow_up_run_in_the_same_session(tmp_path: Path) -> None:
    prompts: list[str] = []
    followup_kwargs: list[dict] = []

    def followup_stream(task: str, *, data_path=None, result_dir=None, max_turns=20, **kwargs):
        prompts.append(task)
        followup_kwargs.append(kwargs)
        run_id = f"run_{len(prompts)}"
        yield {"type": "run_start", "run_id": run_id, "run_dir": str(tmp_path / run_id), "log_path": str(tmp_path / f"{run_id}.log")}
        yield {"type": "final", "content": f"Answer {len(prompts)}", "status": "completed"}
        yield {"type": "run_end", "run_id": run_id, "run_dir": str(tmp_path / run_id), "status": "completed"}

    app = create_app(config=_config(tmp_path), stream_runner=followup_stream)
    client = TestClient(app)
    session_id = client.post("/api/sessions", json={"task": "Run GWAS"}).json()["sessionId"]
    _wait_for_status(client, session_id, "completed")

    response = client.post(
        f"/api/sessions/{session_id}/messages",
        json={"content": "Explain the strongest association."},
    )
    deadline = time.time() + 3
    snapshot = client.get(f"/api/sessions/{session_id}").json()
    while len(snapshot["runs"]) < 2 and time.time() < deadline:
        time.sleep(0.01)
        snapshot = client.get(f"/api/sessions/{session_id}").json()
    snapshot = _wait_for_status(client, session_id, "completed")

    assert response.status_code == 202
    assert len(snapshot["runs"]) == 2
    assert snapshot["sessionStatus"] == "open"
    assert "Explain the strongest association." in prompts[-1]
    assert "Answer 1" in prompts[-1]
    assert followup_kwargs[-1]["session_message"] == "Explain the strongest association."
    assert followup_kwargs[-1]["initial_memory_state"]["task_state"] == {}
    assert followup_kwargs[-1]["prior_run_dirs"] == [str(tmp_path / "run_1")]


def test_session_messages_and_runs_survive_store_reload(tmp_path: Path) -> None:
    config = _config(tmp_path)
    app = create_app(config=config, stream_runner=_fake_stream)
    client = TestClient(app)
    session_id = client.post("/api/sessions", json={"task": "Run demo"}).json()["sessionId"]
    _wait_for_status(client, session_id, "completed")
    client.post(f"/api/sessions/{session_id}/messages", json={"content": "Summarize the result."})
    deadline = time.time() + 3
    snapshot = client.get(f"/api/sessions/{session_id}").json()
    while len(snapshot["runs"]) < 2 and time.time() < deadline:
        time.sleep(0.01)
        snapshot = client.get(f"/api/sessions/{session_id}").json()
    _wait_for_status(client, session_id, "completed")

    restored_client = TestClient(create_app(config=config, stream_runner=_fake_stream))
    restored = restored_client.get(f"/api/sessions/{session_id}").json()

    assert restored["sessionStatus"] == "open"
    assert len(restored["runs"]) == 2
    assert any(message["content"] == "Summarize the result." for message in restored["messages"])


def test_workbench_session_composer_separates_send_pause_and_continue_controls() -> None:
    root = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static"
    markup = (root / "index.html").read_text(encoding="utf-8")
    script = (root / "assets" / "app.js").read_text(encoding="utf-8")

    assert 'id="sendButton"' in markup
    assert 'id="pauseButton"' in markup
    assert 'id="interruptButton"' in markup
    assert 'id="cancelJobButton"' in markup
    assert 'id="continueButton"' in markup
    assert 'id="questionPanel"' in markup
    assert "async function sendSessionMessage" in script
    assert "`/api/sessions/${sessionId}/messages`" in script
    assert 'els.pauseButton.addEventListener("click", pauseTask)' in script
    assert 'els.interruptButton.addEventListener("click", interruptCurrentTool)' in script
    assert 'els.cancelJobButton.addEventListener("click", cancelActiveJob)' in script
    assert 'els.continueButton.addEventListener("click", resumeTask)' in script
    assert '`/api/sessions/${taskId}/interrupt`' in script
    assert '`/api/sessions/${taskId}/jobs/cancel`' in script
    assert "snapshot?.activeQuestion" in script
    assert 'els.taskInput.disabled = false' in script


def _wait_for_status(client: TestClient, task_id: str, status: str) -> dict:
    deadline = time.time() + 3
    while time.time() < deadline:
        payload = client.get(f"/api/tasks/{task_id}").json()
        if payload["status"] == status:
            return payload
        time.sleep(0.02)
    raise AssertionError(f"task {task_id} did not reach {status}")


def test_plan_marks_tool_call_step_completed_after_successful_result() -> None:
    record = TaskRecord(id="task1", task="Run demo", title="Run demo", data_path=None, max_turns=20)

    apply_event_to_record(
        record,
        {
            "type": "tool_call",
            "tool_name": "inspect_workflow_skill",
            "call_id": "call_skill",
            "args": {"skill_id": "scrnaseq-scanpy-core-analysis"},
        },
    )
    apply_event_to_record(
        record,
        {
            "type": "tool_result",
            "tool_name": "inspect_workflow_skill",
            "call_id": "call_skill",
            "ok": True,
            "result": {"scripts": ["scripts/qc_metrics.py"]},
        },
    )

    plan = {item["title"]: item["status"] for item in record.plan}

    assert plan["Call inspect_workflow_skill"] == "completed"
    assert plan["Finish inspect_workflow_skill"] == "completed"


def test_plan_uses_trace_tool_name_when_result_event_omits_tool_name() -> None:
    record = TaskRecord(id="task1", task="Run demo", title="Run demo", data_path=None, max_turns=20)

    apply_event_to_record(
        record,
        {
            "type": "tool_call",
            "tool_name": "inspect_workflow_skill",
            "call_id": "call_skill",
            "args": {"skill_id": "scrnaseq-scanpy-core-analysis"},
        },
    )
    apply_event_to_record(
        record,
        {
            "type": "tool_result",
            "call_id": "call_skill",
            "ok": True,
            "result": {"scripts": ["scripts/qc_metrics.py"]},
        },
    )

    plan = {item["title"]: item["status"] for item in record.plan}

    assert plan["Call inspect_workflow_skill"] == "completed"
    assert "Call tool" not in plan


def test_turn_start_updates_compute_progress() -> None:
    record = TaskRecord(id="task1", task="Run demo", title="Run demo", data_path=None, max_turns=20)

    apply_event_to_record(record, {"type": "turn_start", "turn": 7, "max_turns": 20})

    compute = record.snapshot()["compute"]

    assert compute["turn"] == 7
    assert compute["maxTurns"] == 20


def test_start_events_do_not_overwrite_an_early_pause_request() -> None:
    record = TaskRecord(id="task1", task="Run demo", title="Run demo", data_path=None, max_turns=20, status="pausing")

    apply_event_to_record(record, {"type": "task_started", "status": "running"})
    apply_event_to_record(record, {"type": "run_start", "run_id": "run_test", "status": "running"})

    assert record.status == "pausing"


def test_memory_state_event_is_exposed_and_persisted_in_task_snapshot() -> None:
    memory = {
        "taskState": {
            "task": "Run demo",
            "current_stage": "inspecting_skill",
            "selected_skill": "scrnaseq-scanpy-core-analysis",
        },
        "priorEpisodes": [{"skill_id": "scrnaseq-scanpy-core-analysis", "outcome": "success"}],
        "longTermEnabled": True,
        "namespace": ["bioagent", "default"],
    }
    record = TaskRecord(id="task1", task="Run demo", title="Run demo", data_path=None, max_turns=20)

    apply_event_to_record(record, {"type": "memory_state", "memory": memory})
    restored = record_from_payload(
        {
            "id": "task1",
            "task": "Run demo",
            "title": "Run demo",
            "max_turns": 20,
            "memory": record.memory,
        }
    )

    assert record.snapshot()["memory"] == memory
    assert restored.snapshot()["memory"] == memory


def test_history_load_normalizes_stale_running_tool_plan_steps() -> None:
    record = record_from_payload(
        {
            "id": "task1",
            "task": "Run demo",
            "title": "Run demo",
            "status": "completed",
            "traces": [
                {
                    "id": "call_skill",
                    "toolName": "inspect_workflow_skill",
                    "status": "completed",
                    "input": {},
                    "output": {},
                }
            ],
            "plan": [
                {"id": "p1", "title": "Call inspect_workflow_skill", "status": "running"},
                {"id": "p2", "title": "Finish inspect_workflow_skill", "status": "completed"},
            ],
        }
    )

    plan = {item["title"]: item["status"] for item in record.plan}

    assert plan["Call inspect_workflow_skill"] == "completed"


def test_history_load_uses_explicit_tool_output_failure_status() -> None:
    record = record_from_payload(
        {
            "id": "task1",
            "task": "Run demo",
            "title": "Run demo",
            "status": "failed",
            "traces": [
                {
                    "id": "call_python",
                    "toolName": "execute_python",
                    "status": "completed",
                    "input": {},
                    "output": {"ok": False, "error_reason": "Timed out"},
                }
            ],
        }
    )

    assert record.traces[0]["status"] == "failed"


def test_history_load_marks_orphaned_running_task_as_interrupted() -> None:
    record = record_from_payload(
        {
            "id": "task1",
            "task": "Run demo",
            "title": "Run demo",
            "status": "running",
            "traces": [
                {
                    "id": "call_python",
                    "toolName": "execute_python",
                    "status": "running",
                    "input": {"code": "print('hello')"},
                    "output": None,
                }
            ],
            "compute": {"currentTool": "execute_python"},
        }
    )

    assert record.status == "failed"
    assert record.error == "Run interrupted before completion."
    assert record.final_answer == "Run failed: Run interrupted before completion."
    assert record.traces[0]["status"] == "failed"
    assert record.compute["currentTool"] == ""


def test_run_error_replaces_stale_paused_final_answer() -> None:
    record = TaskRecord(
        id="task1",
        task="Run demo",
        title="Run demo",
        data_path=None,
        max_turns=20,
        status="running",
        final_answer="Task paused. Add instructions to continue.",
    )

    apply_event_to_record(record, {"type": "error", "error": "Context size has been exceeded."})

    assert record.status == "failed"
    assert record.final_answer == "Run failed: Context size has been exceeded."
    assert record.messages[-1]["content"] == record.final_answer


def test_task_api_turns_bioagent_stream_events_into_snapshot(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)

    created = client.post(
        "/api/tasks",
        json={"task": "Run single-cell demo", "data_path": "mas_2/data/bmmc_b_cell.h5ad", "max_turns": 3},
    ).json()
    assert created["taskId"]
    assert created["redirectUrl"] == f"/tasks/{created['taskId']}"

    snapshot = _wait_for_status(client, created["taskId"], "completed")

    assert snapshot["title"] == "Run single-cell demo"
    assert snapshot["dataPath"] == "mas_2/data/bmmc_b_cell.h5ad"
    assert snapshot["maxTurns"] == 3
    assert snapshot["messages"][-1]["content"] == "Analysis completed."
    assert snapshot["traces"][0]["toolName"] == "inspect_workflow_skill"
    assert snapshot["traces"][0]["status"] == "completed"
    assert snapshot["traces"][0]["input"]["skill_id"] == "scrnaseq-scanpy-core-analysis"
    assert snapshot["traces"][0]["output"]["scripts"] == ["scripts/qc_metrics.py"]
    assert snapshot["compute"]["runId"] == "run_web_test"
    assert snapshot["compute"]["resultRoot"].endswith("/outputs")


def test_task_files_api_lists_outputs_under_run_dir(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "completed")

    files = client.get(f"/api/tasks/{task_id}/files").json()["files"]

    assert [node["name"] for node in files] == ["analysis_summary.json", "umap_leiden.png"]


def test_task_file_content_accepts_result_relative_paths(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "completed")

    response = client.get(f"/api/tasks/{task_id}/files/content", params={"path": "analysis_summary.json"})

    assert response.status_code == 200
    assert response.json()["n_cells"] == 10


def test_task_file_content_caps_large_text_previews_but_downloads_full_file(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    snapshot = _wait_for_status(client, task_id, "completed")
    result_root = Path(snapshot["compute"]["resultRoot"])
    full_text = "x" * 220_000
    (result_root / "large.log").write_text(full_text, encoding="utf-8")

    preview = client.get(f"/api/tasks/{task_id}/files/content", params={"path": "large.log"})
    download = client.get(f"/api/tasks/{task_id}/files/download", params={"path": "large.log"})

    assert preview.status_code == 200
    assert len(preview.text) < len(full_text)
    assert preview.text.startswith("x" * 100)
    assert "Preview truncated" in preview.text
    assert download.status_code == 200
    assert download.text == full_text


def test_read_text_preview_reads_only_preview_window() -> None:
    read_sizes: list[int] = []

    class ProbePath:
        def open(self, mode: str):
            assert mode == "rb"

            class ProbeFile(io.BytesIO):
                def read(self, size: int = -1) -> bytes:
                    read_sizes.append(size)
                    return super().read(size)

            return ProbeFile(b"x" * (TEXT_PREVIEW_LIMIT_BYTES + 100))

        def read_bytes(self) -> bytes:
            raise AssertionError("preview should not read the full file")

    preview = _read_text_preview(ProbePath())  # type: ignore[arg-type]

    assert read_sizes == [TEXT_PREVIEW_LIMIT_BYTES + 1]
    assert len(preview) == TEXT_PREVIEW_LIMIT_BYTES + len(TEXT_PREVIEW_TRUNCATED_MESSAGE)
    assert preview.endswith(TEXT_PREVIEW_TRUNCATED_MESSAGE)


def test_task_log_content_supports_logs_outside_run_dir(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "completed")
    log_path = tmp_path / "bioagent.log"
    log_path.write_text("first line\nsecond line\n", encoding="utf-8")
    app.state.task_store.get_task(task_id).log_path = str(log_path)

    preview = client.get(f"/api/tasks/{task_id}/log/content")
    download = client.get(f"/api/tasks/{task_id}/log/download")

    assert preview.status_code == 200
    assert preview.text == "first line\nsecond line\n"
    assert download.status_code == 200
    assert download.text == "first line\nsecond line\n"


def test_result_file_tree_prioritizes_primary_outputs(tmp_path: Path) -> None:
    output_dir = tmp_path / "outputs"
    output_dir.mkdir()
    for name in [
        "bmmc_b_cell_processed.h5ad",
        "cluster_sizes.png",
        "pca_variance.png",
        "qc_metrics.png",
        "summary.json",
        "umap_leiden.png",
    ]:
        (output_dir / name).write_text("x", encoding="utf-8")

    files = build_file_tree(output_dir)

    assert [node["name"] for node in files] == [
        "summary.json",
        "umap_leiden.png",
        "qc_metrics.png",
        "pca_variance.png",
        "cluster_sizes.png",
        "bmmc_b_cell_processed.h5ad",
    ]


def test_task_files_api_prefers_web_result_dir_when_present(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream_with_web_result)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "completed")

    snapshot = client.get(f"/api/tasks/{task_id}").json()

    assert [node["name"] for node in snapshot["resultFiles"]] == ["analysis_summary.json"]
    assert snapshot["compute"]["resultRoot"].endswith(f"/web_{task_id}")


def test_raw_upload_is_size_limited_and_attached_to_a_new_session(tmp_path: Path) -> None:
    received_tasks: list[str] = []

    def capture_stream(task: str, *, result_dir=None, **kwargs):
        received_tasks.append(task)
        run_dir = Path(result_dir).parent / "run_upload"
        yield {"type": "run_start", "run_id": "run_upload", "run_dir": str(run_dir), "log_path": str(run_dir / "run.log")}
        yield {"type": "final", "content": "Upload received.", "status": "completed"}
        yield {"type": "run_end", "run_id": "run_upload", "run_dir": str(run_dir), "status": "completed"}

    config = _config(tmp_path)
    app = create_app(
        config=config,
        stream_runner=capture_stream,
        start_job_watcher=False,
        max_upload_bytes=16,
    )
    client = TestClient(app)

    health = client.get("/api/health").json()
    uploaded = client.post(
        "/api/uploads",
        params={"filename": "cells.csv"},
        content=b"gene,count\nA,1\n",
        headers={"Content-Type": "text/csv"},
    )
    rejected = client.post(
        "/api/uploads",
        params={"filename": "too-large.bin"},
        content=b"x" * 17,
    )

    assert health["maxUploadBytes"] == 16
    assert uploaded.status_code == 200
    attachment = uploaded.json()
    assert attachment["uploadId"].startswith("upload_")
    assert attachment["name"] == "cells.csv"
    assert attachment["size"] == 15
    assert Path(attachment["path"]).read_bytes() == b"gene,count\nA,1\n"
    assert attachment["dockerPath"].startswith("/repo/BioAgent/uploads/")
    assert rejected.status_code == 413
    assert not list((config.project_root / "uploads").rglob("too-large.bin"))

    created = client.post(
        "/api/sessions",
        json={
            "task": "Analyze the uploaded table.",
            "data_path": None,
            "max_turns": 5,
            "attachment_ids": [attachment["uploadId"]],
        },
    ).json()
    snapshot = _wait_for_status(client, created["sessionId"], "completed")

    assert snapshot["dataPath"] == attachment["path"]
    assert snapshot["messages"][0]["attachments"][0]["uploadId"] == attachment["uploadId"]
    assert "<bioagent_attachments>" in received_tasks[0]
    assert attachment["path"] in received_tasks[0]
    assert attachment["dockerPath"] in received_tasks[0]


def test_existing_session_upload_is_bound_to_that_session(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream, start_job_watcher=False)
    client = TestClient(app)
    task_id = client.post("/api/sessions", json={"task": "Run demo"}).json()["sessionId"]
    _wait_for_status(client, task_id, "completed")

    uploaded = client.post(
        "/api/uploads",
        params={"filename": "notes.txt", "session_id": task_id},
        content=b"follow-up evidence",
    )
    missing = client.post(
        "/api/uploads",
        params={"filename": "notes.txt", "session_id": "missing-session"},
        content=b"no",
    )

    assert uploaded.status_code == 200
    assert f"/uploads/{task_id}/" in uploaded.json()["path"]
    assert missing.status_code == 404

    attachment = uploaded.json()
    sent = client.post(
        f"/api/sessions/{task_id}/messages",
        json={
            "content": "Use these notes.",
            "max_turns": 5,
            "attachment_ids": [attachment["uploadId"]],
        },
    )
    messages = client.get(f"/api/sessions/{task_id}").json()["messages"]
    message = next(item for item in reversed(messages) if item["role"] == "user")

    assert sent.status_code == 202
    assert message["role"] == "user"
    assert message["attachments"][0]["uploadId"] == attachment["uploadId"]
    assert "<bioagent_attachments>" in message["content"]


def test_workbench_composer_supports_click_and_drop_uploads() -> None:
    root = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static"
    markup = (root / "index.html").read_text(encoding="utf-8")
    script = (root / "assets" / "app.js").read_text(encoding="utf-8")

    assert 'id="attachmentInput"' in markup
    assert 'id="attachButton"' in markup
    assert 'id="attachmentList"' in markup
    assert 'multiple' in markup
    assert 'addEventListener("drop"' in script
    assert 'addEventListener("dragover"' in script
    assert 'uploadFiles(' in script
    assert 'attachment_ids' in script
    upload = script.index("async function uploadFiles")
    upload_finished = script.index("state.uploadBusy = false", upload)
    composer_restored = script.index("renderComposer(state.snapshot)", upload_finished)
    attachments_restored = script.index("renderAttachmentComposer()", composer_restored)
    assert upload_finished < composer_restored < attachments_restored


def test_session_results_survive_follow_up_run_without_outputs(tmp_path: Path) -> None:
    config = _config(tmp_path)
    first_run = config.runs_dir / "run_first"
    second_run = config.runs_dir / "run_follow_up"
    output = first_run / "outputs" / "umap.png"
    output.parent.mkdir(parents=True)
    output.write_bytes(b"png")
    (second_run / "state").mkdir(parents=True)
    (second_run / "state" / "final_verification.json").write_text("{}", encoding="utf-8")
    record = TaskRecord(
        id="session-results",
        task="Analyze PBMC3K",
        title="Analyze PBMC3K",
        data_path=None,
        max_turns=20,
        status="completed",
        run_dir=str(second_run),
        runs=[
            {"id": "first", "runDir": str(first_run), "status": "failed"},
            {"id": "follow-up", "runDir": str(second_run), "status": "completed"},
        ],
        active_run_id="follow-up",
    )
    store = TaskStore(config=config, start_job_watcher=False)
    store._persist_record(record)
    app = create_app(config=config, start_job_watcher=False)
    client = TestClient(app)

    snapshot = client.get(f"/api/sessions/{record.id}").json()
    preview = client.get(
        f"/api/tasks/{record.id}/files/download",
        params={"path": str(output)},
    )

    assert [node["name"] for node in snapshot["resultFiles"]] == ["umap.png"]
    assert snapshot["compute"]["resultRoot"] == str(first_run / "outputs")
    assert preview.status_code == 200
    assert preview.content == b"png"


def test_runner_error_is_visible_in_task_conversation(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_failing_stream)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]

    snapshot = _wait_for_status(client, task_id, "failed")

    assert snapshot["error"] == "model service is loading"
    assert "model service is loading" in snapshot["finalAnswer"]
    assert "model service is loading" in snapshot["messages"][-1]["content"]


def test_task_resume_continues_same_human_input_run(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_needs_input_stream, resume_runner=_fake_resume)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "needs_user_input")

    resumed = client.post(f"/api/tasks/{task_id}/resume", json={"user_answer": "human", "max_turns": 5}).json()
    snapshot = _wait_for_status(client, task_id, "completed")

    assert resumed["taskId"] == task_id
    assert snapshot["messages"][-2]["role"] == "user"
    assert snapshot["messages"][-2]["content"] == "human"
    assert snapshot["finalAnswer"] == "Resumed with human"
    assert snapshot["messages"][-1]["content"] == "Resumed with human"


def test_task_can_pause_and_resume_same_run_with_user_instruction(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_pausable_stream, resume_runner=_fake_resume)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "running")

    pause_response = client.post(f"/api/tasks/{task_id}/pause")
    paused = _wait_for_status(client, task_id, "paused")
    resumed = client.post(
        f"/api/tasks/{task_id}/resume",
        json={"user_answer": "Use a faster vectorized method", "max_turns": 5},
    ).json()
    completed = _wait_for_status(client, task_id, "completed")

    assert pause_response.status_code == 200
    assert pause_response.json()["status"] in {"pausing", "paused"}
    assert paused["status"] == "paused"
    assert resumed["taskId"] == task_id
    assert completed["messages"][-2]["content"] == "Use a faster vectorized method"


def test_pause_only_pauses_agent_and_does_not_cancel_running_compute(tmp_path: Path, monkeypatch) -> None:
    stopped: list[Path] = []
    monkeypatch.setattr(web_state_module, "stop_active_docker_containers", lambda path: stopped.append(path))
    app = create_app(config=_config(tmp_path), stream_runner=_pausable_stream, resume_runner=_fake_resume)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "running")

    response = client.post(f"/api/sessions/{task_id}/pause")
    _wait_for_status(client, task_id, "paused")

    assert response.status_code == 200
    assert stopped == []


def test_terminal_job_callback_is_persisted_once_without_becoming_user_message(tmp_path: Path) -> None:
    config = _config(tmp_path)
    store = TaskStore(config=config, start_job_watcher=False)
    run_dir = config.runs_dir / "run_callback"
    job_dir = run_dir / "jobs" / "job_20260717_120000_abcdef"
    job_dir.mkdir(parents=True)
    (job_dir / "job.json").write_text(
        json.dumps(
            {
                "job_id": "job_20260717_120000_abcdef",
                "run_id": run_dir.name,
                "run_dir": str(run_dir),
                "status": "completed",
                "callback_state": "pending",
                "exit_code": 0,
                "finished_at": "2026-07-17T12:01:00+00:00",
            }
        ),
        encoding="utf-8",
    )
    record = TaskRecord(
        id="session-callback",
        task="Analyze cells",
        title="Analyze cells",
        data_path=None,
        max_turns=20,
        status="paused",
        run_dir=str(run_dir),
        runs=[{"id": "session-run", "runDir": str(run_dir), "status": "paused"}],
        active_run_id="session-run",
    )
    store._tasks[record.id] = record

    assert store.scan_job_callbacks() == 1
    assert store.scan_job_callbacks() == 0
    callback = next(event for event in record.events if event["type"] == "job_callback")
    assert callback["callbackId"].startswith("callback_")
    assert callback["jobId"] == "job_20260717_120000_abcdef"
    assert callback["status"] == "completed"
    assert record.pending_context_events[0]["type"] == "job_callback"
    assert all(message.get("role") != "user" for message in record.messages)


def test_terminal_job_callback_automatically_continues_an_idle_session(tmp_path: Path) -> None:
    received_context: list[dict] = []

    def callback_stream(task: str, *, result_dir=None, take_pending_messages=None, **kwargs):
        received_context.extend(item for item in take_pending_messages() if isinstance(item, dict))
        run_dir = Path(result_dir).parent / "run_after_callback"
        yield {"type": "run_start", "run_id": "run_after_callback", "run_dir": str(run_dir), "log_path": str(run_dir / "run.log")}
        yield {"type": "final", "content": "Callback handled.", "status": "completed"}
        yield {"type": "run_end", "run_id": "run_after_callback", "run_dir": str(run_dir), "status": "completed"}

    config = _config(tmp_path)
    store = TaskStore(config=config, stream_runner=callback_stream, start_job_watcher=False)
    run_dir = config.runs_dir / "run_callback_idle"
    job_dir = run_dir / "jobs" / "job_20260717_120000_fedcba"
    job_dir.mkdir(parents=True)
    (job_dir / "job.json").write_text(
        json.dumps(
            {
                "job_id": "job_20260717_120000_fedcba",
                "run_id": run_dir.name,
                "run_dir": str(run_dir),
                "status": "completed",
                "callback_state": "pending",
                "exit_code": 0,
            }
        ),
        encoding="utf-8",
    )
    record = TaskRecord(
        id="session-callback-idle",
        task="Analyze cells",
        title="Analyze cells",
        data_path=None,
        max_turns=20,
        status="completed",
        run_dir=str(run_dir),
        runs=[{"id": "first-run", "runDir": str(run_dir), "status": "completed"}],
        active_run_id="first-run",
    )
    store._tasks[record.id] = record

    assert store.scan_job_callbacks() == 1
    deadline = time.time() + 3
    while record.status != "completed" and time.time() < deadline:
        time.sleep(0.01)

    assert len(record.runs) == 2
    assert record.status == "completed"
    assert received_context[0]["type"] == "job_callback"
    assert received_context[0]["jobId"] == "job_20260717_120000_fedcba"
    assert record.final_answer == "Callback handled."


def test_pause_stops_only_containers_named_by_run_cidfiles(tmp_path: Path, monkeypatch) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / ".docker-python.cid").write_text("python-container\n", encoding="utf-8")
    job_dir = run_dir / "jobs" / "job_active"
    job_dir.mkdir(parents=True)
    (job_dir / "container.cid").write_text("job-container\n", encoding="utf-8")
    (job_dir / "job.json").write_text(
        json.dumps({"job_id": "job_active", "status": "running", "container_id": "job-container"}),
        encoding="utf-8",
    )
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        return None

    monkeypatch.setattr(web_state_module.subprocess, "run", fake_run)

    stop_active_docker_containers(run_dir)

    assert calls == [
        ["docker", "rm", "-f", "python-container"],
        ["docker", "rm", "-f", "job-container"],
    ]
    assert json.loads((job_dir / "job.json").read_text(encoding="utf-8"))["status"] == "cancelled"


def test_task_resume_does_not_turn_interaction_lifecycle_into_a_semantic_plan(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_needs_input_stream, resume_runner=_fake_resume)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "needs_user_input")

    client.post(f"/api/tasks/{task_id}/resume", json={"user_answer": "human", "max_turns": 5})
    snapshot = _wait_for_status(client, task_id, "completed")
    assert snapshot["status"] == "completed"
    assert snapshot["plan"] == []


def test_task_resume_streams_updated_memory_state(tmp_path: Path) -> None:
    def memory_resume(resume_id: str, user_answer: str, *, max_turns: int = 20, event_sink=None) -> dict:
        assert event_sink is not None
        event_sink(
            {
                "type": "memory_state",
                "memory": {
                    "taskState": {"current_stage": "resuming", "current_goal": user_answer},
                    "priorEpisodes": [],
                    "longTermEnabled": True,
                    "namespace": ["bioagent", "default"],
                },
            }
        )
        return {
            "final_answer": "Resumed",
            "run_dir": "/tmp/resumed_run",
            "log_path": "/tmp/resumed_run/run.log",
            "turns": 1,
            "status": "completed",
        }

    app = create_app(config=_config(tmp_path), stream_runner=_needs_input_stream, resume_runner=memory_resume)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "needs_user_input")

    client.post(f"/api/tasks/{task_id}/resume", json={"user_answer": "Continue without doublets", "max_turns": 5})
    snapshot = _wait_for_status(client, task_id, "completed")

    assert snapshot["memory"]["taskState"]["current_stage"] == "resuming"
    assert snapshot["memory"]["taskState"]["current_goal"] == "Continue without doublets"


def test_waiting_for_human_input_message_is_not_marked_final(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_needs_input_stream, resume_runner=_fake_resume)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]

    snapshot = _wait_for_status(client, task_id, "needs_user_input")

    assert snapshot["messages"][-1]["content"] == "Which species should I use?"
    assert snapshot["messages"][-1].get("final") is not True


def test_waiting_for_human_input_does_not_duplicate_prompt_message(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_needs_input_stream, resume_runner=_fake_resume)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]

    snapshot = _wait_for_status(client, task_id, "needs_user_input")

    prompts = [message for message in snapshot["messages"] if message["content"] == "Which species should I use?"]
    assert len(prompts) == 1


def test_waiting_for_human_input_does_not_create_a_semantic_plan(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_needs_input_stream, resume_runner=_fake_resume)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]

    snapshot = _wait_for_status(client, task_id, "needs_user_input")
    assert snapshot["status"] == "needs_user_input"
    assert snapshot["plan"] == []


def test_workbench_serves_static_index(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)

    response = client.get("/")

    assert response.status_code == 200
    assert "BioAgent Workbench" in response.text


def test_workbench_static_pages_support_head_requests(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)

    index = client.head("/")
    task_page = client.head("/tasks/example")

    assert index.status_code == 200
    assert task_page.status_code == 200


def test_workbench_favicon_does_not_404(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)

    response = client.get("/favicon.ico")

    assert response.status_code == 204


def test_workbench_serves_static_assets_alias_for_cached_paths(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)

    response = client.get("/static/assets/app.js")

    assert response.status_code == 200
    assert "const state" in response.text


def test_workbench_disables_only_unimplemented_rail_navigation() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert '<button class="rail-button active" id="workbenchNavButton" title="Projects" aria-current="page">P</button>' in markup
    assert '<button class="rail-button" title="Resources (coming soon)" disabled>R</button>' in markup
    assert '<button class="rail-button" title="Skills (coming soon)" disabled>S</button>' in markup
    assert '<button class="rail-button" title="History (coming soon)" disabled>H</button>' in markup
    assert '<button class="rail-button" id="settingsNavButton" title="Model settings">G</button>' in markup
    assert ".rail-button:disabled" in styles
    assert "cursor: default" in styles[styles.index(".rail-button:disabled"):]
    assert "opacity: .45" in styles[styles.index(".rail-button:disabled"):]


def test_workbench_js_closes_terminal_event_stream_and_throttles_refresh() -> None:
    source = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js"
    script = source.read_text(encoding="utf-8")

    assert "function isTerminalStatus(status)" in script
    assert "scheduleSnapshotRefresh(taskId)" in script
    assert "if (shouldCloseSessionStream(state.snapshot))" in script
    assert "closeEvents();" in script


def test_workbench_ignores_stale_realtime_stream_events() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    snapshot_listener = script.index('source.addEventListener("task_snapshot"')
    snapshot_guard = script.index("if (state.taskId !== taskId) return", snapshot_listener)
    snapshot_write = script.index("state.snapshot = JSON.parse(event.data)", snapshot_listener)
    error_handler = script.index("source.onerror", snapshot_write)
    error_guard = script.index("if (state.taskId !== taskId) return", error_handler)
    error_status_check = script.index("if (shouldCloseSessionStream(state.snapshot))", error_handler)

    assert snapshot_guard < snapshot_write
    assert error_guard < error_status_check


def test_workbench_ignores_stale_scheduled_snapshot_refresh_responses() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    refresh_start = script.index("function scheduleSnapshotRefresh")
    fetch_snapshot = script.index("const snapshot = await api", refresh_start)
    stale_guard = script.index("if (state.taskId !== taskId) return", fetch_snapshot)
    set_snapshot = script.index("state.snapshot = snapshot", fetch_snapshot)
    render_call = script.index("render();", set_snapshot)

    assert fetch_snapshot < stale_guard < set_snapshot < render_call


def test_workbench_js_filters_tasks_and_result_files_as_user_types() -> None:
    source = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js"
    script = source.read_text(encoding="utf-8")

    assert 'els.taskSearch.addEventListener("input"' in script
    assert 'els.fileSearch.addEventListener("input"' in script
    assert "state.tasks" in script
    assert "filterFileNodes" in script


def test_workbench_js_supports_resume_mode_for_human_input() -> None:
    source = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js"
    script = source.read_text(encoding="utf-8")

    assert "async function resumeTask" in script
    assert "/resume" in script
    assert "needs_user_input" in script
    assert "Reply" in script


def test_workbench_composer_keeps_terminal_sessions_open_for_follow_up() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    composer = script.index("function renderComposer")
    terminal_check = script.index("isTerminalStatus(snapshot?.status)", composer)
    send_label = script.index('els.sendButton.textContent = "Send"', terminal_check)
    default_label = script.index('els.sendButton.textContent = "Run"', send_label)
    needs_input = script.index('snapshot?.status === "needs_user_input"', composer)
    active_status = script.index("isActiveStatus(snapshot?.status)", needs_input)

    assert needs_input < active_status < terminal_check < send_label < default_label


def test_workbench_js_auto_selects_most_relevant_trace() -> None:
    source = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js"
    script = source.read_text(encoding="utf-8")

    assert "function chooseTraceId" in script
    assert 'trace.status === "failed"' in script
    assert 'trace.status === "running"' in script
    assert "chooseTraceId(state.snapshot" in script


def test_workbench_topbar_shows_current_tool_and_links_to_trace() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert 'id="currentToolPill"' in markup
    assert "function renderCurrentTool" in script
    assert "compute.currentTool" in script
    assert 'if (!trace && !currentTool)' in script
    assert 'els.currentToolPill.textContent = "No tool trace"' in script
    assert "currentToolPill.addEventListener" in script
    assert 'activateTab("trace")' in script
    assert ".tool-pill:disabled" in styles
    assert "background: #f8fafc" in styles[styles.index(".tool-pill:disabled"):]


def test_workbench_humanizes_tool_names_in_trace_surfaces() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    current_tool = script.index("function renderCurrentTool")
    topbar_label = script.index("formatToolName(label)", current_tool)
    trace_browser = script.index("function renderTraceBrowser")
    row_label = script.index("formatToolName(trace.toolName)", trace_browser)
    detail = script.index("function renderTraceDetail")
    detail_title = script.index("formatToolName(trace.toolName)", detail)
    activity = script.index("function renderActivity")
    activity_label = script.index("formatToolName(item.toolName)", activity)
    helper = script.index("function formatToolName")
    humanize = script.index("humanizeToolName", helper)

    assert current_tool < topbar_label < trace_browser < row_label < detail < detail_title < activity < activity_label
    assert helper < humanize


def test_workbench_trace_surfaces_use_human_readable_status_labels() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    trace_browser = script.index("function renderTraceBrowser")
    row_status = script.index("formatTraceMeta(trace, snapshot)", trace_browser)
    row_status_class = script.index("formatTraceStatusClass(trace, snapshot)", row_status)
    row_dot = script.index("status-dot ${escapeHtml(statusClass)}", row_status_class)
    detail = script.index("function renderTraceDetail")
    detail_meta = script.index("els.traceMeta.textContent", detail)
    detail_status = script.index("formatTraceMeta(trace, state.snapshot)", detail_meta)
    meta_helper = script.index("function formatTraceMeta")
    meta_status = script.index("formatTraceStatus(trace, snapshot)", meta_helper)
    helper = script.index("function formatTraceStatus")
    recovered = script.index('return "recovered"', helper)
    fallback = script.index("formatStatusLabel(trace.status)", recovered)
    class_helper = script.index("function formatTraceStatusClass")
    class_recovered = script.index('return "recovered"', class_helper)
    class_fallback = script.index('replaceAll("_", "-")', class_recovered)

    assert trace_browser < row_status < row_status_class < row_dot < detail < detail_meta < detail_status
    assert detail < meta_helper < meta_status < helper
    assert helper < recovered < fallback
    assert helper < class_helper < class_recovered < class_fallback
    assert ".status-dot.recovered" in styles


def test_workbench_trace_rows_and_detail_show_tool_duration() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    trace_browser = script.index("function renderTraceBrowser")
    row_meta = script.index("formatTraceMeta(trace, snapshot)", trace_browser)
    detail = script.index("function renderTraceDetail")
    detail_meta = script.index("formatTraceMeta(trace, state.snapshot)", detail)
    helper = script.index("function formatTraceMeta")
    status = script.index("formatTraceStatus(trace, snapshot)", helper)
    duration = script.index("formatRunDuration(trace)", status)
    return_value = script.index("`${status} · ${duration}`", duration)

    assert trace_browser < row_meta < detail < detail_meta < helper
    assert helper < status < duration < return_value


def test_workbench_trace_tool_names_do_not_push_status_offscreen() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    trace_row = styles.index(".trace-list-copy {")
    tool_name = styles.index(".trace-list-copy strong", trace_row)
    tool_block = styles[tool_name:styles.index("}", tool_name)]

    assert "min-width: 0" in tool_block
    assert "overflow: hidden" in tool_block
    assert "text-overflow: ellipsis" in tool_block
    assert "white-space: nowrap" in tool_block


def test_workbench_topbar_labels_terminal_trace_as_last_tool() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    current_tool = script.index("function renderCurrentTool")
    prefix_const = script.index("const prefix = currentToolLabelPrefix(currentTool, snapshot)", current_tool)
    text_content = script.index("els.currentToolPill.textContent = `${prefix}:", prefix_const)
    helper = script.index("function currentToolLabelPrefix")
    running_case = script.index('return "Running"', helper)
    terminal_case = script.index("isTerminalStatus(snapshot?.status)", running_case)
    last_tool = script.index('return "Last tool"', terminal_case)
    default_tool = script.index('return "Tool"', last_tool)

    assert prefix_const < text_content
    assert helper < running_case < terminal_case < last_tool < default_tool


def test_workbench_topbar_prioritizes_failed_error_summary() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    current_tool = script.index("function renderCurrentTool")
    failure = script.index("const failure = formatSnapshotFailure(snapshot)", current_tool)
    failed_text = script.index('els.currentToolPill.textContent = `Failed: ${failure}`', failure)
    failed_class = script.index('els.currentToolPill.className = "tool-pill failed"', failed_text)
    helper = script.index("function formatSnapshotFailure")
    failed_status = script.index('snapshot?.status !== "failed"', helper)
    error_summary = script.index("formatErrorSummary(snapshot.error || finalAnswer)", failed_status)

    assert current_tool < failure < failed_text < failed_class
    assert helper < failed_status < error_summary


def test_workbench_topbar_surfaces_waiting_for_human_reply() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    current_tool = script.index("function renderCurrentTool")
    waiting_check = script.index('snapshot?.status === "needs_user_input"', current_tool)
    waiting_text = script.index('els.currentToolPill.textContent = "Waiting for reply"', waiting_check)
    waiting_class = script.index('els.currentToolPill.className = "tool-pill needs-input"', waiting_text)
    waiting_enabled = script.index("els.currentToolPill.disabled = false", waiting_class)

    assert current_tool < waiting_check < waiting_text < waiting_class < waiting_enabled
    assert ".tool-pill.needs-input" in styles


def test_workbench_topbar_trace_link_marks_selected_trace_row() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    click_start = script.index("els.currentToolPill.addEventListener")
    click_end = script.index("els.refreshButton.addEventListener", click_start)
    handler = script[click_start:click_end]
    select_trace = handler.index("state.selectedTraceId = trace.id")
    mark_row = handler.index("markSelectedTraceRow()")
    render_detail = handler.index("renderTraceDetail()")
    activate_trace = handler.index('activateTab("trace", { user: true })')

    assert select_trace < mark_row < render_detail < activate_trace


def test_workbench_topbar_waiting_for_reply_focuses_composer() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    click_start = script.index("els.currentToolPill.addEventListener")
    click_end = script.index("els.refreshButton.addEventListener", click_start)
    handler = script[click_start:click_end]
    waiting_check = handler.index('state.snapshot?.status === "needs_user_input"')
    focus_input = handler.index("els.taskInput.focus()", waiting_check)
    trace_lookup = handler.index("const trace = topbarTrace(state.snapshot)", focus_input)

    assert waiting_check < focus_input < trace_lookup


def test_workbench_topbar_does_not_surface_recovered_failed_trace_as_current() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "snapshot?.status === \"failed\"" in script
    assert "const failed = snapshot?.status === \"failed\"" in script
    assert "return running || failed || traces[traces.length - 1] || null" in script


def test_workbench_topbar_marks_completed_failed_trace_as_recovered() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    current_tool = script.index("function renderCurrentTool")
    status_const = script.index('const status = currentTool ? "running" : trace ? formatTraceStatusClass(trace, snapshot) : "idle"', current_tool)
    class_set = script.index("els.currentToolPill.className = `tool-pill ${status}`", status_const)

    assert current_tool < status_const < class_set
    assert ".tool-pill.recovered" in styles


def test_workbench_default_trace_selection_ignores_recovered_failures() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    chooser = script.index("function chooseTraceId")
    running = script.index('trace.status === "running"', chooser)
    failed_gate = script.index('snapshot?.status === "failed"', running)
    failed_lookup = script.index('trace.status === "failed"', failed_gate)
    last_trace = script.index("traces[traces.length - 1].id", failed_lookup)

    assert chooser < running < failed_gate < failed_lookup < last_trace


def test_workbench_marks_selected_trace_row() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert 'trace.id === state.selectedTraceId ? "selected" : ""' in script
    assert "function markSelectedTraceRow" in script
    assert "markSelectedTraceRow();" in script
    assert ".trace-list-row.selected" in styles
    assert ".activity-row.selected" in styles


def test_workbench_result_cards_open_file_preview() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert "openFile(image)" in script
    assert "openFile(summary)" in script
    assert ".plot-card, .summary-card" in styles
    assert "cursor: pointer" in styles


def test_workbench_resets_file_preview_between_tasks() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert "function resetFilePreview" in script
    assert "resetFilePreview();" in script
    assert 'els.downloadLink.classList.add("disabled")' in script
    assert 'els.downloadLink.classList.remove("disabled")' in script
    assert ".icon-button.disabled" in styles


def test_workbench_native_disabled_icon_buttons_match_disabled_links() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    selector = ".icon-button.disabled, .icon-button:disabled"
    block = styles[styles.index(selector):styles.index("}", styles.index(selector))]

    assert "color: var(--muted)" in block
    assert "pointer-events: none" in block
    assert "opacity: .65" in block


def test_workbench_clears_file_search_between_tasks() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "function resetFileSearch" in script
    assert "els.fileSearch.value = \"\"" in script
    assert "syncFileSearchClearButton(\"\")" in script
    select_start = script.index("async function selectTask")
    reset_search = script.index("resetFileSearch();", select_start)
    render_call = script.index("render();", select_start)

    assert reset_search < render_call


def test_workbench_ignores_stale_task_selection_responses() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    select_start = script.index("async function selectTask")
    fetch_snapshot = script.index("const snapshot = await api", select_start)
    stale_guard = script.index("if (state.taskId !== taskId) return", fetch_snapshot)
    set_snapshot = script.index("state.snapshot = snapshot", fetch_snapshot)
    subscribe_call = script.index("subscribe(taskId)", fetch_snapshot)

    assert fetch_snapshot < stale_guard < set_snapshot < subscribe_call


def test_workbench_restores_previous_task_when_task_selection_fails() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    select_start = script.index("async function selectTask")
    previous_task = script.index("const previousTaskId = state.taskId", select_start)
    previous_snapshot = script.index("const previousSnapshot = state.snapshot", previous_task)
    catch_block = script.index("catch (error)", previous_snapshot)
    restore_task = script.index("state.taskId = previousTaskId", catch_block)
    restore_snapshot = script.index("state.snapshot = previousSnapshot", catch_block)
    rerender = script.index("render();", catch_block)
    rethrow = script.index("throw error", catch_block)

    assert previous_task < previous_snapshot < catch_block
    assert restore_task < restore_snapshot < rerender < rethrow


def test_workbench_successful_task_selection_clears_stale_errors() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    select_start = script.index("async function selectTask")
    set_snapshot = script.index("state.snapshot = snapshot", select_start)
    clear_error = script.index('showError("");', set_snapshot)
    render_call = script.index("render();", set_snapshot)

    assert set_snapshot < clear_error < render_call


def test_workbench_ignores_stale_manual_refresh_responses() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    refresh_start = script.index("async function refresh")
    task_capture = script.index("const taskId = state.taskId", refresh_start)
    fetch_snapshot = script.index("const snapshot = await api", task_capture)
    stale_guard = script.index("if (state.taskId !== taskId) return", fetch_snapshot)
    set_snapshot = script.index("state.snapshot = snapshot", fetch_snapshot)
    render_call = script.index("render();", set_snapshot)

    assert task_capture < fetch_snapshot < stale_guard < set_snapshot < render_call


def test_workbench_refresh_buttons_show_errors_in_ui() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert 'els.refreshButton.addEventListener("click", runRefresh)' in script
    assert 'els.filesRefreshButton.addEventListener("click", runRefresh)' in script
    assert "async function runRefresh" in script
    run_refresh = script.index("async function runRefresh")
    clear_error = script.index('showError("")', run_refresh)
    refresh_call = script.index("await refresh()", run_refresh)
    catch_block = script.index("catch (error)", refresh_call)
    show_error = script.index("showError(error.message)", catch_block)

    assert clear_error < refresh_call < catch_block < show_error


def test_workbench_copy_buttons_show_success_and_failure_feedback() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    click_listener = script.index('document.addEventListener("click", async (event)')
    copy_button = script.index('event.target.closest("[data-copy]")', click_listener)
    clipboard_write = script.index("await navigator.clipboard.writeText", copy_button)
    success_feedback = script.index('showCopyFeedback(button, "Copied")', clipboard_write)
    catch_block = script.index("catch (error)", success_feedback)
    failure_error = script.index('showError("Unable to copy to clipboard.")', catch_block)
    helper = script.index("function showCopyFeedback")
    save_label = script.index("button.dataset.copyOriginalLabel || button.textContent", helper)
    clear_timer = script.index("window.clearTimeout(button.copyFeedbackTimer)", save_label)
    set_label = script.index("button.textContent = label", clear_timer)
    set_timer = script.index("button.copyFeedbackTimer = window.setTimeout", set_label)
    restore_label = script.index("button.textContent = button.dataset.copyOriginalLabel", set_timer)
    cleanup_label = script.index("delete button.dataset.copyOriginalLabel", restore_label)
    cleanup_timer = script.index("button.copyFeedbackTimer = null", cleanup_label)

    assert click_listener < copy_button < clipboard_write < success_feedback < catch_block < failure_error
    assert helper < save_label < clear_timer < set_label < set_timer < restore_label < cleanup_label < cleanup_timer


def test_workbench_copy_buttons_can_copy_dynamic_values() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    click_listener = script.index('document.addEventListener("click", async (event)')
    copy_button = script.index('event.target.closest("[data-copy]")', click_listener)
    value_lookup = script.index("button.dataset.copyValue", copy_button)
    target_lookup = script.index("document.getElementById(button.dataset.copy)", value_lookup)
    clipboard_write = script.index("await navigator.clipboard.writeText(copyValue)", target_lookup)
    compute_cards = script.index("function renderComputeCards")
    path_copy_check = script.index("isCopyableComputeFact(label, value)", compute_cards)
    copy_attr = script.index("data-copy-value", path_copy_check)
    helper = script.index("function isCopyableComputeFact")
    run_dir = script.index('"Run Dir"', helper)
    result_root = script.index('"Result Root"', run_dir)

    assert copy_button < value_lookup < target_lookup < clipboard_write
    assert compute_cards < path_copy_check < copy_attr < helper < run_dir < result_root


def test_workbench_formats_fastapi_validation_errors_for_users() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    api_start = script.index("async function api")
    parser_call = script.index("formatApiErrorDetail(payload.detail || detail)", api_start)
    throw_error = script.index("throw new Error(detail)", parser_call)
    formatter = script.index("function formatApiErrorDetail")
    array_case = script.index("Array.isArray(detail)", formatter)
    object_case = script.index('typeof detail === "object"', array_case)
    location_case = script.index("detail.loc.join", object_case)

    assert api_start < parser_call < throw_error
    assert formatter < array_case < object_case < location_case


def test_workbench_supports_browser_back_forward_navigation() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert 'window.addEventListener("popstate", handleLocationChange)' in script
    assert "async function handleLocationChange" in script
    assert 'location.pathname.match(/\\/tasks\\/([^/]+)/)' in script
    assert "await selectTask(match[1])" in script
    assert "showNewTask();" in script


def test_workbench_new_task_clears_active_history_row() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    show_new = script.index("function showNewTask")
    show_new_end = script.index("async function startTask", show_new)
    show_new_body = script[show_new:show_new_end]
    clear_task = show_new_body.index("state.taskId = null")
    refresh_rows = show_new_body.index("renderTasks(state.tasks)")
    render_main = show_new_body.index("render();")

    assert clear_task < refresh_rows < render_main


def test_workbench_recovers_url_when_route_task_selection_fails() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    route_start = script.index("async function handleLocationChange")
    route_catch = script.index("catch (error)", route_start)
    restore_route = script.index("restoreRouteAfterNavigationFailure();", route_catch)
    show_error = script.index("showError(error.message)", route_catch)
    helper_start = script.index("function restoreRouteAfterNavigationFailure")
    visible_path = script.index('state.taskId ? `/tasks/${state.taskId}` : "/"', helper_start)
    replace_state = script.index("history.replaceState(null, \"\", visiblePath)", helper_start)

    assert restore_route < show_error
    assert visible_path < replace_state


def test_workbench_styles_task_status_pills() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert "function statusClass" in script
    assert 'class="pill ${statusClass(task.status)}"' in script
    assert "els.statusPill.className" in script
    assert ".pill.status-completed" in styles
    assert ".pill.status-failed" in styles
    assert ".pill.status-running" in styles
    assert ".pill.status-needs-user-input" in styles


def test_workbench_final_messages_are_visually_distinct() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    render_messages = script.index("function renderMessages")
    final_class = script.index('message.final ? "final" : ""', render_messages)
    final_speaker = script.index('message.final ? "Final"', final_class)

    assert render_messages < final_class < final_speaker
    assert ".message.final" in styles


def test_workbench_status_pills_use_human_readable_labels() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function render()")
    topbar_label = script.index("els.statusPill.textContent = formatStatusLabel(status)", render_start)
    compute_label = script.index("els.computeState.textContent = formatStatusLabel(status)", topbar_label)
    task_rows = script.index("function renderTasks")
    row_label = script.index("formatStatusLabel(task.status)", task_rows)
    helper = script.index("function formatStatusLabel")
    needs_input = script.index('status === "needs_user_input"', helper)
    needs_return = script.index('return "Needs input"', needs_input)
    humanize = script.index("humanizeToolName(status)", needs_return)

    assert render_start < topbar_label < compute_label
    assert task_rows < row_label
    assert helper < needs_input < needs_return < humanize


def test_workbench_result_summary_renders_kpis_from_summary_json() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert "function renderSummaryKpis" in script
    assert 'id="summaryKpiSource"' in script
    assert "summary.json" in script
    assert "n_cells" in script
    assert "n_clusters" in script
    assert "Metrics from" in script
    assert ".kpi-grid" in styles
    assert ".kpi-card" in styles
    assert ".kpi-source" in styles


def test_workbench_result_summary_renders_kpis_without_final_answer() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    summary_start = script.index("function renderResultSummary")
    has_summary = script.index("const summaryFile = findSummaryFile(snapshot)", summary_start)
    empty_check = script.index("if (!finalAnswer && !summaryFile)", has_summary)
    summary_shell = script.index('id="summaryKpis"', empty_check)
    final_pre = script.index('id="resultFinalAnswerText"', summary_shell)
    render_kpis = script.index("renderSummaryKpis(snapshot)", final_pre)

    assert summary_start < has_summary < empty_check < summary_shell < final_pre < render_kpis


def test_workbench_result_summary_omits_empty_kpi_shell_without_summary_json() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    result_summary = script.index("function renderResultSummary")
    summary_markup = script.index("const summaryMarkup = summaryFile ?", result_summary)
    kpi_source = script.index('id="summaryKpiSource"', summary_markup)
    no_summary = script.index(': ""', kpi_source)
    final_markup = script.index("const finalAnswerMarkup = finalAnswer ?", no_summary)
    inner_html = script.index("els.resultSummary.innerHTML = `${summaryMarkup}${finalAnswerMarkup}`", final_markup)

    assert result_summary < summary_markup < kpi_source < no_summary < final_markup < inner_html


def test_workbench_result_summary_can_copy_final_answer() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    result_summary = script.index("function renderResultSummary")
    success_shell = script.index('id="summaryKpis"', result_summary)
    copy_button = script.index('data-copy="resultFinalAnswerText"', success_shell)
    final_pre = script.index('id="resultFinalAnswerText"', copy_button)
    failed_case = script.index('if (snapshot?.status === "failed")', result_summary)
    failed_copy = script.index('data-copy="resultFinalAnswerText"', failed_case)

    assert result_summary < success_shell < copy_button < final_pre
    assert failed_case < failed_copy
    assert ".summary-actions" in styles


def test_workbench_result_summary_accepts_analysis_summary_json() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "function findSummaryFile" in script
    summary_finder = script.index("function findSummaryFile")
    preferred = script.index('file.name === "summary.json"', summary_finder)
    fallback = script.index('file.name.includes("summary") && file.name.endsWith(".json")', preferred)
    render_kpis = script.index("function renderSummaryKpis")
    helper_call = script.index("findSummaryFile(", render_kpis)

    assert preferred < fallback
    assert render_kpis < helper_call


def test_workbench_result_summary_reads_common_summary_schema_variants() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "function firstPresent" in script
    kpi_start = script.index("const cards = [")
    assert 'firstPresent(payload, ["n_cells", "filtered_cells", "cells_final"])' in script[kpi_start:]
    assert 'firstPresent(payload, ["n_hvg", "n_hvgs", "hvg_count"])' in script[kpi_start:]
    assert 'firstPresent(payload, ["pca_n_comps", "n_pcs"])' in script[kpi_start:]
    assert 'firstPresent(qcMetrics, ["median_n_genes", "n_genes_median", "n_genes_by_counts_median"])' in script[kpi_start:]


def test_workbench_result_summary_kpis_show_fetch_errors() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    kpis_start = script.index("function renderSummaryKpis")
    fetch_call = script.index("fetch(`/api/tasks/${taskId}/files/content", kpis_start)
    ok_check = script.index("if (!response.ok) throw new Error", fetch_call)
    json_read = script.index("return response.json()", ok_check)
    catch_block = script.index(".catch(() =>", json_read)
    clear_target = script.index('target.innerHTML = ""', catch_block)
    unavailable_source = script.index('source.textContent = "Metrics unavailable"', clear_target)

    assert fetch_call < ok_check < json_read < catch_block < clear_target < unavailable_source


def test_workbench_result_summary_kpis_ignore_stale_fetch_errors() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    kpis_start = script.index("function renderSummaryKpis")
    catch_block = script.index(".catch(() =>", kpis_start)
    stale_guard = script.index("if (state.taskId !== taskId) return", catch_block)
    clear_target = script.index('target.innerHTML = ""', stale_guard)

    assert catch_block < stale_guard < clear_target


def test_workbench_result_cards_prioritize_primary_outputs() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function renderResultVisuals")
    assert "function sortResultCards" in script
    assert "function resultCardRank" in script
    image_filter = script.index(".filter((file) => file.kind === \"image\")", render_start)
    summary_filter = script.index(".filter((file) => file.name.endsWith(\".json\"))", image_filter)
    sorted_summaries = script.index("const sortedSummaries = summaries.sort(sortResultCards)", summary_filter)
    merged_cards = script.index("const resultCards = [...images, ...sortedSummaries.slice(0, 2)].sort(sortResultCards)", sorted_summaries)
    render_loop = script.index("for (const file of resultCards)", merged_cards)
    summary_branch = script.index("if (file.name.endsWith(\".json\"))", render_loop)
    image_branch = script.index("if (file.kind === \"image\")", summary_branch)
    rank_start = script.index("function resultCardRank")
    umap_rank = script.index('name.includes("umap")', rank_start)
    qc_rank = script.index('name.includes("qc")', umap_rank)
    pca_rank = script.index('name.includes("pca")', qc_rank)
    cluster_rank = script.index('name.includes("cluster")', pca_rank)
    summary_rank = script.index('name === "summary.json"', cluster_rank)

    assert image_filter < summary_filter < sorted_summaries < merged_cards < render_loop < summary_branch < image_branch
    assert umap_rank < qc_rank < pca_rank < cluster_rank < summary_rank


def test_workbench_result_visuals_distinguish_non_previewable_outputs() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function renderResultVisuals")
    files = script.index("const files = flattenFiles", render_start)
    no_cards = script.index("if (!resultCards.length)", files)
    message = script.index('files.length ? "No previewable result cards." : "No result files yet."', no_cards)
    empty_set = script.index("els.resultVisuals.innerHTML", message)

    assert render_start < files < no_cards < message < empty_set


def test_workbench_failed_result_visuals_do_not_repeat_empty_output_state() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function renderResultVisuals")
    clear_visuals = script.index('els.resultVisuals.innerHTML = ""', render_start)
    failed_empty = script.index('if (snapshot?.status === "failed" && !files.length)', clear_visuals)
    early_return = script.index("return", failed_empty)
    generic_empty = script.index("if (!resultCards.length)", early_return)

    assert render_start < clear_visuals < failed_empty < early_return < generic_empty


def test_workbench_result_tab_exposes_complete_output_file_list() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert 'id="resultFileList"' in markup
    render_start = script.index("function render()")
    visuals = script.index("renderResultVisuals(snapshot)", render_start)
    file_list = script.index("renderResultFileList(snapshot)", visuals)
    helper = script.index("function renderResultFileList")
    flatten = script.index("flattenFiles(snapshot?.resultFiles || [])", helper)
    row_class = script.index("output-file-row", flatten)
    open_call = script.index("openFile(file)", row_class)

    assert visuals < file_list < helper < flatten < row_class < open_call
    assert ".output-file-list" in styles
    assert ".output-file-row" in styles


def test_workbench_result_file_list_shows_relative_path_kind_and_size() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    helper = script.index("function renderResultFileList")
    kind = script.index("formatFileKindLabel(file)", helper)
    display_path = script.index("formatFileDisplayPath(file)", kind)
    size = script.index("formatBytes(file.size || 0)", display_path)
    keyboard = script.index("enableKeyboardActivation(row)", size)

    assert helper < kind < display_path < size < keyboard


def test_workbench_file_rows_expose_full_paths_on_hover() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_file = script.index("function renderFileNode")
    tree_title = script.index("row.title = formatFileDisplayPath(node)", render_file)
    tree_click = script.index('row.addEventListener("click", () => openFile(node))', tree_title)
    render_visuals = script.index("function renderResultVisuals")
    summary_title = script.index("card.title = formatFileDisplayPath(summary)", render_visuals)
    image_title = script.index("card.title = formatFileDisplayPath(image)", summary_title)
    render_list = script.index("function renderResultFileList")
    list_title = script.index("row.title = `${formatFileDisplayPath(file)} · ${formatBytes(file.size || 0)}`", render_list)
    list_click = script.index('row.addEventListener("click", () => openFile(file))', list_title)

    assert render_file < tree_title < tree_click
    assert render_visuals < summary_title < image_title < render_list
    assert render_list < list_title < list_click


def test_workbench_result_file_list_marks_selected_file() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    helper = script.index("function renderResultFileList")
    selected_class = script.index('${file.path === state.selectedFilePath ? "selected" : ""}', helper)
    data_path = script.index("row.dataset.filePath = file.path", selected_class)
    marker = script.index("function markSelectedFileRow")
    combined_selector = script.index('".file-node[data-file-path], .output-file-row[data-file-path]"', marker)
    toggle = script.index("row.classList.toggle", combined_selector)

    assert helper < selected_class < data_path
    assert marker < combined_selector < toggle
    assert ".output-file-row.selected" in styles


def test_workbench_result_json_cards_use_compact_preview() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function renderResultVisuals")
    summary_fetch = script.index("fetch(`/api/tasks/${summaryTaskId}/files/content", render_start)
    preview_set = script.index("card.querySelector(\"pre\").textContent = formatSummaryCardPreview(text)", summary_fetch)
    preview_helper = script.index("function formatSummaryCardPreview")
    json_parse = script.index("JSON.parse(text)", preview_helper)
    limiter = script.index("return limitPreview(pretty, 1200)", json_parse)
    limit_helper = script.index("function limitPreview")
    full_preview = script.index("return text", limit_helper)

    assert summary_fetch < preview_set < preview_helper
    assert preview_helper < json_parse < limiter < limit_helper < full_preview


def test_workbench_ignores_stale_summary_card_preview_responses() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "const summaryTaskId = state.taskId" in script
    assert "const summaryPath = summary.path" in script
    assert "if (state.taskId !== summaryTaskId) return" in script
    assert "if (!flattenFiles(state.snapshot?.resultFiles || []).some((file) => file.path === summaryPath)) return" in script


def test_workbench_summary_card_preview_shows_fetch_errors() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function renderResultVisuals")
    summary_fetch = script.index("fetch(`/api/tasks/${summaryTaskId}/files/content", render_start)
    ok_check = script.index("if (!r.ok) throw new Error", summary_fetch)
    text_read = script.index("return r.text()", ok_check)
    catch_block = script.index(".catch(() =>", text_read)
    error_text = script.index("Unable to load preview.", catch_block)

    assert summary_fetch < ok_check < text_read < catch_block < error_text


def test_workbench_file_preview_shows_result_relative_path() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    open_file = script.index("async function openFile")
    meta_set = script.index("els.fileMeta.textContent = formatFileDisplayPath(node)", open_file)
    download_set = script.index("els.downloadLink.href", meta_set)
    helper = script.index("function formatFileDisplayPath")
    result_root = script.index("state.snapshot?.compute?.resultRoot", helper)
    prefix_check = script.index("normalized.startsWith(`${normalizedRoot}/`)", result_root)
    relative_return = script.index("return normalized.slice(normalizedRoot.length + 1)", prefix_check)

    assert open_file < meta_set < download_set
    assert helper < result_root < prefix_check < relative_return


def test_workbench_file_panel_can_copy_selected_file_path() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    file_pane = markup.index('id="filePane"')
    copy_button = markup.index('id="copyFilePathButton"', file_pane)
    download_link = markup.index('id="downloadLink"', copy_button)
    els_list = script.index('"copyFilePathButton"')
    open_file = script.index("async function openFile")
    copy_value = script.index("els.copyFilePathButton.dataset.copyValue = node.path", open_file)
    copy_enabled = script.index('els.copyFilePathButton.classList.remove("disabled")', copy_value)
    reset = script.index("function resetFilePreview")
    copy_disabled = script.index('els.copyFilePathButton.classList.add("disabled")', reset)
    copy_cleared = script.index("delete els.copyFilePathButton.dataset.copyValue", copy_disabled)

    assert file_pane < copy_button < download_link
    assert els_list < open_file < copy_value < copy_enabled
    assert reset < copy_disabled < copy_cleared


def test_workbench_file_panel_guidance_mentions_output_files_and_downloads() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    guidance = "Open an output file to preview text/images or download binary results."
    reset = script.index("function resetFilePreview")
    reset_guidance = script.index(f'els.fileMeta.textContent = "{guidance}"', reset)

    assert guidance in markup
    assert reset < reset_guidance


def test_workbench_merges_live_snapshot_into_task_list() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "function mergeTaskSnapshot" in script
    assert "mergeTaskSnapshot(state.snapshot)" in script
    assert "renderTasks(state.tasks)" in script


def test_workbench_live_task_list_keeps_failure_reason_from_snapshot() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    merge_start = script.index("function mergeTaskSnapshot")
    item_start = script.index("const item = {", merge_start)
    error_field = script.index("error: snapshot.error", item_start)
    status_field = script.index("status: snapshot.status", item_start)

    assert status_field < error_field


def test_workbench_renders_semantic_plan_and_activity_as_separate_surfaces() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function renderPlan")
    semantic_plan = script.index("snapshot?.planState?.steps", render_start)
    activity_start = script.index("function renderActivity")
    activity_source = script.index("snapshot?.activity", activity_start)
    recent_activity = script.index("slice(-8)", activity_source)

    assert 'id="planList"' in markup
    assert 'id="activityList"' in markup
    assert render_start < semantic_plan < activity_start < activity_source < recent_activity


def test_workbench_does_not_repeat_tool_trace_rows_inside_the_conversation() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_messages = script[script.index("function renderMessages"):script.index("function displayMessageContent")]

    assert "traceBlock(" not in render_messages


def test_workbench_trace_tab_contains_compact_plan_when_right_rail_is_hidden() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    trace_pane = markup.index('id="tracePane"')
    trace_plan = markup.index('id="tracePlanList"', trace_pane)
    right_plan = markup.index('id="planList"', trace_plan)
    els_list = script.index('"tracePlanCount", "tracePlanList"')
    render_plan = script.index("function renderPlan")
    trace_count = script.index("els.tracePlanCount.textContent = String(semanticPlan.length)", render_plan)
    trace_render = script.index("renderPlanItems(els.tracePlanList, semanticPlan, snapshot)", trace_count)

    assert trace_pane < trace_plan < right_plan
    assert els_list < render_plan < trace_count < trace_render
    assert ".trace-plan" in styles


def test_workbench_plan_titles_hide_internal_call_prefix() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function renderPlan")
    title_call = script.index("formatPlanTitle(item.title)", render_start)
    step_title = script.index("<strong>", title_call)
    helper = script.index("function formatPlanTitle")
    call_prefix = script.index('title.startsWith("Call ")', helper)
    strip_prefix = script.index('title.slice("Call ".length)', call_prefix)

    assert title_call < step_title
    assert helper < call_prefix < strip_prefix


def test_workbench_plan_titles_humanize_tool_names() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    formatter = script.index("function formatPlanTitle")
    helper_call = script.index("humanizeToolName", formatter)
    helper = script.index("function humanizeToolName")
    underscores = script.index('replaceAll("_", " ")', helper)
    capitalize = script.index("charAt(0).toUpperCase()", underscores)

    assert formatter < helper_call < helper
    assert helper < underscores < capitalize


def test_workbench_plan_marks_intermediate_failures_recovered_for_completed_runs() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    render_start = script.index("function renderPlan")
    status_const = script.index("const displayStatus = formatPlanStatus(item, snapshot)", render_start)
    check_class = script.index("check ${escapeHtml(displayStatus)}", status_const)
    status_text = script.index("${escapeHtml(displayStatus)}", check_class)
    helper = script.index("function formatPlanStatus")
    completed_case = script.index('snapshot?.status === "completed"', helper)
    failed_case = script.index('item.status === "failed"', completed_case)
    recovered = script.index('return "recovered"', failed_case)

    assert status_const < check_class < status_text
    assert helper < completed_case < failed_case < recovered
    assert ".check.recovered" in styles


def test_workbench_status_dots_do_not_render_pending_as_completed() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    trace_dot = styles.index(".status-dot {")
    trace_pending = styles.index("background: #cbd5e1", trace_dot)
    trace_completed = styles.index(".status-dot.completed", trace_pending)
    plan_dot = styles.index(".check {")
    plan_pending = styles.index("background: #cbd5e1", plan_dot)
    plan_completed = styles.index(".check.completed", plan_pending)

    assert trace_dot < trace_pending < trace_completed
    assert plan_dot < plan_pending < plan_completed


def test_workbench_plan_step_titles_do_not_push_status_offscreen() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    step_line = styles.index(".step-line {")
    title = styles.index(".step-line strong", step_line)
    title_block = styles[title:styles.index("}", title)]

    assert "min-width: 0" in title_block
    assert "overflow: hidden" in title_block
    assert "text-overflow: ellipsis" in title_block
    assert "white-space: nowrap" in title_block


def test_workbench_task_rows_show_run_id_for_identical_titles() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "formatTaskMeta(task)" in script
    assert "function formatTaskMeta" in script
    assert "task.id" in script
    assert "updated" in script


def test_workbench_task_rows_expose_full_title_on_hover() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    set_title = script.index("row.title = formatTaskTooltip(task)", render_tasks)
    set_html = script.index("row.innerHTML", set_title)
    helper = script.index("function formatTaskTooltip")
    full_title = script.index("task.title", helper)
    meta = script.index("formatTaskMeta(task)", helper)

    assert set_title < set_html
    assert helper < full_title < meta


def test_workbench_task_rows_search_and_hover_use_full_task_text() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    search_start = script.index("function taskSearchText")
    search_task = script.index("task.task", search_start)
    search_meta = script.index("formatTaskMeta(task)", search_task)
    tooltip_start = script.index("function formatTaskTooltip")
    tooltip_task = script.index("task.task || task.title", tooltip_start)
    tooltip_meta = script.index("formatTaskMeta(task)", tooltip_task)

    assert search_task < search_meta
    assert tooltip_start < tooltip_task < tooltip_meta


def test_workbench_task_rows_render_distinguishable_display_titles() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    display_call = script.index("formatTaskDisplayTitle(task)", render_tasks)
    helper = script.index("function formatTaskDisplayTitle")
    full_task = script.index("task.task || task.title", helper)
    basename = script.index("basename(task.dataPath)", full_task)

    assert display_call < script.index("row.innerHTML", render_tasks)
    assert helper < full_task < basename
    assert ".row-title strong" in styles
    assert ".row-title span" in styles


def test_workbench_topbar_uses_readable_task_display_title() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_start = script.index("function render")
    display_title = script.index("formatTaskDisplayTitle(snapshot)", render_start)
    topbar_title = script.index("els.taskTitle.textContent = displayTitle.title", display_title)
    full_title = script.index('els.taskTitle.title = snapshot?.task || snapshot?.title || "New session"', topbar_title)

    assert render_start < display_title < topbar_title < full_title


def test_workbench_task_rows_keep_visible_meta_short() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    row_sub = script.index('class="row-sub"', render_tasks)
    display_meta_call = script.index("formatTaskDisplayMeta(task, query)", row_sub)
    helper = script.index("function formatTaskDisplayMeta")
    failed_case = script.index('task.status === "failed"', helper)
    updated_case = script.index("updated", failed_case)
    full_meta_helper = script.index("function formatTaskMeta")
    full_data_path = script.index("task.dataPath", full_meta_helper)

    assert row_sub < display_meta_call
    assert helper < failed_case < updated_case
    assert full_meta_helper < full_data_path
    assert ".row-sub" in styles
    assert "text-overflow: ellipsis" in styles[styles.index(".row-sub"):]


def test_workbench_task_rows_show_run_duration_in_visible_and_full_meta() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    visible_meta = script.index("function formatTaskDisplayMeta")
    visible_duration = script.index("const duration = formatRunDuration(task)", visible_meta)
    visible_return = script.index("${duration}", visible_duration)
    full_meta = script.index("function formatTaskMeta")
    full_duration = script.index("const duration = formatRunDuration(task)", full_meta)
    full_return = script.index("${duration}", full_duration)

    assert visible_meta < visible_duration < visible_return
    assert full_meta < full_duration < full_return


def test_workbench_task_rows_show_run_activity_counts_in_meta() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    helper = script.index("function formatTaskActivityCounts")
    messages = script.index("task.messageCount", helper)
    traces = script.index("task.traceCount", messages)
    files = script.index("task.resultFileCount", traces)
    visible_meta = script.index("function formatTaskDisplayMeta")
    visible_counts = script.index("const activity = formatTaskActivityCounts(task)", visible_meta)
    visible_return = script.index("${activity}", visible_counts)
    full_meta = script.index("function formatTaskMeta")
    full_counts = script.index("const activity = formatTaskActivityCounts(task)", full_meta)
    full_return = script.index("${activity}", full_counts)

    assert helper < messages < traces < files
    assert visible_meta < visible_counts < visible_return
    assert full_meta < full_counts < full_return


def test_workbench_task_rows_expose_full_context_to_assistive_tech() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    set_title = script.index("row.title = formatTaskTooltip(task)", render_tasks)
    set_label = script.index('row.setAttribute("aria-label", formatTaskAriaLabel(task))', set_title)
    set_html = script.index("row.innerHTML", set_label)
    helper = script.index("function formatTaskAriaLabel")
    tooltip = script.index("formatTaskTooltip(task)", helper)
    replace_newline = script.index('.replaceAll("\\n", " · ")', tooltip)

    assert set_title < set_label < set_html
    assert helper < tooltip < replace_newline


def test_workbench_task_rows_offer_direct_rerun_action() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    rerun_markup = script.index("data-rerun-task-id", render_tasks)
    click_start = script.index('row.addEventListener("click"', rerun_markup)
    rerun_guard = script.index("event.target.closest(\"[data-rerun-task-id]\")", click_start)
    helper_start = script.index("async function prepareRerunFromList")
    stop_propagation = script.index("event.stopPropagation();", helper_start)
    open_call = script.index("await openTaskFromList(task.id)", stop_propagation)
    focus_call = script.index("els.taskInput.focus()", open_call)

    assert render_tasks < rerun_markup < click_start < rerun_guard
    assert helper_start < stop_propagation < open_call < focus_call
    assert ".task-actions" in styles
    compact_actions = styles[styles.index(".task-log, .task-rerun {"):styles.index("}", styles.index(".task-log, .task-rerun {"))]
    assert "min-height: 24px" in compact_actions
    assert "padding: 2px 7px" in compact_actions
    assert "font-size: 11px" in compact_actions


def test_workbench_task_actions_wrap_without_forcing_sidebar_width() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    actions = styles[styles.index(".task-actions {"):styles.index("}", styles.index(".task-actions {"))]
    row_title = styles[styles.index(".row-title {"):styles.index("}", styles.index(".row-title {"))]
    row_sub_start = styles.index(".row-sub {", styles.index(".row-title span"))
    row_sub = styles[row_sub_start:styles.index("}", row_sub_start)]

    assert "min-width: 0" in actions
    assert "max-width: 100%" in actions
    assert "flex-wrap: wrap" in actions
    assert "min-width: 0" in row_title
    assert "grid-column: 1 / 3" in row_sub


def test_workbench_task_rows_keep_titles_above_actions_in_narrow_sidebar() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    task_row_start = styles.index(".task-row {")
    task_row = styles[task_row_start:styles.index("}", task_row_start)]
    task_actions_start = styles.index(".task-row .task-actions {")
    task_actions = styles[task_actions_start:styles.index("}", task_actions_start)]
    task_meta_start = styles.index(".task-row .row-sub {")
    task_meta = styles[task_meta_start:styles.index("}", task_meta_start)]

    assert "grid-template-columns: minmax(0, 1fr)" in task_row
    assert "justify-content: flex-start" in task_actions
    assert "grid-column: 1" in task_meta


def test_workbench_clicking_active_task_row_preserves_current_context() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    click_start = script.index('row.addEventListener("click"', render_tasks)
    guard = script.index("if (task.id === state.taskId) return", click_start)
    push_state = script.index("history.pushState", click_start)

    assert guard < push_state


def test_workbench_task_row_selection_failure_restores_visible_url() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    click_start = script.index('row.addEventListener("click"', render_tasks)
    push_state = script.index("history.pushState(null, \"\", `/tasks/${task.id}`)", click_start)
    open_call = script.index("openTaskFromList(task.id)", push_state)
    helper_start = script.index("async function openTaskFromList")
    select_call = script.index("await selectTask(taskId)", helper_start)
    catch_block = script.index("catch (error)", select_call)
    restore_route = script.index("restoreRouteAfterNavigationFailure();", catch_block)
    show_error = script.index("showError(error.message)", catch_block)

    assert push_state < open_call
    assert select_call < catch_block < restore_route < show_error


def test_task_list_orders_runs_by_latest_update(tmp_path: Path) -> None:
    config = _config(tmp_path)
    history_dir = config.runs_dir / "web_tasks"
    history_dir.mkdir(parents=True)
    (history_dir / "older_created_recently_updated.json").write_text(
        """{
          "id": "recent_update",
          "task": "Older created",
          "title": "Older created",
          "status": "completed",
          "created_at": "2026-01-01T00:00:00+00:00",
          "updated_at": "2026-01-03T00:00:00+00:00",
          "messages": [],
          "traces": [],
          "plan": [],
          "events": [],
          "compute": {}
        }""",
        encoding="utf-8",
    )
    (history_dir / "newer_created_stale_update.json").write_text(
        """{
          "id": "stale_update",
          "task": "Newer created",
          "title": "Newer created",
          "status": "completed",
          "created_at": "2026-01-02T00:00:00+00:00",
          "updated_at": "2026-01-02T00:00:00+00:00",
          "messages": [],
          "traces": [],
          "plan": [],
          "events": [],
          "compute": {}
        }""",
        encoding="utf-8",
    )

    app = create_app(config=config, stream_runner=_fake_stream)
    client = TestClient(app)

    tasks = client.get("/api/tasks").json()["tasks"]

    assert [task["id"] for task in tasks] == ["recent_update", "stale_update"]


def test_task_list_includes_failure_reason_for_failed_runs(tmp_path: Path) -> None:
    config = _config(tmp_path)
    history_dir = config.runs_dir / "web_tasks"
    history_dir.mkdir(parents=True)
    (history_dir / "failed_model_loading.json").write_text(
        """{
          "id": "failed_model_loading",
          "task": "Failed demo",
          "title": "Failed demo",
          "status": "failed",
          "created_at": "2026-01-01T00:00:00+00:00",
          "updated_at": "2026-01-01T00:00:01+00:00",
          "error": "Loading model",
          "messages": [],
          "traces": [],
          "plan": [],
          "events": [],
          "compute": {}
        }""",
        encoding="utf-8",
    )

    app = create_app(config=config, stream_runner=_fake_stream)
    client = TestClient(app)

    tasks = client.get("/api/tasks").json()["tasks"]

    assert tasks[0]["id"] == "failed_model_loading"
    assert tasks[0]["error"] == "Loading model"


def test_task_list_includes_data_path_for_run_context(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)

    client.post("/api/tasks", json={"task": "Run demo", "data_path": "mas_2/data/bmmc_b_cell.h5ad"})
    tasks = client.get("/api/tasks").json()["tasks"]

    assert tasks[0]["dataPath"] == "mas_2/data/bmmc_b_cell.h5ad"


def test_task_list_includes_run_activity_counts(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "completed")

    tasks = client.get("/api/tasks").json()["tasks"]

    assert tasks[0]["messageCount"] == 3
    assert tasks[0]["traceCount"] == 1
    assert tasks[0]["resultFileCount"] == 2


def test_task_list_includes_log_path_for_direct_log_access(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "completed")

    tasks = client.get("/api/tasks").json()["tasks"]

    assert tasks[0]["logPath"].endswith("/run.log")


def test_task_list_includes_result_file_names_for_search_context(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    task_id = client.post("/api/tasks", json={"task": "Run demo"}).json()["taskId"]
    _wait_for_status(client, task_id, "completed")

    tasks = client.get("/api/tasks").json()["tasks"]

    assert tasks[0]["resultFileNames"] == ["analysis_summary.json", "umap_leiden.png"]


def test_task_list_includes_full_task_for_hover_context(tmp_path: Path) -> None:
    app = create_app(config=_config(tmp_path), stream_runner=_fake_stream)
    client = TestClient(app)
    full_task = "Run demo\nwith full instructions that should not be truncated in task list hover text"

    client.post("/api/tasks", json={"task": full_task})
    tasks = client.get("/api/tasks").json()["tasks"]

    assert tasks[0]["task"] == full_task


def test_task_list_summarizes_provider_error_payloads(tmp_path: Path) -> None:
    config = _config(tmp_path)
    history_dir = config.runs_dir / "web_tasks"
    history_dir.mkdir(parents=True)
    (history_dir / "failed_provider_payload.json").write_text(
        """{
          "id": "failed_provider_payload",
          "task": "Failed demo",
          "title": "Failed demo",
          "status": "failed",
          "created_at": "2026-01-01T00:00:00+00:00",
          "updated_at": "2026-01-01T00:00:01+00:00",
          "error": "Error code: 503 - {'error': {'message': 'Loading model', 'type': 'unavailable_error', 'code': 503}}",
          "messages": [],
          "traces": [],
          "plan": [],
          "events": [],
          "compute": {}
        }""",
        encoding="utf-8",
    )

    app = create_app(config=config, stream_runner=_fake_stream)
    client = TestClient(app)

    tasks = client.get("/api/tasks").json()["tasks"]

    assert tasks[0]["error"] == "Loading model"


def test_workbench_live_task_merge_orders_runs_by_latest_update() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "new Date(right.updatedAt) - new Date(left.updatedAt)" in script


def test_workbench_task_search_matches_id_status_and_metadata() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "function taskSearchText" in script
    assert "taskSearchText(task).includes(query)" in script
    assert "task.status" in script
    assert "task.dataPath" in script
    assert "formatTaskMeta(task)" in script


def test_workbench_task_search_matches_visible_status_labels() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    search_start = script.index("function taskSearchText")
    raw_status = script.index("task.status", search_start)
    visible_status = script.index("formatStatusLabel(task.status)", raw_status)
    meta = script.index("formatTaskMeta(task)", visible_status)

    assert search_start < raw_status < visible_status < meta


def test_workbench_task_search_and_tooltips_include_result_file_names() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    formatter = script.index("function formatTaskMeta")
    outputs = script.index("const outputs = formatTaskOutputNames(task)", formatter)
    output_return = script.index("${outputs}", outputs)
    helper = script.index("function formatTaskOutputNames")
    file_names = script.index("task.resultFileNames", helper)
    join_names = script.index('.join(", ")', file_names)
    tooltip = script.index("function formatTaskTooltip")
    tooltip_meta = script.index("formatTaskMeta(task)", tooltip)
    search = script.index("function taskSearchText")
    search_meta = script.index("formatTaskMeta(task)", search)

    assert formatter < outputs < output_return
    assert helper < file_names < join_names
    assert tooltip < tooltip_meta
    assert search < search_meta


def test_workbench_task_search_shows_matching_output_file_in_visible_meta() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    query = script.index("const query = els.taskSearch.value.trim().toLowerCase()", render_tasks)
    visible_meta = script.index("formatTaskDisplayMeta(task, query)", query)
    formatter = script.index("function formatTaskDisplayMeta")
    match = script.index("const outputMatch = formatTaskOutputMatch(task, query)", formatter)
    return_match = script.index("${outputMatch}", match)
    helper = script.index("function formatTaskOutputMatch")
    file_names = script.index("task.resultFileNames", helper)
    includes = script.index(".includes(normalizedQuery)", file_names)
    label = script.index('" · match: "', includes)

    assert render_tasks < query < visible_meta
    assert formatter < match < return_match
    assert helper < file_names < includes < label


def test_workbench_task_search_can_be_cleared_without_losing_focus() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    search_row = markup.index('class="task-search-row"')
    search_input = markup.index('id="taskSearch"', search_row)
    clear_button = markup.index('id="clearTaskSearchButton"', search_input)
    els_list = script.index('"taskSearch", "clearTaskSearchButton"')
    input_listener = script.index("syncTaskSearchClearButton(query)", script.index("function renderTasks"))
    click_listener = script.index('els.clearTaskSearchButton.addEventListener("click"', els_list)
    clear_value = script.index('els.taskSearch.value = ""', click_listener)
    rerender = script.index("renderTasks(state.tasks)", clear_value)
    focus = script.index("els.taskSearch.focus()", rerender)

    assert search_row < search_input < clear_button
    assert els_list < click_listener < clear_value < rerender < focus
    assert script.index("function renderTasks") < input_listener
    assert ".task-search-row" in styles


def test_workbench_task_rows_offer_direct_log_preview_when_available() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    log_button = script.index('data-open-task-log-id="${escapeHtml(task.id)}"', render_tasks)
    rerun_button = script.index('data-rerun-task-id="${escapeHtml(task.id)}"', log_button)
    handler = script.index("const logButton = row.querySelector(\"[data-open-task-log-id]\")", rerun_button)
    listener = script.index("prepareOpenLogFromList(event, task)", handler)
    row_click_guard = script.index('event.target.closest("[data-open-task-log-id]")', listener)
    helper = script.index("async function prepareOpenLogFromList")
    stop = script.index("event.stopPropagation()", helper)
    open_task = script.index("await openTaskFromList(task.id)", stop)
    open_log = script.index("await openLog()", open_task)

    assert render_tasks < log_button < rerun_button < handler < listener < row_click_guard
    assert helper < stop < open_task < open_log


def test_workbench_task_rows_offer_direct_trace_inspection_when_available() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_tasks = script.index("function renderTasks")
    trace_button = script.index('data-open-task-trace-id="${escapeHtml(task.id)}"', render_tasks)
    log_button = script.index('data-open-task-log-id="${escapeHtml(task.id)}"', trace_button)
    handler = script.index("const traceButton = row.querySelector(\"[data-open-task-trace-id]\")", log_button)
    listener = script.index("prepareOpenTraceFromList(event, task)", handler)
    row_click_guard = script.index('event.target.closest("[data-open-task-trace-id]")', listener)
    helper = script.index("async function prepareOpenTraceFromList")
    stop = script.index("event.stopPropagation()", helper)
    open_task = script.index("await openTaskFromList(task.id)", stop)
    activate_trace = script.index('activateTab("trace", { user: true })', open_task)

    assert render_tasks < trace_button < log_button < handler < listener < row_click_guard
    assert helper < stop < open_task < activate_trace


def test_workbench_task_search_trims_query_whitespace() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    renderer = script.index("function renderTasks")
    query = script.index("els.taskSearch.value.trim().toLowerCase()", renderer)
    filtered = script.index("taskSearchText(task).includes(query)", query)

    assert renderer < query < filtered


def test_workbench_task_rows_show_data_path_when_available() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    formatter = script.index("function formatTaskMeta")
    data_path = script.index("task.dataPath", formatter)
    failed_check = script.index('task.status === "failed"', data_path)
    failed_return = script.index("failed:", failed_check)
    updated = script.index("updated", failed_return)

    assert data_path < failed_check < failed_return < updated


def test_workbench_failed_task_rows_show_failure_reason() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    formatter = script.index("function formatTaskMeta")
    failed_check = script.index('task.status === "failed"', formatter)
    error_text = script.index("failed:", failed_check)
    error_search = script.index("task.error", script.index("function taskSearchText"))

    assert failed_check < error_text
    assert error_search < formatter


def test_workbench_task_tooltips_use_concise_failed_error_summary() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    formatter = script.index("function formatTaskMeta")
    failed_check = script.index('task.status === "failed"', formatter)
    failed_return = script.index("failed:", failed_check)
    concise_error = script.index("formatErrorSummary(task.error)", failed_return)
    tooltip = script.index("function formatTaskTooltip")
    tooltip_meta = script.index("formatTaskMeta(task)", tooltip)

    assert formatter < failed_check < failed_return < concise_error < tooltip < tooltip_meta


def test_workbench_file_search_matches_path_type_and_kind() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "function fileSearchText" in script
    assert "fileSearchText(node).includes(query)" in script
    assert "node.path" in script
    assert "node.type" in script
    assert "node.kind" in script


def test_workbench_file_search_can_be_cleared_without_losing_focus() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    search_row = markup.index('class="file-search-row"')
    search_input = markup.index('id="fileSearch"', search_row)
    clear_button = markup.index('id="clearFileSearchButton"', search_input)
    els_list = script.index('"fileSearch", "clearFileSearchButton"')
    input_listener = script.index("syncFileSearchClearButton(query)", script.index("function renderFiles"))
    click_listener = script.index('els.clearFileSearchButton.addEventListener("click"', els_list)
    clear_value = script.index('els.fileSearch.value = ""', click_listener)
    rerender = script.index("renderFiles(state.snapshot)", clear_value)
    focus = script.index("els.fileSearch.focus()", rerender)

    assert search_row < search_input < clear_button
    assert els_list < click_listener < clear_value < rerender < focus
    assert script.index("function renderFiles") < input_listener
    assert ".file-search-row" in styles


def test_workbench_demo_project_file_fills_data_path() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert 'data-demo-path="mas_2/data/bmmc_b_cell.h5ad"' in markup
    assert "document.querySelectorAll(\"[data-demo-path]\")" in script
    assert "els.dataPathInput.value = row.dataset.demoPath" in script
    assert ".file-row[data-demo-path]" in styles


def test_workbench_clickable_rows_are_keyboard_accessible() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert 'role="button"' in markup
    assert 'tabindex="0"' in markup
    assert "function isActivationKey" in script
    assert 'event.key === "Enter"' in script
    assert 'event.key === " "' in script
    assert 'row.setAttribute("role", "button")' in script
    assert "row.tabIndex = 0" in script
    assert 'row.addEventListener("keydown"' in script
    assert "event.preventDefault()" in script
    assert "row.click()" in script


def test_workbench_keyboard_activation_ignores_nested_control_events() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    helper = script.index("function enableKeyboardActivation")
    keydown = script.index('row.addEventListener("keydown"', helper)
    target_guard = script.index("if (event.target !== row) return", keydown)
    activation_guard = script.index("if (!isActivationKey(event)) return", target_guard)
    row_click = script.index("row.click()", activation_guard)

    assert helper < keydown < target_guard < activation_guard < row_click


def test_workbench_trace_files_and_result_cards_are_keyboard_accessible() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "function enableKeyboardActivation" in script
    render_traces = script.index("function renderTraces")
    trace_click = script.index('row.addEventListener("click"', render_traces)
    trace_keyboard = script.index("enableKeyboardActivation(row);", trace_click)
    render_file = script.index("function renderFileNode")
    file_click = script.index('row.addEventListener("click"', render_file)
    file_keyboard = script.index("enableKeyboardActivation(row);", file_click)
    render_visuals = script.index("function renderResultVisuals")
    first_card_click = script.index('card.addEventListener("click"', render_visuals)
    first_card_keyboard = script.index("enableKeyboardActivation(card);", first_card_click)
    second_card_click = script.index('card.addEventListener("click"', first_card_keyboard)
    second_card_keyboard = script.index("enableKeyboardActivation(card);", second_card_click)

    assert trace_click < trace_keyboard
    assert file_click < file_keyboard
    assert first_card_click < first_card_keyboard < second_card_click < second_card_keyboard


def test_workbench_trace_tab_contains_per_call_browser() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    trace_pane = markup.index('id="tracePane"')
    workbench = markup.index('id="traceWorkbench"', trace_pane)
    call_count = markup.index('id="traceListCount"', workbench)
    trace_list = markup.index('id="traceList"', call_count)
    inspector = markup.index('class="trace-inspector"', trace_list)
    trace_input = markup.index('id="traceInput"', inspector)
    trace_output = markup.index('id="traceOutput"', trace_input)

    assert trace_pane < workbench < call_count < trace_list < inspector < trace_input < trace_output
    assert '"traceListCount", "traceList"' in script
    assert ".trace-workbench" in styles
    assert ".trace-browser" in styles
    assert ".trace-list-row" in styles
    assert ".trace-inspector" in styles


def test_workbench_trace_inspector_prioritizes_input_and_output_before_plan() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")

    trace_pane = markup.index('id="tracePane"')
    inspector = markup.index('class="trace-inspector"', trace_pane)
    trace_input = markup.index('id="traceInput"', inspector)
    trace_output = markup.index('id="traceOutput"', trace_input)
    trace_plan = markup.index('id="tracePlanList"', trace_output)

    assert trace_pane < inspector < trace_input < trace_output < trace_plan


def test_workbench_trace_browser_renders_every_call_and_shared_selection() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_traces = script.index("function renderTraces")
    render_browser_call = script.index("renderTraceBrowser(snapshot)", render_traces)
    browser_helper = script.index("function renderTraceBrowser")
    traces = script.index("const traces = snapshot?.traces || []", browser_helper)
    count = script.index("els.traceListCount.textContent = String(traces.length)", traces)
    map_calls = script.index("traces.map((trace, index)", count)
    trace_meta = script.index("const displayStatus = formatTraceMeta(trace, snapshot)", map_calls)
    call_id = script.index('data-trace-id="${escapeHtml(trace.id)}"', trace_meta)
    tool_name = script.index("formatToolName(trace.toolName)", call_id)
    browser_click = script.index('els.traceList.querySelectorAll("[data-trace-id]")', tool_name)
    shared_selection = script.index("selectTrace(row.dataset.traceId)", browser_click)
    selector_helper = script.index("function selectTrace")
    select_id = script.index("state.selectedTraceId = traceId", selector_helper)
    mark_rows = script.index("markSelectedTraceRow()", select_id)
    render_detail = script.index("renderTraceDetail()", mark_rows)
    activate_tab = script.index('activateTab("trace", { user: true })', render_detail)

    assert render_traces < render_browser_call < browser_helper
    assert browser_helper < traces < count < map_calls < trace_meta < call_id < tool_name < browser_click < shared_selection
    assert selector_helper < select_id < mark_rows < render_detail < activate_tab


def test_workbench_trace_browser_stacks_on_mobile() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    workbench_start = styles.index(".trace-workbench {")
    workbench = styles[workbench_start:styles.index("}", workbench_start)]
    media_start = styles.index("@media (max-width: 980px)")
    media_block = styles[media_start:]

    assert "grid-template-columns: minmax(145px, .42fr) minmax(0, 1fr)" in workbench
    assert ".trace-workbench { grid-template-columns: minmax(0, 1fr); }" in media_block
    assert ".trace-list { max-height: 240px; }" in media_block


def test_workbench_trace_browser_reveals_selected_call_after_tab_becomes_visible() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    browser_helper = script.index("function renderTraceBrowser")
    render_reveal = script.index("revealSelectedTraceRow()", browser_helper)
    reveal_helper = script.index("function revealSelectedTraceRow")
    selected_row = script.index("row.dataset.traceId === state.selectedTraceId", reveal_helper)
    visible_height_guard = script.index("if (!selectedRow || !els.traceList.clientHeight) return", selected_row)
    scroll_position = script.index("els.traceList.scrollTop", visible_height_guard)
    activate_helper = script.index("function activateTab")
    pane_loop = script.index('document.querySelectorAll(".detail-pane").forEach', activate_helper)
    trace_guard = script.index('if (name === "trace")', pane_loop)
    activate_reveal = script.index("revealSelectedTraceRow()", trace_guard)

    assert browser_helper < render_reveal < reveal_helper < selected_row < visible_height_guard < scroll_position
    assert activate_helper < pane_loop < trace_guard < activate_reveal


def test_workbench_trace_detail_distinguishes_missing_output_from_empty_object() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    detail = script.index("function renderTraceDetail")
    waiting = script.index('"Waiting for tool output..."', detail)
    missing = script.index('"No output was received before the call ended."', waiting)
    output = script.index("els.traceOutput.textContent", missing)

    assert detail < waiting < missing < output


def test_workbench_trace_payload_expands_code_and_process_streams() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    detail = script.index("function renderTraceDetail")
    input_format = script.index('formatTracePayload(trace.input, ["code"])', detail)
    output_format = script.index('formatTracePayload(trace.output, ["stdout", "stderr"])', input_format)
    formatter = script.index("function formatTracePayload", output_format)
    raw_string = script.index('if (typeof value === "string") return value;', formatter)
    expanded_section = script.index('`${key}:\\n${text || "(empty)"}`', raw_string)

    assert detail < input_format < output_format < formatter < raw_string < expanded_section


def test_workbench_keyboard_focus_is_visible_on_interactive_rows_and_cards() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    focus_start = styles.index(".task-row:focus-visible")
    assert ".file-row[data-demo-path]:focus-visible" in styles[focus_start:]
    assert ".activity-row:focus-visible" in styles[focus_start:]
    assert ".file-node.is-file:focus-visible" in styles[focus_start:]
    assert ".plot-card:focus-visible" in styles[focus_start:]
    assert ".summary-card:focus-visible" in styles[focus_start:]
    assert "outline: 2px solid var(--blue)" in styles[focus_start:]
    assert "outline-offset: 2px" in styles[focus_start:]
    assert "background: #eef3ff" in styles[focus_start:]


def test_workbench_file_tree_only_files_look_clickable() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    render_file = script.index("function renderFileNode")
    class_marker = script.index('node.type === "file" ? "is-file" : "is-directory"', render_file)
    file_guard = script.index('if (node.type === "file")', class_marker)
    keyboard_call = script.index("enableKeyboardActivation(row);", file_guard)
    file_node_block = styles[styles.index(".file-node {"):styles.index("}", styles.index(".file-node {"))]

    assert class_marker < file_guard < keyboard_call
    assert "cursor: pointer" not in file_node_block
    assert ".file-node.is-file { cursor: pointer; }" in styles
    assert ".file-node.is-file:hover" in styles


def test_workbench_file_tree_labels_binary_outputs_distinctly() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_file = script.index("function renderFileNode")
    label_call = script.index("formatFileKindLabel(node)", render_file)
    helper = script.index("function formatFileKindLabel")
    directory_case = script.index('node.type === "directory"', helper)
    image_case = script.index('node.kind === "image"', directory_case)
    text_case = script.index('node.kind === "text"', image_case)
    binary_return = script.index('return "BIN"', text_case)

    assert render_file < label_call < helper
    assert helper < directory_case < image_case < text_case < binary_return


def test_workbench_file_tree_kind_labels_have_stable_width() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    file_node_block = styles[styles.index(".file-node {"):styles.index("}", styles.index(".file-node {"))]
    label_block = styles[styles.index(".file-node span:first-child {"):styles.index("}", styles.index(".file-node span:first-child {"))]

    assert "grid-template-columns: 34px minmax(0, 1fr) auto" in file_node_block
    assert "font-family: ui-monospace" in label_block
    assert "text-align: center" in label_block
    assert "border-radius: 999px" in label_block


def test_workbench_file_tree_long_names_do_not_force_horizontal_scroll() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    side_panel_block = styles[styles.index(".side-panel {"):styles.index("}", styles.index(".side-panel {"))]
    file_name_block = styles[styles.index(".file-node strong {"):styles.index("}", styles.index(".file-node strong {"))]
    file_size_block = styles[styles.index(".file-node span:last-child {"):styles.index("}", styles.index(".file-node span:last-child {"))]

    assert "overflow-x: hidden" in side_panel_block
    assert "min-width: 0" in file_name_block
    assert "overflow: hidden" in file_name_block
    assert "text-overflow: ellipsis" in file_name_block
    assert "white-space: nowrap" in file_name_block
    assert "white-space: nowrap" in file_size_block


def test_workbench_binary_file_preview_mentions_size_and_download() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    open_file = script.index("async function openFile")
    binary_message = script.index("Binary preview is not available", open_file)
    size_call = script.index("formatBytes(node.size || 0)", binary_message)
    download_word = script.index("Download", binary_message)

    assert binary_message < size_call
    assert binary_message < download_word


def test_workbench_messages_cannot_expand_over_result_panel() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    chat_block = styles[styles.index(".chat {"):styles.index("}", styles.index(".chat {"))]
    stream_block = styles[styles.index(".stream {"):styles.index("}", styles.index(".stream {"))]
    message_block = styles[styles.index(".message {"):styles.index("}", styles.index(".message {"))]
    message_text_block = styles[styles.index(".message p {"):styles.index("}", styles.index(".message p {"))]

    assert "grid-template-columns: minmax(0, 1fr)" in chat_block
    assert "min-width: 0" in stream_block
    assert "width: 100%" in stream_block
    assert "max-width: 100%" in stream_block
    assert "min-width: 0" in message_block
    assert "width: 100%" in message_block
    assert "max-width: 100%" in message_block
    assert "overflow-wrap: anywhere" in message_text_block


def test_workbench_long_messages_do_not_monopolize_chat_stream() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    message_text_block = styles[styles.index(".message p {"):styles.index("}", styles.index(".message p {"))]

    assert "max-height: 340px" in message_text_block
    assert "overflow: auto" in message_text_block


def test_workbench_message_stream_is_vertically_scrollable_inside_chat_column() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    app_block = styles[styles.index(".app {"):styles.index("}", styles.index(".app {"))]
    sidebar_block = styles[styles.index(".sidebar {"):styles.index("}", styles.index(".sidebar {"))]
    chat_block = styles[styles.index(".chat {"):styles.index("}", styles.index(".chat {"))]
    stream_block = styles[styles.index(".stream {"):styles.index("}", styles.index(".stream {"))]
    detail_block = styles[styles.index(".detail {"):styles.index("}", styles.index(".detail {"))]
    right_block = styles[styles.index(".right {"):styles.index("}", styles.index(".right {"))]

    assert "grid-template-rows: 50px minmax(0, 1fr)" in app_block
    assert "min-height: 0" in sidebar_block
    assert "min-height: 0" in chat_block
    assert "min-height: 0" in stream_block
    assert "min-height: 0" in detail_block
    assert "min-height: 0" in right_block
    assert "overflow: auto" in stream_block


def test_workbench_has_single_column_layout_for_narrow_viewports() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    media_start = styles.index("@media (max-width: 980px)")
    media_block = styles[media_start:]

    assert ".app { height: auto; min-height: 100vh; overflow: visible; grid-template-columns: minmax(0, 1fr); grid-template-rows: auto auto auto auto; }" in media_block
    assert ".rail { display: none; }" in media_block
    assert ".topbar, .sidebar, .chat, .detail { grid-column: 1; grid-row: auto; }" in media_block
    assert ".right { display: none; }" in media_block
    assert ".chat { height: clamp(520px, 70vh, 760px); min-height: 0; border-right: 0; border-bottom: 1px solid var(--line); }" in media_block
    assert ".detail { min-height: 360px; }" in media_block


def test_workbench_reserves_duplicate_overview_column_for_ultrawide_screens() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    wide_start = styles.index("@media (max-width: 1500px)")
    mobile_start = styles.index("@media (max-width: 980px)", wide_start)
    wide_block = styles[wide_start:mobile_start]

    assert ".app { grid-template-columns: 52px 248px minmax(0, .9fr) minmax(0, 1.1fr); }" in wide_block
    assert ".right { display: none; }" in wide_block
    assert ".topbar { grid-column: 2 / 5; }" in wide_block


def test_workbench_mobile_chat_keeps_long_runs_inside_scrollable_stream() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    stream_start = styles.index(".stream {")
    stream_block = styles[stream_start:styles.index("}", stream_start)]
    media_start = styles.index("@media (max-width: 980px)")
    media_block = styles[media_start:]

    assert "overflow: auto" in stream_block
    assert "height: clamp(520px, 70vh, 760px)" in media_block
    assert ".chat { height: clamp(520px, 70vh, 760px); min-height: 0;" in media_block


def test_workbench_composer_controls_stay_inside_chat_column() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    composer_tools = styles[styles.index(".composer-tools {"):styles.index("}", styles.index(".composer-tools {"))]
    tool_buttons = styles[styles.index(".tool-buttons {"):styles.index("}", styles.index(".tool-buttons {"))]
    send_block = styles[styles.index(".send {"):styles.index("}", styles.index(".send {"))]

    assert "display: grid" in composer_tools
    assert "grid-template-columns: minmax(0, 1fr) auto" in composer_tools
    assert "min-width: 0" in composer_tools
    assert "min-width: 0" in tool_buttons
    assert "grid-template-columns: minmax(0, 1fr) auto" in tool_buttons
    assert "white-space: nowrap" in send_block


def test_workbench_compute_values_do_not_force_horizontal_scroll() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    compute_grid = styles[styles.index(".compute-grid {"):styles.index("}", styles.index(".compute-grid {"))]
    compute_card = styles[styles.index(".compute-card {"):styles.index("}", styles.index(".compute-card {"))]
    compute_value = styles[styles.index(".compute-card strong {"):styles.index("}", styles.index(".compute-card strong {"))]

    assert "grid-template-columns: repeat(2, minmax(0, 1fr))" in compute_grid
    assert "min-width: 0" in compute_grid
    assert "min-width: 0" in compute_card
    assert "min-width: 0" in compute_value
    assert "overflow-wrap: anywhere" in compute_value


def test_workbench_result_tab_exposes_run_details_when_right_rail_is_hidden() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    result_pane = markup.index('id="resultPane"')
    result_compute = markup.index('id="resultComputeGrid"', result_pane)
    right_compute = markup.index('id="computeGrid"', result_compute)
    render_compute = script.index("function renderCompute")
    facts = script.index("const facts = computeFacts(snapshot)", render_compute)
    right_render = script.index("renderComputeCards(els.computeGrid, facts)", facts)
    result_render = script.index("renderComputeCards(els.resultComputeGrid, facts)", right_render)
    helper = script.index("function computeFacts")
    run_dir = script.index('["Run Dir"', helper)
    result_root = script.index('["Result Root"', run_dir)

    assert result_pane < result_compute < right_compute
    assert render_compute < facts < right_render < result_render < helper < run_dir < result_root
    assert ".result-compute" in styles


def test_workbench_run_details_include_log_path() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    helper = script.index("function computeFacts")
    run_dir = script.index('["Run Dir"', helper)
    log_path = script.index('["Log Path", compute.logPath || "Pending"]', run_dir)
    result_root = script.index('["Result Root"', log_path)

    assert helper < run_dir < log_path < result_root


def test_workbench_run_details_include_input_data_and_turn_limit() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    helper = script.index("function computeFacts")
    run_id = script.index('["Run", compute.runId || "Pending"]', helper)
    data_path = script.index('["Data Path", snapshot?.dataPath || "None"]', run_id)
    max_turns = script.index('["Max Turns", snapshot?.maxTurns || "Pending"]', data_path)
    model = script.index('["Model", compute.model || "Pending"]', max_turns)

    assert helper < run_id < data_path < max_turns < model


def test_workbench_run_details_show_turn_progress() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    helper = script.index("function computeFacts")
    max_turns = script.index('["Max Turns", snapshot?.maxTurns || "Pending"]', helper)
    turns = script.index('["Turns", formatTurnProgress(snapshot)]', max_turns)
    duration = script.index('["Duration", formatRunDuration(snapshot)]', turns)
    formatter = script.index("function formatTurnProgress")
    current_turn = script.index("snapshot?.compute?.turn", formatter)
    limit = script.index("snapshot?.compute?.maxTurns || snapshot?.maxTurns", current_turn)

    assert helper < max_turns < turns < duration
    assert formatter < current_turn < limit


def test_workbench_turn_progress_does_not_guess_missing_history() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    formatter = script.index("function formatTurnProgress")
    missing_turn = script.index("currentTurn === undefined", formatter)
    active_guard = script.index("isActiveStatus(snapshot?.status)", missing_turn)
    active_pending = script.index('return maxTurns ? `0 / ${maxTurns}` : "Pending"', active_guard)
    terminal_guard = script.index('["completed", "failed"].includes(snapshot?.status)', active_pending)
    not_recorded = script.index('return "Not recorded"', terminal_guard)
    idle_pending = script.index('return "Pending"', not_recorded)

    assert formatter < missing_turn < active_guard < active_pending < terminal_guard < not_recorded < idle_pending


def test_workbench_run_details_include_duration_from_run_timestamps() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    duration_helper = script.index("function formatRunDuration")
    start_parse = script.index('Date.parse(snapshot?.createdAt || "")', duration_helper)
    end_parse = script.index('Date.parse(snapshot?.updatedAt || "")', start_parse)
    elapsed = script.index("formatDurationMs(end - start)", end_parse)
    facts_helper = script.index("function computeFacts")
    max_turns = script.index('["Max Turns", snapshot?.maxTurns || "Pending"]', facts_helper)
    duration_fact = script.index('["Duration", formatRunDuration(snapshot)]', max_turns)
    model = script.index('["Model", compute.model || "Pending"]', duration_fact)

    assert duration_helper < start_parse < end_parse < elapsed
    assert facts_helper < max_turns < duration_fact < model


def test_workbench_running_duration_uses_current_time_instead_of_stale_update() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    helper = script.index("function formatRunDuration")
    start_parse = script.index('Date.parse(snapshot?.createdAt || "")', helper)
    running_check = script.index('snapshot?.status === "running"', start_parse)
    current_time = script.index("Date.now()", running_check)
    updated_parse = script.index('Date.parse(snapshot?.updatedAt || "")', current_time)
    elapsed = script.index("formatDurationMs(end - start)", updated_parse)

    assert helper < start_parse < running_check < current_time < updated_parse < elapsed


def test_workbench_log_path_can_open_task_log_preview() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    compute_cards = script.index("function renderComputeCards")
    log_action = script.index('label === "Log Path"', compute_cards)
    open_log = script.index("async function openLog", log_action)
    copy_value = script.index("els.copyFilePathButton.dataset.copyValue = logPath", open_log)
    copy_enabled = script.index('els.copyFilePathButton.classList.remove("disabled")', copy_value)
    content_fetch = script.index("fetch(`/api/tasks/${logTaskId}/log/content`)", open_log)
    download_href = script.index("els.downloadLink.href = `/api/tasks/${state.taskId}/log/download`", open_log)

    assert compute_cards < log_action < open_log
    assert open_log < copy_value < copy_enabled
    assert open_log < content_fetch
    assert open_log < download_href


def test_workbench_topbar_long_task_titles_stay_single_line() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    topbar_title_area = styles[styles.index(".topbar-title {"):styles.index("}", styles.index(".topbar-title {"))]
    crumbs_block = styles[styles.index(".crumbs {"):styles.index("}", styles.index(".crumbs {"))]
    subtle_block = styles[styles.index(".subtle, .meta, .row-sub {"):styles.index("}", styles.index(".subtle, .meta, .row-sub {"))]
    actions_block = styles[styles.index(".top-actions {"):styles.index("}", styles.index(".top-actions {"))]

    assert 'class="topbar-title"' in markup
    assert "min-width: 0" in topbar_title_area
    assert "overflow: hidden" in crumbs_block
    assert "text-overflow: ellipsis" in crumbs_block
    assert "white-space: nowrap" in crumbs_block
    assert "overflow: hidden" in subtle_block
    assert "white-space: nowrap" in subtle_block
    assert "flex-shrink: 0" in actions_block


def test_workbench_result_summary_wraps_long_paths() -> None:
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    summary_pre = styles[styles.index(".result-summary pre {"):styles.index("}", styles.index(".result-summary pre {"))]

    assert "white-space: pre-wrap" in summary_pre
    assert "overflow-wrap: anywhere" in summary_pre


def test_workbench_search_empty_states_distinguish_no_data_from_no_matches() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "No matching sessions." in script
    assert "No sessions yet." in script
    assert "No matching files." in script
    assert "No output files yet." in script
    assert "tasks.length" in script
    assert "sourceFiles.length" in script


def test_workbench_result_summary_highlights_failed_runs() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert 'snapshot?.status === "failed"' in script
    assert "error-summary" in script
    assert "formatErrorSummary(snapshot.error || finalAnswer)" in script
    assert ".error-summary" in styles


def test_workbench_failed_result_summary_offers_log_and_rerun_actions() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    result_summary = script.index("function renderResultSummary")
    failed_case = script.index('if (snapshot?.status === "failed")', result_summary)
    actions = script.index("result-actions", failed_case)
    rerun = script.index("data-rerun-current", actions)
    primary_action = script.index("result-primary", actions)
    open_log = script.index("data-open-log", rerun)
    open_trace = script.index("data-open-failed-trace", open_log)
    trace_label = script.index("Inspect failed trace", open_trace)
    copy_details = script.index("Copy error details", trace_label)
    click_handler = script.index('event.target.closest("[data-rerun-current]")')
    helper = script.index("async function continueFailedRun")
    terminal_guard = script.index('state.snapshot?.status !== "failed"', helper)
    message_call = script.index("sendSessionMessage", terminal_guard)
    continuation = script.index("Continue from the latest failed run", terminal_guard)

    assert failed_case < actions < primary_action < rerun < open_log < open_trace < trace_label < copy_details
    assert click_handler < helper < terminal_guard < message_call < continuation
    assert ".result-actions" in styles
    assert ".result-primary" in styles


def test_workbench_does_not_duplicate_session_task_in_message_stream() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    render_messages = script.index("function renderMessages")
    fallback = script.index("snapshot.messages?.length", render_messages)
    existing_messages = script.index("? snapshot.messages", fallback)
    task_fallback = script.index("snapshot.task", existing_messages)

    assert render_messages < fallback < existing_messages < task_fallback


def test_history_loader_normalizes_legacy_anthropic_content_blocks(tmp_path: Path) -> None:
    payload = {
        "id": "legacy-anthropic",
        "task": "Run analysis",
        "status": "failed",
        "messages": [
            {
                "role": "assistant",
                "content": "[{'text': 'Inspecting the skill.', 'type': 'text'}, {'id': 'toolu_1', 'input': {}, 'name': 'list_files', 'type': 'tool_use'}]",
            },
            {
                "role": "assistant",
                "content": "[{'id': 'toolu_2', 'input': {}, 'name': 'list_jobs', 'type': 'tool_use'}]",
            },
        ],
    }

    record = record_from_payload(payload)

    assert [message["content"] for message in record.messages] == ["Inspecting the skill."]


def test_workbench_failed_result_summary_can_open_failed_trace() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    click_handler = script.index('event.target.closest("[data-open-failed-trace]")')
    trace_lookup = script.index("const trace = topbarTrace(state.snapshot)", click_handler)
    selected_trace = script.index("state.selectedTraceId = trace.id", trace_lookup)
    render_detail = script.index("renderTraceDetail()", selected_trace)
    activate_trace = script.index('activateTab("trace", { user: true })', render_detail)

    assert click_handler < trace_lookup < selected_trace < render_detail < activate_trace


def test_workbench_failed_result_summary_does_not_require_final_answer() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    result_summary = script.index("function renderResultSummary")
    final_answer = script.index("const finalAnswer =", result_summary)
    failed_case = script.index('if (snapshot?.status === "failed")', final_answer)
    empty_case = script.index("if (!finalAnswer && !summaryFile)", failed_case)

    assert final_answer < failed_case < empty_case


def test_workbench_failed_run_surfaces_concise_provider_error_message() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    task_meta = script.index("function formatTaskDisplayMeta")
    task_error = script.index("formatErrorSummary(task.error)", task_meta)
    result_summary = script.index("function renderResultSummary")
    result_error = script.index("formatErrorSummary(snapshot.error || finalAnswer)", result_summary)
    helper = script.index("function formatErrorSummary")
    single_quote_message = script.index("message'\\s*:\\s*'([^']+)'", helper)
    double_quote_message = script.index('"message"\\s*:\\s*"([^"]+)"', single_quote_message)

    assert task_meta < task_error < helper < result_summary < result_error
    assert helper < single_quote_message < double_quote_message


def test_workbench_clears_composer_once_when_waiting_for_human_input() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "replyTaskId" in script
    assert "state.replyTaskId !== snapshot.id" in script
    assert "els.taskInput.value = \"\";" in script


def test_workbench_reply_placeholder_mentions_latest_agent_question() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    composer = script.index("function renderComposer")
    needs_input = script.index('if (snapshot?.status === "needs_user_input")', composer)
    prompt_const = script.index("const replyPrompt = latestAssistantPrompt(snapshot)", needs_input)
    placeholder = script.index('els.taskInput.placeholder = replyPrompt ? `Reply to: ${replyPrompt}`', prompt_const)
    helper = script.index("function latestAssistantPrompt")
    latest_assistant = script.index('message.role === "assistant"', helper)
    first_line = script.index("split(/\\r?\\n/)[0]", latest_assistant)
    truncation = script.index("slice(0, 96)", first_line)

    assert needs_input < prompt_const < placeholder
    assert helper < latest_assistant < first_line < truncation


def test_workbench_hides_new_task_options_while_replying_to_human_input() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert 'id="composerOptions"' in markup
    assert '"composerOptions"' in script
    reply_start = script.index('if (snapshot?.status === "needs_user_input")')
    reply_hide = script.index('els.composerOptions.hidden = true', reply_start)
    active_start = script.index("if (isActiveStatus(snapshot?.status))", reply_start)
    active_hide = script.index('els.composerOptions.hidden = true', active_start)
    idle_show = script.index('els.composerOptions.hidden = false', active_start)

    assert reply_start < reply_hide < active_start < active_hide < idle_show
    assert ".tool-buttons[hidden]" in styles
    assert "state.replyTaskId = snapshot.id" in script
    assert "state.replyTaskId = null" in script


def test_workbench_new_task_restores_default_composer_state() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "defaultTaskText" in script
    assert "defaultDataPath" in script
    assert "defaultMaxTurns" in script
    assert "resetComposerDefaults();" in script
    assert "function resetComposerDefaults" in script
    assert 'activateTab("trace")' in script


def test_workbench_selected_session_prepares_an_empty_message_composer() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    select_start = script.index("async function selectTask")
    select_end = script.index("async function refresh", select_start)
    select_body = script[select_start:select_end]
    set_snapshot = select_body.index("state.snapshot = snapshot")
    sync_call = select_body.index("syncComposerFromSnapshot(state.snapshot);", set_snapshot)
    render_call = select_body.index("render();", sync_call)
    scroll_call = select_body.index("scrollToLatest();", render_call)
    helper_start = script.index("function syncComposerFromSnapshot")
    session_guard = script.index("state.replyTaskId !== snapshot.id", helper_start)
    clear_input = script.index('els.taskInput.value = ""', session_guard)
    path_fill = script.index("els.dataPathInput.value = snapshot.dataPath || \"\"", clear_input)
    turns_fill = script.index("els.maxTurnsInput.value = String(snapshot.maxTurns || state.defaultMaxTurns)", path_fill)

    assert set_snapshot < sync_call < render_call < scroll_call
    assert helper_start < session_guard < clear_input < path_fill < turns_fill


def test_workbench_selected_active_run_leaves_composer_available_for_steering() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    composer = script.index("function renderComposer")
    active = script.index("if (isActiveStatus(snapshot?.status))", composer)
    send = script.index('els.sendButton.textContent = "Send"', active)
    enabled = script.index("els.taskInput.disabled = false", send)
    pause = script.index("els.pauseButton.hidden = false", enabled)

    assert composer < active < send < enabled < pause


def test_workbench_keeps_send_enabled_and_exposes_pause_while_session_is_active() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert "function isActiveStatus" in script
    assert "isActiveStatus(snapshot?.status)" in script
    assert "els.sendButton.disabled = true" in script
    assert "els.sendButton.disabled = false" in script
    assert ".send:disabled" in styles
    assert "els.pauseButton.hidden = false" in script


def test_workbench_disables_run_button_while_task_create_is_in_flight() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    start = script.index("async function startTask")
    begin_submit = script.index('const submitToken = beginSubmit("create")', start)
    create_call = script.index('const created = await api("/api/sessions"', start)
    clear_submitting = script.index("state.isSubmitting = false", create_call)

    assert begin_submit < create_call < clear_submitting
    assert "if (state.isSubmitting)" in script
    assert '"Creating"' in script


def test_workbench_ignores_stale_task_create_responses_after_navigation() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    start = script.index("async function startTask")
    create_call = script.index('const created = await api("/api/sessions"', start)
    stale_guard = script.index("if (state.submitToken !== submitToken) return", create_call)
    add_message = script.index("addUserMessage(payload.task)", create_call)
    push_state = script.index("history.pushState", create_call)

    assert create_call < stale_guard < add_message < push_state
    assert "function cancelSubmit" in script
    assert "cancelSubmit();" in script


def test_workbench_disables_reply_button_while_resume_is_in_flight() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    start = script.index("async function resumeTask")
    begin_submit = script.index('const submitToken = beginSubmit(isPaused ? "continue" : "reply")', start)
    resume_call = script.index("await api(`/api/sessions/${taskId}/resume`", start)
    clear_submitting = script.index("state.isSubmitting = false", resume_call)

    assert begin_submit < resume_call < clear_submitting
    assert 'const labels = { reply: "Replying", continue: "Continuing", pause: "Pausing", message: "Sending" }' in script
    assert "state.submitMode = null" in script


def test_workbench_ignores_stale_resume_responses_after_navigation() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    start = script.index("async function resumeTask")
    task_id_capture = script.index("const taskId = state.taskId", start)
    resume_call = script.index("await api(`/api/sessions/${taskId}/resume`", start)
    stale_guard = script.index("if (state.submitToken !== submitToken) return", resume_call)
    clear_input = script.index('els.taskInput.value = ""', resume_call)
    select_task = script.index("await selectTask(taskId)", resume_call)

    assert task_id_capture < resume_call < stale_guard < clear_input < select_task


def test_workbench_routes_existing_session_submits_to_message_endpoint_before_new_run_setup() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    start = script.index("async function startTask")
    active_guard = script.index("if (state.taskId && state.snapshot)", start)
    send_message = script.index("await sendSessionMessage()", active_guard)
    close_events = script.index("closeEvents();", start)

    assert active_guard < send_message < close_events


def test_workbench_pause_uses_a_separate_control_from_session_send() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    listener = script.index('els.pauseButton.addEventListener("click", pauseTask)')
    endpoint = script.index("async function pauseTask", listener)
    post = script.index("`/api/sessions/${taskId}/pause`", endpoint)

    assert listener < endpoint < post


def test_workbench_paused_composer_can_send_chat_or_explicitly_continue() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    composer = script.index("function renderComposer")
    paused = script.index('snapshot?.status === "paused"', composer)
    send_label = script.index('els.sendButton.textContent = "Send"', paused)
    continue_control = script.index("els.continueButton.hidden = false", send_label)
    pausing = script.index('snapshot?.status === "pausing"', continue_control)
    running = script.index("isActiveStatus(snapshot?.status)", pausing)
    pause_control = script.index("els.pauseButton.hidden = false", running)

    assert composer < paused < send_label < continue_control < pausing < running < pause_control


def test_workbench_pause_statuses_drive_stream_lifecycle() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    terminal = script.index("function isTerminalStatus")
    paused = script.index('"paused"', terminal)
    active = script.index("function isActiveStatus", paused)
    pausing = script.index('"pausing"', active)

    assert terminal < paused < active < pausing


def test_workbench_only_adds_user_message_after_task_create_succeeds() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    start = script.index("async function startTask")
    create_call = script.index('const created = await api("/api/sessions"', start)
    add_message = script.index("addUserMessage(payload.task)", start)

    assert create_call < add_message


def test_workbench_marks_selected_file_row() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert "selectedFilePath" in script
    assert 'node.path === state.selectedFilePath ? "selected" : ""' in script
    assert "state.selectedFilePath = node.path" in script
    assert "markSelectedFileRow();" in script
    assert ".file-node.selected" in styles


def test_workbench_ignores_stale_text_file_preview_responses() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    open_file_start = script.index("async function openFile")
    task_capture = script.index("const previewTaskId = state.taskId", open_file_start)
    path_capture = script.index("const previewPath = node.path", task_capture)
    fetch_call = script.index("await fetch", path_capture)
    stale_guard = script.index("if (state.taskId !== previewTaskId || state.selectedFilePath !== previewPath) return", fetch_call)

    assert task_capture < path_capture < fetch_call < stale_guard
    assert "const previewPath = node.path" in script


def test_workbench_shows_loading_state_for_text_file_preview() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "Loading preview..." in script
    assert "const previewPath = node.path" in script
    assert script.index("Loading preview...") < script.index("await fetch", script.index("const previewPath = node.path"))


def test_workbench_text_file_preview_shows_fetch_errors() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    open_file_start = script.index("async function openFile")
    catch_block = script.index("catch (error)", open_file_start)
    stale_guard = script.index("if (state.taskId !== previewTaskId || state.selectedFilePath !== previewPath) return", catch_block)
    error_preview = script.index("Unable to load preview", stale_guard)

    assert "Unable to load preview" in script
    assert "catch (error)" in script
    assert "escapeHtml(error.message" in script
    assert catch_block < stale_guard < error_preview


def test_workbench_does_not_force_chat_scroll_when_user_reads_history() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert "function shouldStickToBottom" in script
    assert "const stickToBottom = shouldStickToBottom(els.messageStream)" in script
    assert "if (stickToBottom)" in script
    assert "els.messageStream.scrollTop = els.messageStream.scrollHeight" in script


def test_workbench_offers_jump_to_latest_when_chat_has_new_updates() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "styles.css").read_text(encoding="utf-8")

    assert 'id="scrollLatestButton"' in markup
    assert "scrollLatestButton" in script
    assert "function setScrollLatestVisible" in script
    assert "scrollToLatest();" in script
    assert ".scroll-latest" in styles
    assert ".scroll-latest.visible" in styles


def test_workbench_hides_jump_to_latest_after_manual_scroll_to_bottom() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert 'els.messageStream.addEventListener("scroll"' in script
    assert "if (shouldStickToBottom(els.messageStream))" in script
    assert "setScrollLatestVisible(false)" in script


def test_workbench_result_tab_contains_final_summary_region() -> None:
    source = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html"
    markup = source.read_text(encoding="utf-8")

    assert 'id="resultSummary"' in markup


def test_workbench_detail_tabs_show_trace_and_result_counts() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    tabs = markup.index('id="detailTabs"')
    trace_tab = markup.index('id="traceTab"', tabs)
    memory_tab = markup.index('id="memoryTab"', trace_tab)
    result_tab = markup.index('id="resultTab"', memory_tab)
    file_tab = markup.index('id="fileTab"', result_tab)
    els_list = script.index('"traceTab", "memoryTab", "resultTab", "fileTab"')
    render_start = script.index("function render")
    render_counts = script.index("renderTabCounts(snapshot)", render_start)
    helper = script.index("function renderTabCounts")
    trace_count = script.index("snapshot?.traces?.length || 0", helper)
    result_count = script.index("flattenFiles(snapshot?.resultFiles || []).length", trace_count)
    trace_text = script.index('els.traceTab.textContent = traceCount ? `Trace ${traceCount}` : "Trace"', result_count)
    result_text = script.index('els.resultTab.textContent = resultCount ? `Result ${resultCount}` : "Result"', trace_text)
    file_text = script.index('els.fileTab.textContent = "File"', result_text)

    assert tabs < trace_tab < memory_tab < result_tab < file_tab
    assert els_list < render_start < render_counts < helper
    assert helper < trace_count < result_count < trace_text < result_text < file_text


def test_workbench_has_read_only_memory_tab_for_task_state_and_verified_episodes() -> None:
    root = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static"
    markup = (root / "index.html").read_text(encoding="utf-8")
    script = (root / "assets" / "app.js").read_text(encoding="utf-8")
    styles = (root / "assets" / "styles.css").read_text(encoding="utf-8")

    assert 'id="memoryTab"' in markup
    assert 'data-tab="memory"' in markup
    assert 'id="memoryPane"' in markup
    assert 'id="memoryState"' in markup
    assert 'id="memoryEpisodes"' in markup
    assert "function renderMemory(snapshot)" in script
    assert "function formatMemoryTimestamp(value)" in script
    assert "formatMemoryTimestamp(episode.timestamp)" in script
    assert "snapshot?.memory" in script
    assert "priorEpisodes" in script
    assert ".memory-state-grid" in styles
    assert ".memory-episode" in styles


def test_workbench_detail_tab_switch_resets_shared_scroll_only_when_tab_changes() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    detail_section = markup.index('<section class="detail">')
    detail_body = markup.index('id="detailBody"', detail_section)
    els_list = script.index('"detailTabs", "detailBody"')
    helper = script.index("function activateTab")
    active_before = script.index('document.querySelector(".tab.active")?.dataset.tab', helper)
    changed = script.index("const tabChanged = activeTab !== name", active_before)
    scroll_guard = script.index("if (tabChanged)", changed)
    scroll_reset = script.index("els.detailBody.scrollTop = 0", scroll_guard)

    assert detail_section < detail_body
    assert els_list < helper < active_before < changed < scroll_guard < scroll_reset


def test_workbench_detail_tabs_keep_accessible_state_in_sync() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert '<div class="tabs" id="detailTabs" role="tablist"' in markup
    assert 'id="traceTab" role="tab" aria-controls="tracePane" aria-selected="true"' in markup
    assert 'id="resultTab" role="tab" aria-controls="resultPane" aria-selected="false"' in markup
    assert 'id="fileTab" role="tab" aria-controls="filePane" aria-selected="false"' in markup
    assert 'id="tracePane" role="tabpanel" aria-labelledby="traceTab"' in markup
    assert 'id="resultPane" role="tabpanel" aria-labelledby="resultTab"' in markup
    assert 'id="filePane" role="tabpanel" aria-labelledby="fileTab"' in markup

    helper = script.index("function activateTab")
    tab_loop = script.index('document.querySelectorAll(".tab").forEach', helper)
    selected_state = script.index('tab.setAttribute("aria-selected", String(isActive))', tab_loop)
    pane_loop = script.index('document.querySelectorAll(".detail-pane").forEach', selected_state)
    hidden_state = script.index("pane.hidden = !isActive", pane_loop)

    assert helper < tab_loop < selected_state < pane_loop < hidden_state


def test_workbench_detail_tabs_support_arrow_home_and_end_keys() -> None:
    markup = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "index.html").read_text(encoding="utf-8")
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    assert 'id="traceTab" role="tab" aria-controls="tracePane" aria-selected="true" tabindex="0"' in markup
    assert 'id="resultTab" role="tab" aria-controls="resultPane" aria-selected="false" tabindex="-1"' in markup
    assert 'id="fileTab" role="tab" aria-controls="filePane" aria-selected="false" tabindex="-1"' in markup
    assert 'els.detailTabs.addEventListener("keydown", handleDetailTabKeydown)' in script

    helper = script.index("function handleDetailTabKeydown")
    key_guard = script.index('["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)', helper)
    home_case = script.index('event.key === "Home"', key_guard)
    end_case = script.index('event.key === "End"', home_case)
    left_case = script.index('event.key === "ArrowLeft"', end_case)
    activation = script.index('activateTab(nextTab.dataset.tab, { user: true })', left_case)
    focus = script.index("nextTab.focus()", activation)
    activate_helper = script.index("function activateTab")
    roving_tabindex = script.index("tab.tabIndex = isActive ? 0 : -1", activate_helper)

    assert helper < key_guard < home_case < end_case < left_case < activation < focus
    assert activate_helper < roving_tabindex


def test_workbench_terminal_run_opens_result_tab_by_default_without_overriding_user_choice() -> None:
    script = (Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "static" / "assets" / "app.js").read_text(encoding="utf-8")

    state_block = script[script.index("const state = {"):script.index("};", script.index("const state = {"))]
    tab_click = script.index("els.detailTabs.addEventListener")
    user_activation = script.index("activateTab(event.target.dataset.tab, { user: true })", tab_click)
    select_task = script.index("async function selectTask")
    reset_override = script.index("state.userSelectedTab = false", select_task)
    default_call = script.index("activateDefaultTab(state.snapshot)", reset_override)
    helper = script.index("function activateDefaultTab")
    user_guard = script.index("if (state.userSelectedTab) return", helper)
    completed_case = script.index('snapshot?.status === "completed"', user_guard)
    failed_case = script.index('snapshot?.status === "failed"', completed_case)
    result_tab = script.index('activateTab("result")', failed_case)
    trace_tab = script.index('activateTab("trace")', result_tab)

    assert "userSelectedTab: false" in state_block
    assert tab_click < user_activation
    assert select_task < reset_override < default_call
    assert helper < user_guard < completed_case < failed_case < result_tab < trace_tab


def test_task_history_survives_app_restart(tmp_path: Path) -> None:
    config = _config(tmp_path)
    first_app = create_app(config=config, stream_runner=_fake_stream)
    first_client = TestClient(first_app)
    task_id = first_client.post("/api/tasks", json={"task": "Persistent demo"}).json()["taskId"]
    _wait_for_status(first_client, task_id, "completed")

    second_app = create_app(config=config, stream_runner=_fake_stream)
    second_client = TestClient(second_app)

    tasks = second_client.get("/api/tasks").json()["tasks"]
    snapshot = second_client.get(f"/api/tasks/{task_id}").json()

    assert tasks[0]["id"] == task_id
    assert tasks[0]["status"] == "completed"
    assert snapshot["title"] == "Persistent demo"
    assert snapshot["traces"][0]["toolName"] == "inspect_workflow_skill"


def test_task_history_persistence_uses_atomic_replace() -> None:
    source = Path(__file__).resolve().parents[1] / "src" / "bioagent" / "webapp" / "state.py"
    module = source.read_text(encoding="utf-8")

    persist_start = module.index("def _persist_record")
    persist_body = module[persist_start:module.index("\n\ndef record_to_payload", persist_start)]

    assert "tmp_path" in persist_body
    assert "tmp_path.write_text" in persist_body
    assert "tmp_path.replace(path)" in persist_body


def test_failed_history_without_final_answer_gets_error_notice(tmp_path: Path) -> None:
    config = _config(tmp_path)
    history_dir = config.runs_dir / "web_tasks"
    history_dir.mkdir(parents=True)
    (history_dir / "old_failed.json").write_text(
        """{
          "id": "old_failed",
          "task": "Old failed task",
          "title": "Old failed task",
          "status": "failed",
          "created_at": "2026-01-01T00:00:00+00:00",
          "updated_at": "2026-01-01T00:00:01+00:00",
          "error": "Loading model",
          "messages": [],
          "traces": [],
          "plan": [],
          "events": [],
          "compute": {}
        }""",
        encoding="utf-8",
    )

    app = create_app(config=config, stream_runner=_fake_stream)
    client = TestClient(app)

    snapshot = client.get("/api/tasks/old_failed").json()

    assert snapshot["finalAnswer"] == "Run failed: Loading model"
    assert snapshot["messages"][-1]["content"] == "Run failed: Loading model"


def test_completed_history_prefers_last_assistant_summary_over_turn_limit_notice(tmp_path: Path) -> None:
    config = _config(tmp_path)
    history_dir = config.runs_dir / "web_tasks"
    history_dir.mkdir(parents=True)
    (history_dir / "old_completed.json").write_text(
        """{
          "id": "old_completed",
          "task": "Old completed task",
          "title": "Old completed task",
          "status": "completed",
          "created_at": "2026-01-01T00:00:00+00:00",
          "updated_at": "2026-01-01T00:00:01+00:00",
          "final_answer": "Reached max_turns=20 before a final answer.",
          "messages": [
            {"role": "assistant", "content": "real final summary"},
            {"role": "assistant", "content": "Reached max_turns=20 before a final answer."}
          ],
          "traces": [],
          "plan": [],
          "events": [],
          "compute": {}
        }""",
        encoding="utf-8",
    )

    app = create_app(config=config, stream_runner=_fake_stream)
    client = TestClient(app)

    snapshot = client.get("/api/tasks/old_completed").json()

    assert snapshot["finalAnswer"] == "real final summary"
    assert [message["content"] for message in snapshot["messages"]] == ["real final summary"]


def test_completed_history_marks_matching_final_answer_message(tmp_path: Path) -> None:
    config = _config(tmp_path)
    history_dir = config.runs_dir / "web_tasks"
    history_dir.mkdir(parents=True)
    (history_dir / "old_completed.json").write_text(
        """{
          "id": "old_completed",
          "task": "Old completed task",
          "title": "Old completed task",
          "status": "completed",
          "created_at": "2026-01-01T00:00:00+00:00",
          "updated_at": "2026-01-01T00:00:01+00:00",
          "final_answer": "real final summary",
          "messages": [
            {"role": "assistant", "content": "working note"},
            {"role": "assistant", "content": "real final summary"}
          ],
          "traces": [],
          "plan": [],
          "events": [],
          "compute": {}
        }""",
        encoding="utf-8",
    )

    app = create_app(config=config, stream_runner=_fake_stream)
    client = TestClient(app)

    snapshot = client.get("/api/tasks/old_completed").json()

    assert snapshot["messages"][-1]["content"] == "real final summary"
    assert snapshot["messages"][-1]["final"] is True
    assert snapshot["messages"][0].get("final") is not True
