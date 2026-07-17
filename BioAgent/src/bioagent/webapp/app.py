from __future__ import annotations

import asyncio
import time
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable, Literal

from fastapi import FastAPI, HTTPException
from langchain_core.messages import HumanMessage
from fastapi.responses import FileResponse, PlainTextResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from langchain_core.tools import StructuredTool

from bioagent.config import (
    AgentConfig,
    apply_model_settings,
    clear_model_settings,
    load_model_settings,
    save_model_settings,
)
from bioagent.llm import build_llm

from .state import ChatRunner, ResumeRunner, StreamRunner, TaskStore, session_run_dirs as web_state_session_run_dirs, sse_payload


TEXT_PREVIEW_LIMIT_BYTES = 200_000
TEXT_PREVIEW_TRUNCATED_MESSAGE = "\n\n[Preview truncated. Use Download to fetch the full file.]"


class CreateTaskRequest(BaseModel):
    task: str = Field(..., min_length=1)
    data_path: str | None = None
    max_turns: int = Field(20, ge=1, le=100)


class ResumeTaskRequest(BaseModel):
    user_answer: str = Field(..., min_length=1)
    max_turns: int = Field(20, ge=1, le=100)


class SessionMessageRequest(BaseModel):
    content: str = Field(..., min_length=1)
    max_turns: int = Field(20, ge=1, le=100)


class CancelJobRequest(BaseModel):
    job_id: str = ""


class ModelSettingsRequest(BaseModel):
    provider: Literal["openai_compatible", "anthropic"] = "openai_compatible"
    base_url: str = Field(..., min_length=1)
    api_key: str | None = None
    model_name: str = Field(..., min_length=1)
    temperature: float = Field(0.2, ge=0, le=2)
    request_timeout: float = Field(600, ge=1, le=7200)
    mimo_thinking_type: Literal["", "enabled", "disabled"] = ""
    chat_template_enable_thinking: bool | None = None


ModelTestRunner = Callable[[AgentConfig], dict[str, Any]]


def create_app(
    *,
    config: AgentConfig | None = None,
    stream_runner: StreamRunner | None = None,
    resume_runner: ResumeRunner | None = None,
    chat_runner: ChatRunner | None = None,
    model_test_runner: ModelTestRunner | None = None,
    start_job_watcher: bool | None = None,
) -> FastAPI:
    base_config = config or AgentConfig.from_env(include_model_settings=False)
    effective_config = apply_model_settings(base_config)
    store_kwargs: dict[str, Any] = {
        "config": effective_config,
        "start_job_watcher": (config is None and stream_runner is None) if start_job_watcher is None else start_job_watcher,
    }
    if stream_runner:
        store_kwargs["stream_runner"] = stream_runner
    if resume_runner:
        store_kwargs["resume_runner"] = resume_runner
    if chat_runner:
        store_kwargs["chat_runner"] = chat_runner
    store = TaskStore(**store_kwargs)
    model_test_runner = model_test_runner or _test_model_connection
    app = FastAPI(title="BioAgent Workbench")
    app.state.task_store = store
    app.state.base_config = base_config
    static_dir = Path(__file__).with_name("static")

    @app.get("/api/health")
    def health() -> dict[str, Any]:
        return {
            "status": "ok",
            "frontendServed": (static_dir / "index.html").exists(),
            "provider": store.config.provider,
            "model": store.config.model_name,
            "baseUrl": store.config.base_url,
        }

    @app.get("/api/settings/model")
    def get_model_settings() -> dict[str, Any]:
        return _public_model_settings(store.config, source=_model_settings_source(store.config))

    @app.put("/api/settings/model")
    def update_model_settings(payload: ModelSettingsRequest) -> dict[str, Any]:
        updated = _config_from_settings(store.config, payload)
        save_model_settings(updated.project_root, _persisted_model_settings(updated))
        store.config = updated
        return _public_model_settings(updated, source="saved")

    @app.post("/api/settings/model/test")
    def test_model_settings(payload: ModelSettingsRequest) -> dict[str, Any]:
        candidate = _config_from_settings(store.config, payload)
        try:
            result = dict(model_test_runner(candidate))
        except Exception as exc:
            raise HTTPException(status_code=400, detail=f"Model connection failed: {exc}") from exc
        return {
            "ok": bool(result.get("ok", True)),
            "model": candidate.model_name,
            "baseUrl": candidate.base_url,
            **result,
        }

    @app.delete("/api/settings/model")
    def reset_model_settings() -> dict[str, Any]:
        clear_model_settings(base_config.project_root)
        store.config = base_config
        return _public_model_settings(store.config, source="environment")

    @app.get("/favicon.ico", include_in_schema=False)
    def favicon() -> Response:
        return Response(status_code=204)

    @app.post("/api/tasks")
    @app.post("/api/sessions")
    def create_task(payload: CreateTaskRequest) -> dict[str, str]:
        record = store.create_task(task=payload.task, data_path=payload.data_path, max_turns=payload.max_turns)
        return {"taskId": record.id, "sessionId": record.id, "redirectUrl": f"/tasks/{record.id}"}

    @app.post("/api/tasks/{task_id}/resume")
    @app.post("/api/sessions/{task_id}/resume")
    def resume_task(task_id: str, payload: ResumeTaskRequest) -> dict[str, str]:
        try:
            record = store.resume_task(task_id=task_id, user_answer=payload.user_answer, max_turns=payload.max_turns)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        return {"taskId": record.id, "redirectUrl": f"/tasks/{record.id}"}

    @app.post("/api/tasks/{task_id}/pause")
    @app.post("/api/sessions/{task_id}/pause")
    def pause_task(task_id: str) -> dict[str, str]:
        try:
            record = store.pause_task(task_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        return {"taskId": record.id, "status": record.status}

    @app.post("/api/sessions/{task_id}/interrupt")
    def interrupt_current_tool(task_id: str) -> dict[str, str]:
        try:
            record = store.interrupt_current_tool(task_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        return {"sessionId": record.id, "status": "interrupting"}

    @app.post("/api/sessions/{task_id}/jobs/cancel")
    def cancel_session_job(task_id: str, payload: CancelJobRequest) -> dict[str, Any]:
        try:
            result = store.cancel_job(task_id, payload.job_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        if not result.get("ok"):
            raise HTTPException(status_code=409, detail=result.get("error") or "Could not cancel job")
        return result

    @app.post("/api/sessions/{task_id}/messages", status_code=202)
    def send_session_message(task_id: str, payload: SessionMessageRequest) -> dict[str, str]:
        try:
            record = store.send_message(task_id=task_id, content=payload.content, max_turns=payload.max_turns)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        return {
            "sessionId": record.id,
            "runStatus": record.status,
            "interactionStatus": record.interaction_status,
        }

    @app.get("/api/tasks")
    @app.get("/api/sessions")
    def list_tasks() -> dict[str, Any]:
        tasks = store.list_tasks()
        return {"tasks": tasks, "sessions": tasks}

    @app.get("/api/tasks/{task_id}")
    @app.get("/api/sessions/{task_id}")
    def get_task(task_id: str) -> dict[str, Any]:
        try:
            return store.get_task(task_id).snapshot()
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    @app.get("/api/tasks/{task_id}/events")
    @app.get("/api/sessions/{task_id}/events")
    async def task_events(task_id: str) -> StreamingResponse:
        try:
            record = store.get_task(task_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

        async def stream():
            for event in list(record.events):
                yield sse_payload("bioagent_event", event)
            while record.status not in {"completed", "failed", "needs_user_input", "paused"} or record.interaction_status:
                event = await asyncio.to_thread(record.event_queue.get)
                if event is None:
                    break
                yield sse_payload("bioagent_event", event)
            yield sse_payload("task_snapshot", record.snapshot())

        return StreamingResponse(stream(), media_type="text/event-stream")

    @app.get("/api/tasks/{task_id}/files")
    def get_task_files(task_id: str) -> dict[str, Any]:
        try:
            return {"files": store.get_task(task_id).snapshot()["resultFiles"]}
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    @app.get("/api/tasks/{task_id}/files/content")
    def get_task_file_content(task_id: str, path: str) -> PlainTextResponse:
        record = _record_or_404(store, task_id)
        target = _safe_task_path(record, path)
        if target.suffix.lower() not in {".json", ".txt", ".md", ".csv", ".tsv", ".log", ".py", ".r"}:
            raise HTTPException(status_code=400, detail="File is not a text preview type")
        return PlainTextResponse(_read_text_preview(target))

    @app.get("/api/tasks/{task_id}/files/download")
    def download_task_file(task_id: str, path: str) -> FileResponse:
        record = _record_or_404(store, task_id)
        target = _safe_task_path(record, path)
        return FileResponse(target, filename=target.name)

    @app.get("/api/tasks/{task_id}/log/content")
    def get_task_log_content(task_id: str) -> PlainTextResponse:
        record = _record_or_404(store, task_id)
        return PlainTextResponse(_read_text_preview(_log_path_or_404(record)))

    @app.get("/api/tasks/{task_id}/log/download")
    def download_task_log(task_id: str) -> FileResponse:
        target = _log_path_or_404(_record_or_404(store, task_id))
        return FileResponse(target, filename=target.name)

    if static_dir.exists():
        app.mount("/assets", StaticFiles(directory=static_dir / "assets"), name="assets")
        app.mount("/static/assets", StaticFiles(directory=static_dir / "assets"), name="static-assets")

        @app.get("/")
        def index() -> FileResponse:
            return FileResponse(static_dir / "index.html")

        @app.head("/")
        def index_head() -> FileResponse:
            return FileResponse(static_dir / "index.html")

        @app.get("/tasks/{task_id}")
        @app.get("/sessions/{task_id}")
        def task_page(task_id: str) -> FileResponse:
            return FileResponse(static_dir / "index.html")

        @app.head("/tasks/{task_id}")
        @app.head("/sessions/{task_id}")
        def task_page_head(task_id: str) -> FileResponse:
            return FileResponse(static_dir / "index.html")

        @app.get("/settings")
        def settings_page() -> FileResponse:
            return FileResponse(static_dir / "index.html")

        @app.head("/settings")
        def settings_page_head() -> FileResponse:
            return FileResponse(static_dir / "index.html")

    return app


def _record_or_404(store: TaskStore, task_id: str):
    try:
        return store.get_task(task_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


def _safe_task_path(record, requested: str) -> Path:
    run_dirs = [Path(path).resolve() for path in record.result_roots()]
    session_roots = [Path(path).resolve() for path in web_state_session_run_dirs(record)]
    if not session_roots:
        raise HTTPException(status_code=404, detail="Task has no run directory yet")
    requested_path = Path(requested)
    roots: list[Path] = []
    for resolved in [*run_dirs, *session_roots]:
        if resolved not in roots:
            roots.append(resolved)
    candidates = [requested_path.resolve()] if requested_path.is_absolute() else [(base / requested_path).resolve() for base in roots]
    target = next(
        (
            candidate
            for candidate in candidates
            if candidate.exists()
            and any(session_root in candidate.parents or candidate == session_root for session_root in session_roots)
        ),
        None,
    )
    if target is None:
        raise HTTPException(status_code=404, detail="File not found for this task")
    if target.is_dir():
        raise HTTPException(status_code=400, detail="Path is a directory")
    return target


def _log_path_or_404(record) -> Path:
    if not record.log_path:
        raise HTTPException(status_code=404, detail="Task has no log file yet")
    target = Path(record.log_path).resolve()
    if not target.exists():
        raise HTTPException(status_code=404, detail="Log file not found for this task")
    if target.is_dir():
        raise HTTPException(status_code=400, detail="Log path is a directory")
    return target


def _read_text_preview(path: Path) -> str:
    with path.open("rb") as handle:
        data = handle.read(TEXT_PREVIEW_LIMIT_BYTES + 1)
    if len(data) <= TEXT_PREVIEW_LIMIT_BYTES:
        return data.decode("utf-8", errors="replace")
    preview = data[:TEXT_PREVIEW_LIMIT_BYTES].decode("utf-8", errors="replace")
    return f"{preview}{TEXT_PREVIEW_TRUNCATED_MESSAGE}"


def _config_from_settings(current: AgentConfig, payload: ModelSettingsRequest) -> AgentConfig:
    base_url = payload.base_url.strip().rstrip("/")
    if not base_url.startswith(("http://", "https://")):
        raise HTTPException(status_code=422, detail="base_url must start with http:// or https://")
    api_key = (payload.api_key or "").strip() or current.api_key
    return replace(
        current,
        provider=payload.provider,
        base_url=base_url,
        api_key=api_key,
        model_name=payload.model_name.strip(),
        temperature=payload.temperature,
        request_timeout=payload.request_timeout,
        mimo_thinking_type=payload.mimo_thinking_type,
        chat_template_enable_thinking=payload.chat_template_enable_thinking,
    )


def _persisted_model_settings(config: AgentConfig) -> dict[str, Any]:
    return {
        "provider": config.provider,
        "base_url": config.base_url,
        "api_key": config.api_key,
        "model_name": config.model_name,
        "temperature": config.temperature,
        "request_timeout": config.request_timeout,
        "mimo_thinking_type": config.mimo_thinking_type,
        "chat_template_enable_thinking": config.chat_template_enable_thinking,
    }


def _model_settings_source(config: AgentConfig) -> str:
    return "saved" if load_model_settings(config.project_root) else "environment"


def _public_model_settings(config: AgentConfig, *, source: str) -> dict[str, Any]:
    return {
        "provider": config.provider,
        "baseUrl": config.base_url,
        "modelName": config.model_name,
        "apiKeyConfigured": bool(config.api_key.strip()),
        "apiKeyMasked": config.mask_api_key(),
        "temperature": config.temperature,
        "requestTimeout": config.request_timeout,
        "mimoThinkingType": config.mimo_thinking_type,
        "chatTemplateEnableThinking": config.chat_template_enable_thinking,
        "source": source,
    }


def _test_model_connection(config: AgentConfig) -> dict[str, Any]:
    started = time.perf_counter()
    probe = StructuredTool.from_function(
        lambda value: value,
        name="bioagent_connection_probe",
        description="Return the supplied value. Use this tool for the connection test.",
    )
    response = build_llm(config).bind_tools([probe], tool_choice=probe.name).invoke(
        [HumanMessage(content="Call bioagent_connection_probe exactly once with value='OK'.")]
    )
    invalid_tool_calls = getattr(response, "invalid_tool_calls", None) or []
    if invalid_tool_calls:
        detail = str(invalid_tool_calls[0].get("error") or invalid_tool_calls[0].get("args") or "unknown parsing error")
        raise RuntimeError(f"Model returned an invalid tool call: {detail}")
    tool_calls = getattr(response, "tool_calls", None) or []
    if not tool_calls:
        raise RuntimeError("Model responded, but structured tool calling was not verified.")
    elapsed_ms = round((time.perf_counter() - started) * 1000)
    return {
        "ok": True,
        "latencyMs": elapsed_ms,
        "preview": f"Tool calling verified: {tool_calls[0].get('name')}",
    }


app = create_app()
