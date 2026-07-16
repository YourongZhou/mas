from __future__ import annotations

import asyncio
from pathlib import Path
from typing import Any

from fastapi import FastAPI, HTTPException
from fastapi.responses import FileResponse, PlainTextResponse, Response, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from bioagent.config import AgentConfig

from .state import ResumeRunner, StreamRunner, TaskStore, sse_payload


TEXT_PREVIEW_LIMIT_BYTES = 200_000
TEXT_PREVIEW_TRUNCATED_MESSAGE = "\n\n[Preview truncated. Use Download to fetch the full file.]"


class CreateTaskRequest(BaseModel):
    task: str = Field(..., min_length=1)
    data_path: str | None = None
    max_turns: int = Field(20, ge=1, le=100)


class ResumeTaskRequest(BaseModel):
    user_answer: str = Field(..., min_length=1)
    max_turns: int = Field(20, ge=1, le=100)


def create_app(
    *,
    config: AgentConfig | None = None,
    stream_runner: StreamRunner | None = None,
    resume_runner: ResumeRunner | None = None,
) -> FastAPI:
    config = config or AgentConfig.from_env()
    store_kwargs: dict[str, Any] = {"config": config}
    if stream_runner:
        store_kwargs["stream_runner"] = stream_runner
    if resume_runner:
        store_kwargs["resume_runner"] = resume_runner
    store = TaskStore(**store_kwargs)
    app = FastAPI(title="BioAgent Workbench")
    app.state.task_store = store
    static_dir = Path(__file__).with_name("static")

    @app.get("/api/health")
    def health() -> dict[str, Any]:
        return {
            "status": "ok",
            "frontendServed": (static_dir / "index.html").exists(),
            "model": config.model_name,
            "baseUrl": config.base_url,
        }

    @app.get("/favicon.ico", include_in_schema=False)
    def favicon() -> Response:
        return Response(status_code=204)

    @app.post("/api/tasks")
    def create_task(payload: CreateTaskRequest) -> dict[str, str]:
        record = store.create_task(task=payload.task, data_path=payload.data_path, max_turns=payload.max_turns)
        return {"taskId": record.id, "redirectUrl": f"/tasks/{record.id}"}

    @app.post("/api/tasks/{task_id}/resume")
    def resume_task(task_id: str, payload: ResumeTaskRequest) -> dict[str, str]:
        try:
            record = store.resume_task(task_id=task_id, user_answer=payload.user_answer, max_turns=payload.max_turns)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc
        return {"taskId": record.id, "redirectUrl": f"/tasks/{record.id}"}

    @app.get("/api/tasks")
    def list_tasks() -> dict[str, Any]:
        return {"tasks": store.list_tasks()}

    @app.get("/api/tasks/{task_id}")
    def get_task(task_id: str) -> dict[str, Any]:
        try:
            return store.get_task(task_id).snapshot()
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

    @app.get("/api/tasks/{task_id}/events")
    async def task_events(task_id: str) -> StreamingResponse:
        try:
            record = store.get_task(task_id)
        except KeyError as exc:
            raise HTTPException(status_code=404, detail=str(exc)) from exc

        async def stream():
            for event in list(record.events):
                yield sse_payload("bioagent_event", event)
            while record.status not in {"completed", "failed", "needs_user_input"}:
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
        def task_page(task_id: str) -> FileResponse:
            return FileResponse(static_dir / "index.html")

        @app.head("/tasks/{task_id}")
        def task_page_head(task_id: str) -> FileResponse:
            return FileResponse(static_dir / "index.html")

    return app


def _record_or_404(store: TaskStore, task_id: str):
    try:
        return store.get_task(task_id)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


def _safe_task_path(record, requested: str) -> Path:
    if not record.run_dir:
        raise HTTPException(status_code=404, detail="Task has no run directory yet")
    root = Path(record.run_dir).resolve()
    requested_path = Path(requested)
    candidate_roots = [record.result_root(), root]
    roots: list[Path] = []
    for candidate_root in candidate_roots:
        if candidate_root is None:
            continue
        resolved = candidate_root.resolve()
        if resolved not in roots:
            roots.append(resolved)
    candidates = [requested_path.resolve()] if requested_path.is_absolute() else [(base / requested_path).resolve() for base in roots]
    target = next((candidate for candidate in candidates if candidate.exists() and (root in candidate.parents or candidate == root)), None)
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


app = create_app()
