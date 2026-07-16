from __future__ import annotations

import ast
import json
import queue
import threading
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterator

from bioagent.config import AgentConfig
from bioagent.runner import resume_bio_agent, run_bio_agent_stream


StreamRunner = Callable[..., Iterator[dict[str, Any]]]
ResumeRunner = Callable[..., dict[str, Any]]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


@dataclass
class TaskRecord:
    id: str
    task: str
    title: str
    data_path: str | None
    max_turns: int
    status: str = "queued"
    created_at: str = field(default_factory=utc_now)
    updated_at: str = field(default_factory=utc_now)
    messages: list[dict[str, Any]] = field(default_factory=list)
    traces: list[dict[str, Any]] = field(default_factory=list)
    plan: list[dict[str, Any]] = field(default_factory=list)
    events: list[dict[str, Any]] = field(default_factory=list)
    event_queue: queue.Queue[dict[str, Any] | None] = field(default_factory=queue.Queue)
    compute: dict[str, Any] = field(default_factory=dict)
    final_answer: str = ""
    run_dir: str = ""
    log_path: str = ""
    error: str = ""

    def snapshot(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "title": self.title,
            "task": self.task,
            "dataPath": self.data_path,
            "maxTurns": self.max_turns,
            "status": self.status,
            "createdAt": self.created_at,
            "updatedAt": self.updated_at,
            "messages": self.messages,
            "traces": self.traces,
            "plan": self.plan,
            "resultFiles": build_file_tree(self.result_root()),
            "compute": self.compute_snapshot(),
            "finalAnswer": self.final_answer,
            "error": self.error,
        }

    def result_root(self) -> Path | None:
        if not self.run_dir:
            return None
        root = Path(self.run_dir)
        web_result = root / f"web_{self.id}"
        if web_result.is_dir():
            return web_result
        if (root / "outputs").is_dir():
            return root / "outputs"
        return root

    def compute_snapshot(self) -> dict[str, Any]:
        result_root = ""
        root = self.result_root()
        if root:
            result_root = str(root)
        return {
            "runId": self.compute.get("runId", ""),
            "runDir": self.run_dir,
            "logPath": self.log_path,
            "resultRoot": result_root,
            "model": self.compute.get("model", ""),
            "baseUrl": self.compute.get("baseUrl", ""),
            "envProfile": self.compute.get("envProfile", ""),
            "currentTool": self.compute.get("currentTool", ""),
            "turn": self.compute.get("turn"),
            "maxTurns": self.compute.get("maxTurns", self.max_turns),
        }


class TaskStore:
    def __init__(
        self,
        *,
        config: AgentConfig,
        stream_runner: StreamRunner = run_bio_agent_stream,
        resume_runner: ResumeRunner = resume_bio_agent,
    ) -> None:
        self.config = config
        self.stream_runner = stream_runner
        self.resume_runner = resume_runner
        self.history_dir = config.runs_dir / "web_tasks"
        self.history_dir.mkdir(parents=True, exist_ok=True)
        self._tasks: dict[str, TaskRecord] = {}
        self._lock = threading.Lock()
        self._load_history()

    def create_task(self, *, task: str, data_path: str | None = None, max_turns: int = 20) -> TaskRecord:
        task_id = uuid.uuid4().hex[:8]
        title = task.strip().splitlines()[0][:90] or "Untitled BioAgent task"
        record = TaskRecord(id=task_id, task=task, title=title, data_path=data_path, max_turns=max_turns)
        record.compute.update({"model": self.config.model_name, "baseUrl": self.config.base_url})
        with self._lock:
            self._tasks[task_id] = record
            self._persist_record(record)
        thread = threading.Thread(target=self._run_task, args=(record,), name=f"bioagent-web-{task_id}", daemon=True)
        thread.start()
        return record

    def resume_task(self, *, task_id: str, user_answer: str, max_turns: int = 20) -> TaskRecord:
        record = self.get_task(task_id)
        if record.status != "needs_user_input":
            raise ValueError(f"Task {task_id} is not waiting for user input")
        with self._lock:
            record.status = "running"
            record.error = ""
            record.event_queue = queue.Queue()
            record.messages.append({"role": "user", "content": user_answer, "createdAt": utc_now(), "resume": True})
            record.updated_at = utc_now()
            self._persist_record(record)
        thread = threading.Thread(
            target=self._run_resume,
            args=(record, user_answer, max_turns),
            name=f"bioagent-web-resume-{task_id}",
            daemon=True,
        )
        thread.start()
        return record

    def list_tasks(self) -> list[dict[str, Any]]:
        with self._lock:
            records = list(self._tasks.values())
        return [
            {
                "id": record.id,
                "title": record.title,
                "task": record.task,
                "status": record.status,
                "dataPath": record.data_path,
                "createdAt": record.created_at,
                "updatedAt": record.updated_at,
                "error": summarize_error(record.error),
                "logPath": record.log_path,
                "messageCount": len(record.messages),
                "traceCount": len(record.traces),
                "resultFileCount": count_result_files(record.result_root()),
                "resultFileNames": result_file_names(record.result_root()),
            }
            for record in sorted(records, key=lambda item: item.updated_at, reverse=True)
        ]

    def get_task(self, task_id: str) -> TaskRecord:
        with self._lock:
            try:
                return self._tasks[task_id]
            except KeyError as exc:
                raise KeyError(f"Unknown task id: {task_id}") from exc

    def _run_task(self, record: TaskRecord) -> None:
        self._apply_event(record, {"type": "task_started", "status": "running"})
        try:
            result_dir = str(self.config.runs_dir / f"web_{record.id}")
            for event in self.stream_runner(
                record.task,
                data_path=record.data_path,
                result_dir=result_dir,
                max_turns=record.max_turns,
            ):
                self._apply_event(record, event)
        except Exception as exc:
            self._apply_event(record, {"type": "error", "error": str(exc)})
        finally:
            record.event_queue.put(None)

    def _run_resume(self, record: TaskRecord, user_answer: str, max_turns: int) -> None:
        self._apply_event(record, {"type": "resume_started", "status": "running"})
        try:
            resume_id = str(record.compute.get("runId") or "")
            result = self.resume_runner(resume_id, user_answer, max_turns=max_turns)
            self._apply_event(
                record,
                {
                    "type": "final",
                    "content": str(result.get("final_answer") or ""),
                    "status": str(result.get("status") or "completed"),
                    "result": result,
                },
            )
            self._apply_event(
                record,
                {
                    "type": "run_end",
                    "run_id": resume_id or record.compute.get("runId", ""),
                    "run_dir": str(result.get("run_dir") or record.run_dir),
                    "log_path": str(result.get("log_path") or record.log_path),
                    "status": str(result.get("status") or "completed"),
                    "result": result,
                },
            )
        except Exception as exc:
            self._apply_event(record, {"type": "error", "error": str(exc)})
        finally:
            record.event_queue.put(None)

    def _apply_event(self, record: TaskRecord, event: dict[str, Any]) -> None:
        with self._lock:
            apply_event_to_record(record, event)
            record.events.append(event)
            record.updated_at = utc_now()
            self._persist_record(record)
        record.event_queue.put(event)

    def _load_history(self) -> None:
        for path in sorted(self.history_dir.glob("*.json")):
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                record = record_from_payload(payload)
            except Exception:
                continue
            self._tasks[record.id] = record

    def _persist_record(self, record: TaskRecord) -> None:
        path = self.history_dir / f"{record.id}.json"
        tmp_path = path.with_name(f"{path.name}.tmp")
        tmp_path.write_text(json.dumps(record_to_payload(record), ensure_ascii=False, indent=2), encoding="utf-8")
        tmp_path.replace(path)


def record_to_payload(record: TaskRecord) -> dict[str, Any]:
    return {
        "id": record.id,
        "task": record.task,
        "title": record.title,
        "data_path": record.data_path,
        "max_turns": record.max_turns,
        "status": record.status,
        "created_at": record.created_at,
        "updated_at": record.updated_at,
        "messages": record.messages,
        "traces": record.traces,
        "plan": record.plan,
        "events": record.events,
        "compute": record.compute,
        "final_answer": record.final_answer,
        "run_dir": record.run_dir,
        "log_path": record.log_path,
        "error": record.error,
    }


def record_from_payload(payload: dict[str, Any]) -> TaskRecord:
    record = TaskRecord(
        id=str(payload["id"]),
        task=str(payload.get("task") or ""),
        title=str(payload.get("title") or "Untitled BioAgent task"),
        data_path=payload.get("data_path"),
        max_turns=int(payload.get("max_turns") or 20),
        status=str(payload.get("status") or "queued"),
        created_at=str(payload.get("created_at") or utc_now()),
        updated_at=str(payload.get("updated_at") or utc_now()),
        messages=list(payload.get("messages") or []),
        traces=list(payload.get("traces") or []),
        plan=list(payload.get("plan") or []),
        events=list(payload.get("events") or []),
        compute=dict(payload.get("compute") or {}),
        final_answer=str(payload.get("final_answer") or ""),
        run_dir=str(payload.get("run_dir") or ""),
        log_path=str(payload.get("log_path") or ""),
        error=str(payload.get("error") or ""),
    )
    for trace in record.traces:
        output = trace.get("output")
        if isinstance(output, dict) and output.get("ok") is False:
            trace["status"] = "failed"
    if record.status == "running":
        record.status = "failed"
        record.error = "Run interrupted before completion."
        record.compute["currentTool"] = ""
        for trace in record.traces:
            if trace.get("status") == "running":
                trace["status"] = "failed"
    if record.status == "failed" and record.error and not record.final_answer:
        record.final_answer = f"Run failed: {record.error}"
        record.messages.append({"role": "assistant", "content": record.final_answer, "createdAt": record.updated_at, "error": True})
    if record.status == "completed" and record.final_answer.startswith("Reached max_turns="):
        assistant_messages = [
            str(message.get("content") or "").strip()
            for message in record.messages
            if message.get("role") == "assistant" and str(message.get("content") or "").strip()
        ]
        useful_messages = [message for message in assistant_messages if not message.startswith("Reached max_turns=")]
        if useful_messages:
            record.final_answer = useful_messages[-1]
            record.messages = [
                message
                for message in record.messages
                if not str(message.get("content") or "").strip().startswith("Reached max_turns=")
            ]
    if record.status == "completed" and record.final_answer:
        for message in reversed(record.messages):
            if message.get("role") == "assistant" and str(message.get("content") or "").strip() == record.final_answer.strip():
                message["final"] = True
                break
    normalize_plan_from_traces(record)
    return record


def apply_event_to_record(record: TaskRecord, event: dict[str, Any]) -> None:
    event_type = str(event.get("type", ""))
    if event_type == "task_started":
        record.status = "running"
        ensure_plan_step(record, "Understand request", "completed")
        return
    if event_type == "resume_started":
        record.status = "running"
        ensure_plan_step(record, "Wait for user input", "completed")
        ensure_plan_step(record, "Resume with user input", "completed")
        return
    if event_type == "run_start":
        record.status = "running"
        record.run_dir = str(event.get("run_dir") or record.run_dir)
        record.log_path = str(event.get("log_path") or record.log_path)
        record.compute["runId"] = str(event.get("run_id") or record.compute.get("runId") or "")
        ensure_plan_step(record, "Start BioAgent run", "completed")
        return
    if event_type == "turn_start":
        record.compute["turn"] = int(event.get("turn") or 0)
        record.compute["maxTurns"] = int(event.get("max_turns") or record.max_turns)
        return
    if event_type == "assistant_message":
        content = str(event.get("content") or "")
        if content:
            record.messages.append({"role": "assistant", "content": content, "createdAt": utc_now()})
        return
    if event_type == "tool_call":
        call_id = str(event.get("call_id") or uuid.uuid4().hex[:8])
        tool_name = str(event.get("tool_name") or "tool")
        record.compute["currentTool"] = tool_name
        record.traces.append(
            {
                "id": call_id,
                "toolName": tool_name,
                "status": "running",
                "input": event.get("args") or {},
                "output": None,
                "createdAt": utc_now(),
                "updatedAt": utc_now(),
            }
        )
        ensure_plan_step(record, f"Call {tool_name}", "running")
        return
    if event_type == "tool_result":
        call_id = str(event.get("call_id") or "")
        trace = next((item for item in record.traces if item["id"] == call_id), None)
        tool_name = str(event.get("tool_name") or (trace or {}).get("toolName") or "tool")
        status = "completed" if event.get("ok") else "failed"
        if trace is not None:
            trace["status"] = status
            trace["output"] = event.get("result")
            trace["updatedAt"] = utc_now()
        record.compute["currentTool"] = ""
        ensure_plan_step(record, f"Call {tool_name}", status)
        ensure_plan_step(record, f"Finish {tool_name}", status)
        return
    if event_type == "final":
        content = str(event.get("content") or "")
        record.final_answer = content
        record.status = str(event.get("status") or record.status or "completed")
        if content:
            duplicate_prompt = (
                record.status == "needs_user_input"
                and bool(record.messages)
                and record.messages[-1].get("role") == "assistant"
                and str(record.messages[-1].get("content") or "") == content
            )
            if not duplicate_prompt:
                message = {"role": "assistant", "content": content, "createdAt": utc_now()}
                if record.status == "completed":
                    message["final"] = True
                record.messages.append(message)
        if record.status == "needs_user_input":
            ensure_plan_step(record, "Wait for user input", "running")
            return
        ensure_plan_step(record, "Summarize result", "completed")
        return
    if event_type == "run_end":
        record.status = str(event.get("status") or "completed")
        record.run_dir = str(event.get("run_dir") or record.run_dir)
        record.log_path = str(event.get("log_path") or record.log_path)
        record.compute["runId"] = str(event.get("run_id") or record.compute.get("runId") or "")
        result = event.get("result") if isinstance(event.get("result"), dict) else {}
        if not record.final_answer:
            record.final_answer = str(result.get("final_answer") or "")
        if record.status == "needs_user_input":
            return
        ensure_plan_step(record, "Finish run", "completed" if record.status == "completed" else record.status)
        return
    if event_type == "error":
        record.status = "failed"
        record.error = str(event.get("error") or "Unknown error")
        if not record.final_answer:
            record.final_answer = f"Run failed: {record.error}"
        record.messages.append({"role": "assistant", "content": record.final_answer, "createdAt": utc_now(), "error": True})
        ensure_plan_step(record, "Handle error", "failed")


def ensure_plan_step(record: TaskRecord, title: str, status: str) -> None:
    existing = next((item for item in record.plan if item["title"] == title), None)
    if existing:
        existing["status"] = status
        return
    record.plan.append({"id": uuid.uuid4().hex[:8], "title": title, "status": status})


def normalize_plan_from_traces(record: TaskRecord) -> None:
    final_by_tool = {
        str(trace.get("toolName") or "tool"): str(trace.get("status") or "")
        for trace in record.traces
        if str(trace.get("status") or "") in {"completed", "failed"}
    }
    for item in record.plan:
        title = str(item.get("title") or "")
        if not title.startswith("Call ") or item.get("status") != "running":
            continue
        tool_name = title.removeprefix("Call ")
        status = final_by_tool.get(tool_name)
        if status:
            item["status"] = status


def build_file_tree(root: Path | None) -> list[dict[str, Any]]:
    if root is None or not root.exists():
        return []
    children = sorted(root.iterdir(), key=file_tree_sort_key)
    return [_file_node(path) for path in children if path.name != "__pycache__"]


def count_result_files(root: Path | None) -> int:
    if root is None or not root.exists():
        return 0
    if root.is_file():
        return 1
    return sum(count_result_files(path) for path in root.iterdir() if path.name != "__pycache__")


def result_file_names(root: Path | None) -> list[str]:
    return [node["name"] for node in flatten_file_nodes(build_file_tree(root))]


def flatten_file_nodes(nodes: list[dict[str, Any]]) -> list[dict[str, Any]]:
    files: list[dict[str, Any]] = []
    for node in nodes:
        if node.get("type") == "file":
            files.append(node)
        files.extend(flatten_file_nodes(node.get("children", [])))
    return files


def file_tree_sort_key(path: Path) -> tuple[int, int, str]:
    if path.is_dir():
        return (0, 0, path.name.lower())
    return (1, result_file_rank(path.name), path.name.lower())


def result_file_rank(name: str) -> int:
    lower = name.lower()
    if lower == "summary.json":
        return 0
    if "summary" in lower and lower.endswith(".json"):
        return 1
    if "umap" in lower:
        return 2
    if "qc" in lower:
        return 3
    if "pca" in lower:
        return 4
    if "cluster" in lower:
        return 5
    if lower.endswith(".h5ad"):
        return 8
    return 20


def _file_node(path: Path) -> dict[str, Any]:
    if path.is_dir():
        return {
            "name": path.name,
            "path": str(path),
            "type": "directory",
            "children": build_file_tree(path),
        }
    return {
        "name": path.name,
        "path": str(path),
        "type": "file",
        "size": path.stat().st_size,
        "kind": file_kind(path),
    }


def file_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in {".png", ".jpg", ".jpeg", ".gif", ".svg", ".webp"}:
        return "image"
    if suffix in {".json", ".txt", ".md", ".csv", ".tsv", ".log", ".py", ".r"}:
        return "text"
    return "binary"


def summarize_error(error: str) -> str:
    if not error:
        return ""
    marker = " - "
    if marker not in error:
        return error
    _, payload_text = error.split(marker, 1)
    try:
        payload = ast.literal_eval(payload_text)
    except (SyntaxError, ValueError):
        return error
    if not isinstance(payload, dict):
        return error
    provider_error = payload.get("error")
    if not isinstance(provider_error, dict):
        return error
    message = provider_error.get("message")
    return str(message) if message else error


def sse_payload(event: str, data: dict[str, Any]) -> str:
    return f"event: {event}\ndata: {json.dumps(data, ensure_ascii=False)}\n\n"
