from __future__ import annotations

import ast
import inspect
import json
import queue
import subprocess
import threading
import time
import uuid
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Iterator

from bioagent.config import AgentConfig
from bioagent.llm import build_llm, message_text
from bioagent.runner import answer_bio_agent_message, resume_bio_agent, run_bio_agent_stream
from bioagent.skills.registry import list_workflow_skills
from bioagent.tools.jobs import DockerJobManager, TERMINAL_JOB_STATUSES, cancel_run_jobs
from langchain_core.messages import HumanMessage

from .protocol import envelope_event


StreamRunner = Callable[..., Iterator[dict[str, Any]]]
ResumeRunner = Callable[..., dict[str, Any]]
ChatRunner = Callable[[dict[str, Any], str], str]
TitleRunner = Callable[[AgentConfig, str], str]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def fallback_session_title(task: str) -> str:
    first_line = next((line.strip() for line in task.splitlines() if line.strip()), "")
    return " ".join(first_line.split())[:90] or "Untitled BioAgent task"


def clean_session_title(value: str) -> str:
    title = " ".join(str(value or "").strip().splitlines()[0:1]).strip()
    title = title.removeprefix("Title:").removeprefix("title:").removeprefix("标题：").removeprefix("标题:")
    title = title.strip(" `\"'")
    return " ".join(title.split())[:48]


def generate_session_title(config: AgentConfig, task: str) -> str:
    """Ask the configured chat model for a compact, sidebar-friendly session name."""
    prompt = (
        "为下面的 BioAgent 用户任务生成一个简短会话标题。\n"
        "只输出标题本身，不要解释、引号或 Markdown。使用任务的语言，尽量 4-12 个字，最多 48 个字符。\n\n"
        f"任务：\n{task.strip()[:2000]}"
    )
    response = build_llm(config).invoke([HumanMessage(content=prompt)])
    return message_text(response)


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
    plan_state: dict[str, Any] = field(default_factory=dict)
    active_question: dict[str, Any] = field(default_factory=dict)
    events: list[dict[str, Any]] = field(default_factory=list)
    event_queue: queue.Queue[dict[str, Any] | None] = field(default_factory=queue.Queue)
    compute: dict[str, Any] = field(default_factory=dict)
    memory: dict[str, Any] = field(default_factory=dict)
    session_status: str = "open"
    interaction_status: str = ""
    runs: list[dict[str, Any]] = field(default_factory=list)
    active_run_id: str = ""
    pending_user_messages: list[dict[str, str]] = field(default_factory=list)
    pending_context_events: list[dict[str, Any]] = field(default_factory=list)
    delivered_job_callbacks: dict[str, str] = field(default_factory=dict)
    final_answer: str = ""
    run_dir: str = ""
    log_path: str = ""
    error: str = ""
    pause_event: threading.Event = field(default_factory=threading.Event, repr=False)

    def snapshot(self) -> dict[str, Any]:
        semantic_plan = [dict(item) for item in self.plan_state.get("steps") or [] if isinstance(item, dict)]
        return {
            "id": self.id,
            "title": self.title,
            "task": self.task,
            "dataPath": self.data_path,
            "maxTurns": self.max_turns,
            "status": self.status,
            "sessionStatus": self.session_status,
            "runStatus": self.status,
            "interactionStatus": self.interaction_status,
            "runs": self.runs,
            "activeRunId": self.active_run_id,
            "createdAt": self.created_at,
            "updatedAt": self.updated_at,
            "messages": self.messages,
            "traces": self.traces,
            "plan": semantic_plan,
            "planState": self.plan_state,
            "activity": build_tool_activity(self.traces),
            "activeQuestion": self.active_question,
            "events": self.events,
            "resultFiles": self.result_files(),
            "compute": self.compute_snapshot(),
            "memory": self.memory,
            "finalAnswer": self.final_answer,
            "error": self.error,
        }

    def result_root(self) -> Path | None:
        roots = self.result_roots()
        return roots[-1] if roots else None

    def result_roots(self) -> list[Path]:
        run_dirs = [Path(path) for path in session_run_dirs(self)]
        if self.run_dir and Path(self.run_dir) not in run_dirs:
            run_dirs.append(Path(self.run_dir))
        canonical: list[Path] = []
        for root in run_dirs:
            web_result = root / f"web_{self.id}"
            candidate = web_result if web_result.is_dir() else root / "outputs"
            if candidate.is_dir() and candidate not in canonical:
                canonical.append(candidate)
        if canonical:
            return canonical
        if self.run_dir and Path(self.run_dir).exists():
            return [Path(self.run_dir)]
        return []

    def result_files(self) -> list[dict[str, Any]]:
        return build_merged_file_tree(self.result_roots())

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
        chat_runner: ChatRunner = answer_bio_agent_message,
        title_runner: TitleRunner | None = None,
        start_job_watcher: bool = False,
        job_watch_interval_s: float = 2.0,
    ) -> None:
        self.config = config
        self.stream_runner = stream_runner
        self.resume_runner = resume_runner
        self.chat_runner = chat_runner
        self.title_runner = title_runner
        self.history_dir = config.runs_dir / "web_tasks"
        self.history_dir.mkdir(parents=True, exist_ok=True)
        self._tasks: dict[str, TaskRecord] = {}
        self._lock = threading.Lock()
        self._load_history()
        self._job_watch_interval_s = max(0.2, job_watch_interval_s)
        self._job_watcher_stop = threading.Event()
        if start_job_watcher:
            threading.Thread(target=self._job_watch_loop, name="bioagent-job-watcher", daemon=True).start()

    def create_task(
        self,
        *,
        task: str,
        data_path: str | None = None,
        max_turns: int = 20,
        attachments: list[dict[str, Any]] | None = None,
    ) -> TaskRecord:
        task_id = uuid.uuid4().hex[:8]
        title = fallback_session_title(task)
        record = TaskRecord(id=task_id, task=task, title=title, data_path=data_path, max_turns=max_turns)
        record.compute.update({"model": self.config.model_name, "baseUrl": self.config.base_url})
        message = append_session_message(record, "user", task, delivery="consumed", attachments=attachments or [])
        start_session_run(record, kind="agent", trigger_message_id=message["id"])
        with self._lock:
            self._tasks[task_id] = record
            self._persist_record(record)
        thread = threading.Thread(target=self._run_task, args=(record,), name=f"bioagent-web-{task_id}", daemon=True)
        thread.start()
        if self.title_runner is not None:
            threading.Thread(
                target=self._generate_session_title,
                args=(record,),
                name=f"bioagent-title-{task_id}",
                daemon=True,
            ).start()
        return record

    def _generate_session_title(self, record: TaskRecord) -> None:
        """Replace the immediate local fallback without delaying the agent run."""
        try:
            title = clean_session_title(self.title_runner(self.config, record.task)) if self.title_runner else ""
        except Exception:
            return
        if not title or title == record.title:
            return
        with self._lock:
            record.title = title
            event = envelope_event(
                {"type": "session_title", "title": title},
                session_id=record.id,
                current_run_id=str(record.compute.get("runId") or ""),
                prior_events=record.events,
            )
            record.events.append(event)
            record.updated_at = utc_now()
            self._persist_record(record)
        record.event_queue.put(event)

    def send_message(
        self,
        *,
        task_id: str,
        content: str,
        max_turns: int = 20,
        attachments: list[dict[str, Any]] | None = None,
    ) -> TaskRecord:
        record = self.get_task(task_id)
        content = content.strip()
        if not content:
            raise ValueError("Message content is required")
        if record.interaction_status:
            raise ValueError(f"Session {task_id} is already answering a message")
        if record.status == "needs_user_input":
            return self.resume_task(
                task_id=task_id,
                user_answer=content,
                max_turns=max_turns,
                attachments=attachments,
            )
        if record.status in {"queued", "running", "pausing"}:
            with self._lock:
                message = append_session_message(
                    record,
                    "user",
                    content,
                    delivery="queued",
                    attachments=attachments or [],
                )
                record.pending_user_messages.append({"id": message["id"], "content": content})
                event = {"type": "user_message", "message": message, "delivery": "queued"}
                record.events.append(event)
                record.updated_at = utc_now()
                self._persist_record(record)
            record.event_queue.put(event)
            return record
        if record.status == "paused":
            with self._lock:
                message = append_session_message(
                    record,
                    "user",
                    content,
                    delivery="consumed",
                    attachments=attachments or [],
                )
                record.interaction_status = "answering"
                record.event_queue = queue.Queue()
                chat_run = start_session_run(record, kind="chat", trigger_message_id=message["id"], activate=False)
                record.updated_at = utc_now()
                self._persist_record(record)
            thread = threading.Thread(
                target=self._run_chat,
                args=(record, chat_run["id"], content),
                name=f"bioagent-web-chat-{task_id}",
                daemon=True,
            )
            thread.start()
            return record
        if record.status in {"completed", "failed"}:
            with self._lock:
                message = append_session_message(
                    record,
                    "user",
                    content,
                    delivery="consumed",
                    attachments=attachments or [],
                )
                start_session_run(record, kind="agent", trigger_message_id=message["id"])
                record.status = "queued"
                record.error = ""
                record.final_answer = ""
                record.event_queue = queue.Queue()
                record.updated_at = utc_now()
                self._persist_record(record)
            thread = threading.Thread(
                target=self._run_follow_up,
                args=(record,),
                name=f"bioagent-web-followup-{task_id}",
                daemon=True,
            )
            thread.start()
            return record
        raise ValueError(f"Session {task_id} cannot accept messages while status={record.status}")

    def resume_task(
        self,
        *,
        task_id: str,
        user_answer: str,
        max_turns: int = 20,
        attachments: list[dict[str, Any]] | None = None,
    ) -> TaskRecord:
        record = self.get_task(task_id)
        if record.status not in {"needs_user_input", "paused"}:
            raise ValueError(f"Task {task_id} is not waiting for user input or paused")
        resume_kind = record.status
        question_id = str(record.active_question.get("questionId") or "")
        with self._lock:
            record.pause_event.clear()
            record.status = "running"
            record.error = ""
            record.final_answer = ""
            record.event_queue = queue.Queue()
            append_session_message(
                record,
                "user",
                user_answer,
                delivery="consumed",
                resume=True,
                responseToQuestionId=question_id,
                attachments=attachments or [],
            )
            record.updated_at = utc_now()
            self._persist_record(record)
        thread = threading.Thread(
            target=self._run_resume,
            args=(record, user_answer, max_turns, resume_kind),
            name=f"bioagent-web-resume-{task_id}",
            daemon=True,
        )
        thread.start()
        return record

    def pause_task(self, task_id: str) -> TaskRecord:
        record = self.get_task(task_id)
        event = {"type": "pause_requested", "status": "pausing"}
        with self._lock:
            if record.status not in {"queued", "running"}:
                raise ValueError(f"Task {task_id} is not running")
            record.pause_event.set()
            apply_event_to_record(record, event)
            record.events.append(event)
            record.updated_at = utc_now()
            self._persist_record(record)
        record.event_queue.put(event)
        return record

    def interrupt_current_tool(self, task_id: str) -> TaskRecord:
        record = self.get_task(task_id)
        if not record.run_dir:
            raise ValueError(f"Session {task_id} has no active run workspace")
        stop_active_docker_containers(Path(record.run_dir), include_jobs=False)
        self._apply_event(record, {"type": "tool_interrupt_requested", "status": "interrupting"})
        return record

    def cancel_job(self, task_id: str, job_id: str = "") -> dict[str, Any]:
        record = self.get_task(task_id)
        run_dirs = [Path(path) for path in session_run_dirs(record)]
        if not run_dirs:
            raise ValueError(f"Session {task_id} has no job workspace")
        manager = DockerJobManager(self.config, None, run_dirs[-1], prior_run_dirs=run_dirs[:-1])
        result = manager.cancel_job(job_id)
        self._apply_event(
            record,
            {
                "type": "job_cancelled" if result.get("ok") else "job_cancel_failed",
                "jobId": str(result.get("job_id") or job_id),
                "status": str(result.get("status") or "cancelled"),
                "result": result,
            },
        )
        return result

    def take_pending_messages(self, task_id: str) -> list[Any]:
        record = self.get_task(task_id)
        with self._lock:
            pending = list(record.pending_user_messages)
            record.pending_user_messages.clear()
            context_events = list(record.pending_context_events)
            record.pending_context_events.clear()
            pending_ids = {item["id"] for item in pending}
            for message in record.messages:
                if message.get("id") in pending_ids:
                    message["delivery"] = "consumed"
                    message["consumedAt"] = utc_now()
            if pending:
                record.updated_at = utc_now()
                self._persist_record(record)
        return [
            {"type": "user_message", "messageId": item["id"], "content": item["content"]}
            for item in pending
        ] + context_events

    def scan_job_callbacks(self) -> int:
        with self._lock:
            records = list(self._tasks.values())
        delivered = 0
        for record in records:
            run_dirs = [Path(path) for path in session_run_dirs(record)]
            if not run_dirs:
                continue
            manager = DockerJobManager(self.config, None, run_dirs[-1], prior_run_dirs=run_dirs[:-1])
            for job in manager.list_jobs().get("jobs", []):
                job_id = str(job.get("job_id") or "")
                status = str(job.get("status") or "")
                if status == "running":
                    refreshed = manager.get_job(job_id, refresh=True)
                    status = str(refreshed.get("status") or status)
                    job = refreshed
                elif status in TERMINAL_JOB_STATUSES:
                    job = manager.get_job(job_id)
                if not job_id or status not in TERMINAL_JOB_STATUSES or status == "cancelled":
                    continue
                if str(job.get("callback_state") or "") != "pending":
                    continue
                if record.delivered_job_callbacks.get(job_id) == status:
                    continue
                callback_id = f"callback_{uuid.uuid4().hex[:16]}"
                event = {
                    "type": "job_callback",
                    "callbackId": callback_id,
                    "jobId": job_id,
                    "status": status,
                    "exitCode": job.get("exit_code"),
                    "artifacts": job.get("artifacts") or [],
                    "content": f"Job {job_id} reached terminal status {status}. Inspect its exact results and continue.",
                }
                self._apply_event(record, event)
                manager.mark_callback_delivered(job_id, callback_id)
                record.delivered_job_callbacks[job_id] = status
                delivered += 1
                if record.status in {"completed", "failed"} and not record.interaction_status:
                    self._start_callback_follow_up(record)
        return delivered

    def _job_watch_loop(self) -> None:
        while not self._job_watcher_stop.wait(self._job_watch_interval_s):
            try:
                self.scan_job_callbacks()
            except Exception:
                time.sleep(self._job_watch_interval_s)

    def stop_job_watcher(self) -> None:
        self._job_watcher_stop.set()

    def _start_callback_follow_up(self, record: TaskRecord) -> None:
        with self._lock:
            if record.status not in {"completed", "failed"}:
                return
            start_session_run(record, kind="agent")
            record.status = "queued"
            record.error = ""
            record.final_answer = ""
            record.event_queue = queue.Queue()
            record.updated_at = utc_now()
            self._persist_record(record)
        threading.Thread(
            target=self._run_follow_up,
            args=(record, False),
            name=f"bioagent-callback-{record.id}",
            daemon=True,
        ).start()

    def list_tasks(self) -> list[dict[str, Any]]:
        with self._lock:
            records = list(self._tasks.values())
        result: list[dict[str, Any]] = []
        for record in sorted(records, key=lambda item: item.updated_at, reverse=True):
            result_files = record.result_files()
            result.append({
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
                "resultFileCount": len(flatten_file_nodes(result_files)),
                "resultFileNames": [node["name"] for node in flatten_file_nodes(result_files)],
                "sessionStatus": record.session_status,
                "runStatus": record.status,
                "runCount": len(record.runs),
            })
        return result

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
            kwargs: dict[str, Any] = {
                "data_path": record.data_path,
                "result_dir": result_dir,
                "max_turns": record.max_turns,
            }
            if callable_accepts(self.stream_runner, "pause_requested"):
                kwargs["pause_requested"] = record.pause_event.is_set
            if callable_accepts(self.stream_runner, "take_pending_messages"):
                kwargs["take_pending_messages"] = lambda: self.take_pending_messages(record.id)
            for event in self.stream_runner(record.task, **kwargs):
                self._apply_event(record, event)
        except Exception as exc:
            self._apply_event(record, {"type": "error", "error": str(exc)})
        finally:
            record.event_queue.put(None)

    def _run_resume(self, record: TaskRecord, user_answer: str, max_turns: int, resume_kind: str) -> None:
        self._apply_event(record, {"type": "resume_started", "status": "running", "resume_kind": resume_kind})
        try:
            resume_id = str(record.compute.get("runId") or "")
            kwargs: dict[str, Any] = {"max_turns": max_turns}
            terminal_seen = {"final": False, "run_end": False}

            def emit(event: dict[str, Any]) -> None:
                event_type = str(event.get("type") or "")
                if event_type in terminal_seen:
                    terminal_seen[event_type] = True
                self._apply_event(record, event)

            if callable_accepts(self.resume_runner, "pause_requested"):
                kwargs["pause_requested"] = record.pause_event.is_set
            if callable_accepts(self.resume_runner, "event_sink"):
                kwargs["event_sink"] = emit
            if callable_accepts(self.resume_runner, "take_pending_messages"):
                kwargs["take_pending_messages"] = lambda: self.take_pending_messages(record.id)
            result = self.resume_runner(resume_id, user_answer, **kwargs)
            if not terminal_seen["final"]:
                self._apply_event(
                    record,
                    {
                        "type": "final",
                        "content": str(result.get("final_answer") or ""),
                        "status": str(result.get("status") or "completed"),
                        "result": result,
                    },
                )
            if not terminal_seen["run_end"]:
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

    def _run_follow_up(self, record: TaskRecord, include_latest_user_message: bool = True) -> None:
        self._apply_event(record, {"type": "task_started", "status": "running"})
        try:
            kwargs: dict[str, Any] = {
                "data_path": record.data_path,
                "result_dir": str(self.config.runs_dir / f"web_{record.id}"),
                "max_turns": record.max_turns,
            }
            if callable_accepts(self.stream_runner, "pause_requested"):
                kwargs["pause_requested"] = record.pause_event.is_set
            if callable_accepts(self.stream_runner, "take_pending_messages"):
                kwargs["take_pending_messages"] = lambda: self.take_pending_messages(record.id)
            if callable_accepts(self.stream_runner, "initial_memory_state"):
                kwargs["initial_memory_state"] = session_memory_state(record)
            if include_latest_user_message and callable_accepts(self.stream_runner, "session_message"):
                kwargs["session_message"] = latest_user_message(record)
            if callable_accepts(self.stream_runner, "prior_run_dirs"):
                kwargs["prior_run_dirs"] = session_run_dirs(record)
            prompt = build_session_follow_up_prompt(record)
            for event in self.stream_runner(prompt, **kwargs):
                self._apply_event(record, event)
        except Exception as exc:
            self._apply_event(record, {"type": "error", "error": str(exc)})
        finally:
            record.event_queue.put(None)

    def _run_chat(self, record: TaskRecord, chat_run_id: str, content: str) -> None:
        try:
            context = build_session_chat_context(record)
            context["availableSkills"] = [
                {
                    "skillId": skill.skill_id,
                    "name": skill.name,
                    "description": skill.short_description,
                    "runtime": skill.runtime,
                    "environment": skill.env_profile,
                }
                for skill in list_workflow_skills(self.config.workflows_root)
            ]
            answer = self.chat_runner(context, content)
            with self._lock:
                append_session_message(record, "assistant", answer, delivery="consumed", interaction="chat")
                record.interaction_status = ""
                finish_session_run(record, chat_run_id, "completed")
                event = {"type": "chat_message", "content": answer, "status": "completed"}
                record.events.append(event)
                record.updated_at = utc_now()
                self._persist_record(record)
            record.event_queue.put(event)
        except Exception as exc:
            with self._lock:
                record.interaction_status = ""
                finish_session_run(record, chat_run_id, "failed", error=str(exc))
                append_session_message(record, "assistant", f"Unable to answer: {exc}", delivery="consumed", error=True)
                record.updated_at = utc_now()
                self._persist_record(record)
        finally:
            record.event_queue.put(None)

    def _apply_event(self, record: TaskRecord, event: dict[str, Any]) -> None:
        with self._lock:
            event = envelope_event(
                event,
                session_id=record.id,
                current_run_id=str(record.compute.get("runId") or ""),
                prior_events=record.events,
            )
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
        "plan_state": record.plan_state,
        "active_question": record.active_question,
        "events": record.events,
        "compute": record.compute,
        "memory": record.memory,
        "session_status": record.session_status,
        "interaction_status": record.interaction_status,
        "runs": record.runs,
        "active_run_id": record.active_run_id,
        "pending_user_messages": record.pending_user_messages,
        "pending_context_events": record.pending_context_events,
        "delivered_job_callbacks": record.delivered_job_callbacks,
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
        messages=_normalize_stored_messages(payload.get("messages") or []),
        traces=list(payload.get("traces") or []),
        plan=list(payload.get("plan") or []),
        plan_state=dict(payload.get("plan_state") or {}),
        active_question=dict(payload.get("active_question") or {}),
        events=list(payload.get("events") or []),
        compute=dict(payload.get("compute") or {}),
        memory=dict(payload.get("memory") or {}),
        session_status=str(payload.get("session_status") or "open"),
        interaction_status="",
        runs=list(payload.get("runs") or []),
        active_run_id=str(payload.get("active_run_id") or ""),
        pending_user_messages=list(payload.get("pending_user_messages") or []),
        pending_context_events=list(payload.get("pending_context_events") or []),
        delivered_job_callbacks=dict(payload.get("delivered_job_callbacks") or {}),
        final_answer=str(payload.get("final_answer") or ""),
        run_dir=str(payload.get("run_dir") or ""),
        log_path=str(payload.get("log_path") or ""),
        error=str(payload.get("error") or ""),
    )
    if not record.runs and record.compute.get("runId"):
        run = start_session_run(record, kind="agent", activate=True)
        run.update(
            {
                "agentRunId": str(record.compute.get("runId") or ""),
                "status": record.status,
                "runDir": record.run_dir,
                "logPath": record.log_path,
                "endedAt": record.updated_at if record.status in {"completed", "failed", "paused", "needs_user_input"} else "",
            }
        )
    restored_events: list[dict[str, Any]] = []
    for event in record.events:
        restored_events.append(
            envelope_event(
                event,
                session_id=record.id,
                current_run_id=str(record.compute.get("runId") or ""),
                prior_events=restored_events,
            )
        )
    record.events = restored_events
    for trace in record.traces:
        output = trace.get("output")
        if isinstance(output, dict) and output.get("ok") is False:
            trace["status"] = "failed"
    if record.status in {"running", "pausing"}:
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


def _normalize_stored_messages(messages: Any) -> list[dict[str, Any]]:
    normalized: list[dict[str, Any]] = []
    for value in messages if isinstance(messages, list) else []:
        if not isinstance(value, dict):
            continue
        message = dict(value)
        if message.get("role") == "assistant":
            content, recognized = _legacy_anthropic_text(str(message.get("content") or ""))
            if recognized:
                if not content:
                    continue
                message["content"] = content
        normalized.append(message)
    return normalized


def _legacy_anthropic_text(content: str) -> tuple[str, bool]:
    if not content.lstrip().startswith("["):
        return content, False
    try:
        blocks = ast.literal_eval(content)
    except (SyntaxError, ValueError):
        return content, False
    if not isinstance(blocks, list) or not blocks or not all(isinstance(block, dict) for block in blocks):
        return content, False
    block_types = {str(block.get("type") or "") for block in blocks}
    if not block_types or not block_types.issubset({"text", "tool_use"}):
        return content, False
    text = "\n".join(
        str(block.get("text") or "").strip()
        for block in blocks
        if block.get("type") == "text" and str(block.get("text") or "").strip()
    )
    return text, True


def apply_event_to_record(record: TaskRecord, event: dict[str, Any]) -> None:
    event_type = str(event.get("type", ""))
    if event_type == "task_started":
        if record.status != "pausing":
            record.status = "running"
        ensure_plan_step(record, "Understand request", "completed")
        update_active_session_run(record, status="running")
        return
    if event_type == "resume_started":
        record.status = "running"
        record.active_question = {}
        update_active_session_run(record, status="running")
        if event.get("resume_kind") == "paused":
            ensure_plan_step(record, "Pause task", "completed")
            ensure_plan_step(record, "Continue with user instruction", "completed")
        else:
            ensure_plan_step(record, "Wait for user input", "completed")
            ensure_plan_step(record, "Resume with user input", "completed")
        return
    if event_type == "pause_requested":
        record.status = "pausing"
        update_active_session_run(record, status="pausing")
        ensure_plan_step(record, "Pause task", "running")
        return
    if event_type == "run_start":
        if record.status != "pausing":
            record.status = "running"
        record.run_dir = str(event.get("run_dir") or record.run_dir)
        record.log_path = str(event.get("log_path") or record.log_path)
        record.compute["runId"] = str(event.get("run_id") or record.compute.get("runId") or "")
        update_active_session_run(
            record,
            status="running",
            agentRunId=record.compute["runId"],
            runDir=record.run_dir,
            logPath=record.log_path,
        )
        ensure_plan_step(record, "Start BioAgent run", "completed")
        return
    if event_type == "turn_start":
        record.compute["turn"] = int(event.get("turn") or 0)
        record.compute["maxTurns"] = int(event.get("max_turns") or record.max_turns)
        return
    if event_type == "memory_state":
        memory = event.get("memory")
        if isinstance(memory, dict):
            record.memory = memory
        return
    if event_type == "assistant_message":
        content = str(event.get("content") or "")
        if content:
            append_session_message(record, "assistant", content, delivery="consumed")
        return
    if event_type == "steering_message":
        message_id = str(event.get("messageId") or event.get("message_id") or "")
        message = next((item for item in record.messages if item.get("id") == message_id), None)
        if message is not None:
            message["delivery"] = "consumed"
            message["consumedAt"] = str(event.get("createdAt") or utc_now())
            message["consumedByEventId"] = str(event.get("eventId") or "")
            message["consumedByTurnId"] = str(event.get("turnId") or "")
        return
    if event_type == "tool_call":
        call_id = str(event.get("call_id") or uuid.uuid4().hex[:8])
        tool_name = str(event.get("tool_name") or "tool")
        record.compute["currentTool"] = tool_name
        record.traces.append(
            {
                "id": call_id,
                "eventId": str(event.get("eventId") or ""),
                "toolCallId": str(event.get("toolCallId") or call_id),
                "runId": str(event.get("runId") or event.get("run_id") or ""),
                "turnId": str(event.get("turnId") or ""),
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
        status = str(event.get("status") or ("completed" if event.get("ok") else "failed"))
        if trace is not None:
            trace["status"] = status
            trace["output"] = event.get("result")
            trace["resultEventId"] = str(event.get("eventId") or "")
            trace["updatedAt"] = utc_now()
        record.compute["currentTool"] = ""
        result = event.get("result") if isinstance(event.get("result"), dict) else {}
        plan = result.get("plan") if isinstance(result.get("plan"), dict) else None
        if plan is not None:
            record.plan_state = plan
            record.plan = [dict(item) for item in plan.get("steps") or [] if isinstance(item, dict)]
        if result.get("status") == "needs_user_input" and result.get("question"):
            record.active_question = {
                "questionId": str(result.get("question_id") or ""),
                "interactionType": str(result.get("interaction_type") or "clarification"),
                "question": str(result.get("question") or ""),
                "reason": str(result.get("reason") or ""),
                "options": list(result.get("options") or []),
                "allowFreeText": bool(result.get("allow_free_text", True)),
                "inputType": str(result.get("input_type") or "text"),
                "requiredFields": list(result.get("required_fields") or []),
            }
        ensure_plan_step(record, f"Call {tool_name}", status)
        ensure_plan_step(record, f"Finish {tool_name}", status)
        return
    if event_type == "final":
        content = str(event.get("content") or "")
        record.final_answer = content
        record.status = str(event.get("status") or record.status or "completed")
        update_active_session_run(record, status=record.status)
        if content:
            duplicate_prompt = (
                record.status == "needs_user_input"
                and bool(record.messages)
                and record.messages[-1].get("role") == "assistant"
                and str(record.messages[-1].get("content") or "") == content
            )
            if not duplicate_prompt:
                message = append_session_message(record, "assistant", content, delivery="consumed")
                if record.status == "completed":
                    message["final"] = True
        if record.status == "needs_user_input":
            ensure_plan_step(record, "Wait for user input", "running")
            return
        if record.status == "paused":
            ensure_plan_step(record, "Pause task", "completed")
            return
        ensure_plan_step(record, "Summarize result", "completed")
        return
    if event_type == "run_end":
        record.status = str(event.get("status") or "completed")
        record.run_dir = str(event.get("run_dir") or record.run_dir)
        record.log_path = str(event.get("log_path") or record.log_path)
        record.compute["runId"] = str(event.get("run_id") or record.compute.get("runId") or "")
        update_active_session_run(
            record,
            status=record.status,
            agentRunId=record.compute["runId"],
            runDir=record.run_dir,
            logPath=record.log_path,
            endedAt=utc_now() if record.status in {"completed", "failed", "paused", "needs_user_input"} else "",
        )
        result = event.get("result") if isinstance(event.get("result"), dict) else {}
        if not record.final_answer:
            record.final_answer = str(result.get("final_answer") or "")
        if record.status == "needs_user_input":
            return
        if record.status == "paused":
            ensure_plan_step(record, "Pause task", "completed")
            return
        ensure_plan_step(record, "Finish run", "completed" if record.status == "completed" else record.status)
        return
    if event_type == "job_callback":
        job_id = str(event.get("jobId") or event.get("job_id") or "")
        callback_id = str(event.get("callbackId") or event.get("callback_id") or "")
        if job_id:
            record.delivered_job_callbacks[job_id] = str(event.get("status") or "unknown")
        if callback_id and not any(item.get("callbackId") == callback_id for item in record.pending_context_events):
            record.pending_context_events.append(dict(event))
        ensure_plan_step(record, f"Job {job_id} {event.get('status', 'updated')}", "completed")
        return
    if event_type == "error":
        record.status = "failed"
        record.error = str(event.get("error") or "Unknown error")
        record.final_answer = f"Run failed: {record.error}"
        append_session_message(record, "assistant", record.final_answer, delivery="consumed", error=True)
        update_active_session_run(record, status="failed", endedAt=utc_now(), error=record.error)
        ensure_plan_step(record, "Handle error", "failed")


def append_session_message(
    record: TaskRecord,
    role: str,
    content: str,
    *,
    delivery: str,
    **metadata: Any,
) -> dict[str, Any]:
    message = {
        "id": uuid.uuid4().hex[:12],
        "role": role,
        "content": content,
        "createdAt": utc_now(),
        "delivery": delivery,
        **metadata,
    }
    record.messages.append(message)
    return message


def start_session_run(
    record: TaskRecord,
    *,
    kind: str,
    trigger_message_id: str = "",
    activate: bool = True,
) -> dict[str, Any]:
    run = {
        "id": uuid.uuid4().hex[:12],
        "kind": kind,
        "status": "queued" if kind == "agent" else "running",
        "triggerMessageId": trigger_message_id,
        "agentRunId": "",
        "runDir": "",
        "logPath": "",
        "startedAt": utc_now(),
        "endedAt": "",
        "error": "",
    }
    record.runs.append(run)
    if activate:
        record.active_run_id = run["id"]
    return run


def update_active_session_run(record: TaskRecord, **updates: Any) -> None:
    run = next((item for item in record.runs if item.get("id") == record.active_run_id), None)
    if run is not None:
        run.update(updates)


def finish_session_run(record: TaskRecord, run_id: str, status: str, *, error: str = "") -> None:
    run = next((item for item in record.runs if item.get("id") == run_id), None)
    if run is None:
        return
    run.update({"status": status, "endedAt": utc_now(), "error": error})


def build_session_chat_context(record: TaskRecord) -> dict[str, Any]:
    task_state = record.memory.get("taskState") if isinstance(record.memory, dict) else {}
    return {
        "sessionId": record.id,
        "originalTask": record.task,
        "runStatus": record.status,
        "taskState": task_state or {},
        "recentMessages": record.messages[-12:],
        "runId": record.compute.get("runId", ""),
        "resultRoot": record.compute_snapshot().get("resultRoot", ""),
    }


def build_session_follow_up_prompt(record: TaskRecord) -> str:
    recent = record.messages[-12:]
    transcript = "\n".join(
        f"{str(message.get('role') or 'unknown').upper()}: {str(message.get('content') or '')[:3000]}"
        for message in recent
    )
    task_state = record.memory.get("taskState") if isinstance(record.memory, dict) else {}
    failed_trace = next((trace for trace in reversed(record.traces) if trace.get("status") == "failed"), None)
    failure_context = json.dumps(failed_trace or {}, ensure_ascii=False, default=str)[:5000]
    return (
        "Continue an existing BioAgent conversation session. Current observations override prior messages.\n"
        f"Original task:\n{record.task}\n\n"
        f"Compact task state:\n{json.dumps(task_state or {}, ensure_ascii=False, default=str)}\n\n"
        f"Previous run workspace:\n{record.run_dir or '(not available)'}\n"
        f"Previous run log:\n{record.log_path or '(not available)'}\n\n"
        f"Latest failed tool trace:\n{failure_context or '(none)'}\n\n"
        f"Recent conversation:\n{transcript}\n\n"
        "Respond to the latest USER message. Use tools when needed and reuse verified artifacts from the existing session."
    )


def session_memory_state(record: TaskRecord) -> dict[str, Any]:
    memory = record.memory if isinstance(record.memory, dict) else {}
    task_state = memory.get("taskState") if isinstance(memory.get("taskState"), dict) else {}
    return {
        "task_state": task_state,
        "execution_outcome": memory.get("executionOutcome"),
        "execution_tool": str(memory.get("executionTool") or "execution"),
    }


def latest_user_message(record: TaskRecord) -> str:
    for message in reversed(record.messages):
        if message.get("role") == "user":
            return str(message.get("content") or "")
    return ""


def session_run_dirs(record: TaskRecord) -> list[str]:
    result: list[str] = []
    for run in record.runs:
        run_dir = str(run.get("runDir") or "").strip()
        if run_dir and run_dir not in result:
            result.append(run_dir)
    return result


def ensure_plan_step(record: TaskRecord, title: str, status: str) -> None:
    if record.plan_state:
        return
    existing = next((item for item in record.plan if item["title"] == title), None)
    if existing:
        existing["status"] = status
        return
    record.plan.append({"id": uuid.uuid4().hex[:8], "title": title, "status": status})


def callable_accepts(function: Callable[..., Any], parameter: str) -> bool:
    try:
        signature = inspect.signature(function)
    except (TypeError, ValueError):
        return False
    return parameter in signature.parameters or any(
        item.kind == inspect.Parameter.VAR_KEYWORD for item in signature.parameters.values()
    )


def stop_active_docker_containers(run_dir: Path, *, include_jobs: bool = True) -> None:
    for name in (".docker-python.cid", ".docker-r.cid"):
        cidfile = run_dir / name
        if not cidfile.is_file():
            continue
        container_id = cidfile.read_text(encoding="utf-8").strip()
        if not container_id:
            continue
        subprocess.run(
            ["docker", "rm", "-f", container_id],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=30,
            check=False,
        )
    if include_jobs:
        cancel_run_jobs(run_dir)


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


def build_merged_file_tree(roots: list[Path]) -> list[dict[str, Any]]:
    merged: dict[str, dict[str, Any]] = {}
    for root in roots:
        for node in build_file_tree(root):
            name = str(node.get("name") or "")
            existing = merged.get(name)
            if existing and existing.get("type") == node.get("type") == "directory":
                existing["children"] = merge_file_nodes(existing.get("children", []), node.get("children", []))
            else:
                merged[name] = node
    return sorted(merged.values(), key=file_node_sort_key)


def merge_file_nodes(existing: list[dict[str, Any]], incoming: list[dict[str, Any]]) -> list[dict[str, Any]]:
    merged = {str(node.get("name") or ""): node for node in existing}
    for node in incoming:
        name = str(node.get("name") or "")
        current = merged.get(name)
        if current and current.get("type") == node.get("type") == "directory":
            current["children"] = merge_file_nodes(current.get("children", []), node.get("children", []))
        else:
            merged[name] = node
    return sorted(merged.values(), key=file_node_sort_key)


def file_node_sort_key(node: dict[str, Any]) -> tuple[int, int, str]:
    if node.get("type") == "directory":
        return (0, 0, str(node.get("name") or "").lower())
    name = str(node.get("name") or "")
    return (1, result_file_rank(name), name.lower())


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


def build_tool_activity(traces: list[dict[str, Any]]) -> list[dict[str, Any]]:
    order: list[str] = []
    activity: dict[str, dict[str, Any]] = {}
    for index, trace in enumerate(traces):
        key = str(trace.get("toolCallId") or trace.get("id") or trace.get("eventId") or f"trace-{index}")
        if key not in activity:
            order.append(key)
        activity[key] = {
            "id": key,
            "traceId": str(trace.get("id") or key),
            "eventId": str(trace.get("eventId") or ""),
            "toolCallId": str(trace.get("toolCallId") or key),
            "toolName": str(trace.get("toolName") or "tool"),
            "status": str(trace.get("status") or "pending"),
            "runId": str(trace.get("runId") or ""),
            "turnId": str(trace.get("turnId") or ""),
            "createdAt": str(trace.get("createdAt") or ""),
            "updatedAt": str(trace.get("updatedAt") or ""),
        }
    return [activity[key] for key in order]


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
