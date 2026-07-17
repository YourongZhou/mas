from __future__ import annotations

import json
import time
import uuid
from datetime import datetime, timezone
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from langchain_core.messages import HumanMessage, SystemMessage, ToolMessage

from .config import AgentConfig
from .llm import build_llm, message_text, runtime_summary
from .logging_utils import RunLogger
from .memory import build_memory_harness
from .run_state import save_pending_state
from .tools import build_tools
from .tools.jobs import attach_job_provenance
from .tools.planning import PlanStore


CONTEXT_LOOKUP_TOOLS = {
    "list_files",
    "read_file",
    "glob_search",
    "grep_text",
    "read_skill_script",
    "inspect_skill_script_symbols",
    "inspect_skill_function",
}

BUDGET_EXEMPT_TOOLS = {
    *CONTEXT_LOOKUP_TOOLS,
    "list_workflow_skills",
    "inspect_workflow_skill",
    "inspect_image_catalog",
    "validate_script",
    "list_jobs",
    "get_job",
    "poll_job",
    "tail_job",
    "inspect_artifact",
    "finish_task",
}


SYSTEM_PROMPT = """You are a single main-agent loop for bioinformatics analysis.

You own the whole task: understand the request, inspect relevant files or workflow skills, choose a Docker profile, generate code, execute it with tools, repair failures, and finish with a concise report.

Important operating rules:
- Use tools when you need filesystem context, workflow metadata, Docker profile details, code execution, or artifact verification.
- For a complex task with at least three dependent stages, use propose_plan before execution and keep it current with update_plan. Plans describe the model's intended scientific work; they are not a list of tool calls.
- Request plan approval only when ambiguity, cost, or a consequential scientific choice requires the user. Simple or already well-specified tasks should proceed without an approval ceremony.
- When clarification is required, use request_user_input with concise options when possible. Never guess a consequential missing input.
- Choose the lightest sufficient route, similar to a CLI coding agent:
  1. Answer directly for conceptual questions that need no files or execution.
  2. Use file/search tools for inspection-only tasks.
  3. When a mas_2 workflow skill matches, inspect the skill and relevant script symbols/functions, then write the analysis yourself into a persistent run-workspace script.
  4. For ad hoc analysis, inspect files as needed and write the analysis yourself into a persistent run-workspace script.
- For analysis execution, use this default lifecycle: write_run_file -> validate_script -> start_job -> poll_job -> inspect_artifact -> finish_task.
- Keep code in /work/scripts and refer to it by path. Use edit_run_file for minimal repairs after reading the exact error. Do not resend a complete script in later tool calls.
- execute_python and execute_r are compatibility fallbacks for tiny one-off probes only. Do not use them for a full analysis when write_run_file + start_job can do the work.
- start_job returns immediately. For long jobs, call poll_job with wait_s=60..300. Use tail_job only when progress or an error needs inspection; do not burn turns with rapid polling.
- Job ids are issued only by start_job and persisted by the harness. Never invent or derive a job id from a turn number. Leave job_id empty to use the active job, or call list_jobs/get_job when uncertain.
- Jobs remain registered across follow-up runs in the same conversation session. A new run workspace does not mean the prior container died. Continue must inspect or poll the existing session job before considering a restart.
- Only one job may run at a time in a session. Use list_jobs, get_job, poll_job, tail_job, or cancel_job to manage it before starting another.
- poll_job and tail_job return recent stdout/stderr while running; completed and failed logs remain persisted in job metadata. Read the exact logs before repairing code.
- A successful start_job only means the container started. Analysis is successful only when poll_job returns status=completed and exit_code=0.
- Before reporting success, inspect the important output files with inspect_artifact and call finish_task with those evidence_ids. If an artifact changes, inspect it again.
- Prefer existing mas_2 workflow skills and Docker profiles when they match the task, but do not force a skill when the task is simpler without one.
- Gene lookup, enrichment, and cell-type annotation are not hardcoded tools. Use a matching workflow skill when one exists, such as functional-enrichment-from-degs or scrnaseq-scanpy-core-analysis; otherwise write a transparent Python/R script with the execution tools.
- Do not run marker-gene, differential-expression, trajectory, interaction, or report-generation steps unless the user explicitly asks for them. For a core single-cell demo, stop after QC, normalization, HVG, PCA, neighbors, Leiden, UMAP, summary outputs, and processed h5ad.
- Do not call sc.tl.rank_genes_groups or sc.pl.rank_genes_groups unless marker-gene or differential-expression analysis was explicitly requested.
- Manage context progressively: list skills first, inspect only the matching skill, then use read_skill_script, inspect_skill_script_symbols, or inspect_skill_function for the exact scripts/functions you need.
- Do not repeat a context lookup with the same tool and identical arguments immediately; use the existing result or change the lookup scope.
- Do not invent skill script paths or helper function names. Use only exact script paths returned by inspect_workflow_skill or list_files, and verify unknown function signatures before calling them.
- For a primary data path ending in .h5ad, use scanpy's sc.read_h5ad directly unless a verified skill function is specifically needed.
- Do not inspect every script in a workflow skill. Inspect only the few scripts/functions needed for the task, usually 3-6, then write and execute code.
- Hard guardrail: after inspect_workflow_skill for a workflow skill, do not call start_job, execute_python, or execute_r until you have inspected at least one relevant skill script or function with read_skill_script, inspect_skill_script_symbols, or inspect_skill_function.
- For Python single-cell analysis prefer py-scverse-v1. For R/Seurat analysis prefer r-singlecell-v1. For general Python analysis prefer py-general-v1. For Bioconductor/DESeq2/survival prefer r-bioc-v1.
- Do not claim success until execution tools return ok=true or you clearly explain why execution could not be completed.
- Generated scripts must write outputs under /work/outputs when running in Docker.
- /work/outputs is the canonical result directory and the frontend displays it directly. Never copy outputs to a host result_dir after execution.
- The repo is mounted read-only at /repo and the run workspace is mounted at /work during Docker execution.
- If an execution fails, repair based on the exact stdout/stderr and retry with a minimal change. Before guessing a helper function argument, call inspect_skill_function or inspect_skill_script_symbols to verify the signature.
- Execution attempts are globally limited and persisted by the harness, so each start_job, execute_python, or execute_r call should be intentional. validate_script, poll_job, and artifact inspection do not consume attempts.
- max_turns limits state-changing decision steps, not ordinary inspection, validation, artifact verification, or asynchronous job waiting. Do not rush or restart a healthy job because the displayed decision step is near its limit.
- End with final output paths, key findings, and any unresolved blockers.
"""


def messages_for_model(messages: list[Any], config: AgentConfig, *, keep_recent_tool_messages: int = 2) -> list[Any]:
    tool_indices = [idx for idx, message in enumerate(messages) if getattr(message, "type", "") == "tool"]
    keep_indices = set(tool_indices[-keep_recent_tool_messages:])
    tool_names = {
        str(call.get("id", "")): str(call.get("name", ""))
        for message in messages
        if getattr(message, "type", "") == "ai"
        for call in (getattr(message, "tool_calls", None) or [])
    }
    compacted: list[Any] = []
    for idx, message in enumerate(messages):
        if getattr(message, "type", "") != "tool":
            compacted.append(_compact_ai_tool_calls(message))
            continue
        if idx in keep_indices:
            compacted.append(_trim_tool_message(message, max_chars=min(2200, config.max_tool_result_chars)))
            continue
        tool_name = tool_names.get(str(getattr(message, "tool_call_id", "")), "")
        max_chars = 900 if tool_name == "inspect_skill_script_symbols" else 350
        compacted.append(
            _compact_tool_message(
                message,
                max_chars=min(max_chars, config.max_tool_result_chars),
                tool_name=tool_name,
            )
        )
    system_messages = [message for message in compacted if getattr(message, "type", "") == "system"]
    non_system_messages = [message for message in compacted if getattr(message, "type", "") != "system"]
    return [*system_messages, *non_system_messages]


def _compact_ai_tool_calls(message: Any, *, max_code_chars: int = 1400) -> Any:
    if getattr(message, "type", "") != "ai" or not getattr(message, "tool_calls", None):
        return message
    changed = False
    tool_calls: list[dict[str, Any]] = []
    for call in message.tool_calls:
        copied = dict(call)
        args = dict(copied.get("args") or {})
        code = args.get("code")
        if copied.get("name") in {"execute_python", "execute_r"} and isinstance(code, str) and len(code) > max_code_chars:
            marker = (
                "\n[execution code compacted after it was written to scripts/code.py; "
                "use read_file on the script_path from the tool result for exact code]\n"
            )
            head_chars = max(1, (max_code_chars - len(marker)) // 2)
            tail_chars = max(1, max_code_chars - len(marker) - head_chars)
            args["code"] = code[:head_chars].rstrip() + marker + code[-tail_chars:].lstrip()
            copied["args"] = args
            changed = True
        content = args.get("content")
        if copied.get("name") == "write_run_file" and isinstance(content, str) and len(content) > max_code_chars:
            args["content"] = (
                f"[file content compacted after write; {len(content)} characters were written. "
                "Use read_file on the returned path for exact content.]"
            )
            copied["args"] = args
            changed = True
        tool_calls.append(copied)
    if not changed:
        return message
    additional_kwargs = dict(getattr(message, "additional_kwargs", {}) or {})
    additional_kwargs.pop("tool_calls", None)
    return message.model_copy(update={"tool_calls": tool_calls, "additional_kwargs": additional_kwargs})


def messages_with_turn_guidance(
    messages: list[Any],
    config: AgentConfig,
    *,
    turn: int,
    max_turns: int,
    memory_context: str = "",
) -> list[Any]:
    remaining_after_this = max_turns - turn
    guidance = (
        f"Decision budget: this is decision step {turn} of {max_turns}; "
        f"{remaining_after_this} decision step(s) remain after this one. "
        "Filesystem/skill inspection, script validation, job status and waiting tools, artifact inspection, and finish_task "
        "are part of the current decision step and do not consume another step. "
        "When the decision budget is exhausted, only observe an existing job, verify artifacts, finish the task, "
        "or provide the best grounded final answer; do not write, edit, restart, cancel, or launch new execution."
    )
    compacted = messages_for_model(messages, config)
    if memory_context:
        memory_message = SystemMessage(content=memory_context)
        if compacted and getattr(compacted[0], "type", "") == "system":
            compacted = [compacted[0], memory_message, *compacted[1:]]
        else:
            compacted = [memory_message, *compacted]
    return [*compacted, HumanMessage(content=guidance)]


def _trim_tool_message(message: Any, *, max_chars: int) -> ToolMessage:
    content = str(getattr(message, "content", ""))
    if len(content) <= max_chars:
        return message
    tool_call_id = str(getattr(message, "tool_call_id", ""))
    marker = "\n[recent tool output truncated; middle omitted; request a narrower file/skill inspection if needed]\n"
    head_chars = max(1, (max_chars - len(marker)) // 2)
    tail_chars = max(1, max_chars - len(marker) - head_chars)
    return ToolMessage(
        content=content[:head_chars].rstrip() + marker + content[-tail_chars:].lstrip(),
        tool_call_id=tool_call_id,
    )


def _compact_tool_message(message: Any, *, max_chars: int, tool_name: str = "") -> ToolMessage:
    content = str(getattr(message, "content", ""))
    tool_call_id = str(getattr(message, "tool_call_id", ""))
    if tool_name == "inspect_skill_script_symbols":
        try:
            result = json.loads(content)
        except (TypeError, json.JSONDecodeError):
            result = None
        if isinstance(result, dict):
            functions = result.get("functions") or []
            names = [str(item.get("name")) for item in functions if isinstance(item, dict) and item.get("name")]
            header = (
                "Verified skill symbols (preserved from earlier tool output):\n"
                f"skill={result.get('skill_id', '')} script={result.get('script_path', '')}\n"
                f"functions={', '.join(names) or '(none)'}"
            )
            signatures = [
                f"- {item.get('signature')}"
                for item in functions
                if isinstance(item, dict) and item.get("signature")
            ]
            compact = header
            for signature in signatures:
                candidate = f"{compact}\n{signature}"
                if len(candidate) > max_chars:
                    break
                compact = candidate
            return ToolMessage(content=compact[:max_chars], tool_call_id=tool_call_id)
    if tool_name == "inspect_artifact":
        try:
            result = json.loads(content)
        except (TypeError, json.JSONDecodeError):
            result = None
        if isinstance(result, dict):
            compact = {
                "ok": result.get("ok"),
                "evidence_id": result.get("evidence_id"),
                "path": result.get("path"),
                "facts": result.get("facts"),
            }
            return ToolMessage(content=_safe_json(compact, pretty=False)[:max_chars], tool_call_id=tool_call_id)
    preview = content[:max_chars].replace("\n", "\\n")
    return ToolMessage(
        content=(
            "[tool output compressed to keep model context small; "
            "call a file/skill inspection tool again if exact details are needed]"
            f"\npreview: {preview}"
        ),
        tool_call_id=tool_call_id,
    )


def tool_message_text(tool_name: str, result: Any) -> str:
    if tool_name in {"execute_python", "execute_r", "start_job", "poll_job"} and isinstance(result, dict) and not result.get("ok"):
        compact = {
            "ok": False,
            "status": result.get("status"),
            "exit_code": result.get("exit_code"),
            "error_reason": result.get("error_reason") or result.get("error") or "",
            "job_id": result.get("job_id"),
            "active_job_id": result.get("active_job_id"),
            "available_job_ids": result.get("available_job_ids"),
            "stdout_tail": str(result.get("stdout_tail") or result.get("stdout") or "")[-1200:],
            "stderr_tail": str(result.get("stderr_tail") or result.get("stderr") or "")[-1200:],
            "script_path": result.get("script_path"),
        }
        return _safe_json(compact, pretty=False)
    return _safe_json(result, pretty=False)


@dataclass
class AgentRunResult:
    final_answer: str
    run_dir: str
    log_path: str
    turns: int
    status: str = "completed"
    pending_question: str = ""
    resume_id: str = ""
    pending_state_path: str = ""
    pending_interaction: dict[str, Any] | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "final_answer": self.final_answer,
            "run_dir": self.run_dir,
            "log_path": self.log_path,
            "turns": self.turns,
            "status": self.status,
            "pending_question": self.pending_question,
            "resume_id": self.resume_id,
            "pending_state_path": self.pending_state_path,
            "pending_interaction": self.pending_interaction or {},
        }


class BioAgent:
    def __init__(
        self,
        *,
        config: AgentConfig,
        logger: RunLogger,
        run_dir: Path,
        prior_run_dirs: list[Path] | None = None,
    ) -> None:
        self.config = config
        self.logger = logger
        self.run_dir = run_dir
        self.prior_run_dirs = [Path(path) for path in prior_run_dirs or []]
        self.memory_harness = build_memory_harness(config)
        self.plan_store = PlanStore(run_dir)
        self.tools = build_tools(
            config,
            logger,
            run_dir,
            memory_harness=self.memory_harness,
            prior_run_dirs=self.prior_run_dirs,
        )
        self.tool_map = {tool.name: tool for tool in self.tools}
        self.llm = build_llm(config).bind_tools(self.tools)

    def _initial_user_message(self, task: str, data_path: str | None, result_dir: str | None) -> str:
        parts = [
            "User bioinformatics task:",
            task.strip(),
            "",
            f"Host repo root: {self.config.repo_root}",
            f"Run workspace: {self.run_dir}",
            f"Docker /repo mount points to host repo root: {self.config.repo_root}",
            f"Docker /work mount points to this run workspace: {self.run_dir}",
            f"Global execution attempt budget: {self.config.max_execution_attempts}",
        ]
        if data_path:
            normalized_data_path = data_path.replace("\\", "/")
            parts.append(f"Primary data path: {data_path}")
            if not Path(normalized_data_path).is_absolute():
                parts.append(f"Primary data path inside Docker: /repo/{normalized_data_path.lstrip('/')}")
            parts.append("Use this exact primary data path for the analysis; do not search for or invent alternate data paths unless this path is missing or unreadable.")
            parts.append("Do not call list_files only to locate the primary data path; use execute_python with sc.read_h5ad for a lightweight preview if needed.")
            parts.append("Inside Docker, host repo files are visible under /repo. If the data path is inside the repo, translate it to /repo-relative path.")
        if result_dir:
            parts.append(f"Frontend result display hint: {result_dir}")
            parts.append("Do not copy files there. Write all artifacts to /work/outputs; the frontend reads that directory directly.")
        return "\n".join(parts)

    def _emit_event(
        self,
        event_sink: Callable[[dict[str, Any]], None] | None,
        event_type: str,
        **payload: Any,
    ) -> dict[str, Any] | None:
        if event_sink is None:
            return None
        event = {
            "type": event_type,
            "eventId": f"event_{uuid.uuid4().hex[:16]}",
            "createdAt": datetime.now(timezone.utc).isoformat(),
            "run_id": self.logger.run_id,
            "run_dir": str(self.run_dir),
            "log_path": str(self.logger.path),
        }
        event.update(payload)
        event_sink(event)
        return event

    def _emit_memory_state(self, event_sink: Callable[[dict[str, Any]], None] | None) -> None:
        self._emit_event(event_sink, "memory_state", memory=self.memory_harness.public_snapshot())

    def run(
        self,
        task: str,
        *,
        data_path: str | None = None,
        result_dir: str | None = None,
        max_turns: int = 20,
        initial_messages: list[Any] | None = None,
        resume_answer: str | None = None,
        initial_memory_state: dict[str, Any] | None = None,
        event_sink: Callable[[dict[str, Any]], None] | None = None,
        pause_requested: Callable[[], bool] | None = None,
        take_pending_messages: Callable[[], list[Any]] | None = None,
    ) -> AgentRunResult:
        self.memory_harness.start_task(
            task,
            data_path=data_path,
            result_dir=result_dir,
            restored_state=initial_memory_state,
        )
        if resume_answer is not None:
            self.memory_harness.apply_user_input(resume_answer)
        self.logger.node("主 Agent 循环启动", "加载运行配置、工具清单，并开始理解用户任务。")
        self.logger.preview("Runtime", runtime_summary(self.config))
        self.logger.preview("Memory", self.memory_harness.runtime_summary())
        self.logger.progress("可用工具", ", ".join(self.tool_map))
        self._emit_event(
            event_sink,
            "run_start",
            tools=list(self.tool_map),
            runtime=runtime_summary(self.config),
            memory=self.memory_harness.runtime_summary(),
        )
        self._emit_memory_state(event_sink)

        if initial_messages is None:
            messages: list[Any] = [
                SystemMessage(content=SYSTEM_PROMPT),
                HumanMessage(content=self._initial_user_message(task, data_path, result_dir)),
            ]
        else:
            messages = list(initial_messages)
            if resume_answer is not None:
                messages.append(HumanMessage(content=f"User input to continue this run:\n{resume_answer}"))
        final_answer = ""
        status = "completed"
        pending_question = ""
        pending_state_path = ""
        pending_interaction: dict[str, Any] = {}
        workflow_skill_inspected = False
        workflow_skill_detail_seen = False
        last_context_lookup_key: tuple[str, str] | None = None
        persistent_job_started, completed_job_needs_finish = _persistent_job_state(
            [self.run_dir, *self.prior_run_dirs]
        )

        decision_turns = 0
        model_calls = 0
        while True:
            turn = min(max_turns, decision_turns + 1)
            model_calls += 1
            self._append_pending_messages(messages, take_pending_messages, event_sink=event_sink, turn=turn)
            if _pause_is_requested(pause_requested):
                status = "paused"
                pending_question = "Task paused. Add instructions to continue."
                pending_state_path = self._save_pause_checkpoint(messages, pending_question, event_sink=event_sink)
                final_answer = pending_question
                break
            self.logger.node(f"主 Agent 第 {turn} 轮推理", "判断是否需要读取文件、检查 Skill、执行代码或直接回答。")
            self._emit_event(
                event_sink,
                "turn_start",
                turn=turn,
                max_turns=max_turns,
                model_call=model_calls,
                decision_turns_used=decision_turns,
            )
            started = time.perf_counter()
            self.logger.progress("调用主模型", f"turn={turn}")
            self._emit_event(
                event_sink,
                "model_start",
                turn=turn,
                max_turns=max_turns,
                model_call=model_calls,
                decision_turns_used=decision_turns,
            )
            response = self.llm.invoke(
                messages_with_turn_guidance(
                    messages,
                    self.config,
                    turn=turn,
                    max_turns=max_turns,
                    memory_context="\n\n".join(
                        value for value in (self.memory_harness.context_text(), self.plan_store.context_text()) if value
                    ),
                )
            )
            elapsed = time.perf_counter() - started
            messages.append(response)
            response_text = message_text(response)
            self.logger.progress("主模型返回", f"turn={turn} elapsed={elapsed:.2f}s tool_calls={len(response.tool_calls or [])}")
            self.logger.log(f"LLM_RESPONSE turn={turn} elapsed={elapsed:.2f}s tool_calls={len(response.tool_calls or [])}")
            self._emit_event(
                event_sink,
                "model_end",
                turn=turn,
                elapsed_s=elapsed,
                tool_call_count=len(response.tool_calls or []),
            )
            if response_text:
                self.logger.preview("ASSISTANT_CONTENT", response_text, max_chars=6000)
                self._emit_event(event_sink, "assistant_message", turn=turn, content=response_text)

            invalid_tool_calls = getattr(response, "invalid_tool_calls", None) or []
            if invalid_tool_calls:
                invalid_detail = _invalid_tool_call_detail(invalid_tool_calls[0])
                self.logger.preview("INVALID_TOOL_CALL", invalid_detail, max_chars=4000)
                self._emit_event(
                    event_sink,
                    "model_error",
                    turn=turn,
                    error_type="invalid_tool_call",
                    content=invalid_detail,
                )
                decision_turns += 1
                if decision_turns < max_turns:
                    messages.append(
                        HumanMessage(
                            content=(
                                "Your previous tool call could not be parsed. Use the exact tool schema and send one valid JSON "
                                f"arguments object. Parsing error:\n{invalid_detail}"
                            )
                        )
                    )
                    continue
                status = "failed"
                final_answer = f"Model tool call parsing failed on the final turn.\n\n{invalid_detail}"
                break

            tool_calls = response.tool_calls or []
            steering_messages = self._append_pending_messages(
                messages,
                take_pending_messages,
                event_sink=event_sink,
                turn=turn,
            )
            if steering_messages:
                if tool_calls:
                    self._cancel_unstarted_tool_calls(
                        messages,
                        tool_calls,
                        turn=turn,
                        event_sink=event_sink,
                        reason="Cancelled before execution because a new session message changed the active context.",
                    )
                continue
            if _pause_is_requested(pause_requested):
                self._cancel_unstarted_tool_calls(messages, tool_calls, turn=turn, event_sink=event_sink)
                status = "paused"
                pending_question = "Task paused. Add instructions to continue."
                pending_state_path = self._save_pause_checkpoint(messages, pending_question, event_sink=event_sink)
                final_answer = pending_question
                break
            if not tool_calls:
                if completed_job_needs_finish:
                    decision_turns += 1
                    self.logger.progress("作业已完成但产物尚未验证", "要求 inspect_artifact + finish_task 后再收束。")
                    messages.append(
                        HumanMessage(
                            content=(
                                "The analysis job completed, but the task is not grounded yet. Use the exact artifact paths "
                                "returned by poll_job, inspect the important output files with inspect_artifact, then call "
                                "finish_task with the returned evidence_ids. Do not guess filenames or report unverified numbers."
                            )
                        )
                    )
                    continue
                if _is_unusable_final_answer(response_text) or _looks_like_continuation_without_tool(response_text):
                    decision_turns += 1
                    if decision_turns < max_turns:
                        self.logger.progress("主模型未给出可结束回答", "继续要求其修复失败或给出明确最终回答。")
                        messages.append(
                            HumanMessage(
                                content=(
                                    "Your last response did not finish the task and did not call a tool. Continue now: if the latest execution failed, "
                                    "repair the script using the exact error; otherwise provide a concise final answer "
                                    "with output paths and key findings. Do not add new optional analysis steps."
                                )
                            )
                        )
                        continue
                    status = "failed"
                    final_answer = "The model reached the turn limit without a usable final answer or valid tool call."
                    break
                decision_turns += 1
                self.logger.progress("主 Agent 决定结束", "本轮没有工具调用，使用模型文本作为最终回答。")
                final_answer = response_text
                break
            consumes_decision = _tool_calls_consume_decision(tool_calls)
            if consumes_decision and decision_turns >= max_turns:
                self.logger.progress(
                    "决策预算已用尽，拒绝新的有副作用工具调用",
                    f"requested={','.join(str(call.get('name')) for call in tool_calls)}",
                )
                final_answer = response_text or self._final_answer_at_turn_limit(max_turns=max_turns, requested_tool_calls=tool_calls)
                status = "failed"
                break

            stop_after_tools = False
            for call in tool_calls:
                name = call.get("name")
                args = call.get("args") or {}
                call_id = str(call.get("id") or f"toolcall_{uuid.uuid4().hex[:16]}")
                self.logger.node(f"工具调用：{name}", f"call_id={call_id}")
                self.logger.log(f"TOOL_CALL name={name} id={call_id}")
                self.logger.preview("工具入参", _safe_json(args), max_chars=4000)
                self.memory_harness.observe_tool_call(str(name), args)
                self._emit_memory_state(event_sink)
                tool_call_event = self._emit_event(
                    event_sink,
                    "tool_call",
                    turn=turn,
                    tool_name=str(name),
                    call_id=call_id,
                    toolCallId=call_id,
                    args=args,
                )
                context_lookup_key = _context_lookup_call_key(str(name), args)
                duplicate_context_lookup = (
                    context_lookup_key is not None and context_lookup_key == last_context_lookup_key
                )
                last_context_lookup_key = context_lookup_key
                tool = self.tool_map.get(str(name))
                if tool is None:
                    result: Any = {"error": f"Unknown tool: {name}"}
                    self.logger.error_reason(f"Unknown tool: {name}")
                elif duplicate_context_lookup:
                    result = _duplicate_context_lookup(str(name), args)
                    self.logger.progress("跳过连续重复的上下文查阅", f"name={name}")
                    self.logger.log(f"TOOL_GUARD name={name} reason=duplicate_context_lookup")
                elif _should_block_execution_for_skill_detail(str(name), args, workflow_skill_inspected, workflow_skill_detail_seen):
                    result = _workflow_skill_detail_required(str(name))
                    self.logger.progress("执行前需要读取 Skill 脚本", f"name={name}")
                    self.logger.log(f"TOOL_GUARD name={name} reason=workflow_skill_detail_required")
                elif persistent_job_started and str(name) in {"execute_python", "execute_r"}:
                    result = {
                        "ok": False,
                        "status": "persistent_job_lifecycle_required",
                        "error": (
                            "A persistent job lifecycle is already active for this task. Do not switch back to inline code execution. "
                            "Use edit_run_file and start_job for repairs, or inspect_artifact and finish_task for completed outputs."
                        ),
                    }
                    self.logger.progress("拒绝退回整段代码执行", f"name={name}")
                    self.logger.log(f"TOOL_GUARD name={name} reason=persistent_job_lifecycle_required")
                else:
                    try:
                        tool_started = time.perf_counter()
                        self.logger.progress("执行工具", f"name={name}")
                        result = tool.invoke(args)
                        tool_elapsed = time.perf_counter() - tool_started
                        self.logger.progress("工具返回", f"name={name} elapsed={tool_elapsed:.2f}s")
                        self.logger.log(f"TOOL_RESULT name={name} elapsed={tool_elapsed:.2f}s")
                    except Exception as exc:
                        result = {"error": str(exc)}
                        self.logger.error_reason(f"工具 {name} 执行异常：{exc}")
                        self.logger.log(f"TOOL_ERROR name={name} error={exc}")
                result_text = _safe_json(result)
                if str(name) == "start_job" and isinstance(result, dict) and result.get("job_id"):
                    attach_job_provenance(
                        [self.run_dir, *self.prior_run_dirs],
                        job_id=str(result["job_id"]),
                        tool_call_id=call_id,
                        event_id=str((tool_call_event or {}).get("eventId") or ""),
                    )
                    result["created_by_tool_call_id"] = call_id
                    result["created_by_event_id"] = str((tool_call_event or {}).get("eventId") or "")
                message_result_text = tool_message_text(str(name), result)
                self.memory_harness.observe_tool_result(str(name), args, result)
                self._emit_memory_state(event_sink)
                self.logger.preview(f"TOOL_OUTPUT {name}", result_text, max_chars=self.config.max_tool_result_chars)
                messages.append(ToolMessage(content=message_result_text[: self.config.max_tool_result_chars], tool_call_id=call_id))
                self._emit_event(
                    event_sink,
                    "tool_result",
                    turn=turn,
                    tool_name=str(name),
                    call_id=call_id,
                    toolCallId=call_id,
                    parentEventId=str((tool_call_event or {}).get("eventId") or ""),
                    ok=_tool_succeeded(result),
                    status="cancelled" if _pause_is_requested(pause_requested) else None,
                    result=result,
                )
                if _pause_is_requested(pause_requested):
                    status = "paused"
                    pending_question = "Task paused. Add instructions to continue."
                    pending_state_path = self._save_pause_checkpoint(messages, pending_question, event_sink=event_sink)
                    final_answer = pending_question
                    stop_after_tools = True
                    break
                if _tool_succeeded(result):
                    if str(name) == "inspect_workflow_skill":
                        workflow_skill_inspected = True
                    elif str(name) in {"read_skill_script", "inspect_skill_script_symbols", "inspect_skill_function"}:
                        workflow_skill_detail_seen = True
                if str(name) == "poll_job" and isinstance(result, dict) and result.get("status") == "completed":
                    completed_job_needs_finish = True
                elif str(name) == "start_job" and isinstance(result, dict) and result.get("status") == "running":
                    persistent_job_started = True
                elif str(name) == "finish_task" and _tool_succeeded(result):
                    completed_job_needs_finish = False
                if _is_user_input_request(result):
                    pending_question = str(result.get("question") or "")
                    pending_interaction = dict(result)
                    metadata = {
                        "tool_name": name,
                        "reason": result.get("reason", ""),
                        "required_fields": result.get("required_fields", []),
                        "resume_hint": result.get("resume_hint", ""),
                        "pending_interaction": pending_interaction,
                        "memory_state": self.memory_harness.snapshot(),
                    }
                    pending_path = save_pending_state(
                        config=self.config,
                        run_id=self.logger.run_id,
                        run_dir=self.run_dir,
                        log_path=self.logger.path,
                        messages=messages,
                        question=pending_question,
                        metadata=metadata,
                    )
                    pending_state_path = str(pending_path)
                    final_answer = pending_question
                    status = "needs_user_input"
                    self.logger.progress("主 Agent 暂停等待用户输入", f"resume_id={self.logger.run_id} state={pending_path}")
                    stop_after_tools = True
                    break
                tool_final = self._final_answer_from_successful_tool(str(name), result)
                if tool_final:
                    self.logger.progress("工具已形成最终答案", f"name={name}")
                    final_answer = tool_final
                    stop_after_tools = True
                    break
                if _tool_succeeded(result) and str(name) in {"execute_python", "execute_r"} and decision_turns + 1 >= max_turns:
                    self.logger.progress("执行已成功且接近轮数上限", f"name={name} turn={turn} max_turns={max_turns}")
                    final_answer = self._final_answer_from_successful_execution(str(name), result)
                    stop_after_tools = True
                    break
            if consumes_decision:
                decision_turns += 1
            if stop_after_tools:
                break

        self.logger.node("主 Agent 循环结束")
        self.logger.preview("FINAL_ANSWER", final_answer, max_chars=8000)
        self.memory_harness.finish_task(source_run_id=self.logger.run_id, status=status)
        self._emit_memory_state(event_sink)
        result = AgentRunResult(
            final_answer=final_answer,
            run_dir=str(self.run_dir),
            log_path=str(self.logger.path),
            turns=max(1, min(max_turns, decision_turns)),
            status=status,
            pending_question=pending_question,
            resume_id=self.logger.run_id if status in {"needs_user_input", "paused"} else "",
            pending_state_path=pending_state_path,
            pending_interaction=pending_interaction,
        )
        self._emit_event(event_sink, "final", turn=result.turns, content=result.final_answer, status=result.status, result=result.to_dict())
        self._emit_event(event_sink, "run_end", turn=result.turns, status=result.status, result=result.to_dict())
        return result

    def _save_pause_checkpoint(
        self,
        messages: list[Any],
        question: str,
        *,
        event_sink: Callable[[dict[str, Any]], None] | None = None,
    ) -> str:
        self.memory_harness.mark_paused()
        self._emit_memory_state(event_sink)
        path = save_pending_state(
            config=self.config,
            run_id=self.logger.run_id,
            run_dir=self.run_dir,
            log_path=self.logger.path,
            messages=messages,
            question=question,
            metadata={
                "resume_kind": "user_pause",
                "memory_state": self.memory_harness.snapshot(),
            },
        )
        self.logger.progress("主 Agent 已暂停", f"resume_id={self.logger.run_id} state={path}")
        return str(path)

    def _cancel_unstarted_tool_calls(
        self,
        messages: list[Any],
        tool_calls: list[Any],
        *,
        turn: int,
        event_sink: Callable[[dict[str, Any]], None] | None,
        reason: str = "Cancelled before execution because the user paused the task.",
    ) -> None:
        for call in tool_calls:
            name = str(call.get("name") or "tool")
            call_id = str(call.get("id") or f"paused-{turn}-{name}")
            args = call.get("args") or {}
            result = {"ok": False, "status": "cancelled", "error_reason": reason}
            self._emit_event(event_sink, "tool_call", turn=turn, tool_name=name, call_id=call_id, args=args)
            messages.append(ToolMessage(content=_safe_json(result, pretty=False), tool_call_id=call_id))
            self._emit_event(
                event_sink,
                "tool_result",
                turn=turn,
                tool_name=name,
                call_id=call_id,
                ok=False,
                status="cancelled",
                result=result,
            )

    def _append_pending_messages(
        self,
        messages: list[Any],
        callback: Callable[[], list[Any]] | None,
        *,
        event_sink: Callable[[dict[str, Any]], None] | None,
        turn: int | None = None,
    ) -> list[Any]:
        if callback is None:
            return []
        try:
            pending = [item for item in callback() if item]
        except Exception:
            return []
        for item in pending:
            if isinstance(item, dict) and item.get("type") == "job_callback":
                content = str(item.get("content") or f"Job {item.get('jobId')} status={item.get('status')}")
                messages.append(SystemMessage(content=f"Verified asynchronous job callback:\n{content}"))
                self._emit_event(
                    event_sink,
                    "job_callback_consumed",
                    callback_id=item.get("callbackId"),
                    job_id=item.get("jobId"),
                    status=item.get("status"),
                )
                continue
            if isinstance(item, dict) and item.get("type") == "user_message":
                content = str(item.get("content") or "").strip()
                if not content:
                    continue
                messages.append(HumanMessage(content=f"Session steering message:\n{content}"))
                self._emit_event(
                    event_sink,
                    "steering_message",
                    turn=turn,
                    messageId=item.get("messageId"),
                    content=content,
                    status="consumed",
                )
                continue
            content = str(item).strip()
            if not content:
                continue
            messages.append(HumanMessage(content=f"Session steering message:\n{content}"))
            self._emit_event(event_sink, "steering_message", content=content, status="consumed")
        return pending

    def _final_answer_from_successful_tool(self, tool_name: str, result: Any) -> str:
        if tool_name == "finish_task" and isinstance(result, dict) and result.get("ok"):
            return str(result.get("final_answer") or "")
        if tool_name == "run_code_workflow":
            if isinstance(result, dict) and not result.get("ok") and result.get("attempts"):
                return self._final_answer_from_failed_generated_workflow(result, title="通用代码工作流")
            return self._final_answer_from_generated_workflow(result, title="通用代码工作流")
        if tool_name == "run_skill_workflow":
            if isinstance(result, dict) and not result.get("ok") and result.get("attempts"):
                return self._final_answer_from_failed_generated_workflow(result, title="Skill-driven 生信工作流")
            return self._final_answer_from_generated_workflow(result, title="Skill-driven 生信工作流")
        return ""

    def _final_answer_from_successful_execution(self, tool_name: str, result: Any) -> str:
        if not isinstance(result, dict) or not result.get("ok"):
            return ""
        stdout = str(result.get("stdout") or "").strip()
        stderr = str(result.get("stderr") or "").strip()
        output_dir = Path(str(result.get("run_dir") or self.run_dir)) / "outputs"
        output_files: list[str] = []
        if output_dir.exists():
            output_files = sorted(path.name for path in output_dir.iterdir() if path.is_file())
        lines = [
            "代码执行已成功；由于已接近最大轮数，harness 在成功执行后直接收束，避免最后一轮再开启无法修复的新步骤。",
            "",
            f"- Tool：`{tool_name}`",
            f"- Docker profile：`{result.get('env_profile')}`",
            f"- 运行目录：`{result.get('run_dir')}`",
            f"- 输出目录：`{output_dir}`",
            f"- 脚本：`{result.get('script_path')}`",
            f"- 退出码：`{result.get('exit_code')}`",
        ]
        if output_files:
            lines.append(f"- 输出文件：{', '.join(f'`{name}`' for name in output_files)}")
        if stdout:
            lines.extend(["", "stdout 摘要：", stdout[-1600:]])
        if stderr:
            lines.extend(["", "stderr 摘要：", stderr[-800:]])
        return "\n".join(lines)

    def _final_answer_at_turn_limit(self, *, max_turns: int, requested_tool_calls: list[Any]) -> str:
        requested = ", ".join(str(call.get("name")) for call in requested_tool_calls)
        return (
            f"Reached max_turns={max_turns} before a final answer.\n\n"
            f"最后一轮模型仍请求新工具调用（{requested}），harness 已跳过执行，因为没有剩余轮次解释或修复新结果。"
        )

    def _final_answer_from_generated_workflow(self, result: Any, *, title: str) -> str:
        if not isinstance(result, dict) or not result.get("ok"):
            return ""
        execution = result.get("execution") if isinstance(result.get("execution"), dict) else {}
        attempts = result.get("attempts") if isinstance(result.get("attempts"), list) else []
        stdout = str(execution.get("stdout") or "")
        result_line = _extract_result_line(stdout)
        lines = [
            f"{title}已完成。",
            "",
        ]
        if result.get("skill_id"):
            lines.append(f"- 使用 Skill：`{result.get('skill_id')}`")
        lines.extend([
            f"- Runtime：`{result.get('runtime')}`",
            f"- Docker profile：`{result.get('env_profile')}`",
            f"- 输入数据：`{result.get('host_data_path') or '(未指定主输入文件)'}`",
            f"- 输出目录：`{result.get('output_dir')}`",
            f"- 运行目录：`{self.run_dir}`",
            f"- 日志文件：`{self.logger.path}`",
            f"- 生成/执行尝试次数：`{len(attempts)}`",
            f"- 最终脚本：`{execution.get('script_path')}`",
        ])
        if result_line:
            lines.extend(["", "工具返回摘要：", result_line])
        elif stdout.strip():
            lines.extend(["", "stdout 摘要：", stdout.strip()[-1200:]])
        lines.append("")
        if result.get("skill_id"):
            lines.append("这次代码由 Agent 根据 Skill 文档和脚本清单生成，并通过 Docker 执行反馈完成修复闭环。")
        else:
            lines.append("这次代码由 Agent 按任务直接生成，并通过 Docker 执行反馈完成修复闭环。")
        return "\n".join(lines)

    def _final_answer_from_failed_generated_workflow(self, result: dict[str, Any], *, title: str) -> str:
        attempts = result.get("attempts") if isinstance(result.get("attempts"), list) else []
        last = attempts[-1] if attempts else {}
        lines = [
            f"{title}未能在限定尝试次数内完成。",
            "",
        ]
        if result.get("skill_id"):
            lines.append(f"- 使用 Skill：`{result.get('skill_id')}`")
        lines.extend([
            f"- Runtime：`{result.get('runtime')}`",
            f"- Docker profile：`{result.get('env_profile')}`",
            f"- 输入数据：`{result.get('host_data_path') or '(未指定主输入文件)'}`",
            f"- 运行目录：`{self.run_dir}`",
            f"- 日志文件：`{self.logger.path}`",
            f"- 尝试次数：`{len(attempts)}`",
            f"- 最后脚本：`{last.get('script_path') or result.get('last_execution', {}).get('script_path')}`",
            f"- 最后退出码：`{last.get('exit_code')}`",
        ])
        stderr_tail = str(last.get("stderr_tail") or "").strip()
        stdout_tail = str(last.get("stdout_tail") or "").strip()
        if stderr_tail:
            lines.extend(["", "最后 stderr 摘要：", stderr_tail[-1200:]])
        elif stdout_tail:
            lines.extend(["", "最后 stdout 摘要：", stdout_tail[-1200:]])
        lines.append("")
        lines.append("日志中保留了每次代码生成、Docker 执行、失败原因和修复尝试，可据此继续缩小任务范围或增强对应 Skill。")
        return "\n".join(lines)


def _safe_json(value: Any, *, pretty: bool = True) -> str:
    try:
        if pretty:
            return json.dumps(value, ensure_ascii=False, default=str, indent=2)
        return json.dumps(value, ensure_ascii=False, default=str, separators=(",", ":"))
    except Exception:
        return str(value)


def _extract_result_line(stdout: Any) -> str:
    text = str(stdout or "")
    marker = "===RESULT==="
    if marker not in text:
        return ""
    line = text.split(marker, 1)[1].strip().splitlines()[0].strip()
    return line[:-3].strip() if line.endswith("===") else line


def _is_user_input_request(result: Any) -> bool:
    return isinstance(result, dict) and result.get("status") == "needs_user_input" and bool(result.get("question"))


def _tool_succeeded(result: Any) -> bool:
    if not isinstance(result, dict):
        return True
    if result.get("ok") is False:
        return False
    return not bool(result.get("error"))


def _tool_calls_consume_decision(tool_calls: list[Any]) -> bool:
    return any(str(call.get("name") or "") not in BUDGET_EXEMPT_TOOLS for call in tool_calls)


def _pause_is_requested(callback: Callable[[], bool] | None) -> bool:
    if callback is None:
        return False
    try:
        return bool(callback())
    except Exception:
        return False


def _should_block_execution_for_skill_detail(
    tool_name: str,
    args: dict[str, Any],
    skill_inspected: bool,
    skill_detail_seen: bool,
) -> bool:
    if tool_name not in {"start_job", "execute_python", "execute_r"} or not skill_inspected or skill_detail_seen:
        return False
    if tool_name == "execute_python" and _is_lightweight_h5ad_preview(str(args.get("code") or "")):
        return False
    return True


def _is_lightweight_h5ad_preview(code: str) -> bool:
    lowered = code.lower()
    if "read_h5ad" not in lowered or "print" not in lowered:
        return False
    blocked_markers = [
        "/work/outputs",
        ".write",
        "savefig",
        "sc.pp.",
        "sc.tl.",
        "rank_genes_groups",
        "normalize_total",
        "highly_variable_genes",
    ]
    return not any(marker in lowered for marker in blocked_markers)


def _looks_like_continuation_without_tool(text: str) -> bool:
    lowered = text.lower()
    continuation_markers = [
        "let me ",
        "i will ",
        "i'll ",
        "now i will",
        "让我",
        "我将",
        "我会",
        "现在我",
    ]
    completion_markers = [
        "completed",
        "已完成",
        "输出",
        "output",
        "结果",
        "summary",
    ]
    return any(marker in lowered for marker in continuation_markers) and not any(
        marker in lowered for marker in completion_markers
    )


def _is_unusable_final_answer(text: str) -> bool:
    normalized = text.strip()
    return not normalized or normalized in {".", "..", "...", "…"}


def _invalid_tool_call_detail(call: Any) -> str:
    name = str(call.get("name") or "unknown") if isinstance(call, dict) else "unknown"
    args = str(call.get("args") or "") if isinstance(call, dict) else ""
    error = str(call.get("error") or "") if isinstance(call, dict) else str(call)
    return f"tool={name}\narguments={args[:1600]}\nerror={error[:1600]}"


def _workflow_skill_detail_required(tool_name: str) -> dict[str, Any]:
    return {
        "ok": False,
        "status": "workflow_skill_detail_required",
        "error": (
            f"Before calling {tool_name}, inspect at least one relevant workflow skill script or function. "
            "Use inspect_skill_script_symbols for signatures, inspect_skill_function for one function, "
            "or read_skill_script for a script/reference/SKILL.md excerpt."
        ),
        "required_tools": ["inspect_skill_script_symbols", "inspect_skill_function", "read_skill_script"],
    }


def _context_lookup_call_key(tool_name: str, args: dict[str, Any]) -> tuple[str, str] | None:
    if tool_name not in CONTEXT_LOOKUP_TOOLS:
        return None
    return tool_name, json.dumps(args, ensure_ascii=False, default=str, sort_keys=True, separators=(",", ":"))


def _duplicate_context_lookup(tool_name: str, args: dict[str, Any]) -> dict[str, Any]:
    return {
        "ok": False,
        "status": "duplicate_context_lookup",
        "error": (
            f"Skipped consecutive duplicate context lookup {tool_name} with identical arguments. "
            "Use the previous result already in context, change the arguments to inspect a different scope, "
            "or take the next execution/editing action."
        ),
        "tool_name": tool_name,
        "args": args,
    }


def _persistent_job_state(run_dirs: list[Path]) -> tuple[bool, bool]:
    job_files = [path for run_dir in run_dirs for path in (run_dir / "jobs").glob("*/job.json")]
    if not job_files:
        return False, False
    latest = max(job_files, key=lambda path: path.stat().st_mtime)
    try:
        metadata = json.loads(latest.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        metadata = {}
    run_dir = latest.parent.parent.parent
    verified = (run_dir / "state" / "final_verification.json").is_file()
    return True, metadata.get("status") == "completed" and not verified
