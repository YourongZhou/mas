from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

from langchain_core.messages import HumanMessage, SystemMessage, ToolMessage

from .config import AgentConfig
from .llm import build_llm, runtime_summary
from .logging_utils import RunLogger
from .memory import build_memory_harness
from .run_state import save_pending_state
from .tools import build_tools


CONTEXT_LOOKUP_TOOLS = {
    "list_files",
    "read_file",
    "glob_search",
    "grep_text",
    "read_skill_script",
    "inspect_skill_script_symbols",
    "inspect_skill_function",
}
MAX_CONTEXT_LOOKUPS_AFTER_SKILL = 8


SYSTEM_PROMPT = """You are a single main-agent loop for bioinformatics analysis.

You own the whole task: understand the request, inspect relevant files or workflow skills, choose a Docker profile, generate code, execute it with tools, repair failures, and finish with a concise report.

Important operating rules:
- Use tools when you need filesystem context, workflow metadata, Docker profile details, or code execution.
- Choose the lightest sufficient route, similar to a CLI coding agent:
  1. Answer directly for conceptual questions that need no files or execution.
  2. Use file/search tools for inspection-only tasks.
  3. When a mas_2 workflow skill matches, inspect the skill, inspect relevant script symbols/functions, write the code yourself, then execute it with execute_python or execute_r.
  4. For ad hoc analysis, file conversion, summary statistics, plotting, or tasks with no matching skill, inspect files as needed and write the code yourself before using execute_python or execute_r.
- Prefer existing mas_2 workflow skills and Docker profiles when they match the task, but do not force a skill when the task is simpler without one.
- Gene lookup, enrichment, and cell-type annotation are not hardcoded tools. Use a matching workflow skill when one exists, such as functional-enrichment-from-degs or scrnaseq-scanpy-core-analysis; otherwise write a transparent Python/R script with the execution tools.
- Do not run marker-gene, differential-expression, trajectory, interaction, or report-generation steps unless the user explicitly asks for them. For a core single-cell demo, stop after QC, normalization, HVG, PCA, neighbors, Leiden, UMAP, summary outputs, and processed h5ad.
- Do not call sc.tl.rank_genes_groups or sc.pl.rank_genes_groups unless marker-gene or differential-expression analysis was explicitly requested.
- Manage context progressively: list skills first, inspect only the matching skill, then use read_skill_script, inspect_skill_script_symbols, or inspect_skill_function for the exact scripts/functions you need.
- Do not invent skill script paths or helper function names. Use only exact script paths returned by inspect_workflow_skill or list_files, and verify unknown function signatures before calling them.
- For a primary data path ending in .h5ad, use scanpy's sc.read_h5ad directly unless a verified skill function is specifically needed.
- Do not inspect every script in a workflow skill. Inspect only the few scripts/functions needed for the task, usually 3-6, then write and execute code.
- Hard guardrail: after inspect_workflow_skill for a workflow skill, do not call execute_python or execute_r until you have inspected at least one relevant skill script or function with read_skill_script, inspect_skill_script_symbols, or inspect_skill_function.
- For Python single-cell analysis prefer py-scverse-v1. For R/Seurat analysis prefer r-singlecell-v1. For general Python analysis prefer py-general-v1. For Bioconductor/DESeq2/survival prefer r-bioc-v1.
- Do not claim success until execution tools return ok=true or you clearly explain why execution could not be completed.
- Generated scripts must write outputs under /work/outputs when running in Docker.
- The repo is mounted read-only at /repo and the run workspace is mounted at /work during Docker execution.
- If an execution fails, repair based on the exact stdout/stderr and retry with a minimal change. Before guessing a helper function argument, call inspect_skill_function or inspect_skill_script_symbols to verify the signature.
- Execution attempts are globally limited by the harness, so each execute_python or execute_r call should be intentional.
- End with final output paths, key findings, and any unresolved blockers.
"""


def messages_for_model(messages: list[Any], config: AgentConfig, *, keep_recent_tool_messages: int = 2) -> list[Any]:
    tool_indices = [idx for idx, message in enumerate(messages) if getattr(message, "type", "") == "tool"]
    keep_indices = set(tool_indices[-keep_recent_tool_messages:])
    compacted: list[Any] = []
    for idx, message in enumerate(messages):
        if getattr(message, "type", "") != "tool":
            compacted.append(message)
            continue
        if idx in keep_indices:
            compacted.append(_trim_tool_message(message, max_chars=min(2200, config.max_tool_result_chars)))
            continue
        compacted.append(_compact_tool_message(message, max_chars=min(350, config.max_tool_result_chars)))
    return compacted


def messages_with_turn_guidance(messages: list[Any], config: AgentConfig, *, turn: int, max_turns: int) -> list[Any]:
    remaining_after_this = max_turns - turn
    guidance = (
        f"Turn budget: this is turn {turn} of {max_turns}; "
        f"{remaining_after_this} turn(s) remain after this response. "
        "Use tools only when the next result can still be interpreted or repaired within the remaining budget. "
        "On the final turn, do not start a new analysis step; summarize the best verified results and any blockers."
    )
    return [*messages_for_model(messages, config), HumanMessage(content=guidance)]


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


def _compact_tool_message(message: Any, *, max_chars: int) -> ToolMessage:
    content = str(getattr(message, "content", ""))
    preview = content[:max_chars].replace("\n", "\\n")
    tool_call_id = str(getattr(message, "tool_call_id", ""))
    return ToolMessage(
        content=(
            "[tool output compressed to keep model context small; "
            "call a file/skill inspection tool again if exact details are needed]"
            f"\npreview: {preview}"
        ),
        tool_call_id=tool_call_id,
    )


def tool_message_text(tool_name: str, result: Any) -> str:
    if tool_name in {"execute_python", "execute_r"} and isinstance(result, dict) and not result.get("ok"):
        compact = {
            "ok": False,
            "status": result.get("status"),
            "exit_code": result.get("exit_code"),
            "error_reason": result.get("error_reason") or result.get("error") or "",
            "stdout_tail": str(result.get("stdout") or "")[-1200:],
            "stderr_tail": str(result.get("stderr") or "")[-1200:],
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
        }


class BioAgent:
    def __init__(self, *, config: AgentConfig, logger: RunLogger, run_dir: Path) -> None:
        self.config = config
        self.logger = logger
        self.run_dir = run_dir
        self.memory_harness = build_memory_harness(config)
        self.tools = build_tools(config, logger, run_dir, memory_harness=self.memory_harness)
        self.tool_map = {tool.name: tool for tool in self.tools}
        self.llm = build_llm(config).bind_tools(self.tools)

    def _initial_user_message(self, task: str, data_path: str | None, result_dir: str | None) -> str:
        memory_context = self.memory_harness.recall_text(task)
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
            parts.append(f"Requested result directory: {result_dir}")
        if memory_context:
            parts.extend(["", memory_context])
        return "\n".join(parts)

    def _emit_event(self, event_sink: Callable[[dict[str, Any]], None] | None, event_type: str, **payload: Any) -> None:
        if event_sink is None:
            return
        event = {
            "type": event_type,
            "run_id": self.logger.run_id,
            "run_dir": str(self.run_dir),
            "log_path": str(self.logger.path),
        }
        event.update(payload)
        event_sink(event)

    def run(
        self,
        task: str,
        *,
        data_path: str | None = None,
        result_dir: str | None = None,
        max_turns: int = 20,
        initial_messages: list[Any] | None = None,
        resume_answer: str | None = None,
        event_sink: Callable[[dict[str, Any]], None] | None = None,
    ) -> AgentRunResult:
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

        if initial_messages is None:
            messages: list[Any] = [
                SystemMessage(content=SYSTEM_PROMPT),
                HumanMessage(content=self._initial_user_message(task, data_path, result_dir)),
            ]
        else:
            messages = list(initial_messages)
            if resume_answer is not None:
                messages.append(HumanMessage(content=f"User answer for the pending question:\n{resume_answer}"))
        final_answer = ""
        status = "completed"
        pending_question = ""
        pending_state_path = ""
        workflow_skill_inspected = False
        workflow_skill_detail_seen = False
        context_lookups_after_skill = 0

        for turn in range(1, max_turns + 1):
            self.logger.node(f"主 Agent 第 {turn} 轮推理", "判断是否需要读取文件、检查 Skill、执行代码或直接回答。")
            self._emit_event(event_sink, "turn_start", turn=turn, max_turns=max_turns)
            started = time.perf_counter()
            self.logger.progress("调用主模型", f"turn={turn}")
            self._emit_event(event_sink, "model_start", turn=turn, max_turns=max_turns)
            response = self.llm.invoke(messages_with_turn_guidance(messages, self.config, turn=turn, max_turns=max_turns))
            elapsed = time.perf_counter() - started
            messages.append(response)
            self.logger.progress("主模型返回", f"turn={turn} elapsed={elapsed:.2f}s tool_calls={len(response.tool_calls or [])}")
            self.logger.log(f"LLM_RESPONSE turn={turn} elapsed={elapsed:.2f}s tool_calls={len(response.tool_calls or [])}")
            self._emit_event(
                event_sink,
                "model_end",
                turn=turn,
                elapsed_s=elapsed,
                tool_call_count=len(response.tool_calls or []),
            )
            if response.content:
                self.logger.preview("ASSISTANT_CONTENT", str(response.content), max_chars=6000)
                self._emit_event(event_sink, "assistant_message", turn=turn, content=str(response.content))

            tool_calls = response.tool_calls or []
            if not tool_calls:
                response_text = str(response.content or "").strip()
                if (not response_text or _looks_like_continuation_without_tool(response_text)) and turn < max_turns:
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
                self.logger.progress("主 Agent 决定结束", "本轮没有工具调用，使用模型文本作为最终回答。")
                final_answer = response_text
                break
            if turn == max_turns and not any(str(call.get("name")) == "request_user_input" for call in tool_calls):
                self.logger.progress("已到最后一轮，跳过新的工具调用", f"requested={','.join(str(call.get('name')) for call in tool_calls)}")
                response_text = str(response.content or "").strip()
                final_answer = response_text or self._final_answer_at_turn_limit(max_turns=max_turns, requested_tool_calls=tool_calls)
                break

            stop_after_tools = False
            for call in tool_calls:
                name = call.get("name")
                args = call.get("args") or {}
                call_id = call.get("id")
                self.logger.node(f"工具调用：{name}", f"call_id={call_id}")
                self.logger.log(f"TOOL_CALL name={name} id={call_id}")
                self.logger.preview("工具入参", _safe_json(args), max_chars=4000)
                self._emit_event(event_sink, "tool_call", turn=turn, tool_name=str(name), call_id=call_id, args=args)
                tool = self.tool_map.get(str(name))
                if tool is None:
                    result: Any = {"error": f"Unknown tool: {name}"}
                    self.logger.error_reason(f"Unknown tool: {name}")
                elif (
                    workflow_skill_inspected
                    and str(name) in CONTEXT_LOOKUP_TOOLS
                    and context_lookups_after_skill >= MAX_CONTEXT_LOOKUPS_AFTER_SKILL
                ):
                    result = _context_lookup_budget_exhausted(str(name), context_lookups_after_skill)
                    self.logger.progress("上下文查阅次数已达上限", f"name={name} count={context_lookups_after_skill}")
                    self.logger.log(f"TOOL_GUARD name={name} reason=context_lookup_budget_exhausted")
                elif _should_block_execution_for_skill_detail(str(name), args, workflow_skill_inspected, workflow_skill_detail_seen):
                    result = _workflow_skill_detail_required(str(name))
                    self.logger.progress("执行前需要读取 Skill 脚本", f"name={name}")
                    self.logger.log(f"TOOL_GUARD name={name} reason=workflow_skill_detail_required")
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
                message_result_text = tool_message_text(str(name), result)
                self.logger.preview(f"TOOL_OUTPUT {name}", result_text, max_chars=self.config.max_tool_result_chars)
                messages.append(ToolMessage(content=message_result_text[: self.config.max_tool_result_chars], tool_call_id=call_id))
                self._emit_event(
                    event_sink,
                    "tool_result",
                    turn=turn,
                    tool_name=str(name),
                    call_id=call_id,
                    ok=_tool_succeeded(result),
                    result=result,
                )
                if _tool_succeeded(result):
                    if str(name) == "inspect_workflow_skill":
                        workflow_skill_inspected = True
                    elif str(name) in {"read_skill_script", "inspect_skill_script_symbols", "inspect_skill_function"}:
                        workflow_skill_detail_seen = True
                    if workflow_skill_inspected and str(name) in CONTEXT_LOOKUP_TOOLS:
                        context_lookups_after_skill += 1
                if _is_user_input_request(result):
                    pending_question = str(result.get("question") or "")
                    metadata = {
                        "tool_name": name,
                        "reason": result.get("reason", ""),
                        "required_fields": result.get("required_fields", []),
                        "resume_hint": result.get("resume_hint", ""),
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
                if _tool_succeeded(result) and str(name) in {"execute_python", "execute_r"} and turn >= max_turns - 1:
                    self.logger.progress("执行已成功且接近轮数上限", f"name={name} turn={turn} max_turns={max_turns}")
                    final_answer = self._final_answer_from_successful_execution(str(name), result)
                    stop_after_tools = True
                    break
            if stop_after_tools:
                break
        else:
            final_answer = f"Reached max_turns={max_turns} before a final answer."

        self.logger.node("主 Agent 循环结束")
        self.logger.preview("FINAL_ANSWER", final_answer, max_chars=8000)
        result = AgentRunResult(
            final_answer=final_answer,
            run_dir=str(self.run_dir),
            log_path=str(self.logger.path),
            turns=min(max_turns, turn if "turn" in locals() else 0),
            status=status,
            pending_question=pending_question,
            resume_id=self.logger.run_id if status == "needs_user_input" else "",
            pending_state_path=pending_state_path,
        )
        self._emit_event(event_sink, "final", turn=result.turns, content=result.final_answer, status=result.status, result=result.to_dict())
        self._emit_event(event_sink, "run_end", turn=result.turns, status=result.status, result=result.to_dict())
        return result

    def _final_answer_from_successful_tool(self, tool_name: str, result: Any) -> str:
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


def _should_block_execution_for_skill_detail(
    tool_name: str,
    args: dict[str, Any],
    skill_inspected: bool,
    skill_detail_seen: bool,
) -> bool:
    if tool_name not in {"execute_python", "execute_r"} or not skill_inspected or skill_detail_seen:
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


def _context_lookup_budget_exhausted(tool_name: str, count: int) -> dict[str, Any]:
    return {
        "ok": False,
        "status": "context_lookup_budget_exhausted",
        "error": (
            f"Context lookup budget after workflow skill inspection is exhausted "
            f"({count}/{MAX_CONTEXT_LOOKUPS_AFTER_SKILL}); skipped {tool_name}. "
            "Stop inspecting files/scripts. Use the inspected skill context to call execute_python/execute_r now, "
            "or provide a final answer if execution is impossible."
        ),
        "allowed_next_tools": ["execute_python", "execute_r"],
    }
