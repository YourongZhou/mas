from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from langchain_core.messages import HumanMessage, SystemMessage, ToolMessage

from .config import AgentConfig
from .llm import build_llm, runtime_summary
from .logging_utils import RunLogger
from .tools import build_tools


SYSTEM_PROMPT = """You are a single main-agent loop for bioinformatics analysis.

You own the whole task: understand the request, inspect relevant files or workflow skills, choose a Docker profile, generate code, execute it with tools, repair failures, and finish with a concise report.

Important operating rules:
- Use tools when you need filesystem context, workflow metadata, Docker profile details, or code execution.
- Choose the lightest sufficient route, similar to a CLI coding agent:
  1. Answer directly for conceptual questions that need no files or execution.
  2. Use file/search tools for inspection-only tasks.
  3. Use run_skill_workflow only when an existing mas_2 workflow skill clearly matches the task.
  4. Use run_code_workflow for ad hoc analysis, file conversion, summary statistics, plotting, or tasks with no matching skill.
  5. Use execute_python/execute_r for very small one-shot scripts when no repair loop is needed.
- Prefer existing mas_2 workflow skills and Docker profiles when they match the task, but do not force a skill when the task is simpler without one.
- Gene lookup, enrichment, and cell-type annotation are not hardcoded tools. Route them through a matching workflow skill when one exists, such as functional-enrichment-from-degs or scrnaseq-scanpy-core-analysis; otherwise use run_code_workflow so code, dependencies, outputs, and repair attempts are visible in logs.
- Manage context progressively: inspect/list first, load only the matching skill, and rely on run_skill_workflow's compact script/reference selection rather than pasting unrelated workflow material into the main conversation.
- For Python single-cell analysis prefer py-scverse-v1. For R/Seurat analysis prefer r-singlecell-v1. For general Python analysis prefer py-general-v1. For Bioconductor/DESeq2/survival prefer r-bioc-v1.
- Do not claim success until execution tools return ok=true or you clearly explain why execution could not be completed.
- Generated scripts must write outputs under /work/outputs when running in Docker.
- The repo is mounted read-only at /repo and the run workspace is mounted at /work during Docker execution.
- If an execution fails, repair based on the exact stdout/stderr and retry with a minimal change.
- End with final output paths, key findings, and any unresolved blockers.
"""


@dataclass
class AgentRunResult:
    final_answer: str
    run_dir: str
    log_path: str
    turns: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "final_answer": self.final_answer,
            "run_dir": self.run_dir,
            "log_path": self.log_path,
            "turns": self.turns,
        }


class BioAgent:
    def __init__(self, *, config: AgentConfig, logger: RunLogger, run_dir: Path) -> None:
        self.config = config
        self.logger = logger
        self.run_dir = run_dir
        self.tools = build_tools(config, logger, run_dir)
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
        ]
        if data_path:
            parts.append(f"Primary data path: {data_path}")
            parts.append("Inside Docker, host repo files are visible under /repo. If the data path is inside the repo, translate it to /repo-relative path.")
        if result_dir:
            parts.append(f"Requested result directory: {result_dir}")
        return "\n".join(parts)

    def run(self, task: str, *, data_path: str | None = None, result_dir: str | None = None, max_turns: int = 20) -> AgentRunResult:
        self.logger.node("主 Agent 循环启动", "加载运行配置、工具清单，并开始理解用户任务。")
        self.logger.preview("Runtime", runtime_summary(self.config))
        self.logger.progress("可用工具", ", ".join(self.tool_map))

        messages: list[Any] = [
            SystemMessage(content=SYSTEM_PROMPT),
            HumanMessage(content=self._initial_user_message(task, data_path, result_dir)),
        ]
        final_answer = ""

        for turn in range(1, max_turns + 1):
            self.logger.node(f"主 Agent 第 {turn} 轮推理", "判断是否需要读取文件、检查 Skill、执行代码或直接回答。")
            started = time.perf_counter()
            self.logger.progress("调用主模型", f"turn={turn}")
            response = self.llm.invoke(messages)
            elapsed = time.perf_counter() - started
            messages.append(response)
            self.logger.progress("主模型返回", f"turn={turn} elapsed={elapsed:.2f}s tool_calls={len(response.tool_calls or [])}")
            self.logger.log(f"LLM_RESPONSE turn={turn} elapsed={elapsed:.2f}s tool_calls={len(response.tool_calls or [])}")
            if response.content:
                self.logger.preview("ASSISTANT_CONTENT", str(response.content), max_chars=6000)

            tool_calls = response.tool_calls or []
            if not tool_calls:
                self.logger.progress("主 Agent 决定结束", "本轮没有工具调用，使用模型文本作为最终回答。")
                final_answer = str(response.content or "").strip()
                break

            stop_after_tools = False
            for call in tool_calls:
                name = call.get("name")
                args = call.get("args") or {}
                call_id = call.get("id")
                self.logger.node(f"工具调用：{name}", f"call_id={call_id}")
                self.logger.log(f"TOOL_CALL name={name} id={call_id}")
                self.logger.preview("工具入参", _safe_json(args), max_chars=4000)
                tool = self.tool_map.get(str(name))
                if tool is None:
                    result: Any = {"error": f"Unknown tool: {name}"}
                    self.logger.error_reason(f"Unknown tool: {name}")
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
                self.logger.preview(f"TOOL_OUTPUT {name}", result_text, max_chars=self.config.max_tool_result_chars)
                messages.append(ToolMessage(content=result_text[: self.config.max_tool_result_chars], tool_call_id=call_id))
                tool_final = self._final_answer_from_successful_tool(str(name), result)
                if tool_final:
                    self.logger.progress("工具已形成最终答案", f"name={name}")
                    final_answer = tool_final
                    stop_after_tools = True
                    break
            if stop_after_tools:
                break
        else:
            final_answer = f"Reached max_turns={max_turns} before a final answer."

        self.logger.node("主 Agent 循环结束")
        self.logger.preview("FINAL_ANSWER", final_answer, max_chars=8000)
        return AgentRunResult(
            final_answer=final_answer,
            run_dir=str(self.run_dir),
            log_path=str(self.logger.path),
            turns=min(max_turns, turn if "turn" in locals() else 0),
        )

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


def _safe_json(value: Any) -> str:
    try:
        return json.dumps(value, ensure_ascii=False, default=str, indent=2)
    except Exception:
        return str(value)


def _extract_result_line(stdout: Any) -> str:
    text = str(stdout or "")
    marker = "===RESULT==="
    if marker not in text:
        return ""
    line = text.split(marker, 1)[1].strip().splitlines()[0].strip()
    return line[:-3].strip() if line.endswith("===") else line

