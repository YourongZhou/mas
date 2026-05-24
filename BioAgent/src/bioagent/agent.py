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
- Prefer existing mas_2 workflow skills and Docker profiles when they match the task.
- For Python single-cell analysis prefer py-scverse-v1. For R/Seurat analysis prefer r-singlecell-v1. For general Python analysis prefer py-general-v1. For Bioconductor/DESeq2/survival prefer r-bioc-v1.
- When the user asks to analyze an .h5ad single-cell dataset with no custom constraints, first inspect the Scanpy skill/catalog if needed, then call run_scanpy_singlecell_pipeline. This is the preferred cannon-style path for the standard scRNA-seq core workflow.
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
        self.logger.section("AGENT LOOP START")
        self.logger.preview("Runtime", runtime_summary(self.config))
        self.logger.log("TOOLS " + ", ".join(self.tool_map))

        messages: list[Any] = [
            SystemMessage(content=SYSTEM_PROMPT),
            HumanMessage(content=self._initial_user_message(task, data_path, result_dir)),
        ]
        final_answer = ""

        for turn in range(1, max_turns + 1):
            self.logger.section(f"TURN {turn}")
            started = time.perf_counter()
            response = self.llm.invoke(messages)
            elapsed = time.perf_counter() - started
            messages.append(response)
            self.logger.log(f"LLM_RESPONSE turn={turn} elapsed={elapsed:.2f}s tool_calls={len(response.tool_calls or [])}")
            if response.content:
                self.logger.preview("ASSISTANT_CONTENT", str(response.content), max_chars=6000)

            tool_calls = response.tool_calls or []
            if not tool_calls:
                final_answer = str(response.content or "").strip()
                break

            stop_after_tools = False
            for call in tool_calls:
                name = call.get("name")
                args = call.get("args") or {}
                call_id = call.get("id")
                self.logger.log(f"TOOL_CALL name={name} id={call_id} args={_safe_json(args)}")
                tool = self.tool_map.get(str(name))
                if tool is None:
                    result: Any = {"error": f"Unknown tool: {name}"}
                else:
                    try:
                        tool_started = time.perf_counter()
                        result = tool.invoke(args)
                        self.logger.log(f"TOOL_RESULT name={name} elapsed={time.perf_counter() - tool_started:.2f}s")
                    except Exception as exc:
                        result = {"error": str(exc)}
                        self.logger.log(f"TOOL_ERROR name={name} error={exc}")
                result_text = _safe_json(result)
                self.logger.preview(f"TOOL_OUTPUT {name}", result_text, max_chars=self.config.max_tool_result_chars)
                messages.append(ToolMessage(content=result_text[: self.config.max_tool_result_chars], tool_call_id=call_id))
                tool_final = self._final_answer_from_successful_tool(str(name), result)
                if tool_final:
                    final_answer = tool_final
                    stop_after_tools = True
                    break
            if stop_after_tools:
                break
        else:
            final_answer = f"Reached max_turns={max_turns} before a final answer."

        self.logger.section("AGENT LOOP END")
        self.logger.preview("FINAL_ANSWER", final_answer, max_chars=8000)
        return AgentRunResult(
            final_answer=final_answer,
            run_dir=str(self.run_dir),
            log_path=str(self.logger.path),
            turns=min(max_turns, turn if "turn" in locals() else 0),
        )

    def _final_answer_from_successful_tool(self, tool_name: str, result: Any) -> str:
        if tool_name != "run_scanpy_singlecell_pipeline":
            return ""
        if not isinstance(result, dict) or not result.get("ok"):
            return ""

        summary = _extract_scanpy_summary(result.get("stdout", ""))
        output_dir = result.get("expected_output_dir") or str(self.run_dir / "outputs" / "scrnaseq_scanpy_core")
        host_data_path = result.get("host_data_path", "")
        skill_id = result.get("skill_id", "scrnaseq-scanpy-core-analysis")

        lines = [
            "单细胞标准分析流程已完成。",
            "",
            f"- 使用 Skill：`{skill_id}`",
            f"- 输入数据：`{host_data_path}`",
            f"- 输出目录：`{output_dir}`",
            f"- 运行目录：`{self.run_dir}`",
            f"- 日志文件：`{self.logger.path}`",
        ]
        if summary:
            lines.extend(
                [
                    "",
                    "关键结果：",
                    f"- 初始矩阵维度：`{summary.get('initial_shape')}`",
                    f"- 最终矩阵维度：`{summary.get('final_shape')}`",
                    f"- 标准化策略：`{summary.get('normalization')}`",
                    f"- 高变基因数：`{summary.get('n_hvg')}`",
                    f"- PCA 维度数：`{summary.get('n_pcs')}`",
                    f"- 默认聚类列：`{summary.get('cluster_key')}`",
                    f"- 多分辨率聚类：`{summary.get('cluster_status')}`",
                    f"- marker 结果行数：`{summary.get('marker_rows')}`",
                ]
            )
            outputs = summary.get("outputs")
            if isinstance(outputs, list) and outputs:
                lines.extend(["", "主要输出文件："])
                for name in outputs[:16]:
                    lines.append(f"- `{name}`")
        lines.append("")
        lines.append("这次是由单主 Agent 选择工具并调用 Docker 沙箱完成的，不再走旧版多 Agent 的 code_dev/critic 循环。")
        return "\n".join(lines)


def _safe_json(value: Any) -> str:
    try:
        return json.dumps(value, ensure_ascii=False, default=str, indent=2)
    except Exception:
        return str(value)


def _extract_scanpy_summary(stdout: Any) -> dict[str, Any]:
    text = str(stdout or "")
    marker = "===RESULT==="
    if marker not in text:
        return {}
    payload = text.split(marker, 1)[1]
    if marker in payload:
        payload = payload.split(marker, 1)[0]
    else:
        payload = payload.rsplit("===", 1)[0]
    try:
        parsed = json.loads(payload)
    except Exception:
        return {}
    return parsed if isinstance(parsed, dict) else {}

