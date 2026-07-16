from __future__ import annotations

import queue
import threading
from pathlib import Path
from typing import Iterator

from .agent import BioAgent
from .config import AgentConfig
from .logging_utils import RunLogger, now_stamp
from .run_state import clear_pending_state, load_pending_state


def run_bio_agent(
    task: str,
    *,
    data_path: str | None = None,
    result_dir: str | None = None,
    max_turns: int = 20,
) -> dict:
    """Notebook-facing entrypoint."""

    config = AgentConfig.from_env()
    run_id = f"run_{now_stamp()}"
    run_dir = config.runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    with RunLogger(config.logs_dir, run_id=run_id) as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        result = agent.run(
            task,
            data_path=data_path,
            result_dir=result_dir,
            max_turns=max_turns,
        )
        return result.to_dict()


def run_bio_agent_stream(
    task: str,
    *,
    data_path: str | None = None,
    result_dir: str | None = None,
    max_turns: int = 20,
) -> Iterator[dict]:
    """Notebook-facing streaming entrypoint."""

    config = AgentConfig.from_env()
    run_id = f"run_{now_stamp()}"
    run_dir = config.runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    events: queue.Queue[dict | None] = queue.Queue()
    terminal_seen = {"final": False, "run_end": False}

    def emit(event: dict) -> None:
        if event.get("type") == "final":
            terminal_seen["final"] = True
        elif event.get("type") == "run_end":
            terminal_seen["run_end"] = True
        events.put(event)

    def worker() -> None:
        try:
            with RunLogger(config.logs_dir, run_id=run_id) as logger:
                agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
                result = agent.run(
                    task,
                    data_path=data_path,
                    result_dir=result_dir,
                    max_turns=max_turns,
                    event_sink=emit,
                )
                result_dict = result.to_dict()
                if not terminal_seen["final"]:
                    emit(
                        {
                            "type": "final",
                            "run_id": run_id,
                            "run_dir": str(run_dir),
                            "log_path": str(logger.path),
                            "turn": result.turns,
                            "content": result.final_answer,
                            "status": result.status,
                            "result": result_dict,
                        }
                    )
                if not terminal_seen["run_end"]:
                    emit(
                        {
                            "type": "run_end",
                            "run_id": run_id,
                            "run_dir": str(run_dir),
                            "log_path": str(logger.path),
                            "turn": result.turns,
                            "status": result.status,
                            "result": result_dict,
                        }
                    )
        except Exception as exc:
            emit(
                {
                    "type": "error",
                    "run_id": run_id,
                    "run_dir": str(run_dir),
                    "log_path": str(config.logs_dir / f"{run_id}.log"),
                    "error": str(exc),
                }
            )
        finally:
            events.put(None)

    thread = threading.Thread(target=worker, name=f"bioagent-stream-{run_id}", daemon=True)
    thread.start()
    while True:
        event = events.get()
        if event is None:
            thread.join()
            break
        yield event


def format_bio_agent_event(event: dict) -> str:
    """Compact text formatting for notebook streaming display."""

    event_type = event.get("type", "")
    if event_type == "run_start":
        return f"[run_start] {event.get('run_id')} -> {event.get('run_dir')}"
    if event_type == "turn_start":
        return f"[turn] {event.get('turn')}/{event.get('max_turns')}"
    if event_type == "model_start":
        return f"[model] turn={event.get('turn')} calling model"
    if event_type == "model_end":
        return f"[model] turn={event.get('turn')} returned tool_calls={event.get('tool_call_count')}"
    if event_type == "assistant_message":
        return f"[assistant]\n{event.get('content', '')}"
    if event_type == "tool_call":
        return f"[tool_call] {event.get('tool_name')} id={event.get('call_id')}"
    if event_type == "tool_result":
        result = event.get("result") if isinstance(event.get("result"), dict) else {}
        status = "ok" if event.get("ok") else "error"
        reason = result.get("error_reason") or result.get("error") or ""
        suffix = f" reason={reason}" if reason else ""
        return f"[tool_result] {event.get('tool_name')} {status}{suffix}"
    if event_type == "final":
        return f"[final]\n{event.get('content', '')}"
    if event_type == "run_end":
        return f"[run_end] status={event.get('status')} run_id={event.get('run_id')}"
    if event_type == "error":
        return f"[error] {event.get('error')}"
    return str(event)


def display_bio_agent_event(event: dict) -> dict:
    """Display one stream event in a notebook, with a plain print fallback."""

    text = format_bio_agent_event(event)
    try:
        from IPython.display import Markdown, display

        display(Markdown(f"```text\n{text}\n```"))
    except Exception:
        print(text)
    return event


def resume_bio_agent(
    resume_id: str,
    user_answer: str,
    *,
    max_turns: int = 20,
) -> dict:
    """Resume a run that previously returned status='needs_user_input'."""

    config = AgentConfig.from_env()
    pending = load_pending_state(config, resume_id)
    pending.run_dir.mkdir(parents=True, exist_ok=True)

    with RunLogger(config.logs_dir, run_id=pending.run_id, append=True) as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=pending.run_dir)
        result = agent.run(
            task="",
            max_turns=max_turns,
            initial_messages=pending.messages,
            resume_answer=user_answer,
        )
        if result.status != "needs_user_input":
            clear_pending_state(config, resume_id)
        return result.to_dict()


def run_bmmc_singlecell_demo(max_turns: int = 20) -> dict:
    """Run the notebook demo task through the single-agent loop."""

    return run_bio_agent(
        task=(
            "我想用 data 目录下的 bmmc_b_cell.h5ad 单细胞数据进行分析。"
            "请使用新的单主 Agent + 工具调用模式：先确认可用的 Scanpy skill 和 Docker profile，"
            "然后基于 scrnaseq-scanpy-core-analysis Skill 自主阅读 SKILL.md、检查相关 scripts 的符号和函数签名"
            "（优先使用 inspect_skill_script_symbols、inspect_skill_function、read_skill_script），"
            "自己生成分析脚本，并用 execute_python 执行；如果失败，基于 stdout/stderr 继续查看相关函数签名或源码后最小修复，"
            "完成稳健的核心流程：读取数据、QC、标准化、HVG、PCA、邻居图、Leiden、UMAP、"
            "summary JSON、关键图和 processed h5ad。不要额外展开非必要的 DE、marker gene、rank_genes_groups、轨迹、互作或复杂 PDF 报告。"
            "这个 h5ad 可能是 dense numpy AnnData；QC 请手动从 adata.X 兼容 dense/sparse 计算 n_counts、n_genes、pct_mito，"
            "不要依赖 calculate_qc_metrics 生成的列名。保存图片优先用 matplotlib/seaborn 的 fig.savefig 到 /work/outputs，"
            "不要把绝对路径传给 scanpy plotting 的 save 参数。"
        ),
        data_path=r"mas_2\data\bmmc_b_cell.h5ad",
        max_turns=max_turns,
    )
