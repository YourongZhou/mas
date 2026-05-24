from __future__ import annotations

from pathlib import Path

from .agent import BioAgent
from .config import AgentConfig
from .logging_utils import RunLogger, now_stamp


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


def run_bmmc_singlecell_demo(max_turns: int = 12) -> dict:
    """Run the notebook demo task through the single-agent loop."""

    return run_bio_agent(
        task=(
            "我想用 data 目录下的 bmmc_b_cell.h5ad 单细胞数据进行分析。"
            "请使用新的单主 Agent + 工具调用模式：先确认可用的 Scanpy skill 和 Docker profile，"
            "然后调用 run_scanpy_singlecell_pipeline 完成标准单细胞核心分析，最后汇总输出文件。"
        ),
        data_path=r"mas_2\data\bmmc_b_cell.h5ad",
        max_turns=max_turns,
    )

