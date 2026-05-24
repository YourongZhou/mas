from __future__ import annotations

from pathlib import Path

from bioagent.config import AgentConfig
from bioagent.docker_runner import DockerRunner
from bioagent.logging_utils import RunLogger


def execute_python_impl(
    config: AgentConfig,
    logger: RunLogger,
    run_dir: Path,
    code: str,
    env_profile: str = "py-general-v1",
    requirements: str = "",
    timeout_s: int = 900,
) -> dict:
    runner = DockerRunner(config=config, logger=logger, run_dir=run_dir)
    return runner.execute_python(
        code,
        env_profile=env_profile,
        requirements=requirements,
        timeout_s=timeout_s,
    ).to_dict()


def execute_r_impl(
    config: AgentConfig,
    logger: RunLogger,
    run_dir: Path,
    code: str,
    env_profile: str = "r-bioc-v1",
    timeout_s: int = 900,
) -> dict:
    runner = DockerRunner(config=config, logger=logger, run_dir=run_dir)
    return runner.execute_r(code, env_profile=env_profile, timeout_s=timeout_s).to_dict()

