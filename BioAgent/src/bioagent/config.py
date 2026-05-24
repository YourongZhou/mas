from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from dotenv import load_dotenv


def _find_repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _first_nonempty(*names: str) -> str:
    for name in names:
        value = os.getenv(name)
        if value and value.strip():
            return value.strip()
    return ""


@dataclass(frozen=True)
class AgentConfig:
    repo_root: Path
    project_root: Path
    mas2_root: Path
    workflows_root: Path
    docker_root: Path
    logs_dir: Path
    runs_dir: Path
    api_key: str
    base_url: str
    model_name: str
    temperature: float
    request_timeout: float
    mimo_thinking_type: str
    max_tool_result_chars: int = 12000

    @classmethod
    def from_env(cls) -> "AgentConfig":
        repo_root = _find_repo_root()
        project_root = repo_root / "BioAgent"
        mas2_root = repo_root / "mas_2"
        load_dotenv(mas2_root / ".env", override=True)

        logs_dir = project_root / "logs"
        runs_dir = project_root / "runs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        runs_dir.mkdir(parents=True, exist_ok=True)

        return cls(
            repo_root=repo_root,
            project_root=project_root,
            mas2_root=mas2_root,
            workflows_root=mas2_root / "workflows",
            docker_root=mas2_root / "docker",
            logs_dir=logs_dir,
            runs_dir=runs_dir,
            api_key=_first_nonempty("OPENAI_API_KEY", "DASHSCOPE_API_KEY"),
            base_url=os.getenv("BASE_URL", "https://dashscope.aliyuncs.com/compatible-mode/v1"),
            model_name=os.getenv("MODEL_NAME", "qwen-turbo"),
            temperature=float(os.getenv("TEMPERATURE", "0.2")),
            request_timeout=float(os.getenv("LLM_REQUEST_TIMEOUT", "600")),
            mimo_thinking_type=os.getenv("MIMO_THINKING_TYPE", "").strip().lower(),
        )

    def mask_api_key(self) -> str:
        value = self.api_key.strip()
        if not value:
            return "(empty)"
        if len(value) <= 10:
            return value[:2] + "***"
        return f"{value[:6]}...{value[-4:]}"

