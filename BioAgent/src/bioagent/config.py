from __future__ import annotations

import os
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from dotenv import load_dotenv


MODEL_SETTINGS_FILE = "model-settings.json"
MODEL_SETTING_FIELDS = {
    "provider",
    "api_key",
    "base_url",
    "model_name",
    "temperature",
    "request_timeout",
    "mimo_thinking_type",
    "chat_template_enable_thinking",
}


def model_settings_path(project_root: Path) -> Path:
    return project_root / "state" / MODEL_SETTINGS_FILE


def load_model_settings(project_root: Path) -> dict[str, Any]:
    path = model_settings_path(project_root)
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    if not isinstance(payload, dict):
        return {}
    return {key: value for key, value in payload.items() if key in MODEL_SETTING_FIELDS}


def save_model_settings(project_root: Path, settings: dict[str, Any]) -> Path:
    path = model_settings_path(project_root)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {key: settings[key] for key in MODEL_SETTING_FIELDS if key in settings}
    tmp_path = path.with_suffix(".tmp")
    tmp_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    tmp_path.chmod(0o600)
    tmp_path.replace(path)
    path.chmod(0o600)
    return path


def clear_model_settings(project_root: Path) -> None:
    model_settings_path(project_root).unlink(missing_ok=True)


def _find_repo_root() -> Path:
    return Path(__file__).resolve().parents[3]


def _first_nonempty(*names: str) -> str:
    for name in names:
        value = os.getenv(name)
        if value and value.strip():
            return value.strip()
    return ""


def _optional_env_bool(name: str) -> bool | None:
    value = os.getenv(name)
    if value is None or not value.strip():
        return None
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"{name} must be a boolean value, got {value!r}")


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
    chat_template_enable_thinking: bool | None
    memory_enabled: bool = True
    memory_user_id: str = "default"
    memory_namespace: str = "bioagent"
    max_tool_result_chars: int = 6000
    max_execution_attempts: int = 8
    enable_legacy_workflow_tools: bool = False
    provider: str = "openai_compatible"

    @classmethod
    def from_env(cls, *, include_model_settings: bool = True) -> "AgentConfig":
        repo_root = _find_repo_root()
        project_root = repo_root / "BioAgent"
        mas2_root = repo_root / "mas_2"
        load_dotenv(mas2_root / ".env", override=False)

        logs_dir = project_root / "logs"
        runs_dir = project_root / "runs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        runs_dir.mkdir(parents=True, exist_ok=True)

        memory_enabled = _optional_env_bool("BIOAGENT_MEMORY_ENABLED")
        enable_legacy_workflow_tools = _optional_env_bool("BIOAGENT_ENABLE_LEGACY_WORKFLOW_TOOLS")
        config = cls(
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
            chat_template_enable_thinking=_optional_env_bool("LLM_CHAT_TEMPLATE_ENABLE_THINKING"),
            memory_enabled=True if memory_enabled is None else memory_enabled,
            memory_user_id=os.getenv("BIOAGENT_MEMORY_USER_ID", "default").strip() or "default",
            memory_namespace=os.getenv("BIOAGENT_MEMORY_NAMESPACE", "bioagent").strip() or "bioagent",
            max_tool_result_chars=int(os.getenv("BIOAGENT_MAX_TOOL_RESULT_CHARS", "6000")),
            max_execution_attempts=int(os.getenv("BIOAGENT_MAX_ATTEMPTS", "8")),
            enable_legacy_workflow_tools=False if enable_legacy_workflow_tools is None else enable_legacy_workflow_tools,
            provider=os.getenv("MODEL_PROVIDER", "openai_compatible").strip().lower() or "openai_compatible",
        )
        return apply_model_settings(config) if include_model_settings else config

    def mask_api_key(self) -> str:
        value = self.api_key.strip()
        if not value:
            return "(empty)"
        if len(value) <= 10:
            return value[:2] + "***"
        return f"{value[:6]}...{value[-4:]}"


def apply_model_settings(config: AgentConfig) -> AgentConfig:
    settings = load_model_settings(config.project_root)
    if not settings:
        return config
    return replace(config, **settings)
