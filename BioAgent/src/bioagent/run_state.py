from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from langchain_core.messages import BaseMessage, messages_from_dict, messages_to_dict

from .config import AgentConfig


@dataclass(frozen=True)
class PendingRunState:
    run_id: str
    run_dir: Path
    log_path: Path
    messages: list[BaseMessage]
    question: str
    metadata: dict[str, Any]
    path: Path


def pending_state_path(config: AgentConfig, run_id: str) -> Path:
    return config.runs_dir / run_id / "state" / "pending_user_input.json"


def save_pending_state(
    *,
    config: AgentConfig,
    run_id: str,
    run_dir: Path,
    log_path: Path,
    messages: list[BaseMessage],
    question: str,
    metadata: dict[str, Any] | None = None,
) -> Path:
    path = pending_state_path(config, run_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "run_id": run_id,
        "run_dir": str(run_dir),
        "log_path": str(log_path),
        "messages": messages_to_dict(messages),
        "question": question,
        "metadata": metadata or {},
    }
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    return path


def load_pending_state(config: AgentConfig, resume_id: str) -> PendingRunState:
    path = pending_state_path(config, resume_id)
    if not path.is_file():
        raise FileNotFoundError(f"Pending run state not found: {path}")
    payload = json.loads(path.read_text(encoding="utf-8"))
    return PendingRunState(
        run_id=str(payload["run_id"]),
        run_dir=Path(str(payload["run_dir"])),
        log_path=Path(str(payload["log_path"])),
        messages=messages_from_dict(payload["messages"]),
        question=str(payload.get("question") or ""),
        metadata=payload.get("metadata") if isinstance(payload.get("metadata"), dict) else {},
        path=path,
    )


def clear_pending_state(config: AgentConfig, resume_id: str) -> None:
    path = pending_state_path(config, resume_id)
    if path.exists():
        path.unlink()
