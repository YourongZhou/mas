from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from .extraction import extract_verified_episode, task_and_data_signatures
from .long_term import EpisodeStore, JsonEpisodeStore
from .retrieval import EpisodeRetriever
from .schemas import MemoryEpisode, TaskState
from .short_term import EXECUTION_TOOLS, ShortTermMemory


@dataclass
class MemoryHarness:
    enabled: bool
    namespace: tuple[str, ...]
    store: EpisodeStore | None = None
    disabled_reason: str = ""
    short_term: ShortTermMemory | None = None
    prior_episodes: list[MemoryEpisode] = field(default_factory=list)
    execution_outcome: str | None = None
    execution_tool: str = "execution"

    @property
    def tools(self) -> list[Any]:
        return []

    def runtime_summary(self) -> str:
        if not self.enabled:
            return f"short_term=enabled long_term=disabled reason={self.disabled_reason or 'not configured'}"
        location = getattr(self.store, "path", "custom-store")
        return f"short_term=enabled long_term=enabled namespace={'/'.join(self.namespace)} store={location}"

    def start_task(
        self,
        task: str,
        *,
        data_path: str | None = None,
        result_dir: str | None = None,
        restored_state: dict[str, Any] | None = None,
    ) -> TaskState:
        if restored_state and isinstance(restored_state.get("task_state"), dict):
            self.short_term = ShortTermMemory.restore(restored_state["task_state"])
            self.execution_outcome = restored_state.get("execution_outcome")
            self.execution_tool = str(restored_state.get("execution_tool") or "execution")
        else:
            self.short_term = ShortTermMemory.start(task=task, data_path=data_path, result_dir=result_dir)
            self.execution_outcome = None
            self.execution_tool = "execution"

        self.prior_episodes = []
        if self.enabled and self.store is not None:
            task_signature, data_signature = task_and_data_signatures(self.short_term.state)
            self.prior_episodes = EpisodeRetriever(self.store).retrieve(task_signature, data_signature, limit=3)
        return self.short_term.state

    def context_text(self) -> str:
        if self.short_term is None:
            return ""
        sections = [self.short_term.compact_summary()]
        if self.enabled and self.store is not None:
            sections.append(EpisodeRetriever(self.store).format_for_prompt(self.prior_episodes))
        else:
            sections.extend(
                [
                    "Prior verified experience:",
                    "Current observations always override historical experience.",
                    "- Long-term memory is disabled.",
                ]
            )
        return "\n\n".join(sections)

    def observe_tool_call(self, tool_name: str, args: dict[str, Any]) -> None:
        if self.short_term is not None:
            self.short_term.observe_tool_call(tool_name, args)

    def apply_user_input(self, user_input: str) -> None:
        if self.short_term is not None:
            self.short_term.apply_user_input(user_input)

    def observe_tool_result(self, tool_name: str, args: dict[str, Any], result: Any) -> None:
        if self.short_term is None:
            return
        self.short_term.observe_tool_result(tool_name, args, result)
        if tool_name not in EXECUTION_TOOLS or not isinstance(result, dict) or result.get("status") in {"running", "cancelled"}:
            return
        self.execution_tool = tool_name
        self.execution_outcome = "success" if _result_succeeded(result) else "failure"

    def snapshot(self) -> dict[str, Any]:
        return {
            "task_state": self.short_term.state.to_dict() if self.short_term else {},
            "execution_outcome": self.execution_outcome,
            "execution_tool": self.execution_tool,
        }

    def public_snapshot(self) -> dict[str, Any]:
        return {
            "taskState": self.short_term.state.to_dict() if self.short_term else {},
            "priorEpisodes": [episode.to_dict() for episode in self.prior_episodes],
            "longTermEnabled": self.enabled,
            "namespace": list(self.namespace),
            "executionOutcome": self.execution_outcome,
            "executionTool": self.execution_tool,
        }

    def mark_paused(self) -> None:
        if self.short_term is not None:
            self.short_term.mark_paused()

    def finish_task(self, *, source_run_id: str, status: str) -> MemoryEpisode | None:
        if status not in {"completed", "failed"} or self.short_term is None:
            return None
        if status == "completed":
            self.short_term.mark_completed()
        else:
            self.short_term.mark_failed()
        episode = extract_verified_episode(
            self.short_term.state,
            source_run_id=source_run_id,
            execution_outcome=self.execution_outcome,
            execution_tool=self.execution_tool,
        )
        if episode is not None and self.enabled and self.store is not None:
            self.store.add(episode)
        return episode


def build_memory_harness(config: Any) -> MemoryHarness:
    namespace = (config.memory_namespace, config.memory_user_id)
    if not config.memory_enabled:
        return MemoryHarness(False, namespace, disabled_reason="BIOAGENT_MEMORY_ENABLED=false")
    filename = "-".join(_safe_name(part) for part in namespace) + "-episodes.json"
    path = Path(config.project_root) / "memory" / filename
    return MemoryHarness(True, namespace, store=JsonEpisodeStore(path))


def _safe_name(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "-", value.strip()).strip("-.")
    return cleaned or "default"


def _result_succeeded(result: dict[str, Any]) -> bool:
    if "ok" in result:
        return bool(result.get("ok"))
    return not bool(result.get("error") or result.get("error_reason"))


__all__ = [
    "MemoryEpisode",
    "MemoryHarness",
    "TaskState",
    "build_memory_harness",
]
