from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from typing import Any


@dataclass
class TaskState:
    task: str = ""
    current_goal: str = ""
    confirmed_inputs: list[str] = field(default_factory=list)
    expected_outputs: list[str] = field(default_factory=list)
    selected_skill: str = ""
    runtime_environment: str = ""
    active_job_id: str = ""
    job_status: str = ""
    current_stage: str = "understanding"
    active_errors: list[str] = field(default_factory=list)
    resolved_errors: list[str] = field(default_factory=list)
    generated_artifacts: list[str] = field(default_factory=list)
    blockers: list[str] = field(default_factory=list)
    next_action: str = "Understand the task and choose the next verified action."

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "TaskState":
        known = {item.name for item in fields(cls)}
        return cls(**{key: value[key] for key in known if key in value})


@dataclass(frozen=True)
class MemoryEpisode:
    task_signature: str
    data_signature: str
    skill_id: str
    runtime_environment: str
    outcome: str
    verified_root_cause: str
    verified_fix: str
    reusable_lesson: str
    source_run_id: str
    timestamp: str

    def to_dict(self) -> dict[str, str]:
        return asdict(self)

    @classmethod
    def from_dict(cls, value: dict[str, Any]) -> "MemoryEpisode":
        return cls(**{name: str(value.get(name) or "") for name in cls.field_names()})

    @classmethod
    def field_names(cls) -> tuple[str, ...]:
        return tuple(item.name for item in fields(cls))
