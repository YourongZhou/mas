from __future__ import annotations

from typing import Any

from bioagent.config import AgentConfig
from bioagent.skills import inspect_workflow_skill, list_workflow_skills, load_image_catalog


def list_workflow_skills_impl(config: AgentConfig) -> list[dict[str, Any]]:
    return [skill.__dict__ for skill in list_workflow_skills(config.workflows_root)]


def inspect_workflow_skill_impl(config: AgentConfig, skill_id: str, max_chars: int = 20000) -> dict[str, Any]:
    return inspect_workflow_skill(config.workflows_root, skill_id, max_chars=max_chars)


def inspect_image_catalog_impl(config: AgentConfig) -> list[dict[str, Any]]:
    return load_image_catalog(config.docker_root)

