from __future__ import annotations

from typing import Any

from bioagent.config import AgentConfig
from bioagent.skills import (
    inspect_skill_function,
    inspect_skill_script_symbols,
    inspect_workflow_skill,
    list_workflow_skills,
    load_image_catalog,
    read_skill_script,
)


def list_workflow_skills_impl(config: AgentConfig, detail: str = "compact") -> list[dict[str, Any]]:
    skills = list_workflow_skills(config.workflows_root)
    if detail == "full":
        return [skill.__dict__ for skill in skills]
    compact = []
    for skill in skills:
        item = {
            "skill_id": skill.skill_id,
            "name": skill.name,
        }
        if skill.short_description:
            item["short_description"] = _clip(skill.short_description, 80)
        if skill.runtime:
            item["runtime"] = skill.runtime
        if skill.env_profile:
            item["env_profile"] = skill.env_profile
        compact.append(item)
    return compact


def inspect_workflow_skill_impl(config: AgentConfig, skill_id: str, max_chars: int = 3000) -> dict[str, Any]:
    return inspect_workflow_skill(config.workflows_root, skill_id, max_chars=max_chars)


def read_skill_script_impl(
    config: AgentConfig,
    skill_id: str,
    script_path: str,
    line_offset: int = 1,
    max_lines: int = 200,
) -> str:
    return read_skill_script(
        config.workflows_root,
        skill_id,
        script_path,
        line_offset=line_offset,
        max_lines=max_lines,
    )


def inspect_skill_script_symbols_impl(
    config: AgentConfig,
    skill_id: str,
    script_path: str,
    include_docstrings: bool = False,
) -> dict[str, Any]:
    return inspect_skill_script_symbols(
        config.workflows_root,
        skill_id,
        script_path,
        include_docstrings=include_docstrings,
    )


def inspect_skill_function_impl(
    config: AgentConfig,
    skill_id: str,
    function_name: str,
    script_path: str = "",
    max_chars: int = 4000,
) -> dict[str, Any]:
    return inspect_skill_function(
        config.workflows_root,
        skill_id,
        function_name,
        script_path=script_path,
        max_chars=max_chars,
    )


def inspect_image_catalog_impl(config: AgentConfig) -> list[dict[str, Any]]:
    return load_image_catalog(config.docker_root)


def _clip(text: str, max_chars: int) -> str:
    return text if len(text) <= max_chars else text[: max_chars - 3].rstrip() + "..."
