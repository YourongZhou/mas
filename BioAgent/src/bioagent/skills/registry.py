from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


def parse_frontmatter(text: str) -> tuple[dict[str, Any], str]:
    if not text.startswith("---"):
        return {}, text
    parts = text.split("---", 2)
    if len(parts) < 3:
        return {}, text
    raw = parts[1].strip()
    body = parts[2].lstrip()
    try:
        import yaml

        loaded = yaml.safe_load(raw)
        return (loaded if isinstance(loaded, dict) else {}), body
    except Exception:
        meta: dict[str, Any] = {}
        for line in raw.splitlines():
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            meta[key.strip()] = value.strip().strip("\"'")
        return meta, body


@dataclass(frozen=True)
class WorkflowSkill:
    skill_id: str
    name: str
    category: str
    short_description: str
    runtime: str
    env_profile: str
    env_image: str
    path: str


def _skill_from_path(skill_md: Path) -> WorkflowSkill:
    text = skill_md.read_text(encoding="utf-8", errors="replace")
    meta, _body = parse_frontmatter(text)
    skill_id = str(meta.get("id") or skill_md.parent.name).strip()
    return WorkflowSkill(
        skill_id=skill_id,
        name=str(meta.get("name") or skill_id),
        category=str(meta.get("category") or ""),
        short_description=str(meta.get("short-description") or meta.get("short_description") or ""),
        runtime=str(meta.get("runtime") or ""),
        env_profile=str(meta.get("env_profile") or ""),
        env_image=str(meta.get("env_image") or ""),
        path=str(skill_md.parent.resolve()),
    )


def list_workflow_skills(workflows_root: Path) -> list[WorkflowSkill]:
    return [_skill_from_path(path) for path in sorted(workflows_root.glob("*/SKILL.md"))]


def inspect_workflow_skill(workflows_root: Path, skill_id: str, *, max_chars: int = 20000) -> dict[str, Any]:
    skill_md = workflows_root / skill_id / "SKILL.md"
    if not skill_md.is_file():
        return {"error": f"Skill not found: {skill_id}"}
    text = skill_md.read_text(encoding="utf-8", errors="replace")
    meta, body = parse_frontmatter(text)
    scripts = []
    references = []
    for root_name, target in (("scripts", scripts), ("references", references)):
        root = skill_md.parent / root_name
        if root.is_dir():
            target.extend(str(path.relative_to(skill_md.parent)) for path in sorted(root.rglob("*")) if path.is_file())
    return {
        "skill_id": skill_id,
        "path": str(skill_md.parent.resolve()),
        "metadata": meta,
        "body_preview": body[:max_chars],
        "scripts": scripts,
        "references": references,
    }


def load_image_catalog(docker_root: Path) -> list[dict[str, Any]]:
    path = docker_root / "image_catalog.json"
    if not path.is_file():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    return data if isinstance(data, list) else []


def image_for_profile(docker_root: Path, env_profile: str) -> dict[str, Any] | None:
    for entry in load_image_catalog(docker_root):
        if entry.get("env_profile") == env_profile:
            return entry
    return None

