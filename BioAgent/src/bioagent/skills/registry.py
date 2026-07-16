from __future__ import annotations

import ast
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


def read_skill_script(
    workflows_root: Path,
    skill_id: str,
    script_path: str,
    *,
    line_offset: int = 1,
    max_lines: int = 200,
) -> str:
    path = _resolve_skill_file(workflows_root, skill_id, script_path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = max(line_offset, 1)
    end = min(start + max(max_lines, 1) - 1, len(lines))
    return "\n".join(f"{number}: {lines[number - 1]}" for number in range(start, end + 1))


def inspect_skill_script_symbols(
    workflows_root: Path,
    skill_id: str,
    script_path: str,
    *,
    include_docstrings: bool = False,
) -> dict[str, Any]:
    path = _resolve_skill_file(workflows_root, skill_id, script_path)
    source = path.read_text(encoding="utf-8", errors="replace")
    tree = ast.parse(source)
    return {
        "skill_id": skill_id,
        "script_path": _relative_skill_path(workflows_root, skill_id, path),
        "functions": [
            _function_summary(node, include_docstring=include_docstrings)
            for node in tree.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        ],
        "classes": [_class_summary(node, include_docstrings=include_docstrings) for node in tree.body if isinstance(node, ast.ClassDef)],
    }


def inspect_skill_function(
    workflows_root: Path,
    skill_id: str,
    function_name: str,
    *,
    script_path: str = "",
    max_chars: int = 4000,
) -> dict[str, Any]:
    search_paths = [_resolve_skill_file(workflows_root, skill_id, script_path)] if script_path else _python_script_paths(workflows_root, skill_id)
    for path in search_paths:
        source = path.read_text(encoding="utf-8", errors="replace")
        tree = ast.parse(source)
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name == function_name:
                return {
                    "skill_id": skill_id,
                    "script_path": _relative_skill_path(workflows_root, skill_id, path),
                    "name": node.name,
                    "signature": _function_signature(node),
                    "docstring": ast.get_docstring(node) or "",
                    "start_line": node.lineno,
                    "end_line": getattr(node, "end_lineno", node.lineno),
                    "source_preview": _source_for_node(source, node)[:max_chars],
                }
    return {"error": f"Function not found: {function_name}", "skill_id": skill_id, "searched_paths": [str(path) for path in search_paths]}


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


def _resolve_skill_file(workflows_root: Path, skill_id: str, relative_path: str) -> Path:
    if not relative_path.strip():
        raise ValueError("script_path is required")
    if Path(relative_path).is_absolute():
        raise ValueError("script_path must be relative to the workflow skill directory")
    skill_root = (workflows_root / skill_id).resolve()
    path = (skill_root / relative_path).resolve()
    if skill_root not in path.parents and path != skill_root:
        raise ValueError(f"Access denied outside skill: {relative_path}")
    if not path.is_file():
        raise FileNotFoundError(f"Skill file not found: {relative_path}")
    return path


def _relative_skill_path(workflows_root: Path, skill_id: str, path: Path) -> str:
    return str(path.resolve().relative_to((workflows_root / skill_id).resolve()))


def _python_script_paths(workflows_root: Path, skill_id: str) -> list[Path]:
    scripts_root = (workflows_root / skill_id / "scripts").resolve()
    if not scripts_root.is_dir():
        return []
    return sorted(path for path in scripts_root.rglob("*.py") if path.is_file())


def _function_summary(node: ast.FunctionDef | ast.AsyncFunctionDef, *, include_docstring: bool = False) -> dict[str, Any]:
    summary: dict[str, Any] = {
        "name": node.name,
        "signature": _function_signature(node),
        "start_line": node.lineno,
        "end_line": getattr(node, "end_lineno", node.lineno),
    }
    if include_docstring:
        summary["docstring"] = ast.get_docstring(node) or ""
    return summary


def _class_summary(node: ast.ClassDef, *, include_docstrings: bool = False) -> dict[str, Any]:
    methods = [child for child in node.body if isinstance(child, (ast.FunctionDef, ast.AsyncFunctionDef))]
    summary: dict[str, Any] = {
        "name": node.name,
        "start_line": node.lineno,
        "end_line": getattr(node, "end_lineno", node.lineno),
        "methods": [_function_summary(method, include_docstring=include_docstrings) for method in methods],
    }
    if include_docstrings:
        summary["docstring"] = ast.get_docstring(node) or ""
    return summary


def _function_signature(node: ast.FunctionDef | ast.AsyncFunctionDef) -> str:
    prefix = "async def" if isinstance(node, ast.AsyncFunctionDef) else "def"
    returns = f" -> {ast.unparse(node.returns)}" if node.returns is not None else ""
    return f"{prefix} {node.name}({_format_arguments(node.args)}){returns}"


def _format_arguments(args: ast.arguments) -> str:
    positional = list(args.posonlyargs) + list(args.args)
    defaults = [None] * (len(positional) - len(args.defaults)) + list(args.defaults)
    parts = [_format_arg(arg, default) for arg, default in zip(positional, defaults)]
    if args.posonlyargs:
        parts.insert(len(args.posonlyargs), "/")
    if args.vararg:
        parts.append("*" + _format_arg(args.vararg, None))
    elif args.kwonlyargs:
        parts.append("*")
    parts.extend(_format_arg(arg, default) for arg, default in zip(args.kwonlyargs, args.kw_defaults))
    if args.kwarg:
        parts.append("**" + _format_arg(args.kwarg, None))
    return ", ".join(parts)


def _format_arg(arg: ast.arg, default: ast.expr | None) -> str:
    text = arg.arg
    if arg.annotation is not None:
        text += f": {ast.unparse(arg.annotation)}"
    if default is not None:
        text += f"={ast.unparse(default)}"
    return text


def _source_for_node(source: str, node: ast.AST) -> str:
    lines = source.splitlines()
    start = getattr(node, "lineno", 1)
    end = getattr(node, "end_lineno", start)
    return "\n".join(lines[start - 1 : end])
