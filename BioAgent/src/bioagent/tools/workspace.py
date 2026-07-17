from __future__ import annotations

import ast
from pathlib import Path


def resolve_run_path(run_dir: Path, path: str) -> Path:
    root = run_dir.resolve()
    raw = str(path or "").strip()
    if not raw:
        raise ValueError("A run workspace path is required.")
    if raw == "/work":
        candidate = root
    elif raw.startswith("/work/"):
        candidate = root / raw.removeprefix("/work/")
    else:
        value = Path(raw)
        candidate = value if value.is_absolute() else root / value
    resolved = candidate.resolve()
    if resolved != root and root not in resolved.parents:
        raise ValueError(f"Path is outside the run workspace: {path}")
    return resolved


def run_relative_path(run_dir: Path, path: Path) -> str:
    return path.resolve().relative_to(run_dir.resolve()).as_posix()


def write_run_file_impl(run_dir: Path, *, path: str, content: str, overwrite: bool = True) -> dict:
    target = resolve_run_path(run_dir, path)
    if target.exists() and not overwrite:
        return {"ok": False, "error": f"File already exists: {target}"}
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(content, encoding="utf-8")
    return {
        "ok": True,
        "path": str(target),
        "workspace_path": f"/work/{run_relative_path(run_dir, target)}",
        "bytes": target.stat().st_size,
        "lines": content.count("\n") + (1 if content and not content.endswith("\n") else 0),
    }


def edit_run_file_impl(
    run_dir: Path,
    *,
    path: str,
    old_text: str,
    new_text: str,
    replace_all: bool = False,
) -> dict:
    target = resolve_run_path(run_dir, path)
    if not target.is_file():
        return {"ok": False, "error": f"Run workspace file not found: {target}"}
    if not old_text:
        return {"ok": False, "error": "old_text must not be empty."}
    content = target.read_text(encoding="utf-8")
    matches = content.count(old_text)
    if matches == 0:
        return {"ok": False, "error": "old_text was not found; read the current file before editing."}
    if matches > 1 and not replace_all:
        return {
            "ok": False,
            "error": f"old_text matched {matches} locations; provide more context or set replace_all=true.",
        }
    updated = content.replace(old_text, new_text, -1 if replace_all else 1)
    target.write_text(updated, encoding="utf-8")
    return {
        "ok": True,
        "path": str(target),
        "workspace_path": f"/work/{run_relative_path(run_dir, target)}",
        "replacements": matches if replace_all else 1,
        "bytes": target.stat().st_size,
    }


def validate_script_impl(run_dir: Path, *, path: str, runtime: str = "python") -> dict:
    target = resolve_run_path(run_dir, path)
    if not target.is_file():
        return {"ok": False, "error": f"Script not found: {target}"}
    normalized = runtime.strip().lower()
    if normalized == "python":
        try:
            ast.parse(target.read_text(encoding="utf-8"), filename=str(target))
        except SyntaxError as exc:
            return {
                "ok": False,
                "runtime": normalized,
                "path": str(target),
                "line": exc.lineno,
                "offset": exc.offset,
                "error": f"{exc.__class__.__name__}: {exc.msg}",
            }
        return {"ok": True, "runtime": normalized, "path": str(target), "check": "python_ast_parse"}
    if normalized == "r":
        text = target.read_text(encoding="utf-8")
        if not text.strip():
            return {"ok": False, "runtime": normalized, "path": str(target), "error": "R script is empty."}
        return {
            "ok": True,
            "runtime": normalized,
            "path": str(target),
            "check": "nonempty_only",
            "warning": "Full R parsing occurs inside the selected Docker environment when the job starts.",
        }
    return {"ok": False, "error": f"Unsupported runtime: {runtime}. Use python or r."}
