from __future__ import annotations

import glob
import subprocess
from pathlib import Path

from bioagent.config import AgentConfig

from .common import resolve_allowed_path


def list_files_impl(config: AgentConfig, run_dir: Path, path: str = ".", recursive: bool = False, max_entries: int = 200) -> str:
    resolved = resolve_allowed_path(config, run_dir, path)
    if not resolved.exists():
        return f"Path not found: {resolved}"
    if not resolved.is_dir():
        return f"Path is not a directory: {resolved}"
    iterator = resolved.rglob("*") if recursive else resolved.iterdir()
    lines = [f"Directory: {resolved}"]
    for idx, item in enumerate(iterator):
        if idx >= max_entries:
            lines.append("[truncated]")
            break
        kind = "dir" if item.is_dir() else "file"
        size = item.stat().st_size if item.is_file() else 0
        lines.append(f"[{kind}] {item} size={size}")
    return "\n".join(lines)


def read_file_impl(config: AgentConfig, run_dir: Path, path: str, line_offset: int = 1, max_lines: int = 200) -> str:
    resolved = resolve_allowed_path(config, run_dir, path)
    if not resolved.is_file():
        return f"File not found: {resolved}"
    sample = resolved.read_bytes()[:8192]
    if b"\x00" in sample or resolved.suffix.lower() in {".h5ad", ".h5", ".hdf5", ".rds", ".rda", ".loom", ".png", ".jpg", ".jpeg", ".pdf"}:
        return (
            f"File: {resolved}\n"
            f"Binary or structured data file detected; size={resolved.stat().st_size} bytes. "
            "Use a domain execution tool, such as execute_python with scanpy/anndata, to inspect its contents."
        )
    lines = resolved.read_text(encoding="utf-8", errors="replace").splitlines()
    start = max(0, line_offset - 1)
    subset = lines[start : start + max_lines]
    body = "\n".join(f"{idx} | {line}" for idx, line in enumerate(subset, start=start + 1))
    suffix = "\n[truncated]" if start + max_lines < len(lines) else ""
    return f"File: {resolved}\n{body}{suffix}"


def glob_search_impl(config: AgentConfig, run_dir: Path, pattern: str, path: str = ".", max_entries: int = 200) -> str:
    root = resolve_allowed_path(config, run_dir, path)
    matches = sorted(glob.glob(str(root / pattern), recursive=True))
    trimmed = matches[:max_entries]
    return "\n".join(trimmed + (["[truncated]"] if len(matches) > len(trimmed) else [])) or "(no matches)"


def grep_text_impl(
    config: AgentConfig,
    run_dir: Path,
    pattern: str,
    path: str = ".",
    glob_filter: str | None = None,
    case_insensitive: bool = False,
    max_matches: int = 100,
) -> str:
    root = resolve_allowed_path(config, run_dir, path)
    import shutil

    if shutil.which("rg"):
        cmd = ["rg", "-n"]
        if case_insensitive:
            cmd.append("-i")
        if glob_filter:
            cmd.extend(["--glob", glob_filter])
        cmd.extend([pattern, str(root)])
    else:
        cmd = ["findstr", "/S", "/N", "/I" if case_insensitive else "", pattern, str(root / "*")]
        cmd = [part for part in cmd if part]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, check=False)
    text = result.stdout if result.returncode in (0, 1) else (result.stderr or result.stdout)
    lines = text.splitlines()
    return "\n".join(lines[:max_matches]) or "(no matches)"

