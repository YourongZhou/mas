"""LangChain tool bindings for sandbox capabilities."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any, Mapping

from langchain_core.tools import BaseTool, StructuredTool

from src.utils.project_paths import get_mas2_project_root
from src.utils.workflow_skills import get_workflows_root, resolve_workflow_root

from .client import SandboxClient
from .schemas import (
    ExecCommandArgs,
    FileExistsArgs,
    GlobSearchArgs,
    GrepTextArgs,
    ListFilesArgs,
    ReadFileArgs,
)


@dataclass
class SandboxToolBundle:
    client: SandboxClient
    tools: list[BaseTool]
    tool_map: dict[str, BaseTool]
    backend: str
    endpoint: str | None
    allowed_roots: list[str]


def _format_listing(payload: dict) -> str:
    entries = payload.get("entries", [])
    if not entries:
        return f"Directory: {payload.get('path')}\n(no entries)"

    lines = [f"Directory: {payload.get('path')}"]
    for entry in entries:
        lines.append(
            f"[{entry.get('type', 'file')}] {entry.get('path')} "
            f"(size={entry.get('size', 0)})"
        )
    if payload.get("truncated"):
        lines.append("[truncated]")
    return "\n".join(lines)


def _format_read(payload: dict) -> str:
    header = (
        f"File: {payload.get('path')}\n"
        f"Lines: {payload.get('line_offset')}.."
        f"{payload.get('line_offset', 1) + payload.get('line_count', 0) - 1}"
    )
    suffix = "\n[truncated]" if payload.get("truncated") else ""
    return f"{header}\n{payload.get('content', '')}{suffix}".strip()


def _format_glob(payload: dict) -> str:
    matches = payload.get("matches", [])
    lines = [f"Glob root: {payload.get('path')}", f"Pattern: {payload.get('pattern')}"]
    lines.extend(matches or ["(no matches)"])
    if payload.get("truncated"):
        lines.append("[truncated]")
    return "\n".join(lines)


def _format_grep(payload: dict) -> str:
    lines = [
        f"Search path: {payload.get('path')}",
        f"Pattern: {payload.get('pattern')}",
        f"Command: {payload.get('command')}",
        f"Match count: {payload.get('match_count', 0)}",
    ]
    body = payload.get("matches") or "(no matches)"
    lines.append(body)
    if payload.get("truncated"):
        lines.append("[truncated]")
    return "\n".join(lines)


def _format_exists(payload: dict) -> str:
    return (
        f"Path: {payload.get('path')}\n"
        f"exists={payload.get('exists')} "
        f"is_file={payload.get('is_file')} "
        f"is_dir={payload.get('is_dir')}"
    )


def _format_exec(payload: dict) -> str:
    argv = " ".join(payload.get("argv", []))
    stdout = (payload.get("stdout") or "").strip()
    stderr = (payload.get("stderr") or "").strip()
    lines = [
        f"cwd: {payload.get('cwd')}",
        f"argv: {argv}",
        f"exit_code: {payload.get('exit_code')}",
        f"timed_out: {payload.get('timed_out')}",
    ]
    if stdout:
        lines.append("stdout:")
        lines.append(stdout[:6000])
    if stderr:
        lines.append("stderr:")
        lines.append(stderr[:4000])
    return "\n".join(lines)


def build_code_dev_sandbox_tool_bundle(state: Mapping[str, Any]) -> SandboxToolBundle:
    project_root = get_mas2_project_root()
    workflow_root = get_workflows_root()
    skill_id = (state.get("current_step_skill_id") or "").strip()
    skill_root = resolve_workflow_root(skill_id) if skill_id else None

    allowed_roots = [
        str(project_root),
        str(workflow_root),
    ]
    if skill_root:
        allowed_roots.append(str(skill_root))

    for raw_path in (
        state.get("data_path"),
        state.get("result_path"),
    ):
        if raw_path:
            allowed_roots.append(os.path.abspath(str(raw_path)))

    file_paths = state.get("current_step_file_paths") or {}
    for path_group in ("input_files", "output_files"):
        for raw_path in file_paths.get(path_group, []) or []:
            if raw_path:
                allowed_roots.append(os.path.abspath(str(raw_path)))

    deduped_roots = []
    for root in allowed_roots:
        if root and root not in deduped_roots:
            deduped_roots.append(root)

    backend = (os.environ.get("MAS_SANDBOX_BACKEND") or "local").strip().lower()
    endpoint = os.environ.get("MAS_SANDBOX_ENDPOINT")
    token = os.environ.get("MAS_SANDBOX_TOKEN")
    client = SandboxClient(
        workspace_root=str(project_root),
        allowed_roots=deduped_roots,
        writable_roots=[str(project_root / "results")],
        backend=backend,
        endpoint=endpoint,
        token=token,
    )

    def read_file(path: str, line_offset: int = 1, max_lines: int = 200) -> str:
        return _format_read(client.read_file(path, line_offset=line_offset, max_lines=max_lines))

    def list_files(path: str = ".", recursive: bool = False, max_entries: int = 200) -> str:
        return _format_listing(client.list_files(path, recursive=recursive, max_entries=max_entries))

    def glob_search(pattern: str, path: str = ".", max_entries: int = 200) -> str:
        return _format_glob(client.glob_search(pattern, path=path, max_entries=max_entries))

    def grep_text(
        pattern: str,
        path: str = ".",
        glob_filter: str | None = None,
        case_insensitive: bool = False,
        context_lines: int = 0,
        head_limit: int = 100,
    ) -> str:
        return _format_grep(
            client.grep_text(
                pattern,
                path=path,
                glob_filter=glob_filter,
                case_insensitive=case_insensitive,
                context_lines=context_lines,
                head_limit=head_limit,
            )
        )

    def file_exists(path: str) -> str:
        return _format_exists(client.file_exists(path))

    def exec_command(argv: list[str], cwd: str = ".", timeout_s: int = 20) -> str:
        return _format_exec(client.exec_command(argv, cwd=cwd, timeout_s=timeout_s))

    tools = [
        StructuredTool.from_function(
            func=read_file,
            name="read_file",
            description=(
                "Read a UTF-8 text file for code or reference inspection. "
                "Use this after locating a target file."
            ),
            args_schema=ReadFileArgs,
        ),
        StructuredTool.from_function(
            func=list_files,
            name="list_files",
            description=(
                "List files or directories inside the sandbox-visible workspace. "
                "Use this to orient yourself before reading files."
            ),
            args_schema=ListFilesArgs,
        ),
        StructuredTool.from_function(
            func=glob_search,
            name="glob_search",
            description="Find files by glob pattern, such as `**/*.py` or `references/*.md`.",
            args_schema=GlobSearchArgs,
        ),
        StructuredTool.from_function(
            func=grep_text,
            name="grep_text",
            description="Search file contents with ripgrep-style matching.",
            args_schema=GrepTextArgs,
        ),
        StructuredTool.from_function(
            func=file_exists,
            name="file_exists",
            description="Check whether a file or directory exists.",
            args_schema=FileExistsArgs,
        ),
        StructuredTool.from_function(
            func=exec_command,
            name="exec_command",
            description=(
                "Run a read-oriented shell command using structured argv. "
                "Only a small allowlist of inspection commands is permitted."
            ),
            args_schema=ExecCommandArgs,
        ),
    ]

    return SandboxToolBundle(
        client=client,
        tools=tools,
        tool_map={tool.name: tool for tool in tools},
        backend=client.backend,
        endpoint=client.endpoint,
        allowed_roots=deduped_roots,
    )
