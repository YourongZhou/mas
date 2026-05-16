"""Typed request models for sandbox-backed tools."""

from __future__ import annotations

from pydantic import BaseModel, Field


class ReadFileArgs(BaseModel):
    path: str = Field(..., description="Absolute or workspace-relative file path.")
    line_offset: int = Field(
        default=1,
        ge=1,
        description="1-indexed starting line to read from.",
    )
    max_lines: int = Field(
        default=200,
        ge=1,
        le=1000,
        description="Maximum number of lines to return.",
    )


class ListFilesArgs(BaseModel):
    path: str = Field(default=".", description="Directory path to inspect.")
    recursive: bool = Field(default=False, description="Whether to recurse into subdirectories.")
    max_entries: int = Field(
        default=200,
        ge=1,
        le=2000,
        description="Maximum number of entries to return.",
    )


class GlobSearchArgs(BaseModel):
    pattern: str = Field(..., description="Glob pattern such as `**/*.py`.")
    path: str = Field(default=".", description="Directory path used as the glob root.")
    max_entries: int = Field(
        default=200,
        ge=1,
        le=2000,
        description="Maximum number of matches to return.",
    )


class GrepTextArgs(BaseModel):
    pattern: str = Field(..., description="Regular expression or plain-text search pattern.")
    path: str = Field(default=".", description="File or directory path to search.")
    glob_filter: str | None = Field(
        default=None,
        description="Optional glob filter such as `*.py` or `**/*.md`.",
    )
    case_insensitive: bool = Field(default=False, description="Whether matching ignores case.")
    context_lines: int = Field(default=0, ge=0, le=20, description="Context lines around matches.")
    head_limit: int = Field(
        default=100,
        ge=1,
        le=1000,
        description="Maximum number of output lines to return.",
    )


class FileExistsArgs(BaseModel):
    path: str = Field(..., description="Absolute or workspace-relative path.")


class ExecCommandArgs(BaseModel):
    argv: list[str] = Field(
        ...,
        min_length=1,
        description="Command and arguments as a structured argv list.",
    )
    cwd: str = Field(default=".", description="Working directory for the command.")
    timeout_s: int = Field(
        default=20,
        ge=1,
        le=60,
        description="Execution timeout in seconds.",
    )


class WriteFileArgs(BaseModel):
    path: str = Field(..., description="Absolute or workspace-relative file path.")
    content: str = Field(..., description="Full text content to write.")
    create_dirs: bool = Field(default=True, description="Whether missing parent directories are created.")


class ReplaceTextArgs(BaseModel):
    path: str = Field(..., description="Absolute or workspace-relative file path.")
    old_text: str = Field(..., min_length=1, description="First text block to replace.")
    new_text: str = Field(..., description="Replacement text.")
