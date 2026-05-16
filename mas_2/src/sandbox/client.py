"""Sandbox client with pluggable local or HTTP backends."""

from __future__ import annotations

import glob as glob_module
import os
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Iterable

import requests


DEFAULT_EXEC_ALLOWLIST = (
    "ls",
    "find",
    "rg",
    "grep",
    "head",
    "tail",
    "sed",
    "cat",
    "wc",
    "stat",
    "file",
    "pwd",
    "realpath",
)


class SandboxError(RuntimeError):
    """Raised when sandbox operations fail validation or execution."""


class SandboxClient:
    """Tool-facing sandbox client.

    The local backend is used immediately by Code Dev. The HTTP backend is kept
    compatible with the new sandbox server so the same tool schema can later be
    pointed at an in-container sidecar without changing the graph layer.
    """

    def __init__(
        self,
        *,
        workspace_root: str,
        allowed_roots: Iterable[str],
        writable_roots: Iterable[str] | None = None,
        backend: str = "local",
        endpoint: str | None = None,
        token: str | None = None,
        exec_allowlist: Iterable[str] | None = None,
    ) -> None:
        self.workspace_root = Path(workspace_root).resolve()
        self.allowed_roots = self._normalize_roots(allowed_roots)
        self.writable_roots = self._normalize_roots(writable_roots or [])
        self.backend = backend
        self.endpoint = endpoint.rstrip("/") if endpoint else None
        self.token = token
        self.exec_allowlist = tuple(exec_allowlist or DEFAULT_EXEC_ALLOWLIST)

        if self.backend == "http" and not self.endpoint:
            raise SandboxError("HTTP sandbox backend requires an endpoint.")

    @staticmethod
    def _normalize_roots(roots: Iterable[str]) -> tuple[Path, ...]:
        out: list[Path] = []
        for root in roots:
            if not root:
                continue
            path = Path(root).expanduser().resolve()
            if path not in out:
                out.append(path)
        return tuple(out)

    def _resolve_path(self, raw_path: str, *, allow_write: bool = False) -> Path:
        candidate = Path(raw_path).expanduser()
        if not candidate.is_absolute():
            candidate = self.workspace_root / candidate
        resolved = candidate.resolve()

        roots = self.writable_roots if allow_write else self.allowed_roots
        if not roots:
            roots = (self.workspace_root,)
        if not any(resolved == root or str(resolved).startswith(str(root) + os.sep) for root in roots):
            raise SandboxError(
                f"Access denied for path `{resolved}`. "
                f"Allowed roots: {', '.join(str(root) for root in roots)}"
            )
        return resolved

    def _http_headers(self) -> dict[str, str]:
        headers = {"Content-Type": "application/json"}
        if self.token:
            headers["Authorization"] = f"Bearer {self.token}"
        return headers

    def _post(self, route: str, payload: dict) -> dict:
        if self.backend != "http" or not self.endpoint:
            raise SandboxError("HTTP requests are only available on the HTTP backend.")
        response = requests.post(
            f"{self.endpoint}{route}",
            json=payload,
            headers=self._http_headers(),
            timeout=30,
        )
        try:
            response.raise_for_status()
        except requests.HTTPError as exc:
            detail = response.text.strip()
            raise SandboxError(f"Sandbox HTTP error: {detail or exc}") from exc
        return response.json()

    def describe_backend(self) -> dict[str, str | list[str] | None]:
        return {
            "backend": self.backend,
            "endpoint": self.endpoint,
            "allowed_roots": [str(root) for root in self.allowed_roots],
        }

    def read_file(self, path: str, *, line_offset: int = 1, max_lines: int = 200) -> dict:
        if self.backend == "http":
            return self._post(
                "/api/fs/read-text",
                {"path": path, "line_offset": line_offset, "max_lines": max_lines},
            )

        resolved = self._resolve_path(path)
        if not resolved.exists():
            raise SandboxError(f"File not found: {path}")
        if not resolved.is_file():
            raise SandboxError(f"Path is not a file: {path}")

        with resolved.open("r", encoding="utf-8", errors="replace") as handle:
            lines = handle.readlines()

        start = max(0, line_offset - 1)
        end = min(len(lines), start + max_lines)
        subset = lines[start:end]
        numbered = [f"{idx} | {line.rstrip()}" for idx, line in enumerate(subset, start=start + 1)]
        return {
            "path": str(resolved),
            "content": "\n".join(numbered),
            "line_offset": start + 1,
            "line_count": len(subset),
            "truncated": end < len(lines),
            "total_lines": len(lines),
        }

    def list_files(self, path: str = ".", *, recursive: bool = False, max_entries: int = 200) -> dict:
        if self.backend == "http":
            return self._post(
                "/api/fs/list",
                {"path": path, "recursive": recursive, "max_entries": max_entries},
            )

        resolved = self._resolve_path(path)
        if not resolved.exists():
            raise SandboxError(f"Path not found: {path}")
        if not resolved.is_dir():
            raise SandboxError(f"Path is not a directory: {path}")

        entries: list[dict] = []
        iterator = resolved.rglob("*") if recursive else resolved.iterdir()
        for item in iterator:
            if len(entries) >= max_entries:
                break
            try:
                stat = item.stat()
            except OSError:
                continue
            entries.append(
                {
                    "path": str(item),
                    "name": item.name,
                    "type": "directory" if item.is_dir() else "file",
                    "size": stat.st_size,
                }
            )
        return {
            "path": str(resolved),
            "entries": entries,
            "truncated": len(entries) >= max_entries,
        }

    def glob_search(self, pattern: str, *, path: str = ".", max_entries: int = 200) -> dict:
        if self.backend == "http":
            return self._post(
                "/api/fs/glob",
                {"pattern": pattern, "path": path, "max_entries": max_entries},
            )

        resolved = self._resolve_path(path)
        if not resolved.is_dir():
            raise SandboxError(f"Path is not a directory: {path}")

        full_pattern = str(resolved / pattern)
        matches = sorted(glob_module.glob(full_pattern, recursive=True))
        trimmed = matches[:max_entries]
        for match in trimmed:
            self._resolve_path(match)
        return {
            "path": str(resolved),
            "pattern": pattern,
            "matches": trimmed,
            "truncated": len(matches) > len(trimmed),
        }

    def grep_text(
        self,
        pattern: str,
        *,
        path: str = ".",
        glob_filter: str | None = None,
        case_insensitive: bool = False,
        context_lines: int = 0,
        head_limit: int = 100,
    ) -> dict:
        if self.backend == "http":
            return self._post(
                "/api/fs/search",
                {
                    "pattern": pattern,
                    "path": path,
                    "glob_filter": glob_filter,
                    "case_insensitive": case_insensitive,
                    "context_lines": context_lines,
                    "head_limit": head_limit,
                },
            )

        resolved = self._resolve_path(path)
        cmd = ["rg", "-n"]
        if case_insensitive:
            cmd.append("-i")
        if context_lines:
            cmd.extend(["-C", str(context_lines)])
        if glob_filter:
            cmd.extend(["--glob", glob_filter])
        cmd.extend([pattern, str(resolved)])

        if shutil.which("rg") is None:
            cmd = ["grep", "-rn"]
            if case_insensitive:
                cmd.append("-i")
            cmd.extend([pattern, str(resolved)])

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=30, check=False)
        except subprocess.TimeoutExpired as exc:
            raise SandboxError(f"Search timed out after 30s: {exc}") from exc

        output = result.stdout if result.returncode in (0, 1) else (result.stderr or result.stdout)
        lines = [line for line in output.splitlines() if line.strip()]
        trimmed = lines[:head_limit]
        return {
            "path": str(resolved),
            "pattern": pattern,
            "matches": "\n".join(trimmed),
            "match_count": len(lines),
            "truncated": len(lines) > len(trimmed),
            "command": shlex.join(cmd),
        }

    def file_exists(self, path: str) -> dict:
        if self.backend == "http":
            return self._post("/api/fs/exists", {"path": path})

        resolved = self._resolve_path(path)
        exists = resolved.exists()
        return {
            "path": str(resolved),
            "exists": exists,
            "is_file": resolved.is_file() if exists else False,
            "is_dir": resolved.is_dir() if exists else False,
        }

    def exec_command(self, argv: list[str], *, cwd: str = ".", timeout_s: int = 20) -> dict:
        if self.backend == "http":
            return self._post(
                "/api/process/run-command",
                {"argv": argv, "cwd": cwd, "timeout_s": timeout_s},
            )

        if not argv:
            raise SandboxError("Command argv must not be empty.")
        if argv[0] not in self.exec_allowlist:
            raise SandboxError(
                f"Command `{argv[0]}` is not allowed. "
                f"Allowed commands: {', '.join(self.exec_allowlist)}"
            )

        resolved_cwd = self._resolve_path(cwd)
        if not resolved_cwd.is_dir():
            raise SandboxError(f"Working directory is not a directory: {cwd}")

        try:
            result = subprocess.run(
                argv,
                cwd=resolved_cwd,
                capture_output=True,
                text=True,
                timeout=timeout_s,
                check=False,
            )
        except subprocess.TimeoutExpired as exc:
            return {
                "argv": argv,
                "cwd": str(resolved_cwd),
                "stdout": exc.stdout or "",
                "stderr": (exc.stderr or "") + f"\nCommand timed out after {timeout_s} seconds",
                "exit_code": -1,
                "timed_out": True,
            }

        return {
            "argv": argv,
            "cwd": str(resolved_cwd),
            "stdout": result.stdout,
            "stderr": result.stderr,
            "exit_code": result.returncode,
            "timed_out": False,
        }

    def write_file(self, path: str, content: str, *, create_dirs: bool = True) -> dict:
        if self.backend == "http":
            return self._post(
                "/api/fs/write-text",
                {"path": path, "content": content, "create_dirs": create_dirs},
            )

        resolved = self._resolve_path(path, allow_write=True)
        if create_dirs:
            resolved.parent.mkdir(parents=True, exist_ok=True)
        with resolved.open("w", encoding="utf-8") as handle:
            bytes_written = handle.write(content)
        return {"path": str(resolved), "bytes_written": bytes_written}

    def replace_text(self, path: str, old_text: str, new_text: str) -> dict:
        if self.backend == "http":
            return self._post(
                "/api/fs/replace-text",
                {"path": path, "old_text": old_text, "new_text": new_text},
            )

        resolved = self._resolve_path(path, allow_write=True)
        if not resolved.exists():
            raise SandboxError(f"File not found: {path}")
        content = resolved.read_text(encoding="utf-8")
        if old_text not in content:
            raise SandboxError(f"old_text not found in {path}")
        updated = content.replace(old_text, new_text, 1)
        resolved.write_text(updated, encoding="utf-8")
        return {"path": str(resolved), "updated": True}
