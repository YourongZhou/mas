"""FastAPI server for sandbox tools."""

from __future__ import annotations

import argparse
import os
from typing import Iterable

from fastapi import Depends, FastAPI, HTTPException, Request
import uvicorn

from .client import SandboxClient, SandboxError
from .schemas import (
    ExecCommandArgs,
    FileExistsArgs,
    GlobSearchArgs,
    GrepTextArgs,
    ListFilesArgs,
    ReadFileArgs,
    ReplaceTextArgs,
    WriteFileArgs,
)


def _split_roots(raw: str | None, fallback: Iterable[str]) -> list[str]:
    if raw and raw.strip():
        return [part.strip() for part in raw.split(":") if part.strip()]
    return [str(item) for item in fallback if str(item).strip()]


def create_auth_dependency(token: str | None):
    async def verify_token(request: Request) -> None:
        if not token:
            return
        auth = request.headers.get("Authorization", "")
        if auth != f"Bearer {token}":
            raise HTTPException(status_code=401, detail="Unauthorized")

    return verify_token


def create_app(
    *,
    token: str | None,
    workspace_root: str,
    allowed_roots: list[str],
    writable_roots: list[str],
    exec_allowlist: list[str] | None = None,
) -> FastAPI:
    app = FastAPI(title="MAS Sandbox Server", version="2.0.0")
    auth = create_auth_dependency(token)
    client = SandboxClient(
        workspace_root=workspace_root,
        allowed_roots=allowed_roots,
        writable_roots=writable_roots,
        backend="local",
        exec_allowlist=exec_allowlist,
    )

    def _wrap(callable_):
        try:
            return callable_()
        except SandboxError as exc:
            raise HTTPException(status_code=400, detail=str(exc)) from exc

    @app.get("/api/health")
    async def health() -> dict:
        return {"status": "ok", **client.describe_backend()}

    @app.post("/api/fs/read-text", dependencies=[Depends(auth)])
    async def read_text(req: ReadFileArgs) -> dict:
        return _wrap(lambda: client.read_file(req.path, line_offset=req.line_offset, max_lines=req.max_lines))

    @app.post("/api/fs/list", dependencies=[Depends(auth)])
    async def list_files(req: ListFilesArgs) -> dict:
        return _wrap(lambda: client.list_files(req.path, recursive=req.recursive, max_entries=req.max_entries))

    @app.post("/api/fs/glob", dependencies=[Depends(auth)])
    async def glob_files(req: GlobSearchArgs) -> dict:
        return _wrap(lambda: client.glob_search(req.pattern, path=req.path, max_entries=req.max_entries))

    @app.post("/api/fs/search", dependencies=[Depends(auth)])
    async def search_files(req: GrepTextArgs) -> dict:
        return _wrap(
            lambda: client.grep_text(
                req.pattern,
                path=req.path,
                glob_filter=req.glob_filter,
                case_insensitive=req.case_insensitive,
                context_lines=req.context_lines,
                head_limit=req.head_limit,
            )
        )

    @app.post("/api/fs/exists", dependencies=[Depends(auth)])
    async def file_exists(req: FileExistsArgs) -> dict:
        return _wrap(lambda: client.file_exists(req.path))

    @app.post("/api/fs/write-text", dependencies=[Depends(auth)])
    async def write_text(req: WriteFileArgs) -> dict:
        return _wrap(lambda: client.write_file(req.path, req.content, create_dirs=req.create_dirs))

    @app.post("/api/fs/replace-text", dependencies=[Depends(auth)])
    async def replace_text(req: ReplaceTextArgs) -> dict:
        return _wrap(lambda: client.replace_text(req.path, req.old_text, req.new_text))

    @app.post("/api/process/run-command", dependencies=[Depends(auth)])
    async def run_command(req: ExecCommandArgs) -> dict:
        return _wrap(lambda: client.exec_command(req.argv, cwd=req.cwd, timeout_s=req.timeout_s))

    return app


def main() -> None:
    parser = argparse.ArgumentParser(description="MAS sandbox server")
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=int(os.environ.get("MAS_SANDBOX_PORT", "8091")))
    parser.add_argument("--token", default=os.environ.get("MAS_SANDBOX_TOKEN") or "")
    parser.add_argument("--workspace-root", default=os.environ.get("MAS_SANDBOX_WORKSPACE_ROOT", "/workspace"))
    args = parser.parse_args()

    allowed_roots = _split_roots(
        os.environ.get("MAS_SANDBOX_ALLOWED_ROOTS"),
        (args.workspace_root, "/mnt/results", "/mnt/memory"),
    )
    writable_roots = _split_roots(
        os.environ.get("MAS_SANDBOX_WRITABLE_ROOTS"),
        (args.workspace_root, "/mnt/results", "/mnt/memory"),
    )
    exec_allowlist = _split_roots(
        os.environ.get("MAS_SANDBOX_EXEC_ALLOWLIST"),
        (),
    ) or None

    app = create_app(
        token=args.token or None,
        workspace_root=args.workspace_root,
        allowed_roots=allowed_roots,
        writable_roots=writable_roots,
        exec_allowlist=exec_allowlist,
    )
    uvicorn.run(app, host=args.host, port=args.port, log_level="info")


if __name__ == "__main__":
    main()
