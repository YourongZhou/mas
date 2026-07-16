from __future__ import annotations

import os
from pathlib import Path

from bioagent.config import AgentConfig


def allowed_roots(config: AgentConfig, run_dir: Path) -> list[Path]:
    return [
        config.repo_root.resolve(),
        config.mas2_root.resolve(),
        run_dir.resolve(),
    ]


def resolve_allowed_path(config: AgentConfig, run_dir: Path, raw_path: str) -> Path:
    path = Path(raw_path)
    candidates = [path] if path.is_absolute() else [run_dir / path, config.repo_root / path]
    resolved = candidates[0].resolve()
    for candidate in candidates:
        candidate_resolved = candidate.resolve()
        if candidate_resolved.exists():
            resolved = candidate_resolved
            break
    roots = allowed_roots(config, run_dir)
    if not any(resolved == root or str(resolved).startswith(str(root) + os.sep) for root in roots):
        raise ValueError(f"Access denied: {resolved}")
    return resolved


def resolve_data_path(config: AgentConfig, run_dir: Path, raw_path: str) -> Path:
    candidates: list[Path] = []
    p = Path(raw_path)
    if p.is_absolute() and len(p.parts) >= 2 and p.parts[1] == "repo":
        candidates.append(config.repo_root.joinpath(*p.parts[2:]))
    elif p.is_absolute():
        candidates.append(p)
    else:
        candidates.extend(
            [
                config.repo_root / raw_path,
                config.mas2_root / raw_path,
                config.mas2_root / "data" / raw_path,
                config.mas2_root / "data" / p.name,
            ]
        )
    for candidate in candidates:
        try:
            resolved = resolve_allowed_path(config, run_dir, str(candidate))
        except ValueError:
            continue
        if resolved.exists():
            return resolved
    return resolve_allowed_path(config, run_dir, raw_path)


def to_repo_mount_path(config: AgentConfig, host_path: Path) -> str:
    resolved = host_path.resolve()
    repo_root = config.repo_root.resolve()
    rel = resolved.relative_to(repo_root)
    return "/repo/" + rel.as_posix()
