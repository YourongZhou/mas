from __future__ import annotations

from pathlib import Path

import pytest

from bioagent.config import AgentConfig
from bioagent.tools.common import resolve_data_path


def _config(tmp_path: Path) -> AgentConfig:
    repo = tmp_path / "repo"
    mas2 = repo / "mas_2"
    bio = repo / "BioAgent"
    (mas2 / "data").mkdir(parents=True)
    (mas2 / "workflows").mkdir(parents=True)
    (mas2 / "docker").mkdir(parents=True)
    (bio / "logs").mkdir(parents=True)
    (bio / "runs").mkdir(parents=True)
    return AgentConfig(
        repo_root=repo,
        project_root=bio,
        mas2_root=mas2,
        workflows_root=mas2 / "workflows",
        docker_root=mas2 / "docker",
        logs_dir=bio / "logs",
        runs_dir=bio / "runs",
        api_key="test-key",
        base_url="http://example.test/v1",
        model_name="test-model",
        temperature=0.0,
        request_timeout=10.0,
        mimo_thinking_type="",
        chat_template_enable_thinking=None,
    )


def test_resolve_data_path_accepts_repo_mount_path(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    data_path = config.mas2_root / "data" / "bmmc_b_cell.h5ad"
    data_path.write_text("placeholder", encoding="utf-8")

    resolved = resolve_data_path(config, run_dir, "/repo/mas_2/data/bmmc_b_cell.h5ad")

    assert resolved == data_path.resolve()


def test_resolve_data_path_still_rejects_unmounted_absolute_path(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    with pytest.raises(ValueError, match="Access denied"):
        resolve_data_path(config, run_dir, "/etc/passwd")
