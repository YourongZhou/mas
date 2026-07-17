from __future__ import annotations

import subprocess
from pathlib import Path

from bioagent.config import AgentConfig
from bioagent.docker_runner import DockerRunner
from bioagent.logging_utils import RunLogger


def _config(tmp_path: Path) -> AgentConfig:
    repo = tmp_path / "repo"
    mas2 = repo / "mas_2"
    bio = repo / "BioAgent"
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


def test_python_timeout_decodes_output_and_removes_exact_container(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_timeout"
    run_dir.mkdir()
    calls: list[list[str]] = []

    with RunLogger(config.logs_dir, run_id="run_timeout") as logger:
        runner = DockerRunner(config, logger, run_dir)
        monkeypatch.setattr(runner, "_docker_available", lambda: True)
        monkeypatch.setattr(runner, "_image", lambda env_profile, runtime: "mas/py-general:v1")
        monkeypatch.setattr(runner, "_image_exists", lambda image: True)

        def fake_run(command, **kwargs):
            calls.append(command)
            if command[:2] == ["docker", "run"]:
                cidfile = Path(command[command.index("--cidfile") + 1])
                cidfile.write_text("container-timeout\n", encoding="utf-8")
                raise subprocess.TimeoutExpired(
                    command,
                    timeout=1,
                    output=b"partial stdout\n",
                    stderr=b"partial stderr\n",
                )
            return subprocess.CompletedProcess(command, 0, "container-timeout\n", "")

        monkeypatch.setattr(subprocess, "run", fake_run)

        result = runner.execute_python("print('hello')", timeout_s=1)

    assert result.ok is False
    assert result.stdout == "partial stdout\n"
    assert result.stderr == "partial stderr\n\nTimed out"
    assert calls[-1] == ["docker", "rm", "-f", "container-timeout"]
    assert not (run_dir / ".docker-python.cid").exists()
