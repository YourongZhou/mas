from __future__ import annotations

import json
import subprocess
from dataclasses import replace
from pathlib import Path

import pytest

from bioagent.config import AgentConfig
from bioagent.tools.artifacts import ArtifactStore
from bioagent.tools.jobs import DockerJobManager
from bioagent.tools.registry import build_tools
from bioagent.tools.workspace import edit_run_file_impl, write_run_file_impl


def _config(tmp_path: Path) -> AgentConfig:
    repo = tmp_path / "repo"
    mas2 = repo / "mas_2"
    bio = repo / "BioAgent"
    (mas2 / "workflows").mkdir(parents=True)
    (mas2 / "docker").mkdir(parents=True)
    (bio / "logs").mkdir(parents=True)
    (bio / "runs").mkdir(parents=True)
    (mas2 / "docker" / "image_catalog.json").write_text(
        json.dumps(
            [
                {
                    "env_profile": "py-general-v1",
                    "runtime": "python",
                    "image_tag": "mas/py-general:v1",
                }
            ]
        ),
        encoding="utf-8",
    )
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


def test_run_workspace_write_and_exact_edit_are_path_safe(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"

    written = write_run_file_impl(run_dir, path="scripts/analysis.py", content="print('one')\n")
    edited = edit_run_file_impl(
        run_dir,
        path="scripts/analysis.py",
        old_text="one",
        new_text="two",
    )

    assert written["ok"] is True
    assert edited["replacements"] == 1
    assert (run_dir / "scripts" / "analysis.py").read_text(encoding="utf-8") == "print('two')\n"
    with pytest.raises(ValueError, match="outside the run workspace"):
        write_run_file_impl(run_dir, path="../escape.py", content="bad")


def test_execution_attempt_budget_persists_across_tool_rebuilds(tmp_path: Path, monkeypatch) -> None:
    config = replace(_config(tmp_path), max_execution_attempts=1)
    run_dir = config.runs_dir / "run_budget"
    calls = 0

    def fake_execute(*args, **kwargs):
        nonlocal calls
        calls += 1
        return {"ok": True, "stdout": "ok", "stderr": ""}

    monkeypatch.setattr("bioagent.tools.registry.execute_python_impl", fake_execute)
    first_tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=run_dir)}
    second_tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=run_dir)}

    assert first_tools["execute_python"].invoke({"code": "print(1)"})["ok"] is True
    blocked = second_tools["execute_python"].invoke({"code": "print(2)"})

    assert blocked["ok"] is False
    assert blocked["attempts_used"] == 1
    assert calls == 1


def test_async_job_is_persisted_and_can_be_polled_after_manager_rebuild(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_job"
    write_run_file_impl(run_dir, path="scripts/analysis.py", content="print('ok')\n")
    calls: list[list[str]] = []

    def fake_run(command, **kwargs):
        calls.append(command)
        if command[:3] == ["docker", "image", "inspect"]:
            return subprocess.CompletedProcess(command, 0, "[]", "")
        if command[:3] == ["docker", "run", "-d"]:
            return subprocess.CompletedProcess(command, 0, "container-123\n", "")
        if command[:2] == ["docker", "inspect"]:
            state = {"Status": "exited", "Running": False, "ExitCode": 0, "Error": ""}
            return subprocess.CompletedProcess(command, 0, json.dumps(state), "")
        if command[:2] == ["docker", "logs"]:
            return subprocess.CompletedProcess(command, 0, "analysis complete\n", "")
        if command[:3] == ["docker", "rm", "-f"]:
            return subprocess.CompletedProcess(command, 0, "container-123\n", "")
        raise AssertionError(command)

    monkeypatch.setattr("bioagent.tools.jobs.shutil.which", lambda name: "/usr/bin/docker")
    monkeypatch.setattr("bioagent.tools.jobs.subprocess.run", fake_run)

    started = DockerJobManager(config, logger=None, run_dir=run_dir).start_job(
        runtime="python",
        script_path="scripts/analysis.py",
        env_profile="py-general-v1",
    )
    rebuilt = DockerJobManager(config, logger=None, run_dir=run_dir)
    listed = rebuilt.list_jobs()
    polled = rebuilt.poll_job("")

    assert started["status"] == "running"
    assert json.loads((run_dir / "state" / "job_state.json").read_text(encoding="utf-8"))["active_job_id"] == started["job_id"]
    assert listed["active_job_id"] == started["job_id"]
    assert listed["jobs"][0]["job_id"] == started["job_id"]
    assert polled["ok"] is True
    assert polled["status"] == "completed"
    assert polled["exit_code"] == 0
    assert polled["stdout_tail"] == "analysis complete\n"
    assert (run_dir / "jobs" / started["job_id"] / "job.json").is_file()
    assert any("/work/scripts/analysis.py" in part for part in calls[1])


def test_job_manager_rejects_invented_ids_and_exposes_active_job(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_active_job"
    write_run_file_impl(run_dir, path="scripts/analysis.py", content="print('ok')\n")
    docker_runs = 0

    def fake_run(command, **kwargs):
        nonlocal docker_runs
        if command[:3] == ["docker", "image", "inspect"]:
            return subprocess.CompletedProcess(command, 0, "[]", "")
        if command[:3] == ["docker", "run", "-d"]:
            docker_runs += 1
            return subprocess.CompletedProcess(command, 0, "container-running\n", "")
        if command[:2] == ["docker", "inspect"]:
            state = {"Status": "running", "Running": True, "ExitCode": 0, "Error": ""}
            return subprocess.CompletedProcess(command, 0, json.dumps(state), "")
        if command[:2] == ["docker", "logs"]:
            return subprocess.CompletedProcess(command, 0, "step 2 of 5\n", "warning line\n")
        raise AssertionError(command)

    monkeypatch.setattr("bioagent.tools.jobs.shutil.which", lambda name: "/usr/bin/docker")
    monkeypatch.setattr("bioagent.tools.jobs.subprocess.run", fake_run)
    manager = DockerJobManager(config, logger=None, run_dir=run_dir)
    started = manager.start_job(runtime="python", script_path="scripts/analysis.py", env_profile="py-general-v1")

    invented = manager.poll_job("job_17")
    running = manager.poll_job("")
    duplicate = manager.start_job(runtime="python", script_path="scripts/analysis.py", env_profile="py-general-v1")

    assert invented["ok"] is False
    assert invented["status"] == "invalid_job_id"
    assert invented["active_job_id"] == started["job_id"]
    assert running["status"] == "running"
    assert running["stdout_tail"] == "step 2 of 5\n"
    assert running["stderr_tail"] == "warning line\n"
    assert duplicate["ok"] is False
    assert duplicate["status"] == "job_already_running"
    assert duplicate["active_job_id"] == started["job_id"]
    assert docker_runs == 1


def test_job_and_artifacts_remain_accessible_from_a_follow_up_run(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    previous_run = config.runs_dir / "run_previous"
    current_run = config.runs_dir / "run_current"
    write_run_file_impl(previous_run, path="scripts/analysis.py", content="print('ok')\n")
    write_run_file_impl(current_run, path="scripts/analysis.py", content="print('duplicate')\n")
    state = {"running": True}

    def fake_run(command, **kwargs):
        if command[:3] == ["docker", "image", "inspect"]:
            return subprocess.CompletedProcess(command, 0, "[]", "")
        if command[:3] == ["docker", "run", "-d"]:
            return subprocess.CompletedProcess(command, 0, "container-session\n", "")
        if command[:2] == ["docker", "inspect"]:
            docker_state = {
                "Status": "running" if state["running"] else "exited",
                "Running": state["running"],
                "ExitCode": 0,
                "Error": "",
            }
            return subprocess.CompletedProcess(command, 0, json.dumps(docker_state), "")
        if command[:2] == ["docker", "logs"]:
            return subprocess.CompletedProcess(command, 0, "analysis complete\n", "")
        if command[:3] == ["docker", "rm", "-f"]:
            return subprocess.CompletedProcess(command, 0, "container-session\n", "")
        raise AssertionError(command)

    monkeypatch.setattr("bioagent.tools.jobs.shutil.which", lambda name: "/usr/bin/docker")
    monkeypatch.setattr("bioagent.tools.jobs.subprocess.run", fake_run)
    started = DockerJobManager(config, logger=None, run_dir=previous_run).start_job(
        runtime="python",
        script_path="scripts/analysis.py",
        env_profile="py-general-v1",
    )
    follow_up = DockerJobManager(config, logger=None, run_dir=current_run, prior_run_dirs=[previous_run])

    listed = follow_up.list_jobs()
    running = follow_up.poll_job("")
    duplicate = follow_up.start_job(
        runtime="python",
        script_path="scripts/analysis.py",
        env_profile="py-general-v1",
    )

    assert listed["active_job_id"] == started["job_id"]
    assert listed["jobs"][0]["run_id"] == previous_run.name
    assert running["status"] == "running"
    assert duplicate["status"] == "job_already_running"

    output = previous_run / "outputs" / "summary.json"
    output.parent.mkdir(parents=True)
    output.write_text('{"status": "complete"}\n', encoding="utf-8")
    state["running"] = False
    completed = follow_up.poll_job(started["job_id"])
    evidence = ArtifactStore(current_run, config=config, prior_run_dirs=[previous_run]).inspect(str(output))

    assert completed["status"] == "completed"
    assert completed["artifacts"][0]["path"] == str(output)
    assert evidence["ok"] is True
    assert evidence["path"] == str(output)


def test_artifact_evidence_is_required_and_revalidated_at_finish(tmp_path: Path) -> None:
    run_dir = tmp_path / "run"
    outputs = run_dir / "outputs"
    outputs.mkdir(parents=True)
    table = outputs / "summary.csv"
    table.write_text("metric,value\ncells,10\ngenes,20\n", encoding="utf-8")
    store = ArtifactStore(run_dir)

    evidence = store.inspect("outputs/summary.csv")
    finished = store.finish_task(
        summary="Core analysis completed.",
        evidence_ids=[evidence["evidence_id"]],
    )

    assert evidence["facts"]["row_count"] == 2
    assert finished["ok"] is True
    assert "summary.csv" in finished["final_answer"]
    assert (run_dir / "state" / "final_verification.json").is_file()

    unsupported = store.finish_task(
        summary="Core analysis completed with 999 cells.",
        evidence_ids=[evidence["evidence_id"]],
    )
    assert unsupported["ok"] is False
    assert "numeric claims" in unsupported["error"]

    table.write_text("metric,value\ncells,11\n", encoding="utf-8")
    stale = store.finish_task(summary="Done", evidence_ids=[evidence["evidence_id"]])
    assert stale["ok"] is False
    assert "changed since inspection" in stale["error"]


def test_default_tools_expose_workspace_jobs_and_grounded_finish(tmp_path: Path) -> None:
    config = _config(tmp_path)
    names = {tool.name for tool in build_tools(config, logger=None, run_dir=config.runs_dir / "run")}

    assert {
        "write_run_file",
        "edit_run_file",
        "validate_script",
        "start_job",
        "list_jobs",
        "get_job",
        "poll_job",
        "tail_job",
        "cancel_job",
        "inspect_artifact",
        "finish_task",
    } <= names
