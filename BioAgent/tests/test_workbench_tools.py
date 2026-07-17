from __future__ import annotations

import json
import subprocess
from dataclasses import replace
from pathlib import Path

import pytest

from bioagent.config import AgentConfig
from bioagent.tools.artifacts import ArtifactStore
from bioagent.tools.jobs import DockerJobManager, attach_job_provenance
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
    attach_job_provenance(
        [run_dir],
        job_id=started["job_id"],
        tool_call_id="toolcall_start_job",
        event_id="event_start_job",
    )
    output = run_dir / "outputs" / "summary.json"
    output.parent.mkdir(parents=True)
    output.write_text('{"status": "complete"}\n', encoding="utf-8")
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
    assert polled["script_sha256"]
    assert polled["artifacts"][0]["artifact_id"].startswith("artifact_")
    assert polled["artifacts"][0]["created_by_job_id"] == started["job_id"]
    assert polled["artifacts"][0]["created_by_tool_call_id"] == "toolcall_start_job"
    assert polled["artifacts"][0]["created_by_event_id"] == "event_start_job"
    assert polled["execution_manifest_path"].endswith("state/execution_manifest.json")
    assert (run_dir / "jobs" / started["job_id"] / "job.json").is_file()
    evidence = ArtifactStore(run_dir, config=config).inspect(str(output))
    assert evidence["provenance"]["created_by_job_id"] == started["job_id"]
    assert evidence["provenance"]["script_sha256"] == polled["script_sha256"]
    assert any("/work/scripts/analysis.py" in part for part in calls[1])


def test_stale_concurrent_poll_cannot_overwrite_completed_job_as_cancelled(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_poll_race"
    job_id = "job_20260717_130000_abcdef"
    job_dir = run_dir / "jobs" / job_id
    job_dir.mkdir(parents=True)
    metadata = {
        "job_id": job_id,
        "run_id": run_dir.name,
        "run_dir": str(run_dir),
        "container_id": "container-finished",
        "status": "running",
        "callback_state": "pending",
    }
    (job_dir / "job.json").write_text(json.dumps(metadata), encoding="utf-8")
    (run_dir / "outputs").mkdir()
    (run_dir / "outputs" / "summary.json").write_text('{"ok": true}', encoding="utf-8")
    inspect_calls = 0

    def fake_run(command, **kwargs):
        nonlocal inspect_calls
        if command[:2] == ["docker", "inspect"]:
            inspect_calls += 1
            if inspect_calls == 1:
                state = {"Status": "exited", "Running": False, "ExitCode": 0, "Error": ""}
                return subprocess.CompletedProcess(command, 0, json.dumps(state), "")
            return subprocess.CompletedProcess(command, 1, "", "Error: No such object")
        if command[:2] == ["docker", "logs"]:
            return subprocess.CompletedProcess(command, 0, "analysis complete\n", "")
        if command[:3] == ["docker", "rm", "-f"]:
            return subprocess.CompletedProcess(command, 0, "container-finished\n", "")
        raise AssertionError(command)

    monkeypatch.setattr("bioagent.tools.jobs.subprocess.run", fake_run)
    manager = DockerJobManager(config, logger=None, run_dir=run_dir)

    first = manager._poll_once(job_dir, dict(metadata))
    second = manager._poll_once(job_dir, dict(metadata))

    assert first["status"] == "completed"
    assert second["status"] == "completed"
    assert inspect_calls == 1
    assert json.loads((job_dir / "job.json").read_text(encoding="utf-8"))["status"] == "completed"


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


def test_cancel_job_verifies_that_the_container_is_gone(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_cancel_job"
    job_id = "job_20260717_130000_abcdef"
    job_dir = run_dir / "jobs" / job_id
    job_dir.mkdir(parents=True)
    (job_dir / "job.json").write_text(
        json.dumps(
            {
                "job_id": job_id,
                "run_id": run_dir.name,
                "run_dir": str(run_dir),
                "container_id": "container-cancel",
                "status": "running",
                "callback_state": "pending",
            }
        ),
        encoding="utf-8",
    )

    def fake_run(command, **kwargs):
        if command[:3] == ["docker", "rm", "-f"]:
            return subprocess.CompletedProcess(command, 0, "container-cancel\n", "")
        if command[:2] == ["docker", "inspect"]:
            return subprocess.CompletedProcess(command, 1, "", "No such container")
        raise AssertionError(command)

    monkeypatch.setattr("bioagent.tools.jobs.subprocess.run", fake_run)
    result = DockerJobManager(config, logger=None, run_dir=run_dir).cancel_job(job_id)

    assert result["ok"] is True
    assert result["status"] == "cancelled"
    assert result["verified_exit"] is True
    persisted = json.loads((job_dir / "job.json").read_text(encoding="utf-8"))
    assert persisted["cancel_verified"] is True


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
        "propose_plan",
        "update_plan",
    } <= names


def test_plan_tools_persist_semantic_steps_and_block_for_approval(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_plan"
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=run_dir)}

    proposed = tools["propose_plan"].invoke(
        {
            "goal": "Complete a reproducible single-cell core analysis",
            "steps": [
                {
                    "title": "Inspect inputs and skill",
                    "success_criteria": "The dataset and exact skill functions are verified.",
                },
                {
                    "title": "Run and verify analysis",
                    "success_criteria": "Core artifacts pass evidence inspection.",
                },
            ],
            "assumptions": ["Use Scanpy unless the input requires another runtime."],
            "requires_approval": True,
        }
    )

    assert proposed["status"] == "needs_user_input"
    assert proposed["interaction_type"] == "plan_approval"
    assert proposed["question_id"].startswith("question_")
    assert proposed["plan"]["status"] == "awaiting_approval"
    first_step_id = proposed["plan"]["steps"][0]["id"]

    updated = tools["update_plan"].invoke(
        {
            "step_id": first_step_id,
            "step_status": "in_progress",
            "plan_status": "active",
            "note": "Input path confirmed.",
        }
    )

    assert updated["ok"] is True
    assert updated["plan"]["status"] == "active"
    assert updated["plan"]["steps"][0]["status"] == "in_progress"
    assert updated["plan"]["steps"][0]["note"] == "Input path confirmed."
    persisted = json.loads((run_dir / "state" / "plan.json").read_text(encoding="utf-8"))
    assert persisted == updated["plan"]
