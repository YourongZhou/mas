from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path

from langchain_core.messages import AIMessage

import bioagent.agent as agent_module
import bioagent.tools.registry as registry_module
from bioagent.agent import BioAgent
from bioagent.config import AgentConfig
from bioagent.logging_utils import RunLogger
from bioagent.memory import MemoryHarness, build_memory_harness
from bioagent.memory.extraction import extract_verified_episode
from bioagent.memory.long_term import JsonEpisodeStore
from bioagent.memory.retrieval import EpisodeRetriever
from bioagent.memory.schemas import MemoryEpisode, TaskState
from bioagent.memory.short_term import ShortTermMemory
from bioagent.run_state import load_pending_state


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


def _episode(**overrides: str) -> MemoryEpisode:
    values = {
        "task_signature": "single cell scanpy core analysis",
        "data_signature": "format=.h5ad",
        "skill_id": "scrnaseq-scanpy-core-analysis",
        "runtime_environment": "py-scverse-v1",
        "outcome": "success",
        "verified_root_cause": "",
        "verified_fix": "Execution completed successfully in py-scverse-v1.",
        "reusable_lesson": "Use the Scanpy skill with py-scverse-v1 for matching h5ad tasks.",
        "source_run_id": "run_old",
        "timestamp": "2026-07-16T10:00:00+00:00",
    }
    values.update(overrides)
    return MemoryEpisode(**values)


def test_short_term_memory_updates_one_structured_task_state() -> None:
    memory = ShortTermMemory.start(
        task="Run core single-cell analysis",
        data_path="/data/input.h5ad",
        result_dir="/results",
    )
    original_state = memory.state

    memory.observe_tool_call("inspect_workflow_skill", {"skill_id": "scrnaseq-scanpy-core-analysis"})
    memory.observe_tool_call("execute_python", {"env_profile": "py-scverse-v1"})
    memory.observe_tool_result(
        "execute_python",
        {"env_profile": "py-scverse-v1"},
        {"ok": False, "error_reason": "KeyError: pct_mito"},
    )
    memory.observe_tool_result(
        "execute_python",
        {"env_profile": "py-scverse-v1"},
        {"ok": True, "artifacts": ["/work/outputs/processed.h5ad"]},
    )

    assert memory.state is original_state
    assert memory.state.selected_skill == "scrnaseq-scanpy-core-analysis"
    assert memory.state.runtime_environment == "py-scverse-v1"
    assert memory.state.active_errors == []
    assert memory.state.resolved_errors == ["KeyError: pct_mito"]
    assert memory.state.generated_artifacts == ["/work/outputs/processed.h5ad"]
    assert memory.state.next_action == "Review verified outputs and report results."
    summary = memory.compact_summary()
    assert "TaskState" in summary
    assert "selected_skill: scrnaseq-scanpy-core-analysis" in summary
    assert len(summary) < 2200


def test_short_term_memory_overwrites_resume_stage_with_latest_user_input() -> None:
    memory = ShortTermMemory.start(task="Run core analysis", data_path="input.h5ad")
    memory.update(blockers=["Need user decision"], current_stage="waiting_for_user")

    memory.apply_user_input("Skip doublet detection and continue")

    assert memory.state.current_goal == "Skip doublet detection and continue"
    assert memory.state.confirmed_inputs[-1] == "Latest user input: Skip doublet detection and continue"
    assert memory.state.blockers == []
    assert memory.state.current_stage == "resuming"


def test_short_term_memory_persists_harness_issued_active_job() -> None:
    memory = ShortTermMemory.start(task="Run analysis")

    memory.observe_tool_result(
        "start_job",
        {"env_profile": "py-scverse-v1"},
        {"ok": True, "status": "running", "job_id": "job_20260716_120000_abcdef"},
    )
    memory.observe_tool_result(
        "poll_job",
        {"job_id": "job_17"},
        {
            "ok": False,
            "status": "invalid_job_id",
            "active_job_id": "job_20260716_120000_abcdef",
            "error": "invalid job id",
        },
    )

    assert memory.state.active_job_id == "job_20260716_120000_abcdef"
    assert memory.state.job_status == "invalid_job_id"
    assert "active_job_id: job_20260716_120000_abcdef" in memory.compact_summary()


def test_failed_task_keeps_repair_state_for_session_continuation(tmp_path: Path) -> None:
    harness = MemoryHarness(
        enabled=True,
        namespace=("bioagent", "default"),
        store=JsonEpisodeStore(tmp_path / "episodes.json"),
    )
    harness.start_task("Run analysis", data_path="input.h5ad")
    harness.observe_tool_result(
        "poll_job",
        {"job_id": "job_failed"},
        {
            "ok": False,
            "status": "failed",
            "job_id": "job_failed",
            "error_reason": "AttributeError: plotting failed",
        },
    )

    harness.finish_task(source_run_id="run_failed", status="failed")

    state = harness.short_term.state
    assert state.current_stage == "repairing_error"
    assert state.active_errors == ["AttributeError: plotting failed"]
    assert state.next_action == "Use the exact error and verified signatures to make the smallest repair."


def test_json_episode_store_persists_only_episode_schema(tmp_path: Path) -> None:
    path = tmp_path / "memory" / "episodes.json"
    store = JsonEpisodeStore(path)

    store.add(_episode())

    raw = json.loads(path.read_text(encoding="utf-8"))
    assert len(raw) == 1
    assert set(raw[0]) == set(MemoryEpisode.field_names())
    serialized = json.dumps(raw[0])
    assert "stdout" not in serialized
    assert "generated code" not in serialized
    assert "/tmp/" not in serialized
    assert store.list_all() == [_episode()]


def test_extraction_is_conservative_and_sanitizes_temporary_paths() -> None:
    state = TaskState(
        task="Analyze /tmp/private/input.h5ad with Scanpy",
        current_goal="Complete core analysis",
        confirmed_inputs=["Primary data: /tmp/private/input.h5ad"],
        expected_outputs=["processed h5ad"],
        selected_skill="scrnaseq-scanpy-core-analysis",
        runtime_environment="py-scverse-v1",
        current_stage="completed",
        active_errors=[],
        resolved_errors=["FileNotFoundError at /tmp/run/scripts/code.py"],
        generated_artifacts=["/tmp/run/outputs/result.h5ad"],
        blockers=[],
        next_action="Report results",
    )

    assert extract_verified_episode(state, source_run_id="run_none", execution_outcome=None) is None
    episode = extract_verified_episode(state, source_run_id="run_ok", execution_outcome="success")

    assert episode is not None
    assert episode.outcome == "success"
    assert episode.data_signature == "format=.h5ad"
    assert episode.task_signature.startswith("terms=")
    assert "Analyze" not in episode.task_signature
    assert "/tmp/" not in json.dumps(episode.to_dict())
    assert "FileNotFoundError" in episode.verified_root_cause


def test_retrieval_prefers_relevant_verified_episode_and_states_override_rule(tmp_path: Path) -> None:
    store = JsonEpisodeStore(tmp_path / "episodes.json")
    store.add(_episode())
    store.add(
        _episode(
            task_signature="gwas plink association analysis",
            data_signature="format=.bed",
            skill_id="gwas-plink",
            runtime_environment="py-general-v1",
            source_run_id="run_gwas",
        )
    )
    retriever = EpisodeRetriever(store)

    results = retriever.retrieve("Scanpy single-cell analysis", "format=.h5ad", limit=2)
    text = retriever.format_for_prompt(results)

    assert results[0].skill_id == "scrnaseq-scanpy-core-analysis"
    assert "Prior verified experience" in text
    assert "Current observations always override historical experience." in text
    assert "stdout" not in text


class _TwoTurnCaptureLLM:
    def __init__(self) -> None:
        self.calls: list[list] = []

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls.append(messages)
        if len(self.calls) == 1:
            return AIMessage(
                content="",
                tool_calls=[{"name": "list_files", "args": {"path": "."}, "id": "call_list"}],
            )
        return AIMessage(content="Inspection complete")


def test_agent_injects_fresh_task_state_and_prior_experience_into_every_llm_call(
    tmp_path: Path, monkeypatch
) -> None:
    config = _config(tmp_path)
    store = JsonEpisodeStore(tmp_path / "episodes.json")
    store.add(_episode())
    harness = MemoryHarness(enabled=True, namespace=("bioagent", "default"), store=store)
    monkeypatch.setattr(agent_module, "build_memory_harness", lambda config: harness)
    monkeypatch.setattr(registry_module, "list_files_impl", lambda *args, **kwargs: "input.h5ad")
    capture = _TwoTurnCaptureLLM()
    events: list[dict] = []
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = capture
        result = agent.run(
            "Run Scanpy single-cell analysis",
            data_path="input.h5ad",
            max_turns=2,
            event_sink=events.append,
        )

    assert result.final_answer == "Inspection complete"
    assert len(capture.calls) == 2
    first_context = "\n".join(str(message.content) for message in capture.calls[0])
    second_context = "\n".join(str(message.content) for message in capture.calls[1])
    assert "Prior verified experience" in first_context
    assert "Current observations always override historical experience." in first_context
    assert first_context.count("TaskState") == 1
    assert second_context.count("TaskState") == 1
    assert "current_stage: inspecting_files" in second_context
    assert "manage_memory" not in agent.tool_map
    assert "search_memory" not in agent.tool_map
    memory_events = [event["memory"] for event in events if event["type"] == "memory_state"]
    assert memory_events[0]["taskState"]["current_stage"] == "understanding"
    assert any(event["taskState"]["current_stage"] == "inspecting_files" for event in memory_events)
    assert memory_events[0]["priorEpisodes"][0]["skill_id"] == "scrnaseq-scanpy-core-analysis"


def test_paused_run_checkpoints_task_state(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    monkeypatch.setattr(agent_module, "build_memory_harness", lambda config: build_memory_harness(config))
    run_dir = config.runs_dir / "run_pause"
    run_dir.mkdir(parents=True)

    with RunLogger(config.logs_dir, run_id="run_pause") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        result = agent.run(
            "Analyze single-cell data",
            data_path="input.h5ad",
            pause_requested=lambda: True,
        )

    pending = load_pending_state(config, "run_pause")
    task_state = pending.metadata["memory_state"]["task_state"]
    assert result.status == "paused"
    assert task_state["task"] == "Analyze single-cell data"
    assert task_state["confirmed_inputs"] == ["Primary data: input.h5ad"]


class _ExecuteThenFinishLLM:
    def __init__(self) -> None:
        self.calls = 0

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        if self.calls == 1:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "execute_python",
                        "args": {"code": "print('ok')", "env_profile": "py-general-v1"},
                        "id": "call_execute",
                    }
                ],
            )
        return AIMessage(content="Verified execution complete")


def test_agent_writes_verified_episode_only_after_terminal_execution(tmp_path: Path, monkeypatch) -> None:
    config = replace(_config(tmp_path), memory_namespace="test", memory_user_id="user")
    run_dir = config.runs_dir / "run_verified"
    run_dir.mkdir(parents=True)
    monkeypatch.setattr(
        registry_module,
        "execute_python_impl",
        lambda *args, **kwargs: {
            "ok": True,
            "env_profile": "py-general-v1",
            "artifacts": ["/work/outputs/summary.json"],
            "stdout": "ok",
        },
    )

    with RunLogger(config.logs_dir, run_id="run_verified") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _ExecuteThenFinishLLM()
        result = agent.run("Summarize tabular data", data_path="table.csv", max_turns=2)

    episodes = JsonEpisodeStore(config.project_root / "memory" / "test-user-episodes.json").list_all()
    assert result.status == "completed"
    assert len(episodes) == 1
    assert episodes[0].source_run_id == "run_verified"
    assert episodes[0].outcome == "success"
    assert episodes[0].data_signature == "format=.csv"
