from __future__ import annotations

import json
from pathlib import Path

from langchain_core.messages import AIMessage, HumanMessage, SystemMessage, ToolMessage

import bioagent.runner as runner_module
import bioagent.tools.registry as registry_module
from bioagent.agent import AgentRunResult
from bioagent.agent import BioAgent
from bioagent.agent import messages_with_turn_guidance
from bioagent.agent import messages_for_model
from bioagent.agent import tool_message_text
from bioagent.config import AgentConfig
from bioagent.logging_utils import RunLogger
from bioagent.run_state import load_pending_state, save_pending_state
from bioagent.tools.registry import build_tools


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


def test_request_user_input_tool_returns_structured_interrupt(tmp_path: Path) -> None:
    config = _config(tmp_path)
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    result = tools["request_user_input"].invoke(
        {
            "question": "Which species is this dataset?",
            "reason": "The selected single-cell skill needs species-specific QC genes.",
            "required_fields": ["species"],
        }
    )

    assert result["status"] == "needs_user_input"
    assert result["question"] == "Which species is this dataset?"
    assert result["required_fields"] == ["species"]


def test_pending_state_round_trips_messages(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    messages = [SystemMessage("system"), HumanMessage("task")]

    state_path = save_pending_state(
        config=config,
        run_id="run_test",
        run_dir=run_dir,
        log_path=config.logs_dir / "run_test.log",
        messages=messages,
        question="Need species?",
        metadata={
            "required_fields": ["species"],
            "memory_state": {
                "task_state": {"task": "Analyze cells", "selected_skill": "scrnaseq"},
                "execution_outcome": "failure",
                "execution_tool": "execute_python",
            },
        },
    )
    loaded = load_pending_state(config, "run_test")

    assert state_path.exists()
    assert loaded.run_id == "run_test"
    assert loaded.question == "Need species?"
    assert [message.content for message in loaded.messages] == ["system", "task"]


class _FakeLLM:
    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        return AIMessage(
            content="",
            tool_calls=[
                {
                    "name": "request_user_input",
                    "args": {
                        "question": "Which species should I use?",
                        "reason": "Species-specific QC is required.",
                        "required_fields": ["species"],
                    },
                    "id": "call_1",
                }
            ],
        )


class _PauseAfterResponseLLM:
    def __init__(self, pause_state: dict[str, bool]) -> None:
        self.pause_state = pause_state

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.pause_state["requested"] = True
        return AIMessage(
            content="I will run Python.",
            tool_calls=[
                {
                    "name": "execute_python",
                    "args": {"code": "print('should not run')"},
                    "id": "call_pause",
                }
            ],
        )


class _CaptureMessagesLLM:
    def __init__(self) -> None:
        self.messages = []

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.messages = list(messages)
        return AIMessage(content="Steering received.")


def test_agent_injects_pending_session_messages_before_the_next_model_call(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_steering"
    run_dir.mkdir(parents=True)
    pending = ["Use only chromosome 22 and keep the current GWAS goal."]

    with RunLogger(config.logs_dir, run_id="run_steering") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        llm = _CaptureMessagesLLM()
        agent.llm = llm
        result = agent.run(
            "Run the GWAS analysis",
            max_turns=1,
            take_pending_messages=lambda: [pending.pop(0)] if pending else [],
        )

    contents = [str(message.content) for message in llm.messages]
    assert result.status == "completed"
    assert any("Session steering message" in content and "chromosome 22" in content for content in contents)


def test_agent_pauses_after_model_boundary_and_checkpoints_cancelled_tool_call(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_pause"
    run_dir.mkdir(parents=True)
    pause_state = {"requested": False}
    executed = {"count": 0}
    events: list[dict] = []

    def fake_execute(*args, **kwargs):
        executed["count"] += 1
        return {"ok": True}

    monkeypatch.setattr(registry_module, "execute_python_impl", fake_execute)

    with RunLogger(config.logs_dir, run_id="run_pause") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _PauseAfterResponseLLM(pause_state)
        result = agent.run(
            "Analyze data",
            max_turns=3,
            event_sink=events.append,
            pause_requested=lambda: pause_state["requested"],
        )

    assert result.status == "paused"
    assert result.resume_id == "run_pause"
    assert executed["count"] == 0
    loaded = load_pending_state(config, "run_pause")
    assert loaded.metadata["resume_kind"] == "user_pause"
    assert isinstance(loaded.messages[-1], ToolMessage)
    assert "cancelled" in str(loaded.messages[-1].content).lower()
    cancelled = [event for event in events if event.get("type") == "tool_result"]
    assert cancelled[-1]["status"] == "cancelled"


def test_agent_stops_and_persists_state_when_tool_requests_user_input(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _FakeLLM()
        result = agent.run("Analyze this dataset", max_turns=1)

    assert result.status == "needs_user_input"
    assert result.pending_question == "Which species should I use?"
    assert result.resume_id == "run_test"
    assert Path(result.pending_state_path).exists()

    loaded = load_pending_state(config, "run_test")
    assert loaded.question == "Which species should I use?"


def test_resume_bio_agent_restores_pending_state_and_clears_it_after_completion(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    save_pending_state(
        config=config,
        run_id="run_test",
        run_dir=run_dir,
        log_path=config.logs_dir / "run_test.log",
        messages=[SystemMessage("system"), HumanMessage("task")],
        question="Which species?",
        metadata={
            "required_fields": ["species"],
            "memory_state": {
                "task_state": {"task": "Analyze cells", "selected_skill": "scrnaseq"},
                "execution_outcome": "failure",
                "execution_tool": "execute_python",
            },
        },
    )

    captured: dict = {}

    class FakeBioAgent:
        def __init__(self, *, config, logger, run_dir):
            captured["run_dir"] = run_dir

        def run(self, task, *, max_turns, initial_messages, resume_answer, initial_memory_state):
            captured["task"] = task
            captured["max_turns"] = max_turns
            captured["messages"] = initial_messages
            captured["resume_answer"] = resume_answer
            captured["initial_memory_state"] = initial_memory_state
            return AgentRunResult(
                final_answer="completed",
                run_dir=str(run_dir),
                log_path=str(config.logs_dir / "run_test.log"),
                turns=1,
                status="completed",
            )

    monkeypatch.setattr(runner_module.AgentConfig, "from_env", classmethod(lambda cls: config))
    monkeypatch.setattr(runner_module, "BioAgent", FakeBioAgent)

    result = runner_module.resume_bio_agent("run_test", "mouse", max_turns=3)

    assert result["status"] == "completed"
    assert captured["resume_answer"] == "mouse"
    assert captured["initial_memory_state"]["task_state"]["selected_skill"] == "scrnaseq"
    assert [message.content for message in captured["messages"]] == ["system", "task"]
    assert not (run_dir / "state" / "pending_user_input.json").exists()


def test_resume_bio_agent_keeps_checkpoint_when_run_pauses_again(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    save_pending_state(
        config=config,
        run_id="run_test",
        run_dir=run_dir,
        log_path=config.logs_dir / "run_test.log",
        messages=[SystemMessage("system"), HumanMessage("task")],
        question="Paused",
        metadata={"resume_kind": "user_pause"},
    )

    class FakeBioAgent:
        def __init__(self, *, config, logger, run_dir):
            self.config = config
            self.logger = logger
            self.run_dir = run_dir

        def run(self, task, *, max_turns, initial_messages, resume_answer, pause_requested=None):
            save_pending_state(
                config=self.config,
                run_id="run_test",
                run_dir=self.run_dir,
                log_path=self.logger.path,
                messages=initial_messages,
                question="Paused again",
                metadata={"resume_kind": "user_pause"},
            )
            return AgentRunResult(
                final_answer="Paused again",
                run_dir=str(self.run_dir),
                log_path=str(self.logger.path),
                turns=1,
                status="paused",
                resume_id="run_test",
            )

    monkeypatch.setattr(runner_module.AgentConfig, "from_env", classmethod(lambda cls: config))
    monkeypatch.setattr(runner_module, "BioAgent", FakeBioAgent)

    result = runner_module.resume_bio_agent("run_test", "continue", max_turns=3)

    assert result["status"] == "paused"
    assert load_pending_state(config, "run_test").question == "Paused again"


def test_messages_for_model_compacts_old_tool_outputs(tmp_path: Path) -> None:
    config = _config(tmp_path)
    messages = [SystemMessage("system"), HumanMessage("task")]
    for idx in range(8):
        messages.append(ToolMessage(content=f"tool-{idx} " + ("x" * 2000), tool_call_id=f"call_{idx}"))

    compacted = messages_for_model(messages, config)
    tool_messages = [message for message in compacted if getattr(message, "type", "") == "tool"]

    assert len(compacted) == len(messages)
    assert "compressed" in str(tool_messages[0].content)
    assert len(str(tool_messages[0].content)) < 500
    assert "tool-7" in str(tool_messages[-1].content)
    assert len(str(tool_messages[-1].content)) > 1000


def test_messages_for_model_trims_very_large_recent_tool_output(tmp_path: Path) -> None:
    config = _config(tmp_path)
    messages = [
        SystemMessage("system"),
        HumanMessage("task"),
        ToolMessage(content="START" + ("x" * 5000) + "END", tool_call_id="call_1"),
    ]

    compacted = messages_for_model(messages, config)
    tool_message = [message for message in compacted if getattr(message, "type", "") == "tool"][0]

    assert len(str(tool_message.content)) < 2500
    assert "recent tool output truncated" in str(tool_message.content)
    assert str(tool_message.content).startswith("START")
    assert "END" in str(tool_message.content)


def test_messages_for_model_compacts_execution_code_inside_ai_tool_calls(tmp_path: Path) -> None:
    config = _config(tmp_path)
    code = "BEGIN\n" + ("x = 1\n" * 4000) + "END"
    message = AIMessage(
        content="Run the analysis.",
        tool_calls=[
            {
                "name": "execute_python",
                "args": {"code": code, "env_profile": "py-general-v1", "timeout_s": 1800},
                "id": "call_code",
            }
        ],
    )

    compacted = messages_for_model([SystemMessage("system"), HumanMessage("task"), message], config)
    compacted_code = compacted[-1].tool_calls[0]["args"]["code"]

    assert len(compacted_code) < 1800
    assert compacted_code.startswith("BEGIN")
    assert compacted_code.endswith("END")
    assert "execution code compacted" in compacted_code


def test_messages_for_model_preserves_verified_skill_symbols_when_compacting(tmp_path: Path) -> None:
    config = _config(tmp_path)
    symbol_result = {
        "skill_id": "scrnaseq-scanpy-core-analysis",
        "script_path": "scripts/qc_metrics.py",
        "functions": [
            {
                "name": "get_species_mito_pattern",
                "signature": "def get_species_mito_pattern(species: str) -> str",
            },
            {
                "name": "calculate_qc_metrics",
                "signature": "def calculate_qc_metrics(adata, species: str='human')",
            },
        ],
        "classes": [],
    }
    messages = [
        SystemMessage("system"),
        HumanMessage("task"),
        AIMessage(
            content="Inspect symbols.",
            tool_calls=[
                {
                    "name": "inspect_skill_script_symbols",
                    "args": {"skill_id": symbol_result["skill_id"], "script_path": symbol_result["script_path"]},
                    "id": "call_symbols",
                }
            ],
        ),
        ToolMessage(content=json.dumps(symbol_result), tool_call_id="call_symbols"),
        ToolMessage(content="recent one", tool_call_id="call_recent_1"),
        ToolMessage(content="recent two", tool_call_id="call_recent_2"),
    ]

    compacted = messages_for_model(messages, config)
    symbol_message = [message for message in compacted if getattr(message, "tool_call_id", "") == "call_symbols"][0]

    assert "Verified skill symbols" in str(symbol_message.content)
    assert "scripts/qc_metrics.py" in str(symbol_message.content)
    assert "calculate_qc_metrics" in str(symbol_message.content)
    assert "compute_qc_metrics" not in str(symbol_message.content)
    assert len(str(symbol_message.content)) < 1000


def test_messages_with_turn_guidance_adds_remaining_budget(tmp_path: Path) -> None:
    config = _config(tmp_path)
    messages = [SystemMessage("system"), HumanMessage("task")]

    guided = messages_with_turn_guidance(messages, config, turn=19, max_turns=20)

    assert len(guided) == 3
    assert "decision step 19 of 20" in str(guided[-1].content)
    assert "1 decision step(s) remain" in str(guided[-1].content)
    assert "do not write, edit, restart, cancel, or launch new execution" in str(guided[-1].content)


def test_initial_message_tells_agent_to_use_provided_primary_data_path(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        message = agent._initial_user_message("Analyze data", "/home/luting/projects/mas/mas_2/data/bmmc_b_cell.h5ad", None)

    assert "Primary data path: /home/luting/projects/mas/mas_2/data/bmmc_b_cell.h5ad" in message
    assert "Use this exact primary data path" in message
    assert "do not search for or invent alternate data paths" in message

    relative_message = agent._initial_user_message("Analyze data", r"mas_2\data\bmmc_b_cell.h5ad", None)
    assert "Primary data path inside Docker: /repo/mas_2/data/bmmc_b_cell.h5ad" in relative_message
    assert "Do not call list_files only to locate the primary data path" in relative_message


def test_failed_execution_tool_message_surfaces_error_reason_first() -> None:
    text = tool_message_text(
        "execute_python",
        {
            "ok": False,
            "exit_code": 1,
            "error_reason": "TypeError: find_all_cluster_markers() got an unexpected keyword argument 'output_dir'",
            "stdout": "before failure\n",
            "stderr": "font warning\n" * 1000,
        },
    )

    assert text.index("error_reason") < 80
    assert "unexpected keyword argument 'output_dir'" in text[:300]
    assert len(text) < 2500


def test_failed_poll_tool_message_keeps_persisted_job_logs_and_active_id() -> None:
    text = tool_message_text(
        "poll_job",
        {
            "ok": False,
            "status": "failed",
            "job_id": "job_20260716_120000_abcdef",
            "active_job_id": "job_20260716_120000_abcdef",
            "exit_code": 1,
            "error_reason": "Job ended with status=failed",
            "stdout_tail": "loaded 3512 cells\n",
            "stderr_tail": "TypeError: unexpected keyword argument\n",
        },
    )

    assert "job_20260716_120000_abcdef" in text
    assert "loaded 3512 cells" in text
    assert "TypeError: unexpected keyword argument" in text


def _write_guardrail_skill(config: AgentConfig) -> None:
    skill_dir = config.workflows_root / "guardrail-skill"
    scripts_dir = skill_dir / "scripts"
    scripts_dir.mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text(
        """---
id: guardrail-skill
name: Guardrail Skill
runtime: python
env_profile: py-general-v1
---

# Guardrail Skill
""",
        encoding="utf-8",
    )
    (scripts_dir / "helpers.py").write_text("def helper(x):\n    return x\n", encoding="utf-8")


class _InspectThenExecuteLLM:
    def __init__(self) -> None:
        self.messages = []
        self.calls = 0

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        self.messages = list(messages)
        if self.calls == 1:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "inspect_workflow_skill",
                        "args": {"skill_id": "guardrail-skill"},
                        "id": "call_inspect",
                    }
                ],
            )
        if self.calls == 2:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "execute_python",
                        "args": {"code": "print('should not run')", "env_profile": "py-general-v1"},
                        "id": "call_execute",
                    }
                ],
            )
        return AIMessage(content="stopped")


def test_agent_blocks_execution_after_skill_inspection_until_script_is_read(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    _write_guardrail_skill(config)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    def fail_execute(*args, **kwargs):
        raise AssertionError("execute_python should not run before script inspection")

    monkeypatch.setattr(registry_module, "execute_python_impl", fail_execute)
    capture = _InspectThenExecuteLLM()

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = capture
        result = agent.run("Analyze with guardrail skill", max_turns=3)

    assert result.final_answer == "stopped"
    tool_messages = [message for message in capture.messages if getattr(message, "type", "") == "tool"]
    assert any("workflow_skill_detail_required" in str(message.content) for message in tool_messages)
    assert not (run_dir / "scripts" / "code.py").exists()


class _InspectReadThenExecuteLLM:
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
                        "name": "inspect_workflow_skill",
                        "args": {"skill_id": "guardrail-skill"},
                        "id": "call_inspect",
                    }
                ],
            )
        if self.calls == 2:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "inspect_skill_script_symbols",
                        "args": {"skill_id": "guardrail-skill", "script_path": "scripts/helpers.py"},
                        "id": "call_symbols",
                    }
                ],
            )
        if self.calls == 3:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "execute_python",
                        "args": {"code": "print('allowed')", "env_profile": "py-general-v1"},
                        "id": "call_execute",
                    }
                ],
            )
        return AIMessage(content="finished")


class _InspectThenDataPreviewLLM:
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
                        "name": "inspect_workflow_skill",
                        "args": {"skill_id": "guardrail-skill"},
                        "id": "call_inspect",
                    }
                ],
            )
        if self.calls == 2:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "execute_python",
                        "args": {
                            "code": (
                                "import scanpy as sc\n"
                                "adata = sc.read_h5ad('/repo/mas_2/data/example.h5ad')\n"
                                "print(adata.shape)\n"
                                "if 'leiden' in adata.obs.columns:\n"
                                "    print(adata.obs['leiden'].value_counts())"
                            ),
                            "env_profile": "py-scverse-v1",
                        },
                        "id": "call_preview",
                    }
                ],
            )
        return AIMessage(content="finished")


class _SuccessfulExecutionNearBudgetLLM:
    def __init__(self) -> None:
        self.calls = 0

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
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


class _RepeatedContextLookupLLM:
    def __init__(self) -> None:
        self.calls = 0
        self.messages = []

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        self.messages = list(messages)
        if self.calls == 1:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "inspect_workflow_skill",
                        "args": {"skill_id": "guardrail-skill"},
                        "id": "call_inspect",
                    }
                ],
            )
        if self.calls <= 10:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "list_files",
                        "args": {"path": "."},
                        "id": f"call_list_{self.calls}",
                    }
                ],
            )
        return AIMessage(content="stopped")


class _DistinctContextLookupLLM:
    def __init__(self) -> None:
        self.calls = 0
        self.messages = []

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        self.messages = list(messages)
        if self.calls == 1:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "inspect_workflow_skill",
                        "args": {"skill_id": "guardrail-skill"},
                        "id": "call_inspect",
                    }
                ],
            )
        if self.calls <= 10:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "list_files",
                        "args": {"path": f"directory-{self.calls}"},
                        "id": f"call_list_{self.calls}",
                    }
                ],
            )
        return AIMessage(content="stopped")


class _NonConsecutiveContextLookupLLM:
    def __init__(self) -> None:
        self.calls = 0

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        paths = {1: ".", 2: "other", 3: "."}
        if self.calls <= 3:
            return AIMessage(
                content="",
                tool_calls=[
                    {
                        "name": "list_files",
                        "args": {"path": paths[self.calls]},
                        "id": f"call_list_{self.calls}",
                    }
                ],
            )
        return AIMessage(content="stopped")


class _FinalTurnToolLLM:
    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        return AIMessage(
            content="",
            tool_calls=[
                {
                    "name": "execute_python",
                    "args": {"code": "print('should not run')"},
                    "id": "call_execute",
                }
            ],
        )


class _FinalTurnContentAndToolLLM:
    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        return AIMessage(
            content="final summary",
            tool_calls=[
                {
                    "name": "manage_memory",
                    "args": {"action": "create", "content": "remember this"},
                    "id": "call_memory",
                }
            ],
        )


class _EmptyThenFinalLLM:
    def __init__(self) -> None:
        self.calls = 0

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        if self.calls == 1:
            return AIMessage(content="")
        return AIMessage(content="final answer")


class _ContinuationThenFinalLLM:
    def __init__(self) -> None:
        self.calls = 0

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        if self.calls == 1:
            return AIMessage(content="我发现问题了，让我检查数据格式并修复。")
        return AIMessage(content="final answer")


class _WaitInspectAndFinishLLM:
    def __init__(self) -> None:
        self.calls = 0

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        if self.calls <= 2:
            return AIMessage(
                content="Waiting for the existing job.",
                tool_calls=[
                    {
                        "name": "poll_job",
                        "args": {"job_id": "job_20260716_120000_abcdef", "wait_s": 300},
                        "id": f"poll_{self.calls}",
                    }
                ],
            )
        if self.calls == 3:
            return AIMessage(
                content="Verifying the completed output.",
                tool_calls=[
                    {
                        "name": "inspect_artifact",
                        "args": {"path": "/work/outputs/summary.json"},
                        "id": "inspect_1",
                    }
                ],
            )
        tool_message = next(message for message in reversed(messages) if getattr(message, "type", "") == "tool")
        evidence = json.loads(str(tool_message.content))["evidence_id"]
        return AIMessage(
            content="Finishing with verified evidence.",
            tool_calls=[
                {
                    "name": "finish_task",
                    "args": {"summary": "Analysis completed.", "evidence_ids": [evidence]},
                    "id": "finish_1",
                }
            ],
        )


def test_agent_allows_execution_after_skill_script_inspection(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    _write_guardrail_skill(config)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    executed = {"count": 0}

    def fake_execute(*args, **kwargs):
        executed["count"] += 1
        return {"ok": True, "stdout": "allowed", "stderr": ""}

    monkeypatch.setattr(registry_module, "execute_python_impl", fake_execute)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _InspectReadThenExecuteLLM()
        result = agent.run("Analyze with guardrail skill", max_turns=5)

    assert result.final_answer == "finished"
    assert executed["count"] == 1


def test_agent_finalizes_successful_execution_near_turn_limit(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    (run_dir / "outputs").mkdir(parents=True)
    (run_dir / "outputs" / "summary.csv").write_text("metric,value\ncells,10\n", encoding="utf-8")

    def fake_execute(*args, **kwargs):
        return {
            "ok": True,
            "stdout": "analysis complete",
            "stderr": "",
            "env_profile": "py-general-v1",
            "run_dir": str(run_dir),
            "script_path": str(run_dir / "scripts" / "code.py"),
            "exit_code": 0,
        }

    fake_llm = _SuccessfulExecutionNearBudgetLLM()
    monkeypatch.setattr(registry_module, "execute_python_impl", fake_execute)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = fake_llm
        result = agent.run("Analyze data", max_turns=1)

    assert fake_llm.calls == 1
    assert "代码执行已成功" in result.final_answer
    assert "summary.csv" in result.final_answer
    assert result.turns == 1


def test_agent_blocks_only_consecutive_duplicate_context_lookup(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    _write_guardrail_skill(config)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    list_calls = {"count": 0}

    def fake_list_files(config, run_dir, *, path, recursive=False, max_entries=200):
        list_calls["count"] += 1
        return "listing"

    capture = _RepeatedContextLookupLLM()
    monkeypatch.setattr(registry_module, "list_files_impl", fake_list_files)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = capture
        result = agent.run("Analyze with too much context lookup", max_turns=11)

    assert result.final_answer == "stopped"
    assert list_calls["count"] == 1
    tool_messages = [message for message in capture.messages if getattr(message, "type", "") == "tool"]
    assert any("duplicate_context_lookup" in str(message.content) for message in tool_messages)


def test_agent_allows_more_than_eight_distinct_context_lookups(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    _write_guardrail_skill(config)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    paths: list[str] = []

    def fake_list_files(config, run_dir, *, path, recursive=False, max_entries=200):
        paths.append(path)
        return f"listing for {path}"

    capture = _DistinctContextLookupLLM()
    monkeypatch.setattr(registry_module, "list_files_impl", fake_list_files)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = capture
        result = agent.run("Analyze with distinct context lookups", max_turns=11)

    assert result.final_answer == "stopped"
    assert len(paths) == 9
    assert len(set(paths)) == 9
    tool_messages = [message for message in capture.messages if getattr(message, "type", "") == "tool"]
    assert not any("context_lookup_budget_exhausted" in str(message.content) for message in tool_messages)
    assert not any("duplicate_context_lookup" in str(message.content) for message in tool_messages)


def test_agent_allows_same_context_lookup_after_different_call(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    paths: list[str] = []

    def fake_list_files(config, run_dir, *, path, recursive=False, max_entries=200):
        paths.append(path)
        return f"listing for {path}"

    monkeypatch.setattr(registry_module, "list_files_impl", fake_list_files)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _NonConsecutiveContextLookupLLM()
        result = agent.run("Inspect, change scope, then inspect again", max_turns=4)

    assert result.final_answer == "stopped"
    assert paths == [".", "other", "."]


def test_agent_skips_tool_calls_on_final_turn(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    def fail_execute(*args, **kwargs):
        raise AssertionError("final turn tool calls should not execute")

    monkeypatch.setattr(registry_module, "execute_python_impl", fail_execute)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _FinalTurnToolLLM()
        result = agent.run("Analyze data", max_turns=1)

    assert "Reached max_turns=1" in result.final_answer
    assert "已跳过执行" in result.final_answer


def test_agent_keeps_final_turn_content_when_skipping_tool_calls(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _FinalTurnContentAndToolLLM()
        result = agent.run("Analyze data", max_turns=1)

    assert result.final_answer == "final summary"


def test_agent_does_not_finish_on_empty_model_response(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    fake_llm = _EmptyThenFinalLLM()

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = fake_llm
        result = agent.run("Analyze data", max_turns=3)

    assert fake_llm.calls == 2
    assert result.final_answer == "final answer"


def test_agent_does_not_finish_on_continuation_without_tool(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    fake_llm = _ContinuationThenFinalLLM()

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = fake_llm
        result = agent.run("Analyze data", max_turns=3)

    assert fake_llm.calls == 2
    assert result.final_answer == "final answer"


def test_waiting_and_artifact_verification_do_not_consume_decision_budget(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    outputs = run_dir / "outputs"
    outputs.mkdir(parents=True)
    (outputs / "summary.json").write_text('{"status": "complete"}\n', encoding="utf-8")
    poll_count = {"value": 0}

    def fake_poll_job(self, job_id="", *, wait_s=0):
        poll_count["value"] += 1
        status = "running" if poll_count["value"] == 1 else "completed"
        return {
            "ok": True,
            "job_id": job_id,
            "status": status,
            "exit_code": 0 if status == "completed" else None,
            "artifacts": [{"path": str(outputs / "summary.json"), "size_bytes": 23}],
        }

    monkeypatch.setattr("bioagent.tools.jobs.DockerJobManager.poll_job", fake_poll_job)
    fake_llm = _WaitInspectAndFinishLLM()

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = fake_llm
        result = agent.run("Wait for and verify the analysis", max_turns=1)

    assert fake_llm.calls == 4
    assert poll_count["value"] == 2
    assert result.status == "completed"
    assert result.turns == 1
    assert "Verified artifacts" in result.final_answer


def test_agent_allows_lightweight_h5ad_preview_before_skill_script_inspection(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    _write_guardrail_skill(config)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    executed = {"count": 0}

    def fake_execute(*args, **kwargs):
        executed["count"] += 1
        return {"ok": True, "stdout": "(10, 20)", "stderr": ""}

    monkeypatch.setattr(registry_module, "execute_python_impl", fake_execute)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _InspectThenDataPreviewLLM()
        result = agent.run("Preview data with guardrail skill", max_turns=4)

    assert result.final_answer == "finished"
    assert executed["count"] == 1
