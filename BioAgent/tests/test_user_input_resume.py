from __future__ import annotations

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
        metadata={"required_fields": ["species"]},
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
        metadata={"required_fields": ["species"]},
    )

    captured: dict = {}

    class FakeBioAgent:
        def __init__(self, *, config, logger, run_dir):
            captured["run_dir"] = run_dir

        def run(self, task, *, max_turns, initial_messages, resume_answer):
            captured["task"] = task
            captured["max_turns"] = max_turns
            captured["messages"] = initial_messages
            captured["resume_answer"] = resume_answer
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
    assert [message.content for message in captured["messages"]] == ["system", "task"]
    assert not (run_dir / "state" / "pending_user_input.json").exists()


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


def test_messages_with_turn_guidance_adds_remaining_budget(tmp_path: Path) -> None:
    config = _config(tmp_path)
    messages = [SystemMessage("system"), HumanMessage("task")]

    guided = messages_with_turn_guidance(messages, config, turn=19, max_turns=20)

    assert len(guided) == 3
    assert "turn 19 of 20" in str(guided[-1].content)
    assert "1 turn(s) remain" in str(guided[-1].content)
    assert "do not start a new analysis step" in str(guided[-1].content)


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
        result = agent.run("Analyze data", max_turns=2)

    assert fake_llm.calls == 1
    assert "代码执行已成功" in result.final_answer
    assert "summary.csv" in result.final_answer
    assert result.turns == 1


def test_agent_limits_context_lookup_after_skill_inspection(tmp_path: Path, monkeypatch) -> None:
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
    assert list_calls["count"] == 8
    tool_messages = [message for message in capture.messages if getattr(message, "type", "") == "tool"]
    assert any("context_lookup_budget_exhausted" in str(message.content) for message in tool_messages)


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
