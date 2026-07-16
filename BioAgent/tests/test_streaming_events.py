from __future__ import annotations

from pathlib import Path

from langchain_core.messages import AIMessage

import bioagent.runner as runner_module
import bioagent.tools.registry as registry_module
from bioagent import run_bio_agent_stream
from bioagent.agent import AgentRunResult, BioAgent, _tool_succeeded
from bioagent.config import AgentConfig
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


class _ToolThenFinalLLM:
    def __init__(self) -> None:
        self.calls = 0

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.calls += 1
        if self.calls == 1:
            return AIMessage(
                content="I will run Python.",
                tool_calls=[
                    {
                        "name": "execute_python",
                        "args": {"code": "print('hello')", "env_profile": "py-general-v1"},
                        "id": "call_python",
                    }
                ],
            )
        return AIMessage(content="All done.")


def test_agent_run_emits_standard_conversation_events(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    events: list[dict] = []

    def fake_execute(*args, **kwargs):
        return {"ok": True, "stdout": "hello\n", "stderr": "", "exit_code": 0}

    monkeypatch.setattr(registry_module, "execute_python_impl", fake_execute)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = _ToolThenFinalLLM()
        result = agent.run("Analyze data", max_turns=3, event_sink=events.append)

    event_types = [event["type"] for event in events]
    assert event_types[0] == "run_start"
    assert "assistant_message" in event_types
    assert "tool_call" in event_types
    assert "tool_result" in event_types
    assert event_types[-2:] == ["final", "run_end"]
    assert result.final_answer == "All done."
    assert events[-2]["content"] == "All done."
    assert events[-1]["result"]["final_answer"] == "All done."


def test_run_bio_agent_stream_yields_events_from_notebook_entrypoint(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)

    class FakeBioAgent:
        def __init__(self, *, config, logger, run_dir):
            self.logger = logger
            self.run_dir = run_dir

        def run(self, task, *, data_path=None, result_dir=None, max_turns=20, event_sink=None):
            assert event_sink is not None
            event_sink({"type": "assistant_message", "run_id": self.logger.run_id, "content": "streaming"})
            return AgentRunResult(
                final_answer="completed",
                run_dir=str(self.run_dir),
                log_path=str(self.logger.path),
                turns=1,
                status="completed",
            )

    monkeypatch.setattr(runner_module.AgentConfig, "from_env", classmethod(lambda cls: config))
    monkeypatch.setattr(runner_module, "BioAgent", FakeBioAgent)

    events = list(runner_module.run_bio_agent_stream("Analyze data", max_turns=1))

    assert [event["type"] for event in events] == ["assistant_message", "final", "run_end"]
    assert events[0]["content"] == "streaming"
    assert events[-1]["result"]["final_answer"] == "completed"


def test_streaming_entrypoint_is_exported_from_package() -> None:
    assert run_bio_agent_stream is runner_module.run_bio_agent_stream


def test_tool_success_respects_explicit_false_ok_flag() -> None:
    assert _tool_succeeded({"ok": False, "error_reason": "Timed out"}) is False
    assert _tool_succeeded({"ok": True, "stdout": "done"}) is True
