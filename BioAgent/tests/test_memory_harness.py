from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from langchain_core.tools import tool

import bioagent.agent as agent_module
from bioagent.agent import BioAgent
from bioagent.config import AgentConfig
from bioagent.memory import MemoryHarness, build_memory_harness
from bioagent.logging_utils import RunLogger
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


def test_memory_harness_reports_disabled_when_langmem_unavailable(tmp_path: Path) -> None:
    config = _config(tmp_path)

    harness = build_memory_harness(
        config,
        import_langmem_tools=lambda: (_ for _ in ()).throw(ModuleNotFoundError("langmem")),
    )

    assert not harness.enabled
    assert harness.tools == []
    assert "langmem" in harness.disabled_reason.lower()


def test_memory_harness_exposes_langmem_tools_when_available(tmp_path: Path) -> None:
    config = replace(_config(tmp_path), memory_enabled=True, memory_user_id="luting")

    @tool("manage_memory")
    def manage_memory(content: str) -> str:
        """Manage memory."""
        return content

    @tool("search_memory")
    def search_memory(query: str) -> str:
        """Search memory."""
        return query

    def fake_import():
        return lambda **_: manage_memory, lambda **_: search_memory

    harness = build_memory_harness(config, import_langmem_tools=fake_import)

    assert harness.enabled
    assert [tool.name for tool in harness.tools] == ["manage_memory", "search_memory"]
    assert harness.namespace == ("bioagent", "luting")


def test_build_tools_includes_memory_tools(tmp_path: Path) -> None:
    config = _config(tmp_path)

    @tool("search_memory")
    def search_memory(query: str) -> str:
        """Search memory."""
        return query

    harness = MemoryHarness(
        enabled=True,
        namespace=("bioagent", "default"),
        tools=[search_memory],
        disabled_reason="",
    )

    tools = build_tools(config, logger=None, run_dir=tmp_path / "run", memory_harness=harness)

    assert "search_memory" in {tool.name for tool in tools}


def test_memory_harness_recall_formats_search_results() -> None:
    @tool("search_memory")
    def search_memory(query: str, limit: int = 5) -> list[dict]:
        """Search memory."""
        return [
            {
                "value": {
                    "content": "Use group, not cluster, for Scanpy rank_genes_groups output."
                }
            }
        ]

    harness = MemoryHarness(
        enabled=True,
        namespace=("bioagent", "default"),
        tools=[search_memory],
        disabled_reason="",
    )

    text = harness.recall_text("single-cell marker dataframe")

    assert "Relevant BioAgent memories" in text
    assert "Use group, not cluster" in text


class _CaptureLLM:
    def __init__(self) -> None:
        self.messages = []

    def bind_tools(self, tools):
        return self

    def invoke(self, messages):
        self.messages = messages
        from langchain_core.messages import AIMessage

        return AIMessage(content="done")


def test_agent_injects_recalled_memory_into_initial_context(tmp_path: Path, monkeypatch) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    harness = MemoryHarness(
        enabled=True,
        namespace=("bioagent", "default"),
        tools=[],
        disabled_reason="",
    )
    monkeypatch.setattr(harness, "recall_text", lambda query: "Relevant BioAgent memories:\n- prefer py-scverse-v1")
    monkeypatch.setattr(agent_module, "build_memory_harness", lambda config: harness)

    capture = _CaptureLLM()
    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        agent = BioAgent(config=config, logger=logger, run_dir=run_dir)
        agent.llm = capture
        agent.run("Analyze single-cell data", max_turns=1)

    initial_user_message = capture.messages[1].content
    assert "Relevant BioAgent memories" in initial_user_message
    assert "prefer py-scverse-v1" in initial_user_message
