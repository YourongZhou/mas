from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

from langchain_anthropic import ChatAnthropic
from langchain_openai import ChatOpenAI

from bioagent.config import AgentConfig
from bioagent.llm import build_llm, message_text, runtime_summary


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


def test_build_llm_uses_openai_adapter_by_default(tmp_path: Path) -> None:
    llm = build_llm(_config(tmp_path))

    assert isinstance(llm, ChatOpenAI)


def test_build_llm_uses_anthropic_adapter_and_normalizes_v1_base_url(tmp_path: Path) -> None:
    config = replace(
        _config(tmp_path),
        provider="anthropic",
        base_url="https://apiclaude.cc/v1",
        model_name="claude-sonnet-4-6",
    )

    llm = build_llm(config)

    assert isinstance(llm, ChatAnthropic)
    assert llm.anthropic_api_url == "https://apiclaude.cc"
    assert llm.model == "claude-sonnet-4-6"
    assert "provider=anthropic" in runtime_summary(config)


def test_message_text_extracts_only_anthropic_text_blocks() -> None:
    message = SimpleNamespace(
        content=[
            {"type": "text", "text": "I will inspect the skill."},
            {
                "type": "tool_use",
                "id": "toolu_123",
                "name": "inspect_workflow_skill",
                "input": {"skill_id": "scrnaseq-scanpy-core-analysis"},
            },
        ]
    )

    assert message_text(message) == "I will inspect the skill."
