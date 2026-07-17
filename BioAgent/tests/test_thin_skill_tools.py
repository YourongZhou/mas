from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from bioagent.config import AgentConfig
from bioagent.agent import SYSTEM_PROMPT
import bioagent.runner as runner_module
import bioagent.tools.registry as registry_module
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


def _write_scanpy_skill(config: AgentConfig) -> None:
    skill_dir = config.workflows_root / "scrnaseq-scanpy-core-analysis"
    scripts_dir = skill_dir / "scripts"
    scripts_dir.mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text(
        """---
id: scrnaseq-scanpy-core-analysis
name: Scanpy Core
runtime: python
env_profile: py-scverse-v1
---

# Scanpy Core
""",
        encoding="utf-8",
    )
    (scripts_dir / "run_umap.py").write_text(
        '''from __future__ import annotations


def run_umap_reduction(
    adata,
    n_neighbors=None,
    min_dist=0.5,
    spread=1.0,
    random_state=0,
    inplace=True,
):
    """Run UMAP using an existing neighbor graph."""
    return adata
''',
        encoding="utf-8",
    )


def test_default_tools_expose_thin_skill_inspection_not_legacy_workflows(tmp_path: Path) -> None:
    config = _config(tmp_path)

    tool_names = {tool.name for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    assert "inspect_skill_script_symbols" in tool_names
    assert "read_skill_script" in tool_names
    assert "inspect_skill_function" in tool_names
    assert "run_skill_workflow" not in tool_names
    assert "run_code_workflow" not in tool_names


def test_skill_script_tools_read_symbols_and_function_source(tmp_path: Path) -> None:
    config = _config(tmp_path)
    _write_scanpy_skill(config)
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    symbols = tools["inspect_skill_script_symbols"].invoke(
        {
            "skill_id": "scrnaseq-scanpy-core-analysis",
            "script_path": "scripts/run_umap.py",
        }
    )
    assert symbols["script_path"] == "scripts/run_umap.py"
    assert symbols["functions"][0]["name"] == "run_umap_reduction"
    assert "n_neighbors=None" in symbols["functions"][0]["signature"]
    assert "inplace=True" in symbols["functions"][0]["signature"]

    snippet = tools["read_skill_script"].invoke(
        {
            "skill_id": "scrnaseq-scanpy-core-analysis",
            "script_path": "scripts/run_umap.py",
            "line_offset": 4,
            "max_lines": 3,
        }
    )
    assert "4: def run_umap_reduction(" in snippet
    assert "6:     n_neighbors=None," in snippet

    function = tools["inspect_skill_function"].invoke(
        {
            "skill_id": "scrnaseq-scanpy-core-analysis",
            "function_name": "run_umap_reduction",
        }
    )
    assert function["name"] == "run_umap_reduction"
    assert function["script_path"] == "scripts/run_umap.py"
    assert "n_neighbors=None" in function["signature"]
    assert "Run UMAP using an existing neighbor graph." in function["docstring"]
    assert "n_pcs" not in function["signature"]


def test_execution_tools_share_global_attempt_budget(monkeypatch, tmp_path: Path) -> None:
    config = replace(_config(tmp_path), max_execution_attempts=1)
    calls = 0

    def fake_execute_python_impl(config, logger, run_dir, *, code, env_profile, requirements, timeout_s):
        nonlocal calls
        calls += 1
        return {"ok": True, "stdout": "ok", "stderr": ""}

    monkeypatch.setattr(registry_module, "execute_python_impl", fake_execute_python_impl)
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    first = tools["execute_python"].invoke({"code": "print('one')"})
    second = tools["execute_python"].invoke({"code": "print('two')"})

    assert first["ok"] is True
    assert second["ok"] is False
    assert "max_execution_attempts=1" in second["error"]
    assert calls == 1


def test_default_skill_inspection_context_is_small_enough_for_local_models(tmp_path: Path) -> None:
    config = _config(tmp_path)
    skill_dir = config.workflows_root / "large-skill"
    skill_dir.mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text(
        "---\nid: large-skill\n---\n" + ("workflow details\n" * 1000),
        encoding="utf-8",
    )
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    result = tools["inspect_workflow_skill"].invoke({"skill_id": "large-skill"})

    assert len(result["body_preview"]) <= 3000
    assert config.max_tool_result_chars <= 6000


def test_read_file_is_progressively_disclosed_by_default(tmp_path: Path) -> None:
    config = _config(tmp_path)
    large_file = config.repo_root / "large-reference.md"
    large_file.write_text("\n".join(f"line {idx} " + ("x" * 80) for idx in range(1, 301)), encoding="utf-8")
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    default = tools["read_file"].invoke({"path": str(large_file)})
    expanded = tools["read_file"].invoke({"path": str(large_file), "max_lines": 300, "max_chars": 40000})

    assert len(default) <= 3500
    assert "[truncated" in default
    assert len(expanded) > len(default)
    assert "300 | line 300" in expanded


def test_workflow_skill_list_is_compact_by_default_and_full_on_request(tmp_path: Path) -> None:
    config = _config(tmp_path)
    skill_dir = config.workflows_root / "example-skill"
    skill_dir.mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text(
        """---
id: example-skill
name: Example Skill
category: testing
short-description: A deliberately verbose description for routing.
runtime: python
env_profile: py-general-v1
env_image: mas/py-general:v1
---

# Example
""",
        encoding="utf-8",
    )
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    compact = tools["list_workflow_skills"].invoke({})
    full = tools["list_workflow_skills"].invoke({"detail": "full"})

    assert compact == [
        {
            "skill_id": "example-skill",
            "name": "Example Skill",
            "short_description": "A deliberately verbose description for routing.",
            "runtime": "python",
            "env_profile": "py-general-v1",
        }
    ]
    assert full[0]["path"].endswith("example-skill")
    assert full[0]["env_image"] == "mas/py-general:v1"


def test_script_symbols_are_signatures_only_unless_docstrings_requested(tmp_path: Path) -> None:
    config = _config(tmp_path)
    _write_scanpy_skill(config)
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    compact = tools["inspect_skill_script_symbols"].invoke(
        {
            "skill_id": "scrnaseq-scanpy-core-analysis",
            "script_path": "scripts/run_umap.py",
        }
    )
    verbose = tools["inspect_skill_script_symbols"].invoke(
        {
            "skill_id": "scrnaseq-scanpy-core-analysis",
            "script_path": "scripts/run_umap.py",
            "include_docstrings": True,
        }
    )

    assert "docstring" not in compact["functions"][0]
    assert verbose["functions"][0]["docstring"] == "Run UMAP using an existing neighbor graph."


def test_read_skill_script_can_explicitly_return_full_skill_file(tmp_path: Path) -> None:
    config = _config(tmp_path)
    _write_scanpy_skill(config)
    tools = {tool.name: tool for tool in build_tools(config, logger=None, run_dir=tmp_path / "run")}

    text = tools["read_skill_script"].invoke(
        {
            "skill_id": "scrnaseq-scanpy-core-analysis",
            "script_path": "SKILL.md",
            "line_offset": 1,
            "max_lines": 1000,
        }
    )

    assert "1: ---" in text
    assert "# Scanpy Core" in text


def test_prompt_and_demo_do_not_route_to_legacy_workflow_tools(monkeypatch) -> None:
    assert "run_skill_workflow" not in SYSTEM_PROMPT
    assert "run_code_workflow" not in SYSTEM_PROMPT
    assert "Do not invent skill script paths" in SYSTEM_PROMPT
    assert "sc.read_h5ad" in SYSTEM_PROMPT
    assert "Do not inspect every script" in SYSTEM_PROMPT
    assert "Do not call sc.tl.rank_genes_groups" in SYSTEM_PROMPT

    captured: dict = {}

    def fake_run_bio_agent(task, *, data_path=None, result_dir=None, max_turns=20):
        captured["task"] = task
        captured["data_path"] = data_path
        captured["max_turns"] = max_turns
        return {"status": "captured"}

    monkeypatch.setattr(runner_module, "run_bio_agent", fake_run_bio_agent)

    result = runner_module.run_bmmc_singlecell_demo(max_turns=7)

    assert result["status"] == "captured"
    assert "run_skill_workflow" not in captured["task"]
    assert "inspect_skill_function" in captured["task"]
    assert "execute_python" in captured["task"]
    assert "rank_genes_groups" in captured["task"]
    assert "dense numpy AnnData" in captured["task"]
    assert "fig.savefig" in captured["task"]
