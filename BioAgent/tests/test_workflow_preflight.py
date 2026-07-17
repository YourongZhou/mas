from __future__ import annotations

import json
from pathlib import Path

from bioagent.config import AgentConfig
from bioagent.logging_utils import RunLogger
import bioagent.tools.skill_workflow as skill_workflow
from bioagent.tools.skill_workflow import (
    _build_generation_prompt,
    _compact_source_preview,
    _parse_json_object,
    _required_output_paths_from_task,
    run_skill_workflow_impl,
)


def _config(tmp_path: Path) -> AgentConfig:
    repo = tmp_path / "repo"
    mas2 = repo / "mas_2"
    bio = repo / "BioAgent"
    skill_dir = mas2 / "workflows" / "needs-input-skill"
    skill_dir.mkdir(parents=True)
    (mas2 / "docker").mkdir(parents=True)
    (bio / "logs").mkdir(parents=True)
    (bio / "runs").mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text(
        """---
id: needs-input-skill
name: Needs Input Skill
runtime: python
env_profile: py-general-v1
---

# Needs Input Skill

## Clarification Questions

🚨 **ALWAYS ask Question 1 FIRST.**

### 1. Input Files
- Do you have specific data files to analyze?
- Or use example/demo data?

## Standard Workflow

**Step 1 - Load data:**
""",
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


def test_skill_workflow_preflight_blocks_before_code_generation_when_required_input_missing(tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        result = run_skill_workflow_impl(
            config=config,
            logger=logger,
            run_dir=run_dir,
            skill_id="needs-input-skill",
            task="Analyze my data.",
            data_path="",
            max_attempts=1,
        )

    assert result["status"] == "needs_user_input"
    assert result["skill_id"] == "needs-input-skill"
    assert result["required_fields"] == ["data_path"]
    assert not (run_dir / "scripts" / "code.py").exists()


def test_scanpy_core_task_uses_llm_generator_with_skill_and_script_context(monkeypatch, tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    data_path = config.mas2_root / "data" / "bmmc_b_cell.h5ad"
    data_path.parent.mkdir(parents=True)
    data_path.write_text("placeholder", encoding="utf-8")
    skill_dir = config.workflows_root / "scrnaseq-scanpy-core-analysis"
    skill_dir.mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text(
        """---
id: scrnaseq-scanpy-core-analysis
name: Single-Cell RNA-seq Core Analysis (Scanpy)
runtime: python
env_profile: py-scverse-v1
---

# Single-Cell RNA-seq Core Analysis (Scanpy)

## Standard Workflow

**Step 1 — Load and QC**
**Step 2 — Normalize, reduce, cluster and visualize**
""",
        encoding="utf-8",
    )
    scripts_dir = skill_dir / "scripts"
    scripts_dir.mkdir()
    (scripts_dir / "filter_cells.py").write_text(
        "def filter_cells_by_qc(adata, min_genes=200, max_mt_percent=None, inplace=False):\n    return adata\n",
        encoding="utf-8",
    )
    (scripts_dir / "qc_metrics.py").write_text(
        "def calculate_qc_metrics(adata, species='human'):\n    return adata\n",
        encoding="utf-8",
    )

    captured: dict = {}

    def fake_generator(config, prompt, *, system_prompt=skill_workflow.GENERATOR_SYSTEM_PROMPT):
        captured["prompt"] = json.loads(prompt)
        captured["system_prompt"] = system_prompt
        return {
            "runtime": "python",
            "env_profile": "py-scverse-v1",
            "requirements": "",
            "code": "# MODEL GENERATED FROM SKILL SCRIPT CONTEXT\nprint('===RESULT===ok===')",
            "rationale": "Model generated code from the Skill and script previews.",
        }

    monkeypatch.setattr(skill_workflow, "_invoke_code_generator", fake_generator)

    class FakeDocker:
        def __init__(self, config, logger, run_dir):
            self.code = ""

        def execute_python(self, code, *, env_profile, requirements, timeout_s):
            self.code = code
            assert "MODEL GENERATED FROM SKILL SCRIPT CONTEXT" in code

            class Result:
                def to_dict(self_inner):
                    return {
                        "ok": True,
                        "runtime": "python",
                        "env_profile": env_profile,
                        "script_path": str(run_dir / "scripts" / "code.py"),
                        "exit_code": 0,
                        "stdout": "===RESULT===ok===",
                        "stderr": "",
                    }

            return Result()

    monkeypatch.setattr(skill_workflow, "DockerRunner", FakeDocker)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        result = run_skill_workflow_impl(
            config=config,
            logger=logger,
            run_dir=run_dir,
            skill_id="scrnaseq-scanpy-core-analysis",
            task="Complete core single-cell RNA-seq analysis only: QC, normalization, HVG, PCA, neighbors, Leiden, UMAP, summary JSON, key figures, processed h5ad.",
            data_path=str(data_path),
            runtime="python",
            env_profile="py-scverse-v1",
            max_attempts=1,
        )

    assert result["ok"] is True
    assert captured["system_prompt"] == skill_workflow.GENERATOR_SYSTEM_PROMPT
    prompt = captured["prompt"]
    assert prompt["skill_id"] == "scrnaseq-scanpy-core-analysis"
    assert prompt["container_data_path"] == "/repo/mas_2/data/bmmc_b_cell.h5ad"
    assert "Single-Cell RNA-seq Core Analysis" in prompt["skill_body"]
    assert "scripts/filter_cells.py" in prompt["available_script_paths"]
    selected_scripts = {item["path"]: item["preview"] for item in prompt["selected_script_context"]}
    assert "scripts/filter_cells.py" in selected_scripts
    assert "filter_cells_by_qc" in selected_scripts["scripts/filter_cells.py"]
    assert result["attempts"][0]["rationale"] == "Model generated code from the Skill and script previews."


def test_source_preview_preserves_multiline_function_signatures() -> None:
    preview = _compact_source_preview(
        '''import scanpy as sc

def plot_qc_violin(
    adata: "AnnData",
    metrics: Optional[List[str]] = None,
    output_dir: Union[str, Path] = ".",
    figsize: tuple = (12, 4)
) -> None:
    """Create violin plots of QC metrics."""
    body = "not important"
'''
    )

    assert "def plot_qc_violin(" in preview
    assert "metrics: Optional[List[str]] = None" in preview
    assert "output_dir: Union[str, Path] = \".\"" in preview
    assert "Create violin plots of QC metrics." in preview


def test_skill_workflow_retries_when_required_output_is_missing(monkeypatch, tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    data_path = config.mas2_root / "data" / "bmmc_b_cell.h5ad"
    data_path.parent.mkdir(parents=True)
    data_path.write_text("placeholder", encoding="utf-8")
    skill_dir = config.workflows_root / "scrnaseq-scanpy-core-analysis"
    skill_dir.mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text(
        """---
id: scrnaseq-scanpy-core-analysis
name: Single-Cell RNA-seq Core Analysis (Scanpy)
runtime: python
env_profile: py-scverse-v1
---

# Single-Cell RNA-seq Core Analysis (Scanpy)
""",
        encoding="utf-8",
    )

    generated_codes: list[str] = []

    def fake_generator(config, prompt, *, system_prompt=skill_workflow.GENERATOR_SYSTEM_PROMPT):
        payload = json.loads(prompt)
        generated_codes.append(payload.get("repair_context", {}).get("stderr_tail", ""))
        return {
            "runtime": "python",
            "env_profile": "py-scverse-v1",
            "requirements": "",
            "code": f"print('attempt {len(generated_codes)}')",
            "rationale": "test",
        }

    monkeypatch.setattr(skill_workflow, "_invoke_code_generator", fake_generator)

    class FakeDocker:
        calls = 0

        def __init__(self, config, logger, run_dir):
            self.run_dir = run_dir

        def execute_python(self, code, *, env_profile, requirements, timeout_s):
            FakeDocker.calls += 1
            if FakeDocker.calls == 2:
                output_path = self.run_dir / "outputs" / "qc_violin.png"
                output_path.parent.mkdir(parents=True, exist_ok=True)
                output_path.write_bytes(b"fake png")

            class Result:
                def to_dict(self_inner):
                    return {
                        "ok": True,
                        "runtime": "python",
                        "env_profile": env_profile,
                        "script_path": str(run_dir / "scripts" / "code.py"),
                        "exit_code": 0,
                        "stdout": "===RESULT===ok===",
                        "stderr": "",
                    }

            return Result()

    monkeypatch.setattr(skill_workflow, "DockerRunner", FakeDocker)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        result = run_skill_workflow_impl(
            config=config,
            logger=logger,
            run_dir=run_dir,
            skill_id="scrnaseq-scanpy-core-analysis",
            task="Create QC plot: /work/outputs/qc_violin.png",
            data_path=str(data_path),
            runtime="python",
            env_profile="py-scverse-v1",
            max_attempts=2,
        )

    assert result["ok"] is True
    assert FakeDocker.calls == 2
    assert result["attempts"][0]["ok"] is False
    assert "Missing required output files" in result["attempts"][0]["error_reason"]
    assert "qc_violin.png" in generated_codes[1]


def test_required_outputs_include_bare_filenames_when_output_dir_is_named() -> None:
    assert _required_output_paths_from_task(
        "Read input bmmc_b_cell.h5ad. Save all outputs to /work/outputs. "
        "Output files include: processed_bmmc_b_cell.h5ad, umap_clustering.png, "
        "qc_metrics.png, analysis_summary.json."
    ) == [
        "/work/outputs/processed_bmmc_b_cell.h5ad",
        "/work/outputs/umap_clustering.png",
        "/work/outputs/qc_metrics.png",
        "/work/outputs/analysis_summary.json",
    ]


def test_generation_prompt_stays_within_local_context_budget(tmp_path: Path) -> None:
    config = _config(tmp_path)
    skill_dir = config.workflows_root / "large-skill"
    scripts_dir = skill_dir / "scripts"
    refs_dir = skill_dir / "references"
    scripts_dir.mkdir(parents=True)
    refs_dir.mkdir()
    skill_body = (
        "---\nid: large-skill\nruntime: python\nenv_profile: py-scverse-v1\n---\n"
        "# Large Skill\n\n"
        + ("Use scripts for QC, HVG, normalization, PCA, UMAP, plotting, and export.\n" * 400)
    )
    (skill_dir / "SKILL.md").write_text(skill_body, encoding="utf-8")
    for index in range(16):
        (scripts_dir / f"module_{index}.py").write_text(
            "import scanpy as sc\n\n"
            "def workflow_function(\n"
            "    adata,\n"
            "    output_dir='.',\n"
            "    metrics=None,\n"
            "    resolution=0.6,\n"
            ") -> None:\n"
            '    """Run one documented workflow helper."""\n'
            + ("    pass\n" * 300),
            encoding="utf-8",
        )
    for index in range(4):
        (refs_dir / f"ref_{index}.md").write_text("reference\n" * 1000, encoding="utf-8")
    skill = {
        "path": str(skill_dir),
        "metadata": {"runtime": "python", "env_profile": "py-scverse-v1"},
        "body_preview": skill_body,
    }

    prompt = _build_generation_prompt(
        config=config,
        run_dir=config.runs_dir / "run_test",
        skill=skill,
        skill_id="large-skill",
        task="Run QC, HVG, PCA, UMAP and save /work/outputs/summary.json plus /work/outputs/umap.png.",
        runtime="python",
        env_profile="py-scverse-v1",
        host_data_path="",
        container_data_path="",
        previous_code="",
        previous_result=None,
        attempt=1,
    )
    payload = json.loads(prompt)

    assert len(prompt) < 55000
    assert payload["required_output_paths"] == ["/work/outputs/summary.json", "/work/outputs/umap.png"]
    assert payload["selected_script_context"]


def test_repair_prompt_stays_within_local_context_budget(tmp_path: Path) -> None:
    config = _config(tmp_path)
    skill_dir = config.workflows_root / "repair-skill"
    (skill_dir / "scripts").mkdir(parents=True)
    skill_body = "---\nid: repair-skill\nruntime: python\n---\n# Repair Skill\n" + ("workflow details\n" * 1000)
    (skill_dir / "SKILL.md").write_text(skill_body, encoding="utf-8")
    (skill_dir / "scripts" / "plot_qc.py").write_text(
        "def plot_qc_violin(\n"
        "    adata,\n"
        "    metrics=None,\n"
        "    output_dir='.',\n"
        ") -> None:\n"
        '    """Create QC plot."""\n',
        encoding="utf-8",
    )
    skill = {"path": str(skill_dir), "metadata": {}, "body_preview": skill_body}
    previous_result = {
        "exit_code": 1,
        "stdout": "large stdout\n" * 1000,
        "stderr": ("large traceback line\n" * 1000) + "ImportError: cannot import name 'plot_qc_violins'\n",
    }

    prompt = _build_generation_prompt(
        config=config,
        run_dir=config.runs_dir / "run_test",
        skill=skill,
        skill_id="repair-skill",
        task="Run QC and save /work/outputs/qc_violin.png.",
        runtime="python",
        env_profile="py-scverse-v1",
        host_data_path="",
        container_data_path="",
        previous_code="print('bad')\n" * 3000,
        previous_result=previous_result,
        attempt=2,
    )
    payload = json.loads(prompt)

    assert len(prompt) < 55000
    assert payload["repair_context"]["error_reason"] == "ImportError: cannot import name 'plot_qc_violins'"
    assert len(payload["repair_context"]["previous_code"]) < 9000
    assert len(payload["repair_context"]["stderr_tail"]) < 4000


def test_skill_workflow_retries_malformed_json_like_generation(monkeypatch, tmp_path: Path) -> None:
    config = _config(tmp_path)
    run_dir = config.runs_dir / "run_test"
    run_dir.mkdir(parents=True)
    data_path = config.mas2_root / "data" / "input.h5ad"
    data_path.parent.mkdir(parents=True)
    data_path.write_text("placeholder", encoding="utf-8")
    skill_dir = config.workflows_root / "json-skill"
    skill_dir.mkdir(parents=True)
    (skill_dir / "SKILL.md").write_text("---\nid: json-skill\nruntime: python\n---\n# JSON Skill\n", encoding="utf-8")
    calls = {"generator": 0, "docker": 0}

    def fake_generator(config, prompt, *, system_prompt=skill_workflow.GENERATOR_SYSTEM_PROMPT):
        calls["generator"] += 1
        if calls["generator"] == 1:
            return {
                "runtime": "python",
                "env_profile": "py-scverse-v1",
                "requirements": "",
                "code": '{"runtime": "python", "code": "print(1)"',
                "rationale": "malformed json",
            }
        return {
            "runtime": "python",
            "env_profile": "py-scverse-v1",
            "requirements": "",
            "code": "print('===RESULT===ok===')",
            "rationale": "valid code",
        }

    monkeypatch.setattr(skill_workflow, "_invoke_code_generator", fake_generator)

    class FakeDocker:
        def __init__(self, config, logger, run_dir):
            pass

        def execute_python(self, code, *, env_profile, requirements, timeout_s):
            calls["docker"] += 1
            assert not code.lstrip().startswith("{")

            class Result:
                def to_dict(self_inner):
                    return {
                        "ok": True,
                        "runtime": "python",
                        "env_profile": env_profile,
                        "script_path": str(run_dir / "scripts" / "code.py"),
                        "exit_code": 0,
                        "stdout": "===RESULT===ok===",
                        "stderr": "",
                    }

            return Result()

    monkeypatch.setattr(skill_workflow, "DockerRunner", FakeDocker)

    with RunLogger(config.logs_dir, run_id="run_test") as logger:
        result = run_skill_workflow_impl(
            config=config,
            logger=logger,
            run_dir=run_dir,
            skill_id="json-skill",
            task="Run the workflow.",
            data_path=str(data_path),
            runtime="python",
            env_profile="py-scverse-v1",
            max_attempts=2,
        )

    assert result["ok"] is True
    assert calls == {"generator": 2, "docker": 1}
    assert result["attempts"][0]["ok"] is False
    assert "malformed JSON-like content" in result["attempts"][0]["error_reason"]


def test_parse_json_object_recovers_escaped_code_from_json_like_response() -> None:
    parsed = _parse_json_object(
        '{\n'
        '  "runtime": "python",\n'
        '  "env_profile": "py-scverse-v1",\n'
        '  "code": "import sys\\nprint(\\"===RESULT===ok===\\")"\n'
        '  BROKEN_TRAILING_TEXT'
    )

    assert parsed["runtime"] == "python"
    assert parsed["env_profile"] == "py-scverse-v1"
    assert parsed["code"] == 'import sys\nprint("===RESULT===ok===")'
