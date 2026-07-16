from __future__ import annotations

from pydantic import BaseModel, Field


class ReadFileArgs(BaseModel):
    path: str
    line_offset: int = 1
    max_lines: int = 200
    max_chars: int = Field(3000, description="Maximum characters returned by default; set higher with line ranges when more detail is needed.")


class ListFilesArgs(BaseModel):
    path: str = "."
    recursive: bool = False
    max_entries: int = 200


class GlobArgs(BaseModel):
    pattern: str
    path: str = "."
    max_entries: int = 200


class GrepArgs(BaseModel):
    pattern: str
    path: str = "."
    glob_filter: str | None = None
    case_insensitive: bool = False
    max_matches: int = 100


class ListWorkflowSkillsArgs(BaseModel):
    detail: str = Field("compact", description="Use compact for routing fields only; use full when absolute paths/env images are needed.")


class InspectSkillArgs(BaseModel):
    skill_id: str = Field(..., description="Workflow skill id, e.g. scrnaseq-scanpy-core-analysis.")
    max_chars: int = 3000


class ReadSkillScriptArgs(BaseModel):
    skill_id: str = Field(..., description="Workflow skill id, e.g. scrnaseq-scanpy-core-analysis.")
    script_path: str = Field(..., description="Path relative to the skill directory, e.g. scripts/run_umap.py.")
    line_offset: int = 1
    max_lines: int = 200


class InspectSkillScriptSymbolsArgs(BaseModel):
    skill_id: str = Field(..., description="Workflow skill id, e.g. scrnaseq-scanpy-core-analysis.")
    script_path: str = Field(..., description="Python script path relative to the skill directory.")
    include_docstrings: bool = Field(False, description="Keep false for signature-only routing; true adds docstrings when needed.")


class InspectSkillFunctionArgs(BaseModel):
    skill_id: str = Field(..., description="Workflow skill id, e.g. scrnaseq-scanpy-core-analysis.")
    function_name: str = Field(..., description="Function name to inspect in the skill scripts.")
    script_path: str = Field("", description="Optional script path relative to the skill directory. Empty means search scripts/*.py.")
    max_chars: int = 4000


class ExecutePythonArgs(BaseModel):
    code: str = Field(..., description="Complete Python script to execute.")
    env_profile: str = Field("py-general-v1", description="Docker profile, e.g. py-general-v1 or py-scverse-v1.")
    requirements: str = Field("", description="Optional pip requirements not already provided by the profile.")
    timeout_s: int = 900


class ExecuteRArgs(BaseModel):
    code: str = Field(..., description="Complete R script to execute.")
    env_profile: str = Field("r-bioc-v1", description="Docker profile, e.g. r-bioc-v1 or r-singlecell-v1.")
    timeout_s: int = 900


class RequestUserInputArgs(BaseModel):
    question: str = Field(..., description="Concrete question that must be answered before the workflow can proceed.")
    reason: str = Field("", description="Why this answer is required now.")
    required_fields: list[str] = Field(default_factory=list, description="Structured missing fields, e.g. species or data_path.")
    resume_hint: str = Field("", description="Short instruction for how the answer should be used when the run resumes.")


class RunSkillWorkflowArgs(BaseModel):
    skill_id: str = Field(..., description="Workflow skill id from mas_2/workflows.")
    task: str = Field(..., description="Concrete analysis task and expected outputs.")
    data_path: str = Field("", description="Optional host path or repo-relative primary input data path.")
    runtime: str = Field("", description="Optional runtime override: python or r. Empty means use skill metadata.")
    env_profile: str = Field("", description="Optional Docker profile override. Empty means use skill metadata.")
    max_attempts: int = Field(5, description="Maximum generate/execute/repair attempts.")
    timeout_s: int = Field(1800, description="Per-execution timeout in seconds.")


class RunCodeWorkflowArgs(BaseModel):
    task: str = Field(..., description="Concrete analysis or data-processing task and expected outputs.")
    data_path: str = Field("", description="Optional host path or repo-relative primary input data path.")
    runtime: str = Field("python", description="Runtime: python or r.")
    env_profile: str = Field("", description="Optional Docker profile override. Empty means py-general-v1 for Python or r-bioc-v1 for R.")
    max_attempts: int = Field(5, description="Maximum generate/execute/repair attempts.")
    timeout_s: int = Field(900, description="Per-execution timeout in seconds.")
