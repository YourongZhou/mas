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
    timeout_s: int = 1800


class ExecuteRArgs(BaseModel):
    code: str = Field(..., description="Complete R script to execute.")
    env_profile: str = Field("r-bioc-v1", description="Docker profile, e.g. r-bioc-v1 or r-singlecell-v1.")
    timeout_s: int = 1800


class WriteRunFileArgs(BaseModel):
    path: str = Field(..., description="Path inside the run workspace, e.g. scripts/analysis.py or /work/scripts/analysis.py.")
    content: str = Field(..., description="Complete text content to write.")
    overwrite: bool = Field(True, description="Allow replacing an existing run workspace file.")


class EditRunFileArgs(BaseModel):
    path: str = Field(..., description="Existing text file inside the run workspace.")
    old_text: str = Field(..., description="Exact text to replace. Include enough context to make it unique.")
    new_text: str = Field(..., description="Replacement text.")
    replace_all: bool = Field(False, description="Replace all exact matches. False requires a unique match.")


class ValidateScriptArgs(BaseModel):
    path: str = Field(..., description="Script inside the run workspace.")
    runtime: str = Field("python", description="python or r.")


class StartJobArgs(BaseModel):
    runtime: str = Field(..., description="python or r.")
    script_path: str = Field(..., description="Script path inside the run workspace. Write it before starting the job.")
    env_profile: str = Field(..., description="Verified Docker profile from the skill or image catalog.")
    requirements: str = Field("", description="Optional Python pip requirements, one per line.")
    timeout_s: int = Field(1800, description="Hard job timeout in seconds.")


class PollJobArgs(BaseModel):
    job_id: str = Field("", description="Harness-issued job id. Leave empty to use the persisted active job.")
    wait_s: int = Field(0, ge=0, le=300, description="Wait up to this many seconds before returning. Use 60-300 for long jobs.")


class ListJobsArgs(BaseModel):
    status: str = Field("", description="Optional exact status filter, e.g. running, completed, or failed.")


class GetJobArgs(BaseModel):
    job_id: str = Field("", description="Harness-issued job id. Leave empty to use the persisted active job.")
    refresh: bool = Field(False, description="Refresh a running job from Docker before returning.")


class TailJobArgs(BaseModel):
    job_id: str = Field("", description="Harness-issued job id. Leave empty to use the persisted active job.")
    lines: int = Field(200, ge=1, le=2000)


class CancelJobArgs(BaseModel):
    job_id: str = Field("", description="Harness-issued job id. Leave empty to cancel the persisted active job.")


class InspectArtifactArgs(BaseModel):
    path: str = Field(..., description="Artifact path inside the run workspace, usually outputs/...")
    max_rows: int = Field(10, ge=1, le=50)


class FinishTaskArgs(BaseModel):
    summary: str = Field(..., description="Concise conclusion supported by the inspected evidence.")
    evidence_ids: list[str] = Field(..., description="Evidence ids returned by inspect_artifact.")


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
