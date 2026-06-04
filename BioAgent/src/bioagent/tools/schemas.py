from __future__ import annotations

from pydantic import BaseModel, Field


class ReadFileArgs(BaseModel):
    path: str
    line_offset: int = 1
    max_lines: int = 200


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


class InspectSkillArgs(BaseModel):
    skill_id: str = Field(..., description="Workflow skill id, e.g. scrnaseq-scanpy-core-analysis.")
    max_chars: int = 20000


class ExecutePythonArgs(BaseModel):
    code: str = Field(..., description="Complete Python script to execute.")
    env_profile: str = Field("py-general-v1", description="Docker profile, e.g. py-general-v1 or py-scverse-v1.")
    requirements: str = Field("", description="Optional pip requirements not already provided by the profile.")
    timeout_s: int = 900


class ExecuteRArgs(BaseModel):
    code: str = Field(..., description="Complete R script to execute.")
    env_profile: str = Field("r-bioc-v1", description="Docker profile, e.g. r-bioc-v1 or r-singlecell-v1.")
    timeout_s: int = 900


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

