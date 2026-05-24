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


class RunScanpyPipelineArgs(BaseModel):
    data_path: str = Field(..., description="Host path or repo-relative path to an h5ad file.")
    env_profile: str = "py-scverse-v1"
    skill_id: str = "scrnaseq-scanpy-core-analysis"
    timeout_s: int = 1800


class QueryMyGeneArgs(BaseModel):
    gene_symbol: str = Field(..., description="Gene symbol, e.g. TP53 or MS4A1.")


class GeneSetEnrichmentArgs(BaseModel):
    gene_list: list[str] | str = Field(..., description="Gene symbols as a list or comma/space separated string.")
    organism: str = "human"
    databases: list[str] | None = None
    top_k: int = 10


class CellTypeAnnotationArgs(BaseModel):
    gene_list: list[str] | str = Field(..., description="Marker genes as a list or comma/space separated string.")

