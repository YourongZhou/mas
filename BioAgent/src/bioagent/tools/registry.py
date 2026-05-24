from __future__ import annotations

from pathlib import Path

from langchain_core.tools import StructuredTool

from bioagent.config import AgentConfig
from bioagent.logging_utils import RunLogger

from .execution import execute_python_impl, execute_r_impl
from .bio import gene_set_enrichment_impl, query_mygene_impl, run_celltype_annotation_impl
from .filesystem import glob_search_impl, grep_text_impl, list_files_impl, read_file_impl
from .schemas import (
    ExecutePythonArgs,
    ExecuteRArgs,
    CellTypeAnnotationArgs,
    GeneSetEnrichmentArgs,
    GlobArgs,
    GrepArgs,
    InspectSkillArgs,
    ListFilesArgs,
    QueryMyGeneArgs,
    ReadFileArgs,
    RunScanpyPipelineArgs,
)
from .singlecell import run_scanpy_singlecell_pipeline_impl
from .workflow import inspect_image_catalog_impl, inspect_workflow_skill_impl, list_workflow_skills_impl


def build_tools(config: AgentConfig, logger: RunLogger, run_dir: Path) -> list[StructuredTool]:
    def list_files(path: str = ".", recursive: bool = False, max_entries: int = 200) -> str:
        return list_files_impl(config, run_dir, path=path, recursive=recursive, max_entries=max_entries)

    def read_file(path: str, line_offset: int = 1, max_lines: int = 200) -> str:
        return read_file_impl(config, run_dir, path=path, line_offset=line_offset, max_lines=max_lines)

    def glob_search(pattern: str, path: str = ".", max_entries: int = 200) -> str:
        return glob_search_impl(config, run_dir, pattern=pattern, path=path, max_entries=max_entries)

    def grep_text(
        pattern: str,
        path: str = ".",
        glob_filter: str | None = None,
        case_insensitive: bool = False,
        max_matches: int = 100,
    ) -> str:
        return grep_text_impl(
            config,
            run_dir,
            pattern=pattern,
            path=path,
            glob_filter=glob_filter,
            case_insensitive=case_insensitive,
            max_matches=max_matches,
        )

    def list_workflow_skills() -> list[dict]:
        return list_workflow_skills_impl(config)

    def inspect_workflow_skill(skill_id: str, max_chars: int = 20000) -> dict:
        return inspect_workflow_skill_impl(config, skill_id=skill_id, max_chars=max_chars)

    def inspect_image_catalog() -> list[dict]:
        return inspect_image_catalog_impl(config)

    def execute_python(code: str, env_profile: str = "py-general-v1", requirements: str = "", timeout_s: int = 900) -> dict:
        return execute_python_impl(
            config,
            logger,
            run_dir,
            code=code,
            env_profile=env_profile,
            requirements=requirements,
            timeout_s=timeout_s,
        )

    def execute_r(code: str, env_profile: str = "r-bioc-v1", timeout_s: int = 900) -> dict:
        return execute_r_impl(config, logger, run_dir, code=code, env_profile=env_profile, timeout_s=timeout_s)

    def run_scanpy_singlecell_pipeline(
        data_path: str,
        env_profile: str = "py-scverse-v1",
        skill_id: str = "scrnaseq-scanpy-core-analysis",
        timeout_s: int = 1800,
    ) -> dict:
        return run_scanpy_singlecell_pipeline_impl(
            config,
            logger,
            run_dir,
            data_path=data_path,
            env_profile=env_profile,
            skill_id=skill_id,
            timeout_s=timeout_s,
        )

    def query_mygene(gene_symbol: str) -> dict:
        return query_mygene_impl(gene_symbol)

    def gene_set_enrichment(
        gene_list: list[str] | str,
        organism: str = "human",
        databases: list[str] | None = None,
        top_k: int = 10,
    ) -> dict:
        return gene_set_enrichment_impl(gene_list, organism=organism, databases=databases, top_k=top_k)

    def run_celltype_annotation(gene_list: list[str] | str) -> dict:
        return run_celltype_annotation_impl(gene_list)

    return [
        StructuredTool.from_function(list_files, name="list_files", description="列出项目允许范围内的文件和目录。", args_schema=ListFilesArgs),
        StructuredTool.from_function(read_file, name="read_file", description="按行读取文本文件。", args_schema=ReadFileArgs),
        StructuredTool.from_function(glob_search, name="glob_search", description="按 glob 模式查找文件。", args_schema=GlobArgs),
        StructuredTool.from_function(grep_text, name="grep_text", description="搜索文件内容。", args_schema=GrepArgs),
        StructuredTool.from_function(list_workflow_skills, name="list_workflow_skills", description="列出 mas_2/workflows 中迁移可用的 workflow skills。"),
        StructuredTool.from_function(inspect_workflow_skill, name="inspect_workflow_skill", description="读取某个 workflow skill 的元数据、脚本和文档预览。", args_schema=InspectSkillArgs),
        StructuredTool.from_function(inspect_image_catalog, name="inspect_image_catalog", description="查看 mas_2/docker/image_catalog.json 中配置的 Docker profiles。"),
        StructuredTool.from_function(execute_python, name="execute_python", description="在指定 Python Docker profile 中执行生成的 Python 脚本。", args_schema=ExecutePythonArgs),
        StructuredTool.from_function(execute_r, name="execute_r", description="在指定 R Docker profile 中执行生成的 R 脚本。", args_schema=ExecuteRArgs),
        StructuredTool.from_function(
            run_scanpy_singlecell_pipeline,
            name="run_scanpy_singlecell_pipeline",
            description="按照 scrnaseq-scanpy-core-analysis skill 对 .h5ad 单细胞数据执行完整 Scanpy 核心分析炮筒流程。",
            args_schema=RunScanpyPipelineArgs,
        ),
        StructuredTool.from_function(query_mygene, name="query_mygene", description="查询人类基因的 MyGene 信息。", args_schema=QueryMyGeneArgs),
        StructuredTool.from_function(gene_set_enrichment, name="gene_set_enrichment", description="使用 Enrichr/gseapy 对基因列表做 ORA 富集分析。", args_schema=GeneSetEnrichmentArgs),
        StructuredTool.from_function(run_celltype_annotation, name="run_celltype_annotation", description="根据 marker genes 做轻量细胞类型注释。", args_schema=CellTypeAnnotationArgs),
    ]

