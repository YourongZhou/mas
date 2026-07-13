from __future__ import annotations

from pathlib import Path

from langchain_core.tools import StructuredTool

from bioagent.config import AgentConfig
from bioagent.logging_utils import RunLogger

from .execution import execute_python_impl, execute_r_impl
from .filesystem import glob_search_impl, grep_text_impl, list_files_impl, read_file_impl
from .schemas import (
    ExecutePythonArgs,
    ExecuteRArgs,
    GlobArgs,
    GrepArgs,
    InspectSkillArgs,
    ListFilesArgs,
    ReadFileArgs,
    RunCodeWorkflowArgs,
    RunSkillWorkflowArgs,
)
from .skill_workflow import run_code_workflow_impl, run_skill_workflow_impl
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

    def run_skill_workflow(
        skill_id: str,
        task: str,
        data_path: str = "",
        runtime: str = "",
        env_profile: str = "",
        max_attempts: int = 5,
        timeout_s: int = 1800,
    ) -> dict:
        return run_skill_workflow_impl(
            config,
            logger,
            run_dir,
            skill_id=skill_id,
            task=task,
            data_path=data_path,
            runtime=runtime,
            env_profile=env_profile,
            max_attempts=max_attempts,
            timeout_s=timeout_s,
        )

    def run_code_workflow(
        task: str,
        data_path: str = "",
        runtime: str = "python",
        env_profile: str = "",
        max_attempts: int = 5,
        timeout_s: int = 900,
    ) -> dict:
        return run_code_workflow_impl(
            config,
            logger,
            run_dir,
            task=task,
            data_path=data_path,
            runtime=runtime,
            env_profile=env_profile,
            max_attempts=max_attempts,
            timeout_s=timeout_s,
        )

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
            run_skill_workflow,
            name="run_skill_workflow",
            description=(
                "通用 Skill-driven 工作流执行器：读取指定 workflow skill，"
                "由 LLM 根据 Skill 和脚本清单生成代码，使用匹配 Docker profile 执行，失败后基于 stdout/stderr 自动修复重试。"
            ),
            args_schema=RunSkillWorkflowArgs,
        ),
        StructuredTool.from_function(
            run_code_workflow,
            name="run_code_workflow",
            description=(
                "通用无 Skill 代码工作流：当任务不需要 workflow skill 或没有匹配 skill 时，"
                "由 LLM 生成 Python/R 代码，使用 Docker 执行，并基于 stdout/stderr 自动修复重试。"
            ),
            args_schema=RunCodeWorkflowArgs,
        ),
    ]

