from __future__ import annotations

from pathlib import Path

from langchain_core.tools import StructuredTool

from bioagent.config import AgentConfig
from bioagent.logging_utils import RunLogger
from bioagent.memory import MemoryHarness

from .execution import execute_python_impl, execute_r_impl
from .filesystem import glob_search_impl, grep_text_impl, list_files_impl, read_file_impl
from .schemas import (
    ExecutePythonArgs,
    ExecuteRArgs,
    GlobArgs,
    GrepArgs,
    InspectSkillFunctionArgs,
    InspectSkillArgs,
    InspectSkillScriptSymbolsArgs,
    ListFilesArgs,
    ListWorkflowSkillsArgs,
    ReadFileArgs,
    ReadSkillScriptArgs,
    RequestUserInputArgs,
    RunCodeWorkflowArgs,
    RunSkillWorkflowArgs,
)
from .skill_workflow import run_code_workflow_impl, run_skill_workflow_impl
from .workflow import (
    inspect_image_catalog_impl,
    inspect_skill_function_impl,
    inspect_skill_script_symbols_impl,
    inspect_workflow_skill_impl,
    list_workflow_skills_impl,
    read_skill_script_impl,
)


def build_tools(
    config: AgentConfig,
    logger: RunLogger | None,
    run_dir: Path,
    *,
    memory_harness: MemoryHarness | None = None,
) -> list:
    execution_attempts = 0

    def list_files(path: str = ".", recursive: bool = False, max_entries: int = 200) -> str:
        return list_files_impl(config, run_dir, path=path, recursive=recursive, max_entries=max_entries)

    def read_file(path: str, line_offset: int = 1, max_lines: int = 200, max_chars: int = 3000) -> str:
        return read_file_impl(config, run_dir, path=path, line_offset=line_offset, max_lines=max_lines, max_chars=max_chars)

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

    def list_workflow_skills(detail: str = "compact") -> list[dict]:
        return list_workflow_skills_impl(config, detail=detail)

    def inspect_workflow_skill(skill_id: str, max_chars: int = 3000) -> dict:
        return inspect_workflow_skill_impl(config, skill_id=skill_id, max_chars=max_chars)

    def read_skill_script(skill_id: str, script_path: str, line_offset: int = 1, max_lines: int = 200) -> str:
        return read_skill_script_impl(
            config,
            skill_id=skill_id,
            script_path=script_path,
            line_offset=line_offset,
            max_lines=max_lines,
        )

    def inspect_skill_script_symbols(skill_id: str, script_path: str, include_docstrings: bool = False) -> dict:
        return inspect_skill_script_symbols_impl(
            config,
            skill_id=skill_id,
            script_path=script_path,
            include_docstrings=include_docstrings,
        )

    def inspect_skill_function(
        skill_id: str,
        function_name: str,
        script_path: str = "",
        max_chars: int = 4000,
    ) -> dict:
        return inspect_skill_function_impl(
            config,
            skill_id=skill_id,
            function_name=function_name,
            script_path=script_path,
            max_chars=max_chars,
        )

    def inspect_image_catalog() -> list[dict]:
        return inspect_image_catalog_impl(config)

    def execute_python(code: str, env_profile: str = "py-general-v1", requirements: str = "", timeout_s: int = 900) -> dict:
        nonlocal execution_attempts
        if execution_attempts >= config.max_execution_attempts:
            return {
                "ok": False,
                "error": f"Reached global max_execution_attempts={config.max_execution_attempts}.",
                "attempts_used": execution_attempts,
            }
        execution_attempts += 1
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
        nonlocal execution_attempts
        if execution_attempts >= config.max_execution_attempts:
            return {
                "ok": False,
                "error": f"Reached global max_execution_attempts={config.max_execution_attempts}.",
                "attempts_used": execution_attempts,
            }
        execution_attempts += 1
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

    def request_user_input(
        question: str,
        reason: str = "",
        required_fields: list[str] | None = None,
        resume_hint: str = "",
    ) -> dict:
        return {
            "status": "needs_user_input",
            "question": question,
            "reason": reason,
            "required_fields": required_fields or [],
            "resume_hint": resume_hint,
        }

    tools = [
        StructuredTool.from_function(list_files, name="list_files", description="列出项目允许范围内的文件和目录。", args_schema=ListFilesArgs),
        StructuredTool.from_function(read_file, name="read_file", description="按行读取文本文件。", args_schema=ReadFileArgs),
        StructuredTool.from_function(glob_search, name="glob_search", description="按 glob 模式查找文件。", args_schema=GlobArgs),
        StructuredTool.from_function(grep_text, name="grep_text", description="搜索文件内容。", args_schema=GrepArgs),
        StructuredTool.from_function(
            request_user_input,
            name="request_user_input",
            description=(
                "当缺少关键输入且不能安全默认时，请调用此工具暂停当前 run，"
                "向用户提出一个具体问题，并等待 resume 后继续。"
            ),
            args_schema=RequestUserInputArgs,
        ),
        StructuredTool.from_function(list_workflow_skills, name="list_workflow_skills", description="列出 workflow skills。默认 compact 只返回路由字段；需要路径/env image 时用 detail='full'。", args_schema=ListWorkflowSkillsArgs),
        StructuredTool.from_function(inspect_workflow_skill, name="inspect_workflow_skill", description="读取某个 workflow skill 的元数据、scripts/references 清单和短正文预览；完整 SKILL.md 可用 read_skill_script 读取。", args_schema=InspectSkillArgs),
        StructuredTool.from_function(read_skill_script, name="read_skill_script", description="读取 workflow skill 目录内的脚本或参考文件片段，返回带行号文本。", args_schema=ReadSkillScriptArgs),
        StructuredTool.from_function(inspect_skill_script_symbols, name="inspect_skill_script_symbols", description="解析 workflow skill 的 Python 脚本。默认只列函数/类签名；需要 docstring 时设置 include_docstrings=true。", args_schema=InspectSkillScriptSymbolsArgs),
        StructuredTool.from_function(inspect_skill_function, name="inspect_skill_function", description="按函数名查看 workflow skill 脚本中的签名、docstring 和源码预览。", args_schema=InspectSkillFunctionArgs),
        StructuredTool.from_function(inspect_image_catalog, name="inspect_image_catalog", description="查看 mas_2/docker/image_catalog.json 中配置的 Docker profiles。"),
        StructuredTool.from_function(execute_python, name="execute_python", description="在指定 Python Docker profile 中执行生成的 Python 脚本。", args_schema=ExecutePythonArgs),
        StructuredTool.from_function(execute_r, name="execute_r", description="在指定 R Docker profile 中执行生成的 R 脚本。", args_schema=ExecuteRArgs),
    ]
    if config.enable_legacy_workflow_tools:
        tools.extend([
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
        ])
    if memory_harness and memory_harness.enabled:
        tools.extend(memory_harness.tools)
    return tools
