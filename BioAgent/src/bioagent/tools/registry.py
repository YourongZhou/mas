from __future__ import annotations

from pathlib import Path

from langchain_core.tools import StructuredTool

from bioagent.config import AgentConfig
from bioagent.logging_utils import RunLogger
from bioagent.memory import MemoryHarness

from .artifacts import ArtifactStore
from .attempts import ExecutionAttemptBudget
from .execution import execute_python_impl, execute_r_impl
from .filesystem import glob_search_impl, grep_text_impl, list_files_impl, read_file_impl
from .jobs import DockerJobManager
from .schemas import (
    CancelJobArgs,
    EditRunFileArgs,
    ExecutePythonArgs,
    ExecuteRArgs,
    FinishTaskArgs,
    GetJobArgs,
    GlobArgs,
    GrepArgs,
    InspectArtifactArgs,
    InspectSkillFunctionArgs,
    InspectSkillArgs,
    InspectSkillScriptSymbolsArgs,
    ListFilesArgs,
    ListJobsArgs,
    ListWorkflowSkillsArgs,
    PollJobArgs,
    ReadFileArgs,
    ReadSkillScriptArgs,
    RequestUserInputArgs,
    RunCodeWorkflowArgs,
    RunSkillWorkflowArgs,
    StartJobArgs,
    TailJobArgs,
    ValidateScriptArgs,
    WriteRunFileArgs,
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
from .workspace import edit_run_file_impl, validate_script_impl, write_run_file_impl


def build_tools(
    config: AgentConfig,
    logger: RunLogger | None,
    run_dir: Path,
    *,
    memory_harness: MemoryHarness | None = None,
    prior_run_dirs: list[Path] | None = None,
) -> list:
    attempt_budget = ExecutionAttemptBudget(run_dir, config.max_execution_attempts)
    job_manager = DockerJobManager(config, logger, run_dir, prior_run_dirs=prior_run_dirs)
    artifact_store = ArtifactStore(run_dir, config=config, prior_run_dirs=prior_run_dirs)

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

    def execute_python(code: str, env_profile: str = "py-general-v1", requirements: str = "", timeout_s: int = 1800) -> dict:
        attempt = attempt_budget.consume(tool="execute_python")
        if not attempt.get("ok"):
            return attempt
        return execute_python_impl(
            config,
            logger,
            run_dir,
            code=code,
            env_profile=env_profile,
            requirements=requirements,
            timeout_s=timeout_s,
        )

    def execute_r(code: str, env_profile: str = "r-bioc-v1", timeout_s: int = 1800) -> dict:
        attempt = attempt_budget.consume(tool="execute_r")
        if not attempt.get("ok"):
            return attempt
        return execute_r_impl(config, logger, run_dir, code=code, env_profile=env_profile, timeout_s=timeout_s)

    def write_run_file(path: str, content: str, overwrite: bool = True) -> dict:
        return write_run_file_impl(run_dir, path=path, content=content, overwrite=overwrite)

    def edit_run_file(path: str, old_text: str, new_text: str, replace_all: bool = False) -> dict:
        return edit_run_file_impl(
            run_dir,
            path=path,
            old_text=old_text,
            new_text=new_text,
            replace_all=replace_all,
        )

    def validate_script(path: str, runtime: str = "python") -> dict:
        return validate_script_impl(run_dir, path=path, runtime=runtime)

    def start_job(
        runtime: str,
        script_path: str,
        env_profile: str,
        requirements: str = "",
        timeout_s: int = 1800,
    ) -> dict:
        return job_manager.start_job(
            runtime=runtime,
            script_path=script_path,
            env_profile=env_profile,
            requirements=requirements,
            timeout_s=timeout_s,
        )

    def list_jobs(status: str = "") -> dict:
        return job_manager.list_jobs(status=status)

    def get_job(job_id: str = "", refresh: bool = False) -> dict:
        return job_manager.get_job(job_id, refresh=refresh)

    def poll_job(job_id: str = "", wait_s: int = 0) -> dict:
        return job_manager.poll_job(job_id, wait_s=wait_s)

    def tail_job(job_id: str = "", lines: int = 200) -> dict:
        return job_manager.tail_job(job_id, lines=lines)

    def cancel_job(job_id: str = "") -> dict:
        return job_manager.cancel_job(job_id)

    def inspect_artifact(path: str, max_rows: int = 10) -> dict:
        return artifact_store.inspect(path, max_rows=max_rows)

    def finish_task(summary: str, evidence_ids: list[str]) -> dict:
        return artifact_store.finish_task(summary=summary, evidence_ids=evidence_ids)

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
        StructuredTool.from_function(write_run_file, name="write_run_file", description="把持久脚本或文本写入当前 run workspace；返回路径而不是把代码带进后续上下文。", args_schema=WriteRunFileArgs),
        StructuredTool.from_function(edit_run_file, name="edit_run_file", description="对 run workspace 文件做精确文本替换，适合根据报错进行最小修复。", args_schema=EditRunFileArgs),
        StructuredTool.from_function(validate_script, name="validate_script", description="执行不计 attempt 的轻量脚本语法检查。Python 使用 AST parse；R 完整 parse 在 start_job 中执行。", args_schema=ValidateScriptArgs),
        StructuredTool.from_function(start_job, name="start_job", description="按当前 run workspace 中的脚本路径启动持久异步 Docker job，并立即返回 job_id；同一 session 已有运行中 job 时会拒绝重复启动。", args_schema=StartJobArgs),
        StructuredTool.from_function(list_jobs, name="list_jobs", description="列出当前 conversation session 各 run 中由 harness 签发并持久化的 jobs，以及 active_job_id。", args_schema=ListJobsArgs),
        StructuredTool.from_function(get_job, name="get_job", description="读取 session 内一个已登记 job 的持久状态；job_id 留空时优先使用仍在运行的 active job。", args_schema=GetJobArgs),
        StructuredTool.from_function(poll_job, name="poll_job", description="查询异步 job；job_id 留空时使用持久 active job。长任务可设置 wait_s=60..300；返回最新 stdout/stderr。", args_schema=PollJobArgs),
        StructuredTool.from_function(tail_job, name="tail_job", description="按需读取异步 job 的最新 Docker 日志。", args_schema=TailJobArgs),
        StructuredTool.from_function(cancel_job, name="cancel_job", description="终止指定异步 Docker job。", args_schema=CancelJobArgs),
        StructuredTool.from_function(inspect_artifact, name="inspect_artifact", description="通用检查 run 输出文件并登记可复核 evidence；支持表格、JSON、图片、gzip 和文本。", args_schema=InspectArtifactArgs),
        StructuredTool.from_function(finish_task, name="finish_task", description="用 inspect_artifact 返回的 evidence_ids 验证最终产物并形成 grounded final answer。", args_schema=FinishTaskArgs),
        StructuredTool.from_function(execute_python, name="execute_python", description="兼容工具：直接执行完整 Python 代码。优先使用 write_run_file + start_job。", args_schema=ExecutePythonArgs),
        StructuredTool.from_function(execute_r, name="execute_r", description="兼容工具：直接执行完整 R 代码。优先使用 write_run_file + start_job。", args_schema=ExecuteRArgs),
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
    InspectArtifactArgs,
    PollJobArgs,
