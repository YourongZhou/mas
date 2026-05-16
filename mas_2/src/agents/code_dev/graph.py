"""
Code Developer Agent 子图
负责代码生成和执行
"""
import os
import re
import tempfile
import time
import uuid
from datetime import datetime
from pathlib import Path
from langchain_core.messages import SystemMessage, HumanMessage, ToolMessage
from langgraph.graph import StateGraph, START, END
from .state import CodeAgentState
from src.core.llm import get_llm
from src.sandbox.tools import build_code_dev_sandbox_tool_bundle
from .executor import CodeExecutor
from src.utils.docker_log_summary import summarize_docker_stdout
from src.utils.project_paths import get_mas2_project_root
from src.utils.workflow_skills import (
    format_skill_injection_for_code_dev,
    resolve_workflow_root,
    should_mount_workflow_in_docker,
    get_workflows_root,
)
from src.utils.execution_environment import (
    asset_cache_status,
    clear_environment_from_state,
    ensure_asset_cache,
    filter_profiled_requirements,
    format_environment_package_versions_for_prompt,
    resolve_environment_for_state,
    sync_environment_to_state,
)
from ._utils.docker_path import convert_to_docker_path
from ._utils.llm_code_sanitize import sanitize_llm_python_block
from ._utils.base64_support import create_html_with_base64_image
from .prompt import get_code_system_prompt, get_code_user_prompt, get_tool_system_prompt, get_tool_user_prompt, get_final_user_prompt
# 初始化 LLM
llm = get_llm(temperature=0.1)

_MAX_PENDING_ERROR_CHARS = 12000


def _mas2_data_dir_candidates() -> list[str]:
    """与进程 CWD 无关的数据目录候选（优先 mas_2/data，再 ./data）。"""
    roots: list[str] = []
    try:
        roots.append(str(get_mas2_project_root() / "data"))
    except RuntimeError:
        pass
    roots.append(os.path.abspath("./data"))
    return roots


def _env_truthy(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in ("1", "true", "yes", "on")


def _resolve_input_file_host_path(input_ref: str, data_path: str, result_path: str) -> str | None:
    ref = str(input_ref or "").strip()
    if not ref:
        return None

    candidates: list[str] = []
    if os.path.isabs(ref):
        candidates.append(ref)
    else:
        candidates.extend(
            [
                os.path.abspath(ref),
                os.path.join(result_path, ref) if result_path else "",
                os.path.join(result_path, os.path.basename(ref)) if result_path else "",
            ]
        )
        if data_path:
            if os.path.isdir(data_path):
                candidates.extend(
                    [
                        os.path.join(data_path, ref),
                        os.path.join(data_path, os.path.basename(ref)),
                    ]
                )
            elif os.path.isfile(data_path) and os.path.basename(data_path) == os.path.basename(ref):
                candidates.append(data_path)
        for root_candidate in _mas2_data_dir_candidates():
            if root_candidate:
                candidates.append(os.path.join(root_candidate, os.path.basename(ref)))

    seen: set[str] = set()
    for candidate in candidates:
        text = str(candidate or "").strip()
        if not text:
            continue
        abs_path = os.path.abspath(text)
        if abs_path in seen:
            continue
        seen.add(abs_path)
        if os.path.exists(abs_path):
            return abs_path
    return None


def _build_input_mount_plan(data_path: str, result_path: str, input_files: list[str]) -> dict:
    result_root = os.path.abspath(result_path) if result_path else ""
    data_dirs: list[str] = []
    dir_to_index: dict[str, int] = {}
    mappings: list[dict[str, str]] = []
    required_container_paths: list[str] = []
    missing_inputs: list[str] = []

    def _register_data_dir(host_dir: str) -> int:
        abs_dir = os.path.abspath(host_dir)
        if abs_dir not in dir_to_index:
            dir_to_index[abs_dir] = len(data_dirs)
            data_dirs.append(abs_dir)
        return dir_to_index[abs_dir]

    for raw_ref in input_files or []:
        ref = str(raw_ref or "").strip()
        if not ref:
            continue
        resolved = _resolve_input_file_host_path(ref, data_path, result_path)
        container_path = f"/app/data/{os.path.basename(ref)}"

        if resolved:
            abs_resolved = os.path.abspath(resolved)
            if result_root and os.path.commonpath([result_root, abs_resolved]) == result_root:
                rel_path = os.path.relpath(abs_resolved, result_root).replace("\\", "/")
                container_path = f"/app/output/{rel_path}"
            elif os.path.isfile(abs_resolved):
                dir_index = _register_data_dir(os.path.dirname(abs_resolved))
                container_root = "/app/data" if dir_index == 0 else f"/app/data{dir_index}"
                container_path = f"{container_root}/{os.path.basename(abs_resolved)}"
            elif os.path.isdir(abs_resolved):
                dir_index = _register_data_dir(abs_resolved)
                container_path = "/app/data" if dir_index == 0 else f"/app/data{dir_index}"
        else:
            missing_inputs.append(ref)

        mappings.append(
            {
                "input_ref": ref,
                "host_path": resolved or "",
                "container_path": container_path,
            }
        )
        required_container_paths.append(container_path)

    return {
        "data_dirs": data_dirs,
        "mappings": mappings,
        "required_container_paths": required_container_paths,
        "missing_inputs": missing_inputs,
    }


def _exec_output_tail(text: str, n_lines: int) -> str:
    if not text:
        return ""
    lines = text.splitlines()
    if len(lines) <= n_lines:
        return text
    return "\n".join(lines[-n_lines:])


def _extract_informative_error(text: str, max_chars: int = 2000) -> str:
    """从长日志中提取最有信息量的末尾报错片段。"""
    raw = (text or "").strip()
    if not raw:
        return ""

    tb_matches = list(
        re.finditer(
            r"Traceback \(most recent call last\):[\s\S]*?(?=\n\n|\Z)",
            raw,
            re.DOTALL,
        )
    )
    if tb_matches:
        return tb_matches[-1].group(0).strip()[-max_chars:]

    lines = raw.splitlines()
    err_keys = (
        "ModuleNotFoundError",
        "ImportError",
        "TypeError",
        "ValueError",
        "AttributeError",
        "NameError",
        "KeyError",
        "IndexError",
        "RuntimeError",
        "Exception",
        "Error",
        "failed",
        "unable to",
    )
    for i in range(len(lines) - 1, -1, -1):
        if any(k.lower() in lines[i].lower() for k in err_keys):
            start = max(0, i - 6)
            snippet = "\n".join(lines[start : i + 1]).strip()
            return snippet[-max_chars:]

    return _exec_output_tail(raw, 40)[-max_chars:]


def _tail_line_count() -> int:
    try:
        return max(10, int(os.environ.get("MAS_EXEC_OUTPUT_TAIL_LINES", "80")))
    except ValueError:
        return 80


def _maybe_write_full_exec_log(output_str: str) -> str | None:
    """MAS_SAVE_FULL_EXEC_LOG=1 时将完整 stdout 写入 mas_2/results，返回绝对路径。"""
    if not output_str or not _env_truthy("MAS_SAVE_FULL_EXEC_LOG"):
        return None
    root = Path(__file__).resolve().parents[3]
    results_dir = root / "results"
    try:
        results_dir.mkdir(parents=True, exist_ok=True)
        log_name = f"exec_full_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}.log"
        path = results_dir / log_name
        path.write_text(output_str, encoding="utf-8")
        return str(path.resolve())
    except OSError:
        return None


def _build_execute_pending_contribution(
    *,
    code: str,
    requirements: str,
    task: str,
    output_str: str,
    success: bool,
    result_value: str | None = None,
    error_msg: str | None = None,
    output_files: list | None = None,
    timing: dict | None = None,
    phase_timing: dict | None = None,
) -> dict:
    """
    执行后写入 pending_contribution：默认不含完整 output，含 output_display / output_tail；
    MAS_KEEP_FULL_EXEC_OUTPUT_IN_STATE=1 时附带 output；MAS_SAVE_FULL_EXEC_LOG=1 时落盘并写 output_log_path。
    """
    out = output_str or ""
    tail_n = _tail_line_count()
    output_tail = _exec_output_tail(out, tail_n)
    disp = summarize_docker_stdout(out)
    log_path = _maybe_write_full_exec_log(out)
    pending: dict = {
        "code": code,
        "requirements": requirements,
        "task": task,
        "success": success,
        "output_display": disp,
        "output_tail": output_tail,
    }
    if log_path:
        pending["output_log_path"] = log_path
    if _env_truthy("MAS_KEEP_FULL_EXEC_OUTPUT_IN_STATE"):
        pending["output"] = out
    if timing:
        pending["timing"] = timing
    if phase_timing:
        pending["phase_timing"] = phase_timing
    if success:
        pending["result"] = (result_value or "").strip()
        if output_files is not None:
            pending["output_files"] = output_files
    else:
        err = (error_msg or "").strip()
        if len(err) > _MAX_PENDING_ERROR_CHARS:
            err = (
                err[:_MAX_PENDING_ERROR_CHARS]
                + "\n...(pending.error 已截断；见 output_tail / output_log_path 或设置 MAS_KEEP_FULL_EXEC_OUTPUT_IN_STATE)"
            )
        pending["error"] = err
    return pending


def _format_executor_timing(timing: dict | None) -> str:
    if not timing:
        return "无"

    def _fmt(key: str) -> str:
        value = timing.get(key, 0.0)
        try:
            return f"{float(value):.2f}s"
        except (TypeError, ValueError):
            return "n/a"

    install_value = _fmt("pip_install_elapsed_seconds")
    if timing.get("pip_install_skipped", False):
        install_value += " (skipped)"

    parts = [
        f"container_setup={_fmt('container_setup_elapsed_seconds')}",
        f"push_to_container={_fmt('push_to_container_elapsed_seconds')}",
        f"pip_install={install_value}",
        f"python_exec={_fmt('python_exec_elapsed_seconds')}",
        f"collect_outputs={_fmt('collect_outputs_elapsed_seconds')}",
        f"total={_fmt('total_elapsed_seconds')}",
    ]

    install_exit = timing.get("install_exit_code")
    python_exit = timing.get("python_exit_code")
    if install_exit is not None:
        parts.append(f"install_exit_code={install_exit}")
    if python_exit is not None:
        parts.append(f"python_exit_code={python_exit}")
    return ", ".join(parts)


def _fmt_elapsed(value) -> str:
    try:
        return f"{float(value):.2f}s"
    except (TypeError, ValueError):
        return "n/a"


def _ensure_attempt_timing(state: CodeAgentState, attempt_no: int, is_retry: bool) -> dict:
    history = state.get("code_dev_attempt_timings")
    if not isinstance(history, list):
        history = []
        state["code_dev_attempt_timings"] = history

    current = {
        "attempt": attempt_no,
        "is_retry": is_retry,
    }
    history.append(current)
    state["current_attempt_timing"] = current
    return current


def _format_tool_stats(tool_stats: dict[str, dict[str, float | int]]) -> str:
    if not tool_stats:
        return "无"
    parts = []
    for name in sorted(tool_stats.keys()):
        info = tool_stats[name]
        parts.append(
            f"{name}:count={int(info.get('count', 0))},elapsed={_fmt_elapsed(info.get('elapsed_seconds', 0.0))}"
        )
    return "; ".join(parts)


def _format_phase_timing(phase_timing: dict | None) -> str:
    if not isinstance(phase_timing, dict) or not phase_timing:
        return "无"

    parts = [
        f"attempt={phase_timing.get('attempt', '?')}",
        f"is_retry={bool(phase_timing.get('is_retry', False))}",
    ]

    gen = phase_timing.get("generate_code")
    if isinstance(gen, dict):
        parts.append(
            "generate_code="
            + ",".join(
                [
                    f"prompt_prep={_fmt_elapsed(gen.get('prompt_prep_elapsed_seconds'))}",
                    f"tool_loop={_fmt_elapsed(gen.get('tool_loop_elapsed_seconds'))}",
                    f"tool_llm={_fmt_elapsed(gen.get('tool_llm_elapsed_seconds'))}",
                    f"tool_io={_fmt_elapsed(gen.get('tool_io_elapsed_seconds'))}",
                    f"llm_generate={_fmt_elapsed(gen.get('llm_generate_elapsed_seconds'))}",
                    f"parse={_fmt_elapsed(gen.get('parse_elapsed_seconds'))}",
                    f"total={_fmt_elapsed(gen.get('total_elapsed_seconds'))}",
                    f"tool_iterations={gen.get('tool_iterations', 0)}",
                    f"tool_calls={gen.get('tool_calls', 0)}",
                    f"reused_tool_context={bool(gen.get('reused_tool_context', False))}",
                ]
            )
        )

    reflection = phase_timing.get("self_reflection")
    if isinstance(reflection, dict):
        parts.append(
            "self_reflection="
            + ",".join(
                [
                    f"total={_fmt_elapsed(reflection.get('total_elapsed_seconds'))}",
                    f"warnings={reflection.get('warnings_count', 0)}",
                ]
            )
        )

    execute = phase_timing.get("execute_code")
    if isinstance(execute, dict):
        parts.append(
            "execute_code="
            + ",".join(
                [
                    f"prep={_fmt_elapsed(execute.get('prep_elapsed_seconds'))}",
                    f"result_parse={_fmt_elapsed(execute.get('result_parse_elapsed_seconds'))}",
                    f"total={_fmt_elapsed(execute.get('total_elapsed_seconds'))}",
                ]
            )
        )

    retry = phase_timing.get("prepare_retry")
    if isinstance(retry, dict):
        parts.append(
            "prepare_retry="
            + ",".join(
                [
                    f"feedback_extract={_fmt_elapsed(retry.get('feedback_extract_elapsed_seconds'))}",
                    f"total={_fmt_elapsed(retry.get('total_elapsed_seconds'))}",
                ]
            )
        )

    return " | ".join(parts)


def parse_paths_from_query(user_query: str) -> dict:
    """
    从用户查询中解析数据路径和结果路径
    """
    paths = {
        "data_path": "",
        "result_path": ""
    }

    if not user_query:
        return paths

    # 1. 数据路径的正则模式
    # 修复：data_path 模式增加了 ([^\n]+?) 捕获组
    data_patterns = [
        r'数据路径[：:]\s*([^\n]+)',
        r'data[_\s]?path[：:]\s*([^\n]+)',
        r'输入路径[：:]\s*([^\n]+)',
        r'input[_\s]?path[：:]\s*([^\n]+)',
    ]

    # 2. 结果路径的正则模式
    # 优化：去掉了过于激进的关键词 Lookahead 断言，改用更稳健的行结束匹配
    # 这样即使路径里有 /home/data/ 这种词也不会被截断
    result_patterns = [
        r'输出路径[：:]\s*([^\n]+)',
        r'结果路径[：:]\s*([^\n]+)',
        r'output[_\s]?path[：:]\s*([^\n]+)',
        r'result[_\s]?path[：:]\s*([^\n]+)',
        r'保存路径[：:]\s*([^\n]+)',
    ]

    # 尝试匹配数据路径
    for pattern in data_patterns:
        match = re.search(pattern, user_query, re.IGNORECASE)
        if match:
            # 拿到匹配内容后，通过 strip() 移除行尾可能的空白
            path_str = match.group(1).strip()
            # 清理引号和末尾多余的空格
            path_str = path_str.strip('"\' ')
            paths["data_path"] = path_str
            break

    # 兜底识别：从自然语言中提取常见数据文件路径（含组学常见后缀）优先于单纯目录
    if not paths["data_path"]:
        file_matches = re.findall(
            r"([A-Za-z0-9_./:\-\\]+\.(?:vcf\.gz|vcf\.bgz|bcf|vcf\.gz\.tbi|vcf|panel|bed|bim|fam|h5ad|h5|csv|tsv|mtx|loom|txt|json))",
            user_query,
            re.IGNORECASE,
        )
        if file_matches:
            paths["data_path"] = ", ".join([f.strip('"\' ') for f in file_matches])

    # 兜底识别：直接匹配所有长相像绝对路径或文件夹的描述
    if not paths["data_path"]:
        # 匹配位于 ... 下，如 “at E:\... ” 或者 “位于 /home/...”
        loc_match = re.search(r'(?:位于|at|using files at|in|从)\s+([A-Za-z0-9_./:\-\\]+)', user_query, re.IGNORECASE)
        if loc_match:
            paths["data_path"] = loc_match.group(1).strip('"\' ')

    # 兜底识别：直接匹配绝对路径（Linux的 / 或 Windows的 C:\ ）
    if not paths["data_path"]:
        abs_match = re.search(r'([a-zA-Z]:[\\/][A-Za-z0-9_./:\-\\]+|/[A-Za-z0-9_./\-\\]+)', user_query)
        if abs_match:
            paths["data_path"] = abs_match.group(1).strip('"\' ')

    # 兜底识别：如“data/ 目录下”这类描述（如果前面都没匹配到文件）
    if not paths["data_path"]:
        dir_match = re.search(r'([A-Za-z0-9_./:\-\\]+[/\\]?)\s*目录', user_query, re.IGNORECASE)
        if dir_match:
            dir_path = dir_match.group(1).strip('"\' ')
            paths["data_path"] = dir_path

    # 尝试匹配结果路径
    for pattern in result_patterns:
        match = re.search(pattern, user_query, re.IGNORECASE)
        if match:
            path_str = match.group(1).strip()
            path_str = path_str.strip('"\' ')
            paths["result_path"] = path_str
            break

    return paths


def extract_paths_from_state(state: CodeAgentState) -> CodeAgentState:
    """
    从 state 中提取路径信息
    如果 state 中没有路径，则从 user_query 中解析

    优先级：
    1. state 中已有的路径（如果存在）
    2. 从 user_query 中解析的路径
    3. 默认值（result_path 会为每个任务生成唯一的隔离目录，例如 "./result/<task_id>"）
    """
    import uuid
    # 获取或生成任务唯一 ID
    task_id = state.get("task_id")
    if not task_id:
        task_id = uuid.uuid4().hex[:8]
        state["task_id"] = task_id

    # 如果 state 中已经有路径，优先使用（但允许从查询中补充缺失的路径）
    has_data_path = bool(state.get("data_path"))
    has_result_path = bool(state.get("result_path"))

    default_result_path = os.path.join(".", "result", task_id).replace("\\", "/")

    # 从 user_query 中解析路径
    user_query = state.get("user_query", "")
    if user_query:
        parsed_paths = parse_paths_from_query(user_query)
        candidate_result_path = state.get("result_path") or parsed_paths.get("result_path") or default_result_path

        # 更新 data_path（如果 state 中没有）
        if not has_data_path and parsed_paths["data_path"]:
            requested_paths = [p.strip() for p in parsed_paths["data_path"].split(",") if p.strip()]
            missing_paths = []
            for raw_path in requested_paths:
                resolved = _resolve_input_file_host_path(
                    raw_path,
                    parsed_paths["data_path"],
                    candidate_result_path,
                )
                if resolved is None and not os.path.exists(raw_path):
                    missing_paths.append(raw_path)
            if not missing_paths:
                state["data_path"] = parsed_paths["data_path"]
                print(f"  --> 从用户查询中解析到数据路径: {parsed_paths['data_path']}")
            else:
                print(f"  --> 警告：解析到的数据路径部分或全部不存在: {', '.join(missing_paths)}")
                # 仍然设置路径，让后续代码处理
                state["data_path"] = parsed_paths["data_path"]
        elif not has_data_path:
            # 如果没有解析到且 state 中也没有，保持为空字符串
            state["data_path"] = ""

        # 更新 result_path（如果 state 中没有）
        if not has_result_path and parsed_paths["result_path"]:
            state["result_path"] = parsed_paths["result_path"]
            print(f"  --> 从用户查询中解析到结果路径: {parsed_paths['result_path']}")
        elif not has_result_path:
            # 如果没有解析到且 state 中也没有，使用默认值
            state["result_path"] = default_result_path
            print(f"  --> 为当前任务分配默认结果隔离路径: {state['result_path']}")

    # 确保 result_path 有默认值
    if not state.get("result_path"):
        state["result_path"] = default_result_path
        print(f"  --> 为当前任务分配默认结果隔离路径: {state['result_path']}")

    return state


def generate_code(state: CodeAgentState) -> CodeAgentState:
    """
    生成代码节点
    调用 LLM 生成代码，写入 pending_contribution
    """
    attempt_no = state.get('internal_iteration_count', 0) + 1
    generate_start = time.perf_counter()
    print(f"--- [Code Dev] 正在生成代码 (迭代 {attempt_no}) ---")

    critic_feedback = state.get("critique_feedback", "")
    internal_feedback = state.get("feedback", "")
    
    final_feedback = ""
    if critic_feedback:
        print(f"  --> 收到 Critic 的驳回意见: {critic_feedback[:500]}...")
        final_feedback = critic_feedback
    elif internal_feedback:
        print(f"  --> 收到内部执行的错误反馈: {internal_feedback[:500]}...")
        final_feedback = internal_feedback
    is_retry_attempt = bool(critic_feedback or internal_feedback or state.get("internal_iteration_count", 0) > 0)
    attempt_timing = _ensure_attempt_timing(state, attempt_no=attempt_no, is_retry=is_retry_attempt)
    attempt_timing["retry_reason"] = "critic_feedback" if critic_feedback else ("internal_feedback" if internal_feedback else "")

    # 首先从 user_query 中提取路径（如果 state 中没有）
    prompt_prep_start = time.perf_counter()
    state = extract_paths_from_state(state)

    # 构建 Prompt：如果有反馈，说明是修正模式
    context_instruction = ""
    previous_code = state.get("scanpy_code", "")
    previous_requirements = state.get("requirements_txt", "")

    if state.get("feedback"):
        context_instruction = f"""
        【重要！代码执行遇到错误，请务必修改代码或依赖】
        上一次执行遇到了以下报错：
        {final_feedback}

        上一次使用的代码如下：
        ```python
        {previous_code}
        ```
        上一次使用的 requirements.txt 如下：
        ```txt
        {previous_requirements}
        ```
        请仔细分析上述报错。如果报错指向缺少某个 Python 包（如 ModuleNotFoundError, ImportError 等），你**必须**在生成的 requirements.txt 代码块中补充该包。
        如果报错逻辑错误，请确保你输出的新的 python 代码块已经修复了该问题。
        请不要输出多余解释，务必包含新的 ```python 和 ```txt 块！
        """
    # 获取当前步骤的输入、预期输出和文件路径
    current_step_input = state.get('current_step_input', '')
    current_step_expected_output = state.get('current_step_expected_output', '')
    current_step_file_paths = state.get('current_step_file_paths', {})
    
    # 从计划中获取文件路径，如果存在则优先使用
    input_files = current_step_file_paths.get('input_files', []) if current_step_file_paths else []
    output_files = current_step_file_paths.get('output_files', []) if current_step_file_paths else []
    
    # 获取数据路径和结果路径
    # 优先使用计划中指定的输入文件，否则使用原有的data_path
    if input_files:
        # 使用计划中指定的第一个输入文件作为数据路径
        data_path = input_files[0]
    else:
        data_path = state.get('data_path', '')
    
    result_path = state.get('result_path', './result')
    
    # 如果计划中指定了输出文件路径，使用第一个作为参考路径
    if output_files:
        # 从输出文件路径中提取目录路径
        first_output = output_files[0]
        if os.path.isabs(first_output):
            result_path = os.path.dirname(first_output)
        else:
            result_path = os.path.dirname(first_output) if os.path.dirname(first_output) else result_path

    # 转换为 Docker 路径
    # 注意：如果 data_path 是文件，convert_to_docker_path 会返回 /app/data/filename.h5ad
    # 如果 data_path 是目录，会返回 /app/data
    if data_path:
        # 检查路径是否存在，如果不存在，尝试判断是文件还是目录
        if os.path.exists(data_path):
            docker_data_path = convert_to_docker_path(data_path, 'data')
        else:
            # 如果路径不存在，根据是否有扩展名判断（.h5ad 等通常是文件）
            if os.path.splitext(data_path)[1]:
                # 有扩展名，可能是文件
                filename = os.path.basename(data_path)
                docker_data_path = f"/app/data/{filename}"
            else:
                # 没有扩展名，可能是目录
                docker_data_path = "/app/data"
    else:
        docker_data_path = '/app/data'

    docker_output_path = '/app/output'

    skill_id = state.get("current_step_skill_id")
    resolved_env = resolve_environment_for_state(state, skill_id)
    # 统一使用较大的注入限制以兼容复杂 SKILL
    inj_limit = 12000
    skill_block = format_skill_injection_for_code_dev(skill_id, max_chars=inj_limit)
    env_package_versions_block = format_environment_package_versions_for_prompt(resolved_env)
    if env_package_versions_block:
        skill_block = (skill_block + "\n" + env_package_versions_block).strip()

    req_instruction = "requirements.txt 须且仅需列出代码中实际调用的第三方依赖环境包，若是标准库可留空。"
    if resolved_env:
        core_names = ", ".join(resolved_env.core_dependencies) if resolved_env.core_dependencies else "(none)"
        allowlist_names = ", ".join(resolved_env.allowed_runtime_extras) if resolved_env.allowed_runtime_extras else "(none)"
        req_instruction = (
            f"当前步骤已绑定预配置环境 profile={resolved_env.env_profile}。requirements.txt 仅允许声明运行时额外依赖；"
            f"禁止覆盖或重装核心依赖 {core_names}。允许的额外依赖白名单仅限：{allowlist_names}。"
            "若无需额外依赖，请保持 requirements.txt 为空。"
        )
    if state.get("global_requirements"):
        req_instruction = "Supervisor 已经提前为你生成并安装了该任务全局所需的 global requirements。因此，在此步骤生成代码时，请【留空 requirements.txt】！**除非**你明确在修复之前的环境报错（如 `ModuleNotFoundError`）并需要额外补充缺失包，才在 requirements.txt 中列出。"

    system_prompt = get_code_system_prompt(skill_block, docker_data_path, docker_output_path, req_instruction)

    # 构建任务描述，优先使用当前步骤的输入
    plan = state.get("plan", [])
    current_step_index = state.get("current_step_index", 0)
    if plan and 0 <= current_step_index < len(plan):
        step = plan[current_step_index]
        s_id = getattr(step, "step_id", current_step_index + 1)
        s_name = getattr(step, "name", "")
        s_desc = getattr(step, "description", "")
        s_acc = getattr(step, "acceptance_criteria", "")
        s_in = getattr(step, "input_files", [])
        s_out = getattr(step, "output_files", [])
        task_description = f"[步骤 {s_id}] {s_name}\n- 描述: {s_desc}\n- 验收: {s_acc}\n- 输入: {s_in}\n- 输出: {s_out}"
    else:
        task_description = current_step_input if current_step_input else state.get('task', state.get('user_query', ''))
    
    # 构建预期输出说明
    expected_output_note = ""
    if current_step_expected_output:
        expected_output_note = f"\n\n【预期输出要求】\n{current_step_expected_output}\n请确保生成的代码能够满足上述要求。"
    
    # 构建文件路径说明
    input_mount_plan = _build_input_mount_plan(data_path, result_path, input_files)
    canonical_mapping_lines = [
        f"- {item['input_ref']} -> {item['container_path']}"
        + (f" (host: {item['host_path']})" if item.get("host_path") else "")
        for item in input_mount_plan["mappings"]
    ]

    file_paths_note = ""
    if input_files:
        file_paths_note += f"\n【输入文件】\n" + "\n".join([f"- {f}" for f in input_files])
        if canonical_mapping_lines:
            file_paths_note += "\n【已验证的容器内输入路径映射】\n" + "\n".join(canonical_mapping_lines)
        file_paths_note += "\n（必须严格使用上面的容器内路径读取输入，禁止臆造其它文件名或路径。）"
        if input_mount_plan["missing_inputs"]:
            file_paths_note += "\n【警告】以下输入文件暂未在宿主机解析到，若继续使用它们，执行前挂载预检会直接失败：\n"
            file_paths_note += "\n".join([f"- {item}" for item in input_mount_plan["missing_inputs"]])
    if output_files:
        file_paths_note += f"\n【必须生成的输出文件】\n" + "\n".join([f"- {f}" for f in output_files])
        # 确保输出文件保存到指定路径
        docker_output_paths = [convert_to_docker_path(f, 'output') for f in output_files]
        file_paths_note += f"\n【Docker容器内输出路径】\n" + "\n".join([f"- {p}" for p in docker_output_paths])
    
    exploration_context = ""
    # 如果已存在探查结果，则强制传入供代码开发参考（不再局限于探查计划已全部完成后才传）
    if state.get("data_exploration_results"):
        res_str = "\n\n".join(state["data_exploration_results"])
        exploration_context = f"\n【前期数据探查参考】\n已获取的数据统计或结构特征如下：\n{res_str}\n请在编写代码时充分考虑上述数据特征（如针对真实列名、数据格式、特定字段等进行操作）。\n"

    historical_outputs_context = ""
    completed_outputs = state.get("completed_steps_outputs", [])
    if completed_outputs:
        # 取最近的3-4个步骤防止溢出
        recent_outputs = completed_outputs[-4:]
        outputs_str = "\n\n".join(recent_outputs)
        historical_outputs_context = f"\n【前期分析步骤的执行上下文摘要】\n以下是前面环节产生的日志与输出产物报告。如果你的代码需要载入前面环节生成的中间文件（如结果表或对象数据），请务必要参考其中的输出路径以保证挂载文件获取正确：\n{outputs_str}\n"

    user_prompt = get_code_user_prompt(historical_outputs_context, exploration_context, task_description, expected_output_note, file_paths_note, context_instruction)

    # 替换 prompt 中的本地路径为 Docker 路径，防止大模型嵌套创建目录
    if result_path and result_path != './result':
        base_result = os.path.basename(result_path)
        if base_result:
            # 用户 Prompt
            user_prompt = user_prompt.replace(f"./result/{base_result}", "/app/output")
            user_prompt = user_prompt.replace(f"result/{base_result}", "/app/output")
            user_prompt = user_prompt.replace(f"/app/output/{base_result}", "/app/output")
            # 系统 Prompt
            system_prompt = system_prompt.replace(f"./result/{base_result}", "/app/output")
            system_prompt = system_prompt.replace(f"result/{base_result}", "/app/output")
            system_prompt = system_prompt.replace(f"/app/output/{base_result}", "/app/output")
    prompt_prep_elapsed = time.perf_counter() - prompt_prep_start

    try:
        sandbox_bundle = build_code_dev_sandbox_tool_bundle(state)
        state["sandbox_backend"] = sandbox_bundle.backend
        state["sandbox_endpoint"] = sandbox_bundle.endpoint
        state["sandbox_allowed_roots"] = sandbox_bundle.allowed_roots
        llm_with_tools = llm.bind_tools(sandbox_bundle.tools)
        
        # =============== 新增：工具调用专用的信息收集阶段 ===============
        is_exploration = not state.get("data_exploration_done", True)
        is_retry = state.get("internal_iteration_count", 0) > 0
        tool_context_str = state.get("step_tool_context", "")
        tool_loop_elapsed = 0.0
        tool_llm_elapsed = 0.0
        tool_io_elapsed = 0.0
        tool_iterations = 0
        tool_calls_count = 0
        tool_stats: dict[str, dict[str, float | int]] = {}
        reused_tool_context = False

        if is_exploration:
            print(f"🔍 [Code Dev] 当前处于数据探查阶段，直接生成探查代码（跳过 SKILL 工具调研）...")
        elif is_retry and tool_context_str:
            print(f"🔄 [Code Dev] 当前步骤为重试迭代，使用已有的 Tool Context 进行重试，不再重新调研...")
            reused_tool_context = True
        elif not skill_id:
            print(f"⏩ [Code Dev] Supervisor 未指派具体 SKILL (skill_id 为空)，跳过工具调研阶段直接生成代码...")
        else:
            tool_loop_start = time.perf_counter()
            print("\n" + "="*60)
            print(f"🔍 [Code Dev Tool Loop] 开始为执行任务查阅 SKILL 目录...")
            print(
                f"🔐 [Code Dev Tool Loop] sandbox backend={sandbox_bundle.backend} "
                f"allowed_roots={sandbox_bundle.allowed_roots}"
            )
            wf_host_path = resolve_workflow_root(skill_id) if skill_id else ""
            
            tool_system_prompt = get_tool_system_prompt(skill_id, wf_host_path)
            tool_messages = [
                SystemMessage(content=tool_system_prompt),
                HumanMessage(content=get_tool_user_prompt(task_description))
            ]
            
            tool_iterations = 0
            max_tool_iterations = 8
            collected_info = []

            while tool_iterations < max_tool_iterations:
                tool_llm_start = time.perf_counter()
                response = llm_with_tools.invoke(tool_messages)
                tool_llm_elapsed += time.perf_counter() - tool_llm_start
                tool_messages.append(response)
                
                if not response.tool_calls:
                    break

                print(f"  [Code Dev Tool Loop] 正在执行第 {tool_iterations+1} 轮工具调用...")
                for tool_call in response.tool_calls:
                    tool_calls_count += 1
                    print(f"  [Code Dev Tool Loop] 调用工具: {tool_call['name']}({tool_call['args']})")
                    tool_call_start = time.perf_counter()
                    try:
                        tool_instance = sandbox_bundle.tool_map.get(tool_call["name"])
                        if tool_instance is None:
                            tool_result = f"Error: Unknown tool {tool_call['name']}"
                        else:
                            tool_result = tool_instance.invoke(tool_call["args"])
                    except Exception as e:
                        tool_result = f"Error: 工具执行异常: {str(e)}"
                        print(f"  [Code Dev Tool Loop] 异常: {str(e)}")
                    tool_call_elapsed = time.perf_counter() - tool_call_start
                    tool_io_elapsed += tool_call_elapsed
                    stat = tool_stats.setdefault(tool_call["name"], {"count": 0, "elapsed_seconds": 0.0})
                    stat["count"] = int(stat.get("count", 0)) + 1
                    stat["elapsed_seconds"] = float(stat.get("elapsed_seconds", 0.0)) + tool_call_elapsed
                    
                    preview = str(tool_result)[:500].replace('\n', ' ')
                    print(f"  [Code Dev Tool Loop] 结果预览: {preview}...")
                    
                    tool_msg = ToolMessage(content=str(tool_result), tool_call_id=tool_call["id"], name=tool_call["name"])
                    tool_messages.append(tool_msg)
                    
                    # 截取一部分结果供后续 Code 生成使用（避免超长提示词卡死）
                    collected_info.append(f"---- 工具 {tool_call['name']} 目标 {tool_call['args']} 返回：\n{str(tool_result)[:3000]}")
                
                tool_iterations += 1
                
            if collected_info:
                tool_context_str = "\n【来自工具调研获取的最佳背景内容】\n" + "\n".join(collected_info) + "\n\n请结合上方参考资料和 SKILL 要求，完成任务代码编写。\n"
                
                # 兼容不同系统路径分隔符及 Python dict/JSON 序列化后的转义反斜杠
                host_wf_root_os = get_workflows_root()
                host_wf_root_posix = host_wf_root_os.replace('\\', '/')
                host_wf_root_escaped = host_wf_root_os.replace('\\', '\\\\')
                
                # 依次替换转义格式、原生 OS 格式、POSIX 格式
                tool_context_str = tool_context_str.replace(host_wf_root_escaped, "/app/workflows")
                tool_context_str = tool_context_str.replace(host_wf_root_os, "/app/workflows")
                tool_context_str = tool_context_str.replace(host_wf_root_posix, "/app/workflows")

                tool_context_str = tool_context_str.replace('\\', '/')
                # 修复可能出现的双斜杠，但不影响 URL 中的 ://
                tool_context_str = re.sub(r'(?<!:)//+', '/', tool_context_str)
                state["step_tool_context"] = tool_context_str

            print(f"🔍 [Code Dev Tool Loop] 调研结束，进入正式代码生成...")
            tool_loop_elapsed = time.perf_counter() - tool_loop_start
            print(
                "【Code Dev Tool Loop 耗时】"
                f"total={_fmt_elapsed(tool_loop_elapsed)}, "
                f"llm={_fmt_elapsed(tool_llm_elapsed)}, "
                f"tool_io={_fmt_elapsed(tool_io_elapsed)}, "
                f"iterations={tool_iterations}, "
                f"tool_calls={tool_calls_count}, "
                f"stats={_format_tool_stats(tool_stats)}"
            )
            print("="*60 + "\n")
        # =============== 信息收集阶段结束 ===============

        # 将工具内容合并进 user_prompt 作为上下文，并添加 SKILL_BLOCK 要求
        final_user_prompt = get_final_user_prompt(tool_context_str, skill_block, user_prompt)

        messages = [
            SystemMessage(content=system_prompt),
            HumanMessage(content=final_user_prompt)
        ]

        # 打印最终提示词
        print("\n" + "*"*60)
        print(f"🔍 [Code Dev Debug] 发送给 LLM 的代码生成提示词 (迭代 {state.get('internal_iteration_count', 0) + 1})")
        print("*"*60)
        print("【System Prompt】:")
        print(system_prompt.strip())
        print("\n【User Prompt】:")
        print(final_user_prompt.strip())
        print("*"*60 + "\n")

        # 生成正式代码
        llm_generate_start = time.perf_counter()
        response = llm.invoke(messages)
        llm_generate_elapsed = time.perf_counter() - llm_generate_start
        text = response.content

        # 提取代码块和 requirements
        parse_start = time.perf_counter()
        # requirements 支持多种标签，避免模型输出轻微变化导致解析失败
        # 允许 ```python / ```py、大小写及首尾空白；非贪婪匹配至第一个闭合 ```
        python_pattern = r"```(?:python|py)\s*\n(.*?)```"
        requirements_patterns = [
            r'```requirements(?:\.txt)?\n(.*?)\n```',
            r'```txt\n(.*?)\n```',
            r'```text\n(.*?)\n```',
        ]

        python_match = re.search(python_pattern, text, re.DOTALL | re.IGNORECASE)
        requirements_match = None
        for pattern in requirements_patterns:
            m = re.search(pattern, text, re.DOTALL | re.IGNORECASE)
            if m:
                requirements_match = m
                break

        # 兜底：允许“requirements.txt:”后直接给列表（无代码块）
        if not requirements_match:
            plain_req_match = re.search(
                r'requirements(?:\.txt)?\s*[:：]\s*\n([\s\S]+)$',
                text,
                re.IGNORECASE,
            )
            if plain_req_match:
                plain_req_text = plain_req_match.group(1).strip()
                if plain_req_text:
                    requirements_match = plain_req_match

        if python_match:
            code = sanitize_llm_python_block(python_match.group(1).strip())
        else:
            # 如果没有代码块，尝试提取整个响应
            print("没有获取到 code 代码块，作为备用提取整个响应")
            code = sanitize_llm_python_block(text.strip())

        if requirements_match:
            requirements = requirements_match.group(1).strip()
        else:
            print("未获取到 requirements.txt 代码块，使用默认值")
            # 默认 requirements 置空
            requirements = ""
        parse_elapsed = time.perf_counter() - parse_start

        # 更新状态
        state["scanpy_code"] = code
        state["requirements_txt"] = requirements
        state["internal_iteration_count"] = state.get("internal_iteration_count", 0) + 1

        # 将代码写入 pending_contribution
        state["pending_contribution"] = {
            "code": code,
            "requirements": requirements,
            "task": state.get("task", "")
        }

        attempt_timing["generate_code"] = {
            "prompt_prep_elapsed_seconds": round(prompt_prep_elapsed, 3),
            "tool_loop_elapsed_seconds": round(tool_loop_elapsed, 3),
            "tool_llm_elapsed_seconds": round(tool_llm_elapsed, 3),
            "tool_io_elapsed_seconds": round(tool_io_elapsed, 3),
            "llm_generate_elapsed_seconds": round(llm_generate_elapsed, 3),
            "parse_elapsed_seconds": round(parse_elapsed, 3),
            "total_elapsed_seconds": round(time.perf_counter() - generate_start, 3),
            "tool_iterations": tool_iterations,
            "tool_calls": tool_calls_count,
            "tool_stats": {
                name: {
                    "count": int(info.get("count", 0)),
                    "elapsed_seconds": round(float(info.get("elapsed_seconds", 0.0)), 3),
                }
                for name, info in tool_stats.items()
            },
            "reused_tool_context": reused_tool_context,
        }
        print(f"【Code Dev生成耗时】{_format_phase_timing(attempt_timing)}")
        print(f"  --> 代码生成成功，代码长度: {len(code)} 字符")

    except Exception as e:
        error_msg = f"代码生成失败: {str(e)}"
        print(f"  --> {error_msg}")
        state["scanpy_code"] = f"代码生成失败: {e}"
        state["requirements_txt"] = f"requirements.txt 生成失败：{str(e)}"
        state["pending_contribution"] = {"error": error_msg}
        attempt_timing["generate_code"] = {
            "prompt_prep_elapsed_seconds": round(prompt_prep_elapsed, 3),
            "total_elapsed_seconds": round(time.perf_counter() - generate_start, 3),
            "failed": True,
        }
        print(f"【Code Dev生成耗时】{_format_phase_timing(attempt_timing)}")
        print(f"模型调用失败：{e}")

    return state


def self_reflection(state: CodeAgentState) -> CodeAgentState:
    """
    自我检查节点
    简单的代码质量检查
    """
    code = state.get("scanpy_code", "")
    if not code or code.startswith("# Error"):
        return state

    reflection_start = time.perf_counter()
    print("--- [Code Dev] 进行自我检查 ---")

    # 简单的安全检查
    dangerous_patterns = ["eval(", "exec(", "__import__", "open("]
    warnings = []
    for pattern in dangerous_patterns:
        if pattern in code:
            warnings.append(f"检测到潜在风险: {pattern}")

    if warnings:
        print(f"  --> 警告: {', '.join(warnings)}")

    attempt_timing = state.get("current_attempt_timing")
    if isinstance(attempt_timing, dict):
        attempt_timing["self_reflection"] = {
            "warnings_count": len(warnings),
            "total_elapsed_seconds": round(time.perf_counter() - reflection_start, 3),
        }
        print(f"【Code Dev自检耗时】attempt={attempt_timing.get('attempt', '?')}, total={_fmt_elapsed(attempt_timing['self_reflection']['total_elapsed_seconds'])}, warnings={len(warnings)}")

    return state


def execute_code(state: CodeAgentState) -> CodeAgentState:
    """
    执行代码节点
    在 Docker 容器中执行生成的代码
    """
    execute_start = time.perf_counter()
    prep_start = execute_start
    print("--- [Code Dev] 正在执行代码 ---")

    code = state.get("scanpy_code", "")
    requirements = state.get("requirements_txt", "")
    attempt_timing = state.get("current_attempt_timing")

    if not code or code.startswith("# Error"):
        state["success"] = False
        state["analysis_result"] = "代码生成失败，无法执行"
        return state

    # 确保结果目录存在
    result_path = state.get("result_path", "./result")
    
    # 如果计划中指定了输出文件路径，确保对应的目录存在
    current_step_file_paths = state.get("current_step_file_paths", {})
    output_files = current_step_file_paths.get("output_files", []) if current_step_file_paths else []
    if output_files:
        for output_file in output_files:
            output_dir = os.path.dirname(output_file) if os.path.dirname(output_file) else result_path
            os.makedirs(output_dir, exist_ok=True)
    
    os.makedirs(result_path, exist_ok=True)

    # 构建完整的可执行代码（参考 umap_langgraph.py 的改进）
    _skill_for_exec = state.get("current_step_skill_id")
    
    _skill_path_docker = f"/app/workflows/{_skill_for_exec}" if _skill_for_exec else "/app/workflows"

    header = f"""
# 基础库导入（确保代码独立运行）
import sys
import os
sys.path.append(os.getcwd())
if '{_skill_path_docker}' not in sys.path:
    sys.path.insert(0, '{_skill_path_docker}')
# --- DEBUG START: 检查挂载情况 ---
print("DEBUG: Checking /app/data contents...")
try:
    if os.path.exists('/app/data'):
        print(f"DEBUG: Files in /app/data: {{os.listdir('/app/data')}}")
    else:
        print("DEBUG: /app/data does not exist!")
    
    if os.path.exists('/app/output'):
        print(f"DEBUG: Files in /app/output: {{os.listdir('/app/output')}}")
    if os.path.exists('{_skill_path_docker}'):
        print(f"DEBUG: Files in {_skill_path_docker}: {{os.listdir('{_skill_path_docker}')}}")
except Exception as e:
    print(f"DEBUG: Error checking directories: {{e}}")
# --- DEBUG END ---
print("===MAS_EXEC_START===")
"""
# 核心分析代码（来自大模型生成）
    llm_code = state.get("scanpy_code", "")

    footer = f"""
import os

# --- 补充输出文件信息 ---
try:
    if 'analysis_summary' in locals():
        print(f"===RESULT==={{analysis_summary}}===")
    else:
        out_files = os.listdir('/app/output') if os.path.exists('/app/output') else []
        if out_files:
            print(f"===RESULT===Step executed. Generated files: {{', '.join(out_files)}}===")
        else:
            print(f"===RESULT===Step executed (No specific result string).===")
except Exception as e:
    print(f"===RESULT===Execution finished, but result extraction failed: {{str(e)}}===")
"""
    
        

    full_code = header + "\n# --- LLM Generated Code ---\n" + llm_code + "\n" + footer

    data_path = state.get("data_path", "")
    current_step_file_paths = state.get("current_step_file_paths", {})
    input_files = current_step_file_paths.get("input_files", []) if current_step_file_paths else []
    input_mount_plan = _build_input_mount_plan(data_path, result_path, input_files)
    actual_data_dirs = list(input_mount_plan["data_dirs"])
    actual_data_path = actual_data_dirs[0] if actual_data_dirs else data_path

    for mapping in input_mount_plan["mappings"]:
        if mapping.get("host_path"):
            print(
                f"  --> 输入映射: {mapping['input_ref']} -> {mapping['container_path']} "
                f"(host: {mapping['host_path']})"
            )
    for missing in input_mount_plan["missing_inputs"]:
        print(f"  --> 警告：未在宿主机解析到输入文件: {missing}")

    if not actual_data_dirs and (not actual_data_path or not os.path.exists(actual_data_path)):
        for default_data_dir in _mas2_data_dir_candidates():
            if default_data_dir and os.path.isdir(default_data_dir):
                actual_data_path = default_data_dir
                actual_data_dirs = [default_data_dir]
                print(f"  --> 数据路径无效，使用默认 data 目录: {actual_data_path}")
                break

    # 在 Docker 容器中执行代码
    # 传递 input_files 以便 executor 智能确定需要挂载的目录
    from src.utils.workflow_skills import get_workflows_root
    workflow_host = get_workflows_root()
    resolved_env = resolve_environment_for_state(state, _skill_for_exec)
    executor_container_id = state.get("docker_container_id")
    asset_cache_host_path = None

    if resolved_env:
        print(
            f"ENV RESOLVE profile={resolved_env.env_profile} "
            f"runtime={resolved_env.runtime} reason={resolved_env.skill_id}"
        )
        prev_sig = state.get("docker_env_signature")
        prev_profile = state.get("docker_env_profile")
        if prev_sig and prev_sig != resolved_env.env_signature:
            print(
                f"ENV SWITCH from={prev_profile or 'legacy'} to={resolved_env.env_profile}"
            )
            executor_container_id = None
        for asset_name, exists, asset_path in asset_cache_status(resolved_env):
            print(
                f"ASSET CACHE {'HIT' if exists else 'MISS'} "
                f"asset={asset_name} path={asset_path}"
            )
        asset_cache_root = ensure_asset_cache(resolved_env)
        if asset_cache_root is not None:
            asset_cache_host_path = str(asset_cache_root)

        if not resolved_env.is_python:
            err = (
                f"当前 workflow `{_skill_for_exec}` 绑定 runtime={resolved_env.runtime} / "
                f"image={resolved_env.env_image}，但当前 Code Dev 执行器仅支持 Python。"
            )
            print(f"  --> [ENV UNSUPPORTED] {err}")
            state["success"] = False
            state["analysis_result"] = err
            state["pending_contribution"] = _build_execute_pending_contribution(
                code=full_code,
                requirements=requirements,
                task=state.get("task", ""),
                output_str="",
                success=False,
                error_msg=err,
            )
            return state

        filtered_requirements = filter_profiled_requirements(resolved_env, requirements)
        if filtered_requirements.blocked_lines:
            print(
                "  --> [PROFILE REQUIREMENTS FILTER] 已忽略以下 requirements: "
                + ", ".join(filtered_requirements.blocked_lines)
            )
        if filtered_requirements.blocked_core_lines:
            print(
                "  --> [PROFILE REQUIREMENTS FILTER] 已阻止核心依赖覆盖: "
                + ", ".join(filtered_requirements.blocked_core_lines)
            )
        sync_environment_to_state(
            state,
            resolved_env,
            env_cache_key=filtered_requirements.env_cache_key,
            env_cache_volume=filtered_requirements.env_cache_volume,
        )

        if workflow_host and not executor_container_id:
            print(f"  --> Workflow 目录将挂载到容器 /app/workflows: {workflow_host}")

        shared_cache_volume = (
            filtered_requirements.env_cache_volume
            if resolved_env.is_shared
            else None
        )
        executor = CodeExecutor(
            docker_path=None,
            container_id=executor_container_id,
            data_dir=actual_data_path if actual_data_path and os.path.exists(actual_data_path) else None,
            data_dirs=actual_data_dirs if actual_data_dirs else None,
            input_files=input_files if input_files else None,
            output_dir=result_path,
            workflow_host_path=workflow_host,
            docker_image=resolved_env.env_image,
            runtime=resolved_env.runtime,
            env_profile=resolved_env.env_profile,
            env_signature=resolved_env.env_signature,
            env_mode=resolved_env.env_mode,
            asset_cache_host_path=asset_cache_host_path,
            env_cache_key=filtered_requirements.env_cache_key,
            env_cache_volume=shared_cache_volume,
            required_input_container_paths=input_mount_plan["required_container_paths"],
        )

        if resolved_env.is_shared:
            cache_result = executor.ensure_env_cache_ready(
                requirements_str=filtered_requirements.allowed_text,
                timeout=900,
            )
            cache_timing = cache_result.get("timing", {}) if isinstance(cache_result, dict) else {}
            if cache_timing:
                print(
                    "  --> [环境缓存预热] "
                    f"elapsed={float(cache_timing.get('cache_prepare_elapsed_seconds', 0.0)):.2f}s, "
                    f"cache_hit={bool(cache_timing.get('cache_hit', False))}, "
                    f"install={bool(cache_timing.get('cache_install_performed', False))}"
                )
            if not cache_result.get("success"):
                err = cache_result.get("error", "环境缓存准备失败")
                print(f"  --> [ENV CACHE FAILED] {err}")
                state["success"] = False
                state["analysis_result"] = err
                state["pending_contribution"] = _build_execute_pending_contribution(
                    code=full_code,
                    requirements=filtered_requirements.allowed_text,
                    task=state.get("task", ""),
                    output_str=cache_result.get("output", "") if isinstance(cache_result, dict) else "",
                    success=False,
                    error_msg=err,
                )
                return state
            requirements = ""
        else:
            requirements = filtered_requirements.allowed_text
    else:
        if state.get("docker_env_signature"):
            print(
                f"ENV SWITCH from={state.get('docker_env_profile') or 'profiled'} to=legacy"
            )
            clear_environment_from_state(state)
        if workflow_host and not state.get("docker_container_id"):
            print(f"  --> Workflow 目录将挂载到容器 /app/workflows: {workflow_host}")

        executor = CodeExecutor(
            docker_path=None,
            container_id=None if state.get("docker_env_signature") else state.get('docker_container_id'),
            data_dir=actual_data_path if actual_data_path and os.path.exists(actual_data_path) else None,
            data_dirs=actual_data_dirs if actual_data_dirs else None,
            input_files=input_files if input_files else None,
            output_dir=result_path,
            workflow_host_path=workflow_host,
            required_input_container_paths=input_mount_plan["required_container_paths"],
        )

        # 仅 legacy/unprofiled workflow 使用 global_requirements fallback
        if not state.get('docker_container_id') and state.get('global_requirements'):
            print(f"  --> [LEGACY ENV] 首次执行，组合并预安装 global_requirements")
            global_reqs = state['global_requirements'].strip()
            if global_reqs:
                if requirements.strip():
                    requirements = global_reqs + "\n" + requirements.strip()
                else:
                    requirements = global_reqs

    print("\n" + "="*50)
    print(f"🚀 [Code Dev] 准备推送到 Docker 执行的代码与环境 (迭代 {state.get('internal_iteration_count', 0)})")
    print("="*50)
    print("【requirements.txt】")
    print(requirements.strip() if requirements.strip() else "(空)")
    print("-" * 50)
    print("【code.py】")
    print(full_code.strip())
    print("="*50 + "\n")
    prep_elapsed = time.perf_counter() - prep_start

    try:
        step_name = state.get("current_step_input") or _skill_for_exec or f"step_{state.get('current_step_index', 0)}"
        print(f"EXEC START step={step_name}")
        executor_start = time.perf_counter()
        result = executor.execute(code_str=full_code, requirements_str=requirements, timeout=600)
        executor_call_elapsed = time.perf_counter() - executor_start
        print(f"EXEC END step={step_name} elapsed={executor_call_elapsed:.2f}s")
        exec_timing = result.get("timing") if isinstance(result, dict) else None

        if result.get('container_id'):
            new_cid = result['container_id']
            if new_cid != state.get('docker_container_id'):
                print(f"  --> [容器建立成功] 已成功构建持久化容器，容器ID: {new_cid[:12]}")
            state['docker_container_id'] = new_cid
            if resolved_env:
                sync_environment_to_state(state, resolved_env, container_id=new_cid)

        if exec_timing:
            print(f"【Docker阶段耗时】{_format_executor_timing(exec_timing)}")
            if resolved_env and not exec_timing.get("pip_install_skipped", True):
                print(
                    f"ENV EXTRA INSTALL profile={resolved_env.env_profile} "
                    f"elapsed={float(exec_timing.get('pip_install_elapsed_seconds', 0.0)):.2f}s"
                )

        # 控制台仅输出关键摘要，避免安装依赖等噪声日志淹没真实错误。
        result_parse_start = time.perf_counter()
        output_str = result.get('output', '')
        output_str = output_str.split("===MAS_EXEC_START===")[-1] if "===MAS_EXEC_START===" in output_str else output_str

        # 检查执行是否成功（从executor返回的success字段）
        executor_success = result.get('success', True)
        
        # 检查output中是否包含错误信息（即使executor返回success=True，代码执行也可能失败）
        has_error_in_output = any(keyword in output_str for keyword in [
            'Traceback (most recent call last):', 'Error:', 'Exception:', 'TypeError:', 'ValueError:',
            'AttributeError:', 'NameError:', 'KeyError:', 'IndexError:',
            'SyntaxError:', 'IndentationError:', 'AssertionError:'
        ])
        result_part = ""
        result_looks_failed = False
        if "===RESULT===" in output_str:
            parts = output_str.rsplit("===RESULT===", 1)
            result_part = parts[-1]
            if result_part.endswith("==="):
                result_part = result_part[:-3]
            elif "===" in result_part:
                result_part = result_part.rsplit("===", 1)[0]
                
            # 去除 pip 等安装包的前置外部日志，只保留真实执行部分
            stdout_part = parts[0].strip()
                
            if stdout_part:
                clean_lines = []
                for line in stdout_part.split('\n'):
                    if line.startswith("DEBUG: ") or line.startswith("# --- DEBUG ") or line.strip() == "":
                        continue
                    clean_lines.append(line)
                clean_stdout = "\n".join(clean_lines).strip()
                
                if clean_stdout:
                    if result_part.strip() in ["Step completed successfully (No specific return value).", "Step executed (No specific result string).", "Step executed.", ""]:
                        result_part = clean_stdout
                    else:
                        result_part = f"{clean_stdout}\n{result_part}"

            _res_low = result_part.strip().lower()
            result_looks_failed = any(
                hint in _res_low
                for hint in (
                    "analysis failed",
                    "execution failed",
                    "fatal error",
                    "traceback (most recent call last)",
                    "unable to execute",
                    "执行失败",
                )
            )

        if executor_success and "===RESULT===" in output_str and not has_error_in_output and not result_looks_failed:
            _succ = (result_part or "Execution successful").strip()
            print(f"【Docker执行摘要】SUCCESS: {_succ}")
        else:
            _err_summary = str(result.get('error', '') or "").strip() or _extract_informative_error(output_str)
            if not _err_summary:
                _err_summary = "代码执行失败，但未提取到明确错误"
            print(f"【Docker执行摘要】FAILED: {_err_summary}")
        result_parse_elapsed = time.perf_counter() - result_parse_start

        if isinstance(attempt_timing, dict):
            attempt_timing["execute_code"] = {
                "prep_elapsed_seconds": round(prep_elapsed, 3),
                "executor_call_elapsed_seconds": round(executor_call_elapsed, 3),
                "result_parse_elapsed_seconds": round(result_parse_elapsed, 3),
                "total_elapsed_seconds": round(time.perf_counter() - execute_start, 3),
            }
            print(f"【Code Dev执行耗时】{_format_phase_timing(attempt_timing)}")
        
        # 提取结果（参考 umap_langgraph.py 的改进）
        if executor_success and "===RESULT===" in output_str and not has_error_in_output and not result_looks_failed:
            # 提取结果部分
            state["analysis_result"] = result_part.strip()
            state["success"] = True  # 标记为成功
            print("  --> 代码执行成功！已提取分析结果")

            state["pending_contribution"] = _build_execute_pending_contribution(
                code=code,
                requirements=requirements,
                task=state.get("task", ""),
                output_str=output_str,
                success=True,
                result_value=state["analysis_result"],
                output_files=result.get("files", []),
                timing=exec_timing,
                phase_timing=attempt_timing if isinstance(attempt_timing, dict) else None,
            )
        else:
            # 执行失败或没有找到结果标记
            # 优先使用executor返回的error字段，否则从output中提取错误信息
            error_msg = result.get('error', '')
            
            # 如果没有error字段，尝试从output中提取错误信息
            if not error_msg and output_str:
                # 优先提取完整的Traceback信息
                if 'Traceback' in output_str:
                    # 提取从Traceback到最后一个错误行的内容
                    traceback_match = re.search(
                        r'(Traceback \(most recent call last\):.*?)(?=\n\n|\n[A-Z][a-z]+:|\Z)',
                        output_str,
                        re.DOTALL
                    )
                    if traceback_match:
                        error_msg = traceback_match.group(1).strip()
                    else:
                        # 如果没找到完整traceback，提取包含错误的行
                        lines = output_str.split('\n')
                        error_lines = []
                        for i, line in enumerate(lines):
                            if any(keyword in line for keyword in [
                                'TypeError', 'ValueError', 'AttributeError', 
                                'NameError', 'KeyError', 'IndexError', 'Error:'
                            ]):
                                # 包含这一行和前面几行上下文
                                start = max(0, i - 3)
                                error_lines = lines[start:i+1]
                                break
                        if error_lines:
                            error_msg = '\n'.join(error_lines).strip()
                
                # 如果还是没找到，使用完整输出作为错误信息（便于排查）
                if not error_msg:
                    error_msg = _extract_informative_error(output_str)
            
            # 如果仍然没有错误信息，使用默认值
            if not error_msg:
                if result_looks_failed and result_part.strip():
                    error_msg = result_part.strip()
                else:
                    error_msg = '代码执行失败，但未找到具体错误信息'
            
            # 构建错误信息
            if "===RESULT===" not in output_str:
                state["analysis_result"] = f"代码执行失败\\n错误信息：{error_msg}"
            else:
                state["analysis_result"] = f"代码执行完成，但未找到结果标记\\n错误日志：{error_msg}"
            
            state["success"] = False
            _fail_preview = error_msg if len(error_msg) <= 500 else error_msg[:500] + "..."
            print(f"  --> 代码执行失败: {_fail_preview}")

            state["pending_contribution"] = _build_execute_pending_contribution(
                code=code,
                requirements=requirements,
                task=state.get("task", ""),
                output_str=output_str or "",
                success=False,
                error_msg=error_msg,
                timing=exec_timing,
                phase_timing=attempt_timing if isinstance(attempt_timing, dict) else None,
            )

    except Exception as e:
        # 处理其他运行时错误（参考 umap_langgraph.py 的改进）
        error_msg = f"Docker代码运行失败：{str(e)}"
        _r = locals().get("result")
        _exec_log = (_r.get("output") or "") if isinstance(_r, dict) else ""
        concise_log = _extract_informative_error(_exec_log) if _exec_log else "无"
        state["analysis_result"] = f"{error_msg}\\n错误日志：{concise_log}"
        state["success"] = False
        print(f"  --> {error_msg}")

        state["pending_contribution"] = _build_execute_pending_contribution(
            code=code,
            requirements=requirements,
            task=state.get("task", ""),
            output_str=_exec_log or "",
            success=False,
            error_msg=error_msg,
            timing=_r.get("timing") if isinstance(_r, dict) else None,
            phase_timing=attempt_timing if isinstance(attempt_timing, dict) else None,
        )

    return state

def display_result(state: CodeAgentState) -> CodeAgentState:
    """
    功能：展示分析结果文本+UMAP聚类图（增加容错，避免解码崩溃）
    """
    is_exploration = not state.get("data_exploration_done", True)

    if state.get("success", False):
        if is_exploration:
            print("\n=== [Code Dev] 探查执行结果 ===")
            print("-" * 30)
            print(str(state.get("analysis_result", "")).strip())
            print("-" * 30)
        else:
            print("\n=== [Code Dev] 展现分析结果 ===")
            # 显示文本结果（优先保证文本能看到）
            print("单细胞分析结果：")
            print("-"*30)
            print(str(state.get("analysis_result", "")).strip())
            
            # 显示所有PNG图片（增加容错）
            print("正在处理PNG图片并创建对应的HTML文件：")
            print("-"*30)

            import os
            
            # Ensure result_path exists before accessing
            if "result_path" in state and os.path.exists(state["result_path"]):
                # Get all PNG files in the result directory
                result_dir = state['result_path']
                png_files = [f for f in os.listdir(result_dir) if f.lower().endswith('.png')]

                # Check if we have write permissions to the result directory
                if not os.access(result_dir, os.W_OK):
                    # If no write access to result directory, create a single temporary directory for all outputs
                    with tempfile.TemporaryDirectory() as temp_dir:
                        for png_file in png_files:
                            png_path = os.path.join(result_dir, png_file)

                            # Generate corresponding HTML filename (e.g., leiden.png -> leiden_decoded.html)
                            base_name = os.path.splitext(png_file)[0]
                            html_filename = f"{base_name}_decoded.html"
                            output_html_path = os.path.join(temp_dir, html_filename)

                            # Call the function to create HTML with base64 image
                            create_html_with_base64_image(png_path, output_html_path)

                            # Inform user of location
                            print(f"HTML file saved at: {output_html_path}")
                            print(f"Please copy the file manually to result directory if needed.")
                else:
                    # We have write access, proceed normally for each PNG file
                    for png_file in png_files:
                        png_path = os.path.join(result_dir, png_file)

                        # Generate corresponding HTML filename (e.g., leiden.png -> leiden_decoded.html)
                        base_name = os.path.splitext(png_file)[0]
                        html_filename = f"{base_name}_decoded.html"
                        output_html_path = os.path.join(result_dir, html_filename)

                        create_html_with_base64_image(png_path, output_html_path)
                        print(f"Created HTML file: {html_filename} for {png_file}")
    else:
        # 超过迭代次数仍失败，给出准确精简的向用户解释的错误原因
        print("\n=== [Code Dev] 任务执行失败 ===")
        print("运行失败详情：")
        print("-"*30)
        error_info = str(state.get("analysis_result", "未知错误"))[:1000]
        explanation = f"经过 {state.get('internal_iteration_count', 0)} 次尝试，代码仍然无法成功执行。\n主要原因是:\n{error_info}\n请检查数据路径、文件格式或依赖是否正确。"
        print(explanation)
        
        # 将结果覆盖为对用户的精简解释
        state["analysis_result"] = explanation

    # 保留 Docker 容器，生命周期由 main.py 中 finalize_step 管理

    state["last_worker"] = "code_dev"
    return state


def should_retry(state: CodeAgentState) -> str:
    """
    判断是否应该重试代码生成
    如果执行失败且没有达到最大重试次数，则重试
    """
    max_retries = 3
    iteration_count = state.get("internal_iteration_count", 0)
    success = state.get("success", False)

    if success:
        print(f"【Code Dev重试决策】attempt={iteration_count}, success=True, next=end")
        return "end"
    elif iteration_count < max_retries:
        print(f"  --> 执行失败，准备重试 (第 {iteration_count + 1}/{max_retries} 次)")
        print(f"【Code Dev重试决策】attempt={iteration_count}, success=False, next=retry, remaining={max_retries - iteration_count}")
        return "retry"
    else:
        print(f"  --> 已达到最大重试次数 ({max_retries})，停止重试")
        print(f"【Code Dev重试决策】attempt={iteration_count}, success=False, next=end")
        return "end"


def prepare_retry(state: CodeAgentState) -> CodeAgentState:
    """
    准备重试节点
    将错误信息作为反馈，用于下一次代码生成
    """
    if not state.get("success", False):
        retry_start = time.perf_counter()
        pending = state.get("pending_contribution") or {}
        feedback = ""

        if isinstance(pending, dict):
            feedback = _extract_informative_error(str(pending.get("error", "") or ""))
            if not feedback:
                feedback = _extract_informative_error(str(pending.get("output_tail", "") or ""))
            if not feedback:
                feedback = _extract_informative_error(str(pending.get("output_display", "") or ""))

        if not feedback:
            feedback = _extract_informative_error(str(state.get("analysis_result", "") or ""))
        if not feedback:
            feedback = "代码执行失败"

        state["feedback"] = feedback
        _preview = feedback if len(feedback) <= 500 else feedback[:500] + "..."
        print(f"  --> 设置反馈信息: {_preview}")
        attempt_timing = state.get("current_attempt_timing")
        if isinstance(attempt_timing, dict):
            attempt_timing["prepare_retry"] = {
                "feedback_extract_elapsed_seconds": round(time.perf_counter() - retry_start, 3),
                "feedback_length": len(feedback),
                "total_elapsed_seconds": round(time.perf_counter() - retry_start, 3),
            }
            print(f"【Code Dev重试准备耗时】{_format_phase_timing(attempt_timing)}")
    return state


# 构建子图
workflow = StateGraph(CodeAgentState)

# 添加节点
workflow.add_node("generate_code", generate_code)
workflow.add_node("self_reflection", self_reflection)
workflow.add_node("execute_code", execute_code)
workflow.add_node("prepare_retry", prepare_retry)
workflow.add_node("display_result", display_result)

# 定义边
workflow.add_edge(START, "generate_code")
workflow.add_edge("generate_code", "self_reflection")
workflow.add_edge("self_reflection", "execute_code")
# 条件边：根据执行结果决定是否重试
workflow.add_conditional_edges(
    "execute_code",
    should_retry,
    {
        "retry": "prepare_retry",  # 准备重试
        "end": "display_result", # 结束
    }
)

# 准备重试后，返回生成代码节点
workflow.add_edge("prepare_retry", "generate_code")
workflow.add_edge("display_result", END)

# 编译子图
code_agent_graph = workflow.compile()
