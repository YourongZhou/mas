"""
Code Developer Agent 子图
负责代码生成和执行
"""
import os
import re
import tempfile
import uuid
from datetime import datetime
from pathlib import Path
from langchain_core.messages import SystemMessage, HumanMessage, ToolMessage
from langgraph.graph import StateGraph, START, END
from .state import CodeAgentState
from src.core.llm import get_llm
from src.tools import read_local_file, list_directory
from .executor import CodeExecutor
from src.utils.docker_log_summary import summarize_docker_stdout
from src.utils.project_paths import get_mas2_project_root
from src.utils.workflow_skills import (
    format_skill_injection_for_code_dev,
    resolve_workflow_root,
    should_mount_workflow_in_docker,
)
from ._utils.docker_path import convert_to_docker_path
from ._utils.llm_code_sanitize import sanitize_llm_python_block
from ._utils.base64_support import create_html_with_base64_image
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

    # 兜底识别：如“data/ 目录下”这类描述（优先于文件名识别）
    if not paths["data_path"]:
        dir_match = re.search(r'([A-Za-z0-9_./\\-]+[/\\])\s*目录', user_query, re.IGNORECASE)
        if dir_match:
            dir_path = dir_match.group(1).strip('"\' ')
            paths["data_path"] = dir_path

    # 兜底识别：从自然语言中提取常见数据文件路径（含组学常见后缀）
    if not paths["data_path"]:
        file_match = re.search(
            r"([A-Za-z0-9_./\\-]+\.(?:vcf\.gz|vcf\.bgz|bcf|vcf\.gz\.tbi|vcf|panel|bed|bim|fam|h5ad|h5|csv|tsv|mtx|loom))",
            user_query,
            re.IGNORECASE,
        )
        if file_match:
            file_path = file_match.group(1).strip('"\' ')
            paths["data_path"] = file_path

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

        # 更新 data_path（如果 state 中没有）
        if not has_data_path and parsed_paths["data_path"]:
            # 验证路径是否存在
            if os.path.exists(parsed_paths["data_path"]):
                state["data_path"] = parsed_paths["data_path"]
                print(f"  --> 从用户查询中解析到数据路径: {parsed_paths['data_path']}")
            else:
                print(f"  --> 警告：解析到的数据路径不存在: {parsed_paths['data_path']}")
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
    print(f"--- [Code Dev] 正在生成代码 (迭代 {state.get('internal_iteration_count', 0) + 1}) ---")

    critic_feedback = state.get("critique_feedback", "")
    internal_feedback = state.get("feedback", "")
    
    final_feedback = ""
    if critic_feedback:
        print(f"  --> 收到 Critic 的驳回意见: {critic_feedback[:500]}...")
        final_feedback = critic_feedback
    elif internal_feedback:
        print(f"  --> 收到内部执行的错误反馈: {internal_feedback[:500]}...")
        final_feedback = internal_feedback

    # 首先从 user_query 中提取路径（如果 state 中没有）
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
    # 统一使用较大的注入限制以兼容复杂 SKILL
    inj_limit = 12000
    skill_block = format_skill_injection_for_code_dev(skill_id, max_chars=inj_limit)

    req_instruction = "requirements.txt 须且仅需列出代码中实际调用的第三方依赖环境包，若是标准库可留空。"
    if state.get("global_requirements"):
        req_instruction = "Supervisor 已经提前为你生成并安装了该任务全局所需的 global requirements。因此，在此步骤生成代码时，请【留空 requirements.txt】！**除非**你明确在修复之前的环境报错（如 `ModuleNotFoundError`）并需要额外补充缺失包，才在 requirements.txt 中列出。"

    system_prompt = f"""
你是专业的 Python工程师（专注于生物信息数据分析、科学计算及绘图代码生成）。请仅返回 Python 代码和 requirements.txt 包列表（无额外解释），严格按下方代码块格式输出。

【Workflow 技能上下文】
{skill_block}

【代码生成硬性约束】
1. 路径规范：代码在 Docker 内运行。由 {docker_data_path} 读取数据，将产出写入 {docker_output_path}。针对纯计算任务无文件 IO 则无需强行读盘。
2. 勿默认套用单细胞流程：除非任务/SKILL 明确要求（如 h5ad、Leiden、UMAP），否则不要随意使用 scanpy 库或假设存在 AnnData。
3. 脚本复用与流程：若存在 /app/workflows，请加入 sys.path 并复用其中 `scripts/` 的已有逻辑。API、容错与绘图均以 SKILL 为准，不要随意重写已有脚本核心逻辑。
4. 输出契约：必须在脚本末尾使用 print(f"===RESULT==={{analysis_summary}}===") 以向系统反馈摘要（analysis_summary 应为一段表示执行概括的字符串）。
5. 依赖声明：{req_instruction}
6. 绘图准则：仅当明确需要绘图任务时才引入 matplotlib（必须选用 Agg 后端），以 plt.savefig 保存至 {docker_output_path} 下即可，严禁 show=True。纯数值任务免引入绘图库。
7. 文件校验：若用户提供了输入文件列表，必须使用其确切路径或挂载名，禁止臆造其它文件名。
8. 读表容错：遇到读取 .panel / 样本列表等文本时，要求使用 pandas.read_csv(..., sep=r"\\\\s+", comment='#', header=None) 此类容错方法，绝对禁止假定表格中存在名为 ID 的列。

格式要求：
python 代码全部被包括在 ```python 和 ``` 之间。
requirements.txt 内容全部被包括在 ```txt 和 ``` 之间。
【禁止】在 ```python 代码块内部再写任何 ``` 或 ```python 行；代码块内必须且只能为可执行的 Python 源码。
    """

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
    file_paths_note = ""
    if input_files:
        file_paths_note += f"\n【输入文件】\n" + "\n".join([f"- {f}" for f in input_files])
        file_paths_note += (
            "\n（Docker 内数据目录为 /app/data/，请用 os.path.join('/app/data', '<basename>') 或下列 basename 读取；"
            "勿使用未出现在上述列表中的文件名。）"
        )
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

    user_prompt = f"""{historical_outputs_context}{exploration_context}
【当前生成任务需求】
{task_description}

{expected_output_note}
{file_paths_note}

{context_instruction}
    """

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

    try:
        tools = [list_directory, read_local_file]
        llm_with_tools = llm.bind_tools(tools)
        
        # =============== 新增：工具调用专用的信息收集阶段 ===============
        is_exploration = not state.get("data_exploration_done", True)
        is_retry = state.get("internal_iteration_count", 0) > 0
        tool_context_str = state.get("step_tool_context", "")

        if is_exploration:
            print(f"🔍 [Code Dev] 当前处于数据探查阶段，直接生成探查代码（跳过 SKILL 工具调研）...")
        elif is_retry and tool_context_str:
            print(f"🔄 [Code Dev] 当前步骤为重试迭代，使用已有的 Tool Context 进行重试，不再重新调研...")
        else:
            print("\n" + "="*60)
            print(f"🔍 [Code Dev Tool Loop] 开始为执行任务查阅 SKILL 目录...")
            wf_host_path = resolve_workflow_root(skill_id) if skill_id else ""
            
            tool_system_prompt = f"""
你是代码生成前的问题拆解与调研专家。
即将执行的环节对应技能(SKILL): `{skill_id}`。
你的工作任务：先使用 `list_directory` 浏览该技能对应的文件目录结构，必须重点关注 `{wf_host_path}` 等相关路径；随后请从目录中找出与接下来任务最相关的极少数核心 `scripts` / `references` 脚本或示例，调用 `read_local_file` 查看 1-5 个最相关文件内容。切忌盲目读取所有脚本，完整的 SKILL.md 信息将在此调研结束后自动添加给实际编写代码的模型！

【要求】
- 必须且仅允许使用 Tool Calling （如 list_directory / read_local_file）发起工具调用请求。
- 严禁在此刻生成任何以 ```python 格式代表最终代码的代码块，不要做模拟总结。
- 当你认为已查阅到所需的所有上下文时，请直接回复字符串文本 "TOOL_CALLS_DONE" 以标志文档查阅彻底结束。
"""
            tool_messages = [
                SystemMessage(content=tool_system_prompt),
                HumanMessage(content=f"当前生成任务需求：\n{task_description}\n\n请立刻开始使用工具收集信息，完毕后回复 'TOOL_CALLS_DONE'。")
            ]
            
            tool_iterations = 0
            max_tool_iterations = 10
            collected_info = []

            while tool_iterations < max_tool_iterations:
                response = llm_with_tools.invoke(tool_messages)
                tool_messages.append(response)
                
                if not response.tool_calls:
                    break

                print(f"  [Code Dev Tool Loop] 正在执行第 {tool_iterations+1} 轮推理...")
                for tool_call in response.tool_calls:
                    print(f"  [Code Dev Tool Loop] 调用工具: {tool_call['name']}({tool_call['args']})")
                    try:
                        if tool_call["name"] == "list_directory":
                            tool_result = list_directory.invoke(tool_call['args'])
                        elif tool_call["name"] == "read_local_file":
                            tool_result = read_local_file.invoke(tool_call['args'])
                        else:
                            tool_result = f"Error: Unknown tool {tool_call['name']}"
                    except Exception as e:
                        tool_result = f"Error: 工具执行异常: {str(e)}"
                        print(f"  [Code Dev Tool Loop] 异常: {str(e)}")
                    
                    preview = str(tool_result)[:500].replace('\n', ' ')
                    print(f"  [Code Dev Tool Loop] 结果预览: {preview}...")
                    
                    tool_msg = ToolMessage(content=str(tool_result), tool_call_id=tool_call["id"], name=tool_call["name"])
                    tool_messages.append(tool_msg)
                    
                    # 截取一部分结果供后续 Code 生成使用（避免超长提示词卡死）
                    collected_info.append(f"---- 工具 {tool_call['name']} 目标 {tool_call['args']} 返回：\n{str(tool_result)[:3000]}")
                
                tool_iterations += 1
                
            if collected_info:
                tool_context_str = "\n【来自工具调研获取的最佳背景内容】\n" + "\n".join(collected_info) + "\n\n请结合上方参考资料和 SKILL 要求，完成任务代码编写。\n"
                state["step_tool_context"] = tool_context_str

            print(f"🔍 [Code Dev Tool Loop] 调研结束，进入正式代码生成...")
            print("="*60 + "\n")
        # =============== 信息收集阶段结束 ===============

        # 将工具内容合并进 user_prompt 作为上下文，并添加 SKILL_BLOCK 要求
        final_user_prompt = f"{tool_context_str}\n【Workflow技能参考】：\n{skill_block}\n\n【正式任务：所需功能要求】\n{user_prompt}"

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
        response = llm.invoke(messages)
        text = response.content

        # 提取代码块和 requirements
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
            # 默认 requirements
            requirements = "scanpy>=1.9.0\nmatplotlib>=3.4.0\nnumpy>=1.21.0\npandas>=1.3.0\nscipy>=1.7.0\nanndata>=0.8.0\nigraph"

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

        print(f"  --> 代码生成成功，代码长度: {len(code)} 字符")

    except Exception as e:
        error_msg = f"代码生成失败: {str(e)}"
        print(f"  --> {error_msg}")
        state["scanpy_code"] = f"代码生成失败: {e}"
        state["requirements_txt"] = f"requirements.txt 生成失败：{str(e)}"
        state["pending_contribution"] = {"error": error_msg}
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

    print("--- [Code Dev] 进行自我检查 ---")

    # 简单的安全检查
    dangerous_patterns = ["eval(", "exec(", "__import__", "open("]
    warnings = []
    for pattern in dangerous_patterns:
        if pattern in code:
            warnings.append(f"检测到潜在风险: {pattern}")

    if warnings:
        print(f"  --> 警告: {', '.join(warnings)}")

    return state


def execute_code(state: CodeAgentState) -> CodeAgentState:
    """
    执行代码节点
    在 Docker 容器中执行生成的代码
    """
    print("--- [Code Dev] 正在执行代码 ---")

    code = state.get("scanpy_code", "")
    requirements = state.get("requirements_txt", "")

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
import scanpy as sc
import matplotlib.pyplot as plt
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

# 关键配置
plt.switch_backend('Agg')  # 关闭matplotlib弹窗
sc.settings.verbosity = 3  # 显示Scanpy详细日志
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

    # 智能确定数据路径：检查输入文件的实际位置
    data_path = state.get("data_path", "")
    current_step_file_paths = state.get("current_step_file_paths", {})
    input_files = current_step_file_paths.get("input_files", []) if current_step_file_paths else []
    
    # 如果计划中指定了输入文件，检查它们实际存在的位置
    actual_data_path = data_path
    if input_files:
        # 检查第一个输入文件的实际位置
        first_input = input_files[0]
        
        # 如果输入文件路径是绝对路径且存在，使用其所在目录
        if os.path.isabs(first_input) and os.path.exists(first_input):
            actual_data_path = os.path.dirname(first_input) if os.path.isfile(first_input) else first_input
            print(f"  --> 检测到输入文件: {first_input}")
            print(f"  --> 使用输入文件所在目录作为数据路径: {actual_data_path}")
        # 如果输入文件路径是相对路径，尝试在 result_path 中查找
        elif not os.path.isabs(first_input):
            # 尝试在 result_path 中查找（可能是上一轮的输出）
            candidate_paths = [
                os.path.join(result_path, first_input),  # result_path/input_file
                os.path.join(result_path, os.path.basename(first_input)),  # result_path/filename
                first_input  # 直接使用相对路径
            ]
            
            for candidate in candidate_paths:
                if os.path.exists(candidate):
                    actual_data_path = os.path.dirname(candidate) if os.path.isfile(candidate) else candidate
                    print(f"  --> 在候选路径中找到输入文件: {candidate}")
                    print(f"  --> 使用该输入文件所在目录作为数据路径: {actual_data_path}")
                    break
            else:
                # 仅文件名时：在 mas_2/data 与 ./data 下按 basename 查找（不依赖 Notebook CWD）
                base = os.path.basename(first_input)
                found_in_project_data = False
                for root_candidate in _mas2_data_dir_candidates():
                    if not root_candidate or not os.path.isdir(root_candidate):
                        continue
                    cand = os.path.join(root_candidate, base)
                    if os.path.exists(cand):
                        actual_data_path = root_candidate
                        print(f"  --> 在数据目录 {root_candidate} 中找到输入文件: {cand}")
                        found_in_project_data = True
                        break
                if not found_in_project_data and data_path and os.path.exists(data_path):
                    candidate_in_data = os.path.join(data_path, first_input) if os.path.isdir(data_path) else None
                    if candidate_in_data and os.path.exists(candidate_in_data):
                        actual_data_path = data_path
                        print(f"  --> 在原始 data_path 中找到输入文件: {candidate_in_data}")
                    else:
                        # 如果还是找不到，使用 result_path（因为可能是上一轮的输出）
                        actual_data_path = result_path
                        print(f"  --> 未找到输入文件，使用 result_path 作为数据路径: {actual_data_path}")
        # 如果输入文件路径是绝对路径但不存在，检查是否在 result_path 中
        elif os.path.isabs(first_input) and not os.path.exists(first_input):
            filename = os.path.basename(first_input)
            candidate_in_result = os.path.join(result_path, filename)
            if os.path.exists(candidate_in_result):
                actual_data_path = result_path
                print(f"  --> 输入文件不存在，但在 result_path 中找到同名文件: {candidate_in_result}")
                print(f"  --> 使用 result_path 作为数据路径: {actual_data_path}")
            else:
                for root_candidate in _mas2_data_dir_candidates():
                    if not root_candidate or not os.path.isdir(root_candidate):
                        continue
                    cand = os.path.join(root_candidate, filename)
                    if os.path.exists(cand):
                        actual_data_path = root_candidate
                        print(f"  --> 绝对路径无效，在 {root_candidate} 中找到同名文件: {cand}")
                        break
    
    # 如果 actual_data_path 为空或不存在，优先尝试 mas_2/data 与 ./data（不依赖 Notebook CWD）
    if not actual_data_path or not os.path.exists(actual_data_path):
        for default_data_dir in _mas2_data_dir_candidates():
            if default_data_dir and os.path.isdir(default_data_dir):
                actual_data_path = default_data_dir
                print(f"  --> 数据路径无效，使用默认 data 目录: {actual_data_path}")
                break

    # 如果仍无有效数据路径，再回退到 result_path
    if not actual_data_path or not os.path.exists(actual_data_path):
        if result_path and os.path.exists(result_path):
            actual_data_path = result_path
            print(f"  --> 数据路径无效，使用 result_path 作为数据路径: {actual_data_path}")
    
    # 确保数据路径存在
    if actual_data_path and not os.path.exists(actual_data_path):
        print(f"  --> 警告：数据路径不存在: {actual_data_path}，将尝试创建")
        try:
            os.makedirs(actual_data_path, exist_ok=True)
        except Exception as e:
            print(f"  --> 无法创建数据路径: {e}")

    # 在 Docker 容器中执行代码
    # 传递 input_files 以便 executor 智能确定需要挂载的目录
    from src.utils.workflow_skills import get_workflows_root
    workflow_host = get_workflows_root()
    if workflow_host and not state.get("docker_container_id"):
        print(f"  --> Workflow 目录将挂载到容器 /app/workflows: {workflow_host}")

    executor = CodeExecutor(
        docker_path=None,
        container_id=state.get('docker_container_id'),
        data_dir=actual_data_path if actual_data_path and os.path.exists(actual_data_path) else None,
        input_files=input_files if input_files else None,
        output_dir=result_path,
        workflow_host_path=workflow_host,
    )

    # 首次启动容器才合并 global_requirements 进行集中安装
    if not state.get('docker_container_id') and state.get('global_requirements'):
        print(f"  --> [初次执行] 组合并预安装 global_requirements")
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

    try:
        result = executor.execute(code_str=full_code, requirements_str=requirements, timeout=600)

        if result.get('container_id'):
            new_cid = result['container_id']
            if new_cid != state.get('docker_container_id'):
                print(f"  --> [容器建立成功] 已成功构建持久化容器，容器ID: {new_cid[:12]}")
            state['docker_container_id'] = new_cid

        # 控制台仅输出关键摘要，避免安装依赖等噪声日志淹没真实错误。
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
    max_retries = 10
    iteration_count = state.get("internal_iteration_count", 0)
    success = state.get("success", False)

    if success:
        return "end"
    elif iteration_count < max_retries:
        print(f"  --> 执行失败，准备重试 (第 {iteration_count + 1}/{max_retries} 次)")
        return "retry"
    else:
        print(f"  --> 已达到最大重试次数 ({max_retries})，停止重试")
        return "end"


def prepare_retry(state: CodeAgentState) -> CodeAgentState:
    """
    准备重试节点
    将错误信息作为反馈，用于下一次代码生成
    """
    if not state.get("success", False):
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

