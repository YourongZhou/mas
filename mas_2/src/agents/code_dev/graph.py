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
from langchain_core.messages import SystemMessage, HumanMessage, AIMessage
from langgraph.graph import StateGraph, START, END
from .state import CodeAgentState
from src.core.llm import get_llm
from .executor import CodeExecutor
from src.utils.docker_log_summary import summarize_docker_stdout
from src.utils.project_paths import get_mas2_project_root
from src.utils.workflow_skills import (
    format_skill_injection_for_code_dev,
    resolve_workflow_root,
    should_mount_workflow_in_docker,
    use_scanpy_code_style,
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

    docker_output_path = convert_to_docker_path(result_path, 'output') if result_path else '/app/output'

    skill_id = state.get("current_step_skill_id")
    inj_limit = 12000 if use_scanpy_code_style(skill_id) else 4500
    skill_block = format_skill_injection_for_code_dev(skill_id, max_chars=inj_limit)

    if use_scanpy_code_style(skill_id):
        system_prompt = f"""
你是专业的单细胞数据分析工程师，请仅返回 Python 代码和 requirement.txt 包列表（无额外解释），并严格按下方代码块格式输出。

【Docker 路径（不可改）】读取数据：{docker_data_path}；写入产出：{docker_output_path}。
【输出契约】必须使用 print(f"===RESULT==={{analysis_summary}}===")；requirements.txt 须列出代码中所有第三方依赖。
【流程与 API】Leiden/UMAP/邻域图顺序、防御性检查、容错与绘图要求以注入的 SKILL 中 **「MAS Code Agent contract (Docker / generated code)」** 为准，并与本技能 `scripts/` 用法、全文说明一致。

{skill_block}

格式：
python代码全部被包括在```python 和```之间
requirement.txt内容全部被包括在```txt 和 ```之间
注意：请严格按照上述格式返回内容，确保代码和requirements.txt清晰分隔。
【禁止】在 ```python 代码块内部再写任何 ``` 或 ```python 行；代码块内只能是可执行的 Python 源码。
    """
    else:
        system_prompt = f"""
你是 Python 工程师（任务可能是科学计算、组学分析或其它领域）。请仅返回 Python 代码和 requirements.txt 包列表（无额外解释），严格按下方代码块格式输出。

【勿默认套用单细胞 Scanpy 流程】除非任务描述或下方 Workflow/SKILL 明确要求（如 h5ad、Leiden、UMAP），否则不要使用 scanpy，也不要假设存在 AnnData。

【Workflow 技能上下文】
{skill_block}

【硬性约束】
1. 代码在 Docker 内运行。若任务需要读数据/写文件：读取使用 {docker_data_path}，写入使用 {docker_output_path}。若任务为纯计算且无文件 IO，不必强行读盘。
2. 若存在 /app/workflow，可加入 sys.path 并复用其中 `scripts/`，与 SKILL 一致；不要随意重写 SKILL 已有脚本逻辑。
3. 必须：print(f"===RESULT==={{analysis_summary}}===")，analysis_summary 为字符串，概括本步结果。
4. requirements.txt 仅列出代码中实际 import 的第三方包；无第三方依赖时可留空。
5. 【仅当任务需要绘图时】再使用 matplotlib（Agg 后端）、plt.savefig 到 {docker_output_path}，禁用 show=True。纯数值或无图任务不要 import 绘图库。
6. 【禁止】在 ```python 代码块内部再写 ``` 或 ```python；块内仅含可执行 Python。
7. 若用户提供了输入文件列表，必须使用其中的确切路径或 basename（挂载在 /app/data/），禁止臆造其它文件名。
8. 读取 .panel / 样本列表等文本表时：先用 pandas.read_csv(..., sep=r"\\s+", comment='#', header=None) 或打印列名再选列，禁止假设存在名为 ID 的列。

格式：
python 代码在 ```python 与 ``` 之间；requirements 在 ```txt 与 ``` 之间。
    """
        
    # 3. 输出 analysis_summary 变量时需进行安全检查：
    # - 必须包含：细胞总数 (adata.n_obs)、基因总数 (adata.n_vars)
    # - 【仅当存在聚类结果时】才包含：聚类数量 (len(adata.obs['leiden'].cat.categories))
    # - 示例代码：
    #     n_clusters = len(adata.obs['leiden'].cat.categories) if 'leiden' in adata.obs else 0
    #     analysis_summary = f"细胞总数：{{adata.n_obs}}，基因总数：{{adata.n_vars}}，聚类数量：{{n_clusters}}"

    # 构建任务描述，优先使用当前步骤的输入
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
    
    user_prompt = f"""
    任务：{task_description}
    数据路径（Docker容器内）：{docker_data_path}
    结果路径（Docker容器内）：{docker_output_path}
    {file_paths_note}
    {expected_output_note}
    {context_instruction}
    """

    messages = [
        SystemMessage(content=system_prompt),
        HumanMessage(content=user_prompt)
    ]

    # 每次都打印发送给大模型的完整 Prompt (User Prompt, 包含可能带有的修正指令)，便于排错
    print("\n" + "*"*60)
    print(f"🔍 [Code Dev Debug] 发送给 LLM 的提示词 (迭代 {state.get('internal_iteration_count', 0) + 1})")
    print("*"*60)
    print(user_prompt.strip())
    print("*"*60 + "\n")

    try:
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
    _use_scanpy_style = use_scanpy_code_style(_skill_for_exec)

    if _use_scanpy_style:
        header = f"""
# 基础库导入（确保代码独立运行）
import sys
import os
sys.path.append(os.getcwd())
if '/app/workflow' not in sys.path:
    sys.path.insert(0, '/app/workflow')
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
    if os.path.exists('/app/workflow'):
        print(f"DEBUG: Files in /app/workflow: {{os.listdir('/app/workflow')}}")
except Exception as e:
    print(f"DEBUG: Error checking directories: {{e}}")
# --- DEBUG END ---

# 关键配置
plt.switch_backend('Agg')  # 关闭matplotlib弹窗
sc.settings.verbosity = 3  # 显示Scanpy详细日志
"""
    else:
        header = f"""
import sys
import os
if '/app/workflow' not in sys.path:
    sys.path.insert(0, '/app/workflow')
sys.path.append(os.getcwd())

# --- DEBUG START ---
print("DEBUG: Checking /app/data contents...")
try:
    if os.path.exists('/app/data'):
        print(f"DEBUG: Files in /app/data: {{os.listdir('/app/data')}}")
    else:
        print("DEBUG: /app/data does not exist!")
    if os.path.exists('/app/output'):
        print(f"DEBUG: Files in /app/output: {{os.listdir('/app/output')}}")
    if os.path.exists('/app/workflow'):
        print(f"DEBUG: Files in /app/workflow: {{os.listdir('/app/workflow')}}")
except Exception as e:
    print(f"DEBUG: Error checking directories: {{e}}")
# --- DEBUG END ---
"""
# 核心分析代码（来自大模型生成）
    llm_code = state.get("scanpy_code", "")

    footer = f"""
import os

# --- 智能结果输出 ---
try:
    # 1. 优先检查代码中是否定义了自定义摘要
    if 'analysis_summary' in locals():
        print(f"===RESULT==={{analysis_summary}}===")
    
    # 2. 如果没有摘要，但有 adata 对象，输出基础维度信息
    elif 'adata' in locals():
        # 基础信息
        res_str = f"Execution successful. Data shape: {{adata.n_obs}} cells x {{adata.n_vars}} genes."
        
        # 尝试通过 adata.uns/obs 推断做了什么，作为补充信息
        infos = []
        if 'pca' in adata.uns: infos.append("PCA done")
        if 'neighbors' in adata.uns: infos.append("Neighbors computed")
        if 'leiden' in adata.obs: 
            n_clust = len(adata.obs['leiden'].unique())
            infos.append(f"Leiden clusters: {{n_clust}}")
        if 'louvain' in adata.obs:
            n_clust = len(adata.obs['louvain'].unique())
            infos.append(f"Louvain clusters: {{n_clust}}")
            
        if infos:
            res_str += f" (Progress: {{', '.join(infos)}})"
            
        print(f"===RESULT==={{res_str}}===")
        
    # 3. 如果只是普通脚本（没有adata），检查是否有文件生成
    else:
        # 检查 output 目录下的新文件
        out_files = os.listdir('/app/output')
        if out_files:
            print(f"===RESULT===Step completed. Generated files: {{', '.join(out_files)}}===")
        else:
            print(f"===RESULT===Step completed successfully (No specific return value).===")

except Exception as e:
    # 最后的安全网
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
    wf_sid = state.get("current_step_skill_id")
    workflow_host = None
    if should_mount_workflow_in_docker(wf_sid):
        workflow_host = resolve_workflow_root(wf_sid)
        if workflow_host and not state.get("docker_container_id"):
            print(f"  --> Workflow 目录将挂载到容器 /app/workflow: {workflow_host}")

    executor = CodeExecutor(
        docker_path=None,
        container_id=state.get('docker_container_id'),
        data_dir=actual_data_path if actual_data_path and os.path.exists(actual_data_path) else None,
        input_files=input_files if input_files else None,
        output_dir=result_path,
        workflow_host_path=workflow_host,
    )

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
            state['docker_container_id'] = result['container_id']

        # 控制台仅输出关键摘要，避免安装依赖等噪声日志淹没真实错误。
        output_str = result.get('output', '')

        # 检查执行是否成功（从executor返回的success字段）
        executor_success = result.get('success', True)
        
        # 检查output中是否包含错误信息（即使executor返回success=True，代码执行也可能失败）
        has_error_in_output = any(keyword in output_str for keyword in [
            'Traceback', 'Error:', 'Exception:', 'TypeError', 'ValueError', 
            'AttributeError', 'NameError', 'KeyError', 'IndexError'
        ])

        result_part = ""
        result_looks_failed = False
        if "===RESULT===" in output_str:
            result_part = output_str.split("===RESULT===", 1)[1]
            if "===" in result_part:
                result_part = result_part.split("===", 1)[0]
            _res_low = result_part.strip().lower()
            result_looks_failed = any(
                hint in _res_low
                for hint in (
                    "analysis failed",
                    "failed",
                    "error",
                    "exception",
                    "traceback",
                    "unable to",
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

    if state.get("success", False):
        print("\n=== [Code Dev] 展现分析结果 ===")
        # 显示文本结果（优先保证文本能看到）
        print("单细胞分析结果：")
        print("-"*30)
        print(state["analysis_result"])

        # 显示所有PNG图片（增加容错）
        print("正在处理PNG图片并创建对应的HTML文件：")
        print("-"*30)

        import os
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

    # 清理复用的 Docker 容器
    if state.get("docker_container_id"):
        try:
            import docker
            client = docker.from_env()
            container = client.containers.get(state["docker_container_id"])
            container.remove(force=True)
            print(f"  --> [Code Dev] 已清理复用的 Docker 容器: {state['docker_container_id'][:12]}")
        except Exception as e:
            print(f"  --> [Code Dev] 重置/清理复用 Docker 容器失败: {e}")
        # 清除记录，避免下次任务复用
        state["docker_container_id"] = None

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

