"""
Code Developer Agent 子图
负责代码生成和执行
"""
import os
import re
import tempfile
from langchain_core.messages import SystemMessage, HumanMessage, AIMessage
from langgraph.graph import StateGraph, START, END
from .state import CodeAgentState
from src.core.llm import get_llm
from .executor import CodeExecutor
from ._utils.docker_path import convert_to_docker_path
from ._utils.base64_support import create_html_with_base64_image
# 初始化 LLM
llm = get_llm(model_name="qwen-plus", temperature=0.1)


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
    3. 默认值（result_path 默认为 "./result"）
    """
    # 如果 state 中已经有路径，优先使用（但允许从查询中补充缺失的路径）
    has_data_path = bool(state.get("data_path"))
    has_result_path = bool(state.get("result_path"))

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
            state["result_path"] = "./result"

    # 确保 result_path 有默认值
    if not state.get("result_path"):
        state["result_path"] = "./result"

    return state


def generate_code(state: CodeAgentState) -> CodeAgentState:
    """
    生成代码节点
    调用 LLM 生成代码，写入 pending_contribution
    """
    print(f"--- [Code Dev] 正在生成代码 (迭代 {state.get('internal_iteration_count', 0) + 1}) ---")

    # 首先从 user_query 中提取路径（如果 state 中没有）
    state = extract_paths_from_state(state)

    # 构建 Prompt：如果有反馈，说明是修正模式
    context_instruction = ""
    if state.get("feedback"):
        context_instruction = f"""
        【重要！这是修改重试】
        上一次生成的代码或结果被审核员驳回。
        驳回意见/错误信息：{state['feedback']}
        请根据上述意见修改代码。
        """

    # 获取数据路径和结果路径
    data_path = state.get('data_path', '')
    result_path = state.get('result_path', './result')

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

    system_prompt = f"""
你是专业的单细胞数据分析工程师，仅返回【纯Python代码】（无解释、无注释、无markdown）和 【纯requirement.txt包列表】（无解释、无注释），必须严格遵守：
1. 只使用Leiden聚类（sc.tl.leiden），禁止使用Louvain聚类（sc.tl.louvain）；
2. 完整导入所有依赖，确保代码能独立运行；
3. 必须输出analysis_summary变量，格式为：analysis_summary = f"细胞总数：{{adata.n_obs}}，基因总数：{{adata.n_vars}}，聚类数量：{{len(adata.obs['leiden'].cat.categories)}}"
4. UMAP图标题固定为'Clustering UMAP'，无特殊字符；
5. 读取数据时cache=False，anndata设置allow_write_nullable_strings=True。
6. 生成给docker环境的requirements.txt，确保包含所有代码中用到的包。
7. 代码将在 Docker 容器中运行，数据路径为 {docker_data_path}，输出路径为 {docker_output_path}
8. 必须使用 print(f"===RESULT==={{analysis_summary}}===") 输出结果标记
9. 如果生成图片，请保存到 {docker_output_path} 目录

格式：
python代码全部被包括在```python 和```之间
requirement.txt内容全部被包括在```md 和 ```之间
注意：请严格按照上述格式返回内容，确保代码和requirements.txt清晰分隔。
{context_instruction}
    """

    user_prompt = f"""
    任务：{state.get('task', state.get('user_query', ''))}
    数据路径（Docker容器内）：{docker_data_path}
    结果路径（Docker容器内）：{docker_output_path}

    请生成完整的单细胞数据分析代码，包括：
    1. 数据读取（使用 {docker_data_path}）
    2. 数据预处理和质控
    3. 标准化和降维
    4. Leiden 聚类（使用 sc.tl.leiden，不要使用 Louvain）
    5. UMAP 可视化
    6. 结果输出（保存图片到 {docker_output_path}，输出 analysis_summary）
    """

    messages = [
        SystemMessage(content=system_prompt),
        HumanMessage(content=user_prompt)
    ]

    try:
        response = llm.invoke(messages)
        text = response.content

        # 提取代码块和 requirements
        # 支持两种格式：```python 和 ```md（参考 umap_langgraph.py）
        python_pattern = r'```python\n(.*?)\n```'
        requirements_pattern_md = r'```md\n(.*?)\n```'
        requirements_pattern_requirements = r'```requirements\n(.*?)\n```'

        python_match = re.search(python_pattern, text, re.DOTALL)
        # 优先尝试 ```md 格式，如果没有则尝试 ```requirements 格式
        requirements_match = re.search(requirements_pattern_md, text, re.DOTALL)
        if not requirements_match:
            requirements_match = re.search(requirements_pattern_requirements, text, re.DOTALL)

        if python_match:
            code = python_match.group(1).strip()
            print("获取到了 code，前面部分内容：")
            print(code[:300])
        else:
            # 如果没有代码块，尝试提取整个响应
            print("没有获取到 code，尝试提取整个响应")
            code = text.strip()

        if requirements_match:
            requirements = requirements_match.group(1).strip()
            print("获取到了 requirements.txt，内容：")
            print(requirements)
        else:
            print("没有获取到 requirements.txt，使用默认值")
            # 默认 requirements
            requirements = "scanpy>=1.9.0\nmatplotlib>=3.4.0\nnumpy>=1.21.0\npandas>=1.3.0\nscipy>=1.7.0\nanndata>=0.8.0\nigraph\nleidenalg"

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
    os.makedirs(result_path, exist_ok=True)

    # 构建完整的可执行代码（参考 umap_langgraph.py 的改进）
    header = f"""
# 基础库导入（确保代码独立运行）
import sys
import os
sys.path.append(os.getcwd())
import scanpy as sc
import matplotlib.pyplot as plt

# 关键配置
plt.switch_backend('Agg')  # 关闭matplotlib弹窗
sc.settings.verbosity = 3  # 显示Scanpy详细日志
"""
# 核心分析代码（来自大模型生成）
    llm_code = state.get("scanpy_code", "")

    footer = f"""
# 强制输出结果标记（关键：解决list index out of range）
try:
    # 确保即使中间步骤有警告，也能输出结果标记
    if 'analysis_summary' in locals():
        print(f"===RESULT==={{analysis_summary}}===")
    else:
        # 如果没有 analysis_summary，尝试从 adata 提取基本信息
        if 'adata' in locals():
            fallback_result = f"细胞总数：{{adata.n_obs}}，基因总数：{{adata.n_vars}}"
            if 'leiden' in adata.obs.columns:
                fallback_result += f"，聚类数量：{{len(adata.obs['leiden'].cat.categories)}}"
            print(f"===RESULT==={{fallback_result}}===")
        else:
            print(f"===RESULT===代码执行完成，但未找到结果变量===")
except Exception as e:
    # 即使analysis_summary未定义，也输出基础细胞/基因数
    try:
        if 'adata' in locals():
            fallback_result = f"细胞总数：{{adata.n_obs}}，基因总数：{{adata.n_vars}}（聚类步骤失败：{{str(e)[:50]}}）"
            print(f"===RESULT==={{fallback_result}}===")
        else:
            print(f"===RESULT===代码执行失败，无法提取结果：{{str(e)}}===")
    except:
        # 终极兜底：输出空标记，避免index out of range
        print(f"===RESULT===代码执行失败，无法提取结果===")
    """

    # 创建临时目录运行代码
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_script_path = os.path.join(temp_dir, "code.py")
        temp_requirements_path = os.path.join(temp_dir, "requirements.txt")
        full_code = header + "\n# --- LLM Generated Code ---\n" + llm_code + "\n" + footer

        # 将代码和 requirements 写入临时文件
        with open(temp_script_path, "w", encoding="utf-8") as f:
            f.write(full_code)

        with open(temp_requirements_path, "w", encoding="utf-8") as f:
            f.write(requirements)

        # 获取数据路径
        data_path = state.get("data_path", "")

        # 在 Docker 容器中执行代码
        executor = CodeExecutor(
            docker_path=temp_dir,
            data_dir=data_path if data_path and os.path.exists(data_path) else None,
            output_dir=result_path
        )

        try:
            result = executor.execute(timeout=600)  # 10分钟超时

            # 打印执行日志
            print(f"【Docker代码执行日志】: {result.get('output', '')[:500]}...")

            # 提取结果（参考 umap_langgraph.py 的改进）
            output_str = result.get('output', '')
            if "===RESULT===" in output_str:
                # 提取结果部分
                result_part = output_str.split("===RESULT===")[1]
                if "===" in result_part:
                    result_part = result_part.split("===")[0]
                state["analysis_result"] = result_part.strip()
                state["success"] = True  # 标记为成功
                print("  --> 代码执行成功！已提取分析结果")

                # 更新 pending_contribution
                state["pending_contribution"] = {
                    "code": code,
                    "requirements": requirements,
                    "task": state.get("task", ""),
                    "result": state["analysis_result"],
                    "success": True,
                    "output_files": result.get('files', [])
                }
            else:
                # 如果没有找到结果标记，返回错误信息
                error_msg = result.get('error', '未知错误')
                state["analysis_result"] = f"代码执行完成，但未找到结果标记\\n错误日志摘要：{output_str[:500]}"
                state["success"] = False
                print(f"  --> 代码执行失败: {error_msg}")

                # 更新 pending_contribution
                state["pending_contribution"] = {
                    "code": code,
                    "requirements": requirements,
                    "task": state.get("task", ""),
                    "error": error_msg,
                    "success": False
                }

        except Exception as e:
            # 处理其他运行时错误（参考 umap_langgraph.py 的改进）
            error_msg = f"Docker代码运行失败：{str(e)}"
            state["analysis_result"] = f"{error_msg}\\n错误日志：{result if 'result' in locals() else '无'}"
            state["success"] = False
            print(f"  --> {error_msg}")

            # 更新 pending_contribution
            state["pending_contribution"] = {
                "code": code,
                "requirements": requirements,
                "task": state.get("task", ""),
                "error": error_msg,
                "success": False
            }

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
            print(f"Warning: No write access to {result_dir}, using temporary directory: {temp_dir}")
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
        # 显示失败原因（增加调试信息）
        print("运行失败详情：")
        print("-"*30)
        print(state["analysis_result"])
        # 打印提取的结果标记，帮助排查
        print(f"调试信息：")
        print(f"提取的analysis_result原始值：{state.get('analysis_result', '空')[:100]}")
        # print(f"提取的umap_base64原始值：{state.get('umap_base64', '空')[:50]}")
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
        state["feedback"] = state.get("analysis_result", "代码执行失败")
        print(f"  --> 设置反馈信息: {state['feedback'][:100]}...")
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

