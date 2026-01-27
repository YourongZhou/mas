import dashscope
import os
import tempfile
import base64
import re
from typing import TypedDict  # 用于定义状态类

# ========== LangGraph核心（构建Agent工作流） ==========
from langgraph.graph import StateGraph, END

# ========== DashScope（调用通义千问模型） ==========
import dashscope
from dashscope import Generation

# ========== Scanpy+可视化（单细胞分析核心） ==========
import scanpy as sc
import matplotlib.pyplot as plt

# ========== 环境适配（避免弹窗/设置样式） ==========
plt.switch_backend('Agg')  # 关闭matplotlib弹窗
sc.settings.verbosity = 3  # 显示Scanpy详细日志
sc.set_figure_params(dpi=80, facecolor="white")  # 设置图片样式


def convert_to_docker_path(local_path, mode = "data"):
    """
    将本地文件路径转换为Docker容器内的路径

    Args:
        local_path (str): 本地文件路径，例如 '/home/zwhuang/Data/scHetero/adata_raw_new.h5ad'
        base_docker_path (str): Docker容器内的基础路径，默认为 '/app/data'

    Returns:
        str: Docker容器内的路径
    """
    if mode == "data":
        base_docker_path = "/app/data"
    elif mode == "output":
        base_docker_path = "/app/output"
    else:
        raise ValueError("mode参数必须为'data'或'output'")
    # Check if the path is a directory
    if os.path.isdir(local_path):
        # If it's a directory, the whole directory is mounted to base_docker_path
        return base_docker_path
    else:
        # If it's a file, append the filename to base_docker_path
        filename = os.path.basename(local_path)
        docker_path = os.path.join(base_docker_path, filename)
        return docker_path


def create_html_with_base64_image(image_path, output_html_path):
    """
    Reads a PNG image, converts it to base64, and creates an HTML file with the embedded image.

    Args:
        image_path (str): Path to the input PNG image
        output_html_path (str): Path to the output HTML file
    """
    # Check if the image file exists
    if not os.path.exists(image_path):
        print(f"Error: Image file '{image_path}' does not exist.")
        return

    # Read the image file in binary mode
    with open(image_path, 'rb') as image_file:
        # Encode the image as base64
        base64_encoded = base64.b64encode(image_file.read()).decode('utf-8')

    # Create HTML content with the base64-encoded image
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>UMAP Clustering</title>
</head>
<body>
    <h1>Leiden Clustering UMAP</h1>
    <img src="data:image/png;base64,{base64_encoded}" alt="UMAP Visualization">
</body>
</html>"""

    # Write the HTML content to the output file
    with open(output_html_path, 'w') as html_file:
        html_file.write(html_content)

    print(f"Successfully created HTML file with base64-encoded image: {output_html_path}")


class ScanpyAgentState(TypedDict):
    """
    单细胞分析Agent的状态类（替换教程的代码修复AgentState）
    每个字段存储Agent运行的关键信息，全程可追溯
    """
    task: str               # 你的分析任务（比如"PBMC3k做UMAP聚类"）
    data_path: str          # PBMC3k数据的本地路径
    result_path: str        # 结果保存路径
    scanpy_code: str        # 大模型生成的Scanpy分析代码
    requirements_txt: str   # 大模型生成的requirements.txt内容
    analysis_result: str    # 分析结果（细胞数、聚类数等）
    # umap_base64: str        # UMAP图的base64编码（用于显示）
    success: bool           # 代码运行是否成功（True/False)


# ====================== 节点1：生成Scanpy分析代码 ======================
def generate_scanpy_code(state: ScanpyAgentState) -> ScanpyAgentState:
    """
    功能：调用通义千问模型，根据你的任务+数据路径生成完整的Scanpy代码
    关键：强制模型生成analysis_summary变量（解决没结果的问题）
    """
    # 1. 配置DashScope API Key（替换为你的有效Key）
    dashscope.api_key = "sk-f5f0fa472cb545ea9faa4ed6c2e45c42"

    # 2. System Prompt
    system_prompt = """
你是专业的单细胞数据分析工程师，仅返回【纯Python代码】（无解释、无注释、无markdown）和 【纯requirement.txt包列表】（无解释、无注释），必须严格遵守：
1. 只使用Leiden聚类（sc.tl.leiden），禁止使用Louvain聚类（sc.tl.louvain）；
2. 完整导入所有依赖，确保代码能独立运行；
3. 必须输出analysis_summary;
4. UMAP图标题固定为'Clustering UMAP'，无特殊字符；
5. 读取数据时cache=False，anndata设置allow_write_nullable_strings=True。
6. 生成给docker环境的requirements.txt，确保包含所有必要包。

格式：
python代码全部被包括在```python 和```之间
requirement.txt内容全部被包括在```md 和 ```之间
注意：请严格按照上述格式返回内容，确保代码和requirements.txt清晰分隔。
    """

    # 3. User Prompt
    user_prompt = f"""
请生成单细胞数据分析的Scanpy代码，数据路径为：{convert_to_docker_path(state['data_path'], 'data')}
代码内容：
# 完整导入
import anndata
import scanpy as sc
import matplotlib.pyplot as plt
import igraph
import leidenalg

# 关键配置
anndata.settings.allow_write_nullable_strings = True

# 读取数据
adata = sc.read_10x_mtx('{convert_to_docker_path(state['data_path'], 'data')}', var_names='gene_symbols', cache=False)

# 质控
adata.var_names_make_unique()
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
adata = adata[(adata.obs.n_genes_by_counts < 2500) & (adata.obs.pct_counts_mt < 5)]

# 标准化+降维
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata)
sc.pp.neighbors(adata)

# 核心：Leiden聚类（仅用这个，禁用Louvain）
sc.tl.leiden(adata, resolution=1.0)
sc.tl.umap(adata)

# 结果定义
analysis_summary = f"细胞总数：{{adata.n_obs}}，基因总数：{{adata.n_vars}}，聚类数量：{{len(adata.obs['leiden'].cat.categories)}}"

# 生成UMAP图
sc.pl.umap(adata, color='leiden', title='Leiden Clustering UMAP', show=False)
plt.savefig(f"{convert_to_docker_path(state['result_path'], 'output')}/leiden.png")
plt.close()

# 强制输出结果标记（避免index out of range）
print(f"===RESULT==={{analysis_summary}}===")
# print(f"===IMAGE==={{umap_base64}}===")


2. 请生成对应的requirements.txt内容，确保代码能在docker容器中成功运行，至少包含以下包：
scanpy>=1.9.0
anndata>=0.8.0
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0
igraph
leidenalg
    """

    # 4. 构造通义千问所需的messages格式（System+User）
    messages = [
        {"role": "system", "content": system_prompt.strip()},
        {"role": "user", "content": user_prompt.strip()}
    ]

    # 5. 调用通义千问生成代码
    try:
        response = Generation.call(
            model="qwen-plus",  # 模型名称
            messages=messages,  # 传入拆分后的System+User Prompt
            temperature=0.1     # 低温度，保证代码稳定
        )
        # 提取生成的代码并清理格式
        text = response.output["text"]

        python_pattern = r'```python\n(.*?)\n```'
        md_pattern = r'```md\n(.*?)\n```'
        python_match = re.search(python_pattern, text, re.DOTALL)
        md_match = re.search(md_pattern, text, re.DOTALL)
        if python_match:
            state["scanpy_code"] = python_match.group(1).strip()
            print("生成Scanpy代码：\n" + "-"*60 + "\n", state["scanpy_code"], "\n" + "-"*60)

        if md_match:
            state["requirements_txt"] = md_match.group(1).strip()
            print("生成requirements代码：\n" + "-"*60 + "\n", state["requirements_txt"], "\n" + "-"*60)


    except Exception as e:
        state["scanpy_code"] = f"模型调用失败：{str(e)}"
        state["requirements_txt"] = f"模型调用失败：{str(e)}"
        print(f"模型调用失败：{e}")

    return state

# ====================== 节点2：安全运行Scanpy代码 ======================
# 导入Docker执行器
from executor import CodeExecutor

def run_scanpy_code(state: ScanpyAgentState) -> ScanpyAgentState:
    """
    运行生成的Scanpy代码，提取分析结果和UMAP图base64编码
    核心特性：
        1. 容错处理：代码执行失败也不会崩溃，保留完整错误日志
        2. 强制提取：优先识别===RESULT===/===IMAGE===标记，解决list index out of range
        3. 调试友好：打印stdout/stderr，方便排查代码执行问题
        4. Docker执行：在Docker容器中安全执行代码
    """
    # 1. 拼接完整的可执行代码（补充基础导入，确保代码独立运行）
    full_scanpy_code = f"""
# 基础库导入（确保代码独立运行）
import sys
import os
sys.path.append(os.getcwd())

import scanpy as sc
import matplotlib.pyplot as plt
import anndata
import igraph
import leidenalg

# 核心分析代码（来自大模型生成）
{state['scanpy_code']}

# 强制输出结果标记（关键：解决list index out of range）
try:
    # 确保即使中间步骤有警告，也能输出结果标记
    print(f"===RESULT==={{analysis_summary}}===")
    # 确保UMAP图base64存在，无则赋值为空字符串
    if 'umap_base64' in locals():
        print(f"===IMAGE==={{umap_base64}}===")
    else:
        print(f"===IMAGE======")
except Exception as e:
    # 即使analysis_summary未定义，也输出基础细胞/基因数
    try:
        fallback_result = f"细胞总数：{{adata.n_obs}}，基因总数：{{adata.n_vars}}（聚类步骤失败：{{str(e)[:50]}}）"
        print(f"===RESULT==={{fallback_result}}===")
        print(f"===IMAGE======")
    except:
        # 终极兜底：输出空标记，避免index out of range
        print(f"===RESULT===代码执行失败，无法提取结果===")
        print(f"===IMAGE======")
    """

    # 2. 创建临时目录运行代码（避免污染当前目录）
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_script_path = os.path.join(temp_dir, "code.py")
        temp_txt_path = os.path.join(temp_dir, "requirements.txt")

        # 将拼接好的代码写入临时文件
        with open(temp_script_path, "w", encoding="utf-8") as f:
            f.write(full_scanpy_code)

        with open(temp_txt_path, "w", encoding="utf-8") as f:
            f.write(state['requirements_txt'])

        # # Create requirements.txt file for the Docker container
        # requirements_path = os.path.join(temp_dir, "requirements.txt")
        # with open(requirements_path, "w", encoding="utf-8") as f:
        #     f.write("scanpy>=1.11.0\n")
        #     f.write("matplotlib>=3.5.0\n")
        #     f.write("anndata>=0.11.0\n")
        #     f.write("igraph\n")
        #     f.write("leidenalg\n")
        #     f.write("numpy\n")
        #     f.write("pandas\n")

        # 3. Execute code in Docker container
        executor = CodeExecutor(temp_dir, data_dir=state['data_path'], output_dir = f"{os.path.dirname(os.path.abspath(__file__))}/result")
        try:
            result = executor.execute()

            # 4. Print execution logs for debugging
            print(f"【Docker代码执行日志】: {result}")

            # 5. Extract results from Docker execution output
            output_str = result.get('output', '')
            if "===RESULT===" in output_str:
                # Extract the result part
                result_part = output_str.split("===RESULT===")[1]
                if "===" in result_part:
                    result_part = result_part.split("===")[0]
                state["analysis_result"] = result_part.strip()
                state["success"] = True  # Mark as successful
                print("Scanpy代码在Docker中运行成功，已提取分析结果！")
            else:
                # Fallback: if no marker is found, return basic error info
                state["analysis_result"] = f"代码执行完成，但未找到结果标记\\错误日志摘要：{output_str[:500]}"
                state["success"] = False

            # # Extract UMAP image base64 encoding
            # if "===IMAGE===" in result:
            #     image_part = result.split("===IMAGE===")[1]
            #     if "===" in image_part:
            #         image_part = image_part.split("===")[0]
            #     state["umap_base64"] = image_part.strip()
            # else:
            #     state["umap_base64"] = ""  # Assign empty string if no image, to avoid decode errors

        except Exception as e:
            # Handle other runtime errors (such as file permissions, path issues)
            state["analysis_result"] = f"Docker代码运行失败：{str(e)}\\错误日志：{result if 'result' in locals() else '无'}"
            state["success"] = False
            print(f"Scanpy代码在Docker中运行异常：{e}")

    # Return the state object filled with results
    return state

# ====================== 节点3：显示结果 ======================
def display_result(state: ScanpyAgentState) -> ScanpyAgentState:
    """
    功能：展示分析结果文本+UMAP聚类图（增加容错，避免解码崩溃）
    """

    if state["success"]:
        # 显示文本结果（优先保证文本能看到）
        print("单细胞分析结果：")
        print("-"*30)
        print(state["analysis_result"])

        # 显示UMAP图（增加容错）
        print("UMAP聚类图：")
        print("-"*30)

        # Create output path in a location where we have write permissions
        output_html_path = f"{state['result_path']}/leiden_decoded.html"

        # Check if we have write permissions to the result directory
        import os
        result_dir = state['result_path']
        if not os.access(result_dir, os.W_OK):
            # If no write access to result directory, create a temporary directory for output
            import tempfile
            with tempfile.TemporaryDirectory() as temp_dir:
                output_html_path = os.path.join(temp_dir, "leiden_decoded.html")
                print(f"Warning: No write access to {result_dir}, using temporary directory: {temp_dir}")

                # Call the function to create HTML with base64 image
                create_html_with_base64_image(f"{state['result_path']}/leiden.png", output_html_path)

                # Inform user of location
                print(f"HTML file saved at: {output_html_path}")
                print(f"Please copy the file manually to result directory if needed.")
        else:
            # We have write access, proceed normally
            create_html_with_base64_image(f"{state['result_path']}/leiden.png", output_html_path)
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


# 1. 创建状态图
workflow = StateGraph(ScanpyAgentState)

# 2. 添加节点（对应步骤4的3个函数）
workflow.add_node("generate_code", generate_scanpy_code)  # 生成代码节点
workflow.add_node("run_code", run_scanpy_code)            # 运行代码节点
workflow.add_node("display_result", display_result)       # 显示结果节点

# 3. 设置执行顺序（线性流程，无循环）
workflow.set_entry_point("generate_code")                # 入口：先生成代码
workflow.add_edge("generate_code", "run_code")           # 生成代码后→运行代码
workflow.add_edge("run_code", "display_result")          # 运行代码后→显示结果
workflow.add_edge("display_result", END)                 # 显示结果后→结束

# 4. 编译工作流（生成可执行的Agent）
app = workflow.compile()

# # （可选）可视化工作流（需安装graphviz：pip install graphviz）
# try:
#     display(Image(app.get_graph().draw_mermaid_png()))
#     print("工作流可视化完成！")
# except:
#     print("提示：安装graphviz后可可视化工作流，不影响核心功能。")


# 1. 定义初始状态（仅需修改data_path！）
initial_state = {
    "task": "完成PBMC3k数据集的完整单细胞分析，生成UMAP聚类图并提取关键统计信息",
    "data_path": "/home/zwhuang/Data/pbmc3k/filtered_gene_bc_matrices/hg19",  # 替换为你的PBMC3k数据路径
    "result_path": f"{os.path.dirname(os.path.abspath(__file__))}/result",  # 结果保存路径
    "scanpy_code": "",
    "requirements_txt": "",
    "analysis_result": "",
    # "umap_base64": "",
    "success": False
}

# 确保输出目录存在
output_dir = f"{initial_state['result_path']}"
os.makedirs(output_dir, exist_ok=True)
print(f"确保输出目录存在{output_dir}")

# 2. 启动Agent（核心操作！）
print("开始运行单细胞分析Agent...")
final_state = app.invoke(initial_state)
print("Agent运行结束！")