import os
import re
import base64
from typing import TypedDict, Optional, List
from langgraph.graph import StateGraph, END
from dashscope import Generation
import dashscope

# 引入你的 GlobalState 定义，用于类型提示
from src.schema import GlobalState 
from src.executor import CodeExecutor # 引入你写好的 Executor

# ================= 1. 定义子图私有状态 =================
class ScanpyAgentState(TypedDict):
    """子图内部使用的状态，父图不可见"""
    task: str               # 当前任务
    feedback: str           # 来自 Critic 的反馈（如果是驳回重做）
    data_path: str          # 数据路径
    result_path: str        # 结果路径
    scanpy_code: str        # 生成的代码
    requirements_txt: str   # 依赖
    analysis_result: str    # 执行结果文本
    success: bool           # 是否执行成功

# ================= 2. 改造生成函数 (支持 Feedback) =================
def generate_scanpy_code(state: ScanpyAgentState) -> ScanpyAgentState:
    dashscope.api_key = os.getenv("DASHSCOPE_API_KEY") # 建议从环境变量取

    # 构建 Prompt：如果有反馈，说明是修正模式
    context_instruction = ""
    if state.get("feedback"):
        context_instruction = f"""
        【重要！这是修改重试】
        上一次生成的代码或结果被审核员驳回。
        驳回意见/错误信息：{state['feedback']}
        请根据上述意见修改代码。
        """
    
    system_prompt = f"""
    你是专业的单细胞数据分析工程师。
    {context_instruction}
    ... (这里放入你原本的 System Prompt) ...
    """
    
    # 简化：这里直接用 state['task'] 作为 user input
    # 实际项目中可以将 path 动态注入 prompt
    user_prompt = f"任务：{state['task']}\n数据路径：{state['data_path']}\n..."

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": user_prompt}
    ]

    try:
        # ... (保持你原有的调用逻辑) ...
        response = Generation.call(model="qwen-plus", messages=messages, temperature=0.1)
        text = response.output["text"]
        
        # ... (保持你原本的正则提取逻辑) ...
        # 假设提取到了 code 和 req
        python_pattern = r'```python\n(.*?)\n```'
        match = re.search(python_pattern, text, re.DOTALL)
        if match:
            state["scanpy_code"] = match.group(1).strip()
            state["requirements_txt"] = "scanpy\nmatplotlib\n..." # 简化示例
        else:
            # 兜底：如果没提取到，可能模型直接返回了代码
            state["scanpy_code"] = text 

    except Exception as e:
        state["scanpy_code"] = f"# Error generating code: {e}"
    
    return state

# ================= 3. 运行函数 (保持不变，直接复用) =================
def run_scanpy_code(state: ScanpyAgentState) -> ScanpyAgentState:
    # ... (这里完全复制你原本的 run_scanpy_code 逻辑) ...
    # 只要确保 input/output 都是 ScanpyAgentState 即可
    
    # 为了演示，假设 executor 逻辑在这里调用
    # executor = CodeExecutor(...)
    # state["analysis_result"] = "模拟运行结果：细胞聚类完成..."
    # state["success"] = True
    return state

# ================= 4. 结果处理 (可选) =================
# 你原来的 display_result 主要是打印，现在我们需要它作为最后一步
# 确保结果写入 state["analysis_result"] 即可

# ================= 5. 【关键】适配器逻辑 =================

def input_adapter(global_state: GlobalState) -> ScanpyAgentState:
    """
    【入参映射】：GlobalState -> ScanpyAgentState
    这是父图进入子图的第一步。
    """
    # 1. 获取任务：取最后一条 User 消息，或者是 Supervisor 指派的具体内容
    last_msg = global_state["messages"][-1].content
    
    # 2. 获取反馈：如果 Critic 驳回，应该在 GlobalState 里有记录
    # 假设 GlobalState 有一个 'critic_feedback' 字段，或者从 messages 里找上一条 Critic 的消息
    feedback = global_state.get("critic_feedback", "")
    
    # 3. 动态确定路径（可以是写死的，也可以是配置的）
    data_path = "/home/zwhuang/Data/pbmc3k/..." 
    
    return {
        "task": last_msg,
        "feedback": feedback,
        "data_path": data_path,
        "result_path": "./result",
        "scanpy_code": "",
        "requirements_txt": "",
        "analysis_result": "",
        "success": False
    }

def output_adapter(subgraph_state: ScanpyAgentState) -> dict:
    """
    【出参映射】：ScanpyAgentState -> GlobalState 更新量
    子图跑完后，要把结果转换回父图能懂的格式（Message）。
    """
    result_text = subgraph_state["analysis_result"]
    if not subgraph_state["success"]:
        result_text = f"代码执行失败：{result_text}"
    
    from langchain_core.messages import AIMessage
    
    # 构造一条 AI Message 返回给父图
    # 这里我们清空 critic_feedback，因为任务已经重新执行了
    return {
        "messages": [AIMessage(content=f"【Code Agent 报告】\n{result_text}")],
        "code_result": result_text,
        "critic_feedback": "" # 清空反馈，等待下一次 Critic
    }

# ================= 6. 构建并编译子图 =================

workflow = StateGraph(ScanpyAgentState)

workflow.add_node("generate", generate_scanpy_code)
workflow.add_node("run", run_scanpy_code)

workflow.set_entry_point("generate")
workflow.add_edge("generate", "run")
workflow.add_edge("run", END)

# 编译子图 app
scanpy_app = workflow.compile()

# ================= 7. 导出给父图使用的节点函数 =================

def code_worker_node(state: GlobalState):
    """
    这是挂载到父图上的实际节点函数。
    它负责：Global -> Local (Invoke) -> Global
    """
    # 1. 转换输入
    sub_input = input_adapter(state)
    
    # 2. 调用子图
    sub_output = scanpy_app.invoke(sub_input)
    
    # 3. 转换输出
    return output_adapter(sub_output)

# 导出这个 scanpy_app 或者封装好的 code_worker_node 都可以
# 推荐导出 code_worker_node，这样父图写法最简单
code_app = code_worker_node