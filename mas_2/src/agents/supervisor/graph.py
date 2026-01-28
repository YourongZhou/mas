"""
Supervisor Agent 子图
负责调度决策，决定下一个执行的 worker
"""
from langchain_core.messages import SystemMessage, HumanMessage
from langgraph.graph import StateGraph, START, END
from pydantic import BaseModel, Field
from typing import Literal
from .state import SupervisorAgentState
from src.core.llm import get_llm

# 初始化 LLM
llm = get_llm(model_name="qwen-plus", temperature=0.5)


class RouteDecision(BaseModel):
    """路由决策模型"""
    next_worker: Literal["rag_researcher", "code_dev", "tool_caller", "critic", "FINISH"] = Field(
        ...,
        description="下一个要执行的 worker"
    )
    reasoning: str = Field(..., description="决策理由")


def make_decision(state: SupervisorAgentState) -> SupervisorAgentState:
    """
    决策节点
    根据当前状态决定下一个执行的 worker
    """
    print("--- [Supervisor] 正在做调度决策 ---")
    
    # 获取当前状态
    user_query = state.get("user_query", "")
    rag_context = state.get("rag_context", "")
    code_solution = state.get("code_solution", "")
    final_report = state.get("final_report", "")
    is_approved = state.get("is_approved", False)
    last_worker = state.get("last_worker", "")
    pending_contribution = state.get("pending_contribution")
    
    # 构建决策 Prompt
    prompt = f"""
    你是项目经理，负责协调多个 AI 代理完成用户任务。
    
    当前项目状态：
    - 用户问题: {user_query}
    - RAG 上下文: {"已获取" if rag_context else "未获取"}
    - 代码解决方案: {"已生成" if code_solution else "未生成"}
    - 最终报告: {"已生成" if final_report else "未生成"}
    - 上一个 Worker: {last_worker}
    - 待审核内容: {"有" if pending_contribution else "无"}
    - 审核状态: {"已通过" if is_approved else "未通过/待审核"}
    
    可用的 Worker：
    1. rag_researcher: 检索相关文献和文档
    2. code_dev: 生成和执行代码
    3. tool_caller: 分析用户问题，调用工具返回结果（比如给定DE 基因鉴定细胞型种类）
    4. critic: 审核工作成果
    5. FINISH: 任务完成
    
    请根据当前状态，决定下一个要执行的 worker。
    决策原则：
    - 如果还没有 RAG 上下文，优先调用 rag_researcher；如果没有需要检索的文献和文档，优先调用 code_dev；
    - 如果有待审核内容，必须调用 critic
    - 如果需要生成代码，调用 code_dev
    - 如果需要调用现有工具，调用 tool_caller
    - 如果所有工作都已完成，返回 FINISH
    
    请给出你的决策和理由。
    """
    
    try:
        chain = llm.with_structured_output(RouteDecision)
        decision = chain.invoke(prompt)
        
        print(f"  --> 决策: {decision.next_worker}")
        print(f"  --> 理由: {decision.reasoning}")
        
        # 更新状态
        if decision.next_worker == "FINISH":
            state["next_worker"] = "FINISH"
        else:
            state["next_worker"] = decision.next_worker
        
    except Exception as e:
        print(f"  --> 决策失败: {e}，默认选择 critic")
        # 如果有待审核内容，默认选择 critic
        if pending_contribution:
            state["next_worker"] = "critic"
        else:
            state["next_worker"] = "rag_researcher"
    
    return state


# 构建子图
workflow = StateGraph(SupervisorAgentState)

# 添加节点
workflow.add_node("make_decision", make_decision)

# 定义边
workflow.add_edge(START, "make_decision")
workflow.add_edge("make_decision", END)

# 编译子图
supervisor_agent_graph = workflow.compile()

