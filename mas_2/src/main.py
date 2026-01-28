"""
主图编排模块
将所有 SubGraph Agent 连接成完整的工作流
"""
from langgraph.graph import StateGraph, START, END
from src.core.state import GlobalState

# 导入所有编译好的子图
from src.agents.supervisor.graph import supervisor_agent_graph
from src.agents.critic.graph import critic_agent_graph
from src.agents.code_dev.graph import code_agent_graph
from src.agents.rag_researcher.graph import rag_agent_graph
from src.agents.tool_caller.graph import tool_caller_agent_graph


# ==================== Wrapper 节点 ====================

def wrap_rag_researcher(state: GlobalState) -> GlobalState:
    """包装 RAG Researcher，更新 last_worker"""
    result = rag_agent_graph.invoke(state)
    return {
        **result,
        "last_worker": "rag_researcher"
    }


def wrap_code_dev(state: GlobalState) -> GlobalState:
    """包装 Code Dev，更新 last_worker"""
    result = code_agent_graph.invoke(state)
    return {
        **result,
        "last_worker": "code_dev"
    }


def wrap_tool_caller(state: GlobalState) -> GlobalState:
    """包装 Tool Caller，更新 last_worker"""
    result = tool_caller_agent_graph.invoke(state)
    return {
        **result,
        "last_worker": "data_analyst"  # 映射为 data_analyst
    }


# ==================== 路由逻辑 ====================

def supervisor_router(state: GlobalState) -> str:
    """
    Supervisor 路由：根据 next_worker 字段决定去哪个 Agent
    
    Returns:
        目标节点的名称
    """
    next_worker = state.get("next_worker", "FINISH")
    
    if next_worker == "FINISH":
        return "finalize"
    elif next_worker == "rag_researcher":
        return "rag_researcher"
    elif next_worker == "code_dev":
        return "code_dev"
    elif next_worker == "data_analyst":
        # data_analyst 使用 tool_caller 实现
        return "tool_caller"
    elif next_worker == "critic":
        return "critic"
    else:
        # 默认结束
        return "finalize"


def critic_router(state: GlobalState) -> str:
    """
    Critic 路由：根据 is_approved 字段决定下一步
    
    - 如果通过审核：返回 supervisor 进行下一轮决策
    - 如果驳回：返回 last_worker 重新执行
    
    Returns:
        目标节点的名称
    """
    is_approved = state.get("is_approved", False)
    last_worker = state.get("last_worker", "")
    
    if is_approved:
        # 审核通过 -> 回 Supervisor 进行下一轮决策
        print("  --> 审核通过，返回 Supervisor")
        return "supervisor"
    else:
        # 审核驳回 -> 返回上一个 Worker 重做
        print(f"  --> 审核驳回，返回 {last_worker} 重做")
        
        # 映射 worker 名称到节点名称
        worker_to_node = {
            "rag_researcher": "rag_researcher",
            "code_dev": "code_dev",
            "data_analyst": "tool_caller",
            "tool_caller": "tool_caller"
        }
        
        return worker_to_node.get(last_worker, "supervisor")


def finalize_step(state: GlobalState) -> GlobalState:
    """
    最终化节点：整合所有结果，生成最终答案
    """
    print("\n=== [Finalize] 生成最终答案 ===")
    
    # 整合所有结果
    rag_context = state.get("rag_context", "")
    code_solution = state.get("code_solution", "")
    final_report = state.get("final_report", "")
    
    # 构建最终答案
    parts = []
    if rag_context:
        parts.append(f"【相关文献】\n{rag_context}")
    if code_solution:
        parts.append(f"【代码方案】\n{code_solution}")
    if final_report:
        parts.append(f"【分析报告】\n{final_report}")
    
    if not parts:
        final_answer = state.get("user_query", "任务已完成")
    else:
        final_answer = "\n\n".join(parts)
    
    return {
        **state,
        "final_answer": final_answer,
        "next_worker": "FINISH"
    }


# ==================== 构建主图 ====================

workflow = StateGraph(GlobalState)

# 1. 添加所有节点
# Supervisor 和 Critic 直接使用子图
workflow.add_node("supervisor", supervisor_agent_graph)
workflow.add_node("critic", critic_agent_graph)

# Worker 节点使用 wrapper 函数，确保更新 last_worker
workflow.add_node("rag_researcher", wrap_rag_researcher)
workflow.add_node("code_dev", wrap_code_dev)
workflow.add_node("tool_caller", wrap_tool_caller)  # data_analyst 使用 tool_caller

# Finalize 节点
workflow.add_node("finalize", finalize_step)

# 2. 定义流程

# START -> supervisor
workflow.add_edge(START, "supervisor")

# supervisor -> (路由) -> [rag_researcher, code_dev, tool_caller, critic, finalize]
workflow.add_conditional_edges(
    "supervisor",
    supervisor_router,
    {
        "rag_researcher": "rag_researcher",
        "code_dev": "code_dev",
        "tool_caller": "tool_caller",
        "critic": "critic",
        "finalize": "finalize"
    }
)

# 所有 Worker -> critic (所有产出必须经过审核)
workflow.add_edge("rag_researcher", "critic")
workflow.add_edge("code_dev", "critic")
workflow.add_edge("tool_caller", "critic")

# critic -> (路由) -> [supervisor OR Worker]
workflow.add_conditional_edges(
    "critic",
    critic_router,
    {
        "supervisor": "supervisor",
        "rag_researcher": "rag_researcher",
        "code_dev": "code_dev",
        "tool_caller": "tool_caller"
    }
)

# finalize -> END
workflow.add_edge("finalize", END)

# 3. 编译主图
graph = workflow.compile()

