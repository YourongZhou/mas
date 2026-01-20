from langgraph.graph import StateGraph, START, END
from src.schema import GlobalState

# 导入所有节点
from src.agents.supervisor import supervisor_node
from src.agents.critic import critic_node
from src.agents.rag import rag_app as rag_worker
from src.agents.code import code_app as code_worker
from src.agents.tool import react_graph as tool_worker

# --- 路由逻辑 ---

# 1. 经理路由：决定是否结束，或者去哪个 Worker
def supervisor_router(state: GlobalState):
    if state["next_action"] == "FINISH":
        return "finalize"
    return state["active_worker"] # 去指定的 Worker

# 2. 审核路由：决定是回经理汇报，还是打回给工人
def critic_router(state: GlobalState):
    if state["is_approved"]:
        # 通过 -> 回 Supervisor 进行下一轮决策
        return "supervisor"
    else:
        # 驳回 -> 谁干的就回谁那里去重做
        return state["active_worker"]

def finalize_step(state: GlobalState):
    return {"final_answer": state.get("code_result", "任务结束")}

# --- 构建图 ---
workflow = StateGraph(GlobalState)

# 添加节点
workflow.add_node("supervisor", supervisor_node)
workflow.add_node("critic", critic_node)
workflow.add_node("rag_agent", rag_worker)
workflow.add_node("code_agent", code_worker)
workflow.add_node("tool_agent", tool_worker)
workflow.add_node("finalize", finalize_step)

# 连线：启动 -> 经理
workflow.add_edge(START, "supervisor")

# 连线：经理 -> (选择工人)
workflow.add_conditional_edges(
    "supervisor",
    supervisor_router,
    {
        "rag_agent": "rag_agent",
        "code_agent": "code_agent",
        "tool_agent": "tool_agent",
        "finalize": "finalize"
    }
)

# 连线：工人 -> 审核 (所有产出必须过审核)
workflow.add_edge("rag_agent", "critic")
workflow.add_edge("code_agent", "critic")
workflow.add_edge("tool_agent", "critic")

# 连线：审核 -> (通过回经理 / 驳回回工人)
workflow.add_conditional_edges(
    "critic",
    critic_router,
    {
        "supervisor": "supervisor", # 只要通过了，就让经理看下一步怎么走
        "rag_agent": "rag_agent",
        "code_agent": "code_agent",
        "tool_agent": "tool_agent"
    }
)

workflow.add_edge("finalize", END)

graph = workflow.compile()