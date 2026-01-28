"""
Tool Caller Agent 子图
负责工具调用和结果解释

这个子图可以直接接收 GlobalState，因为 ToolCallerAgentState 继承自 GlobalState。
子图内部会使用 tool_name, tool_args, tool_result 等字段进行工具调用，
最终将结果写入 pending_contribution 字段，供 Critic 审核。
"""
from langgraph.graph import StateGraph, START, END
from .state import ToolCallerAgentState
from .nodes import decision_node, tool_execution_node, interpret_node


def route_logic(state: ToolCallerAgentState) -> str:
    """
    条件边逻辑：
    - 如果 decision 节点决定了 tool_name，则去 'execute_tool'
    - 如果没有 tool_name (null)，则直接去 'interpret'
    """
    if state.get("tool_name"):
        return "execute_tool"
    else:
        return "interpret"


# === 构建子图 ===
builder = StateGraph(ToolCallerAgentState)

# 1. 添加节点
builder.add_node("decision", decision_node)
builder.add_node("execute_tool", tool_execution_node)
builder.add_node("interpret", interpret_node)

# 2. 定义流程
builder.add_edge(START, "decision")

# 3. 添加条件分支
builder.add_conditional_edges(
    "decision",
    route_logic,
    {
        "execute_tool": "execute_tool",
        "interpret": "interpret"
    }
)

# 4. 闭环
builder.add_edge("execute_tool", "interpret")
builder.add_edge("interpret", END)

# 5. 编译子图
tool_caller_agent_graph = builder.compile()

