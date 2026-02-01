"""
Tool Caller Agent 的节点定义
"""
import json
from langchain_core.messages import HumanMessage
from .state import ToolCallerAgentState
from src.core.llm import get_llm
from .prompts import get_decision_system_prompt, get_interpret_system_prompt
from .tools import TOOL_MAP # 确保这导入的是包含 BaseTool 对象的字典

# 初始化 LLM
llm = get_llm(model_name="qwen-plus", temperature=0.1)


def decision_node(state: ToolCallerAgentState) -> ToolCallerAgentState:
    """
    [节点 1] 决策节点：决定调用哪个工具，准备什么参数
    """
    print("\n=== [Tool Caller - Decision Node] ===")
    user_input = state.get("user_query", "")
    
    # 回退策略：如果 user_query 为空，尝试从 messages 历史中取最后一条
    if not user_input and state.get("messages"):
        last_msg = state["messages"][-1]
        user_input = last_msg.content if hasattr(last_msg, "content") else str(last_msg)
    
    # 1. 获取动态生成的 Prompt
    system_msg = get_decision_system_prompt()
    human_msg = HumanMessage(content=user_input)

    # 2. 调用 LLM
    response = llm.invoke([system_msg, human_msg])
    
    # 3. 解析 JSON
    try:
        content = response.content.replace("```json", "").replace("```", "").strip()
        decision = json.loads(content)
    except Exception as e:
        print(f"JSON Parse Error: {e}")
        # 出错时回退到不做任何操作
        decision = {"tool_name": None, "tool_args": {}}

    print(f"Decision: Use {decision.get('tool_name')} with args {decision.get('tool_args')}")

    return {
        **state,
        "tool_name": decision.get("tool_name"),
        "tool_args": decision.get("tool_args", {})
    }


def tool_execution_node(state: ToolCallerAgentState) -> ToolCallerAgentState:
    """
    [节点 2] 通用工具执行节点
    根据 tool_name 动态查找并运行函数 (适配 LangChain BaseTool)
    """
    print("\n=== [Tool Caller - Tool Execution Node] ===")
    name = state.get("tool_name")
    args = state.get("tool_args", {})

    # 1. 从注册表中查找工具对象
    tool_instance = TOOL_MAP.get(name)
    
    if tool_instance:
        try:
            # 使用 LangChain 标准 invoke 方法执行
            # invoke 会自动处理参数验证
            result = tool_instance.invoke(args)
            print(f"Tool '{name}' executed successfully")
        except Exception as e:
            result = f"Tool Execution Error: {str(e)}"
            print(f"Tool execution failed: {e}")
    else:
        result = f"Error: Tool '{name}' not found."
        print(f"Tool '{name}' not found in TOOL_MAP")

    return {
        **state,
        "tool_result": result
    }


def interpret_node(state: ToolCallerAgentState) -> ToolCallerAgentState:
    """
    [节点 3] 结果解释节点
    """
    print("\n=== [Tool Caller - Interpret Node] ===")
    
    # 如果没有调用工具（比如只是打招呼），直接回答
    if not state.get("tool_name"):
        answer = "I can help you search biological databases. Please ask about a gene."
        return {
            **state,
            "tool_final_answer": answer,
            "pending_contribution": answer
        }

    system_msg = get_interpret_system_prompt()
    
    user_query = state.get("user_query", "")
    context_content = (
        f"User Input: {user_query}\n"
        f"Tool Name: {state.get('tool_name')}\n"
        f"Tool Result: {state.get('tool_result')}"
    )
    
    response = llm.invoke([system_msg, HumanMessage(content=context_content)])
    answer = response.content

    return {
        **state,
        "tool_final_answer": answer,
        "pending_contribution": answer
    }  