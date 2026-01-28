"""
Tool Caller Agent 的 Prompt 定义
"""
from langchain_core.messages import SystemMessage
from .tools import AVAILABLE_TOOLS


def get_decision_system_prompt() -> SystemMessage:
    """
    动态构建决策节点的 System Prompt。
    它会自动读取 tools 里定义的所有工具描述。
    """
    
    # 自动生成工具描述文本
    tools_desc = "\n".join([
        f"- {t['name']}: {t['description']}" 
        for t in AVAILABLE_TOOLS
    ])

    content = (
        "You are an intelligent router agent.\n"
        "Analyze the user input and decide which tool to use.\n\n"
        "You are forbidden from interpret the user input by yourself. You can only select proper tools for it.\n"
        "Available Tools:\n"
        f"{tools_desc}\n\n"
        "OUTPUT FORMAT (Strict JSON):\n"
        "{\n"
        '  "tool_name": "name of the tool OR null if no tool needed",\n'
        '  "tool_args": { "key": "value" } \n'
        "}\n\n"
        "EXAMPLES:\n"
        'User: "Check xxx gene" -> {"tool_name": "mygene", "tool_args": {"gene_symbol": "xxx"}}\n'
        'User: "Hi" -> {"tool_name": null, "tool_args": {}}'
    )
    
    return SystemMessage(content=content)


def get_interpret_system_prompt() -> SystemMessage:
    """结果解释节点的 Prompt"""
    return SystemMessage(
        content=(
            "You are a helpful assistant.\n"
            "You will receive a user question and the raw result from a tool.\n"
            "Please summarize the tool result to answer the user's question professionally.\n"
            "Please use CHINESE in your interpretation."
        )
    )

