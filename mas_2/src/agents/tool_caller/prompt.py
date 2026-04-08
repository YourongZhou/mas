"""
Tool Caller Agent 的 Prompt 定义
"""
from langchain_core.messages import SystemMessage
from .tools import AVAILABLE_TOOLS 

def get_tool_schema(tool) -> str:
    """辅助函数：获取工具的名称、描述和参数结构"""
    # 获取参数 schema (例如: {'gene_symbol': {'title': 'Gene Symbol', 'type': 'string'}})
    args_schema = tool.args
    
    # 格式化参数描述，让 LLM 知道具体的参数名
    args_str = ", ".join([f"{name}: {info.get('type', 'any')}" for name, info in args_schema.items()])
    
    return f"- {tool.name}: {tool.description} (Args: {args_str})"

def get_decision_system_prompt() -> SystemMessage:
    """
    动态构建决策节点的 System Prompt。
    自动读取 LangChain Tool 对象的元数据。
    """
    
    # 自动生成工具描述文本 (适配 BaseTool 对象)
    tools_desc = "\n".join([get_tool_schema(t) for t in AVAILABLE_TOOLS])

    content = (
        "You are an intelligent router agent.\n"
        "Analyze the user input and decide which tool to use.\n\n"
        "You are forbidden from interpret the user input by yourself. You can only select proper tools for it.\n"
        "Available Tools:\n"
        f"{tools_desc}\n\n"
        "OUTPUT FORMAT (Strict JSON):\n"
        "{\n"
        '  "tool_name": "name of the tool OR null if no tool needed",\n'
        '  "tool_args": { "argument_name": "value" } \n'
        "}\n\n"
        "EXAMPLES:\n"
        'User: "Check TP53 gene" -> {"tool_name": "query_mygene", "tool_args": {"gene_symbol": "TP53"}}\n'
        'User: "Hi" -> {"tool_name": null, "tool_args": {}}'
    )
    
    return SystemMessage(content=content)


def get_interpret_system_prompt() -> SystemMessage:
    """结果解释节点的 Prompt """
    return SystemMessage(
        content=(
            "You are a helpful assistant.\n"
            "You will receive a user question and the raw result from a tool.\n"
            "Please summarize the tool result to answer the user's question professionally.\n"
            "Please use CHINESE in your interpretation."
        )
    )