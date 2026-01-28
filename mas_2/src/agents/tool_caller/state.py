"""
Tool Caller Agent 状态定义
定义工具调用相关的状态字段
"""
from typing import TypedDict, Optional, Dict, Any
from src.core.state import GlobalState


class ToolCallerAgentState(GlobalState):
    """
    Tool Caller Agent 的状态
    继承自 GlobalState，添加工具调用相关的私有字段
    """
    
    # === 工具调用相关字段（内部使用）===
    # 决定使用的工具名称，如果为 None 则表示不使用工具
    tool_name: Optional[str]
    # 传递给工具的参数字典
    tool_args: Dict[str, Any]
    # 工具运行后的原始结果
    tool_result: Any
    # 工具调用的最终答案（解释后的结果）
    tool_final_answer: str

