from typing import TypedDict, Any, Dict, Optional

class AgentState(TypedDict):
    """
    定义 Agent 的状态流转数据结构。
    使用 TypedDict 确保类型安全。
    """
    
    # 用户的原始输入
    user_input: str

    # === 决策层输出 ===
    # 决定使用的工具名称，如果为 None 则表示不使用工具
    tool_name: Optional[str]
    
    # 传递给工具的参数字典 (例如: {"gene_symbol": "TP53"} 或 {"a": 1, "b": 2})
    # 这样设计比把具体参数写死在 State 里更灵活，方便扩展
    tool_args: Dict[str, Any]

    # === 执行层输出 ===
    # 工具运行后的原始结果 (可能是 dict, str, int 等)
    tool_result: Any

    # === 最终输出 ===
    # 给用户的最终自然语言回复
    final_answer: str