"""
Supervisor Agent 状态定义
定义调度决策相关的状态字段
"""
from typing import TypedDict
from src.core.state import GlobalState


class SupervisorAgentState(GlobalState):
    """
    Supervisor Agent 的状态
    继承自 GlobalState，用于调度决策
    """
    # Supervisor 不需要额外的私有字段
    # 直接使用 GlobalState 的字段即可
    pass

