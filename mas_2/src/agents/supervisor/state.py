"""
Supervisor Agent 状态定义
定义调度决策相关的状态字段
"""
from typing import TypedDict, Optional
from src.core.state import GlobalState


class SupervisorAgentState(GlobalState):
    """
    Supervisor Agent 的状态
    继承自 GlobalState，用于调度决策
    """
    # 从初步探查节点提取出的适配系统当前任务的技能ID
    selected_skill_id: Optional[str]

