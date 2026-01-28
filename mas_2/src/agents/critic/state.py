"""
Critic Agent 状态定义
定义审核相关的状态字段
"""
from typing import TypedDict, Optional
from src.core.state import GlobalState


class CriticAgentState(GlobalState):
    """
    Critic Agent 的状态
    继承自 GlobalState，添加审核相关的字段
    """
    
    # === 审核相关字段 ===
    # 待审核的内容类型
    content_type: str  # "code", "docs", "db_result", "image"
    # 审核结果详情
    review_details: Optional[str]

