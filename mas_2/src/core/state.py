"""
全局状态定义模块
定义 LangGraph 的状态结构
"""
from typing import TypedDict, List, Optional, Literal, Any, Annotated
try:
    from langchain_core.messages import BaseMessage
    from langgraph.graph.message import add_messages
except ImportError:
    # 兼容性处理：如果导入失败，使用备用方案
    from typing import Any as BaseMessage
    def add_messages(left: List, right: List) -> List:
        """简单的消息合并函数"""
        return left + right


class GlobalState(TypedDict):
    """
    全局状态定义
    
    包含消息历史、用户查询、调度控制、产出结果和交互反馈等字段
    """
    
    # === 基础字段 ===
    # 消息历史（使用 add_messages reducer）
    messages: Annotated[List[BaseMessage], add_messages]
    # 用户的原始查询
    user_query: str
    
    # === 调度控制字段 ===
    # 下一个要执行的 worker（由 Supervisor 决定）
    next_worker: Literal["rag_researcher", "code_dev", "data_analyst", "critic", "FINISH"]
    # 上一个执行的 worker（用于追踪）
    last_worker: str
    
    # === 产出字段 ===
    # 最终报告（整合所有结果）
    final_report: str
    # 代码解决方案
    code_solution: str
    # RAG 检索到的上下文
    rag_context: str
    
    # === 交互字段 ===
    # Worker 提交的待审核草稿（临时缓冲区）
    pending_contribution: Any
    # Critic 的审核反馈
    critique_feedback: Optional[str]
    # 是否通过审核
    is_approved: bool

