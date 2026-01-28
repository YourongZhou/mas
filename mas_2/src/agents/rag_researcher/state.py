"""
RAG Researcher Agent 状态定义
定义 RAG 检索相关的状态字段
"""
from typing import TypedDict, List
from src.core.state import GlobalState


class RAGAgentState(GlobalState):
    """
    RAG Researcher Agent 的状态
    继承自 GlobalState，添加 RAG 检索相关的字段
    """
    
    # === RAG 检索相关字段 ===
    # 检索查询
    search_query: str
    # 检索到的文档列表
    retrieved_docs: List[str]
    # 检索到的文档数量
    doc_count: int

