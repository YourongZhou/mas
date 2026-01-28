"""
RAG Researcher Agent 集成测试
使用真实的向量库进行端到端测试

注意：此测试需要配置 ChromaDB 和向量库
如果向量库未配置，测试可能会返回默认消息
"""
import pytest
import os
from src.agents.rag_researcher.graph import rag_agent_graph
from langchain_core.messages import HumanMessage


@pytest.mark.skipif(
    not os.getenv("CHROMA_PERSIST_PATH") and not os.path.exists("./chroma_db"),
    reason="需要配置 ChromaDB 向量库"
)
def test_rag_researcher_search(rag_researcher_state):
    """
    测试 RAG Researcher 检索文档
    
    前提条件：
    - ChromaDB 向量库已配置
    - 向量库中已有文档数据
    """
    state = rag_researcher_state.copy()
    state["user_query"] = "基因功能"
    state["search_query"] = "基因功能"
    
    # 执行 RAG Researcher Agent 子图
    result = rag_agent_graph.invoke(state)
    
    # 断言：retrieved_docs 字段存在
    assert "retrieved_docs" in result, "结果中应包含 retrieved_docs 字段"
    
    # 断言：doc_count 字段存在
    assert "doc_count" in result, "结果中应包含 doc_count 字段"
    
    # 断言：pending_contribution 已更新
    assert "pending_contribution" in result, "结果中应包含 pending_contribution 字段"
    
    # 断言：rag_context 已更新
    assert "rag_context" in result, "结果中应包含 rag_context 字段"
    
    # 检查检索结果
    docs = result.get("retrieved_docs", [])
    doc_count = result.get("doc_count", 0)
    
    assert isinstance(docs, list), "retrieved_docs 应该是列表"
    assert doc_count == len(docs), "doc_count 应该等于 retrieved_docs 的长度"
    
    # 如果向量库有数据，应该有检索结果
    # 如果没有数据，可能会返回默认消息
    if len(docs) > 0:
        assert all(isinstance(doc, str) for doc in docs), "所有文档应该是字符串"
        assert result.get("rag_context"), "如果有文档，rag_context 不应为空"


def test_rag_researcher_empty_query(rag_researcher_state):
    """
    测试 RAG Researcher 处理空查询
    """
    state = rag_researcher_state.copy()
    state["user_query"] = ""
    state["search_query"] = ""
    
    result = rag_agent_graph.invoke(state)
    
    # 断言：应该处理空查询
    assert "retrieved_docs" in result
    assert "doc_count" in result
    assert "pending_contribution" in result


def test_rag_researcher_without_vector_db(rag_researcher_state):
    """
    测试 RAG Researcher 在没有向量库时的行为
    
    这个测试验证错误处理逻辑
    """
    # 临时修改环境变量，指向不存在的路径
    original_path = os.environ.get("CHROMA_PERSIST_PATH")
    os.environ["CHROMA_PERSIST_PATH"] = "/tmp/nonexistent_chroma_db_test"
    
    try:
        state = rag_researcher_state.copy()
        state["user_query"] = "测试查询"
        state["search_query"] = "测试查询"
        
        result = rag_agent_graph.invoke(state)
        
        # 断言：即使向量库不存在，也应该返回结果
        assert "retrieved_docs" in result
        assert "pending_contribution" in result
        
        # 可能会返回错误消息或空列表
        docs = result.get("retrieved_docs", [])
        assert isinstance(docs, list)
        
    finally:
        # 恢复原始环境变量
        if original_path:
            os.environ["CHROMA_PERSIST_PATH"] = original_path
        elif "CHROMA_PERSIST_PATH" in os.environ:
            del os.environ["CHROMA_PERSIST_PATH"]


def test_rag_researcher_specific_query(rag_researcher_state):
    """
    测试 RAG Researcher 处理特定查询
    """
    state = rag_researcher_state.copy()
    state["user_query"] = "TP53 基因突变"
    state["search_query"] = "TP53 基因突变"
    
    result = rag_agent_graph.invoke(state)
    
    # 断言：基本字段存在
    assert "retrieved_docs" in result
    assert "doc_count" in result
    assert "pending_contribution" in result
    assert "rag_context" in result
    
    # 验证结果格式
    assert isinstance(result.get("retrieved_docs"), list)
    assert isinstance(result.get("doc_count"), int)

