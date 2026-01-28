"""
Critic Agent 集成测试
使用真实的 LLM 进行端到端测试
"""
import pytest
from src.agents.critic.graph import critic_agent_graph
from langchain_core.messages import HumanMessage


def test_critic_review_code(critic_state):
    """
    测试 Critic 审核代码
    """
    state = critic_state.copy()
    state["last_worker"] = "code_dev"
    state["pending_contribution"] = {
        "code": "def hello():\n    print('Hello, World!')"
    }
    state["user_query"] = "写一个打印 hello 的函数"
    
    # 执行 Critic Agent 子图
    result = critic_agent_graph.invoke(state)
    
    # 断言：is_approved 字段存在
    assert "is_approved" in result, "结果中应包含 is_approved 字段"
    assert isinstance(result.get("is_approved"), bool), "is_approved 应该是布尔值"
    
    # 断言：review_details 字段存在
    assert "review_details" in result, "结果中应包含 review_details 字段"
    
    # 断言：content_type 字段存在
    assert "content_type" in result, "结果中应包含 content_type 字段"
    
    # 如果审核通过，应该归档到 code_solution
    if result.get("is_approved"):
        assert result.get("code_solution") is not None or result.get("pending_contribution") is None


def test_critic_review_docs(critic_state):
    """
    测试 Critic 审核文档列表
    """
    state = critic_state.copy()
    state["last_worker"] = "rag_researcher"
    state["pending_contribution"] = [
        "文档1：关于基因的研究",
        "文档2：基因功能分析",
        "文档3：相关文献摘要"
    ]
    state["user_query"] = "查找基因相关文献"
    
    result = critic_agent_graph.invoke(state)
    
    # 断言：基本字段存在
    assert "is_approved" in result
    assert "review_details" in result
    assert "content_type" in result
    
    # 如果审核通过，应该归档到 rag_context
    if result.get("is_approved"):
        assert result.get("rag_context") or result.get("pending_contribution") is None


def test_critic_review_tool_result(critic_state):
    """
    测试 Critic 审核工具调用结果
    """
    state = critic_state.copy()
    state["last_worker"] = "data_analyst"
    state["pending_contribution"] = "TP53 基因信息：位于17号染色体，编码肿瘤抑制蛋白"
    state["user_query"] = "查询 TP53 基因信息"
    
    result = critic_agent_graph.invoke(state)
    
    # 断言：基本字段存在
    assert "is_approved" in result
    assert "review_details" in result
    assert "content_type" in result


def test_critic_reject_with_feedback(critic_state):
    """
    测试 Critic 驳回并给出反馈
    """
    state = critic_state.copy()
    state["last_worker"] = "code_dev"
    state["pending_contribution"] = {
        "code": "print hello"  # 不完整的代码
    }
    state["user_query"] = "写一个完整的函数"
    
    result = critic_agent_graph.invoke(state)
    
    # 断言：如果驳回，应该有反馈
    if not result.get("is_approved"):
        assert result.get("critique_feedback") is not None, \
            "驳回时应该提供反馈"
        assert len(result.get("critique_feedback", "")) > 0, \
            "反馈不应为空"


def test_critic_review_empty_contribution(critic_state):
    """
    测试 Critic 审核空内容
    """
    state = critic_state.copy()
    state["last_worker"] = "code_dev"
    state["pending_contribution"] = None
    state["user_query"] = "测试查询"
    
    result = critic_agent_graph.invoke(state)
    
    # 断言：应该处理空内容
    assert "is_approved" in result
    # 空内容应该被驳回
    assert result.get("is_approved") == False or result.get("critique_feedback") is not None

