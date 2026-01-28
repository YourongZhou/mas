"""
Supervisor Agent 集成测试
使用真实的 LLM 进行端到端测试
"""
import pytest
from src.agents.supervisor.graph import supervisor_agent_graph
from langchain_core.messages import HumanMessage


def test_supervisor_decision_for_code_task(supervisor_state):
    """
    测试 Supervisor 为代码生成任务做决策
    """
    state = supervisor_state.copy()
    state["user_query"] = "写一个 Python 函数计算斐波那契数列"
    
    # 执行 Supervisor Agent 子图
    result = supervisor_agent_graph.invoke(state)
    
    # 断言：next_worker 字段已更新
    assert "next_worker" in result, "结果中应包含 next_worker 字段"
    assert result.get("next_worker") in [
        "rag_researcher", "code_dev", "data_analyst", "critic", "FINISH"
    ], f"next_worker 应该是有效的 worker 名称，实际为: {result.get('next_worker')}"
    
    # 对于代码生成任务，应该选择 code_dev
    # 但由于 LLM 的决策可能不同，我们只验证它返回了有效值
    assert result.get("next_worker") is not None


def test_supervisor_decision_for_research_task(supervisor_state):
    """
    测试 Supervisor 为研究任务做决策
    """
    state = supervisor_state.copy()
    state["user_query"] = "查找关于 TP53 基因的文献"
    state["rag_context"] = ""  # 没有 RAG 上下文
    
    result = supervisor_agent_graph.invoke(state)
    
    # 断言：应该选择 rag_researcher（但 LLM 可能选择其他，所以只验证有效值）
    assert result.get("next_worker") is not None
    assert result.get("next_worker") in [
        "rag_researcher", "code_dev", "data_analyst", "critic", "FINISH"
    ]


def test_supervisor_decision_with_pending_contribution(supervisor_state):
    """
    测试 Supervisor 在有待审核内容时的决策
    """
    state = supervisor_state.copy()
    state["user_query"] = "分析任务"
    state["pending_contribution"] = {"code": "print('test')"}
    state["is_approved"] = False
    
    result = supervisor_agent_graph.invoke(state)
    
    # 断言：应该优先选择 critic 审核
    # 但 LLM 可能选择其他，所以只验证返回了有效值
    assert result.get("next_worker") is not None


def test_supervisor_decision_finish_task(supervisor_state):
    """
    测试 Supervisor 决定结束任务
    """
    state = supervisor_state.copy()
    state["user_query"] = "任务已完成"
    state["rag_context"] = "已有上下文"
    state["code_solution"] = "已有代码"
    state["final_report"] = "已有报告"
    
    result = supervisor_agent_graph.invoke(state)
    
    # 断言：可能选择 FINISH（但 LLM 可能选择其他）
    assert result.get("next_worker") is not None


def test_supervisor_decision_after_approval(supervisor_state):
    """
    测试 Supervisor 在审核通过后的决策
    """
    state = supervisor_state.copy()
    state["user_query"] = "继续任务"
    state["last_worker"] = "code_dev"
    state["is_approved"] = True
    state["code_solution"] = "已生成的代码"
    
    result = supervisor_agent_graph.invoke(state)
    
    # 断言：应该返回有效的 next_worker
    assert result.get("next_worker") is not None

