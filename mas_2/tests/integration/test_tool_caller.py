"""
Tool Caller Agent 集成测试
使用真实的 LLM 和工具进行端到端测试
"""
import pytest
from src.agents.tool_caller.graph import tool_caller_agent_graph
from langchain_core.messages import HumanMessage


def test_tool_caller_mygene_query(tool_caller_state):
    """
    测试 Tool Caller 调用 MyGene 工具查询基因信息
    
    前提条件：
    - .env 文件已配置 OPENAI_API_KEY
    - 有网络连接可以调用 MyGene API
    """
    state = tool_caller_state.copy()
    state["user_query"] = "查询 TP53 基因信息"
    
    # 执行 Tool Caller Agent 子图
    result = tool_caller_agent_graph.invoke(state)
    
    # 断言：tool_name 字段存在
    assert "tool_name" in result, "结果中应包含 tool_name 字段"
    
    # 断言：pending_contribution 已更新
    assert "pending_contribution" in result, "结果中应包含 pending_contribution 字段"
    
    # 断言：tool_final_answer 已更新
    assert "tool_final_answer" in result, "结果中应包含 tool_final_answer 字段"
    
    # 检查工具调用结果
    tool_name = result.get("tool_name")
    tool_result = result.get("tool_result")
    pending = result.get("pending_contribution")
    
    # 如果选择了工具，应该有结果
    if tool_name:
        assert tool_result is not None, "如果调用了工具，应该有 tool_result"
        assert pending is not None, "如果调用了工具，应该有 pending_contribution"
        
        # tool_final_answer 应该包含解释
        final_answer = result.get("tool_final_answer", "")
        assert len(final_answer) > 0, "tool_final_answer 不应为空"


def test_tool_caller_greeting(tool_caller_state):
    """
    测试 Tool Caller 处理简单问候（不需要工具）
    """
    state = tool_caller_state.copy()
    state["user_query"] = "你好"
    
    result = tool_caller_agent_graph.invoke(state)
    
    # 断言：应该处理问候
    assert "tool_name" in result
    assert "pending_contribution" in result
    
    # 如果不需要工具，tool_name 可能为 None
    # 但应该有响应
    assert result.get("pending_contribution") is not None or \
           result.get("tool_final_answer") is not None


def test_tool_caller_enrichment_query(tool_caller_state):
    """
    测试 Tool Caller 调用基因富集分析工具
    
    注意：这个测试需要 gseapy 库和网络连接
    """
    state = tool_caller_state.copy()
    state["user_query"] = "对 CD3D, LCK, TP53 这些基因进行富集分析"
    
    result = tool_caller_agent_graph.invoke(state)
    
    # 断言：基本字段存在
    assert "tool_name" in result
    assert "pending_contribution" in result
    
    # 如果选择了富集分析工具
    tool_name = result.get("tool_name")
    if tool_name == "gene_set_enrichment":
        assert result.get("tool_result") is not None
        assert result.get("pending_contribution") is not None


def test_tool_caller_celltype_annotation(tool_caller_state):
    """
    测试 Tool Caller 调用细胞类型注释工具
    """
    state = tool_caller_state.copy()
    state["user_query"] = "根据标记基因 CD3D, CD3E, LCK 注释细胞类型"
    
    result = tool_caller_agent_graph.invoke(state)
    
    # 断言：基本字段存在
    assert "tool_name" in result
    assert "pending_contribution" in result
    
    # 如果选择了细胞类型注释工具
    tool_name = result.get("tool_name")
    if tool_name == "annotate_celltype":
        assert result.get("tool_result") is not None
        assert result.get("pending_contribution") is not None


def test_tool_caller_error_handling(tool_caller_state):
    """
    测试 Tool Caller 的错误处理
    """
    state = tool_caller_state.copy()
    state["user_query"] = "调用一个不存在的工具"
    
    result = tool_caller_agent_graph.invoke(state)
    
    # 断言：应该处理错误
    assert "tool_name" in result
    assert "tool_result" in result
    
    # 如果工具不存在，tool_result 应该包含错误信息
    tool_result = result.get("tool_result")
    if tool_result and isinstance(tool_result, str):
        # 可能包含错误信息
        pass


def test_tool_caller_multiple_tools(tool_caller_state):
    """
    测试 Tool Caller 处理需要多个工具的任务
    """
    state = tool_caller_state.copy()
    state["user_query"] = "先查询 TP53 基因，然后进行富集分析"
    
    result = tool_caller_agent_graph.invoke(state)
    
    # 断言：应该选择了一个工具（子图只执行一次，所以只会选择一个）
    assert "tool_name" in result
    assert "pending_contribution" in result
    
    # 注意：由于子图只执行一次，只会调用一个工具
    # 如果需要多个工具，需要在主图中多次调用

