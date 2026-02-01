"""
Tool Caller Agent 单元测试
使用 Mock 模拟 LLM 调用和工具执行
"""
import pytest
import json
from unittest.mock import patch, MagicMock
from src.agents.tool_caller.graph import tool_caller_agent_graph
from src.agents.tool_caller.state import ToolCallerAgentState


@pytest.fixture
def base_state():
    """基础测试状态"""
    return {
        "user_query": "",
        "messages": [],
        "tool_name": None,
        "tool_args": {},
        "tool_result": None,
        "tool_final_answer": "",
        "pending_contribution": ""
    }


def test_tool_caller_flow_success(base_state):
    """
    测试完整流程：Decision -> Execute Tool -> Interpret
    """
    state = base_state.copy()
    state["user_query"] = "查一下 TP53 基因"

    # 1. 模拟 Decision 节点的 LLM 响应
    decision_content = json.dumps({
        "tool_name": "mock_tool",
        "tool_args": {"gene": "TP53"}
    })
    
    # 2. 模拟 Interpret 节点的 LLM 响应
    interpret_content = "TP53 是一个著名的抑癌基因。"

    # 3. 模拟工具函数的执行结果
    mock_tool_result = {"symbol": "TP53", "type": "coding"}

    # 创建一个模拟 LangChain BaseTool 的对象
    mock_tool_instance = MagicMock()
    mock_tool_instance.name = "mock_tool"
    mock_tool_instance.description = "测试工具"
    # 配置 invoke 方法的返回值
    mock_tool_instance.invoke.return_value = mock_tool_result

    # 构造 Mock 的 TOOL_MAP，里面放的是对象
    mock_tool_map = {
        "mock_tool": mock_tool_instance
    }

    # 开始 Patch
    with patch('src.agents.tool_caller.nodes.llm') as mock_llm, \
         patch('src.agents.tool_caller.nodes.TOOL_MAP', mock_tool_map):
        
        mock_msg_decision = MagicMock()
        mock_msg_decision.content = decision_content
        
        mock_msg_interpret = MagicMock()
        mock_msg_interpret.content = interpret_content
        
        mock_llm.invoke.side_effect = [mock_msg_decision, mock_msg_interpret]

        # 执行图
        result = tool_caller_agent_graph.invoke(state)

        # === 断言验证 ===
        assert result["tool_name"] == "mock_tool"
        assert result["tool_args"] == {"gene": "TP53"}
        
        # 验证结果正确
        assert result["tool_result"] == mock_tool_result
        
        # 验证 invoke 被调用了 (替代了以前的 func.called)
        mock_tool_instance.invoke.assert_called_once_with({"gene": "TP53"})
        
        assert result["tool_final_answer"] == interpret_content


def test_tool_caller_no_tool_needed(base_state):
    """
    测试不需要工具的流程：Decision -> Interpret (Direct Answer)
    """
    state = base_state.copy()
    state["user_query"] = "Hi"

    decision_content = json.dumps({
        "tool_name": None,
        "tool_args": {}
    })

    with patch('src.agents.tool_caller.nodes.llm') as mock_llm:
        mock_msg = MagicMock()
        mock_msg.content = decision_content
        mock_llm.invoke.return_value = mock_msg

        result = tool_caller_agent_graph.invoke(state)

        assert result["tool_name"] is None
        expected_answer = "I can help you search biological databases. Please ask about a gene."
        assert result["tool_final_answer"] == expected_answer
        assert mock_llm.invoke.call_count == 1


def test_decision_node_json_error(base_state):
    """
    测试 Decision 节点处理 JSON 解析错误的情况
    """
    state = base_state.copy()
    state["user_query"] = "Invalid JSON input"

    bad_content = "Thinking: I should use a tool but here is plain text."

    with patch('src.agents.tool_caller.nodes.llm') as mock_llm:
        mock_msg = MagicMock()
        mock_msg.content = bad_content
        mock_llm.invoke.return_value = mock_msg

        result = tool_caller_agent_graph.invoke(state)

        assert result["tool_name"] is None
        assert result["tool_args"] == {}
        assert "I can help you search" in result["tool_final_answer"]


def test_tool_execution_error(base_state):
    """
    测试工具执行抛出异常的情况
    """
    state = base_state.copy()
    state["user_query"] = "Break the tool"

    decision_content = json.dumps({
        "tool_name": "broken_tool",
        "tool_args": {}
    })
    
    interpret_content = "工具执行出错，请检查输入。"

    # 模拟一个会抛出异常的 Tool 对象
    mock_broken_tool = MagicMock()
    mock_broken_tool.name = "broken_tool"
    # 配置 invoke 方法抛出异常
    mock_broken_tool.invoke.side_effect = ValueError("Simulated API Crash")

    mock_tool_map = {
        "broken_tool": mock_broken_tool
    }

    with patch('src.agents.tool_caller.nodes.llm') as mock_llm, \
         patch('src.agents.tool_caller.nodes.TOOL_MAP', mock_tool_map):
        
        mock_msg_decision = MagicMock()
        mock_msg_decision.content = decision_content
        mock_msg_interpret = MagicMock()
        mock_msg_interpret.content = interpret_content
        
        mock_llm.invoke.side_effect = [mock_msg_decision, mock_msg_interpret]

        result = tool_caller_agent_graph.invoke(state)

        # 断言
        assert "Tool Execution Error" in str(result["tool_result"])
        assert "Simulated API Crash" in str(result["tool_result"])
        assert result["tool_final_answer"] == interpret_content


def test_tool_not_found(base_state):
    """
    测试 LLM 幻觉生成了不存在的 tool_name
    """
    state = base_state.copy()
    state["user_query"] = "Use magic tool"

    decision_content = json.dumps({
        "tool_name": "non_existent_tool",
        "tool_args": {}
    })
    
    interpret_content = "抱歉，找不到该工具。"

    # 空的工具表，模拟找不到工具的情况
    with patch('src.agents.tool_caller.nodes.llm') as mock_llm, \
         patch('src.agents.tool_caller.nodes.TOOL_MAP', {}): 
        
        mock_msg_decision = MagicMock()
        mock_msg_decision.content = decision_content
        mock_msg_interpret = MagicMock()
        mock_msg_interpret.content = interpret_content
        
        mock_llm.invoke.side_effect = [mock_msg_decision, mock_msg_interpret]

        result = tool_caller_agent_graph.invoke(state)

        assert "Error: Tool 'non_existent_tool' not found" in result["tool_result"]