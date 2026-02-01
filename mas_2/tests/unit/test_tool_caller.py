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
    模拟场景：用户查询基因，LLM 决定调用工具，工具返回结果，LLM 解释结果。
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

    # 定义 Mock 工具函数
    def mock_func(args):
        return mock_tool_result

    # 构造 Mock 的 TOOL_MAP
    mock_tool_map = {
        "mock_tool": {
            "name": "mock_tool",
            "description": "测试工具",
            "func": mock_func
        }
    }

    # 开始 Patch
    # 注意：我们要 Patch nodes.py 里的 llm 和 TOOL_MAP
    with patch('src.agents.tool_caller.nodes.llm') as mock_llm, \
         patch('src.agents.tool_caller.nodes.TOOL_MAP', mock_tool_map):
        
        # 配置 LLM Mock 的 side_effect (顺序返回)
        # 第一次调用是 Decision 节点，第二次调用是 Interpret 节点
        mock_msg_decision = MagicMock()
        mock_msg_decision.content = decision_content
        
        mock_msg_interpret = MagicMock()
        mock_msg_interpret.content = interpret_content
        
        mock_llm.invoke.side_effect = [mock_msg_decision, mock_msg_interpret]

        # 执行图
        result = tool_caller_agent_graph.invoke(state)

        # === 断言验证 ===
        
        # 1. 验证 Decision 节点是否正确写入了 tool_name 和 args
        assert result["tool_name"] == "mock_tool"
        assert result["tool_args"] == {"gene": "TP53"}
        
        # 2. 验证 Execute 节点是否正确运行了 Mock 工具并写入 result
        assert result["tool_result"] == mock_tool_result
        
        # 3. 验证 Interpret 节点是否生成了最终答案并写入 pending_contribution
        assert result["tool_final_answer"] == interpret_content
        assert result["pending_contribution"] == interpret_content
        
        # 4. 验证 LLM 被调用了两次
        assert mock_llm.invoke.call_count == 2


def test_tool_caller_no_tool_needed(base_state):
    """
    测试不需要工具的流程：Decision -> Interpret (Direct Answer)
    模拟场景：用户说 "Hi"，LLM 决定 tool_name: null
    """
    state = base_state.copy()
    state["user_query"] = "Hi"

    # 模拟 Decision 响应：不使用工具
    decision_content = json.dumps({
        "tool_name": None,
        "tool_args": {}
    })

    with patch('src.agents.tool_caller.nodes.llm') as mock_llm:
        mock_msg = MagicMock()
        mock_msg.content = decision_content
        mock_llm.invoke.return_value = mock_msg

        result = tool_caller_agent_graph.invoke(state)

        # 验证 Decision 结果
        assert result["tool_name"] is None
        
        # 验证 Interpret 节点逻辑：
        # 根据 nodes.py 中的 interpret_node，如果没有 tool_name，
        # 会直接返回固定话术，并且不调用 LLM 进行解释。
        expected_answer = "I can help you search biological databases. Please ask about a gene."
        assert result["tool_final_answer"] == expected_answer
        assert result["pending_contribution"] == expected_answer
        
        # 验证 LLM 只被调用了一次（Decision 节点）
        assert mock_llm.invoke.call_count == 1


def test_decision_node_json_error(base_state):
    """
    测试 Decision 节点处理 JSON 解析错误的情况
    """
    state = base_state.copy()
    state["user_query"] = "Invalid JSON input"

    # 模拟 LLM 返回非 JSON 内容
    bad_content = "Thinking: I should use a tool but here is plain text."

    with patch('src.agents.tool_caller.nodes.llm') as mock_llm:
        mock_msg = MagicMock()
        mock_msg.content = bad_content
        mock_llm.invoke.return_value = mock_msg

        result = tool_caller_agent_graph.invoke(state)

        # 断言：发生 JSON 错误时，nodes.py 中的 except 块应该捕获异常，
        # 并返回 tool_name: None
        assert result["tool_name"] is None
        assert result["tool_args"] == {}
        
        # 此时应该走入 interpret 的默认分支
        assert "I can help you search" in result["tool_final_answer"]


def test_tool_execution_error(base_state):
    """
    测试工具执行抛出异常的情况
    """
    state = base_state.copy()
    state["user_query"] = "Break the tool"

    # 1. Decision 决定调用工具
    decision_content = json.dumps({
        "tool_name": "broken_tool",
        "tool_args": {}
    })
    
    # 2. Interpret 尝试解释错误信息
    interpret_content = "工具执行出错，请检查输入。"

    # 定义一个会抛出异常的 Mock 工具
    def broken_func(args):
        raise ValueError("Simulated API Crash")

    mock_tool_map = {
        "broken_tool": {
            "name": "broken_tool",
            "description": "Broken Tool",
            "func": broken_func
        }
    }

    with patch('src.agents.tool_caller.nodes.llm') as mock_llm, \
         patch('src.agents.tool_caller.nodes.TOOL_MAP', mock_tool_map):
        
        mock_msg_decision = MagicMock()
        mock_msg_decision.content = decision_content
        
        mock_msg_interpret = MagicMock()
        mock_msg_interpret.content = interpret_content
        
        mock_llm.invoke.side_effect = [mock_msg_decision, mock_msg_interpret]

        result = tool_caller_agent_graph.invoke(state)

        # 断言：工具执行结果应该包含错误信息
        assert "Tool Execution Error" in str(result["tool_result"])
        assert "Simulated API Crash" in str(result["tool_result"])
        
        # 断言：Interpret 节点接收到了错误结果并生成了回答
        assert result["tool_final_answer"] == interpret_content


def test_tool_not_found(base_state):
    """
    测试 LLM 幻觉生成了不存在的 tool_name
    """
    state = base_state.copy()
    state["user_query"] = "Use magic tool"

    # Decision 决定调用一个不存在的工具
    decision_content = json.dumps({
        "tool_name": "non_existent_tool",
        "tool_args": {}
    })
    
    interpret_content = "抱歉，找不到该工具。"

    with patch('src.agents.tool_caller.nodes.llm') as mock_llm, \
         patch('src.agents.tool_caller.nodes.TOOL_MAP', {}): # 空的工具表
        
        mock_msg_decision = MagicMock()
        mock_msg_decision.content = decision_content
        
        mock_msg_interpret = MagicMock()
        mock_msg_interpret.content = interpret_content
        
        mock_llm.invoke.side_effect = [mock_msg_decision, mock_msg_interpret]

        result = tool_caller_agent_graph.invoke(state)

        # 断言：Execute 节点处理了找不到工具的情况
        assert "Error: Tool 'non_existent_tool' not found" in result["tool_result"]