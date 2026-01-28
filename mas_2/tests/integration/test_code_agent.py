"""
Code Dev Agent 集成测试
使用真实的 LLM 进行端到端测试
"""
import pytest
from src.agents.code_dev.graph import code_agent_graph
from langchain_core.messages import HumanMessage


def test_code_agent_execution(code_dev_state):
    """
    测试 Code Dev Agent 的完整执行流程
    
    前提条件：
    - .env 文件已配置 OPENAI_API_KEY 等环境变量
    - 有网络连接可以调用 LLM API
    """
    # 准备测试状态
    state = code_dev_state.copy()
    state["user_query"] = "写一个函数 print hello"
    state["task"] = "写一个函数 print hello"
    
    # 执行 Code Dev Agent 子图
    result = code_agent_graph.invoke(state)
    
    # 断言：pending_contribution 字段存在且包含代码相关内容
    assert "pending_contribution" in result, "结果中应包含 pending_contribution 字段"
    
    pending = result.get("pending_contribution")
    assert pending is not None, "pending_contribution 不应为空"
    
    # 检查 pending_contribution 的内容
    # 它可能是一个字典（包含 code 字段）或字符串
    if isinstance(pending, dict):
        code = pending.get("code", "")
        assert "print" in code.lower() or "hello" in code.lower(), \
            f"生成的代码应包含 'print' 或 'hello'，实际内容: {code[:200]}"
    else:
        # 如果是字符串，直接检查
        pending_str = str(pending)
        assert "print" in pending_str.lower() or "hello" in pending_str.lower(), \
            f"pending_contribution 应包含 'print' 或 'hello'，实际内容: {pending_str[:200]}"
    
    # 断言：scanpy_code 字段已更新
    assert "scanpy_code" in result, "结果中应包含 scanpy_code 字段"
    assert result.get("scanpy_code"), "scanpy_code 不应为空"
    
    # 断言：internal_iteration_count 已增加
    assert result.get("internal_iteration_count", 0) > 0, \
        "internal_iteration_count 应该增加"


def test_code_agent_with_feedback(code_dev_state):
    """
    测试 Code Dev Agent 在收到反馈后的重试逻辑
    """
    state = code_dev_state.copy()
    state["user_query"] = "写一个函数计算两数之和"
    state["task"] = "写一个函数计算两数之和"
    state["feedback"] = "请添加类型注解"
    state["internal_iteration_count"] = 1
    
    # 执行 Code Dev Agent
    result = code_agent_graph.invoke(state)
    
    # 断言：应该生成了代码
    assert result.get("pending_contribution") is not None
    assert result.get("internal_iteration_count", 0) > 1

