"""
Pytest 配置文件
定义全局 fixtures
"""
import pytest
from langchain_core.messages import HumanMessage
from typing import Dict, Any


@pytest.fixture
def basic_global_state() -> Dict[str, Any]:
    """
    返回一个初始化好的 GlobalState 字典
    
    这个 fixture 可以在所有测试中使用，避免重复定义空字段
    """
    return {
        # === 基础字段 ===
        "messages": [HumanMessage(content="测试查询")],
        "user_query": "测试查询",
        
        # === 调度控制字段 ===
        "next_worker": "rag_researcher",
        "last_worker": "",
        
        # === 产出字段 ===
        "final_report": "",
        "code_solution": "",
        "rag_context": "",
        
        # === 交互字段 ===
        "pending_contribution": None,
        "critique_feedback": None,
        "is_approved": False,
    }


@pytest.fixture
def code_dev_state() -> Dict[str, Any]:
    """
    专门用于 Code Dev Agent 测试的状态
    """
    return {
        "messages": [HumanMessage(content="写一个函数 print hello")],
        "user_query": "写一个函数 print hello",
        "next_worker": "code_dev",
        "last_worker": "",
        "final_report": "",
        "code_solution": "",
        "rag_context": "",
        "pending_contribution": None,
        "critique_feedback": None,
        "is_approved": False,
        # Code Dev 特有字段
        "task": "写一个函数 print hello",
        "feedback": None,
        "data_path": "",
        "result_path": "./result",
        "scanpy_code": "",
        "requirements_txt": "",
        "analysis_result": "",
        "success": False,
        "internal_iteration_count": 0,
    }


@pytest.fixture
def supervisor_state() -> Dict[str, Any]:
    """
    专门用于 Supervisor Agent 测试的状态
    """
    return {
        "messages": [HumanMessage(content="分析这个任务")],
        "user_query": "分析这个任务",
        "next_worker": "rag_researcher",
        "last_worker": "",
        "final_report": "",
        "code_solution": "",
        "rag_context": "",
        "pending_contribution": None,
        "critique_feedback": None,
        "is_approved": False,
    }


@pytest.fixture
def rag_researcher_state() -> Dict[str, Any]:
    """
    专门用于 RAG Researcher Agent 测试的状态
    """
    return {
        "messages": [HumanMessage(content="查询基因相关信息")],
        "user_query": "查询基因相关信息",
        "search_query": "查询基因相关信息",
        "next_worker": "rag_researcher",
        "last_worker": "",
        "final_report": "",
        "code_solution": "",
        "rag_context": "",
        "pending_contribution": None,
        "critique_feedback": None,
        "is_approved": False,
        # RAG Researcher 特有字段
        "retrieved_docs": [],
        "doc_count": 0,
    }


@pytest.fixture
def critic_state() -> Dict[str, Any]:
    """
    专门用于 Critic Agent 测试的状态
    """
    return {
        "messages": [HumanMessage(content="审核工作成果")],
        "user_query": "测试查询",
        "next_worker": "critic",
        "last_worker": "code_dev",
        "final_report": "",
        "code_solution": "",
        "rag_context": "",
        "pending_contribution": {"code": "print('hello')"},
        "critique_feedback": None,
        "is_approved": False,
        # Critic 特有字段
        "content_type": "code",
        "review_details": None,
    }


@pytest.fixture
def tool_caller_state() -> Dict[str, Any]:
    """
    专门用于 Tool Caller Agent 测试的状态
    """
    return {
        "messages": [HumanMessage(content="查询 TP53 基因信息")],
        "user_query": "查询 TP53 基因信息",
        "next_worker": "data_analyst",
        "last_worker": "",
        "final_report": "",
        "code_solution": "",
        "rag_context": "",
        "pending_contribution": None,
        "critique_feedback": None,
        "is_approved": False,
        # Tool Caller 特有字段
        "tool_name": None,
        "tool_args": {},
        "tool_result": None,
        "tool_final_answer": "",
    }

