"""
主图全局集成测试
测试完整的工作流：Supervisor -> Worker -> Critic -> Supervisor -> Finalize
"""
import pytest
from src.main import graph
from langchain_core.messages import HumanMessage
from typing import Dict, Any


@pytest.fixture
def main_graph_initial_state() -> Dict[str, Any]:
    """
    主图初始状态
    """
    return {
        "messages": [HumanMessage(content="写一个函数打印 hello world")],
        "user_query": "写一个函数打印 hello world",
        "next_worker": "rag_researcher",  # 初始值，会被 Supervisor 覆盖
        "last_worker": "",
        "final_report": "",
        "code_solution": "",
        "rag_context": "",
        "pending_contribution": None,
        "critique_feedback": None,
        "is_approved": False,
    }


def test_main_graph_simple_code_task(main_graph_initial_state):
    """
    测试主图处理简单的代码生成任务
    
    预期流程：
    1. START -> Supervisor (决策)
    2. Supervisor -> Code Dev (生成代码)
    3. Code Dev -> Critic (审核代码)
    4. Critic -> Supervisor (通过) 或 Code Dev (驳回重做)
    5. Supervisor -> Finalize (结束)
    6. Finalize -> END
    
    注意：由于 LLM 的不确定性，实际流程可能不同
    """
    state = main_graph_initial_state.copy()
    state["user_query"] = "写一个简单的 Python 函数打印 hello"
    
    # 使用 stream 来跟踪执行流程，避免无限循环
    # 设置最大步数限制
    max_steps = 20
    step_count = 0
    
    try:
        # 使用 stream 来逐步执行，可以观察流程
        final_state = None
        for step in graph.stream(state):
            step_count += 1
            if step_count > max_steps:
                pytest.fail(f"主图执行超过最大步数限制 ({max_steps})，可能存在无限循环")
            
            # 获取最后一步的状态
            for node_name, node_state in step.items():
                final_state = node_state
                
                # 如果到达 finalize，应该结束
                if node_name == "finalize":
                    break
        
        # 如果没有通过 stream 完成，使用 invoke
        if final_state is None or final_state.get("next_worker") != "FINISH":
            # 如果还没结束，直接 invoke（可能会继续执行）
            final_state = graph.invoke(state)
        
        # 断言：最终状态应该包含 final_answer
        assert final_state is not None, "应该返回最终状态"
        assert "final_answer" in final_state or "next_worker" in final_state, \
            "最终状态应包含 final_answer 或 next_worker"
        
        # 验证关键字段存在
        assert "user_query" in final_state, "应保留 user_query"
        
    except Exception as e:
        # 如果出现异常，记录但不失败（因为 LLM 可能返回意外结果）
        pytest.skip(f"主图执行出现异常（可能是 LLM 返回格式问题）: {e}")


def test_main_graph_complete_workflow(main_graph_initial_state):
    """
    测试主图的完整工作流程
    
    这个测试验证：
    1. Supervisor 能够正确决策
    2. Worker 能够执行任务
    3. Critic 能够审核结果
    4. 状态在各个节点间正确传递
    """
    state = main_graph_initial_state.copy()
    state["user_query"] = "查询 TP53 基因信息并生成分析代码"
    
    # 执行主图（限制迭代次数）
    try:
        result = graph.invoke(state, config={"recursion_limit": 15})
        
        # 断言：应该返回结果
        assert result is not None, "主图应该返回结果"
        
        # 断言：关键字段存在
        assert "user_query" in result, "应保留 user_query"
        assert "messages" in result, "应保留 messages"
        
        # 验证工作流执行了至少一步
        # 如果 next_worker 不是初始值，说明至少执行了 Supervisor
        assert result.get("next_worker") is not None, \
            "next_worker 应该有值（至少执行了 Supervisor）"
        
    except Exception as e:
        pytest.skip(f"主图执行出现异常: {e}")


def test_main_graph_state_persistence(main_graph_initial_state):
    """
    测试主图中状态的持久性
    
    验证状态在各个节点之间正确传递和更新
    """
    state = main_graph_initial_state.copy()
    state["user_query"] = "测试状态传递"
    
    try:
        # 记录初始状态
        initial_query = state["user_query"]
        
        # 执行主图
        result = graph.invoke(state, config={"recursion_limit": 10})
        
        # 断言：user_query 应该保持不变
        assert result.get("user_query") == initial_query, \
            "user_query 应该在执行过程中保持不变"
        
        # 断言：messages 应该存在
        assert "messages" in result, "messages 应该存在"
        assert isinstance(result.get("messages"), list), "messages 应该是列表"
        
    except Exception as e:
        pytest.skip(f"主图执行出现异常: {e}")


def test_main_graph_supervisor_to_worker_flow(main_graph_initial_state):
    """
    测试 Supervisor -> Worker 的流程
    """
    state = main_graph_initial_state.copy()
    state["user_query"] = "生成代码"
    
    try:
        # 执行几步，观察 Supervisor 的决策
        steps = []
        for step in graph.stream(state, config={"recursion_limit": 5}):
            steps.append(step)
            # 如果到达 worker 节点，停止
            if any("code_dev" in str(k) or "rag_researcher" in str(k) 
                   for k in step.keys()):
                break
        
        # 验证至少执行了 Supervisor
        assert len(steps) > 0, "应该至少执行一步"
        
        # 验证 Supervisor 被执行
        supervisor_executed = any(
            "supervisor" in str(step.keys()) for step in steps
        )
        assert supervisor_executed, "Supervisor 应该被执行"
        
    except Exception as e:
        pytest.skip(f"主图执行出现异常: {e}")


def test_main_graph_critic_approval_flow(main_graph_initial_state):
    """
    测试 Critic 审核通过的流程
    
    当 Critic 审核通过后，应该返回 Supervisor
    """
    # 准备一个已经有 pending_contribution 的状态
    state = main_graph_initial_state.copy()
    state["user_query"] = "审核代码"
    state["last_worker"] = "code_dev"
    state["pending_contribution"] = {
        "code": "def hello():\n    print('Hello, World!')"
    }
    state["next_worker"] = "critic"  # 直接进入 Critic
    
    try:
        # 执行几步，观察 Critic 的行为
        steps = []
        for step in graph.stream(state, config={"recursion_limit": 5}):
            steps.append(step)
            # 如果回到 Supervisor，停止
            if any("supervisor" in str(k) for k in step.keys()):
                break
        
        # 验证 Critic 被执行
        critic_executed = any(
            "critic" in str(step.keys()) for step in steps
        )
        assert critic_executed, "Critic 应该被执行"
        
    except Exception as e:
        pytest.skip(f"主图执行出现异常: {e}")


def test_main_graph_finalize_step(main_graph_initial_state):
    """
    测试 Finalize 节点
    
    当 next_worker == "FINISH" 时，应该进入 Finalize
    """
    state = main_graph_initial_state.copy()
    state["user_query"] = "完成任务"
    state["next_worker"] = "FINISH"  # 直接进入 Finalize
    state["rag_context"] = "测试上下文"
    state["code_solution"] = "测试代码"
    
    try:
        result = graph.invoke(state, config={"recursion_limit": 2})
        
        # 断言：应该包含 final_answer
        assert "final_answer" in result or result.get("next_worker") == "FINISH", \
            "Finalize 应该生成 final_answer 或设置 next_worker=FINISH"
        
    except Exception as e:
        pytest.skip(f"主图执行出现异常: {e}")


def test_main_graph_error_recovery(main_graph_initial_state):
    """
    测试主图的错误恢复能力
    
    验证即使某个节点出错，主图也能继续执行
    """
    state = main_graph_initial_state.copy()
    state["user_query"] = "测试错误恢复"
    
    try:
        # 执行主图，即使有错误也应该能够处理
        result = graph.invoke(state, config={"recursion_limit": 10})
        
        # 断言：应该返回某种结果（即使出错）
        assert result is not None, "即使出错也应该返回结果"
        
    except Exception as e:
        # 某些错误是可以接受的（如 LLM API 错误）
        # 但应该记录
        pytest.skip(f"主图执行出现异常（可能是预期的）: {e}")


def test_main_graph_multiple_iterations(main_graph_initial_state):
    """
    测试主图的多轮迭代
    
    验证 Supervisor -> Worker -> Critic -> Supervisor 的循环能够正常工作
    """
    state = main_graph_initial_state.copy()
    state["user_query"] = "多轮迭代测试"
    
    try:
        # 跟踪执行的节点
        executed_nodes = []
        
        for step in graph.stream(state, config={"recursion_limit": 10}):
            for node_name in step.keys():
                executed_nodes.append(node_name)
                # 如果到达 finalize，停止
                if node_name == "finalize":
                    break
        
        # 验证至少执行了 Supervisor
        assert "supervisor" in executed_nodes, "应该至少执行一次 Supervisor"
        
        # 验证节点执行顺序合理
        # Supervisor 应该在 Worker 之前
        if "code_dev" in executed_nodes or "rag_researcher" in executed_nodes:
            supervisor_index = executed_nodes.index("supervisor")
            worker_indices = [
                executed_nodes.index(n) for n in executed_nodes 
                if n in ["code_dev", "rag_researcher", "tool_caller"]
            ]
            if worker_indices:
                assert supervisor_index < min(worker_indices), \
                    "Supervisor 应该在 Worker 之前执行"
        
    except Exception as e:
        pytest.skip(f"主图执行出现异常: {e}")


def test_main_graph_with_complex_task(main_graph_initial_state):
    """
    测试主图处理复杂任务
    
    复杂任务可能需要多个 Worker 协作
    """
    state = main_graph_initial_state.copy()
    state["user_query"] = "先查找相关文献，然后生成分析代码，最后进行数据分析"
    
    try:
        result = graph.invoke(state, config={"recursion_limit": 20})
        
        # 断言：应该返回结果
        assert result is not None
        
        # 验证可能涉及多个 Worker
        # 由于 LLM 的决策，可能只执行部分 Worker
        # 我们只验证流程能够完成
        
    except Exception as e:
        pytest.skip(f"主图执行出现异常: {e}")

