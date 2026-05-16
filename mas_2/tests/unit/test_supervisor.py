"""
Supervisor Agent 单元测试
使用 Mock 模拟 LLM 调用
"""
import pytest
from unittest.mock import patch, MagicMock
from src.agents.supervisor.graph import RouteDecision, generate_plan, make_decision, supervisor_agent_graph
from src.core.state import PlanStep


def _build_single_step_plan() -> list[PlanStep]:
    """构造最小可执行计划，避免单测触发真实计划生成。"""
    return [
        PlanStep(
            step_id=1,
            name="生成代码",
            description="生成并执行一段测试代码",
            input_files=[],
            output_files=["./result/step_1_test.txt"],
            acceptance_criteria="输出测试结果",
        )
    ]


def test_supervisor_decision_with_mock(supervisor_state):
    """
    测试 Supervisor Agent 的决策逻辑
    
    使用 @patch 模拟 LLM 的返回，验证路由决策是否正确
    """
    # 准备测试状态
    state = supervisor_state.copy()
    state["user_query"] = "需要生成代码"
    state["plan"] = _build_single_step_plan()
    state["current_step_index"] = 0
    
    # 创建模拟的 RouteDecision 对象
    mock_decision = RouteDecision(
        next_worker="code_dev",
        reasoning="用户需要生成代码，应该调用 code_dev worker"
    )
    
    # Mock LLM 的 with_structured_output 方法
    with patch('src.agents.supervisor.graph.llm') as mock_llm:
        # 创建模拟的 chain 对象
        mock_chain = MagicMock()
        mock_chain.invoke.return_value = mock_decision
        
        # 设置 with_structured_output 返回 mock_chain
        mock_llm.with_structured_output.return_value = mock_chain
        
        # 执行 Supervisor Agent 子图
        result = supervisor_agent_graph.invoke(state)
        
        # 断言：next_worker 字段确实变成了 "code_dev"
        assert result.get("next_worker") == "code_dev", \
            f"期望 next_worker='code_dev'，实际为 '{result.get('next_worker')}'"

        # 断言：当前步骤上下文已写入状态
        assert result.get("current_step_input") == "生成并执行一段测试代码"
        assert result.get("current_step_expected_output") == "输出测试结果"
        assert result.get("current_step_file_paths") == {
            "input_files": [],
            "output_files": ["./result/step_1_test.txt"],
        }
        assert result.get("current_step_skill_id") is None
        
        # 验证 LLM 被调用
        mock_llm.with_structured_output.assert_called_once()
        mock_chain.invoke.assert_called_once()


def test_supervisor_decision_finish(supervisor_state):
    """
    测试 Supervisor 决定结束任务的情况
    """
    state = supervisor_state.copy()
    state["user_query"] = "任务已完成"
    state["rag_context"] = "已有上下文"
    state["code_solution"] = "已有代码"
    state["plan"] = _build_single_step_plan()
    state["current_step_index"] = 0
    
    # 创建模拟的 RouteDecision 对象（决定结束）
    mock_decision = RouteDecision(
        next_worker="FINISH",
        reasoning="所有任务已完成，可以结束"
    )
    
    with patch('src.agents.supervisor.graph.llm') as mock_llm:
        mock_chain = MagicMock()
        mock_chain.invoke.return_value = mock_decision
        mock_llm.with_structured_output.return_value = mock_chain
        
        result = supervisor_agent_graph.invoke(state)
        
        # 断言：next_worker 应该是 "FINISH"
        assert result.get("next_worker") == "FINISH", \
            f"期望 next_worker='FINISH'，实际为 '{result.get('next_worker')}'"


def test_supervisor_decision_rag_researcher(supervisor_state):
    """
    测试 Supervisor 决定调用 RAG Researcher 的情况
    """
    state = supervisor_state.copy()
    state["user_query"] = "查询相关文献"
    state["rag_context"] = ""  # 没有 RAG 上下文
    state["plan"] = _build_single_step_plan()
    state["current_step_index"] = 0
    
    mock_decision = RouteDecision(
        next_worker="rag_researcher",
        reasoning="需要检索相关文献"
    )
    
    with patch('src.agents.supervisor.graph.llm') as mock_llm:
        mock_chain = MagicMock()
        mock_chain.invoke.return_value = mock_decision
        mock_llm.with_structured_output.return_value = mock_chain
        
        result = supervisor_agent_graph.invoke(state)
        
        assert result.get("next_worker") == "rag_researcher"


def test_supervisor_decision_error_handling(supervisor_state):
    """
    测试 Supervisor 在 LLM 调用失败时的错误处理
    """
    state = supervisor_state.copy()
    state["user_query"] = "测试错误处理"
    state["pending_contribution"] = "待审核内容"
    state["plan"] = _build_single_step_plan()
    state["current_step_index"] = 0
    
    # Mock LLM 抛出异常
    with patch('src.agents.supervisor.graph.llm') as mock_llm:
        mock_chain = MagicMock()
        mock_chain.invoke.side_effect = Exception("LLM API 调用失败")
        mock_llm.with_structured_output.return_value = mock_chain
        
        result = supervisor_agent_graph.invoke(state)
        
        # 断言：当有 pending_contribution 时，应该默认选择 critic
        assert result.get("next_worker") == "critic", \
            "当 LLM 调用失败且有待审核内容时，应默认选择 critic"


def test_supervisor_decision_error_handling_no_pending(supervisor_state):
    """
    测试 Supervisor 在 LLM 调用失败且没有待审核内容时的错误处理
    """
    state = supervisor_state.copy()
    state["user_query"] = "测试错误处理"
    state["pending_contribution"] = None
    state["plan"] = _build_single_step_plan()
    state["current_step_index"] = 0
    
    with patch('src.agents.supervisor.graph.llm') as mock_llm:
        mock_chain = MagicMock()
        mock_chain.invoke.side_effect = Exception("LLM API 调用失败")
        mock_llm.with_structured_output.return_value = mock_chain
        
        result = supervisor_agent_graph.invoke(state)
        
        # 断言：当没有 pending_contribution 时，应该默认选择 rag_researcher
        assert result.get("next_worker") == "rag_researcher", \
            "当 LLM 调用失败且没有待审核内容时，应默认选择 rag_researcher"


def test_supervisor_direct_finish_when_plan_done(supervisor_state):
    """测试计划完成且审核通过时，不依赖 LLM 直接 FINISH。"""
    state = supervisor_state.copy()
    state["plan"] = _build_single_step_plan()
    state["current_step_index"] = 1  # 已超出最后一步
    state["is_approved"] = True
    state["pending_contribution"] = None

    with patch('src.agents.supervisor.graph.llm') as mock_llm:
        result = supervisor_agent_graph.invoke(state)

        assert result.get("next_worker") == "FINISH"
        mock_llm.with_structured_output.assert_not_called()


def test_supervisor_rag_pass_does_not_finish_single_code_step(supervisor_state):
    """RAG 通过后，单步且含产出路径的计划不得直接 FINISH，应继续当前步并走动态决策。"""
    state = supervisor_state.copy()
    state["plan"] = _build_single_step_plan()
    state["current_step_index"] = 0
    state["is_approved"] = True
    state["last_worker"] = "rag_researcher"
    state["rag_context"] = "已检索"
    state["pending_contribution"] = None

    mock_decision = RouteDecision(
        next_worker="code_dev",
        reasoning="已有 RAG，应执行代码生成产出文件",
    )
    with patch("src.agents.supervisor.graph.llm") as mock_llm:
        mock_chain = MagicMock()
        mock_chain.invoke.return_value = mock_decision
        mock_llm.with_structured_output.return_value = mock_chain

        result = supervisor_agent_graph.invoke(state)

    assert result.get("current_step_index") == 0
    assert result.get("next_worker") == "code_dev"


def test_supervisor_rag_pass_advances_pure_retrieval_step(supervisor_state):
    """纯检索步（无产出文件、无代码类关键词）RAG 通过后应推进并直接 FINISH。"""
    state = supervisor_state.copy()
    state["plan"] = [
        PlanStep(
            step_id=1,
            name="文献检索",
            description="从向量库检索相关文献并列出要点",
            input_files=[],
            output_files=[],
            acceptance_criteria="列出至少三条相关摘要",
        )
    ]
    state["current_step_index"] = 0
    state["is_approved"] = True
    state["last_worker"] = "rag_researcher"
    state["rag_context"] = "已检索"
    state["pending_contribution"] = None

    with patch("src.agents.supervisor.graph.llm") as mock_llm:
        result = supervisor_agent_graph.invoke(state)
        mock_llm.with_structured_output.assert_not_called()

    assert result.get("next_worker") == "FINISH"


def test_supervisor_preselects_skill_before_exploration(supervisor_state, monkeypatch):
    """有数据且尚未 exploration 时，应先预选 skill 再生成探查计划。"""
    state = supervisor_state.copy()
    state["user_query"] = "我想用data目录下的bmmc_b_cell.h5ad单细胞数据进行分析"
    state["plan"] = []
    state["data_exploration_done"] = False
    state["selected_skill_id"] = None

    call_order: list[str] = []

    def _fake_select_skill(current_state):
        call_order.append("select")
        current_state["selected_skill_id"] = "scrnaseq-scanpy-core-analysis"
        return current_state

    def _fake_generate_exploration_plan(current_state, data_path):
        call_order.append("explore")
        assert data_path == "bmmc_b_cell.h5ad"
        assert current_state.get("selected_skill_id") == "scrnaseq-scanpy-core-analysis"
        current_state["plan"] = _build_single_step_plan()
        current_state["current_step_index"] = 0
        return current_state

    mock_decision = RouteDecision(
        next_worker="code_dev",
        reasoning="已有 exploration 计划，继续执行代码",
    )

    monkeypatch.setattr("src.agents.supervisor.graph.select_skill_for_plan", _fake_select_skill)
    monkeypatch.setattr("src.agents.supervisor.graph.generate_exploration_plan", _fake_generate_exploration_plan)
    monkeypatch.setattr(
        "src.agents.code_dev.graph.parse_paths_from_query",
        lambda query: {"data_path": "bmmc_b_cell.h5ad", "input_files": ["bmmc_b_cell.h5ad"]},
    )

    with patch("src.agents.supervisor.graph.llm") as mock_llm:
        mock_chain = MagicMock()
        mock_chain.invoke.return_value = mock_decision
        mock_llm.with_structured_output.return_value = mock_chain

        result = make_decision(state)

    assert call_order == ["select", "explore"]
    assert result.get("selected_skill_id") == "scrnaseq-scanpy-core-analysis"
    assert result.get("next_worker") == "code_dev"


def test_generate_plan_profile_prewarm_does_not_persist_task_container(supervisor_state):
    state = supervisor_state.copy()
    state["selected_skill_id"] = "scrnaseq-scanpy-core-analysis"
    state["docker_container_id"] = None

    plan_response = MagicMock()
    plan_response.plan = _build_single_step_plan()
    plan_response.global_requirements = None

    with patch("src.agents.supervisor.graph.llm_plan") as mock_llm, patch(
        "src.agents.supervisor.graph.CodeExecutor"
    ) as mock_executor_cls:
        mock_chain = MagicMock()
        mock_chain.invoke.return_value = plan_response
        mock_llm.with_structured_output.return_value = mock_chain

        mock_executor = MagicMock()
        mock_executor.ensure_env_cache_ready.return_value = {
            "success": True,
            "timing": {
                "cache_prepare_elapsed_seconds": 0.01,
                "cache_hit": True,
                "cache_install_performed": False,
            },
        }
        mock_executor_cls.return_value = mock_executor

        result = generate_plan(state)

    assert result["docker_container_id"] is None
    assert result["docker_env_cache_key"]
    assert result["docker_env_cache_volume"]
    mock_executor.ensure_env_cache_ready.assert_called_once()
