"""
主图编排模块
将所有 SubGraph Agent 连接成完整的工作流
"""
import os
from pathlib import Path

from langgraph.graph import StateGraph, START, END
from langchain_core.messages import SystemMessage, HumanMessage
from src.core.state import GlobalState
from src.core.llm import get_llm, format_llm_runtime_banner

# 导入所有编译好的子图
from src.agents.supervisor.graph import supervisor_agent_graph
from src.agents.critic.graph import critic_agent_graph
from src.agents.code_dev.graph import code_agent_graph
from src.agents.rag_researcher.graph import rag_agent_graph
from src.agents.tool_caller.graph import tool_caller_agent_graph
from src.utils.docker_log_summary import summarize_docker_stdout


# ==================== Wrapper 节点 ====================

def runtime_banner_step(state: GlobalState) -> GlobalState:
    """在每轮运行最开始打印实际注册的 LLM 运行配置摘要。"""
    print(format_llm_runtime_banner())
    return state

def wrap_rag_researcher(state: GlobalState) -> GlobalState:
    """包装 RAG Researcher，更新 last_worker"""
    result = rag_agent_graph.invoke(state)
    return {
        **result,
        "last_worker": "rag_researcher"
    }


def wrap_tool_caller(state: GlobalState) -> GlobalState:
    """包装 Tool Caller，更新 last_worker"""
    result = tool_caller_agent_graph.invoke(state)
    return {
        **result,
        "last_worker": "data_analyst"  # 映射为 data_analyst
    }


# ==================== 路由逻辑 ====================

def supervisor_router(state: GlobalState) -> str:
    """
    Supervisor 路由：根据 next_worker 字段决定去哪个 Agent
    
    Returns:
        目标节点的名称
    """
    next_worker = state.get("next_worker", "FINISH")
    
    if next_worker == "FINISH":
        return "finalize"
    elif next_worker == "rag_researcher":
        return "rag_researcher"
    elif next_worker == "code_dev":
        return "code_dev"
    elif next_worker == "tool_caller":
        return "tool_caller"
    elif next_worker == "data_analyst":
        # data_analyst 使用 tool_caller 实现
        return "tool_caller"
    elif next_worker == "critic":
        return "critic"
    else:
        # 默认结束
        return "finalize"


def critic_router(state: GlobalState) -> str:
    """
    Critic 路由：根据 is_approved 字段决定下一步
    
    - 如果通过审核：返回 supervisor 进行下一轮决策
    - 如果驳回：返回 last_worker 重新执行
    
    Returns:
        目标节点的名称
    """
    is_approved = state.get("is_approved", False)
    last_worker = state.get("last_worker", "")
    
    if is_approved:
        # 审核通过 -> 回 Supervisor 进行下一轮决策
        print("  --> 审核通过，返回 Supervisor")
        return "supervisor"
    else:
        if state.get("step_blocked"):
            print("  --> 当前步骤达到最大驳回轮次，返回 Supervisor 做终止决策")
            return "supervisor"
        # 审核驳回 -> 返回上一个 Worker 重做
        print(f"  --> 审核驳回，返回 {last_worker} 重做")
        
        # 映射 worker 名称到节点名称
        worker_to_node = {
            "rag_researcher": "rag_researcher",
            "code_dev": "code_dev",
            "data_analyst": "tool_caller",
            "tool_caller": "tool_caller"
        }
        
        return worker_to_node.get(last_worker, "supervisor")


def finalize_step(state: GlobalState) -> GlobalState:
    """
    最终化节点：整合所有结果，生成最终答案
    """
    print("\n=== [Finalize] 生成最终答案 ===")

    def _truncate(text: str, max_len: int = 8000) -> str:
        text = str(text or "").strip()
        if len(text) <= max_len:
            return text
        return text[:max_len] + "\n...(内容过长已截断)"

    def _format_plan(plan) -> str:
        if not plan:
            return "无"

        lines = []
        for idx, step in enumerate(plan, start=1):
            if hasattr(step, "step_id"):
                step_id = getattr(step, "step_id", idx)
                name = getattr(step, "name", "")
                description = getattr(step, "description", "")
                acceptance = getattr(step, "acceptance_criteria", "")
                input_files = getattr(step, "input_files", []) or []
                output_files = getattr(step, "output_files", []) or []
            elif isinstance(step, dict):
                step_id = step.get("step_id", idx)
                name = step.get("name", "")
                description = step.get("description", "")
                acceptance = step.get("acceptance_criteria", "")
                input_files = step.get("input_files", []) or []
                output_files = step.get("output_files", []) or []
            else:
                lines.append(f"{idx}. {step}")
                continue

            lines.append(f"{step_id}. {name}")
            if description:
                lines.append(f"   - 说明: {description}")
            if input_files:
                lines.append(f"   - 输入: {', '.join(map(str, input_files))}")
            if output_files:
                lines.append(f"   - 输出: {', '.join(map(str, output_files))}")
            if acceptance:
                lines.append(f"   - 验收: {acceptance}")
            if isinstance(step, dict):
                skill_id = step.get("skill_id")
            else:
                skill_id = getattr(step, "skill_id", None)
            if skill_id:
                lines.append(f"   - skill_id: {skill_id}")

        return "\n".join(lines)

    def _format_code_solution(code_solution) -> str:
        if not code_solution:
            return ""
        if not isinstance(code_solution, dict):
            return str(code_solution)

        def _format_timing(timing) -> str:
            if not isinstance(timing, dict) or not timing:
                return ""

            def _fmt(key: str) -> str:
                value = timing.get(key, 0.0)
                try:
                    return f"{float(value):.2f}s"
                except (TypeError, ValueError):
                    return "n/a"

            install_value = _fmt("pip_install_elapsed_seconds")
            if timing.get("pip_install_skipped", False):
                install_value += " (skipped)"

            parts = [
                f"container_setup={_fmt('container_setup_elapsed_seconds')}",
                f"push_to_container={_fmt('push_to_container_elapsed_seconds')}",
                f"pip_install={install_value}",
                f"python_exec={_fmt('python_exec_elapsed_seconds')}",
                f"collect_outputs={_fmt('collect_outputs_elapsed_seconds')}",
                f"total={_fmt('total_elapsed_seconds')}",
            ]
            if timing.get("install_exit_code") is not None:
                parts.append(f"install_exit_code={timing.get('install_exit_code')}")
            if timing.get("python_exit_code") is not None:
                parts.append(f"python_exit_code={timing.get('python_exit_code')}")
            return ", ".join(parts)

        def _format_phase_timing(phase_timing) -> str:
            if not isinstance(phase_timing, dict) or not phase_timing:
                return ""

            def _fmt(key: str, mapping: dict) -> str:
                value = mapping.get(key, 0.0)
                try:
                    return f"{float(value):.2f}s"
                except (TypeError, ValueError):
                    return "n/a"

            parts = [
                f"attempt={phase_timing.get('attempt', '?')}",
                f"is_retry={bool(phase_timing.get('is_retry', False))}",
            ]

            gen = phase_timing.get("generate_code")
            if isinstance(gen, dict):
                parts.append(
                    "generate_code="
                    + ", ".join(
                        [
                            f"prompt_prep={_fmt('prompt_prep_elapsed_seconds', gen)}",
                            f"tool_loop={_fmt('tool_loop_elapsed_seconds', gen)}",
                            f"tool_llm={_fmt('tool_llm_elapsed_seconds', gen)}",
                            f"tool_io={_fmt('tool_io_elapsed_seconds', gen)}",
                            f"llm_generate={_fmt('llm_generate_elapsed_seconds', gen)}",
                            f"parse={_fmt('parse_elapsed_seconds', gen)}",
                            f"total={_fmt('total_elapsed_seconds', gen)}",
                            f"tool_iterations={gen.get('tool_iterations', 0)}",
                            f"tool_calls={gen.get('tool_calls', 0)}",
                        ]
                    )
                )

            reflection = phase_timing.get("self_reflection")
            if isinstance(reflection, dict):
                parts.append(
                    "self_reflection="
                    + ", ".join(
                        [
                            f"total={_fmt('total_elapsed_seconds', reflection)}",
                            f"warnings={reflection.get('warnings_count', 0)}",
                        ]
                    )
                )

            execute = phase_timing.get("execute_code")
            if isinstance(execute, dict):
                parts.append(
                    "execute_code="
                    + ", ".join(
                        [
                            f"prep={_fmt('prep_elapsed_seconds', execute)}",
                            f"executor_call={_fmt('executor_call_elapsed_seconds', execute)}",
                            f"result_parse={_fmt('result_parse_elapsed_seconds', execute)}",
                            f"total={_fmt('total_elapsed_seconds', execute)}",
                        ]
                    )
                )

            retry = phase_timing.get("prepare_retry")
            if isinstance(retry, dict):
                parts.append(
                    "prepare_retry="
                    + ", ".join(
                        [
                            f"feedback_extract={_fmt('feedback_extract_elapsed_seconds', retry)}",
                            f"total={_fmt('total_elapsed_seconds', retry)}",
                        ]
                    )
                )
            return " | ".join(parts)

        code_text = str(code_solution.get("code", "")).strip()
        exec_result = str(code_solution.get("result", "")).strip()
        raw_output = str(code_solution.get("output", "")).strip()
        display_output = str(code_solution.get("output_display", "")).strip()
        tail_output = str(code_solution.get("output_tail", "")).strip()
        log_disk = str(code_solution.get("output_log_path", "")).strip()
        timing_text = _format_timing(code_solution.get("timing"))
        phase_timing_text = _format_phase_timing(code_solution.get("phase_timing"))
        if not display_output and raw_output:
            display_output = summarize_docker_stdout(raw_output)
        if not display_output and tail_output:
            display_output = summarize_docker_stdout(tail_output)
        use_full_log = os.environ.get("MAS_FULL_EXEC_LOG_IN_FINALIZE", "").strip().lower() in (
            "1",
            "true",
            "yes",
        )
        if use_full_log:
            log_text = raw_output
            if not log_text and log_disk:
                try:
                    log_text = Path(log_disk).read_text(encoding="utf-8")
                except OSError:
                    log_text = display_output or tail_output
            log_label = "完整执行日志"
        else:
            log_text = display_output or tail_output or raw_output
            log_label = "压缩执行日志"
        output_files = code_solution.get("output_files", []) or []

        lines = []
        if code_text:
            lines.append("代码内容:")
            lines.append(code_text)
        if exec_result:
            lines.append("执行结果（===RESULT=== 提取）:")
            lines.append(exec_result)
        if timing_text:
            lines.append("Docker 阶段耗时:")
            lines.append(timing_text)
        if phase_timing_text:
            lines.append("Code Dev 阶段耗时:")
            lines.append(phase_timing_text)
        if log_text:
            lines.append(f"{log_label}:")
            lines.append(log_text)
        if output_files:
            lines.append("产出文件:")
            lines.extend([f"- {f}" for f in output_files])

        return "\n".join(lines)

    # 整合所有结果
    user_query = state.get("user_query", "")
    plan = state.get("plan", [])
    current_step_index = state.get("current_step_index", 0)
    rag_context = state.get("rag_context", "")
    code_solution = state.get("code_solution", "")
    final_report = state.get("final_report", "")
    critique_feedback = state.get("critique_feedback", "")
    is_approved = state.get("is_approved", False)
    pending_contribution = state.get("pending_contribution", "")

    plan_text = _format_plan(plan)
    code_solution_text = _format_code_solution(code_solution)

    # 先准备兜底答案，确保任何情况下都有输出
    fallback_parts = []
    if rag_context:
        fallback_parts.append(f"【相关文献】\n{rag_context}")
    if code_solution_text:
        fallback_parts.append(f"【代码方案】\n{code_solution_text}")
    if final_report:
        fallback_parts.append(f"【分析报告】\n{final_report}")
    if critique_feedback:
        fallback_parts.append(f"【审核反馈】\n{critique_feedback}")

    fallback_answer = "\n\n".join(fallback_parts) if fallback_parts else user_query or "任务已完成"

    # 用 LLM 进行最终整合
    system_prompt = (
        "你是多智能体系统的最终答复整合器。"
        "请基于用户原始问题、执行计划和各环节中间结果，"
        "产出一份结构清晰、结论明确、可直接给用户查看的最终回答。"
        "回答要求："
        "1) 先给结论，再给关键依据；"
        "2) 对代码执行结果与产出文件要明确点出；"
        "3) 若信息不足，明确说明缺口，不要编造；"
        "4) 使用中文，简洁专业。"
    )

    user_prompt = f"""
请整合以下上下文，生成最终回答。

[用户原始问题]
{_truncate(user_query, 2000)}

[执行计划]
{_truncate(plan_text, 5000)}

[当前步骤索引]
{current_step_index}

[RAG 检索上下文]
{_truncate(rag_context, 8000)}

[代码执行产出]
{_truncate(code_solution_text, 256_000)}

[分析报告]
{_truncate(final_report, 8000)}

[Critic 审核信息]
是否通过: {is_approved}
反馈: {_truncate(str(critique_feedback), 2000)}

[待审核贡献（如有）]
{_truncate(str(pending_contribution), 4000)}
""".strip()

    final_answer = ""
    try:
        llm = get_llm(temperature=0.2, streaming=True)
        response = llm.invoke([
            SystemMessage(content=system_prompt),
            HumanMessage(content=user_prompt)
        ])
        llm_content = getattr(response, "content", "")

        if isinstance(llm_content, list):
            llm_content = "\n".join(str(item) for item in llm_content)

        final_answer = str(llm_content).strip()
        if not final_answer:
            raise ValueError("LLM 返回空内容")
    except Exception as e:
        print(f"[Finalize] LLM 汇总失败，使用兜底答案: {e}")
        final_answer = fallback_answer

    # 任务结束后统一清理当前实验专属 task container
    if state.get("docker_container_id"):
        try:
            import docker
            client = docker.from_env()
            container = client.containers.get(state["docker_container_id"])
            container.remove(force=True)
            print(f"\n=== [Finalize] 已清理实验 Docker task container: {state['docker_container_id'][:12]} ===")
        except Exception as e:
            print(f"\n=== [Finalize] 清理实验 Docker task container 失败: {e} ===")
        state["docker_container_id"] = None

    return {
        **state,
        "final_answer": final_answer,
        "next_worker": "FINISH"
    }


# ==================== 构建主图 ====================

workflow = StateGraph(GlobalState)

# 1. 添加所有节点
workflow.add_node("runtime_banner", runtime_banner_step)
# Supervisor 和 Critic 直接使用子图
workflow.add_node("supervisor", supervisor_agent_graph)
workflow.add_node("critic", critic_agent_graph)

# Worker 节点使用 wrapper 函数，确保更新 last_worker
workflow.add_node("rag_researcher", wrap_rag_researcher)
# 子图直接挂入主图，便于 astream_events 暴露 generate_code / execute_code 等内部节点
workflow.add_node("code_dev", code_agent_graph)
workflow.add_node("tool_caller", wrap_tool_caller)  # data_analyst 使用 tool_caller

# Finalize 节点
workflow.add_node("finalize", finalize_step)

# 2. 定义流程

# START -> runtime_banner -> supervisor
workflow.add_edge(START, "runtime_banner")
workflow.add_edge("runtime_banner", "supervisor")

# supervisor -> (路由) -> [rag_researcher, code_dev, tool_caller, critic, finalize]
workflow.add_conditional_edges(
    "supervisor",
    supervisor_router,
    {
        "rag_researcher": "rag_researcher",
        "code_dev": "code_dev",
        "tool_caller": "tool_caller",
        "critic": "critic",
        "finalize": "finalize"
    }
)

# 所有 Worker -> critic (所有产出必须经过审核)
workflow.add_edge("rag_researcher", "critic")
workflow.add_edge("code_dev", "critic")
workflow.add_edge("tool_caller", "critic")

# critic -> (路由) -> [supervisor OR Worker]
workflow.add_conditional_edges(
    "critic",
    critic_router,
    {
        "supervisor": "supervisor",
        "rag_researcher": "rag_researcher",
        "code_dev": "code_dev",
        "tool_caller": "tool_caller"
    }
)

# finalize -> END
workflow.add_edge("finalize", END)

# 3. 编译主图
graph = workflow.compile()
