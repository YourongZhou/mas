"""
Supervisor Agent 子图
负责调度决策，决定下一个执行的 worker
"""
from langchain_core.messages import SystemMessage, HumanMessage
from langgraph.graph import StateGraph, START, END
from pydantic import BaseModel, Field
from typing import Literal, List
from .state import SupervisorAgentState
from src.core.llm import get_llm
from src.core.state import PlanStep

# 初始化 LLM
llm = get_llm(temperature=0.5)
llm_plan = get_llm(temperature=0.3)  # 用于计划生成，温度稍低以保证稳定性


class RouteDecision(BaseModel):
    """路由决策模型"""
    next_worker: Literal["rag_researcher", "code_dev", "tool_caller", "critic", "FINISH"] = Field(
        ...,
        description="下一个要执行的 worker"
    )
    reasoning: str = Field(..., description="决策理由")


class PlanResponse(BaseModel):
    """计划生成响应模型"""
    plan: List[PlanStep] = Field(..., description="完整的执行计划列表")


def generate_plan(state: SupervisorAgentState, retry_count: int = 0, max_retries: int = 3) -> SupervisorAgentState:
    """
    生成计划节点
    根据用户查询生成完整的执行计划列表
    如果失败会自动重试，最多重试 max_retries 次
    """
    if retry_count == 0:
        print("--- [Supervisor] 正在生成执行计划 ---")
    else:
        print(f"--- [Supervisor] 正在重新生成执行计划 (重试 {retry_count}/{max_retries}) ---")
    
    user_query = state.get("user_query", "")
    result_path = state.get("result_path", "./result")
    
    # 如果是重试，在prompt中添加提示
    retry_hint = ""
    if retry_count > 0:
        retry_hint = f"\n\n【重要】这是第 {retry_count + 1} 次尝试生成计划。请确保：\n"
        retry_hint += "1. 返回的JSON格式完全正确\n"
        retry_hint += "2. 每个步骤的字段都完整填写\n"
        retry_hint += "3. step_id 从1开始，连续递增\n"
        retry_hint += "4. 所有必填字段（step_id, name, description, acceptance_criteria）都有值\n"
    
    system_prompt = """你是一个专业的项目规划师，负责将复杂的任务分解为可执行的步骤。

请根据用户的任务需求，生成一个详细的执行计划。每个步骤应该：
1. 有明确的输入和输出
2. 指定需要的输入文件路径（如果有）
3. 指定必须生成的输出文件路径（如果有）
4. 包含明确的验收标准，用于判断任务是否成功

对于代码开发任务，输出文件路径格式应为：{result_path}/step_{{step_id}}_{{filename}}
例如：./result/step_1_umap.png, ./result/step_2_clustering.png

请以 JSON 格式返回计划列表。确保返回的格式完全符合 PlanResponse 模型的要求。"""
    
    user_prompt = f"""
用户任务：{user_query}
结果保存路径：{result_path}
{retry_hint}
请生成详细的执行计划，包含以下信息：
- step_id: 步骤序号（从1开始，必须连续递增）
- name: 步骤名称（简短的动词短语，必填）
- description: 详细的任务描述，包含具体的参数要求（必填）
- input_files: 本步骤需要的输入文件路径列表（如果没有则为空列表 []）
- output_files: 本步骤必须生成的输出文件路径列表（格式：{result_path}/step_{{step_id}}_{{filename}}，如果没有则为空列表 []）
- acceptance_criteria: 验收标准，明确说明如何判断任务成功（必填）

请确保计划覆盖用户任务的所有要求，并且所有字段都正确填写。
"""
    
    try:
        chain = llm_plan.with_structured_output(PlanResponse)
        messages = [
            SystemMessage(content=system_prompt),
            HumanMessage(content=user_prompt)
        ]
        response = chain.invoke(messages)
        
        plan = response.plan
        
        # 验证计划的有效性
        if not plan or len(plan) == 0:
            raise ValueError("生成的计划为空")
        
        # 验证步骤ID的连续性
        step_ids = [step.step_id for step in plan]
        if step_ids != list(range(1, len(plan) + 1)):
            raise ValueError(f"步骤ID不连续: {step_ids}")
        
        # 验证必填字段
        for i, step in enumerate(plan):
            if not step.name or not step.description or not step.acceptance_criteria:
                raise ValueError(f"步骤 {i+1} 缺少必填字段")
        
        print(f"  --> 计划生成成功，共 {len(plan)} 个步骤")
        for step in plan:
            print(f"     步骤 {step.step_id}: {step.name}")
        
        # 更新状态
        state["plan"] = plan
        state["current_step_index"] = 0
        
    except Exception as e:
        error_msg = str(e)
        print(f"  --> 计划生成失败: {error_msg}")
        
        # 如果还有重试机会，则重试
        if retry_count < max_retries:
            print(f"  --> 将在 {retry_count + 1} 秒后重试...")
            import time
            time.sleep(1)  # 短暂延迟，避免过快重试
            return generate_plan(state, retry_count=retry_count + 1, max_retries=max_retries)
        else:
            print(f"  --> 已达到最大重试次数 ({max_retries})，将使用动态决策模式")
            # 如果计划生成失败且重试次数用完，保持 plan 为空，使用原有的动态决策模式
            state["plan"] = []
            state["current_step_index"] = 0
    
    return state


def make_decision(state: SupervisorAgentState) -> SupervisorAgentState:
    """
    决策节点
    根据当前状态决定下一个执行的 worker
    如果计划存在，则根据计划派遣任务；否则使用动态决策
    """
    print("--- [Supervisor] 正在做调度决策 ---")
    
    # 获取当前状态
    user_query = state.get("user_query", "")
    plan = state.get("plan", [])
    current_step_index = state.get("current_step_index", 0)
    rag_context = state.get("rag_context", "")
    code_solution = state.get("code_solution", "")
    final_report = state.get("final_report", "")
    is_approved = state.get("is_approved", False)
    last_worker = state.get("last_worker", "")
    pending_contribution = state.get("pending_contribution")
    result_path = state.get("result_path", "./result")
    
    # 如果计划为空，先生成计划
    if not plan:
        state = generate_plan(state)
        plan = state.get("plan", [])
        current_step_index = state.get("current_step_index", 0)
    
    # 步骤推进逻辑：如果上一步审核通过，且不是critic执行的，则推进步骤索引
    if plan and is_approved and last_worker != "critic" and last_worker:
        # 上一步审核通过，推进到下一步
        new_index = current_step_index + 1
        if new_index < len(plan):
            state["current_step_index"] = new_index
            current_step_index = new_index
            print(f"  --> 步骤 {current_step_index} 审核通过，推进到步骤 {current_step_index + 1}")
        else:
            # 所有步骤已完成
            state["current_step_index"] = new_index
            current_step_index = new_index
    
    # 如果有待审核内容，必须先审核
    if pending_contribution:
        state["next_worker"] = "critic"
        return state
    
    # 如果计划存在且未完成，根据计划派遣任务
    if plan and current_step_index < len(plan):
        current_step = plan[current_step_index]
        print(f"  --> 当前执行步骤 {current_step_index + 1}/{len(plan)}: {current_step.name}")
        
        # 更新当前步骤的输入、预期输出和文件路径
        state["current_step_input"] = current_step.description
        state["current_step_expected_output"] = current_step.acceptance_criteria
        state["current_step_file_paths"] = {
            "input_files": current_step.input_files,
            "output_files": current_step.output_files
        }
        
        # 根据步骤描述决定派遣哪个 worker
        step_desc_lower = current_step.description.lower()
        step_name_lower = current_step.name.lower()
        
        # 判断任务类型
        if any(keyword in step_desc_lower or keyword in step_name_lower 
               for keyword in ["代码", "code", "编程", "执行", "运行", "分析", "处理数据"]):
            state["next_worker"] = "code_dev"
            print(f"  --> 派遣到 code_dev")
        elif any(keyword in step_desc_lower or keyword in step_name_lower 
                 for keyword in ["检索", "搜索", "文献", "文档", "查询", "search", "retrieve"]):
            state["next_worker"] = "rag_researcher"
            print(f"  --> 派遣到 rag_researcher")
        elif any(keyword in step_desc_lower or keyword in step_name_lower 
                 for keyword in ["工具", "数据库", "调用", "tool", "database", "查询基因", "鉴定"]):
            state["next_worker"] = "tool_caller"
            print(f"  --> 派遣到 tool_caller")
        else:
            # 默认使用动态决策
            state = _make_dynamic_decision(state, user_query, rag_context, code_solution, 
                                          final_report, is_approved, last_worker, pending_contribution)
        
        return state
    
    # 如果所有步骤已完成，检查是否还有待处理的任务
    if plan and current_step_index >= len(plan):
        print("  --> 所有计划步骤已完成")
        state["next_worker"] = "FINISH"
        return state
    
    # 如果没有计划或计划已完成，使用动态决策
    return _make_dynamic_decision(state, user_query, rag_context, code_solution, 
                                  final_report, is_approved, last_worker, pending_contribution)


def _make_dynamic_decision(state: SupervisorAgentState, user_query: str, rag_context: str,
                           code_solution: str, final_report: str, is_approved: bool,
                           last_worker: str, pending_contribution) -> SupervisorAgentState:
    """
    动态决策函数（原有逻辑）
    当没有计划或需要动态调整时使用
    """
    
    # 构建决策 Prompt
    # 注意：DashScope API 要求当使用 response_format: json_object 时，消息中必须包含 "json" 字样
    system_prompt = """你是项目经理，负责协调多个 AI 代理完成用户任务。

请根据当前状态，决定下一个要执行的 worker，并以 JSON 格式返回结果。

返回的 JSON 格式必须严格遵循以下结构：
{
  "next_worker": "rag_researcher" | "code_dev" | "tool_caller" | "critic" | "FINISH",
  "reasoning": "你的决策理由"
}

字段说明：
- next_worker: 下一个要执行的 worker，必须是以下之一：rag_researcher, code_dev, tool_caller, critic, FINISH
- reasoning: 决策理由的详细说明"""
    
    user_prompt = f"""
当前项目状态：
- 用户问题: {user_query}
- RAG 上下文: {"已获取" if rag_context else "未获取"}
- 代码解决方案: {"已生成" if code_solution else "未生成"}
- 最终报告: {"已生成" if final_report else "未生成"}
- 上一个 Worker: {last_worker}
- 待审核内容: {"有" if pending_contribution else "无"}
- 审核状态: {"已通过" if is_approved else "未通过/待审核"}

可用的 Worker：
1. rag_researcher: 检索相关文献和文档
2. code_dev: 生成和执行代码
3. tool_caller: 分析用户问题，调用工具返回结果（比如给定DE 基因鉴定细胞型种类）
4. critic: 审核工作成果
5. FINISH: 任务完成

决策原则：
- 如果有必须要搜索的内容，调用 rag_researcher
- 如果有待审核内容，必须调用 critic
- 如果需要生成代码，调用 code_dev
- 如果需要调用现有工具，调用 tool_caller
- 如果所有工作都已完成，返回 FINISH

请返回 JSON 格式的决策结果，包含 next_worker 和 reasoning 两个字段。
"""
    
    try:
        chain = llm.with_structured_output(RouteDecision)
        messages = [
            SystemMessage(content=system_prompt),
            HumanMessage(content=user_prompt)
        ]
        decision = chain.invoke(messages)
        
        print(f"  --> 决策: {decision.next_worker}")
        print(f"  --> 理由: {decision.reasoning}")
        
        # 更新状态
        if decision.next_worker == "FINISH":
            state["next_worker"] = "FINISH"
        else:
            state["next_worker"] = decision.next_worker
        
    except Exception as e:
        print(f"  --> 决策失败: {e}，默认选择 critic")
        # 如果有待审核内容，默认选择 critic
        if pending_contribution:
            state["next_worker"] = "critic"
        else:
            state["next_worker"] = "rag_researcher"
    
    return state


# 构建子图
workflow = StateGraph(SupervisorAgentState)

# 添加节点
workflow.add_node("make_decision", make_decision)

# 定义边
workflow.add_edge(START, "make_decision")
workflow.add_edge("make_decision", END)

# 编译子图
supervisor_agent_graph = workflow.compile()

