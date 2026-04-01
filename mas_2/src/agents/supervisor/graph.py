"""
Supervisor Agent 子图
负责调度决策，决定下一个执行的 worker
"""
from langchain_core.messages import SystemMessage, HumanMessage
from langgraph.graph import StateGraph, START, END
from pydantic import BaseModel, Field
from typing import Literal, List, Optional
from .state import SupervisorAgentState
from src.core.llm import get_llm
from src.core.state import PlanStep
from src.utils.workflow_skills import format_skills_catalog_for_prompt
from .exploration_plan import generate_exploration_plan

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
        retry_hint += "5. skill_id 可为空：仅当步骤明确对应下方目录中的某一技能时填写\n"

    skills_catalog = format_skills_catalog_for_prompt()

    system_prompt = f"""你是一个专业的项目规划师，负责将复杂的任务分解为可执行的步骤。

请根据用户的任务需求，生成一个详细的执行计划。每个步骤应该：
1. 有明确的输入和输出
2. 指定需要的输入文件路径（如果有）
3. 指定必须生成的输出文件路径（如果有）
4. 包含明确的验收标准，用于判断任务是否成功

【执行前提与背景】
（如果存在前期的初步数据检查运行结果，请参考上下文信息。接下来只需根据具体要求对核心任务进行结构拆解映射。）

【可用 Workflow 技能目录】
这里列出了本地目录下的所有可供参考的分析流程（对应 skill_id）。
{skills_catalog}

【skill_id 规则与分析流程规划】
- 在生成分析流程时，请先阅览上述的技能目录。若发现某个技能的名称或描述非标契合用户的分析任务（例如用户提出要进行基于 Scanpy 的单细胞分析，理应匹配 `scrnaseq-scanpy-core-analysis`），此时你应该为其设置此配套的 `skill_id`。
- **一旦选定可适配的技能并绑定了 skill_id**，那么对应的“正式分析流程部分”的计划粒度请务必设计的粗粒度一点，通常只需要整合为少数 1 步或 2 步核心动作即可即可（例如：“遵循 Scanpy core analysis 的推荐流程执行质控、降维、聚类并绘制图像”）。过于琐碎的拆分不利于对应组件遵从其原本最佳实践的完整工作流。
- 如果不确定，或无任何技能相匹配，不需要勉强，可不填写 skill_id（填 null）。

【数据路径与 input_files】
- 用户若在任务中说明了数据位置（例如「数据路径：…」「输入为 xxx.h5ad」），必须在依赖该数据的步骤的 input_files 中填入具体路径（与原文一致）；不要将路径只写在 description 而遗漏 input_files。
- 多文件输入时，按依赖顺序列出全部需要的文件路径。

【验收标准 acceptance_criteria】
- 若步骤填写了 skill_id，验收标准应对齐该技能文档中期望产生的高层面结果。
- 验收标准须可判定：明确打印了什么文本维度指标，或生成了什么图像才能算成功。

【计划粒度要求】
- 带有 skill_id 的主体分析流程应只包含最小数量的核心调度模块（如一个端到端总控运行步骤）。

【RAG步骤规划原则】
- 可以自主决定是否加入 RAG 检索步骤。
- 当任务存在方法不确定性、参数不明确、需要领域依据时，建议先加 1 个 RAG 步骤。
- 当任务非常基础且需求明确（例如常规 Scanpy 基础流程），可不加 RAG 步骤。

对于代码开发任务，输出文件路径格式应为：{result_path}/step_{{step_id}}_{{filename}}
例如：./result/step_1_umap.png, ./result/step_2_clustering.png

请以 JSON 格式返回计划列表。确保返回的格式完全符合 PlanResponse 模型的要求。"""
    
    exploration_context = ""
    if state.get("data_exploration_results"):
        res_str = "\n\n".join(state["data_exploration_results"])
        exploration_context = f"\n【前期数据探查结果】\n在正式规划分析流程前，我们已经对输入数据进行了摸底，探查结果如下：\n{res_str}\n请根据上述数据特征规划分析流程。\n"

    user_prompt = f"""{exploration_context}
用户任务：{user_query}
结果保存路径：{result_path}
{retry_hint}
请生成详细的执行计划，包含以下信息：
- step_id: 步骤序号（从1开始，必须连续递增）
- name: 步骤名称（简短的动词短语，必填）
- description: 详细的任务描述，包含具体的参数要求（必填）
- input_files: 本步骤需要的输入文件路径列表（用户已提供数据路径时必须填入；无则 []）
- output_files: 本步骤必须生成的输出文件路径列表（格式：{result_path}/step_{{step_id}}_{{filename}}，如果没有则为空列表 []）
- acceptance_criteria: 验收标准，明确说明如何判断任务成功；若绑定 skill_id 则与 SKILL 文档要点一致（必填）
- skill_id: 可选，字符串或 null。若本步骤绑定上表中的某一技能则填写其 id，否则为 null。

请确保计划覆盖用户任务的所有要求，并且所有字段都正确填写。若用户给出了结果/输出目录，计划中 output_files 的路径前缀应与用户期望一致。
"""
    
    try:
        chain = llm_plan.with_structured_output(PlanResponse)
        messages = [
            SystemMessage(content=system_prompt),
            HumanMessage(content=user_prompt)
        ]
        
        print("\n" + "="*60)
        print("🔍 [Supervisor Debug] 发送给 LLM 的生成最终计划提示词：")
        print("[System Prompt]:")
        print(system_prompt)
        print("-" * 60)
        print("[User Prompt]:")
        print(user_prompt)
        print("="*60 + "\n")
        
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
        print("\n" + "="*30 + " 📃 [正式分析计划生成结果] " + "="*30)
        for step in plan:
            print(f"  [步骤 {step.step_id}] {step.name}")
            print(f"  - 描述: {step.description}")
            print(f"  - 验收: {step.acceptance_criteria}")
            print(f"  - 输入: {step.input_files}")
            print(f"  - 输出: {step.output_files}")
            print(f"  - Skill: {step.skill_id}")
            print(f"  " + "-"*40)
        print("="*80 + "\n")
        
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


def _step_needs_execution_beyond_rag(step: Optional[PlanStep]) -> bool:
    """
    当前计划步骤在 RAG 之后是否还需要 code_dev / tool 等交付产物。
    若为 True，则 RAG 经 Critic 通过后不推进 current_step_index，避免单步「绘图」类任务在只做 RAG 后就 FINISH。
    """
    if step is None:
        return False
    outs = list(getattr(step, "output_files", None) or [])
    if outs:
        return True
    blob = " ".join(
        [
            str(getattr(step, "name", "") or ""),
            str(getattr(step, "description", "") or ""),
            str(getattr(step, "acceptance_criteria", "") or ""),
        ]
    ).lower()
    hints = (
        "python",
        "代码",
        "执行",
        "docker",
        "绘图",
        "绘制",
        "可视化",
        "figure",
        "plot",
        ".png",
        ".pdf",
        ".jpg",
        "matplotlib",
        "scanpy",
        "seaborn",
        "qqman",
        "plotly",
        "输出文件",
        "脚本",
        "保存图",
        "manhattan",
    )
    return any(h in blob for h in hints)


def make_decision(state: SupervisorAgentState) -> SupervisorAgentState:
    """
    决策节点
    根据当前状态决定下一个执行的 worker
    如果计划存在，则根据计划派遣任务；否则使用动态决策
    """
    print("--- [Supervisor] 正在做调度决策 ---")
    
    # 获取当前状态
    user_query = state.get("user_query", "")

    # 判断并初始化 result_path 和 task_id
    result_path = state.get("result_path")
    if not result_path:
        # 如果尚未提取路径，则尝试解析并生成结果路径目录
        try:
            from src.agents.code_dev.graph import extract_paths_from_state
            # typed dict 转换为包含 CodeAgentState 特有字段的假 dict，或者直接传因为 python dict 是鸭子类型
            state = extract_paths_from_state(state)
        except Exception as e:
            # 万一模块导入出错，提供后备处理
            import os
            import uuid
            task_id = state.get("task_id", uuid.uuid4().hex[:8])
            state["task_id"] = task_id
            state["result_path"] = os.path.join(".", "result", task_id).replace("\\", "/")
            
    # 从中读出最新的 plan 和由于上方代码可能更新出的 result_path
    result_path = state.get("result_path", "./result")
    plan = state.get("plan", [])
    current_step_index = state.get("current_step_index", 0)
    rag_context = state.get("rag_context", "")
    code_solution = state.get("code_solution", "")
    final_report = state.get("final_report", "")
    is_approved = state.get("is_approved", False)
    last_worker = state.get("last_worker", "")
    pending_contribution = state.get("pending_contribution")
    
    # 如果计划为空，先生成计划
    if not plan:
        data_exploration_done = state.get("data_exploration_done", False)
        try:
            from src.agents.code_dev.graph import parse_paths_from_query
            parsed_paths = parse_paths_from_query(user_query)
            data_path = parsed_paths.get("data_path")
            has_data = bool(data_path)
            
            if not has_data:
                # 兼容旧版的逻辑如果有的回退
                input_files = parsed_paths.get("input_files", [])
                has_data = len(input_files) > 0
                data_path = input_files[0] if has_data else None
                
        except Exception:
            has_data = False
            data_path = None
            
        if has_data and not data_exploration_done:
            state = generate_exploration_plan(state, data_path)
        else:
            state = generate_plan(state)
            
        plan = state.get("plan", [])
        current_step_index = state.get("current_step_index", 0)
    
    # 步骤推进逻辑：如果上一步审核通过，且不是 critic 执行的，则推进步骤索引
    # 注意：仅 RAG 通过不能算作「代码/交付类」步骤完成，否则单步计划会在 RAG 后直接 FINISH。
    if plan and is_approved and last_worker != "critic" and last_worker:
        skip_advance = (
            last_worker == "rag_researcher"
            and 0 <= current_step_index < len(plan)
            and _step_needs_execution_beyond_rag(plan[current_step_index])
        )
        if skip_advance:
            print(
                "  --> RAG 已通过，本步仍需代码/工具交付，暂不推进步骤索引，继续当前步"
            )
        else:
            new_index = current_step_index + 1
            if new_index < len(plan):
                state["current_step_index"] = new_index
                current_step_index = new_index
                print(
                    f"  --> 步骤 {current_step_index} 审核通过，推进到步骤 {current_step_index + 1}"
                )
            else:
                state["current_step_index"] = new_index
                current_step_index = new_index
    
    # 若有计划，则只负责提供当前步骤上下文给大模型，不做关键词硬路由
    if plan and current_step_index < len(plan):
        current_step = plan[current_step_index]
        print(f"  --> 当前执行步骤 {current_step_index + 1}/{len(plan)}: {current_step.name}")

        state["current_step_input"] = current_step.description
        state["current_step_expected_output"] = current_step.acceptance_criteria
        state["current_step_file_paths"] = {
            "input_files": current_step.input_files,
            "output_files": current_step.output_files
        }
        sid = getattr(current_step, "skill_id", None)
        state["current_step_skill_id"] = sid if sid else None
    elif plan and current_step_index >= len(plan):
        # 关键收敛条件：计划已完成 + 无待审核内容 + 最近一步已审核通过 -> 直接结束
        if is_approved and not pending_contribution:
            if not state.get("data_exploration_done", True):
                print("  --> 探查计划已全部完成且审核通过，进入正式分析计划阶段")
                state["data_exploration_done"] = True
                state["plan"] = []
                state["current_step_index"] = 0
                state["is_approved"] = False
                state["last_worker"] = "supervisor"
                return make_decision(state)  # 重新走决策逻辑，这次会生成正式计划
                
            print("  --> 所有计划步骤已完成且审核通过，直接 FINISH")
            state["next_worker"] = "FINISH"
            return state

        print("  --> 所有计划步骤已完成，交由大模型决定是否 FINISH")
        state["current_step_input"] = ""
        state["current_step_expected_output"] = ""
        state["current_step_file_paths"] = {"input_files": [], "output_files": []}
        state["current_step_skill_id"] = None

    # 统一使用动态决策（包括是否 FINISH）
    return _make_dynamic_decision(state, user_query, rag_context, code_solution, 
                                  final_report, is_approved, last_worker, pending_contribution)


def _make_dynamic_decision(state: SupervisorAgentState, user_query: str, rag_context: str,
                           code_solution: str, final_report: str, is_approved: bool,
                           last_worker: str, pending_contribution) -> SupervisorAgentState:
    """
    动态决策函数（原有逻辑）
    当没有计划或需要动态调整时使用
    """
    
    plan = state.get("plan", [])
    current_step_index = state.get("current_step_index", 0)
    current_step_input = state.get("current_step_input", "")
    current_step_expected_output = state.get("current_step_expected_output", "")

    plan_progress = f"{current_step_index}/{len(plan)}" if plan else "无计划"
    current_step_name = ""
    if plan and 0 <= current_step_index < len(plan):
        current_step_name = plan[current_step_index].name

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
- 计划进度: {plan_progress}
- 当前步骤名称: {current_step_name if current_step_name else "无"}
- 当前步骤输入: {current_step_input if current_step_input else "无"}
- 当前步骤验收标准: {current_step_expected_output if current_step_expected_output else "无"}
- RAG 上下文: {"已获取" if rag_context else "未获取"}
- 代码解决方案: {"已生成" if code_solution else "未生成"}
- 最终报告: {"已生成" if final_report else "未生成"}
- 上一个 Worker: {last_worker}
- 待审核内容: {"有" if pending_contribution else "无"}
- 审核状态: {"已通过" if is_approved else "未通过/待审核"}

可用的 Worker：
1. rag_researcher: 检索相关文献和文档
2. code_dev: 生成和执行代码（包含数据读取、预处理、聚类、UMAP绘图、结果文件生成）
3. tool_caller: 仅用于调用现有内置工具（run_celltype_annotation / gene_set_enrichment / query_mygene）
4. critic: 审核工作成果
5. FINISH: 任务完成

决策原则：
- 可自主决定是否先调用 rag_researcher：
    - 当任务存在方法不确定性、参数不明确、需要领域依据时，优先考虑先检索
    - 当任务非常基础且目标明确时，可直接进入 code_dev
- 如果有待审核内容，必须调用 critic
- 若任务需要编写或执行 Python 代码（例如读取 h5ad、运行 scanpy、绘制 UMAP、保存图片/文件），必须调用 code_dev
- 只有当任务本身是“调用内置工具接口”时才调用 tool_caller，不要把常规分析/绘图任务派给 tool_caller
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

