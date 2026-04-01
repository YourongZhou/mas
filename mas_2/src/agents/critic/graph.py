"""
Critic Agent 子图
负责审核 Worker 的工作成果
"""
import base64
import re
from langchain_core.messages import SystemMessage, HumanMessage
from langgraph.graph import StateGraph, START, END
from .state import CriticAgentState
from src.core.llm import get_llm
from src.utils.workflow_skills import format_skill_for_critic

# 初始化 LLM
llm = get_llm(temperature=0.1) # 建议降低温度，让审核更死板、更守规矩
llm_vision = get_llm(model_name="qwen-vl-plus", temperature=0.1)

# --- 全局 System Prompt (增强版) --- 
CRITIC_SYSTEM_PROMPT = """
Role: Senior Bioinformatics Reviewer (Nature/Cell Standard) & Technical Auditor
Profile:
You are a rigorous AI auditor. Your goal is to verify if the CURRENT STEP has been completed according to its SPECIFIC acceptance criteria.

*** CRITICAL RULES FOR MULTI-STEP TASKS ***
1. SCOPE IS LIMITED: You are reviewing ONE STEP of a larger plan (e.g., Step 1 of 7).
   - DO NOT reject the work because it hasn't finished the *entire* project yet.
   - IF Step 1 is "Load Data", and the code loads data: PASS IT. DO NOT ask for UMAP/Clustering/Annotation if that is in Step 4.
   - ONLY judge based on the "Current Step Acceptance Criteria".

2. DOCKER ENVIRONMENT AWARENESS:
   - The code runs in a Docker container.
   - Path Mismatch is EXPECTED: The user says `/home/user/data/file.h5ad`, but the code uses `/app/data/file.h5ad`.
   - THIS IS CORRECT BEHAVIOR (Volume Mounting).
   - NEVER reject code solely because the file path looks different from the user's prompt, AS LONG AS the code executed successfully.

3. EXECUTION LOG IS KING:
   - If the `Execution Result` shows "SUCCESS" or produced the expected output files, you MUST trust the code works, even if the paths look weird.
   - Do not hallucinate errors if the log says it worked.

Output Protocol:
- If the work meets the *current step's* criteria, reply with exactly: "PASS"
- If the work is flawed, reply in the following format:
  [FAIL]
  CRITICAL ISSUE: <Describe the scientific or technical error>
  SUGGESTION: <Actionable advice to fix it>
- Reply in Chinese.
"""


def _normalize_base64_image(image_b64: str, default_mime: str = "image/png") -> str:
    """规范化 base64 图片数据"""
    if not image_b64:
        raise ValueError("empty image base64")

    b64 = image_b64.strip()
    if b64.startswith("data:image/"):
        return b64

    if "base64," in b64:
        b64 = b64.split("base64,", 1)[1].strip()

    b64 = re.sub(r"\s+", "", b64)
    try:
        base64.b64decode(b64, validate=True)
    except Exception as exc:
        raise ValueError("invalid base64 image data") from exc

    return f"data:{default_mime};base64,{b64}"


def check_umap_image(image_base64: str, query: str, expected_output: str = None, 
                     step_context: dict = None) -> str:
    """审核 UMAP 图片质量"""
    image_system_prompt = """
    --- Visualization Review Task ---
    Task: Evaluate the scientific visualization quality and relevance to the User Question.
    
    [Universal Criteria] (Must Have)
    1. Labels: Axis labels must be visible.
    2. Clarity: No severe blurring.
    
    [Step-Specific Context]
    If this is an intermediate step, do not demand final publication polish.
    """

    # 完整 System Prompt
    full_system_prompt = f"{CRITIC_SYSTEM_PROMPT}\n{image_system_prompt}"
    
    step_context_note = ""
    if step_context:
        step_num = step_context.get("step_num", "")
        total_steps = step_context.get("total_steps", "")
        step_context_note = f"\n\n【步骤上下文】当前是步骤 {step_num}/{total_steps}。请只检查本步骤要求生成的图片。\n"
    
    expected_output_note = ""
    if expected_output:
        expected_output_note = f"\n\n【验收标准】\n{expected_output}"
    
    user_prompt = f"User question: {query}{step_context_note}{expected_output_note}"
    
    try:
        data_url = _normalize_base64_image(image_base64)
    except ValueError as exc:
        return f"INVALID_IMAGE: {exc}"

    response = llm_vision.invoke([
        SystemMessage(content=full_system_prompt),
        HumanMessage(
            content=[
                {"type": "text", "text": user_prompt},
                {"type": "image_url", "image_url": {"url": data_url}},
            ]
        ),
    ])
    return response.content


def check_code(content: str, query: str, execution_result: str = None,
               expected_output: str = None, step_context: dict = None,
               skill_note: str = None) -> str:
    """审核代码"""
    # --- 针对代码审核的强化 Prompt ---
    code_system_prompt = """
    你是一个资深代码审查员。请按以下优先级进行检查：

    1. 【最高优先级】执行结果检查 (Execution Check):
       - 查看提供的【代码执行结果/日志】。
       - 如果日志显示 "EXECUTION SUCCESS" 或成功输出了结果标记（如 ===RESULT===），则代码**通过**。
       - 只要运行成功，**绝对不要**因为文件路径与用户输入不同而驳回（这是Docker映射的正常现象）。
       - 只有在日志显示 "Traceback", "Error", "Exception" 时才判定为失败。

    2. 步骤范围检查 (Scope Check):
       - 当前是分步执行模式。
       - **严禁**要求代码包含当前步骤未提及的功能。
       - 例子：如果当前步骤是"读取数据"，代码只要读取并保存了数据就是 PASS。**不要**抱怨"未进行聚类"或"未画图"。

    3. 代码逻辑检查:
       - 只有在没有执行日志的情况下，才深度检查逻辑漏洞。
    """
    
    # 构建步骤上下文信息
    step_context_note = ""
    if step_context:
        step_name = step_context.get("step_name", "")
        step_num = step_context.get("step_num", "")
        total_steps = step_context.get("total_steps", "")
        step_context_note = f"\n\n【重要！当前步骤上下文】\n"
        step_context_note += f"- 进度：步骤 {step_num} / {total_steps}\n"
        step_context_note += f"- 本步骤任务：{step_name}\n"
        step_context_note += f"- **审核红线**：只要完成了'{step_name}'，无论后面还有多少步骤没做，都必须给 PASS。\n"
        
    expected_output_note = ""
    if expected_output:
        expected_output_note = f"\n\n【本步骤验收标准】\n{expected_output}\n(只要满足此标准即可，不要自行加码)"

    execution_note = ""
    if execution_result:
        execution_note = f"\n\n【代码执行结果/日志】\n{execution_result}\n\n"
        # 增加提示引导 Critic 信任日志
        if "EXECUTION SUCCESS" in execution_result or "===RESULT===" in execution_result:
            execution_note += "提示：日志显示执行成功。请忽略路径差异，直接通过。\n"
        else:
            execution_note += "提示：日志显示执行可能存在问题，请仔细检查报错信息。\n"

    skill_extra = skill_note or ""

    user_prompt = f"""
    用户问题: {query}
    待审核代码: 
    ```python
    {content}
    ```
    {execution_note}
    {step_context_note}{expected_output_note}{skill_extra}
    """
    
    merged_system_prompt = f"{CRITIC_SYSTEM_PROMPT}\n\n{code_system_prompt}"
    response = llm.invoke([
        SystemMessage(content=merged_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content


def check_docs(content: list, query: str, expected_output: str = None, 
              step_context: dict = None) -> str:
    """审核文献"""
    docs_str = "\n".join(content) if isinstance(content, list) else str(content)
    
    docs_system_prompt = """
    你是一个科研审稿人。
    请检查文献是否与问题相关且包含足够信息。
    对于分步任务，只要满足当前步骤的检索要求即可。
    """
    
    step_context_note = ""
    if step_context:
        step_name = step_context.get("step_name", "")
        step_num = step_context.get("step_num", "")
        total_steps = step_context.get("total_steps", "")
        step_context_note = f"\n\n【步骤上下文】步骤 {step_num}/{total_steps}: {step_name}。\n"
    
    expected_output_note = ""
    if expected_output:
        expected_output_note = f"\n\n【验收标准】\n{expected_output}"
    
    user_prompt = f"""
    用户问题: {query}
    检索到的文献: {docs_str}
    {step_context_note}{expected_output_note}
    """
    
    merged_system_prompt = f"{CRITIC_SYSTEM_PROMPT}\n\n{docs_system_prompt}"
    response = llm.invoke([
        SystemMessage(content=merged_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content


def check_db(content: str, query: str, expected_output: str = None, 
            step_context: dict = None) -> str:
    """审核数据库结果"""
    db_system_prompt = """
    你是一个数据分析师。
    请检查数据查询结果是否为空，以及是否符合当前步骤的要求。
    """
    
    step_context_note = ""
    if step_context:
        step_name = step_context.get("step_name", "")
        step_num = step_context.get("step_num", "")
        total_steps = step_context.get("total_steps", "")
        step_context_note = f"\n\n【步骤上下文】步骤 {step_num}/{total_steps}: {step_name}。\n"
    
    expected_output_note = ""
    if expected_output:
        expected_output_note = f"\n\n【验收标准】\n{expected_output}"
    
    user_prompt = f"""
    用户问题: {query}
    数据库查询结果: {content}
    {step_context_note}{expected_output_note}
    """
    
    merged_system_prompt = f"{CRITIC_SYSTEM_PROMPT}\n\n{db_system_prompt}"
    response = llm.invoke([
        SystemMessage(content=merged_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content


def review_contribution(state: CriticAgentState) -> CriticAgentState:
    """
    审核节点
    """
    pending = state.get("pending_contribution")
    query = state.get("user_query", "")
    last_worker = state.get("last_worker", "")
    
    expected_output = state.get("current_step_expected_output")
    plan = state.get("plan", [])
    current_step_index = state.get("current_step_index", 0)
    
    step_context = None
    if plan and current_step_index < len(plan):
        current_step = plan[current_step_index]
        step_context = {
            "step_name": current_step.name,
            "step_num": str(current_step_index + 1),
            "total_steps": str(len(plan)),
            "step_description": current_step.description
        }
    
    print(f"--- [Critic] 正在审核 {last_worker} 的产出 ---")
    if step_context:
        print(f"  --> 步骤 {step_context['step_num']}/{step_context['total_steps']}: {step_context['step_name']}")
    
    if pending is None:
        state["is_approved"] = False
        state["critique_feedback"] = "未找到待审核内容"
        return state
    
    feedback = ""
    
    # 检查是否是图片
    if isinstance(pending, dict) and (
        "umap_base64" in pending or "image_base64" in pending
    ):
        image_b64 = pending.get("umap_base64") or pending.get("image_base64")
        feedback = check_umap_image(image_b64, query, expected_output, step_context)
        state["content_type"] = "image"
    
    # 检查是否是代码 (增强版提取逻辑)
    elif isinstance(pending, dict) and "code" in pending:
        code = pending.get("code", "")
        
        execution_result = ""
        
        # 1. 优先检查 success 字段
        is_exec_success = pending.get("success", False)
        
        # 2. 强制提取错误信息
        if not is_exec_success:
            error_msg = pending.get("error", "")
            if not error_msg:
                error_msg = (
                    pending.get("output_tail")
                    or pending.get("output_display")
                    or pending.get("output")
                    or "Unknown execution error"
                )
            execution_result = f"EXECUTION FAILED (CRITICAL):\n{error_msg}"
        else:
            result_msg = pending.get("result", "")
            if not result_msg:
                result_msg = (
                    pending.get("output_display")
                    or pending.get("output_tail")
                    or pending.get("output")
                    or "Execution successful but no output."
                )
            execution_result = f"EXECUTION SUCCESS:\n{result_msg}"
            
        # 3. 打印调试信息，确保我们知道传给 LLM 的是什么
        print(f"  --> [Debug] 传给 Critic 的执行结果片段: {execution_result[:100]}...")
        
        skill_note = format_skill_for_critic(state.get("current_step_skill_id"))
        feedback = check_code(
            code, query, execution_result, expected_output, step_context, skill_note=skill_note
        )
        state["content_type"] = "code"
    
    # 其他类型检查...
    elif isinstance(pending, list) or (isinstance(pending, dict) and "docs" in str(pending)):
        docs = pending if isinstance(pending, list) else pending.get("docs", [])
        feedback = check_docs(docs, query, expected_output, step_context)
        state["content_type"] = "docs"
    
    elif isinstance(pending, str) or (isinstance(pending, dict) and "result" in str(pending)):
        content = pending if isinstance(pending, str) else str(pending.get("result", ""))
        feedback = check_db(content, query, expected_output, step_context)
        state["content_type"] = "db_result"
    
    else:
        content_str = str(pending)
        skill_note = format_skill_for_critic(state.get("current_step_skill_id"))
        feedback = check_code(
            content_str, query, None, expected_output, step_context, skill_note=skill_note
        )
        state["content_type"] = "code"
    
    # 判断是否通过
    is_pass = "PASS" in feedback.upper() or "通过" in feedback
    
    if is_pass:
        print(f"  --> 审核通过！")
        state["is_approved"] = True
        state["critique_feedback"] = None
        
        if last_worker == "code_dev":
            # 保留 dict，便于下游展示完整 output / result（勿 str() 丢失结构）
            state["code_solution"] = pending if isinstance(pending, dict) else str(pending)
            
            # 探索阶段结束前，收集探索结果
            if not state.get("data_exploration_done", True):
                res_val = pending.get("result", "") if isinstance(pending, dict) else str(pending)
                if res_val:
                    exp_results = state.get("data_exploration_results", [])
                    state["data_exploration_results"] = exp_results + [res_val]

        elif last_worker == "rag_researcher":
            if isinstance(pending, list):
                state["rag_context"] = "\n\n".join(pending)
            else:
                state["rag_context"] = str(pending)
        elif last_worker == "data_analyst":
            state["final_report"] = str(pending)
        
        state["pending_contribution"] = None
    else:
        print(f"  --> 审核驳回！意见: {feedback}")
        state["is_approved"] = False
        state["critique_feedback"] = feedback
    
    state["review_details"] = feedback
    return state


# 构建子图
workflow = StateGraph(CriticAgentState)
workflow.add_node("review_contribution", review_contribution)
workflow.add_edge(START, "review_contribution")
workflow.add_edge("review_contribution", END)
critic_agent_graph = workflow.compile()