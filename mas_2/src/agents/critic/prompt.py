def get_umap_image_system_prompt() -> str:
    image_system_prompt = """
    --- Visualization Review Task ---
    Task: Evaluate the scientific visualization quality and relevance to the User Question.
    
    [Universal Criteria] (Must Have)
    1. Labels: Axis labels must be visible.
    2. Clarity: No severe blurring.
    
    [Step-Specific Context]
    If this is an intermediate step, do not demand final publication polish.
    """
    return f"{CRITIC_SYSTEM_PROMPT}\n{image_system_prompt}"

def get_umap_image_user_prompt(query: str, step_context_note: str, expected_output_note: str) -> str:
    return f"User question: {query}{step_context_note}{expected_output_note}"

def get_code_system_prompt() -> str:
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
    return f"{CRITIC_SYSTEM_PROMPT}\n\n{code_system_prompt}"

def get_code_user_prompt(query: str, content: str, execution_note: str, step_context_note: str, expected_output_note: str, skill_extra: str) -> str:
    parts = [f"用户问题: {query}", f"待审核代码: \n```python\n{content}\n```"]
    if execution_note:
        parts.append(execution_note.strip())
    if step_context_note:
        parts.append(step_context_note.strip())
    if expected_output_note:
        parts.append(expected_output_note.strip())
    if skill_extra:
        parts.append(skill_extra.strip())
    return "\n\n".join(parts)

def get_docs_system_prompt() -> str:
    docs_system_prompt = """
    你是一个科研审稿人。
    请检查文献是否与问题相关且包含足够信息。
    对于分步任务，只要满足当前步骤的检索要求即可。
    """
    return f"{CRITIC_SYSTEM_PROMPT}\n\n{docs_system_prompt}"

def get_docs_user_prompt(query: str, docs_str: str, step_context_note: str, expected_output_note: str) -> str:
    parts = [f"用户问题: {query}", f"检索到的文献: {docs_str}"]
    if step_context_note:
        parts.append(step_context_note.strip())
    if expected_output_note:
        parts.append(expected_output_note.strip())
    return "\n\n".join(parts)

def get_db_system_prompt() -> str:
    db_system_prompt = """
    你是一个数据分析师。
    请检查数据查询结果是否为空，以及是否符合当前步骤的要求。
    """
    return f"{CRITIC_SYSTEM_PROMPT}\n\n{db_system_prompt}"

def get_db_user_prompt(query: str, content: str, step_context_note: str, expected_output_note: str) -> str:
    parts = [f"用户问题: {query}", f"数据库查询结果: {content}"]
    if step_context_note:
        parts.append(step_context_note.strip())
    if expected_output_note:
        parts.append(expected_output_note.strip())
    return "\n\n".join(parts)

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
