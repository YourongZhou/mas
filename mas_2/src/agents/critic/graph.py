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

# 初始化 LLM
llm = get_llm(model_name="qwen-plus", temperature=0.5)
llm_vision = get_llm(model_name="qwen-vl-plus", temperature=0.2)

# --- 全局 System Prompt --- 
CRITIC_SYSTEM_PROMPT = """
你是一个资深的 AI 工作成果审核专家。
你的任务是对其他 AI 代理的工作成果进行严格审核。
审核标准必须严格、一致、可复现。
如果审核通过，请只回复 "PASS"。
如果审核不通过，请给出具体、可操作的修改建议。
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


def check_umap_image(image_base64: str, query: str) -> str:
    """审核 UMAP 图片质量"""
    umap_system_prompt = """
    You are an expert in embedding visualization review.
    Task: Determine if the UMAP plot shows good clustering.

    Criteria for PASS:
    - clusters are compact and visually separable
    - limited overlap between clusters
    - reasonable outliers (not dominating the plot)
    - not a single indistinct blob
    """
    
    user_prompt = f"User question: {query}"
    
    try:
        data_url = _normalize_base64_image(image_base64)
    except ValueError as exc:
        return f"INVALID_IMAGE: {exc}"

    message = HumanMessage(
        content=[
            {"type": "text", "text": CRITIC_SYSTEM_PROMPT},
            {"type": "text", "text": umap_system_prompt},
            {"type": "text", "text": user_prompt},
            {"type": "image_url", "image_url": {"url": data_url}},
        ]
    )
    response = llm_vision.invoke([message])
    return response.content


def check_code(content: str, query: str) -> str:
    """审核代码"""
    code_system_prompt = """
    你是一个资深代码审查员。
    请检查：
    1. 代码是否安全？
    2. 是否直接回答了问题？
    3. 是否有死循环风险？
    4. 代码逻辑是否合理？
    """
    
    user_prompt = f"""
    用户问题: {query}
    待审核代码: {content}
    """
    
    response = llm.invoke([
        SystemMessage(content=CRITIC_SYSTEM_PROMPT),
        SystemMessage(content=code_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content


def check_docs(content: list, query: str) -> str:
    """审核文献"""
    docs_str = "\n".join(content) if isinstance(content, list) else str(content)
    
    docs_system_prompt = """
    你是一个科研审稿人。
    请检查：
    1. 文献是否与问题强相关？
    2. 是否包含足够的信息？
    3. 文献质量是否可靠？
    """
    
    user_prompt = f"""
    用户问题: {query}
    检索到的文献: {docs_str}
    """
    
    response = llm.invoke([
        SystemMessage(content=CRITIC_SYSTEM_PROMPT),
        SystemMessage(content=docs_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content


def check_db_result(content: str, query: str) -> str:
    """审核数据库结果"""
    db_system_prompt = """
    你是一个数据分析师。
    请检查：
    1. 数据格式是否正确？
    2. 是否为空结果？
    3. 数据是否与问题相关？
    """
    
    user_prompt = f"""
    用户问题: {query}
    数据库查询结果: {content}
    """
    
    response = llm.invoke([
        SystemMessage(content=CRITIC_SYSTEM_PROMPT),
        SystemMessage(content=db_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content


def review_contribution(state: CriticAgentState) -> CriticAgentState:
    """
    审核节点
    根据内容类型调用相应的审核函数
    """
    pending = state.get("pending_contribution")
    query = state.get("user_query", "")
    last_worker = state.get("last_worker", "")
    
    print(f"--- [Critic] 正在审核 {last_worker} 的产出 ---")
    
    if pending is None:
        feedback = "未找到待审核内容"
        state["is_approved"] = False
        state["critique_feedback"] = feedback
        return state
    
    # 判断内容类型并分发审核
    feedback = ""
    
    # 检查是否是图片
    if isinstance(pending, dict) and (
        "umap_base64" in pending or "image_base64" in pending
    ):
        image_b64 = pending.get("umap_base64") or pending.get("image_base64")
        feedback = check_umap_image(image_b64, query)
        state["content_type"] = "image"
    
    # 检查是否是代码
    elif isinstance(pending, dict) and "code" in pending:
        code = pending.get("code", "")
        feedback = check_code(code, query)
        state["content_type"] = "code"
    
    # 检查是否是文档列表
    elif isinstance(pending, list) or (isinstance(pending, dict) and "docs" in str(pending)):
        docs = pending if isinstance(pending, list) else pending.get("docs", [])
        feedback = check_docs(docs, query)
        state["content_type"] = "docs"
    
    # 检查是否是数据库结果
    elif isinstance(pending, str) or (isinstance(pending, dict) and "result" in str(pending)):
        content = pending if isinstance(pending, str) else str(pending.get("result", ""))
        feedback = check_db_result(content, query)
        state["content_type"] = "db_result"
    
    else:
        # 通用审核
        content_str = str(pending)
        if "code" in content_str.lower() or "def " in content_str or "import " in content_str:
            feedback = check_code(content_str, query)
            state["content_type"] = "code"
        else:
            feedback = check_docs([content_str], query)
            state["content_type"] = "docs"
    
    # 判断是否通过
    is_pass = "PASS" in feedback.upper() or "通过" in feedback
    
    if is_pass:
        print(f"  --> 审核通过！")
        state["is_approved"] = True
        state["critique_feedback"] = None
        
        # 根据 worker 类型归档数据
        if last_worker == "code_dev":
            state["code_solution"] = str(pending)
        elif last_worker == "rag_researcher":
            # 如果是列表，转换为字符串
            if isinstance(pending, list):
                state["rag_context"] = "\n\n".join(pending)
            else:
                state["rag_context"] = str(pending)
        elif last_worker == "data_analyst":
            state["final_report"] = str(pending)
        
        # 清空待审核区
        state["pending_contribution"] = None
    else:
        print(f"  --> 审核驳回！意见: {feedback}")
        state["is_approved"] = False
        state["critique_feedback"] = feedback
    
    state["review_details"] = feedback
    return state


# 构建子图
workflow = StateGraph(CriticAgentState)

# 添加节点
workflow.add_node("review_contribution", review_contribution)

# 定义边
workflow.add_edge(START, "review_contribution")
workflow.add_edge("review_contribution", END)

# 编译子图
critic_agent_graph = workflow.compile()

