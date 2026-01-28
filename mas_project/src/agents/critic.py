import base64
import re
import os

os.environ["LANGCHAIN_TRACING_V2"] = "true"

from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage, HumanMessage
from src.schema import GlobalState

# --- 全局 System Prompt --- 
CRITIC_SYSTEM_PROMPT = """
你是一个资深的 AI 工作成果审核专家。
你的任务是对其他 AI 代理的工作成果进行严格审核。
审核标准必须严格、一致、可复现。
如果审核通过，请只回复 "PASS"。
如果审核不通过，请给出具体、可操作的修改建议。
"""

llm = ChatOpenAI(
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
    model="qwen-plus",
    temperature=0.5
)

# Vision model for UMAP image review
llm_vision = ChatOpenAI(
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
    model="qwen-vl-plus",
    temperature=0.2
)

def _normalize_base64_image(image_b64: str, default_mime: str = "image/png") -> str:
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

# --- Sub-feature 4: Review UMAP image quality ---
def check_umap_image(image_base64: str, query: str):
    # 子角色 System Prompt
    umap_system_prompt = """
    You are an expert in embedding visualization review.
    Task: Determine if the UMAP plot shows good clustering.

    Criteria for PASS:
    - clusters are compact and visually separable
    - limited overlap between clusters
    - reasonable outliers (not dominating the plot)
    - not a single indistinct blob
    """
    
    # User Prompt
    user_prompt = f"""
    User question: {query}
    """
    
    try:
        data_url = _normalize_base64_image(image_base64)
    except ValueError as exc:
        return f"INVALID_IMAGE: {exc}"

    # 发送给模型
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

# --- 子功能 1: 审核代码 ---
def check_code(content: str, query: str):
    # 子角色 System Prompt
    code_system_prompt = """
    你是一个资深代码审查员。
    请检查：
    1. 代码是否安全？
    2. 是否直接回答了问题？
    3. 是否有死循环风险？
    """
    
    # User Prompt
    user_prompt = f"""
    用户问题: {query}
    待审核代码: {content}
    """
    
    # 发送给模型
    response = llm.invoke([
        SystemMessage(content=CRITIC_SYSTEM_PROMPT),
        SystemMessage(content=code_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content

# --- 子功能 2: 审核文献 ---
def check_docs(content: list, query: str):
    # content 是 list，转成字符串给 LLM 看
    docs_str = "\n".join(content)
    
    # 子角色 System Prompt
    docs_system_prompt = """
    你是一个科研审稿人。
    请检查：
    1. 文献是否与问题强相关？
    2. 是否包含足够的信息？
    """
    
    # User Prompt
    user_prompt = f"""
    用户问题: {query}
    检索到的文献: {docs_str}
    """
    
    # 发送给模型
    response = llm.invoke([
        SystemMessage(content=CRITIC_SYSTEM_PROMPT),
        SystemMessage(content=docs_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content

# --- 子功能 3: 审核数据库结果 ---
def check_db(content: str, query: str):
    # 子角色 System Prompt
    db_system_prompt = """
    你是一个数据分析师。
    请检查：
    1. 数据格式是否正确？
    2. 是否为空结果？
    """
    
    # User Prompt
    user_prompt = f"""
    用户问题: {query}
    数据库查询结果: {content}
    """
    
    # 发送给模型
    response = llm.invoke([
        SystemMessage(content=CRITIC_SYSTEM_PROMPT),
        SystemMessage(content=db_system_prompt),
        HumanMessage(content=user_prompt)
    ])
    return response.content

# --- 主入口节点 ---
def critic_node(state: GlobalState):
    worker = state["active_worker"]
    pending = state["pending_result"]
    query = state["user_query"]
    
    print(f"--- [Critic] 正在审核 {worker} 的产出 ---")
    
    # 1. 分发逻辑 (Dispatcher)
    if isinstance(pending, dict) and (
        "umap_base64" in pending or "image_base64" in pending
    ):
        image_b64 = pending.get("umap_base64") or pending.get("image_base64")
        feedback = check_umap_image(image_b64, query)
    elif worker == "code_agent":
        feedback = check_code(pending, query)
    elif worker == "rag_agent":
        feedback = check_docs(pending, query)
    elif worker == "tool_agent":
        feedback = check_db(pending, query)
    else:
        feedback = "未知 Worker类型，无法审核"

    # 2. 判断是否通过
    if "PASS" in feedback:
        print(f"  --> 审核通过！正在归档数据...")
        
        # 【关键步骤】数据转正 (Draft -> Final)
        updates = {
            "is_approved": True,
            "critique_feedback": None, # 清空旧的批评
            "pending_result": None     # 清空草稿区
        }
        
        # 根据是谁干的活，把数据搬运到对应的正式字段
        if worker == "code_agent":
            updates["code_result"] = pending
        elif worker == "rag_agent":
            updates["rag_docs"] = pending
        elif worker == "tool_agent":
            updates["db_results"] = pending
            
        return updates
        
    else:
        print(f"  --> 审核驳回！意见: {feedback}")
        return {
            "is_approved": False,
            "critique_feedback": feedback
            # 注意：pending_result 保留不动，或者清空都可以，取决于Worker是否需要参照旧草稿
        }
