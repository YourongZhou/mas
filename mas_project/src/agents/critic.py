from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage
from src.schema import GlobalState

llm = ChatOpenAI(
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
    model="qwen-plus",
    temperature=0.5
)

# --- 子功能 1: 审核代码 ---
def check_code(content: str, query: str):
    prompt = f"""
    你是一个资深代码审查员。
    用户问题: {query}
    待审核代码: {content}
    
    请检查：
    1. 代码是否安全？
    2. 是否直接回答了问题？
    3. 是否有死循环风险？
    
    如果不通过，请给出具体修改建议。
    如果通过，请只回复 "PASS"。
    """
    response = llm.invoke([HumanMessage(content=prompt)])
    return response.content

# --- 子功能 2: 审核文献 ---
def check_docs(content: list, query: str):
    # content 是 list，转成字符串给 LLM 看
    docs_str = "\n".join(content)
    prompt = f"""
    你是一个科研审稿人。
    用户问题: {query}
    检索到的文献: {docs_str}
    
    请检查：
    1. 文献是否与问题强相关？
    2. 是否包含足够的信息？
    
    如果不通过，请建议关键词修改。
    如果通过，请只回复 "PASS"。
    """
    response = llm.invoke([HumanMessage(content=prompt)])
    return response.content

# --- 子功能 3: 审核数据库结果 ---
def check_db(content: str, query: str):
    prompt = f"""
    你是一个数据分析师。
    用户问题: {query}
    数据库查询结果: {content}
    
    请检查：
    1. 数据格式是否正确？
    2. 是否为空结果？
    
    如果不通过，请说明原因。
    如果通过，请只回复 "PASS"。
    """
    response = llm.invoke([HumanMessage(content=prompt)])
    return response.content

# --- 主入口节点 ---
def critic_node(state: GlobalState):
    worker = state["active_worker"]
    pending = state["pending_result"]
    query = state["user_query"]
    
    print(f"--- [Critic] 正在审核 {worker} 的产出 ---")
    
    # 1. 分发逻辑 (Dispatcher)
    if worker == "code_agent":
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