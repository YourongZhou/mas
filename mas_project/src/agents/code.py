from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage, HumanMessage
from langgraph.graph import StateGraph, START, END
from src.schema import GlobalState

llm = ChatOpenAI(
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
    model="qwen-plus",
    temperature=0.5
)

def coding_step(state: GlobalState):
    print("--- [Code组] 正在生成解决方案 ---")
    
    # 获取上下文
    docs = "\n".join(state.get("rag_docs", []))
    feedback = state.get("critique_feedback", "无")
    
    # 构建 Prompt：包含 RAG 资料 + 之前的审核意见
    prompt = f"""
    任务：基于参考资料回答用户问题。
    参考资料：{docs}
    上次审核意见：{feedback} (如果是"无"则忽略)
    
    用户问题：{state['user_query']}
    """
    
    response = llm.invoke([HumanMessage(content=prompt)])
    
    # 更新 code_result
    return {"pending_result": response.content}

# 构建子图
workflow = StateGraph(GlobalState)
workflow.add_node("generate_code", coding_step)
workflow.add_edge(START, "generate_code")
workflow.add_edge("generate_code", END)

code_app = workflow.compile()