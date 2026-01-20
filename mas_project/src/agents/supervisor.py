from langchain_openai import ChatOpenAI
from pydantic import BaseModel, Field
from typing import Literal
from src.schema import GlobalState

class RouteDecision(BaseModel):
    next_worker: Literal["rag_agent", "code_agent", "tool_agent", "FINISH"] = Field(...)
    reasoning: str = Field(...)

llm = ChatOpenAI(
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
    model="qwen-plus",
    temperature=0.5
)

def supervisor_node(state: GlobalState):
    # 构建 Prompt
    prompt = f"""
    你是项目经理。当前项目状态：
    - RAG文献数: {len(state.get('rag_docs') or [])}
    - DB结果: {state.get('db_results', '无')}
    - 代码/最终结果: {state.get('code_result', '无')}
    
    用户问题: {state['user_query']}
    
    请停止。
    """
    
    chain = llm.with_structured_output(RouteDecision)
    decision = chain.invoke(prompt)
    
    print(f"\n[Supervisor] 决策: {decision.next_worker} ({decision.reasoning})")
    
    # 更新状态
    if decision.next_worker == "FINISH":
        return {"next_action": "FINISH"}
    else:
        return {
            "next_action": "continue",
            "active_worker": decision.next_worker # 设定当前的活跃工人
        }