import os
from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage
from langgraph.graph import StateGraph, MessagesState, START, END

# llm setup
# 建议：将 BASE_URL 也放入环境变量，或者像这样保留
llm = ChatOpenAI(
    # 不要硬编码 API Key，会自动读取环境变量 OPENAI_API_KEY
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
    model="qwen-plus",
    temperature=0.5
)

# 算数 agent
def math_agent(state: MessagesState):
    system_prompt = SystemMessage(content="""
    你是一个普通的计算员，会做基本的数学运算。
    请根据用户的问题或者审核员的建议，直接给出答案。不要给出计算过程。
    """)
    messages = [system_prompt] + state["messages"]
    response = llm.invoke(messages)
    response.name = "math_agent"
    return {"messages": [response]}

# 检查 agent
def check_agent(state: MessagesState):
    system_prompt = SystemMessage(content="""
    你是一个数学家，擅长检查各种数学问题。
    请检查最近一次 math_agent 给出的答案是否正确。
    如果正确，请给出 "correct"，并确保回答中没有"incorrect"这个词。
    如果错误，请给出 "incorrect"，并给出新的计算建议。
    """)
    messages = [system_prompt] + state["messages"]
    response = llm.invoke(messages)
    response.name = "check_agent"
    return {"messages": [response]}

def should_continue(state: MessagesState):
    last_message = state["messages"][-1]
    # 稍微放宽判断条件，防止大小写问题
    content = last_message.content.lower()
    if 'incorrect' in content:
        return "math_agent"
    else:
        return "END"

# --- 构建图 ---

workflow = StateGraph(MessagesState)

workflow.add_node("math_agent", math_agent)
workflow.add_node("check_agent", check_agent)

workflow.add_edge(START, "math_agent")
workflow.add_edge("math_agent", "check_agent")

workflow.add_conditional_edges(
    "check_agent",
    should_continue,
    {
        "math_agent": "math_agent",
        "END": END,
    },
)

# 【关键】这里必须赋值给一个变量，供 Studio 导入
# 通常命名为 graph 或 app
graph = workflow.compile()