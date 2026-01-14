from langchain_openai import ChatOpenAI
from langchain_core.messages import HumanMessage, SystemMessage
from langgraph.graph import StateGraph, MessagesState, START, END
import os

# langsmith setup
os.environ["LANGCHAIN_TRACING_V2"] = "true"
os.environ["LANGCHAIN_API_KEY"] = "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"
os.environ["LANGCHAIN_PROJECT"] = "chatbot_test"

# llm setup
llm = ChatOpenAI(
    # 根据自己的需求配置，可以是环境变量，也可以是文本内容
    api_key="sk-3e43aba7e80343bc96fb5e7d549837ac", 
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
    model="qwen-plus"
)

# 算数 agent
def math_agent(state: MessagesState):
    system_prompt = SystemMessage("""
    你是一个普通的计算员，会做基本的数学运算。
    请根据用户的问题或者审核员的建议，直接给出答案。不要给出计算过程。
    """)
    messages = [system_prompt] + state["messages"]
    response = llm.invoke(messages)

    # 标记 response 名称
    response.name = "math_agent"
    return {"messages": [response]}


# 检查 agent
def check_agent(state: MessagesState):
    system_prompt = SystemMessage("""
    你是一个数学家，擅长检查各种数学问题。
    请检查最近一次 math_agent 给出的答案是否正确。
    如果正确，请给出 "correct"，并确保回答中没有"incorrect"这个词。
    如果错误，请给出新的计算建议，请给出 "incorrect"。
    """)
    messages = [system_prompt] + state["messages"]
    response = llm.invoke(messages)

    # 标记 response 名称
    response.name = "check_agent"
    return {"messages": [response]}

def should_continue(state: MessagesState):
    """
    根据审核结果，判断是否重新计算。
    返回下一个节点的名称。
    """
    last_message = state["messages"][-1]
    if 'incorrect' in last_message.content:
        print("审核不通过，打回中...")
        return "math_agent"
    else:
        print("审核通过，返回结果...")
        return "END"

# 构建图
workflow = StateGraph(MessagesState)

# 添加节点
workflow.add_node("math_agent", math_agent)
workflow.add_node("check_agent", check_agent)

# 添加边
workflow.add_edge(START, "math_agent")
workflow.add_edge("math_agent", "check_agent")
# 复杂判断
workflow.add_conditional_edges(
    "check_agent",
    should_continue,
    {
        "math_agent": "math_agent",
        "END": END,
    },
)
workflow = workflow.compile()

workflow.invoke({"messages": [{"role": "user", "content": "what is 3^{12345} (mod 100)?"}]})

# 生成图的流程图片
qa_async_ascii = workflow.get_graph().print_ascii()
print(qa_async_ascii)