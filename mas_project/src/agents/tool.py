# src/agents/tool.py

from langchain_openai import ChatOpenAI
from langchain_core.messages import SystemMessage, HumanMessage
from langchain_core.tools import tool
from langgraph.graph import StateGraph, MessagesState, START, END
from langgraph.prebuilt import ToolNode, tools_condition
from src.schema import GlobalState

# --- 1. 定义工具 (Mock 数据库) ---
# 真正的开发中，这里替换为 SQL 查询或 API 调用
@tool
def query_gene_database(gene_name: str):
    """根据基因名称查询 Gene Ontology 和相关信息。"""
    print(f"  --> [Tool] 正在查询数据库: {gene_name}")
    # 模拟返回结果
    if "p53" in gene_name.lower():
        return "Gene: TP53; Function: Tumor suppressor; Pathway: Apoptosis."
    elif "egfr" in gene_name.lower():
        return "Gene: EGFR; Function: Epidermal growth factor receptor."
    else:
        return "No record found."

@tool
def search_web(query: str):
    """搜索互联网以获取补充信息。"""
    print(f"  --> [Tool] 正在搜索: {query}")
    return f"Search results for {query}: Some relevant bio-papers..."

# 工具列表
tools = [query_gene_database, search_web]

# --- 2. 设置 LLM ---
llm = ChatOpenAI(
    base_url="https://dashscope.aliyuncs.com/compatible-mode/v1",
    model="qwen-plus",
    temperature=0.5
)
# 关键：把工具绑定到模型上
llm_with_tools = llm.bind_tools(tools)

# --- 3. 定义子图逻辑 ---

# 3.1 适配器：从 GlobalState 进入 Subgraph 的 MessagesState
def agent_node(state: MessagesState):
    """
    Agent 节点：负责思考和决定调用哪个工具。
    注意：这里的 input state 是 MessagesState (列表)，
    而不是 GlobalState (字典)。
    """
    messages = state["messages"]
    
    # 如果是第一次进入，messages 可能为空，或者我们手动注入了 System Prompt
    # 但为了简单，我们在下面的 entry_builder 里处理初始化
    
    response = llm_with_tools.invoke(messages)
    return {"messages": [response]}

# 3.2 工具节点：LangGraph 预置的，非常方便
tool_node = ToolNode(tools)

# --- 4. 构建子图 ---
# 注意：这个子图内部使用 MessagesState，只有 Messages
workflow = StateGraph(MessagesState)

workflow.add_node("agent", agent_node)
workflow.add_node("tools", tool_node)

workflow.add_edge(START, "agent")

# 条件边：如果 LLM 想调工具 -> tools；否则 -> END
workflow.add_conditional_edges(
    "agent",
    tools_condition,
)
workflow.add_edge("tools", "agent") # 工具跑完回 Agent 继续思考

# 编译内部子图
react_graph = workflow.compile()

# --- 5. 【核心】包装成主图可用的节点 ---
# 这一步是为了把 GlobalState (字典) 转换成 MessagesState (列表)，
# 跑完后再转回 GlobalState 更新 db_results。

def run_tool_agent(state: GlobalState):
    print("--- [Tool组] 启动数据库查询 Agent ---")
    
    # A. 准备输入：把用户的 Query 包装成第一条消息
    initial_msg = [
        SystemMessage(content="你是一个生物数据库专家。请利用工具查询必要信息。查到后请总结。"),
        HumanMessage(content=state["user_query"])
    ]
    
    # B. 调用上面的 ReAct 子图
    # invoke 的输入必须符合 MessagesState 的结构
    result = react_graph.invoke({"messages": initial_msg})
    
    # C. 提取结果
    # 拿到最后一条 AI 的回复（也就是工具跑完后的总结）
    final_message = result["messages"][-1]
    
    # D. 返回给主图
    return {"pending_result": final_message.content}