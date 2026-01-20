from langgraph.graph import StateGraph, START, END
from src.schema import GlobalState

# 模拟 RAG 搜索步骤
def search_step(state: GlobalState):
    query = state["user_query"]
    print(f"--- [RAG组] 正在搜索: {query} ---")
    
    # TODO: 这里替换为真实的 Vector DB 搜索逻辑
    mock_docs = [
        f"文献1:文献1中将量子纠缠定义为一种典型的非经典相关性，是量子系统中多个子系统之间所形成的整体态性质。当系统处于纠缠态时，其整体波函数无法表示为各子系统态的直积形式，即整体信息不能由局部信息完全决定。这意味着，即便各子系统在空间上被分离，对其中一个系统的测量结果仍会瞬时影响另一个系统的测量统计分布，这种现象突破了经典物理中局域实在论的直觉。爱因斯坦曾将其称为“鬼魅般的超距作用”，以表达其对经典因果观的挑战。从数学角度看，纠缠态通常通过密度矩阵的不可分解性、冯·诺依曼熵或纠缠熵等指标进行刻画。量子纠缠不仅是量子力学区别于经典力学的核心特征之一，也是量子信息科学的基础资源，在量子通信、量子计算与量子精密测量中发挥着不可替代的作用。",
        f"文献2:文献2系统总结了近年来量子纠缠在理论与实验层面的重要进展。在理论方面，研究者提出了多种新的纠缠度量方法，并在多体系统、拓扑量子态和开放量子系统中深入探讨了纠缠的动力学行为。这些工作有助于揭示复杂量子系统中的相变机制和信息传播规律。在实验层面，高维纠缠态、长距离纠缠分发以及卫星量子通信取得了突破性成果，使得纠缠从实验室走向实际应用成为可能。此外，纠缠还被广泛应用于量子误差校正、量子隐形传态和分布式量子计算等关键任务中，显著提升了量子信息处理的鲁棒性和效率。这些进展不仅推动了量子技术的工程化进程，也进一步加深了人类对量子世界基本规律的理解。",
    ]
    
    # 只更新 rag_docs 字段
    return {"pending_result": mock_docs}

# 构建子图
workflow = StateGraph(GlobalState)
workflow.add_node("search", search_step)
workflow.add_edge(START, "search")
workflow.add_edge("search", END)

# 导出编译好的子图
rag_app = workflow.compile()