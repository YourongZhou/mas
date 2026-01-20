from typing import TypedDict, List, Optional, Literal, Any

class GlobalState(TypedDict):
    messages: List[str]
    user_query: str
    
    # --- 1. 正式归档区 (Supervisor 做决策的依据) ---
    # 只有审核通过的数据才会出现在这里
    rag_docs: List[str]
    db_results: str
    code_result: str
    
    # --- 2. 调度控制区 ---
    # 当前是谁在干活？(rag_agent / code_agent / tool_agent)
    active_worker: Literal["rag_agent", "code_agent", "tool_agent"]
    # Supervisor 的下一步指令
    next_action: Literal["continue", "FINISH"]
    
    # --- 3. 待审核草稿区 (Worker 的输出) ---
    # 这是一个临时缓冲区，审核通过前，这里的数据是“脏”的
    pending_result: Any 
    
    # --- 4. 审核反馈 ---
    is_approved: bool
    critique_feedback: Optional[str]
    final_answer: str