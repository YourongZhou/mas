"""
RAG Researcher Agent 子图
负责向量检索和文档查询
"""
import os
from functools import lru_cache
from typing import List
from langgraph.graph import StateGraph, START, END
from .state import RAGAgentState

try:
    import chromadb
    from sentence_transformers import SentenceTransformer
    CHROMA_AVAILABLE = True
except ImportError:
    CHROMA_AVAILABLE = False
    print("警告: chromadb 或 sentence-transformers 未安装，RAG 功能将受限")

# ---------------------------
# 配置
# ---------------------------

CHROMA_PERSIST_PATH = os.getenv("CHROMA_PERSIST_PATH", "./chroma_db")
CHROMA_COLLECTION = os.getenv("CHROMA_COLLECTION", "docs")
EMBED_MODEL_NAME = os.getenv("EMBED_MODEL_NAME", "sentence-transformers/all-MiniLM-L6-v2")
TOP_K = int(os.getenv("RAG_TOP_K", "5"))


# ---------------------------
# 嵌入与向量库
# ---------------------------

@lru_cache(maxsize=1)
def _get_embedder():
    if not CHROMA_AVAILABLE:
        return None
    return SentenceTransformer(EMBED_MODEL_NAME)


class _STEmbeddingFunction:
    def __init__(self, model):
        self.model = model

    def __call__(self, texts: List[str]):
        return self.model.encode(texts, normalize_embeddings=True).tolist()


@lru_cache(maxsize=1)
def _get_collection():
    if not CHROMA_AVAILABLE:
        return None
    
    embedder = _get_embedder()
    if embedder is None:
        return None
    
    client = chromadb.PersistentClient(path=CHROMA_PERSIST_PATH)
    embedding_fn = _STEmbeddingFunction(embedder)
    return client.get_or_create_collection(
        name=CHROMA_COLLECTION,
        embedding_function=embedding_fn,
    )


def _vector_search(query: str, top_k: int = TOP_K) -> List[str]:
    """向量检索函数"""
    collection = _get_collection()
    if collection is None:
        return ["向量库未配置，请检查 chromadb 和 sentence-transformers 是否已安装"]
    
    try:
        result = collection.query(query_texts=[query], n_results=top_k)
        docs = result.get("documents") or []
        flat = docs[0] if docs else []
        if not flat:
            return ["未找到相关文献片段"]
        return flat
    except Exception as e:
        return [f"检索出错: {str(e)}"]


def search_documents(state: RAGAgentState) -> RAGAgentState:
    """
    检索文档节点
    从向量库中检索相关文档
    """
    # 优先使用当前步骤的输入作为查询内容
    current_step_input = state.get("current_step_input", "")
    current_step_expected_output = state.get("current_step_expected_output", "")
    
    # 构建查询：如果有当前步骤输入，使用它；否则使用search_query或user_query
    if current_step_input:
        query = current_step_input
        print(f"--- [RAG Researcher] 正在查询向量库（根据计划步骤）: {query[:100]}... ---")
    else:
        query = state.get("search_query") or state.get("user_query", "")
        print(f"--- [RAG Researcher] 正在查询向量库: {query} ---")
    
    # 如果有预期输出，在日志中显示
    if current_step_expected_output:
        print(f"  --> 预期输出要求: {current_step_expected_output[:100]}...")
    
    docs = _vector_search(query, top_k=TOP_K)
    
    # 更新状态
    state["retrieved_docs"] = docs
    state["doc_count"] = len(docs)
    
    # 将结果写入 pending_contribution
    state["pending_contribution"] = docs
    
    # 同时更新 rag_context
    state["rag_context"] = "\n\n".join(docs)
    
    print(f"  --> 检索到 {len(docs)} 条文档")
    
    return state


# 构建子图
workflow = StateGraph(RAGAgentState)

# 添加节点
workflow.add_node("search_documents", search_documents)

# 定义边
workflow.add_edge(START, "search_documents")
workflow.add_edge("search_documents", END)

# 编译子图
rag_agent_graph = workflow.compile()

