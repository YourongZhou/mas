"""RAG Worker：接入外部向量库（Chroma 持久化），返回 top-k 相关片段。"""

import os
from functools import lru_cache
from typing import List

import chromadb
from langgraph.graph import StateGraph, START, END
from sentence_transformers import SentenceTransformer

from src.schema import GlobalState

# ---------------------------
# 配置
# ---------------------------

# 持久化向量库路径与集合名称，可通过环境变量覆盖
CHROMA_PERSIST_PATH = os.getenv("CHROMA_PERSIST_PATH", "./chroma_db")
CHROMA_COLLECTION = os.getenv("CHROMA_COLLECTION", "docs")

# 嵌入模型名称；需与向量库中入库时使用的模型保持一致
EMBED_MODEL_NAME = os.getenv("EMBED_MODEL_NAME", "sentence-transformers/all-MiniLM-L6-v2")
TOP_K = int(os.getenv("RAG_TOP_K", "5"))


# ---------------------------
# 嵌入与向量库
# ---------------------------

@lru_cache(maxsize=1)
def _get_embedder():
    return SentenceTransformer(EMBED_MODEL_NAME)


class _STEmbeddingFunction:
    def __init__(self, model):
        self.model = model

    def __call__(self, texts: List[str]):
        return self.model.encode(texts, normalize_embeddings=True).tolist()


@lru_cache(maxsize=1)
def _get_collection():
    embedder = _get_embedder()
    client = chromadb.PersistentClient(path=CHROMA_PERSIST_PATH)
    embedding_fn = _STEmbeddingFunction(embedder)
    return client.get_or_create_collection(
        name=CHROMA_COLLECTION,
        embedding_function=embedding_fn,
    )


def _vector_search(query: str, top_k: int = TOP_K) -> List[str]:
    collection = _get_collection()
    result = collection.query(query_texts=[query], n_results=top_k)

    docs = result.get("documents") or []
    # chromadb 返回 list[list[str]] 结构
    flat = docs[0] if docs else []
    if not flat:
        return ["未找到相关文献片段"]
    return flat


# ---------------------------
# LangGraph 节点
# ---------------------------

def search_step(state: GlobalState):
    query = state["user_query"]
    print(f"--- [RAG组] 正在查询向量库: {query} ---")

    docs = _vector_search(query, top_k=TOP_K)
    return {"pending_result": docs}


# 构建子图
workflow = StateGraph(GlobalState)
workflow.add_node("search", search_step)
workflow.add_edge(START, "search")
workflow.add_edge("search", END)

# 导出编译好的子图
rag_app = workflow.compile()