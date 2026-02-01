"""
celltype.py

基于 Enrichr 和 LLM Expert 的双引擎细胞类型注释工具。
集成逻辑：
1. Enrichr (PanglaoDB) 查库
2. LLM 专家推理
3. 归一化 + 冲突仲裁 (LCA/Lineage Check)
"""

import json
import re
import requests
from dataclasses import dataclass, asdict
from typing import Any, Dict, List, Literal, Tuple, Union
from langchain_core.tools import tool


from langchain_core.messages import SystemMessage, HumanMessage
from src.core.llm import get_llm  # 使用新的 LLM 工厂

# 初始化 LLM 实例
llm = get_llm(model_name="qwen-plus", temperature=0.1)

# =========================
# 0. 数据结构与辅助函数
# =========================

@dataclass
class EngineResult:
    engine: Literal["enrichr", "llm_expert"]
    label: str
    confidence: float
    normalized_label: str = "unknown"
    lineage: str = "unknown"
    meta: Dict[str, Any] = None

def _clean_genes(gene_list: List[str]) -> List[str]:
    """清洗基因列表，去除非法字符"""
    out = []
    for g in gene_list:
        if not isinstance(g, str) or not g:
            continue
        g = g.strip().upper()
        # 只保留字母数字
        g = re.sub(r"[^A-Z0-9_\-\.]", "", g)
        if g:
            out.append(g)
    # 去重保持顺序
    seen = set()
    dedup = []
    for g in out:
        if g not in seen:
            seen.add(g)
            dedup.append(g)
    return dedup

def _http(method: str, url: str, **kwargs) -> Any:
    """带重试和错误处理的 HTTP 请求"""
    try:
        r = requests.request(method, url, timeout=30, **kwargs)
        r.raise_for_status()
        return r.json()
    except Exception as e:
        print(f"[DiffEngine] HTTP Error: {e}")
        return {}

# =========================
# 1. 引擎 A: Enrichr
# =========================

def enrichr_vote(genes: List[str]) -> EngineResult:
    """调用 Enrichr API 进行注释"""
    ENRICHR_URL = "https://maayanlab.cloud/Enrichr"
    
    clean_genes = _clean_genes(genes)
    if not clean_genes:
        return EngineResult("enrichr", "unknown", 0.0, meta={"reason": "empty_genes"})

    # 1. 上传基因列表
    try:
        payload = {"description": "marker_genes"}
        files = {"list": ("genes.txt", ("\n".join(clean_genes)).encode("utf-8"))}
        r = requests.post(f"{ENRICHR_URL}/addList", data=payload, files=files, timeout=30)
        r.raise_for_status()
        user_list_id = r.json()["userListId"]
    except Exception as e:
        return EngineResult("enrichr", "unknown", 0.0, meta={"error": str(e)})

    # 2. 获取富集结果 (PanglaoDB)
    library = "PanglaoDB_Augmented_2021"
    res = _http("GET", f"{ENRICHR_URL}/enrich", params={"userListId": str(user_list_id), "backgroundType": library})
    
    rows = res.get(library, [])
    if not rows:
        return EngineResult("enrichr", "unknown", 0.0, meta={"reason": "no_hits"})

    # 3. 计算置信度 (基于前两名的分离度)
    best = rows[0]
    term = str(best[1]).strip() # Cell Type Name
    
    def get_score(row):
        try: return float(row[4]) # Combined Score
        except: return 0.0

    best_s = get_score(rows[0])
    second_s = get_score(rows[1]) if len(rows) > 1 else 0.0
    
    conf = 0.0
    if best_s > 0:
        sep = (best_s - second_s) / best_s
        # 简单的 sigmoid-like 映射
        conf = max(0.0, min(1.0, 0.4 + 0.6 * sep))

    return EngineResult("enrichr", term, conf, meta={"library": library, "top_score": best_s})

# =========================
# 2. 引擎 B: LLM Expert
# =========================

def llm_expert_vote(genes: List[str]) -> EngineResult:
    """让 LLM 扮演专家进行注释"""
    prompt = SystemMessage(content="""
You are a single-cell annotation expert.
Input: A list of marker genes.
Output: JSON only.
{
  "label": "most probable cell type",
  "confidence": 0.0 to 1.0,
  "evidence": "brief reason why"
}
Do not output extra text.
""")
    
    user_msg = HumanMessage(content=json.dumps({"genes": genes}, ensure_ascii=False))
    
    try:
        resp = llm.invoke([prompt, user_msg])
        content = resp.content.replace("```json", "").replace("```", "").strip()
        data = json.loads(content)
        return EngineResult(
            "llm_expert", 
            data.get("label", "unknown"), 
            float(data.get("confidence", 0.0)),
            meta={"evidence": data.get("evidence")}
        )
    except Exception as e:
        return EngineResult("llm_expert", "unknown", 0.0, meta={"error": str(e)})

# =========================
# 3. 仲裁与归一化
# =========================

def llm_normalize_and_lineage(label: str) -> Tuple[str, str]:
    """归一化名称并判断谱系"""
    if label.lower() == "unknown":
        return "unknown", "unknown"
        
    sys_msg = SystemMessage(content="""
Normalize the cell type label and identify its lineage.
Output JSON only:
{
  "normalized_label": "standard name (e.g., T cell, Macrophage)",
  "lineage": "immune | epithelial | stromal | neural | other | unknown"
}
""")
    try:
        resp = llm.invoke([sys_msg, HumanMessage(content=label)])
        data = json.loads(resp.content.replace("```json", "").replace("```", "").strip())
        return (
            data.get("normalized_label", "unknown"), 
            data.get("lineage", "unknown").lower()
        )
    except:
        return label, "unknown"

def llm_lca_check(labels: List[str]) -> Tuple[str, bool]:
    """查找最小公共祖先 (LCA)"""
    sys_msg = SystemMessage(content="""
Find the Lowest Common Ancestor (LCA) for these cell types.
Output JSON only:
{
  "lca_label": "broad category",
  "same_lineage": true/false
}
""")
    try:
        resp = llm.invoke([sys_msg, HumanMessage(content=json.dumps(labels))])
        data = json.loads(resp.content.replace("```json", "").replace("```", "").strip())
        return data.get("lca_label", "unknown"), data.get("same_lineage", False)
    except:
        return "unknown", False

# =========================
# 4. 主执行逻辑 (工具入口)
# =========================

@tool
def run_celltype_annotation(args: Dict[str, Any]) -> Dict[str, Any]:
    """
    [Tool Function] 执行细胞类型注释
    Args:
        args: 包含 'gene_list' (list or str)
    """
    raw_input = args.get("gene_list")
    
    # 1. 解析输入
    gene_list = []
    if isinstance(raw_input, list):
        gene_list = raw_input
    elif isinstance(raw_input, str):
        # 简单分割，如果更复杂可以用 enrichment.py 里的 flatten
        gene_list = [x.strip() for x in raw_input.split(",")]
    
    gene_list = _clean_genes(gene_list)
    print(f"[CellType Tool] Annotating {len(gene_list)} genes: {gene_list[:5]}...")

    if len(gene_list) < 2:
        return {"error": "Too few valid genes provided. Need at least 2."}

    # 2. 双引擎运行
    res_enrichr = enrichr_vote(gene_list)
    res_llm = llm_expert_vote(gene_list)

    # 3. 归一化
    res_enrichr.normalized_label, res_enrichr.lineage = llm_normalize_and_lineage(res_enrichr.label)
    res_llm.normalized_label, res_llm.lineage = llm_normalize_and_lineage(res_llm.label)

    print(f"  -> Enrichr says: {res_enrichr.normalized_label} ({res_enrichr.confidence:.2f})")
    print(f"  -> LLM says:     {res_llm.normalized_label} ({res_llm.confidence:.2f})")

    # 4. 决策逻辑
    final_label = "unknown"
    decision_path = "unresolved"

    # A. 完全一致
    if res_enrichr.normalized_label.lower() == res_llm.normalized_label.lower():
        final_label = res_enrichr.normalized_label
        decision_path = "consensus"
    
    # B. 置信度压制
    elif res_enrichr.confidence >= 0.75 and res_llm.confidence <= 0.6:
        final_label = res_enrichr.normalized_label
        decision_path = "enrichr_high_conf"
    elif res_llm.confidence >= 0.8 and res_enrichr.confidence <= 0.6:
        final_label = res_llm.normalized_label
        decision_path = "llm_high_conf"
    
    # C. 高置信度冲突 -> LCA 降级
    elif res_enrichr.confidence > 0.6 and res_llm.confidence > 0.6:
        lca, same = llm_lca_check([res_enrichr.normalized_label, res_llm.normalized_label])
        if same and lca != "unknown":
            final_label = lca
            decision_path = "conflict_resolved_by_LCA"
        else:
            final_label = "ambiguous_conflict"
            decision_path = "cross_lineage_conflict"
    else:
        # 低置信度回退
        if res_enrichr.confidence > res_llm.confidence:
            final_label = f"{res_enrichr.normalized_label} (Low Conf)"
        else:
            final_label = f"{res_llm.normalized_label} (Low Conf)"
        decision_path = "low_conf_guess"

    return {
        "final_prediction": final_label,
        "decision_logic": decision_path,
        "details": {
            "enrichr": asdict(res_enrichr),
            "llm_expert": asdict(res_llm)
        }
    }
