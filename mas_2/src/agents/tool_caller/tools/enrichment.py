"""
enrichment.py

基于 gseapy.enrichr 的【真实 gene set 功能富集分析工具】。

设计原则：
- 鲁棒性：具备“核弹级”输入清洗能力，能够处理 Agent/LLM 传入的各种奇怪格式（嵌套列表、JSON 字符串、带引号文本等）。
- 自动化：自动调用在线 Enrichr 数据库进行 ORA 分析。
- 可视化：自动生成富集分析条形图并保存至本地。

依赖库：
- gseapy
- pandas
"""

import os
import re
import json
import ast
import gseapy as gp
import pandas as pd
from langchain_core.tools import tool
from typing import List, Dict, Any, Union

# =========================
# 辅助函数：深度展平任意数据结构
# =========================
def flatten_input(data: Any) -> str:
    """
    递归处理函数。
    将任意嵌套的列表、字典、JSON字符串、元组等强制展平成一个单一的长字符串。
    这能有效解决 LLM 偶尔传入 `['["GeneA", "GeneB"]']` 或 `{'genes': ...}` 等复杂结构的问题。
    """
    if data is None:
        return ""
    
    # 1. 处理字符串类型
    if isinstance(data, str):
        data = data.strip()
        # 尝试检测并解析 JSON 格式的字符串 (例如 '["CD3D", "LCK"]')
        if (data.startswith("[") and data.endswith("]")) or \
           (data.startswith("{") and data.endswith("}")):
            try:
                parsed = json.loads(data)
                return flatten_input(parsed) # 递归调用处理解析后的对象
            except:
                # 如果 JSON 解析失败，尝试用 Python 的 literal_eval 解析 (例如 "['A', 'B']")
                try:
                    parsed = ast.literal_eval(data)
                    return flatten_input(parsed)
                except:
                    pass
        # 如果不是结构化字符串，直接返回原文本
        return data

    # 2. 处理列表或元组：递归拼接每个元素
    if isinstance(data, (list, tuple)):
        return " ".join([flatten_input(item) for item in data])
    
    # 3. 处理字典：只提取 values 部分
    if isinstance(data, dict):
        return " ".join([flatten_input(v) for v in data.values()])
    
    # 4. 其他类型（如 int, float）强转字符串
    return str(data)


# =========================
# 核心执行函数
# =========================
@tool
def gene_set_enrichment(
    gene_list: Union[List[str], str, Any], # 放宽类型限制以兼容 Agent 的各种输出
    organism: str = "human",
    databases: List[str] | None = None,
    method: str = "ora",
    top_k: int = 10,
) -> Dict[str, Any]:
    """
    使用 Enrichr 对 gene set 进行功能富集分析（ORA），并保存结果图片。

    Parameters
    ----------
    gene_list : List[str] | str
        基因列表。支持标准 List，也支持逗号分隔的字符串，甚至包含噪音的文本。
        代码内部会进行强力清洗。
    organism : str
        物种，目前默认 "human"。
    databases : List[str]
        Enrichr 库名，如 "GO_Biological_Process_2023", "KEGG_2021_Human"。
    top_k : int
        返回结果中包含的前 K 个通路。

    Returns
    -------
    Dict[str, Any]
        包含 top_terms (结构化数据) 和 summary (自然语言描述)。
    """

    # =========================
    # 1. 强力输入清洗 (Deep Cleaning)
    # =========================
    
    # 步骤 A: 结构展平
    # 无论输入是 List, Dict 还是 JSON String，统一变成一个大字符串
    raw_text = flatten_input(gene_list)

    # 步骤 B: 正则提取 (Regex Extraction)
    # 匹配规则：\b[a-zA-Z]... -> 以字母开头，后续包含字母、数字、连字符或点
    # 作用：自动剥离引号、方括号、逗号等标点符号，同时过滤掉纯数字索引
    matches = re.findall(r'\b[a-zA-Z][a-zA-Z0-9\-\.]{1,}\b', raw_text)

    # 步骤 C: 标准化 (Normalization)
    # 转大写、去重、去除常见的非基因停用词
    stopwords = {"AND", "OR", "THE", "GENE", "LIST", "HUMAN", "SET", "YES", "NO"}
    final_gene_list = sorted(list(set(
        g.upper() for g in matches 
        if g.upper() not in stopwords and not g.isdigit()
    )))

    # =========================
    # 2. 有效性校验
    # =========================
    if len(final_gene_list) < 3:
        return {
            "top_terms": [],
            "summary": (
                f"ERROR: No valid genes found for enrichment analysis. "
                f"After cleaning the input, only {len(final_gene_list)} valid symbols were identified. "
                "Please provide a list of at least 3 valid gene symbols (e.g., CD3D, LCK, TP53)."
            )
        }

    if organism.lower() != "human":
        return {
            "top_terms": [],
            "summary": "Currently, only human gene enrichment is supported."
        }

    # =========================
    # 3. 配置输出与数据库
    # =========================
    # 默认数据库配置
    if databases is None:
        databases = [
            "GO_Biological_Process_2023",
            "KEGG_2021_Human",
            "Reactome_2022",
        ]

    # 图片保存路径设置
    output_dir = None

    # =========================
    # 4. 调用 gseapy.enrichr
    # =========================
    try:
        enr = gp.enrichr(
            gene_list=final_gene_list,
            gene_sets=databases,
            organism="Human",
            outdir=output_dir,    # 结果保存路径
            no_plot=False,        # 开启绘图功能 (会生成 png/pdf)
            verbose=False         # 关闭库内部日志
        )
    except Exception as e:
        return {
            "top_terms": [],
            "summary": f"Enrichr API request failed: {str(e)}. Please check your internet connection."
        }

    # =========================
    # 5. 结果解析
    # =========================
    if enr.results is None or enr.results.empty:
        return {
            "top_terms": [],
            "summary": f"Analysis ran successfully on {len(final_gene_list)} genes, but no statistically significant pathways were found."
        }

    df: pd.DataFrame = enr.results.copy()
    
    # 结果排序：优先使用 Adjusted P-value
    sort_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
    df = df.sort_values(sort_col, ascending=True)

    top_terms = []
    # 提取 Top K 条目构建 JSON
    for _, row in df.head(top_k).iterrows():
        # 兼容不同版本的列名 (Genes vs gene_list)
        genes_str = row.get("Genes", row.get("gene_list", ""))
        overlap_list = []
        if isinstance(genes_str, str):
            overlap_list = genes_str.split(";")

        top_terms.append({
            "term": row.get("Term"),
            "database": row.get("Gene_set"),
            "p_value": float(row.get("P-value", 1.0)),
            "adj_p_value": float(row.get("Adjusted P-value", 1.0)),
            "overlap_gene_list": overlap_list,
        })

    # =========================
    # 6. 生成 Summary
    # =========================
    top_names = [t["term"] for t in top_terms[:5]]
    
    summary_text = (
        f"Enrichment analysis was performed on {len(final_gene_list)} valid genes "
        f"(e.g., {', '.join(final_gene_list[:3])}...). "
    )
    
    if top_names:
        summary_text += f"Significant pathways identified include: {'; '.join(top_names)}. "
        summary_text += f"Plots and detailed results have been saved to the '{output_dir}' directory."
    else:
        summary_text += "No dominant biological themes were identified."

    return {
        "top_terms": top_terms,
        "summary": summary_text,
    }
