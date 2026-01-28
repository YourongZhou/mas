import requests
from typing import Dict

def query_mygene(args: Dict) -> Dict:
    """
    MyGene 工具的具体执行逻辑
    args: 包含 'gene_symbol' 的字典
    """
    # 从通用参数字典中提取特定参数
    gene = args.get("gene_symbol")
    if not gene:
        return {"error": "Missing gene_symbol"}

    print(f"\n[Tool Execution] Querying MyGene for: {gene}")

    url = "https://mygene.info/v3/query"
    params = {
        "q": gene,
        "species": "human",
        "size": 1
    }

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        # 返回第一个匹配结果
        return data["hits"][0] if data.get("hits") else {"error": "No hits found"}
    except Exception as e:
        return {"error": str(e)}

# 定义工具元数据
MYGENE_DEF = {
    "name": "mygene",
    "description": (
        "Retrieve gene information. "
        "Use this tool when the user asks about gene functions, summary, or biology. "
        "Input should contain 'gene_symbol'."
    ),
    "func": query_mygene
}