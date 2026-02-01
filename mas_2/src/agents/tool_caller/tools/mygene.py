from langchain_core.tools import tool
import requests
from typing import Dict, Union

@tool
def query_mygene(gene_symbol: str) -> Union[Dict, str]:
    """
    Retrieve gene information (function, summary, genomic location).
    Use this tool when the user asks about specific gene biology.
    
    Args:
        gene_symbol: The standard symbol of the gene (e.g., "TP53", "CD4").
    """
    if not gene_symbol:
        return {"error": "Missing gene_symbol"}

    print(f"\n[Tool Execution] Querying MyGene for: {gene_symbol}")

    url = "https://mygene.info/v3/query"
    params = {
        "q": gene_symbol,  
        "species": "human",
        "size": 1
    }

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        return data["hits"][0] if data.get("hits") else {"error": "No hits found"}
    except Exception as e:
        return {"error": str(e)}

