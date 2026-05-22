from __future__ import annotations

import re
from typing import Any

import requests


def _clean_gene_list(gene_list: list[str] | str | Any) -> list[str]:
    if isinstance(gene_list, str):
        text = gene_list
    elif isinstance(gene_list, (list, tuple, set)):
        text = " ".join(str(item) for item in gene_list)
    else:
        text = str(gene_list)
    genes = re.findall(r"\b[A-Za-z][A-Za-z0-9_.-]{1,}\b", text)
    stop = {"GENE", "GENES", "LIST", "AND", "OR", "THE"}
    out: list[str] = []
    seen = set()
    for gene in genes:
        symbol = gene.upper()
        if symbol in stop or symbol in seen:
            continue
        seen.add(symbol)
        out.append(symbol)
    return out


def query_mygene_impl(gene_symbol: str) -> dict[str, Any]:
    symbol = str(gene_symbol or "").strip()
    if not symbol:
        return {"error": "Missing gene_symbol"}
    try:
        response = requests.get(
            "https://mygene.info/v3/query",
            params={"q": symbol, "species": "human", "size": 1},
            timeout=30,
        )
        response.raise_for_status()
        data = response.json()
        hits = data.get("hits") or []
        return hits[0] if hits else {"error": "No hits found", "gene_symbol": symbol}
    except Exception as exc:
        return {"error": str(exc), "gene_symbol": symbol}


def gene_set_enrichment_impl(
    gene_list: list[str] | str,
    organism: str = "human",
    databases: list[str] | None = None,
    top_k: int = 10,
) -> dict[str, Any]:
    genes = _clean_gene_list(gene_list)
    if len(genes) < 3:
        return {"error": "Need at least 3 valid gene symbols.", "genes": genes}
    if organism.lower() != "human":
        return {"error": "Only human enrichment is currently supported.", "organism": organism}
    databases = databases or ["GO_Biological_Process_2023", "KEGG_2021_Human", "Reactome_2022"]
    try:
        import gseapy as gp

        enr = gp.enrichr(gene_list=genes, gene_sets=databases, organism="Human", outdir=None, no_plot=True, verbose=False)
        if enr.results is None or enr.results.empty:
            return {"genes": genes, "top_terms": [], "summary": "No significant terms returned."}
        df = enr.results.copy()
        sort_col = "Adjusted P-value" if "Adjusted P-value" in df.columns else "P-value"
        df = df.sort_values(sort_col, ascending=True).head(top_k)
        return {
            "genes": genes,
            "databases": databases,
            "top_terms": df.to_dict(orient="records"),
        }
    except Exception as exc:
        return {"error": str(exc), "genes": genes, "databases": databases}


def run_celltype_annotation_impl(gene_list: list[str] | str) -> dict[str, Any]:
    genes = set(_clean_gene_list(gene_list))
    marker_sets = {
        "B cell": {"MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA"},
        "Plasma cell": {"MZB1", "XBP1", "JCHAIN", "IGHG1", "SDC1"},
        "T cell": {"CD3D", "CD3E", "CD2", "TRAC", "IL7R"},
        "NK cell": {"NKG7", "GNLY", "PRF1", "KLRD1", "GZMB"},
        "Monocyte": {"LYZ", "S100A8", "S100A9", "FCN1", "CTSS"},
        "Dendritic cell": {"FCER1A", "CST3", "CLEC10A", "LILRA4"},
        "Erythroid": {"HBB", "HBA1", "HBA2", "ALAS2"},
        "Platelet": {"PPBP", "PF4", "NRGN"},
    }
    scores = []
    for label, markers in marker_sets.items():
        overlap = sorted(genes & markers)
        scores.append({"label": label, "score": len(overlap), "overlap": overlap})
    scores.sort(key=lambda item: item["score"], reverse=True)
    best = scores[0] if scores else {"label": "unknown", "score": 0, "overlap": []}
    return {
        "final_prediction": best["label"] if best["score"] > 0 else "unknown",
        "confidence": min(1.0, best["score"] / 4.0),
        "details": scores,
        "input_genes": sorted(genes),
    }

