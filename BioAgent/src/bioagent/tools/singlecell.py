from __future__ import annotations

import json
from pathlib import Path

from bioagent.config import AgentConfig
from bioagent.docker_runner import DockerRunner
from bioagent.logging_utils import RunLogger

from .common import resolve_data_path, to_repo_mount_path


def build_scanpy_pipeline_code(container_data_path: str) -> str:
    payload = json.dumps(container_data_path)
    return f'''
import json
import os
import warnings

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

INPUT_PATH = {payload}
OUTDIR = "/work/outputs/scrnaseq_scanpy_core"
os.makedirs(OUTDIR, exist_ok=True)

summary = {{}}
adata = sc.read_h5ad(INPUT_PATH)
adata.var_names_make_unique()
summary["input_path"] = INPUT_PATH
summary["initial_shape"] = list(adata.shape)
summary["obs_columns"] = list(map(str, adata.obs.columns))
summary["var_columns"] = list(map(str, adata.var.columns))
summary["obsm_keys"] = list(map(str, adata.obsm.keys()))
summary["x_is_sparse"] = bool(sparse.issparse(adata.X))

adata.var["mt"] = adata.var_names.astype(str).str.upper().str.startswith(("MT-", "MT."))
qc_vars = ["mt"] if bool(adata.var["mt"].any()) else None
sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, log1p=False, inplace=True)
qc_cols = [c for c in ["n_genes_by_counts", "total_counts", "pct_counts_mt"] if c in adata.obs]
if qc_cols:
    adata.obs[qc_cols].describe().to_csv(os.path.join(OUTDIR, "qc_metrics_summary.csv"))
    adata.obs[qc_cols].to_csv(os.path.join(OUTDIR, "cell_qc_metrics.csv"))

if qc_cols:
    fig, axes = plt.subplots(1, len(qc_cols), figsize=(4 * len(qc_cols), 3))
    if len(qc_cols) == 1:
        axes = [axes]
    for ax, col in zip(axes, qc_cols):
        values = pd.to_numeric(adata.obs[col], errors="coerce")
        ax.hist(values.dropna(), bins=40, color="#4C78A8", alpha=0.85)
        ax.set_title(col)
    fig.tight_layout()
    fig.savefig(os.path.join(OUTDIR, "qc_histograms.png"), dpi=180)
    plt.close(fig)

if adata.raw is None:
    adata.raw = adata

try:
    max_x = float(adata.X.max())
except Exception:
    max_x = 0.0
summary["max_x_before_norm"] = max_x
if max_x > 30:
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    summary["normalization"] = "normalize_total_1e4_log1p"
else:
    summary["normalization"] = "input_appears_already_normalized_or_log_scaled"

n_top = min(2000, max(1, adata.n_vars - 1))
try:
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat", subset=False)
except Exception as exc:
    summary["hvg_warning"] = str(exc)
    adata.var["highly_variable"] = True
summary["n_hvg"] = int(adata.var.get("highly_variable", pd.Series(False, index=adata.var_names)).sum())

try:
    sc.pp.scale(adata, max_value=10)
except Exception as exc:
    summary["scale_warning"] = str(exc)

n_comps = max(2, min(50, adata.n_obs - 1, adata.n_vars - 1))
try:
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=True)
except Exception:
    sc.tl.pca(adata, n_comps=n_comps, use_highly_variable=False)
summary["n_pcs"] = n_comps

sc.pp.neighbors(adata, n_neighbors=min(15, max(2, adata.n_obs - 1)), n_pcs=min(30, n_comps))
cluster_key = "leiden_0.8"
cluster_status = {{}}
for resolution in [0.4, 0.6, 0.8, 1.0]:
    key = f"leiden_{{resolution}}"
    try:
        sc.tl.leiden(adata, resolution=resolution, key_added=key)
        cluster_status[key] = int(adata.obs[key].nunique())
    except Exception as exc:
        cluster_status[key] = f"failed: {{exc}}"
summary["cluster_status"] = cluster_status

if cluster_key not in adata.obs:
    fallback_cols = [c for c in adata.obs.columns if str(c).startswith("leiden")]
    if fallback_cols:
        cluster_key = str(fallback_cols[0])

if "X_umap" not in adata.obsm:
    sc.tl.umap(adata)

plot_colors = [cluster_key]
for candidate in ["cell type", "cell_type", "CellType", "annotation"]:
    if candidate in adata.obs:
        plot_colors.append(candidate)
        break
for color in plot_colors:
    try:
        sc.pl.umap(adata, color=color, show=False)
        safe = str(color).replace(" ", "_").replace("/", "_")
        plt.title(f"UMAP by {{color}}")
        plt.savefig(os.path.join(OUTDIR, f"umap_{{safe}}.png"), dpi=220, bbox_inches="tight")
        plt.close()
    except Exception as exc:
        summary.setdefault("plot_warnings", []).append(f"{{color}}: {{exc}}")

markers_path = os.path.join(OUTDIR, "cluster_markers.csv")
try:
    if cluster_key in adata.obs and adata.obs[cluster_key].nunique() > 1:
        sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon")
        marker_frames = []
        for group in adata.obs[cluster_key].cat.categories if hasattr(adata.obs[cluster_key], "cat") else sorted(adata.obs[cluster_key].unique()):
            df = sc.get.rank_genes_groups_df(adata, group=group)
            df.insert(0, "cluster", group)
            marker_frames.append(df)
        markers = pd.concat(marker_frames, ignore_index=True)
        markers.to_csv(markers_path, index=False)
        summary["marker_rows"] = int(len(markers))
except Exception as exc:
    summary["marker_warning"] = str(exc)

if "X_umap" in adata.obsm:
    pd.DataFrame(adata.obsm["X_umap"], index=adata.obs_names, columns=["UMAP1", "UMAP2"]).to_csv(os.path.join(OUTDIR, "umap_coordinates.csv"))
if "X_pca" in adata.obsm:
    pcs = adata.obsm["X_pca"]
    pd.DataFrame(pcs, index=adata.obs_names, columns=[f"PC{{i+1}}" for i in range(pcs.shape[1])]).to_csv(os.path.join(OUTDIR, "pca_coordinates.csv"))
adata.obs.to_csv(os.path.join(OUTDIR, "cell_metadata.csv"))
processed_path = os.path.join(OUTDIR, "adata_processed.h5ad")
adata.write_h5ad(processed_path, compression="gzip")

summary["final_shape"] = list(adata.shape)
summary["cluster_key"] = cluster_key
summary["outputs"] = sorted(os.listdir(OUTDIR))
pdf_path = os.path.join(OUTDIR, "scrna_analysis_report.pdf")
with PdfPages(pdf_path) as pdf:
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.08, 0.94, "Single-cell RNA-seq Scanpy Core Analysis", fontsize=16, weight="bold")
    fig.text(0.08, 0.89, f"Input: {{INPUT_PATH}}", fontsize=9)
    fig.text(0.08, 0.86, f"Initial shape: {{summary['initial_shape']}}", fontsize=10)
    fig.text(0.08, 0.83, f"Final shape: {{summary['final_shape']}}", fontsize=10)
    fig.text(0.08, 0.80, f"Cluster key: {{cluster_key}}", fontsize=10)
    fig.text(0.08, 0.77, f"Cluster status: {{cluster_status}}", fontsize=9)
    fig.text(0.08, 0.72, "Outputs:\\n" + "\\n".join(summary["outputs"][:24]), fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)
    for image_name in ["qc_histograms.png", f"umap_{{cluster_key.replace(' ', '_').replace('/', '_')}}.png"]:
        path = os.path.join(OUTDIR, image_name)
        if os.path.exists(path):
            img = plt.imread(path)
            fig, ax = plt.subplots(figsize=(8.5, 7))
            ax.imshow(img)
            ax.axis("off")
            ax.set_title(image_name)
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

summary["outputs"] = sorted(os.listdir(OUTDIR))
with open(os.path.join(OUTDIR, "analysis_summary.json"), "w", encoding="utf-8") as handle:
    json.dump(summary, handle, ensure_ascii=False, indent=2)
with open(os.path.join(OUTDIR, "analysis_summary.txt"), "w", encoding="utf-8") as handle:
    handle.write(json.dumps(summary, ensure_ascii=False, indent=2))

print("===RESULT===" + json.dumps(summary, ensure_ascii=False) + "===")
'''


def run_scanpy_singlecell_pipeline_impl(
    config: AgentConfig,
    logger: RunLogger,
    run_dir: Path,
    data_path: str,
    env_profile: str = "py-scverse-v1",
    skill_id: str = "scrnaseq-scanpy-core-analysis",
    timeout_s: int = 1800,
) -> dict:
    host_data_path = resolve_data_path(config, run_dir, data_path)
    if host_data_path.suffix.lower() != ".h5ad":
        return {"ok": False, "error": f"Expected .h5ad input, got: {host_data_path}"}
    try:
        container_path = to_repo_mount_path(config, host_data_path)
    except ValueError as exc:
        return {
            "ok": False,
            "error": f"Data path must be inside repo root so Docker can see it via /repo: {exc}",
            "host_data_path": str(host_data_path),
        }
    code = build_scanpy_pipeline_code(container_path)
    runner = DockerRunner(config=config, logger=logger, run_dir=run_dir)
    result = runner.execute_python(code, env_profile=env_profile, requirements="", timeout_s=timeout_s).to_dict()
    result["skill_id"] = skill_id
    result["host_data_path"] = str(host_data_path)
    result["expected_output_dir"] = str(run_dir / "outputs" / "scrnaseq_scanpy_core")
    return result

