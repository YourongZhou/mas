# Code segment 1 (from code_dev)
# ============================================================
import scanpy as sc
import pandas as pd
import os

# 设置路径
input_path = "/app/data/pbmc3k.h5ad"
output_path = "/app/output"

# 读取数据
adata = sc.read_hdf5(input_path, cache=False)

# 质控
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# 标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# PCA降维
sc.tl.pca(adata, n_comps=50)

# UMAP可视化
sc.tl.umap(adata)

# Leiden聚类
sc.tl.leiden(adata, resolution=0.5)

# 保存UMAP图
sc.pl.umap(adata, color='leiden', title='Clustering UMAP', show=False)
plt.savefig(os.path.join(output_path, "clustering_umap.png"))

# 差异表达分析
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

# 获取每个聚类的top10差异基因
cluster_results = {}
for cluster in adata.obs['leiden'].unique():
    genes = adata.uns['rank_genes_groups']['names'][cluster]
    cluster_results[cluster] = genes[:10]

# 生成分析摘要
analysis_summary = f"细胞总数：{adata.n_obs}，基因总数：{adata.n_vars}，聚类数量：{len(adata.obs['leiden'].cat.categories)}"

# 输出结果标记
print(f"===RESULT==={analysis_summary}===")

# Code segment 2 (from code_solution)
# ============================================================
{'code': 'import scanpy as sc\nimport pandas as pd\nimport os\n\n# 设置路径\ninput_path = "/app/data/pbmc3k.h5ad"\noutput_path = "/app/output"\n\n# 读取数据\nadata = sc.read_hdf5(input_path, cache=False)\n\n# 质控\nsc.pp.filter_cells(adata, min_genes=200)\nsc.pp.filter_genes(adata, min_cells=3)\n\n# 标准化\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# PCA降维\nsc.tl.pca(adata, n_comps=50)\n\n# UMAP可视化\nsc.tl.umap(adata)\n\n# Leiden聚类\nsc.tl.leiden(adata, resolution=0.5)\n\n# 保存UMAP图\nsc.pl.umap(adata, color=\'leiden\', title=\'Clustering UMAP\', show=False)\nplt.savefig(os.path.join(output_path, "clustering_umap.png"))\n\n# 差异表达分析\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\')\n\n# 获取每个聚类的top10差异基因\ncluster_results = {}\nfor cluster in adata.obs[\'leiden\'].unique():\n    genes = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    cluster_results[cluster] = genes[:10]\n\n# 生成分析摘要\nanalysis_summary = f"细胞总数：{adata.n_obs}，基因总数：{adata.n_vars}，聚类数量：{len(adata.obs[\'leiden\'].cat.categories)}"\n\n# 输出结果标记\nprint(f"===RESULT==={analysis_summary}===")', 'requirements': 'scanpy>=1.9.0\nmatplotlib>=3.4.0\nnumpy>=1.21.0\npandas>=1.3.0\nscipy>=1.7.0\nanndata>=0.8.0\nigraph\nleidenalg', 'task': '', 'error': '未知错误', 'success': False}

# Code segment 3 (from code_solution)
# ============================================================
{'code': 'import scanpy as sc\nimport pandas as pd\nimport os\n\n# 设置路径\ninput_path = "/app/data/pbmc3k.h5ad"\noutput_path = "/app/output"\n\n# 读取数据\nadata = sc.read_hdf5(input_path, cache=False)\n\n# 质控\nsc.pp.filter_cells(adata, min_genes=200)\nsc.pp.filter_genes(adata, min_cells=3)\n\n# 标准化\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# PCA降维\nsc.tl.pca(adata, n_comps=50)\n\n# UMAP可视化\nsc.tl.umap(adata)\n\n# Leiden聚类\nsc.tl.leiden(adata, resolution=0.5)\n\n# 保存UMAP图\nsc.pl.umap(adata, color=\'leiden\', title=\'Clustering UMAP\', show=False)\nplt.savefig(os.path.join(output_path, "clustering_umap.png"))\n\n# 差异表达分析\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\')\n\n# 获取每个聚类的top10差异基因\ncluster_results = {}\nfor cluster in adata.obs[\'leiden\'].unique():\n    genes = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    cluster_results[cluster] = genes[:10]\n\n# 生成分析摘要\nanalysis_summary = f"细胞总数：{adata.n_obs}，基因总数：{adata.n_vars}，聚类数量：{len(adata.obs[\'leiden\'].cat.categories)}"\n\n# 输出结果标记\nprint(f"===RESULT==={analysis_summary}===")', 'requirements': 'scanpy>=1.9.0\nmatplotlib>=3.4.0\nnumpy>=1.21.0\npandas>=1.3.0\nscipy>=1.7.0\nanndata>=0.8.0\nigraph\nleidenalg', 'task': '', 'error': '未知错误', 'success': False}

# Code segment 4 (from code_solution)
# ============================================================
{'code': 'import scanpy as sc\nimport pandas as pd\nimport os\n\n# 设置路径\ninput_path = "/app/data/pbmc3k.h5ad"\noutput_path = "/app/output"\n\n# 读取数据\nadata = sc.read_hdf5(input_path, cache=False)\n\n# 质控\nsc.pp.filter_cells(adata, min_genes=200)\nsc.pp.filter_genes(adata, min_cells=3)\n\n# 标准化\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# PCA降维\nsc.tl.pca(adata, n_comps=50)\n\n# UMAP可视化\nsc.tl.umap(adata)\n\n# Leiden聚类\nsc.tl.leiden(adata, resolution=0.5)\n\n# 保存UMAP图\nsc.pl.umap(adata, color=\'leiden\', title=\'Clustering UMAP\', show=False)\nplt.savefig(os.path.join(output_path, "clustering_umap.png"))\n\n# 差异表达分析\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\')\n\n# 获取每个聚类的top10差异基因\ncluster_results = {}\nfor cluster in adata.obs[\'leiden\'].unique():\n    genes = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    cluster_results[cluster] = genes[:10]\n\n# 生成分析摘要\nanalysis_summary = f"细胞总数：{adata.n_obs}，基因总数：{adata.n_vars}，聚类数量：{len(adata.obs[\'leiden\'].cat.categories)}"\n\n# 输出结果标记\nprint(f"===RESULT==={analysis_summary}===")', 'requirements': 'scanpy>=1.9.0\nmatplotlib>=3.4.0\nnumpy>=1.21.0\npandas>=1.3.0\nscipy>=1.7.0\nanndata>=0.8.0\nigraph\nleidenalg', 'task': '', 'error': '未知错误', 'success': False}

