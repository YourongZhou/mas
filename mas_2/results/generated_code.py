# Code segment 1 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# 基础质控：计算线粒体基因比例（假设mt基因以'MT-'或'mt-'开头）
adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 过滤低质量细胞和基因（典型阈值，可根据实际调整）
adata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限
adata = adata[adata.obs.pct_counts_mt < 15, :]        # 线粒体比例上限
adata = adata[adata.obs.total_counts > 500, :]       # 总UMI数下限

# 过滤低表达基因（至少在3个细胞中表达）
sc.pp.filter_genes(adata, min_cells=3)

# 标准化 & 对数变换
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 2. PCA降维
print("Step 2: Performing PCA...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# 4. Leiden聚类
print("Step 4: Clustering with Leiden...")
sc.tl.leiden(adata, key_added='leiden', resolution=1.0)  # 可调resolution控制簇数

# 5. 差异表达分析：每个簇 vs 所有其他簇（'t-test' method for speed & robustness）
print("Step 5: Differential expression analysis per cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', use_raw=False, n_genes=1000)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes per cluster...")
clusters = adata.obs['leiden'].unique().tolist()
deg_dict = {}

for cluster in clusters:
    # 获取该cluster的top 10基因名（rank_genes_groups结果按logfoldchange排序，但scanpy默认按'names'字段返回）
    # 使用rank_genes_groups结果中的'names'字段（按统计显著性排序）
    gene_names = adata.uns['rank_genes_groups']['names'][cluster][:10]
    deg_dict[cluster] = gene_names.tolist()

# 保存结果为CSV（每个簇一行，top 10基因用逗号分隔）
deg_df = pd.DataFrame.from_dict(deg_dict, orient='index', columns=[f'top_{i+1}' for i in range(10)])
deg_csv_path = os.path.join(result_dir, "leiden_top10_de_genes.csv")
deg_df.to_csv(deg_csv_path)
print(f"✅ Top 10 DE genes per cluster saved to: {deg_csv_path}")

# 可选：保存带注释的AnnData对象（含UMAP、leiden、DE结果）
adata.write(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))
print(f"✅ Processed AnnData saved to: {os.path.join(result_dir, 'pbmc3k_processed_with_de.h5ad')}")

# 返回每个聚类的差异表达基因列表（字典格式）
deg_dict

# Code segment 2 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础质控：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（典型阈值，可根据实际调整）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]       # 总UMI数下限\n\n# 过滤低表达基因（至少在3个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=3)\n\n# 标准化 & 对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)  # 可调resolution控制簇数\n\n# 5. 差异表达分析：每个簇 vs 所有其他簇（\'t-test\' method for speed & robustness）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'t-test\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\ndeg_dict = {}\n\nfor cluster in clusters:\n    # 获取该cluster的top 10基因名（rank_genes_groups结果按logfoldchange排序，但scanpy默认按\'names\'字段返回）\n    # 使用rank_genes_groups结果中的\'names\'字段（按统计显著性排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    deg_dict[cluster] = gene_names.tolist()\n\n# 保存结果为CSV（每个簇一行，top 10基因用逗号分隔）\ndeg_df = pd.DataFrame.from_dict(deg_dict, orient=\'index\', columns=[f\'top_{i+1}\' for i in range(10)])\ndeg_csv_path = os.path.join(result_dir, "leiden_top10_de_genes.csv")\ndeg_df.to_csv(deg_csv_path)\nprint(f"✅ Top 10 DE genes per cluster saved to: {deg_csv_path}")\n\n# 可选：保存带注释的AnnData对象（含UMAP、leiden、DE结果）\nadata.write(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\nprint(f"✅ Processed AnnData saved to: {os.path.join(result_dir, \'pbmc3k_processed_with_de.h5ad\')}")\n\n# 返回每个聚类的差异表达基因列表（字典格式）\ndeg_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 3 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础质控：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（典型阈值，可根据实际调整）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]       # 总UMI数下限\n\n# 过滤低表达基因（至少在3个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=3)\n\n# 标准化 & 对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)  # 可调resolution控制簇数\n\n# 5. 差异表达分析：每个簇 vs 所有其他簇（\'t-test\' method for speed & robustness）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'t-test\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\ndeg_dict = {}\n\nfor cluster in clusters:\n    # 获取该cluster的top 10基因名（rank_genes_groups结果按logfoldchange排序，但scanpy默认按\'names\'字段返回）\n    # 使用rank_genes_groups结果中的\'names\'字段（按统计显著性排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    deg_dict[cluster] = gene_names.tolist()\n\n# 保存结果为CSV（每个簇一行，top 10基因用逗号分隔）\ndeg_df = pd.DataFrame.from_dict(deg_dict, orient=\'index\', columns=[f\'top_{i+1}\' for i in range(10)])\ndeg_csv_path = os.path.join(result_dir, "leiden_top10_de_genes.csv")\ndeg_df.to_csv(deg_csv_path)\nprint(f"✅ Top 10 DE genes per cluster saved to: {deg_csv_path}")\n\n# 可选：保存带注释的AnnData对象（含UMAP、leiden、DE结果）\nadata.write(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\nprint(f"✅ Processed AnnData saved to: {os.path.join(result_dir, \'pbmc3k_processed_with_de.h5ad\')}")\n\n# 返回每个聚类的差异表达基因列表（字典格式）\ndeg_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 4 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# 基础质控：计算线粒体基因比例（假设mt基因以'MT-'或'mt-'开头）
adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 过滤低质量细胞和基因（典型阈值，可根据实际调整）
adata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限
adata = adata[adata.obs.pct_counts_mt < 15, :]        # 线粒体比例上限
adata = adata[adata.obs.total_counts > 500, :]       # 总UMI数下限

# 过滤低表达基因（至少在10个细胞中表达）
sc.pp.filter_genes(adata, min_cells=10)

# 标准化 + 对数变换
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 2. PCA 降维
print("Step 2: Performing PCA...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP 可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# 4. Leiden 聚类
print("Step 4: Clustering with Leiden...")
sc.tl.leiden(adata, key_added='leiden', resolution=1.0)  # 可调 resolution 控制聚类粒度

# 5. 差异表达分析（每个 cluster vs 所有其他 cells）
print("Step 5: Performing differential expression analysis per cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=adata.n_vars, use_raw=False)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes per cluster...")
clusters = adata.obs['leiden'].unique().tolist()
de_genes_per_cluster = {}

for cluster in clusters:
    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，兼顾统计显著性）
    # rank_genes_groups 输出为 pandas DataFrame，每列一个 group
    try:
        # 提取 'names' 列中对应 cluster 的 top 10 基因名
        gene_names = adata.uns['rank_genes_groups']['names'][cluster][:10]
        # 转为 Python list（避免 numpy array）
        de_genes_per_cluster[cluster] = gene_names.tolist() if hasattr(gene_names, 'tolist') else list(gene_names)
    except (KeyError, IndexError):
        de_genes_per_cluster[cluster] = []

# 保存结果到 CSV 文件（每个 cluster 一行，基因用逗号分隔）
result_df = pd.DataFrame([
    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} 
    for cl, genes in de_genes_per_cluster.items()
])
result_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)

# 同时保存为更易读的 JSON 格式（可选）
import json
with open(os.path.join(result_dir, "leiden_top10_de_genes.json"), "w") as f:
    json.dump(de_genes_per_cluster, f, indent=2)

print("✅ Analysis completed.")
print("📁 Results saved to:", result_dir)
print("📋 Top 10 DE genes per cluster:")
for cl, genes in de_genes_per_cluster.items():
    print(f"  Cluster {cl}: {genes}")

# 返回字典形式的结果（符合任务要求）
de_genes_per_cluster

# Code segment 5 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础质控：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（典型阈值，可根据实际调整）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]       # 总UMI数下限\n\n# 过滤低表达基因（至少在3个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=3)\n\n# 标准化 & 对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)  # 可调resolution控制簇数\n\n# 5. 差异表达分析：每个簇 vs 所有其他簇（\'t-test\' method for speed & robustness）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'t-test\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\ndeg_dict = {}\n\nfor cluster in clusters:\n    # 获取该cluster的top 10基因名（rank_genes_groups结果按logfoldchange排序，但scanpy默认按\'names\'字段返回）\n    # 使用rank_genes_groups结果中的\'names\'字段（按统计显著性排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    deg_dict[cluster] = gene_names.tolist()\n\n# 保存结果为CSV（每个簇一行，top 10基因用逗号分隔）\ndeg_df = pd.DataFrame.from_dict(deg_dict, orient=\'index\', columns=[f\'top_{i+1}\' for i in range(10)])\ndeg_csv_path = os.path.join(result_dir, "leiden_top10_de_genes.csv")\ndeg_df.to_csv(deg_csv_path)\nprint(f"✅ Top 10 DE genes per cluster saved to: {deg_csv_path}")\n\n# 可选：保存带注释的AnnData对象（含UMAP、leiden、DE结果）\nadata.write(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\nprint(f"✅ Processed AnnData saved to: {os.path.join(result_dir, \'pbmc3k_processed_with_de.h5ad\')}")\n\n# 返回每个聚类的差异表达基因列表（字典格式）\ndeg_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 6 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础质控：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（典型阈值，可根据实际调整）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]       # 总UMI数下限\n\n# 过滤低表达基因（至少在10个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=10)\n\n# 标准化 + 对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)  # 可调 resolution 控制聚类粒度\n\n# 5. 差异表达分析（每个 cluster vs 所有其他 cells）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_per_cluster = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，兼顾统计显著性）\n    # rank_genes_groups 输出为 pandas DataFrame，每列一个 group\n    try:\n        # 提取 \'names\' 列中对应 cluster 的 top 10 基因名\n        gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n        # 转为 Python list（避免 numpy array）\n        de_genes_per_cluster[cluster] = gene_names.tolist() if hasattr(gene_names, \'tolist\') else list(gene_names)\n    except (KeyError, IndexError):\n        de_genes_per_cluster[cluster] = []\n\n# 保存结果到 CSV 文件（每个 cluster 一行，基因用逗号分隔）\nresult_df = pd.DataFrame([\n    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} \n    for cl, genes in de_genes_per_cluster.items()\n])\nresult_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 同时保存为更易读的 JSON 格式（可选）\nimport json\nwith open(os.path.join(result_dir, "leiden_top10_de_genes.json"), "w") as f:\n    json.dump(de_genes_per_cluster, f, indent=2)\n\nprint("✅ Analysis completed.")\nprint("📁 Results saved to:", result_dir)\nprint("📋 Top 10 DE genes per cluster:")\nfor cl, genes in de_genes_per_cluster.items():\n    print(f"  Cluster {cl}: {genes}")\n\n# 返回字典形式的结果（符合任务要求）\nde_genes_per_cluster', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 7 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础质控：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（典型阈值，可根据实际调整）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]       # 总UMI数下限\n\n# 过滤低表达基因（至少在10个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=10)\n\n# 标准化 + 对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)  # 可调 resolution 控制聚类粒度\n\n# 5. 差异表达分析（每个 cluster vs 所有其他 cells）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_per_cluster = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，兼顾统计显著性）\n    # rank_genes_groups 输出为 pandas DataFrame，每列一个 group\n    try:\n        # 提取 \'names\' 列中对应 cluster 的 top 10 基因名\n        gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n        # 转为 Python list（避免 numpy array）\n        de_genes_per_cluster[cluster] = gene_names.tolist() if hasattr(gene_names, \'tolist\') else list(gene_names)\n    except (KeyError, IndexError):\n        de_genes_per_cluster[cluster] = []\n\n# 保存结果到 CSV 文件（每个 cluster 一行，基因用逗号分隔）\nresult_df = pd.DataFrame([\n    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} \n    for cl, genes in de_genes_per_cluster.items()\n])\nresult_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 同时保存为更易读的 JSON 格式（可选）\nimport json\nwith open(os.path.join(result_dir, "leiden_top10_de_genes.json"), "w") as f:\n    json.dump(de_genes_per_cluster, f, indent=2)\n\nprint("✅ Analysis completed.")\nprint("📁 Results saved to:", result_dir)\nprint("📋 Top 10 DE genes per cluster:")\nfor cl, genes in de_genes_per_cluster.items():\n    print(f"  Cluster {cl}: {genes}")\n\n# 返回字典形式的结果（符合任务要求）\nde_genes_per_cluster', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 8 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# 基础 QC：计算线粒体基因比例（假设mt基因以'MT-'或'mt-'开头）
adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 过滤低质量细胞和基因（可根据实际需求调整阈值）
# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%
adata = adata[adata.obs.n_genes_by_counts >= 500]
adata = adata[adata.obs.pct_counts_mt < 20]

# 过滤低表达基因（在至少 10 个细胞中表达）
sc.pp.filter_genes(adata, min_cells=10)

# 标准化与对数变换
sc.pp.normalize_total(adata, target_sum=1e4)  # 归一化到每细胞10000 counts
sc.pp.log1p(adata)  # log(1+x)

# 2. PCA 降维
print("Step 2: Performing PCA...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP 可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# 4. Leiden 聚类
print("Step 4: Clustering with Leiden...")
sc.tl.leiden(adata, key_added='leiden', resolution=1.0)

# 5. 差异表达分析（每个 cluster vs rest）
print("Step 5: Performing differential expression analysis per cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=adata.n_vars, use_raw=False)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes for each cluster...")

# 初始化结果字典
de_genes_dict = {}

# 获取 rank_genes_groups 结果（返回的是 AnnData.uns['rank_genes_groups'] 中的 structured arrays）
groups = adata.uns['rank_genes_groups']['names'].dtype.names
for group in groups:
    # 获取该 cluster 的 top 10 基因名（按 wilcoxon score 排序）
    gene_names = adata.uns['rank_genes_groups']['names'][group][:10]
    # 可选：同时获取 logfoldchanges 和 pvals
    logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges'][group][:10]
    pvals = adata.uns['rank_genes_groups']['pvals_adj'][group][:10]
    
    # 构建 DataFrame 并存入字典
    de_df = pd.DataFrame({
        'gene': gene_names,
        'logfoldchange': logfoldchanges,
        'pval_adj': pvals
    })
    de_genes_dict[group] = de_df

# 7. 保存结果
print("Saving results...")
# 保存 UMAP + Leiden 可视化图
sc.pl.umap(adata, color='leiden', save='_leiden_umap.png', show=False)
sc.pl.umap(adata, color='leiden', legend_loc='on data', save='_leiden_umap_ondat.png', show=False)

# 保存每个 cluster 的 top10 DE 基因为 CSV 文件
for cluster_id, df in de_genes_dict.items():
    csv_path = os.path.join(result_dir, f"cluster_{cluster_id}_top10_de_genes.csv")
    df.to_csv(csv_path, index=False)

# 同时汇总为一个 Excel 文件（可选）
summary_excel = os.path.join(result_dir, "all_clusters_top10_de_genes.xlsx")
with pd.ExcelWriter(summary_excel) as writer:
    for cluster_id, df in de_genes_dict.items():
        df.to_excel(writer, sheet_name=f"Cluster_{cluster_id}", index=False)

# 打印各 cluster 的 top10 基因（简要输出）
print("\nTop 10 DE genes per cluster:")
for cluster_id, df in de_genes_dict.items():
    print(f"\nCluster {cluster_id}:")
    print(df[['gene', 'logfoldchange', 'pval_adj']].to_string(index=False))

# 返回每个聚类的差异表达基因列表（仅 gene 名称列表，便于后续使用）
result_lists = {cluster_id: df['gene'].tolist() for cluster_id, df in de_genes_dict.items()}
print(f"\n✅ Analysis completed. Results saved to '{result_dir}'.")
print(f"✅ Returned DE gene lists (dict): keys = cluster IDs, values = list of top 10 gene names.")

# Code segment 9 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础质控：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（典型阈值，可根据实际调整）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]       # 总UMI数下限\n\n# 过滤低表达基因（至少在10个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=10)\n\n# 标准化 + 对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)  # 可调 resolution 控制聚类粒度\n\n# 5. 差异表达分析（每个 cluster vs 所有其他 cells）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_per_cluster = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，兼顾统计显著性）\n    # rank_genes_groups 输出为 pandas DataFrame，每列一个 group\n    try:\n        # 提取 \'names\' 列中对应 cluster 的 top 10 基因名\n        gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n        # 转为 Python list（避免 numpy array）\n        de_genes_per_cluster[cluster] = gene_names.tolist() if hasattr(gene_names, \'tolist\') else list(gene_names)\n    except (KeyError, IndexError):\n        de_genes_per_cluster[cluster] = []\n\n# 保存结果到 CSV 文件（每个 cluster 一行，基因用逗号分隔）\nresult_df = pd.DataFrame([\n    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} \n    for cl, genes in de_genes_per_cluster.items()\n])\nresult_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 同时保存为更易读的 JSON 格式（可选）\nimport json\nwith open(os.path.join(result_dir, "leiden_top10_de_genes.json"), "w") as f:\n    json.dump(de_genes_per_cluster, f, indent=2)\n\nprint("✅ Analysis completed.")\nprint("📁 Results saved to:", result_dir)\nprint("📋 Top 10 DE genes per cluster:")\nfor cl, genes in de_genes_per_cluster.items():\n    print(f"  Cluster {cl}: {genes}")\n\n# 返回字典形式的结果（符合任务要求）\nde_genes_per_cluster', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 10 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础 QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际需求调整阈值）\n# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%\nadata = adata[adata.obs.n_genes_by_counts >= 500]\nadata = adata[adata.obs.pct_counts_mt < 20]\n\n# 过滤低表达基因（在至少 10 个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=10)\n\n# 标准化与对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)  # 归一化到每细胞10000 counts\nsc.pp.log1p(adata)  # log(1+x)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)\n\n# 5. 差异表达分析（每个 cluster vs rest）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes for each cluster...")\n\n# 初始化结果字典\nde_genes_dict = {}\n\n# 获取 rank_genes_groups 结果（返回的是 AnnData.uns[\'rank_genes_groups\'] 中的 structured arrays）\ngroups = adata.uns[\'rank_genes_groups\'][\'names\'].dtype.names\nfor group in groups:\n    # 获取该 cluster 的 top 10 基因名（按 wilcoxon score 排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][group][:10]\n    # 可选：同时获取 logfoldchanges 和 pvals\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][group][:10]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][group][:10]\n    \n    # 构建 DataFrame 并存入字典\n    de_df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval_adj\': pvals\n    })\n    de_genes_dict[group] = de_df\n\n# 7. 保存结果\nprint("Saving results...")\n# 保存 UMAP + Leiden 可视化图\nsc.pl.umap(adata, color=\'leiden\', save=\'_leiden_umap.png\', show=False)\nsc.pl.umap(adata, color=\'leiden\', legend_loc=\'on data\', save=\'_leiden_umap_ondat.png\', show=False)\n\n# 保存每个 cluster 的 top10 DE 基因为 CSV 文件\nfor cluster_id, df in de_genes_dict.items():\n    csv_path = os.path.join(result_dir, f"cluster_{cluster_id}_top10_de_genes.csv")\n    df.to_csv(csv_path, index=False)\n\n# 同时汇总为一个 Excel 文件（可选）\nsummary_excel = os.path.join(result_dir, "all_clusters_top10_de_genes.xlsx")\nwith pd.ExcelWriter(summary_excel) as writer:\n    for cluster_id, df in de_genes_dict.items():\n        df.to_excel(writer, sheet_name=f"Cluster_{cluster_id}", index=False)\n\n# 打印各 cluster 的 top10 基因（简要输出）\nprint("\\nTop 10 DE genes per cluster:")\nfor cluster_id, df in de_genes_dict.items():\n    print(f"\\nCluster {cluster_id}:")\n    print(df[[\'gene\', \'logfoldchange\', \'pval_adj\']].to_string(index=False))\n\n# 返回每个聚类的差异表达基因列表（仅 gene 名称列表，便于后续使用）\nresult_lists = {cluster_id: df[\'gene\'].tolist() for cluster_id, df in de_genes_dict.items()}\nprint(f"\\n✅ Analysis completed. Results saved to \'{result_dir}\'.")\nprint(f"✅ Returned DE gene lists (dict): keys = cluster IDs, values = list of top 10 gene names.")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 11 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础 QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际需求调整阈值）\n# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%\nadata = adata[adata.obs.n_genes_by_counts >= 500]\nadata = adata[adata.obs.pct_counts_mt < 20]\n\n# 过滤低表达基因（在至少 10 个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=10)\n\n# 标准化与对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)  # 归一化到每细胞10000 counts\nsc.pp.log1p(adata)  # log(1+x)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)\n\n# 5. 差异表达分析（每个 cluster vs rest）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes for each cluster...")\n\n# 初始化结果字典\nde_genes_dict = {}\n\n# 获取 rank_genes_groups 结果（返回的是 AnnData.uns[\'rank_genes_groups\'] 中的 structured arrays）\ngroups = adata.uns[\'rank_genes_groups\'][\'names\'].dtype.names\nfor group in groups:\n    # 获取该 cluster 的 top 10 基因名（按 wilcoxon score 排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][group][:10]\n    # 可选：同时获取 logfoldchanges 和 pvals\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][group][:10]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][group][:10]\n    \n    # 构建 DataFrame 并存入字典\n    de_df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval_adj\': pvals\n    })\n    de_genes_dict[group] = de_df\n\n# 7. 保存结果\nprint("Saving results...")\n# 保存 UMAP + Leiden 可视化图\nsc.pl.umap(adata, color=\'leiden\', save=\'_leiden_umap.png\', show=False)\nsc.pl.umap(adata, color=\'leiden\', legend_loc=\'on data\', save=\'_leiden_umap_ondat.png\', show=False)\n\n# 保存每个 cluster 的 top10 DE 基因为 CSV 文件\nfor cluster_id, df in de_genes_dict.items():\n    csv_path = os.path.join(result_dir, f"cluster_{cluster_id}_top10_de_genes.csv")\n    df.to_csv(csv_path, index=False)\n\n# 同时汇总为一个 Excel 文件（可选）\nsummary_excel = os.path.join(result_dir, "all_clusters_top10_de_genes.xlsx")\nwith pd.ExcelWriter(summary_excel) as writer:\n    for cluster_id, df in de_genes_dict.items():\n        df.to_excel(writer, sheet_name=f"Cluster_{cluster_id}", index=False)\n\n# 打印各 cluster 的 top10 基因（简要输出）\nprint("\\nTop 10 DE genes per cluster:")\nfor cluster_id, df in de_genes_dict.items():\n    print(f"\\nCluster {cluster_id}:")\n    print(df[[\'gene\', \'logfoldchange\', \'pval_adj\']].to_string(index=False))\n\n# 返回每个聚类的差异表达基因列表（仅 gene 名称列表，便于后续使用）\nresult_lists = {cluster_id: df[\'gene\'].tolist() for cluster_id, df in de_genes_dict.items()}\nprint(f"\\n✅ Analysis completed. Results saved to \'{result_dir}\'.")\nprint(f"✅ Returned DE gene lists (dict): keys = cluster IDs, values = list of top 10 gene names.")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 12 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# QC：计算线粒体基因比例（假设为 human，mt 基因前缀为 'MT-')
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 过滤低质量细胞和基因（可根据实际数据调整阈值）
adata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 去除高基因数异常细胞（如doublets）
adata = adata[adata.obs.pct_counts_mt < 15, :]        # 去除线粒体比例过高细胞（凋亡/低质量）
adata = adata[adata.obs.n_genes_by_counts > 500, :]   # 去除基因数过少的细胞
adata = adata[:, adata.var.n_cells_by_counts >= 10]   # 保留至少在10个细胞中表达的基因

# 标准化与对数转换
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 2. PCA 降维
print("Step 2: Performing PCA...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP 可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# 4. Leiden 聚类
print("Step 4: Clustering with Leiden...")
sc.tl.leiden(adata, resolution=1.0, key_added='leiden')

# 5. 差异表达分析（每个 cluster vs rest）
print("Step 5: Differential expression analysis per cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=adata.n_vars, use_raw=False)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes per cluster...")
clusters = adata.obs['leiden'].unique().tolist()
de_genes_dict = {}

for cluster in clusters:
    # 获取该 cluster 的 top genes（wilcoxon 检验的 logfoldchanges 排序）
    # 注意：rank_genes_groups 输出结构为字典，含 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj'
    gene_names = adata.uns['rank_genes_groups']['names'][cluster]
    logfc = adata.uns['rank_genes_groups']['logfoldchanges'][cluster]
    pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][cluster]
    
    # 构建 DataFrame 并按 logfoldchange 降序取 top 10（确保有足够基因）
    top10_df = pd.DataFrame({
        'gene': gene_names[:10],
        'logfoldchange': logfc[:10],
        'pval_adj': pvals_adj[:10]
    }).sort_values('logfoldchange', ascending=False).head(10)
    
    de_genes_dict[cluster] = top10_df['gene'].tolist()

# 保存结果到 CSV 文件（每个 cluster 一行，top10 基因用逗号分隔）
summary_df = pd.DataFrame([
    {'cluster': cl, 'top10_de_genes': ', '.join(genes)} 
    for cl, genes in de_genes_dict.items()
])
summary_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes_summary.csv"), index=False)

# 可选：保存更详细的 DE 表（每个 cluster 单独文件）
for cluster, genes in de_genes_dict.items():
    # 重新提取完整 top 10 行详细信息
    idx = list(clusters).index(cluster)
    detailed_df = pd.DataFrame({
        'gene': adata.uns['rank_genes_groups']['names'][idx][:10],
        'logfoldchange': adata.uns['rank_genes_groups']['logfoldchanges'][idx][:10],
        'pval_adj': adata.uns['rank_genes_groups']['pvals_adj'][idx][:10],
        'score': adata.uns['rank_genes_groups']['scores'][idx][:10]
    }).sort_values('logfoldchange', ascending=False)
    detailed_df.to_csv(
        os.path.join(result_dir, f"cluster_{cluster}_top10_de_genes_detailed.csv"),
        index=False
    )

# 打印结果概览
print("\n✅ Analysis completed.")
print(f"Found {len(clusters)} clusters: {clusters}")
for cl, genes in de_genes_dict.items():
    print(f"Cluster {cl}: {genes}")

# 可选：保存带注释的 AnnData 用于后续探索
adata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))

# 返回每个聚类的差异表达基因列表（按要求返回字典格式）
de_genes_dict

# Code segment 13 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# 基础 QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际需求调整阈值）\n# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%\nadata = adata[adata.obs.n_genes_by_counts >= 500]\nadata = adata[adata.obs.pct_counts_mt < 20]\n\n# 过滤低表达基因（在至少 10 个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=10)\n\n# 标准化与对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)  # 归一化到每细胞10000 counts\nsc.pp.log1p(adata)  # log(1+x)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, key_added=\'leiden\', resolution=1.0)\n\n# 5. 差异表达分析（每个 cluster vs rest）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes for each cluster...")\n\n# 初始化结果字典\nde_genes_dict = {}\n\n# 获取 rank_genes_groups 结果（返回的是 AnnData.uns[\'rank_genes_groups\'] 中的 structured arrays）\ngroups = adata.uns[\'rank_genes_groups\'][\'names\'].dtype.names\nfor group in groups:\n    # 获取该 cluster 的 top 10 基因名（按 wilcoxon score 排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][group][:10]\n    # 可选：同时获取 logfoldchanges 和 pvals\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][group][:10]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][group][:10]\n    \n    # 构建 DataFrame 并存入字典\n    de_df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval_adj\': pvals\n    })\n    de_genes_dict[group] = de_df\n\n# 7. 保存结果\nprint("Saving results...")\n# 保存 UMAP + Leiden 可视化图\nsc.pl.umap(adata, color=\'leiden\', save=\'_leiden_umap.png\', show=False)\nsc.pl.umap(adata, color=\'leiden\', legend_loc=\'on data\', save=\'_leiden_umap_ondat.png\', show=False)\n\n# 保存每个 cluster 的 top10 DE 基因为 CSV 文件\nfor cluster_id, df in de_genes_dict.items():\n    csv_path = os.path.join(result_dir, f"cluster_{cluster_id}_top10_de_genes.csv")\n    df.to_csv(csv_path, index=False)\n\n# 同时汇总为一个 Excel 文件（可选）\nsummary_excel = os.path.join(result_dir, "all_clusters_top10_de_genes.xlsx")\nwith pd.ExcelWriter(summary_excel) as writer:\n    for cluster_id, df in de_genes_dict.items():\n        df.to_excel(writer, sheet_name=f"Cluster_{cluster_id}", index=False)\n\n# 打印各 cluster 的 top10 基因（简要输出）\nprint("\\nTop 10 DE genes per cluster:")\nfor cluster_id, df in de_genes_dict.items():\n    print(f"\\nCluster {cluster_id}:")\n    print(df[[\'gene\', \'logfoldchange\', \'pval_adj\']].to_string(index=False))\n\n# 返回每个聚类的差异表达基因列表（仅 gene 名称列表，便于后续使用）\nresult_lists = {cluster_id: df[\'gene\'].tolist() for cluster_id, df in de_genes_dict.items()}\nprint(f"\\n✅ Analysis completed. Results saved to \'{result_dir}\'.")\nprint(f"✅ Returned DE gene lists (dict): keys = cluster IDs, values = list of top 10 gene names.")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 14 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设为 human，mt 基因前缀为 \'MT-\')\nadata.var["mt"] = adata.var_names.str.startswith("MT-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 去除高基因数异常细胞（如doublets）\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 去除线粒体比例过高细胞（凋亡/低质量）\nadata = adata[adata.obs.n_genes_by_counts > 500, :]   # 去除基因数过少的细胞\nadata = adata[:, adata.var.n_cells_by_counts >= 10]   # 保留至少在10个细胞中表达的基因\n\n# 标准化与对数转换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 5. 差异表达分析（每个 cluster vs rest）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top genes（wilcoxon 检验的 logfoldchanges 排序）\n    # 注意：rank_genes_groups 输出结构为字典，含 \'names\', \'scores\', \'logfoldchanges\', \'pvals\', \'pvals_adj\'\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    logfc = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster]\n    pvals_adj = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][cluster]\n    \n    # 构建 DataFrame 并按 logfoldchange 降序取 top 10（确保有足够基因）\n    top10_df = pd.DataFrame({\n        \'gene\': gene_names[:10],\n        \'logfoldchange\': logfc[:10],\n        \'pval_adj\': pvals_adj[:10]\n    }).sort_values(\'logfoldchange\', ascending=False).head(10)\n    \n    de_genes_dict[cluster] = top10_df[\'gene\'].tolist()\n\n# 保存结果到 CSV 文件（每个 cluster 一行，top10 基因用逗号分隔）\nsummary_df = pd.DataFrame([\n    {\'cluster\': cl, \'top10_de_genes\': \', \'.join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nsummary_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes_summary.csv"), index=False)\n\n# 可选：保存更详细的 DE 表（每个 cluster 单独文件）\nfor cluster, genes in de_genes_dict.items():\n    # 重新提取完整 top 10 行详细信息\n    idx = list(clusters).index(cluster)\n    detailed_df = pd.DataFrame({\n        \'gene\': adata.uns[\'rank_genes_groups\'][\'names\'][idx][:10],\n        \'logfoldchange\': adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][idx][:10],\n        \'pval_adj\': adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][idx][:10],\n        \'score\': adata.uns[\'rank_genes_groups\'][\'scores\'][idx][:10]\n    }).sort_values(\'logfoldchange\', ascending=False)\n    detailed_df.to_csv(\n        os.path.join(result_dir, f"cluster_{cluster}_top10_de_genes_detailed.csv"),\n        index=False\n    )\n\n# 打印结果概览\nprint("\\n✅ Analysis completed.")\nprint(f"Found {len(clusters)} clusters: {clusters}")\nfor cl, genes in de_genes_dict.items():\n    print(f"Cluster {cl}: {genes}")\n\n# 可选：保存带注释的 AnnData 用于后续探索\nadata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\n\n# 返回每个聚类的差异表达基因列表（按要求返回字典格式）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 15 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设为 human，mt 基因前缀为 \'MT-\')\nadata.var["mt"] = adata.var_names.str.startswith("MT-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 去除高基因数异常细胞（如doublets）\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 去除线粒体比例过高细胞（凋亡/低质量）\nadata = adata[adata.obs.n_genes_by_counts > 500, :]   # 去除基因数过少的细胞\nadata = adata[:, adata.var.n_cells_by_counts >= 10]   # 保留至少在10个细胞中表达的基因\n\n# 标准化与对数转换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 5. 差异表达分析（每个 cluster vs rest）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top genes（wilcoxon 检验的 logfoldchanges 排序）\n    # 注意：rank_genes_groups 输出结构为字典，含 \'names\', \'scores\', \'logfoldchanges\', \'pvals\', \'pvals_adj\'\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    logfc = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster]\n    pvals_adj = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][cluster]\n    \n    # 构建 DataFrame 并按 logfoldchange 降序取 top 10（确保有足够基因）\n    top10_df = pd.DataFrame({\n        \'gene\': gene_names[:10],\n        \'logfoldchange\': logfc[:10],\n        \'pval_adj\': pvals_adj[:10]\n    }).sort_values(\'logfoldchange\', ascending=False).head(10)\n    \n    de_genes_dict[cluster] = top10_df[\'gene\'].tolist()\n\n# 保存结果到 CSV 文件（每个 cluster 一行，top10 基因用逗号分隔）\nsummary_df = pd.DataFrame([\n    {\'cluster\': cl, \'top10_de_genes\': \', \'.join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nsummary_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes_summary.csv"), index=False)\n\n# 可选：保存更详细的 DE 表（每个 cluster 单独文件）\nfor cluster, genes in de_genes_dict.items():\n    # 重新提取完整 top 10 行详细信息\n    idx = list(clusters).index(cluster)\n    detailed_df = pd.DataFrame({\n        \'gene\': adata.uns[\'rank_genes_groups\'][\'names\'][idx][:10],\n        \'logfoldchange\': adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][idx][:10],\n        \'pval_adj\': adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][idx][:10],\n        \'score\': adata.uns[\'rank_genes_groups\'][\'scores\'][idx][:10]\n    }).sort_values(\'logfoldchange\', ascending=False)\n    detailed_df.to_csv(\n        os.path.join(result_dir, f"cluster_{cluster}_top10_de_genes_detailed.csv"),\n        index=False\n    )\n\n# 打印结果概览\nprint("\\n✅ Analysis completed.")\nprint(f"Found {len(clusters)} clusters: {clusters}")\nfor cl, genes in de_genes_dict.items():\n    print(f"Cluster {cl}: {genes}")\n\n# 可选：保存带注释的 AnnData 用于后续探索\nadata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\n\n# 返回每个聚类的差异表达基因列表（按要求返回字典格式）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 16 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据（使用用户提供的路径）
adata_path = "/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad"
if not os.path.exists(adata_path):
    raise FileNotFoundError(f"数据文件不存在: {adata_path}")

print("✅ 正在读取 AnnData 文件...")
adata = sc.read_h5ad(adata_path)

# 2. 预处理：QC、过滤、标准化
print("✅ 开始质量控制（QC）...")
# 添加 QC 指标（如果尚未计算）
if 'n_genes_by_counts' not in adata.obs.columns:
    adata.var['n_cells_by_counts'] = (adata.X > 0).sum(axis=0).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=0)
    adata.var['mean_counts'] = adata.X.mean(axis=0).A1 if hasattr(adata.X, 'A1') else adata.X.mean(axis=0)
    adata.obs['n_genes_by_counts'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, 'A1') else (adata.X > 0).sum(axis=1)
    adata.obs['total_counts'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, 'A1') else adata.X.sum(axis=1)

# 过滤低质量细胞（示例阈值，可根据实际调整）
sc.pp.filter_cells(adata, min_genes=200)   # 至少表达200个基因
sc.pp.filter_genes(adata, min_cells=3)      # 至少在3个细胞中表达

# 标准化与对数变换
print("✅ 开始标准化和对数变换...")
sc.pp.normalize_total(adata, target_sum=1e4)  # 归一化至10000 counts per cell
sc.pp.log1p(adata)                            # log1p transformation

# 3. PCA 降维
print("✅ 计算 PCA...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 4. UMAP 可视化
print("✅ 计算 UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# 5. Leiden 聚类
print("✅ 执行 Leiden 聚类...")
sc.tl.leiden(adata, resolution=1.0, key_added='leiden')

# 6. 差异表达分析（每个 cluster vs all others）
print("✅ 开始差异表达分析（rank test）...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=1000, use_raw=False)

# 提取每个聚类的 top 10 差异表达基因
print("✅ 提取各聚类 top 10 差异表达基因...")
cluster_de_genes = {}
for cluster_id in adata.obs['leiden'].cat.categories:
    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，兼顾 pvals_adj 显著性）
    # 使用 rank_genes_groups 输出的 DataFrame 格式提取
    try:
        # 构建结果 DataFrame（scanpy >= 1.9 兼容方式）
        gene_names = adata.uns['rank_genes_groups']['names'][cluster_id]
        logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges'][cluster_id]
        pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][cluster_id]
        
        # 合并为 DataFrame 并排序（优先按 logfoldchange，其次 pval_adj）
        df = pd.DataFrame({
            'gene': gene_names,
            'logfoldchange': logfoldchanges,
            'pval_adj': pvals_adj
        }).dropna().sort_values(['logfoldchange', 'pval_adj'], ascending=[False, True])
        
        top10_genes = df.head(10)['gene'].tolist()
        cluster_de_genes[cluster_id] = top10_genes
    except Exception as e:
        print(f"⚠️  警告：无法提取 cluster {cluster_id} 的 DE 基因：{e}")
        cluster_de_genes[cluster_id] = []

# 7. 保存结果
print("✅ 保存结果...")
# 保存聚类注释和 UMAP 坐标
adata.obs[['leiden']].to_csv(os.path.join(result_dir, "clusters.csv"))
umap_df = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP1', 'UMAP2'], index=adata.obs.index)
umap_df.to_csv(os.path.join(result_dir, "umap_coordinates.csv"))

# 保存每个 cluster 的 top 10 DE 基因
de_summary = []
for clust, genes in cluster_de_genes.items():
    for rank, gene in enumerate(genes, start=1):
        de_summary.append({'cluster': clust, 'rank': rank, 'gene': gene})
de_df = pd.DataFrame(de_summary)
de_df.to_csv(os.path.join(result_dir, "top10_de_genes_per_cluster.csv"), index=False)

# 可选：保存完整 DE 结果（前50）
try:
    sc.write_results_to_file(
        adata, 
        os.path.join(result_dir, "rank_genes_groups.xlsx"),
        n_genes=50
    )
except:
    # 兼容旧版 scanpy（无 write_results_to_file）
    pass

print("✅ 分析完成！结果已保存至:", result_dir)
print("\n📊 各聚类 top 10 差异表达基因汇总：")
for clust, genes in cluster_de_genes.items():
    print(f"Cluster {clust}: {genes}")

# 返回每个聚类的差异表达基因列表（按要求）
cluster_de_genes

# Code segment 17 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设为 human，mt 基因前缀为 \'MT-\')\nadata.var["mt"] = adata.var_names.str.startswith("MT-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 去除高基因数异常细胞（如doublets）\nadata = adata[adata.obs.pct_counts_mt < 15, :]        # 去除线粒体比例过高细胞（凋亡/低质量）\nadata = adata[adata.obs.n_genes_by_counts > 500, :]   # 去除基因数过少的细胞\nadata = adata[:, adata.var.n_cells_by_counts >= 10]   # 保留至少在10个细胞中表达的基因\n\n# 标准化与对数转换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 5. 差异表达分析（每个 cluster vs rest）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top genes（wilcoxon 检验的 logfoldchanges 排序）\n    # 注意：rank_genes_groups 输出结构为字典，含 \'names\', \'scores\', \'logfoldchanges\', \'pvals\', \'pvals_adj\'\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    logfc = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster]\n    pvals_adj = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][cluster]\n    \n    # 构建 DataFrame 并按 logfoldchange 降序取 top 10（确保有足够基因）\n    top10_df = pd.DataFrame({\n        \'gene\': gene_names[:10],\n        \'logfoldchange\': logfc[:10],\n        \'pval_adj\': pvals_adj[:10]\n    }).sort_values(\'logfoldchange\', ascending=False).head(10)\n    \n    de_genes_dict[cluster] = top10_df[\'gene\'].tolist()\n\n# 保存结果到 CSV 文件（每个 cluster 一行，top10 基因用逗号分隔）\nsummary_df = pd.DataFrame([\n    {\'cluster\': cl, \'top10_de_genes\': \', \'.join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nsummary_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes_summary.csv"), index=False)\n\n# 可选：保存更详细的 DE 表（每个 cluster 单独文件）\nfor cluster, genes in de_genes_dict.items():\n    # 重新提取完整 top 10 行详细信息\n    idx = list(clusters).index(cluster)\n    detailed_df = pd.DataFrame({\n        \'gene\': adata.uns[\'rank_genes_groups\'][\'names\'][idx][:10],\n        \'logfoldchange\': adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][idx][:10],\n        \'pval_adj\': adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][idx][:10],\n        \'score\': adata.uns[\'rank_genes_groups\'][\'scores\'][idx][:10]\n    }).sort_values(\'logfoldchange\', ascending=False)\n    detailed_df.to_csv(\n        os.path.join(result_dir, f"cluster_{cluster}_top10_de_genes_detailed.csv"),\n        index=False\n    )\n\n# 打印结果概览\nprint("\\n✅ Analysis completed.")\nprint(f"Found {len(clusters)} clusters: {clusters}")\nfor cl, genes in de_genes_dict.items():\n    print(f"Cluster {cl}: {genes}")\n\n# 可选：保存带注释的 AnnData 用于后续探索\nadata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\n\n# 返回每个聚类的差异表达基因列表（按要求返回字典格式）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 18 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据（使用用户提供的路径）\nadata_path = "/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad"\nif not os.path.exists(adata_path):\n    raise FileNotFoundError(f"数据文件不存在: {adata_path}")\n\nprint("✅ 正在读取 AnnData 文件...")\nadata = sc.read_h5ad(adata_path)\n\n# 2. 预处理：QC、过滤、标准化\nprint("✅ 开始质量控制（QC）...")\n# 添加 QC 指标（如果尚未计算）\nif \'n_genes_by_counts\' not in adata.obs.columns:\n    adata.var[\'n_cells_by_counts\'] = (adata.X > 0).sum(axis=0).A1 if hasattr(adata.X, \'A1\') else (adata.X > 0).sum(axis=0)\n    adata.var[\'mean_counts\'] = adata.X.mean(axis=0).A1 if hasattr(adata.X, \'A1\') else adata.X.mean(axis=0)\n    adata.obs[\'n_genes_by_counts\'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, \'A1\') else (adata.X > 0).sum(axis=1)\n    adata.obs[\'total_counts\'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, \'A1\') else adata.X.sum(axis=1)\n\n# 过滤低质量细胞（示例阈值，可根据实际调整）\nsc.pp.filter_cells(adata, min_genes=200)   # 至少表达200个基因\nsc.pp.filter_genes(adata, min_cells=3)      # 至少在3个细胞中表达\n\n# 标准化与对数变换\nprint("✅ 开始标准化和对数变换...")\nsc.pp.normalize_total(adata, target_sum=1e4)  # 归一化至10000 counts per cell\nsc.pp.log1p(adata)                            # log1p transformation\n\n# 3. PCA 降维\nprint("✅ 计算 PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 4. UMAP 可视化\nprint("✅ 计算 UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 5. Leiden 聚类\nprint("✅ 执行 Leiden 聚类...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 6. 差异表达分析（每个 cluster vs all others）\nprint("✅ 开始差异表达分析（rank test）...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=1000, use_raw=False)\n\n# 提取每个聚类的 top 10 差异表达基因\nprint("✅ 提取各聚类 top 10 差异表达基因...")\ncluster_de_genes = {}\nfor cluster_id in adata.obs[\'leiden\'].cat.categories:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，兼顾 pvals_adj 显著性）\n    # 使用 rank_genes_groups 输出的 DataFrame 格式提取\n    try:\n        # 构建结果 DataFrame（scanpy >= 1.9 兼容方式）\n        gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster_id]\n        logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster_id]\n        pvals_adj = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][cluster_id]\n        \n        # 合并为 DataFrame 并排序（优先按 logfoldchange，其次 pval_adj）\n        df = pd.DataFrame({\n            \'gene\': gene_names,\n            \'logfoldchange\': logfoldchanges,\n            \'pval_adj\': pvals_adj\n        }).dropna().sort_values([\'logfoldchange\', \'pval_adj\'], ascending=[False, True])\n        \n        top10_genes = df.head(10)[\'gene\'].tolist()\n        cluster_de_genes[cluster_id] = top10_genes\n    except Exception as e:\n        print(f"⚠️  警告：无法提取 cluster {cluster_id} 的 DE 基因：{e}")\n        cluster_de_genes[cluster_id] = []\n\n# 7. 保存结果\nprint("✅ 保存结果...")\n# 保存聚类注释和 UMAP 坐标\nadata.obs[[\'leiden\']].to_csv(os.path.join(result_dir, "clusters.csv"))\numap_df = pd.DataFrame(adata.obsm[\'X_umap\'], columns=[\'UMAP1\', \'UMAP2\'], index=adata.obs.index)\numap_df.to_csv(os.path.join(result_dir, "umap_coordinates.csv"))\n\n# 保存每个 cluster 的 top 10 DE 基因\nde_summary = []\nfor clust, genes in cluster_de_genes.items():\n    for rank, gene in enumerate(genes, start=1):\n        de_summary.append({\'cluster\': clust, \'rank\': rank, \'gene\': gene})\nde_df = pd.DataFrame(de_summary)\nde_df.to_csv(os.path.join(result_dir, "top10_de_genes_per_cluster.csv"), index=False)\n\n# 可选：保存完整 DE 结果（前50）\ntry:\n    sc.write_results_to_file(\n        adata, \n        os.path.join(result_dir, "rank_genes_groups.xlsx"),\n        n_genes=50\n    )\nexcept:\n    # 兼容旧版 scanpy（无 write_results_to_file）\n    pass\n\nprint("✅ 分析完成！结果已保存至:", result_dir)\nprint("\\n📊 各聚类 top 10 差异表达基因汇总：")\nfor clust, genes in cluster_de_genes.items():\n    print(f"Cluster {clust}: {genes}")\n\n# 返回每个聚类的差异表达基因列表（按要求）\ncluster_de_genes', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 19 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据（使用用户提供的路径）\nadata_path = "/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad"\nif not os.path.exists(adata_path):\n    raise FileNotFoundError(f"数据文件不存在: {adata_path}")\n\nprint("✅ 正在读取 AnnData 文件...")\nadata = sc.read_h5ad(adata_path)\n\n# 2. 预处理：QC、过滤、标准化\nprint("✅ 开始质量控制（QC）...")\n# 添加 QC 指标（如果尚未计算）\nif \'n_genes_by_counts\' not in adata.obs.columns:\n    adata.var[\'n_cells_by_counts\'] = (adata.X > 0).sum(axis=0).A1 if hasattr(adata.X, \'A1\') else (adata.X > 0).sum(axis=0)\n    adata.var[\'mean_counts\'] = adata.X.mean(axis=0).A1 if hasattr(adata.X, \'A1\') else adata.X.mean(axis=0)\n    adata.obs[\'n_genes_by_counts\'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, \'A1\') else (adata.X > 0).sum(axis=1)\n    adata.obs[\'total_counts\'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, \'A1\') else adata.X.sum(axis=1)\n\n# 过滤低质量细胞（示例阈值，可根据实际调整）\nsc.pp.filter_cells(adata, min_genes=200)   # 至少表达200个基因\nsc.pp.filter_genes(adata, min_cells=3)      # 至少在3个细胞中表达\n\n# 标准化与对数变换\nprint("✅ 开始标准化和对数变换...")\nsc.pp.normalize_total(adata, target_sum=1e4)  # 归一化至10000 counts per cell\nsc.pp.log1p(adata)                            # log1p transformation\n\n# 3. PCA 降维\nprint("✅ 计算 PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 4. UMAP 可视化\nprint("✅ 计算 UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 5. Leiden 聚类\nprint("✅ 执行 Leiden 聚类...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 6. 差异表达分析（每个 cluster vs all others）\nprint("✅ 开始差异表达分析（rank test）...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=1000, use_raw=False)\n\n# 提取每个聚类的 top 10 差异表达基因\nprint("✅ 提取各聚类 top 10 差异表达基因...")\ncluster_de_genes = {}\nfor cluster_id in adata.obs[\'leiden\'].cat.categories:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，兼顾 pvals_adj 显著性）\n    # 使用 rank_genes_groups 输出的 DataFrame 格式提取\n    try:\n        # 构建结果 DataFrame（scanpy >= 1.9 兼容方式）\n        gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster_id]\n        logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster_id]\n        pvals_adj = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][cluster_id]\n        \n        # 合并为 DataFrame 并排序（优先按 logfoldchange，其次 pval_adj）\n        df = pd.DataFrame({\n            \'gene\': gene_names,\n            \'logfoldchange\': logfoldchanges,\n            \'pval_adj\': pvals_adj\n        }).dropna().sort_values([\'logfoldchange\', \'pval_adj\'], ascending=[False, True])\n        \n        top10_genes = df.head(10)[\'gene\'].tolist()\n        cluster_de_genes[cluster_id] = top10_genes\n    except Exception as e:\n        print(f"⚠️  警告：无法提取 cluster {cluster_id} 的 DE 基因：{e}")\n        cluster_de_genes[cluster_id] = []\n\n# 7. 保存结果\nprint("✅ 保存结果...")\n# 保存聚类注释和 UMAP 坐标\nadata.obs[[\'leiden\']].to_csv(os.path.join(result_dir, "clusters.csv"))\numap_df = pd.DataFrame(adata.obsm[\'X_umap\'], columns=[\'UMAP1\', \'UMAP2\'], index=adata.obs.index)\numap_df.to_csv(os.path.join(result_dir, "umap_coordinates.csv"))\n\n# 保存每个 cluster 的 top 10 DE 基因\nde_summary = []\nfor clust, genes in cluster_de_genes.items():\n    for rank, gene in enumerate(genes, start=1):\n        de_summary.append({\'cluster\': clust, \'rank\': rank, \'gene\': gene})\nde_df = pd.DataFrame(de_summary)\nde_df.to_csv(os.path.join(result_dir, "top10_de_genes_per_cluster.csv"), index=False)\n\n# 可选：保存完整 DE 结果（前50）\ntry:\n    sc.write_results_to_file(\n        adata, \n        os.path.join(result_dir, "rank_genes_groups.xlsx"),\n        n_genes=50\n    )\nexcept:\n    # 兼容旧版 scanpy（无 write_results_to_file）\n    pass\n\nprint("✅ 分析完成！结果已保存至:", result_dir)\nprint("\\n📊 各聚类 top 10 差异表达基因汇总：")\nfor clust, genes in cluster_de_genes.items():\n    print(f"Cluster {clust}: {genes}")\n\n# 返回每个聚类的差异表达基因列表（按要求）\ncluster_de_genes', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 20 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# QC：计算线粒体基因比例（假设mt基因以'MT-'或'mt-'开头）
adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 过滤低质量细胞和基因（可根据实际数据调整阈值）
adata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限
adata = adata[adata.obs.pct_counts_mt < 20, :]        # 线粒体比例上限
adata = adata[adata.obs.total_counts > 500, :]         # 总UMI数下限

# 保留高变基因（用于后续降维）
sc.pp.normalize_total(adata, target_sum=1e4)  # 标准化至10000 UMI/cell
sc.pp.log1p(adata)                            # 对数变换
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# 2. PCA 降维
print("Step 2: Performing PCA...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP 可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# 4. Leiden 聚类
print("Step 4: Clustering with Leiden...")
sc.tl.leiden(adata, resolution=1.0, key_added='leiden')

# 5. 差异表达分析：每个聚类 vs 所有其他细胞（'t-test' 或 'wilcoxon' 更稳健）
print("Step 5: Performing differential expression analysis per cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', use_raw=False, n_genes=1000)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes per cluster...")
clusters = adata.obs['leiden'].unique().tolist()
de_genes_dict = {}

for cluster in clusters:
    # 获取该cluster的top N基因（按logfoldchange排序，或scores，这里用scores更稳定）
    gene_names = adata.uns['rank_genes_groups']['names'][cluster]
    logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges'][cluster]
    pvals = adata.uns['rank_genes_groups']['pvals'][cluster]
    
    # 构建DataFrame便于排序和筛选
    df = pd.DataFrame({
        'gene': gene_names,
        'logfoldchange': logfoldchanges,
        'pval': pvals
    }).head(10)  # rank_genes_groups已按score排序，前10即top10
    
    de_genes_dict[cluster] = df['gene'].tolist()

# 保存结果到CSV（每个cluster一行，top10基因用逗号分隔）
summary_df = pd.DataFrame([
    {'cluster': cl, 'top_10_de_genes': ', '.join(genes)} 
    for cl, genes in de_genes_dict.items()
])
summary_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)

# 可选：保存完整DE结果（所有基因）
pd.DataFrame(adata.uns['rank_genes_groups']['names']).to_csv(
    os.path.join(result_dir, "all_ranked_genes_per_cluster.csv"), index=False
)

# 输出每个聚类的top10基因（供用户直接查看）
print("\nTop 10 DE genes per Leiden cluster:")
for cluster, genes in de_genes_dict.items():
    print(f"Cluster {cluster}: {genes}")

# 可选：保存AnnData对象（含UMAP、聚类、DE结果）
adata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))

print(f"\n✅ Analysis completed. Results saved to '{result_dir}'.")

# Code segment 21 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据（使用用户提供的路径）\nadata_path = "/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad"\nif not os.path.exists(adata_path):\n    raise FileNotFoundError(f"数据文件不存在: {adata_path}")\n\nprint("✅ 正在读取 AnnData 文件...")\nadata = sc.read_h5ad(adata_path)\n\n# 2. 预处理：QC、过滤、标准化\nprint("✅ 开始质量控制（QC）...")\n# 添加 QC 指标（如果尚未计算）\nif \'n_genes_by_counts\' not in adata.obs.columns:\n    adata.var[\'n_cells_by_counts\'] = (adata.X > 0).sum(axis=0).A1 if hasattr(adata.X, \'A1\') else (adata.X > 0).sum(axis=0)\n    adata.var[\'mean_counts\'] = adata.X.mean(axis=0).A1 if hasattr(adata.X, \'A1\') else adata.X.mean(axis=0)\n    adata.obs[\'n_genes_by_counts\'] = (adata.X > 0).sum(axis=1).A1 if hasattr(adata.X, \'A1\') else (adata.X > 0).sum(axis=1)\n    adata.obs[\'total_counts\'] = adata.X.sum(axis=1).A1 if hasattr(adata.X, \'A1\') else adata.X.sum(axis=1)\n\n# 过滤低质量细胞（示例阈值，可根据实际调整）\nsc.pp.filter_cells(adata, min_genes=200)   # 至少表达200个基因\nsc.pp.filter_genes(adata, min_cells=3)      # 至少在3个细胞中表达\n\n# 标准化与对数变换\nprint("✅ 开始标准化和对数变换...")\nsc.pp.normalize_total(adata, target_sum=1e4)  # 归一化至10000 counts per cell\nsc.pp.log1p(adata)                            # log1p transformation\n\n# 3. PCA 降维\nprint("✅ 计算 PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 4. UMAP 可视化\nprint("✅ 计算 UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 5. Leiden 聚类\nprint("✅ 执行 Leiden 聚类...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 6. 差异表达分析（每个 cluster vs all others）\nprint("✅ 开始差异表达分析（rank test）...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=1000, use_raw=False)\n\n# 提取每个聚类的 top 10 差异表达基因\nprint("✅ 提取各聚类 top 10 差异表达基因...")\ncluster_de_genes = {}\nfor cluster_id in adata.obs[\'leiden\'].cat.categories:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，兼顾 pvals_adj 显著性）\n    # 使用 rank_genes_groups 输出的 DataFrame 格式提取\n    try:\n        # 构建结果 DataFrame（scanpy >= 1.9 兼容方式）\n        gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster_id]\n        logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster_id]\n        pvals_adj = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][cluster_id]\n        \n        # 合并为 DataFrame 并排序（优先按 logfoldchange，其次 pval_adj）\n        df = pd.DataFrame({\n            \'gene\': gene_names,\n            \'logfoldchange\': logfoldchanges,\n            \'pval_adj\': pvals_adj\n        }).dropna().sort_values([\'logfoldchange\', \'pval_adj\'], ascending=[False, True])\n        \n        top10_genes = df.head(10)[\'gene\'].tolist()\n        cluster_de_genes[cluster_id] = top10_genes\n    except Exception as e:\n        print(f"⚠️  警告：无法提取 cluster {cluster_id} 的 DE 基因：{e}")\n        cluster_de_genes[cluster_id] = []\n\n# 7. 保存结果\nprint("✅ 保存结果...")\n# 保存聚类注释和 UMAP 坐标\nadata.obs[[\'leiden\']].to_csv(os.path.join(result_dir, "clusters.csv"))\numap_df = pd.DataFrame(adata.obsm[\'X_umap\'], columns=[\'UMAP1\', \'UMAP2\'], index=adata.obs.index)\numap_df.to_csv(os.path.join(result_dir, "umap_coordinates.csv"))\n\n# 保存每个 cluster 的 top 10 DE 基因\nde_summary = []\nfor clust, genes in cluster_de_genes.items():\n    for rank, gene in enumerate(genes, start=1):\n        de_summary.append({\'cluster\': clust, \'rank\': rank, \'gene\': gene})\nde_df = pd.DataFrame(de_summary)\nde_df.to_csv(os.path.join(result_dir, "top10_de_genes_per_cluster.csv"), index=False)\n\n# 可选：保存完整 DE 结果（前50）\ntry:\n    sc.write_results_to_file(\n        adata, \n        os.path.join(result_dir, "rank_genes_groups.xlsx"),\n        n_genes=50\n    )\nexcept:\n    # 兼容旧版 scanpy（无 write_results_to_file）\n    pass\n\nprint("✅ 分析完成！结果已保存至:", result_dir)\nprint("\\n📊 各聚类 top 10 差异表达基因汇总：")\nfor clust, genes in cluster_de_genes.items():\n    print(f"Cluster {clust}: {genes}")\n\n# 返回每个聚类的差异表达基因列表（按要求）\ncluster_de_genes', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 22 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 20, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]         # 总UMI数下限\n\n# 保留高变基因（用于后续降维）\nsc.pp.normalize_total(adata, target_sum=1e4)  # 标准化至10000 UMI/cell\nsc.pp.log1p(adata)                            # 对数变换\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 5. 差异表达分析：每个聚类 vs 所有其他细胞（\'t-test\' 或 \'wilcoxon\' 更稳健）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该cluster的top N基因（按logfoldchange排序，或scores，这里用scores更稳定）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals\'][cluster]\n    \n    # 构建DataFrame便于排序和筛选\n    df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval\': pvals\n    }).head(10)  # rank_genes_groups已按score排序，前10即top10\n    \n    de_genes_dict[cluster] = df[\'gene\'].tolist()\n\n# 保存结果到CSV（每个cluster一行，top10基因用逗号分隔）\nsummary_df = pd.DataFrame([\n    {\'cluster\': cl, \'top_10_de_genes\': \', \'.join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nsummary_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 可选：保存完整DE结果（所有基因）\npd.DataFrame(adata.uns[\'rank_genes_groups\'][\'names\']).to_csv(\n    os.path.join(result_dir, "all_ranked_genes_per_cluster.csv"), index=False\n)\n\n# 输出每个聚类的top10基因（供用户直接查看）\nprint("\\nTop 10 DE genes per Leiden cluster:")\nfor cluster, genes in de_genes_dict.items():\n    print(f"Cluster {cluster}: {genes}")\n\n# 可选：保存AnnData对象（含UMAP、聚类、DE结果）\nadata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\n\nprint(f"\\n✅ Analysis completed. Results saved to \'{result_dir}\'.")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 23 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 20, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]         # 总UMI数下限\n\n# 保留高变基因（用于后续降维）\nsc.pp.normalize_total(adata, target_sum=1e4)  # 标准化至10000 UMI/cell\nsc.pp.log1p(adata)                            # 对数变换\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 5. 差异表达分析：每个聚类 vs 所有其他细胞（\'t-test\' 或 \'wilcoxon\' 更稳健）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该cluster的top N基因（按logfoldchange排序，或scores，这里用scores更稳定）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals\'][cluster]\n    \n    # 构建DataFrame便于排序和筛选\n    df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval\': pvals\n    }).head(10)  # rank_genes_groups已按score排序，前10即top10\n    \n    de_genes_dict[cluster] = df[\'gene\'].tolist()\n\n# 保存结果到CSV（每个cluster一行，top10基因用逗号分隔）\nsummary_df = pd.DataFrame([\n    {\'cluster\': cl, \'top_10_de_genes\': \', \'.join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nsummary_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 可选：保存完整DE结果（所有基因）\npd.DataFrame(adata.uns[\'rank_genes_groups\'][\'names\']).to_csv(\n    os.path.join(result_dir, "all_ranked_genes_per_cluster.csv"), index=False\n)\n\n# 输出每个聚类的top10基因（供用户直接查看）\nprint("\\nTop 10 DE genes per Leiden cluster:")\nfor cluster, genes in de_genes_dict.items():\n    print(f"Cluster {cluster}: {genes}")\n\n# 可选：保存AnnData对象（含UMAP、聚类、DE结果）\nadata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\n\nprint(f"\\n✅ Analysis completed. Results saved to \'{result_dir}\'.")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 24 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# QC：计算线粒体基因比例、基因数、UMI数等
adata.var["mt"] = adata.var_names.str.startswith("MT-")  # human mitochondrial genes
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
# 过滤低质量细胞和基因（可根据实际数据调整阈值）
adata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限
adata = adata[adata.obs.pct_counts_mt < 10, :]         # 线粒体比例上限
adata = adata[adata.obs.total_counts > 500, :]          # 总UMI数下限

# 标准化与对数变换
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 2. PCA 降维
print("Step 2: Performing PCA...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP 可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# 4. Leiden 聚类
print("Step 4: Clustering with Leiden...")
sc.tl.leiden(adata, resolution=1.0)  # 可根据需要调整 resolution

# 5. 差异表达分析（每个 cluster vs all others）
print("Step 5: Differential expression analysis per cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=adata.n_vars, use_raw=False)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes per cluster...")
clusters = adata.obs['leiden'].cat.categories.tolist()
de_genes_dict = {}

for cluster in clusters:
    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，或 scores；这里用 scores 更稳健）
    # rank_genes_groups 结果存储在 adata.uns['rank_genes_groups']
    gene_names = adata.uns['rank_genes_groups']['names'][cluster][:10]
    logfoldchanges = adata.uns['rank_genes_groups']['logfoldchanges'][cluster][:10]
    pvals = adata.uns['rank_genes_groups']['pvals'][cluster][:10]
    
    # 构建 DataFrame 并转为列表字典
    de_genes_df = pd.DataFrame({
        'gene': gene_names,
        'logfoldchange': logfoldchanges,
        'pval': pvals
    })
    de_genes_dict[cluster] = de_genes_df.to_dict('records')

# 保存结果到 JSON 或 CSV（每个 cluster 一个文件）
print("Saving results...")
for cluster, genes in de_genes_dict.items():
    df = pd.DataFrame(genes)
    df.to_csv(os.path.join(result_dir, f"cluster_{cluster}_top10_de_genes.csv"), index=False)

# 同时保存汇总的字典（可选）
import json
with open(os.path.join(result_dir, "all_clusters_top10_de_genes.json"), "w") as f:
    json.dump(de_genes_dict, f, indent=2, default=str)

print(f"✅ Analysis completed. Results saved to {result_dir}/")
print("Cluster-wise top 10 DE genes summary:")
for cluster in clusters:
    print(f"  Cluster {cluster}: {len(de_genes_dict[cluster])} genes")

# Code segment 25 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 20, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]         # 总UMI数下限\n\n# 保留高变基因（用于后续降维）\nsc.pp.normalize_total(adata, target_sum=1e4)  # 标准化至10000 UMI/cell\nsc.pp.log1p(adata)                            # 对数变换\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added=\'leiden\')\n\n# 5. 差异表达分析：每个聚类 vs 所有其他细胞（\'t-test\' 或 \'wilcoxon\' 更稳健）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该cluster的top N基因（按logfoldchange排序，或scores，这里用scores更稳定）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals\'][cluster]\n    \n    # 构建DataFrame便于排序和筛选\n    df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval\': pvals\n    }).head(10)  # rank_genes_groups已按score排序，前10即top10\n    \n    de_genes_dict[cluster] = df[\'gene\'].tolist()\n\n# 保存结果到CSV（每个cluster一行，top10基因用逗号分隔）\nsummary_df = pd.DataFrame([\n    {\'cluster\': cl, \'top_10_de_genes\': \', \'.join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nsummary_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 可选：保存完整DE结果（所有基因）\npd.DataFrame(adata.uns[\'rank_genes_groups\'][\'names\']).to_csv(\n    os.path.join(result_dir, "all_ranked_genes_per_cluster.csv"), index=False\n)\n\n# 输出每个聚类的top10基因（供用户直接查看）\nprint("\\nTop 10 DE genes per Leiden cluster:")\nfor cluster, genes in de_genes_dict.items():\n    print(f"Cluster {cluster}: {genes}")\n\n# 可选：保存AnnData对象（含UMAP、聚类、DE结果）\nadata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\n\nprint(f"\\n✅ Analysis completed. Results saved to \'{result_dir}\'.")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 26 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例、基因数、UMI数等\nadata.var["mt"] = adata.var_names.str.startswith("MT-")  # human mitochondrial genes\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 10, :]         # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]          # 总UMI数下限\n\n# 标准化与对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0)  # 可根据需要调整 resolution\n\n# 5. 差异表达分析（每个 cluster vs all others）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].cat.categories.tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，或 scores；这里用 scores 更稳健）\n    # rank_genes_groups 结果存储在 adata.uns[\'rank_genes_groups\']\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster][:10]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals\'][cluster][:10]\n    \n    # 构建 DataFrame 并转为列表字典\n    de_genes_df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval\': pvals\n    })\n    de_genes_dict[cluster] = de_genes_df.to_dict(\'records\')\n\n# 保存结果到 JSON 或 CSV（每个 cluster 一个文件）\nprint("Saving results...")\nfor cluster, genes in de_genes_dict.items():\n    df = pd.DataFrame(genes)\n    df.to_csv(os.path.join(result_dir, f"cluster_{cluster}_top10_de_genes.csv"), index=False)\n\n# 同时保存汇总的字典（可选）\nimport json\nwith open(os.path.join(result_dir, "all_clusters_top10_de_genes.json"), "w") as f:\n    json.dump(de_genes_dict, f, indent=2, default=str)\n\nprint(f"✅ Analysis completed. Results saved to {result_dir}/")\nprint("Cluster-wise top 10 DE genes summary:")\nfor cluster in clusters:\n    print(f"  Cluster {cluster}: {len(de_genes_dict[cluster])} genes")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 27 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例、基因数、UMI数等\nadata.var["mt"] = adata.var_names.str.startswith("MT-")  # human mitochondrial genes\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 10, :]         # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]          # 总UMI数下限\n\n# 标准化与对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0)  # 可根据需要调整 resolution\n\n# 5. 差异表达分析（每个 cluster vs all others）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].cat.categories.tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，或 scores；这里用 scores 更稳健）\n    # rank_genes_groups 结果存储在 adata.uns[\'rank_genes_groups\']\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster][:10]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals\'][cluster][:10]\n    \n    # 构建 DataFrame 并转为列表字典\n    de_genes_df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval\': pvals\n    })\n    de_genes_dict[cluster] = de_genes_df.to_dict(\'records\')\n\n# 保存结果到 JSON 或 CSV（每个 cluster 一个文件）\nprint("Saving results...")\nfor cluster, genes in de_genes_dict.items():\n    df = pd.DataFrame(genes)\n    df.to_csv(os.path.join(result_dir, f"cluster_{cluster}_top10_de_genes.csv"), index=False)\n\n# 同时保存汇总的字典（可选）\nimport json\nwith open(os.path.join(result_dir, "all_clusters_top10_de_genes.json"), "w") as f:\n    json.dump(de_genes_dict, f, indent=2, default=str)\n\nprint(f"✅ Analysis completed. Results saved to {result_dir}/")\nprint("Cluster-wise top 10 DE genes summary:")\nfor cluster in clusters:\n    print(f"  Cluster {cluster}: {len(de_genes_dict[cluster])} genes")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 28 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# QC：计算线粒体基因比例（假设mt基因以'MT-'或'mt-'开头）
adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 过滤低质量细胞和基因（可根据实际需求调整阈值）
# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%
adata = adata[adata.obs.n_genes_by_counts >= 500, :]
adata = adata[adata.obs.pct_counts_mt < 20, :]

# 标准化与对数转换
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# 高变基因筛选（用于后续降维）
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# 2. PCA 降维
print("Step 2: Performing PCA...")
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP 可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, min_dist=0.4, spread=1.0)

# 4. Leiden 聚类
print("Step 4: Clustering with Leiden...")
sc.tl.leiden(adata, resolution=1.0, key_added="leiden")

# 5. 差异表达分析（每个聚类 vs 所有其他细胞）
print("Step 5: Performing differential expression analysis per cluster...")
# 使用 't-test' 方法（快速且稳定），也可替换为 'wilcoxon' 或 'logreg'
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test', use_raw=False, n_genes=1000)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes for each cluster...")
clusters = adata.obs['leiden'].unique().tolist()
de_genes_dict = {}

for cluster in clusters:
    # 获取该 cluster 的 top 10 基因名（rank_genes_groups 结果按 group 排序）
    gene_names = adata.uns['rank_genes_groups']['names'][cluster][:10]
    # 确保返回为 Python list（避免 numpy array）
    de_genes_dict[cluster] = gene_names.tolist() if hasattr(gene_names, 'tolist') else list(gene_names)

# 保存结果到 CSV 文件（每个 cluster 一行，top 10 基因用逗号分隔）
result_df = pd.DataFrame.from_dict(de_genes_dict, orient='index', columns=[f'top_{i+1}' for i in range(10)])
result_csv = os.path.join(result_dir, "leiden_de_genes_top10.csv")
result_df.to_csv(result_csv)
print(f"✅ Top 10 DE genes per cluster saved to: {result_csv}")

# 可选：保存带注释的 AnnData 对象（含 UMAP、leiden、DE 结果）
adata.write(os.path.join(result_dir, "pbmc3k_processed.h5ad"))
print(f"✅ Processed AnnData saved to: {os.path.join(result_dir, 'pbmc3k_processed.h5ad')}")

# 返回每个聚类的差异表达基因列表（字典格式）
de_genes_dict

# Code segment 29 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例、基因数、UMI数等\nadata.var["mt"] = adata.var_names.str.startswith("MT-")  # human mitochondrial genes\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 10, :]         # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]          # 总UMI数下限\n\n# 标准化与对数变换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0)  # 可根据需要调整 resolution\n\n# 5. 差异表达分析（每个 cluster vs all others）\nprint("Step 5: Differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].cat.categories.tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top N 基因（按 logfoldchange 排序，或 scores；这里用 scores 更稳健）\n    # rank_genes_groups 结果存储在 adata.uns[\'rank_genes_groups\']\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    logfoldchanges = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster][:10]\n    pvals = adata.uns[\'rank_genes_groups\'][\'pvals\'][cluster][:10]\n    \n    # 构建 DataFrame 并转为列表字典\n    de_genes_df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfoldchanges,\n        \'pval\': pvals\n    })\n    de_genes_dict[cluster] = de_genes_df.to_dict(\'records\')\n\n# 保存结果到 JSON 或 CSV（每个 cluster 一个文件）\nprint("Saving results...")\nfor cluster, genes in de_genes_dict.items():\n    df = pd.DataFrame(genes)\n    df.to_csv(os.path.join(result_dir, f"cluster_{cluster}_top10_de_genes.csv"), index=False)\n\n# 同时保存汇总的字典（可选）\nimport json\nwith open(os.path.join(result_dir, "all_clusters_top10_de_genes.json"), "w") as f:\n    json.dump(de_genes_dict, f, indent=2, default=str)\n\nprint(f"✅ Analysis completed. Results saved to {result_dir}/")\nprint("Cluster-wise top 10 DE genes summary:")\nfor cluster in clusters:\n    print(f"  Cluster {cluster}: {len(de_genes_dict[cluster])} genes")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 30 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际需求调整阈值）\n# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%\nadata = adata[adata.obs.n_genes_by_counts >= 500, :]\nadata = adata[adata.obs.pct_counts_mt < 20, :]\n\n# 标准化与对数转换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 高变基因筛选（用于后续降维）\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.4, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added="leiden")\n\n# 5. 差异表达分析（每个聚类 vs 所有其他细胞）\nprint("Step 5: Performing differential expression analysis per cluster...")\n# 使用 \'t-test\' 方法（快速且稳定），也可替换为 \'wilcoxon\' 或 \'logreg\'\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'t-test\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes for each cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top 10 基因名（rank_genes_groups 结果按 group 排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    # 确保返回为 Python list（避免 numpy array）\n    de_genes_dict[cluster] = gene_names.tolist() if hasattr(gene_names, \'tolist\') else list(gene_names)\n\n# 保存结果到 CSV 文件（每个 cluster 一行，top 10 基因用逗号分隔）\nresult_df = pd.DataFrame.from_dict(de_genes_dict, orient=\'index\', columns=[f\'top_{i+1}\' for i in range(10)])\nresult_csv = os.path.join(result_dir, "leiden_de_genes_top10.csv")\nresult_df.to_csv(result_csv)\nprint(f"✅ Top 10 DE genes per cluster saved to: {result_csv}")\n\n# 可选：保存带注释的 AnnData 对象（含 UMAP、leiden、DE 结果）\nadata.write(os.path.join(result_dir, "pbmc3k_processed.h5ad"))\nprint(f"✅ Processed AnnData saved to: {os.path.join(result_dir, \'pbmc3k_processed.h5ad\')}")\n\n# 返回每个聚类的差异表达基因列表（字典格式）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 31 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际需求调整阈值）\n# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%\nadata = adata[adata.obs.n_genes_by_counts >= 500, :]\nadata = adata[adata.obs.pct_counts_mt < 20, :]\n\n# 标准化与对数转换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 高变基因筛选（用于后续降维）\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.4, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added="leiden")\n\n# 5. 差异表达分析（每个聚类 vs 所有其他细胞）\nprint("Step 5: Performing differential expression analysis per cluster...")\n# 使用 \'t-test\' 方法（快速且稳定），也可替换为 \'wilcoxon\' 或 \'logreg\'\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'t-test\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes for each cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top 10 基因名（rank_genes_groups 结果按 group 排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    # 确保返回为 Python list（避免 numpy array）\n    de_genes_dict[cluster] = gene_names.tolist() if hasattr(gene_names, \'tolist\') else list(gene_names)\n\n# 保存结果到 CSV 文件（每个 cluster 一行，top 10 基因用逗号分隔）\nresult_df = pd.DataFrame.from_dict(de_genes_dict, orient=\'index\', columns=[f\'top_{i+1}\' for i in range(10)])\nresult_csv = os.path.join(result_dir, "leiden_de_genes_top10.csv")\nresult_df.to_csv(result_csv)\nprint(f"✅ Top 10 DE genes per cluster saved to: {result_csv}")\n\n# 可选：保存带注释的 AnnData 对象（含 UMAP、leiden、DE 结果）\nadata.write(os.path.join(result_dir, "pbmc3k_processed.h5ad"))\nprint(f"✅ Processed AnnData saved to: {os.path.join(result_dir, \'pbmc3k_processed.h5ad\')}")\n\n# 返回每个聚类的差异表达基因列表（字典格式）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 32 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# QC：计算线粒体基因比例（假设mt基因以'MT-'或'mt-'开头）
adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 过滤低质量细胞和基因（可根据实际数据调整阈值）
adata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限
adata = adata[adata.obs.pct_counts_mt < 20, :]        # 线粒体比例上限
adata = adata[adata.obs.total_counts > 500, :]         # 总UMI数下限

# 保留高变基因（用于后续降维）
sc.pp.normalize_total(adata, target_sum=1e4)  # 标准化至10000 UMI/cell
sc.pp.log1p(adata)                            # 对数变换
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]

# 2. PCA 降维
print("Step 2: Performing PCA...")
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP 可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# 4. Leiden 聚类
print("Step 4: Clustering with Leiden...")
sc.tl.leiden(adata, resolution=0.6, key_added="leiden")

# 5. 差异表达分析（每个聚类 vs 所有其他细胞）
print("Step 5: Performing differential expression analysis per cluster...")
# 使用 scanpy 内置的 rank_genes_groups（Wilcoxon test，推荐用于单细胞）
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=1000, use_raw=False)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes per cluster...")
clusters = adata.obs['leiden'].unique().tolist()
de_genes_dict = {}

for cluster in clusters:
    # 获取该 cluster 的 top 10 基因名（rank_genes_groups 输出为 pandas DataFrame）
    gene_df = sc.get.rank_genes_groups_df(adata, group=cluster, key='rank_genes_groups')
    top10_genes = gene_df.head(10)['names'].tolist()
    de_genes_dict[cluster] = top10_genes

# 保存结果到 CSV 文件（每个 cluster 一行，genes 用逗号分隔）
result_df = pd.DataFrame([
    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} 
    for cl, genes in de_genes_dict.items()
])
result_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)

# 可选：保存完整 rank_genes_groups 结果（h5ad 元数据中已包含，也可导出为 CSV）
# sc.write_results_to_csv(adata, os.path.join(result_dir, "rank_genes_groups"))

print(f"✅ Analysis completed. Top 10 DE genes per cluster saved to {os.path.join(result_dir, 'leiden_top10_de_genes.csv')}")
print("Cluster-wise top 10 DE genes:")
for cl, genes in de_genes_dict.items():
    print(f"  Cluster {cl}: {genes}")

# 返回每个聚类的差异表达基因列表（按要求返回字典）
de_genes_dict

# Code segment 33 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际需求调整阈值）\n# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%\nadata = adata[adata.obs.n_genes_by_counts >= 500, :]\nadata = adata[adata.obs.pct_counts_mt < 20, :]\n\n# 标准化与对数转换\nsc.pp.normalize_total(adata, target_sum=1e4)\nsc.pp.log1p(adata)\n\n# 高变基因筛选（用于后续降维）\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.4, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=1.0, key_added="leiden")\n\n# 5. 差异表达分析（每个聚类 vs 所有其他细胞）\nprint("Step 5: Performing differential expression analysis per cluster...")\n# 使用 \'t-test\' 方法（快速且稳定），也可替换为 \'wilcoxon\' 或 \'logreg\'\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'t-test\', use_raw=False, n_genes=1000)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes for each cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top 10 基因名（rank_genes_groups 结果按 group 排序）\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster][:10]\n    # 确保返回为 Python list（避免 numpy array）\n    de_genes_dict[cluster] = gene_names.tolist() if hasattr(gene_names, \'tolist\') else list(gene_names)\n\n# 保存结果到 CSV 文件（每个 cluster 一行，top 10 基因用逗号分隔）\nresult_df = pd.DataFrame.from_dict(de_genes_dict, orient=\'index\', columns=[f\'top_{i+1}\' for i in range(10)])\nresult_csv = os.path.join(result_dir, "leiden_de_genes_top10.csv")\nresult_df.to_csv(result_csv)\nprint(f"✅ Top 10 DE genes per cluster saved to: {result_csv}")\n\n# 可选：保存带注释的 AnnData 对象（含 UMAP、leiden、DE 结果）\nadata.write(os.path.join(result_dir, "pbmc3k_processed.h5ad"))\nprint(f"✅ Processed AnnData saved to: {os.path.join(result_dir, \'pbmc3k_processed.h5ad\')}")\n\n# 返回每个聚类的差异表达基因列表（字典格式）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 34 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 20, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]         # 总UMI数下限\n\n# 保留高变基因（用于后续降维）\nsc.pp.normalize_total(adata, target_sum=1e4)  # 标准化至10000 UMI/cell\nsc.pp.log1p(adata)                            # 对数变换\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=0.6, key_added="leiden")\n\n# 5. 差异表达分析（每个聚类 vs 所有其他细胞）\nprint("Step 5: Performing differential expression analysis per cluster...")\n# 使用 scanpy 内置的 rank_genes_groups（Wilcoxon test，推荐用于单细胞）\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=1000, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top 10 基因名（rank_genes_groups 输出为 pandas DataFrame）\n    gene_df = sc.get.rank_genes_groups_df(adata, group=cluster, key=\'rank_genes_groups\')\n    top10_genes = gene_df.head(10)[\'names\'].tolist()\n    de_genes_dict[cluster] = top10_genes\n\n# 保存结果到 CSV 文件（每个 cluster 一行，genes 用逗号分隔）\nresult_df = pd.DataFrame([\n    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nresult_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 可选：保存完整 rank_genes_groups 结果（h5ad 元数据中已包含，也可导出为 CSV）\n# sc.write_results_to_csv(adata, os.path.join(result_dir, "rank_genes_groups"))\n\nprint(f"✅ Analysis completed. Top 10 DE genes per cluster saved to {os.path.join(result_dir, \'leiden_top10_de_genes.csv\')}")\nprint("Cluster-wise top 10 DE genes:")\nfor cl, genes in de_genes_dict.items():\n    print(f"  Cluster {cl}: {genes}")\n\n# 返回每个聚类的差异表达基因列表（按要求返回字典）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 35 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 20, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]         # 总UMI数下限\n\n# 保留高变基因（用于后续降维）\nsc.pp.normalize_total(adata, target_sum=1e4)  # 标准化至10000 UMI/cell\nsc.pp.log1p(adata)                            # 对数变换\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=0.6, key_added="leiden")\n\n# 5. 差异表达分析（每个聚类 vs 所有其他细胞）\nprint("Step 5: Performing differential expression analysis per cluster...")\n# 使用 scanpy 内置的 rank_genes_groups（Wilcoxon test，推荐用于单细胞）\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=1000, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top 10 基因名（rank_genes_groups 输出为 pandas DataFrame）\n    gene_df = sc.get.rank_genes_groups_df(adata, group=cluster, key=\'rank_genes_groups\')\n    top10_genes = gene_df.head(10)[\'names\'].tolist()\n    de_genes_dict[cluster] = top10_genes\n\n# 保存结果到 CSV 文件（每个 cluster 一行，genes 用逗号分隔）\nresult_df = pd.DataFrame([\n    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nresult_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 可选：保存完整 rank_genes_groups 结果（h5ad 元数据中已包含，也可导出为 CSV）\n# sc.write_results_to_csv(adata, os.path.join(result_dir, "rank_genes_groups"))\n\nprint(f"✅ Analysis completed. Top 10 DE genes per cluster saved to {os.path.join(result_dir, \'leiden_top10_de_genes.csv\')}")\nprint("Cluster-wise top 10 DE genes:")\nfor cl, genes in de_genes_dict.items():\n    print(f"  Cluster {cl}: {genes}")\n\n# 返回每个聚类的差异表达基因列表（按要求返回字典）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 36 (from code_dev)
# ============================================================
import os
import scanpy as sc
import pandas as pd
import numpy as np

# 设置结果保存路径
result_dir = "./result"
os.makedirs(result_dir, exist_ok=True)

# 1. 读取数据并进行预处理（QC、标准化）
print("Step 1: Loading and preprocessing data...")
adata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")

# QC：计算线粒体基因比例（假设mt基因以'MT-'或'mt-'开头）
adata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# 过滤低质量细胞和基因（可根据实际数据调整阈值）
# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%
adata = adata[adata.obs.n_genes_by_counts >= 500]
adata = adata[adata.obs.pct_counts_mt < 20]

# 过滤低表达基因（至少在10个细胞中表达）
sc.pp.filter_genes(adata, min_cells=10)

# 标准化与对数转换
sc.pp.normalize_total(adata, target_sum=1e4)  # 每个细胞总UMI数归一化至10000
sc.pp.log1p(adata)  # log1p transformation

# 2. PCA 降维
print("Step 2: Performing PCA...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]  # 仅使用高变基因进行下游分析
sc.pp.scale(adata, max_value=10)  # 标准化（z-score），限制最大值防异常值影响
sc.tl.pca(adata, svd_solver='arpack')

# 3. UMAP 可视化
print("Step 3: Computing UMAP...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata, min_dist=0.3, spread=1.0)

# 4. Leiden 聚类
print("Step 4: Clustering with Leiden algorithm...")
sc.tl.leiden(adata, resolution=1.0, key_added="leiden")

# 5. 差异表达分析（每个聚类 vs 所有其他细胞）
print("Step 5: Performing differential expression analysis per cluster...")
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon', n_genes=adata.n_vars, use_raw=False)

# 6. 提取每个聚类的 top 10 差异表达基因
print("Step 6: Extracting top 10 DE genes for each cluster...")
clusters = adata.obs['leiden'].unique()
de_genes_dict = {}

for cluster in sorted(clusters):
    # 获取该 cluster 的 top genes（wilcoxon 方法下按 logfoldchange 排序）
    # 注意：rank_genes_groups 输出结构为字典，含 'names', 'scores', 'logfoldchanges', 'pvals', 'pvals_adj'
    gene_names = adata.uns['rank_genes_groups']['names'][cluster]
    logfc = adata.uns['rank_genes_groups']['logfoldchanges'][cluster]
    pvals_adj = adata.uns['rank_genes_groups']['pvals_adj'][cluster]
    
    # 构建 DataFrame 并筛选显著（adj p < 0.05）且 logFC > 0 的前10个（更稳健）
    df = pd.DataFrame({
        'gene': gene_names,
        'logfoldchange': logfc,
        'pval_adj': pvals_adj
    }).head(100)  # 先取前100避免索引越界
    
    # 筛选显著且上调的基因（可选：也可放宽为所有top10，不强制显著性）
    sig_df = df[df['pval_adj'] < 0.05].nlargest(10, 'logfoldchange')
    if len(sig_df) < 10:
        # 若显著基因不足10个，则补上前10个（按logfoldchange排序）
        sig_df = df.nlargest(10, 'logfoldchange')
    
    de_genes_dict[cluster] = sig_df['gene'].tolist()

# 保存结果到 CSV（每个 cluster 一行，基因用逗号分隔）
print("Saving results...")
result_df = pd.DataFrame([
    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} 
    for cl, genes in de_genes_dict.items()
])
result_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)

# 可选：保存完整 rank_genes_groups 结果（便于后续查看）
import pickle
with open(os.path.join(result_dir, "rank_genes_groups_results.pkl"), "wb") as f:
    pickle.dump(adata.uns['rank_genes_groups'], f)

# 输出每个聚类的 top 10 基因（供用户直接查看）
print("\n=== Top 10 DE genes per Leiden cluster ===")
for cluster, genes in de_genes_dict.items():
    print(f"Cluster {cluster}: {genes}")

# 可选：保存 AnnData 对象（含所有分析结果）
adata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))

print(f"\n✅ Analysis completed. Results saved to '{result_dir}'.")

# Code segment 37 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\nadata = adata[adata.obs.n_genes_by_counts < 2500, :]  # 基因数上限\nadata = adata[adata.obs.pct_counts_mt < 20, :]        # 线粒体比例上限\nadata = adata[adata.obs.total_counts > 500, :]         # 总UMI数下限\n\n# 保留高变基因（用于后续降维）\nsc.pp.normalize_total(adata, target_sum=1e4)  # 标准化至10000 UMI/cell\nsc.pp.log1p(adata)                            # 对数变换\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.scale(adata, max_value=10)\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden...")\nsc.tl.leiden(adata, resolution=0.6, key_added="leiden")\n\n# 5. 差异表达分析（每个聚类 vs 所有其他细胞）\nprint("Step 5: Performing differential expression analysis per cluster...")\n# 使用 scanpy 内置的 rank_genes_groups（Wilcoxon test，推荐用于单细胞）\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=1000, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes per cluster...")\nclusters = adata.obs[\'leiden\'].unique().tolist()\nde_genes_dict = {}\n\nfor cluster in clusters:\n    # 获取该 cluster 的 top 10 基因名（rank_genes_groups 输出为 pandas DataFrame）\n    gene_df = sc.get.rank_genes_groups_df(adata, group=cluster, key=\'rank_genes_groups\')\n    top10_genes = gene_df.head(10)[\'names\'].tolist()\n    de_genes_dict[cluster] = top10_genes\n\n# 保存结果到 CSV 文件（每个 cluster 一行，genes 用逗号分隔）\nresult_df = pd.DataFrame([\n    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nresult_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 可选：保存完整 rank_genes_groups 结果（h5ad 元数据中已包含，也可导出为 CSV）\n# sc.write_results_to_csv(adata, os.path.join(result_dir, "rank_genes_groups"))\n\nprint(f"✅ Analysis completed. Top 10 DE genes per cluster saved to {os.path.join(result_dir, \'leiden_top10_de_genes.csv\')}")\nprint("Cluster-wise top 10 DE genes:")\nfor cl, genes in de_genes_dict.items():\n    print(f"  Cluster {cl}: {genes}")\n\n# 返回每个聚类的差异表达基因列表（按要求返回字典）\nde_genes_dict', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

# Code segment 38 (from code_solution)
# ============================================================
{'code': 'import os\nimport scanpy as sc\nimport pandas as pd\nimport numpy as np\n\n# 设置结果保存路径\nresult_dir = "./result"\nos.makedirs(result_dir, exist_ok=True)\n\n# 1. 读取数据并进行预处理（QC、标准化）\nprint("Step 1: Loading and preprocessing data...")\nadata = sc.read_h5ad("/home/luting/projects/mas/mas_2/data/pbmc3k.h5ad")\n\n# QC：计算线粒体基因比例（假设mt基因以\'MT-\'或\'mt-\'开头）\nadata.var["mt"] = adata.var_names.str.startswith("MT-") | adata.var_names.str.startswith("mt-")\nsc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)\n\n# 过滤低质量细胞和基因（可根据实际数据调整阈值）\n# 示例：保留至少 500 个基因表达的细胞，且线粒体比例 < 20%\nadata = adata[adata.obs.n_genes_by_counts >= 500]\nadata = adata[adata.obs.pct_counts_mt < 20]\n\n# 过滤低表达基因（至少在10个细胞中表达）\nsc.pp.filter_genes(adata, min_cells=10)\n\n# 标准化与对数转换\nsc.pp.normalize_total(adata, target_sum=1e4)  # 每个细胞总UMI数归一化至10000\nsc.pp.log1p(adata)  # log1p transformation\n\n# 2. PCA 降维\nprint("Step 2: Performing PCA...")\nsc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)\nadata = adata[:, adata.var.highly_variable]  # 仅使用高变基因进行下游分析\nsc.pp.scale(adata, max_value=10)  # 标准化（z-score），限制最大值防异常值影响\nsc.tl.pca(adata, svd_solver=\'arpack\')\n\n# 3. UMAP 可视化\nprint("Step 3: Computing UMAP...")\nsc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)\nsc.tl.umap(adata, min_dist=0.3, spread=1.0)\n\n# 4. Leiden 聚类\nprint("Step 4: Clustering with Leiden algorithm...")\nsc.tl.leiden(adata, resolution=1.0, key_added="leiden")\n\n# 5. 差异表达分析（每个聚类 vs 所有其他细胞）\nprint("Step 5: Performing differential expression analysis per cluster...")\nsc.tl.rank_genes_groups(adata, \'leiden\', method=\'wilcoxon\', n_genes=adata.n_vars, use_raw=False)\n\n# 6. 提取每个聚类的 top 10 差异表达基因\nprint("Step 6: Extracting top 10 DE genes for each cluster...")\nclusters = adata.obs[\'leiden\'].unique()\nde_genes_dict = {}\n\nfor cluster in sorted(clusters):\n    # 获取该 cluster 的 top genes（wilcoxon 方法下按 logfoldchange 排序）\n    # 注意：rank_genes_groups 输出结构为字典，含 \'names\', \'scores\', \'logfoldchanges\', \'pvals\', \'pvals_adj\'\n    gene_names = adata.uns[\'rank_genes_groups\'][\'names\'][cluster]\n    logfc = adata.uns[\'rank_genes_groups\'][\'logfoldchanges\'][cluster]\n    pvals_adj = adata.uns[\'rank_genes_groups\'][\'pvals_adj\'][cluster]\n    \n    # 构建 DataFrame 并筛选显著（adj p < 0.05）且 logFC > 0 的前10个（更稳健）\n    df = pd.DataFrame({\n        \'gene\': gene_names,\n        \'logfoldchange\': logfc,\n        \'pval_adj\': pvals_adj\n    }).head(100)  # 先取前100避免索引越界\n    \n    # 筛选显著且上调的基因（可选：也可放宽为所有top10，不强制显著性）\n    sig_df = df[df[\'pval_adj\'] < 0.05].nlargest(10, \'logfoldchange\')\n    if len(sig_df) < 10:\n        # 若显著基因不足10个，则补上前10个（按logfoldchange排序）\n        sig_df = df.nlargest(10, \'logfoldchange\')\n    \n    de_genes_dict[cluster] = sig_df[\'gene\'].tolist()\n\n# 保存结果到 CSV（每个 cluster 一行，基因用逗号分隔）\nprint("Saving results...")\nresult_df = pd.DataFrame([\n    {"cluster": cl, "top_10_de_genes": ", ".join(genes)} \n    for cl, genes in de_genes_dict.items()\n])\nresult_df.to_csv(os.path.join(result_dir, "leiden_top10_de_genes.csv"), index=False)\n\n# 可选：保存完整 rank_genes_groups 结果（便于后续查看）\nimport pickle\nwith open(os.path.join(result_dir, "rank_genes_groups_results.pkl"), "wb") as f:\n    pickle.dump(adata.uns[\'rank_genes_groups\'], f)\n\n# 输出每个聚类的 top 10 基因（供用户直接查看）\nprint("\\n=== Top 10 DE genes per Leiden cluster ===")\nfor cluster, genes in de_genes_dict.items():\n    print(f"Cluster {cluster}: {genes}")\n\n# 可选：保存 AnnData 对象（含所有分析结果）\nadata.write_h5ad(os.path.join(result_dir, "pbmc3k_processed_with_de.h5ad"))\n\nprint(f"\\n✅ Analysis completed. Results saved to \'{result_dir}\'.")', 'requirements': 'scanpy\nmatplotlib\nnumpy\npandas\nseaborn', 'task': ''}

