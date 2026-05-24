# BioAgent：面向生信分析的单主 Agent 系统

`bioagent` 是对原 `mas_2` MultiAgent System 的一次结构性改造实验。它不再使用：

```text
supervisor -> code_dev -> critic -> supervisor
```

这类多 Agent 图循环，而是改成更接近 `agent_cli_example` 的模式：

```text
用户任务
-> 单个主 Agent 维护完整上下文
-> 主 Agent 判断是否需要工具
-> 工具执行：读文件 / 查 Skill / 查 Docker 环境 / 执行代码 / 跑生信流程
-> 工具结果回灌给同一个主 Agent
-> 主 Agent 继续推理或给出最终结果
```

目标是让整体链路更清晰：**只有一个 Agent 决策者，所有外部能力都以 Tool 形式显式暴露，所有执行结果都回到同一条消息历史中。**

---

## 1. 与 `mas_2` 的关系

本目录与 `mas_2` 并列：

```text
mas/
  mas_2/       原 MultiAgent System，保留不动
  BioAgent/    新单主 Agent Loop 实现
```

新系统复用 `mas_2` 中已经沉淀好的资源：

| 资源 | 来源 | 新系统中的使用方式 |
|---|---|---|
| LLM API 配置 | `mas_2/.env` | 自动加载，复用 `OPENAI_API_KEY`、`BASE_URL`、`MODEL_NAME`、`MIMO_THINKING_TYPE` |
| Workflow Skills | `mas_2/workflows/*/SKILL.md` | 通过 `skills/registry.py` 统一扫描和读取 |
| Docker Profiles | `mas_2/docker/image_catalog.json` | 通过 `inspect_image_catalog` 和 Docker 执行工具使用 |
| 单细胞数据示例 | `mas_2/data/bmmc_b_cell.h5ad` | Notebook 示例直接调用 |
| 预装环境理念 | `py-scverse-v1`、`r-singlecell-v1` 等 | 工具执行时显式指定 profile，不静默 fallback |

---

## 2. 当前迁移状态

这不是简单复制 `mas_2/src`，而是把旧系统的能力重新映射到“单主 Agent + Tools + Skills + Sandbox”的结构中。

| `mas_2` 旧能力 | 新系统迁移状态 | 说明 |
|---|---:|---|
| Supervisor 任务规划 | 已改造 | 不再单独存在，由主 Agent 在同一循环内完成计划和工具选择 |
| Code Dev 代码生成 | 已改造 | 主 Agent 可直接生成代码并调用 `execute_python` / `execute_r` |
| Critic 审核 | 设计上移除 | 通过工具返回值、退出码、stdout/stderr 和最终输出检查形成闭环，不再额外启一个 Critic Agent |
| Workflow Skill 选择和读取 | 已迁移 | `list_workflow_skills`、`inspect_workflow_skill` |
| Docker image catalog 解析 | 已迁移 | `inspect_image_catalog`，执行工具也会检查 profile/image |
| Python Docker 执行 | 已迁移 | `execute_python` |
| R Docker 执行 | 已迁移 | `execute_r` |
| 单细胞 Scanpy 核心分析 | 已增强迁移 | `run_scanpy_singlecell_pipeline` 提供炮筒式完整流程 |
| MyGene 查询 | 已迁移 | `query_mygene` |
| 基因集富集 | 已迁移 | `gene_set_enrichment`，基于 `gseapy/Enrichr` |
| 细胞类型注释工具 | 已迁移为轻量版 | `run_celltype_annotation`，基于 marker overlap |
| RAG Researcher | 未作为独立 Agent 迁移 | 当前通过 Skill 文档读取替代；如需文献 RAG，可后续作为 Tool 加入 |
| Chainlit 前端 | 明确不迁移 | 你的要求是不需要前端，执行方式仍然是 Notebook |
| Notebook 执行 | 已支持 | `notebooks/run_agent_loop.ipynb` |
| 运行日志 | 已支持 | 每次运行写入 `BioAgent/logs/run_*.log` |

---

## 3. 新目录结构

```text
BioAgent/
  README.md
  .gitignore
  notebooks/
    run_agent_loop.ipynb
  src/
    bioagent/
      agent.py              单主 Agent 循环
      config.py             加载 mas_2/.env 和路径配置
      docker_runner.py      Docker 沙箱执行
      llm.py                ChatOpenAI 创建与运行时摘要
      logging_utils.py      运行日志
      runner.py             Notebook 入口函数
      skills/
        registry.py         扫描/读取 mas_2 workflow skills 与 Docker catalog
      tools/
        registry.py         LangChain Tool 注册
        filesystem.py       文件读取、搜索、glob
        workflow.py         workflow skill / image catalog 工具
        execution.py        Python/R 代码执行工具
        singlecell.py       Scanpy 单细胞炮筒流程
        bio.py              MyGene、富集、细胞类型注释
        schemas.py          Tool 参数 schema
```

---

## 4. Agent Loop 工作方式

主 Agent 每轮只做一件事：**要么调用工具，要么输出最终答案。**

```text
Turn 1:
  LLM 读取用户任务
  决定调用 list_workflow_skills / inspect_image_catalog / inspect_workflow_skill

Turn 2:
  LLM 看到工具结果
  对 h5ad 单细胞任务调用 run_scanpy_singlecell_pipeline

Turn 3:
  LLM 看到执行结果
  如果成功，汇总输出文件
  如果失败，基于 stderr/stdout 决定修复、重试或说明阻塞
```

相比旧的 `supervisor + code_dev + critic`，新的链路更短：

- 没有跨 Agent 状态同步问题
- 没有 Critic 驳回后上下文丢失问题
- 工具调用和结果都保留在同一条消息历史里
- 日志更容易追踪每一次 LLM 请求和每一个工具调用

---

## 5. 当前可用 Tools

### 文件与搜索

| Tool | 作用 |
|---|---|
| `list_files` | 列出项目允许范围内的文件 |
| `read_file` | 按行读取文本文件 |
| `glob_search` | 按 glob 查找文件 |
| `grep_text` | 搜索文件内容 |

### Workflow / Skill

| Tool | 作用 |
|---|---|
| `list_workflow_skills` | 列出 `mas_2/workflows` 中所有 Skill |
| `inspect_workflow_skill` | 读取某个 Skill 的 frontmatter、正文预览、scripts、references |
| `inspect_image_catalog` | 读取 `mas_2/docker/image_catalog.json` |

### 沙箱执行

| Tool | 作用 |
|---|---|
| `execute_python` | 在指定 Python Docker profile 中执行 Python 脚本 |
| `execute_r` | 在指定 R Docker profile 中执行 R 脚本 |

### 生信工具

| Tool | 作用 |
|---|---|
| `run_scanpy_singlecell_pipeline` | 针对 `.h5ad` 执行完整 Scanpy 单细胞核心分析 |
| `query_mygene` | 查询基因信息 |
| `gene_set_enrichment` | 基因集富集分析 |
| `run_celltype_annotation` | 根据 marker genes 做轻量细胞类型注释 |

---

## 6. 单细胞炮筒流程

对于这类任务：

```text
我想用 data 目录下的 bmmc_b_cell.h5ad 单细胞数据进行分析
```

主 Agent 应优先调用：

```text
run_scanpy_singlecell_pipeline
```

该工具对应 `scrnaseq-scanpy-core-analysis` Skill，默认使用：

```text
env_profile = py-scverse-v1
image       = mas/py-scverse:v1
runtime     = python
```

它会在 Docker 中执行一条完整流程：

1. 读取 `.h5ad`
2. 输出数据结构摘要
3. 计算 QC 指标
4. 生成 QC 分布图
5. 标准化 / log1p
6. 高变基因选择
7. scale + PCA
8. neighbors
9. Leiden 多分辨率聚类
10. UMAP
11. marker gene 计算
12. 导出 processed h5ad、metadata、PCA/UMAP 坐标、marker CSV
13. 生成 PDF 报告和 summary JSON/TXT

预期输出目录：

```text
BioAgent/runs/run_YYYYMMDD_HHMMSS/outputs/scrnaseq_scanpy_core/
```

典型输出文件：

```text
adata_processed.h5ad
analysis_summary.json
analysis_summary.txt
cell_metadata.csv
cell_qc_metrics.csv
cluster_markers.csv
pca_coordinates.csv
qc_histograms.png
scrna_analysis_report.pdf
umap_coordinates.csv
umap_leiden_0.8.png
```

---

## 7. Notebook 使用方式

打开：

```text
BioAgent/notebooks/run_agent_loop.ipynb
```

核心调用：

```python
import sys
from pathlib import Path

repo_root = Path(r"E:\master\novellab\Agent\code\mas")
sys.path.insert(0, str(repo_root / "BioAgent" / "src"))

from bioagent.runner import run_bmmc_singlecell_demo

result = run_bmmc_singlecell_demo(max_turns=12)
result
```

也可以自定义任务：

```python
from bioagent.runner import run_bio_agent

result = run_bio_agent(
    task="我想用 data 目录下的 bmmc_b_cell.h5ad 单细胞数据进行分析",
    data_path=r"mas_2\data\bmmc_b_cell.h5ad",
    max_turns=12,
)
result
```

---

## 8. Docker 环境要求

新系统不会静默 fallback 到错误环境。

如果调用 `run_scanpy_singlecell_pipeline`，本地必须存在：

```text
mas/py-scverse:v1
```

如果镜像不存在，工具会明确返回：

```text
Docker image not found locally: mas/py-scverse:v1
```

这比旧系统更直接：不会把 `python:3.11` legacy 容器误当成 `py-scverse-v1` 来用。

查看 Docker profiles：

```python
from bioagent.config import AgentConfig
from bioagent.skills import load_image_catalog

config = AgentConfig.from_env()
load_image_catalog(config.docker_root)
```

---

## 9. 日志

每次运行会生成：

```text
BioAgent/logs/run_YYYYMMDD_HHMMSS.log
```

日志包含：

- 实际模型配置摘要
- tools 列表
- 每轮 LLM 请求耗时
- 每次工具调用参数
- 每次工具输出预览
- Docker command
- Docker exit code
- stdout/stderr 预览
- final answer

---

## 10. 验证命令

推荐使用 `mas_agent` Conda 环境：

```powershell
$env:PYTHONDONTWRITEBYTECODE='1'
$env:PYTHONPATH='BioAgent/src'
C:\Users\WYX\.conda\envs\mas_agent\python.exe -c "from bioagent.config import AgentConfig; from bioagent.skills import list_workflow_skills, load_image_catalog; c=AgentConfig.from_env(); print(c.model_name, len(list_workflow_skills(c.workflows_root)), len(load_image_catalog(c.docker_root)))"
```

语法检查：

```powershell
C:\Users\WYX\.conda\envs\mas_agent\python.exe -c "import ast, pathlib; files=list(pathlib.Path('BioAgent/src').rglob('*.py')); [ast.parse(p.read_text(encoding='utf-8'), filename=str(p)) for p in files]; print('ast_ok', len(files))"
```

---

## 11. 当前边界

目前已经迁移的是主流程能力和常用生信工具能力。以下能力还没有作为独立 Tool 完整恢复：

- 文献 RAG 检索链路
- 原 Critic Agent 的独立 LLM 审核
- Chainlit 前端
- 对每个 Skill 的专用 deterministic runner

不过所有 `mas_2/workflows/*/SKILL.md` 都可以通过 Skill 工具读取，所有 Docker profile 都可以通过 catalog 工具识别，主 Agent 也可以基于这些 Skill 生成并执行 Python/R 代码。后续如果要继续增强，最建议逐个把高频 Skill 做成类似 `run_scanpy_singlecell_pipeline` 的专用炮筒工具。


