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
-> 工具执行：读文件 / 查 Skill / 查 Docker 环境 / 通用代码执行 / Skill 工作流
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
| Skill-driven 生信工作流 | 已增强迁移 | `run_skill_workflow` 读取 Skill、生成代码、Docker 执行并按报错修复 |
| 无 Skill 通用任务 | 已支持 | `run_code_workflow` 直接按任务生成 Python/R 代码、Docker 执行并修复 |
| MyGene 查询 | 已通用化 | 不再写死专用 Tool；无匹配 Skill 时通过 `run_code_workflow` 生成代码执行 |
| 基因集富集 | 已通用化 | 优先走 `functional-enrichment-from-degs` Skill；无匹配场景时走 `run_code_workflow` |
| 细胞类型注释 | 已通用化 | 优先走单细胞相关 Skill；无匹配场景时走 `run_code_workflow` |
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
        skill_workflow.py   Skill-driven 代码生成、执行与修复闭环
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
  如果有匹配 Skill，调用 run_skill_workflow
  如果没有匹配 Skill，调用 run_code_workflow 或 execute_python / execute_r

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
| `run_skill_workflow` | 读取指定 Skill，生成代码，使用 Docker 执行，失败后按 stdout/stderr 修复重试 |
| `run_code_workflow` | 无需 Skill 的通用代码工作流，适合统计、转换、绘图、轻量分析和临时数据处理 |

### 生信任务路由

| 任务类型 | 推荐路径 |
|---|---|
| MyGene / 外部基因信息查询 | `run_code_workflow` 生成可审计脚本执行 |
| DEG 后功能富集 | `run_skill_workflow(skill_id="functional-enrichment-from-degs")` |
| 单细胞 marker / 细胞类型注释 | `run_skill_workflow` 选择 Scanpy/Seurat 单细胞 Skill |
| 没有明确 Skill 的临时生信分析 | `run_code_workflow` |

---

## 6. 工作流路由

BioAgent 不会强制所有任务都走 Skill。它的路由策略更接近 `agent_cli_example` 的主循环设计：

```text
概念解释 / 普通问答
-> 直接回答

文件检查 / 日志分析 / 搜索
-> list_files / read_file / grep_text / glob_search

有明确匹配的生信工作流
-> inspect_workflow_skill
-> run_skill_workflow

没有匹配 Skill 的通用数据分析
-> run_code_workflow

非常小的一次性脚本
-> execute_python / execute_r
```

这样单细胞、bulk RNA-seq、Seurat、生存分析等任务可以借助对应 Skill；而简单表格汇总、文件转换、画图、日志解析、临时统计等任务不需要强行套 Skill。

## 7. Skill-driven 工作流

对于这类任务：

```text
我想用 data 目录下的 bmmc_b_cell.h5ad 单细胞数据进行分析
```

主 Agent 应优先调用：

```text
run_skill_workflow
```

对单细胞 `.h5ad` 示例，Agent 会选择 `scrnaseq-scanpy-core-analysis` Skill，默认使用：

```text
env_profile = py-scverse-v1
image       = mas/py-scverse:v1
runtime     = python
```

它不是预先写死的单细胞模板，而是按以下闭环执行：

1. 读取 `mas_2/workflows/<skill_id>/SKILL.md`
2. 读取该 Skill 的 scripts / references 清单和关键预览
3. 根据用户任务、输入数据、Skill 约束和 Docker profile 生成完整 Python/R 脚本
4. 在 `/repo` 只读挂载和 `/work` 运行目录中执行
5. 所有结果写入 `/work/outputs`
6. 如果执行失败，把上一版代码、exit code、stdout/stderr 回灌给代码生成器修复
7. 成功后输出 `===RESULT===...` 摘要，并在主 Agent 中形成最终报告

预期输出目录：

```text
BioAgent/runs/run_YYYYMMDD_HHMMSS/outputs/
```

具体输出文件由 Skill 和生成代码决定；单细胞任务通常会包含 processed `.h5ad`、QC 表、UMAP/PCA 坐标、marker 结果、图像和报告。

---

## 8. Notebook 使用方式

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

## 9. Docker 环境要求

新系统不会静默 fallback 到错误环境。

如果执行单细胞 Scanpy Skill，本地必须存在：

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

## 10. 日志

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
- Skill workflow 阶段进展：选定 Skill、生成代码、执行、修复重试、完成/耗尽尝试
- Code workflow 阶段进展：生成代码、执行、修复重试、完成/耗尽尝试

---

## 11. 验证命令

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

## 12. 当前边界

目前已经迁移的是主流程能力和常用生信工具能力。以下能力还没有作为独立 Tool 完整恢复：

- 文献 RAG 检索链路
- 原 Critic Agent 的独立 LLM 审核
- Chainlit 前端
- 对每个 Skill 的专用 deterministic runner

不过所有 `mas_2/workflows/*/SKILL.md` 都可以通过 Skill 工具读取，所有 Docker profile 都可以通过 catalog 工具识别，`run_skill_workflow` 会基于这些 Skill 生成并执行 Python/R 代码。后续如果要继续增强，最建议补充每个 Skill 的输出 schema、QC gate 和结果验证器，而不是回到硬编码单流程。


