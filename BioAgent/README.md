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
| Skill-driven 生信工作流 | 已改为薄工具模式 | 主 Agent 自主读取 Skill/scripts/function 签名，生成代码后用 `execute_python` / `execute_r` 执行和修复 |
| 无 Skill 通用任务 | 已支持 | 主 Agent 直接生成 Python/R 代码并用执行工具运行 |
| MyGene 查询 | 已通用化 | 不再写死专用 Tool；无匹配 Skill 时由主 Agent 生成可审计脚本执行 |
| 基因集富集 | 已通用化 | 优先读取 `functional-enrichment-from-degs` Skill；无匹配场景时由主 Agent 生成脚本执行 |
| 细胞类型注释 | 已通用化 | 优先读取单细胞相关 Skill；无匹配场景时由主 Agent 生成脚本执行 |
| RAG Researcher | 未作为独立 Agent 迁移 | 当前通过 Skill 文档读取替代；如需文献 RAG，可后续作为 Tool 加入 |
| Chainlit 前端 | 不迁移 | 旧 Chainlit UI 不再沿用 |
| 本地 Workbench 前端 | 已支持 | Biomni-style 本地页面，展示实时事件、工具调用、结果文件和历史 run |
| Notebook 执行 | 已支持 | `notebooks/run_agent_loop.ipynb` |
| 运行日志 | 已支持 | 每次运行写入 `BioAgent/logs/run_*.log` |
| LangMem Memory Harness | 已支持 | 主 Agent 暴露 `manage_memory` / `search_memory`，当前默认使用进程内 LangGraph store |
| 中途等待用户输入 | 已支持 | `request_user_input` 和 Skill preflight 可返回 `status="needs_user_input"`，通过 `resume_bio_agent()` 继续 |

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
      webapp/
        app.py              FastAPI 本地 Workbench 后端与静态文件服务
        state.py            Web 任务状态、事件映射和结果文件树
        static/             无构建依赖的前端页面
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
  如果有匹配 Skill，读取 scripts/function 签名并自己写代码
  如果没有匹配 Skill，直接写 Python/R 脚本

Turn 3:
  LLM 调用 execute_python / execute_r 并看到执行结果
  如果成功，汇总输出文件
  如果失败，基于 stderr/stdout 再查函数签名或源码，最小修复后重试
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
| `read_skill_script` | 读取某个 Skill 内脚本或 reference 的带行号片段 |
| `inspect_skill_script_symbols` | 解析 Skill 内 Python 脚本的函数、类和方法签名 |
| `inspect_skill_function` | 查看某个 Skill 函数的签名、docstring 和源码预览 |
| `inspect_image_catalog` | 读取 `mas_2/docker/image_catalog.json` |

### 沙箱执行

| Tool | 作用 |
|---|---|
| `execute_python` | 在指定 Python Docker profile 中执行 Python 脚本 |
| `execute_r` | 在指定 R Docker profile 中执行 R 脚本 |
| `request_user_input` | 缺少关键输入且不能安全默认时暂停 run，保存状态并等待用户回答后 resume |
| `manage_memory` | LangMem 记忆写入/更新/删除工具；启用 LangMem 后自动暴露 |
| `search_memory` | LangMem 记忆检索工具；启用 LangMem 后自动暴露 |

### Memory Harness

BioAgent 通过 LangMem 的 hot-path memory tools 为主 Agent 增加记忆能力：

```text
manage_memory  -> 记录用户偏好、项目约定、环境事实、可复用失败经验
search_memory  -> 在后续 run 中检索相关记忆
```

每次新 run 开始时，BioAgent 会用当前任务文本检索最多 5 条相关记忆，并把压缩后的 `Relevant BioAgent memories` 注入初始上下文。主 Agent 仍然可以在推理中主动调用 `search_memory` 或 `manage_memory`。

配置项：

| 环境变量 | 默认值 | 说明 |
|---|---|---|
| `BIOAGENT_MEMORY_ENABLED` | `true` | 是否启用 memory harness |
| `BIOAGENT_MEMORY_USER_ID` | `default` | LangMem namespace 中的用户维度 |
| `BIOAGENT_MEMORY_NAMESPACE` | `bioagent` | LangMem namespace 前缀 |

当前默认 store 是进程内 `InMemoryStore`。这是真实 LangMem 工具接入，但不是跨 Python 进程的长期持久化数据库；如果需要生产级长期记忆，下一步应把 store 替换成 Postgres/Mongo 等 LangGraph `BaseStore` 后端。

依赖安装：

```bash
python -m pip install -r BioAgent/requirements.txt
```

### 本地 Workbench 前端

BioAgent 提供一个轻量本地 Web UI，路径是 `bioagent.webapp`。它不复用旧 Chainlit，也不把分析逻辑包成 deterministic workflow；前端直接消费 `run_bio_agent_stream()` 的标准事件。

启动：

```bash
cd /home/luting/projects/mas/BioAgent
./scripts/run_workbench.sh
```

默认地址：

```text
http://127.0.0.1:8013/
```

主要接口：

| API | 作用 |
|---|---|
| `POST /api/tasks` | 创建一个 BioAgent run |
| `GET /api/tasks` | 查看历史 run 列表 |
| `GET /api/tasks/{task_id}` | 获取当前 snapshot |
| `GET /api/tasks/{task_id}/events` | SSE 实时事件流 |
| `GET /api/tasks/{task_id}/files` | 查看 run 输出文件树 |
| `GET /api/tasks/{task_id}/files/content` | 预览文本/JSON/CSV/日志 |
| `GET /api/tasks/{task_id}/files/download` | 下载或显示结果文件 |

### 中途等待用户输入与 Resume

当主 Agent 或 Skill preflight 发现缺少关键输入时，run 会返回：

```python
{
    "status": "needs_user_input",
    "pending_question": "...",
    "resume_id": "run_YYYYMMDD_HHMMSS",
    "pending_state_path": ".../pending_user_input.json",
}
```

用户回答后，在 notebook 中继续：

```python
from bioagent.runner import resume_bio_agent

result = resume_bio_agent(
    resume_id="run_YYYYMMDD_HHMMSS",
    user_answer="species is mouse; use /path/to/data.h5ad",
)
```

`resume_bio_agent()` 会恢复同一个 run 的消息历史、追加用户回答，并继续原来的单主 Agent loop。

### 生信任务路由

| 任务类型 | 推荐路径 |
|---|---|
| MyGene / 外部基因信息查询 | 主 Agent 生成可审计脚本并用 `execute_python` 执行 |
| DEG 后功能富集 | 读取 `functional-enrichment-from-degs` Skill 和 scripts 后生成脚本执行 |
| 单细胞 marker / 细胞类型注释 | 读取 Scanpy/Seurat 单细胞 Skill、scripts 和函数签名后生成脚本执行 |
| 没有明确 Skill 的临时生信分析 | 主 Agent 直接生成 Python/R 脚本并执行 |

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
-> inspect_skill_script_symbols / inspect_skill_function / read_skill_script
-> execute_python / execute_r

没有匹配 Skill 的通用数据分析
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
list_workflow_skills
inspect_workflow_skill
inspect_skill_script_symbols / inspect_skill_function / read_skill_script
execute_python
```

对单细胞 `.h5ad` 示例，Agent 会选择 `scrnaseq-scanpy-core-analysis` Skill，默认使用：

```text
env_profile = py-scverse-v1
image       = mas/py-scverse:v1
runtime     = python
```

它不是预先写死的单细胞模板，也不再通过一个胖 workflow tool 包办，而是按以下闭环执行：

1. 读取 `mas_2/workflows/<skill_id>/SKILL.md`
2. 按需读取该 Skill 的 scripts / references、函数签名和源码片段
3. 主 Agent 根据用户任务、输入数据、Skill 约束和 Docker profile 生成完整 Python/R 脚本
4. 在 `/repo` 只读挂载和 `/work` 运行目录中执行
5. 所有结果写入 `/work/outputs`
6. 如果执行失败，主 Agent 根据 exit code、stdout/stderr 再查相关函数签名或源码后最小修复
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
- Skill 读取阶段进展：选定 Skill、读取 scripts、检查函数签名和源码片段
- 执行阶段进展：生成代码、执行、修复重试、完成/耗尽全局 attempt 预算

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

不过所有 `mas_2/workflows/*/SKILL.md` 都可以通过 Skill 工具读取，所有 Docker profile 都可以通过 catalog 工具识别。默认路径是主 Agent 自主读取 Skill/scripts/function 签名，生成 Python/R 代码并用执行工具迭代。后续如果要继续增强，最建议补充每个 Skill 的输出 schema、QC gate 和结果验证器，而不是回到硬编码单流程。
