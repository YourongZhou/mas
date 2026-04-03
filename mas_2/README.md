# MAS-2：单细胞分析多智能体系统

基于 **LangGraph** 和 **Chainlit** 构建的多智能体系统，能够根据用户的自然语言指令，自动规划、编写代码、执行分析并汇总报告，专注于单细胞 RNA 测序（scRNA-seq）数据分析场景。

---

## 目录结构

```
mas_2/
├── app.py                        # Chainlit 前端入口
├── chainlit.md                   # 聊天界面欢迎语
├── requirements.txt              # Python 依赖
├── langgraph.json                # LangGraph Studio 配置
├── .env / .env.example           # 环境变量
├── .chainlit/
│   └── config.toml               # Chainlit 服务配置
├── chroma_db/                    # ChromaDB 向量数据库（持久化）
├── data/                         # 输入数据（如 10x Genomics hg19 矩阵）
├── notebooks/
│   └── test_single_cell_analysis.ipynb  # 流程测试 Notebook
├── results/                      # 执行结果（JSON、图片等）
├── src/
│   ├── main.py                   # 主图编排（LangGraph 主入口）
│   ├── core/
│   │   ├── config.py             # 全局配置（模型、API Key 等）
│   │   ├── state.py              # 全局状态定义（GlobalState）
│   │   └── llm.py                # LLM 工厂函数
│   └── agents/
│       ├── supervisor/           # 调度智能体
│       ├── code_dev/             # 代码开发与执行智能体
│       ├── critic/               # 代码与结果审核智能体
│       ├── rag_researcher/       # 文献检索智能体
│       └── tool_caller/          # 生物信息工具调用智能体
└── tests/
    ├── unit/
    └── integration/
```

---

## 系统架构

### 整体流程

```text
用户输入
   │
   ▼
Supervisor（规划与调度）
   │ 1. 前期数据探查（可选）：生成探查代码，经 Code Dev 执行后汇总数据特征
   │ 2. 匹配 Workflow Skill：根据用户需求和数据特征匹配内置领域技能流程
   │ 3. 全局环境与容器初始化：提前提取依赖（requirements.txt），拉起可复用的 Docker 容器
   │ 4. 制定核心步骤计划（PlanStep 列表）
   │
   ├── 根据当前计划步骤决定下一个 Worker
   │      ├──► RAG Researcher（文献/知识检索）
   │      ├──► Code Dev（代码生成与执行，利用复用的 Docker 容器）
   │      └──► Tool Caller（生物信息工具调用）
   │
   └──► Critic（产出审核）
           │
           ├── 通过 ──► 记录当前步骤输出上下文，返回 Supervisor（继续下一步）
           └── 不通过 ──► 提取反馈，返回对应 Worker（重试）
                                    │
                              Supervisor 判断全部步骤完成
                                    │
                                    ▼
                               Finalize（LLM 根据全局结果汇总最终答复与报告）
```

### 状态流转（GlobalState）

所有 Agent 共享同一个 `GlobalState`，这是系统的核心数据总线。关键字段包括：

**基础与调度控制：**
- `messages`: 对话历史轨迹
- `user_query`: 用户原始问题的记录
- `task_id`: 本次任务的唯一标识，用于结果与环境隔离
- `next_worker` / `last_worker`: 记录当前调度的流转状态

**任务规划与执行上下文：**
- `plan`: 完整的执行计划列表（包含步骤 ID、名称、描述、输入/输出文件列表以及具体的验收标准）
- `current_step_index`: 当前执行到的步骤序号
- `current_step_input` / `current_step_expected_output` / `current_step_file_paths`: 提取出来的当前执行步骤的具体上下文参数，约束 Code Dev 产出的边界
- `current_step_skill_id`: 绑定该步骤对应的领域级标准化 SOP（Workflow Skill）

**数据分析特有状态：**
- `data_exploration_done`: 标志位，指示前期对庞大数据集的自动摸底动作是否已完成
- `data_exploration_results`: 储存从原始数据中嗅探到的（如形态大小、稀疏度等）客观条件日志，指导后续正式的运算编排

**全局环境依赖（Docker 复用）：**
- `global_requirements`: 从选定的技能描述里萃取的一揽子统一依赖项配置
- `docker_container_id`: 持久化的底层服务环境，一次初始化供所有分支步骤重用，避免每次 step 反复安装重复包

**结果缓存与交互流转：**
- `completed_steps_outputs`: 将每次审核通过的正确数据/结果存入环形缓冲区，给下一次串联任务提供承上启下的输入上下文
- `pending_contribution`: Worker 即将提交的待审结果，供 Critic 审核
- `critique_feedback`: Critic 打回请求时的具象化纠偏建议
- `is_approved`: 当前阶段的准入印可，决定了能否跃迁至下一步

---

## 各智能体说明

### Supervisor（`src/agents/supervisor/graph.py`）

**核心职责**：超级大脑，负责解析需求、前期摸地、确立方案、管控上下文历史并调度后续的具体干活角色。

**核心特性**：
- **两阶段策略**：首先（可选地）触发 `data analyst` 性质的探索型策略，先通过小片段代码了解数据底层样貌，而不是盲目直接编排主流程。
- **Workflow Skill 模板匹配**：如果用户意图明确，智能匹配到本地库预置的专业流程（基于 `workflows/*/SKILL.md`）。
- **执行计划精细化（PlanStep）**：计划包含 `description`（需做什么）、`input/output_files`（路径强约束）、`acceptance_criteria`（验收标尺）。
- **提前锁定底层基础设施**：一经定调，便在此层级通过 `CodeExecutor` 提前初始化复用的持久化 Docker 容器，并执行一揽子依赖配置，此后只将 `container_id` 传递，不再折腾环境。

---

### Code Dev（`src/agents/code_dev/graph.py`）

**核心职责**：听令于前置计划约束，在此框架下针对某一步进行代码补全、生成与沙箱执行。

**核心特性**：
- **精确注入 Prompt 上下文**：不再使用模糊要求，Prompt 中会被严密穿插“步骤说明/规范输入输出/验收标准”，迫使大模型输出对齐步骤任务。
- **安全检查与容错重启**：自带 `self_reflection` 防止恶意调用（如拦截 `eval` / `exec`），若执行存在抛错（精准匹配 `Error:`、`Exception:` 后缀痕迹而过滤掉无关的 `warning`），自我拉回历史堆栈连环修复，至多尝试重试 3 次。
- **复用持久级容器与 Volume 桥接**：接收 Supervisor 移交的 `docker_container_id`，利用重用环境执行 Python 脚本，并通过 base64 桥接生成图表供链上前台渲染。

---

### Critic（`src/agents/critic/graph.py`）

**核心职责**：代码与结果验证的客观审核官。拿着每一步的验收标准 `acceptance_criteria` 来量化检验 Workers 产出。

**核心特性**：
- 它只评估“当前步骤”的标准是否达到，从不为未要求的事项横生枝节。
- 能够理解 Docker 沙箱和物理路径上的虚实映射（如 `/app/data`），防止低级路径认知带来的误判。
- 若执行日志中报告流程大体通畅收尾，默认尊重执行结果实体；遇到明确错误或者逻辑明显缺项时写入 `critique_feedback` 打回并驳回状态标识。
- **视觉能力**：对于带画图渲染类型的步骤，内部挂载 `qwen-vl-plus` 通过像素层面的特征来做质实验收（如对 UMAP 点阵或聚类簇形态做把控）。

---

### RAG Researcher（`src/agents/rag_researcher/graph.py`）

**核心职责**：充当外部文献/代码规范“备忘录”，当当前步骤脱离常规编程直觉，涉及深层生物学机制或冷门库调参时前去拉取检索。

**核心特性**：
- 直接连接并受制于本项目定制化的 Chroma 向量库资源。
- `CHROMA_PERSIST_PATH` 可通过 `src/utils/project_paths.py` 自动化相对定位，无视当前运行工作路径变动引起的文件错乱。
- 知识落库将同步刷新至 `pending_contribution` 以供调用层消化及审核层判定用。

---

### Tool Caller（`src/agents/tool_caller/graph.py`）

**核心职责**：绕过不可靠的黑盒大模型泛泛回答，通过确定性工具层级调用固化的生物信息学 API/库，再通过 LLM 重组反馈并解析专业术语。

**核心流程**：
- decision（工具路由） → execute_tool（入参洗流并执行） → interpret（大白话释义及专业推论） → END

---

## 包含的工具集（Tools）

项目中工具分为 **调用型专有工具 (Tool Caller 专属)** 以及 **系统级辅助工具** 两层架构。

### 1. 专有分析工具 (`src/agents/tool_caller/tools/`)

这些工具具有极强的生物信息领域的专有属性，专为 Tool Caller 调度设计。

*   **`run_celltype_annotation` ( `celltype.py` )**: 基于 Enrichr 数据库（PanglaoDB 等图谱）检索叠加 LLM 专家级判断逻辑组成的“双引擎”细胞类型靶向推断。具备底层基因列表的清洗、抗干扰以及打分仲裁能力。
*   **`gene_set_enrichment` ( `enrichment.py` )**: 底层利用 `gseapy` 封装而来的真实基因富集操作接口。其具备出色的正则展平清洗能力能防止大模型的幻觉杂乱 JSON 数据源格式污染，直达 Enrichr 在线数据库执行 ORA，并可在本地完成条形图表绘制。
*   **`query_mygene` ( `mygene.py` )**: 一个直连 MyGene.info 云端服务端口的信息查询桥梁，方便对于模糊靶点基因为底板提供功能概述和详细身份映射。

### 2. 系统核心底层辅助工具 (`src/tools/`)

这些工具是剥离了具体业务的高层能力切片，能跨 Agent 复用或者用于物理系统干预。

*   **`read_local_file` ( `file_tools.py` )**: 一款健壮的本地路径大文件文本态分页获取工具。它解决了物理系统长篇日志、巨大 Workflow 附属 Python 文本文献跨界加载打爆上下文视窗的痛点（内置了 `offset` / `limit` 控制与容错纠偏）。

---

## Docker 隔离环境及高级复用机制

## Docker 隔离环境及高级复用配置机制

Code Dev 利用 Docker 搭建执行沙箱已是传统，但在 `mas_2` 体系下，Docker 环境被重构为**基于全局规划期拉起（Eager Initialization）、通过 Container ID 跨步复用**的现代构型。

1. **镜像基础升级机制**：选用纯正的 `python:3.11` 为底座，而抛弃了精简的 slim 镜像。这一举措从根本上确保了 C++ 编译套件（gcc、g++）的在场，完美跨过了单细胞分析中某些刁钻包（如 `annoy` ）对原生编译底座的强依赖障碍。
2. **状态共享与环境固化（复用生命周期）**：在 Supervisor 开始分发 `plan` 前的黄金空窗期，提取 Workflow 的 `global_requirements` 依赖图谱并只行使其一次，以此实例化出唯一、持久态的 `docker_container_id` 注入到 `GlobalState` 中。之后的每一步 `execute_code`，Code Dev 将略过初始化耗时且不稳定的安装步骤直接调用该复用容器内建的 Python 解释器。
3. **挂载卷全态复原（Volume Mount策略）**：该核心容器由 `CodeExecutor` 全权接管，并默认支持无限制多数据仓的路由挂载机制（根据入参解析为 `/app/data`、`/app/data1` 等）和用于落图的 `/app/output`。更关键的是，它将系统架构层级的本地 `/workflows` 只读式软链至内建环境的 `/app/workflows` 目录下，使得大模型生成的代码可无障碍导入系统内置的高质分析函数。

---

## 前端（Chainlit）

### 位置

`mas_2/app.py`

### 启动方式

```bash
cd mas_2
chainlit run app.py -w
```

默认访问地址：`http://localhost:8000`

### 核心流程

1. 用户在聊天框输入分析需求（如"请对提供的 10x 数据进行聚类和细胞类型鉴定"）
2. `build_initial_state(user_query)` 构建初始状态
3. `graph.stream()` 流式运行主图，实时展示各 Agent 输出：
      - Supervisor：输出执行计划和当前步骤
      - Code Dev：展示生成的代码、执行日志、输出图片
      - Tool Caller：展示工具调用报告
      - Critic：展示审核结论与反馈
      - Finalize：展示 LLM 整合后的最终答复（失败时回退兜底内容）
4. 结果保存到 `results/chainlit_result_<timestamp>.json`

### 配置文件

`.chainlit/config.toml` 控制会话超时、文件上传、MCP 集成等参数。

---

## 后端（LangGraph）

### 主图入口

`src/main.py:graph`

LangGraph Studio 配置（`langgraph.json`）：

```json
{
  "graphs": {
    "agent": "src/main.py:graph"
  }
}
```

### 直接调用（不经过 Chainlit）

```python
from src.main import graph
from langchain_core.messages import HumanMessage

initial_state = {
    "messages": [HumanMessage(content="分析 data/ 目录下的单细胞数据")],
    "user_query": "分析 data/ 目录下的单细胞数据",
    "plan": [],
    "current_step_index": 0,
    "next_worker": None,
    "last_worker": None,
    "pending_contribution": None,
    "critique_feedback": None,
    "is_approved": False,
}

config = {"configurable": {"thread_id": "test-001"}}
result = graph.invoke(initial_state, config=config)
print(result["final_answer"])
```

### 使用 LangGraph Studio

```bash
cd mas_2
langgraph up
```

---

## 环境配置

### 依赖安装

```bash
cd mas_2
pip install -r requirements.txt
```

### 环境变量（`.env`）

```env
# LLM 服务（DashScope 兼容 OpenAI 接口）
OPENAI_API_KEY=your_dashscope_api_key
BASE_URL=https://dashscope.aliyuncs.com/compatible-mode/v1
MODEL_NAME=qwen-turbo-latest
TEMPERATURE=0.5

# LangSmith 追踪（可选）
LANGCHAIN_TRACING_V2=false
LANGCHAIN_API_KEY=your_langsmith_api_key
LANGCHAIN_PROJECT=mas_2
```

### Docker 环境（代码执行沙箱）

Code Dev Agent 依赖 Docker 运行代码沙箱，请确保本机已安装并启动 Docker：

```bash
docker info  # 验证 Docker 可用
```

---

## 快速开始

```bash
# 1. 安装依赖
cd mas_2
pip install -r requirements.txt

# 2. 配置环境变量
cp .env.example .env
# 编辑 .env，填入 API Key

# 3. 准备数据
# 将 10x Genomics 数据放入 data/ 目录

# 4. 启动前端
chainlit run app.py -w

# 5. 浏览器访问 http://localhost:8000，输入分析需求
```

---

## 技术栈

| 类别 | 技术 |
|------|------|
| 多智能体框架 | [LangGraph](https://github.com/langchain-ai/langgraph) |
| LLM | 阿里云百炼（Qwen-Turbo / Qwen-VL-Plus），兼容 OpenAI 接口 |
| 前端 | [Chainlit](https://github.com/Chainlit/chainlit) |
| 向量数据库 | [ChromaDB](https://www.trychroma.com/) + SentenceTransformer |
| 代码执行沙箱 | Docker |
| 生物信息库 | Scanpy, gseapy, MyGene.info API |
| 追踪监控 | LangSmith |
