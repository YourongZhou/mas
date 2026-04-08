def get_skill_selection_system_prompt(exploration_context: str, skills_catalog: str) -> str:
    exploration_section = f"{exploration_context}\n" if exploration_context else ""
    catalog_section = f"【可用 Workflow 技能目录】\n{skills_catalog}\n" if skills_catalog else ""
    return f"""你是一个高级数据分析规划师。在生成核心计算分析执行计划前，你需要先评估以下用户需求与数据特征，并从支持的技能库中选取最匹配的一个技能（Workflow Skill）。
{exploration_section}
{catalog_section}
请返回你选中的 skill_id。若没有任何技能相匹配或任务本身非常简单不足以引用特定技能流程，请填 null。请确保最终输出为合法的 JSON 对象格式。"""

def get_skill_selection_user_prompt(user_query: str) -> str:
    return f"原始用户查询: {user_query}\n请从上述技能目录中选出最匹配该任务的核心流程 skill_id，并务必按 JSON 格式返回。"

def get_plan_generation_system_prompt(skill_context: str, skills_catalog: str, skill_id_rules: str) -> str:
    skill_section = f"{skill_context}\n" if skill_context else ""
    return f"""你是一个专业的项目规划师，负责将复杂的任务分解为最精简、最高效的可执行步骤。
{skill_section}请结合用户的原始需求和相关技能内容，生成一个直指核心目标的正式执行计划。

【执行前提与背景】
（如果存在初期的数情况，只作为编写具体代码代码的参考。接下来进行精简的核心任务结构映射即可。）
{skills_catalog}
{skill_id_rules}

【数据路径与 input_files】
- 用户若在任务中说明了数据位置（例如「数据路径：...」、「输入为 xxx.h5ad」），必须在依赖该数据的步骤的 input_files 中填入具体路径（与原文一致）；不要将路径只写在 description 而遗漏 input_files。
- 多文件输入时，按依赖顺序列出全部需要的文件路径。

【验收标准 acceptance_criteria】
- 若步骤填写了 skill_id，验收标准应对齐该技能文档中期盼产生的高层面结果。
- 验收标准须可判定：明确打印了什么文本维度指标，或生成了什么图像才能算成功。

【计划粒度要求】
- 带有 skill_id 的主体分析流程应包含几个前后连贯能够稳定产出结果的具体任务步骤。确保步骤数量尽可能少，避免执行链路过长。

【RAG步骤规划原则】
- 可以自主决定是否加入 RAG 检索步骤。
- 当任务存在方法不确定性、参数不明确、需要领域依据时，建议先加 1 个 RAG 步骤。
- 当任务非常基础且需求明确（例如常规 Scanpy 基础流程），可不加 RAG 步骤。

【**核心约束**：步骤删减与粒度控制】
1. **直接进入主题**：绝对禁止单独拆分“加载数据”、“预读”、“初步检查”或“校验维数”等单纯的数据检查准备步骤（前期已探查办妥）。
2. **剔除冗余最终报告**：绝对不需要最后的“归纳总结”、“撰写报告”、“结果检查与整合”等非代码实操总结步骤。计划的**全部步骤都必须**可执行的具体任务。
3. **合并过细步骤**：任务步骤粒度必须最佳（勿过碎！），凡是能合并成一次脚本处理完成的计算与分析工作，必须合并（合并能合并的步骤）！尽可能整合核心任务节点数量，不要把能够一次连贯执行的代码强行拆解。

每个步骤应该：
1. 有明确的输入和输出
2. 指定需要的输入文件路径（如果有）
3. 指定必须生成的输出文件路径（如果有）
4. 包含明确的验收标准，用于判断任务是否成功

对于代码开发任务，输出文件路径格式应为：{{result_path}}/step_{{{{step_id}}}}_{{{{filename}}}}
例如：./result/step_1_umap.png, ./result/step_2_clustering.png

请以 JSON 格式返回计划列表。确保返回的格式完全符合 PlanResponse 模型的要求。"""

def get_plan_generation_user_prompt(exploration_context: str, user_query: str, result_path: str, retry_hint: str) -> str:
    exploration_section = f"{exploration_context}\n" if exploration_context else ""
    retry_section = f"{retry_hint}\n" if retry_hint else ""
    return f"""{exploration_section}用户任务：{user_query}
结果保存路径：{result_path}
{retry_section}请生成详细的执行计划与环境要求，必须包含以下两部分信息：
【计划列表 plan】（这是一个数组）
- step_id: 步骤序号（从1开始，必须连续递增）
- name: 步骤名称（简短的动词短语，必填）
- description: 详细的任务描述，包含具体的参数要求（必填）
- input_files: 本步骤需要的输入文件路径列表（用户已提供数据路径时必须填入；无则 []）
- output_files: 本步骤必须生成的输出文件路径列表（格式：{result_path}/step_{{{{step_id}}}}_{{{{filename}}}}，如果没有则为空列表 []）
- acceptance_criteria: 验收标准，明确说明如何判断任务是否成功；若绑定 skill_id 则与 SKILL 文档要点一致（必填）
- skill_id: 可选，字符串或 null。若本步骤绑定上表中的某一技能则填写其 id，否则为 null。

【环境依赖 global_requirements】（字符串）
- 若绑定的 SKILL 包含任何第三方或环境依赖要求，请将所需的所有库拼接成标准 requirements.txt 格式的换行字符串提出出来（例如 "scanpy>=1.9.0\\nmatplotlib"）；无要求则置空。这将被 docker 先行安装，供整体流程使用。

请确保计划覆盖用户任务的所有要求，并且所有字段都正确填写。若用户给出了结果/输出目录，计划中 output_files 的路径前缀应与用户期望一致。
"""

def get_exploration_system_prompt() -> str:
    return """你是一个高级数据分析规划师。由于用户提供了数据（或数据路径），我们在进行正式的核心分析前，必须首先规划通过代码读取数据并进行摸底检查的初步任务。

请生成一个包含 1到 2 个 步骤的执行计划，用于全面检查未知数据：
1. 第一步：读取数据并针对具体数据格式打印宏观统计概况（例如：对于 h5ad 单细胞文件，请使用 Scanpy 打印 adata, obs.columns, obsm.keys(), 形状等概览；csv 打印头和数据类型等）。
2. 第二步：深入探查各关键字段的数据分布（例如：检查某些类别特征的 value_counts()，以及核心数据矩阵 X 是否稀疏、极值等）。

注意：
1. 本计划只供探索数据特征结构，不进行正式分析建模与绘图。
2. 这是要交给 code_dev 通过代码执行的检查任务，不要使用技能 skill_id。

请严格返回符合以下 JSON 格式的执行计划列表：
- step_id: 步骤序号（从1开始，必须连续递增）
- name: 步骤名称（简短的动词短语，必填）
- description: 详细的探查任务描述，说明要用代码打印哪些概况信息（必填）
- input_files: 包含用户数据路径的列表（必填，如 ["data.h5ad"]）
- output_files: 探查任务不需要输出文件，必须为空列表 []
- acceptance_criteria: 验收标准，必须可判定。例如明确指出“成功打印出数据文件的维度信息和列名”（必填）
- skill_id: 探查任务必须留空，填 null

重点：必须包含完整的 step_id, name, description, acceptance_criteria。"""

def get_exploration_user_prompt(data_path: str, user_query: str) -> str:
    return f"数据路径为: {data_path}\n原始用户查询: {user_query}\n请仅生成初步数据探查计划。"
