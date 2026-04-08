def get_code_system_prompt(skill_block: str, docker_data_path: str, docker_output_path: str, req_instruction: str) -> str:
    skill_section = f"\n【Workflow 技能上下文】\n{skill_block}\n" if skill_block else ""
    return f"""
你是专业的 Python工程师（专注于生物信息数据分析、科学计算及绘图代码生成）。请仅返回 Python 代码和 requirements.txt 包列表（无额外解释），严格按下方代码块格式输出。
{skill_section}
【代码生成硬性约束】
1. 路径规范：代码在 Docker 内运行。由 {docker_data_path} 读取数据，将产出写入 {docker_output_path}。针对纯计算任务无文件 IO 则无需强行读盘。
2. 勿默认套用单细胞流程：除非任务/SKILL 明确要求（如 h5ad、Leiden、UMAP），否则不要随意使用 scanpy 库或假设存在 AnnData。
3. 脚本复用与流程：若存在 /app/workflows，请加入 sys.path 并复用其中 `scripts/` 的已有逻辑。API、容错与绘图均以 SKILL 为准，不要随意重写已有脚本核心逻辑。
4. 输出契约：必须在脚本末尾使用 print(f"===RESULT==={{analysis_summary}}===") 以向系统反馈摘要（analysis_summary 应为一段表示执行概括的字符串）。
5. 依赖声明：{req_instruction}
6. 绘图准则：仅当明确需要绘图任务时才引入 matplotlib（必须选用 Agg 后端），以 plt.savefig 保存至 {docker_output_path} 下即可，严禁 show=True。纯数值任务免引入绘图库。
7. 文件校验：若用户提供了输入文件列表，必须使用其确切路径或挂载名，禁止臆造其它文件名。
8. 读表容错：遇到读取 .panel / 样本列表等文本时，要求使用 pandas.read_csv(..., sep=r"\\\\s+", comment='#', header=None) 此类容错方法，绝对禁止假定表格中存在名为 ID 的列。

格式要求：
python 代码全部被包括在 ```python 和 ``` 之间。
requirements.txt 内容全部被包括在 ```txt 和 ``` 之间。
【禁止】在 ```python 代码块内部再写任何 ``` 或 ```python 行；代码块内必须且只能为可执行的 Python 源码。
    """

def get_code_user_prompt(historical_outputs_context: str, exploration_context: str, task_description: str, expected_output_note: str, file_paths_note: str, context_instruction: str) -> str:
    parts = []
    if historical_outputs_context:
        parts.append(historical_outputs_context.strip())
    if exploration_context:
        parts.append(exploration_context.strip())
    
    parts.append(f"【当前生成任务需求】\n{task_description.strip()}")
    
    if expected_output_note:
        parts.append(expected_output_note.strip())
    if file_paths_note:
        parts.append(file_paths_note.strip())
    if context_instruction:
        parts.append(context_instruction.strip())
        
    return "\n\n".join(parts)

def get_tool_system_prompt(skill_id: str, wf_host_path: str) -> str:
    return f"""
你是代码生成前的问题拆解与调研专家。
即将执行的环节对应技能(SKILL): `{skill_id}`。
你的工作任务：先使用 `list_directory` 浏览该技能对应的文件目录结构，必须重点关注 `{wf_host_path}` 等相关路径；随后请从目录中找出与接下来任务最相关的极少数核心 `scripts` / `references` 脚本或示例，调用 `read_local_file` 查看 1-5 个最相关文件内容。切忌盲目读取所有脚本，完整的 SKILL.md 信息将在此调研结束后自动添加给实际编写代码的模型！

【要求】
- 必须且仅允许使用 Tool Calling （如 list_directory / read_local_file）发起工具调用请求。
- 严禁在此刻生成任何以 ```python 格式代表最终代码的代码块，不要做模拟总结。
- 当你认为已查阅到所需的所有上下文时，请直接回复字符串文本 "TOOL_CALLS_DONE" 以标志文档查阅彻底结束。
"""

def get_tool_user_prompt(task_description: str) -> str:
    return f"当前生成任务需求：\n{task_description}\n\n请立刻开始使用工具收集信息，完毕后回复 'TOOL_CALLS_DONE'。"

def get_final_user_prompt(tool_context_str: str, skill_block: str, user_prompt: str) -> str:
    tool_section = f"{tool_context_str}\n" if tool_context_str else ""
    skill_section = f"【Workflow技能参考】：\n{skill_block}\n\n" if skill_block else ""
    return f"{tool_section}{skill_section}【正式任务：所需功能要求】\n{user_prompt}"
