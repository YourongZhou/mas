from .state import SupervisorAgentState
from src.core.llm import get_llm
from src.core.state import PlanStep
from langchain_core.messages import SystemMessage, HumanMessage
import time
from pydantic import BaseModel, Field, model_validator
from typing import List

class PlanResponse(BaseModel):
    plan: List[PlanStep] = Field(..., description="完整的执行计划列表")

    @model_validator(mode="before")
    @classmethod
    def allow_list(cls, data):
        if isinstance(data, list):
            return {"plan": data}
        return data

def generate_exploration_plan(state: SupervisorAgentState, data_path: str, retry_count: int = 0, max_retries: int = 3) -> SupervisorAgentState:
    print("--- [Supervisor] 正在生成初步数据探查计划 ---")
    
    system_prompt = """你是一个高级数据分析规划师。由于用户提供了数据（或数据路径），我们在进行正式的核心分析前，必须首先规划通过代码读取数据并进行摸底检查的初步任务。

请生成一个包含 1到2个 步骤的执行计划，用于全面检查未知数据：
1. 第一步：读取数据并打印宏观统计概况（例如：对于 h5ad 单细胞文件，请使用 Scanpy 打印 adata, obs.columns, obsm.keys(), 形状等概览；csv 打印头和数据类型等）。
2. 第二步：深入探查各关键字段的数据分布（例如：检查某些类别特征的 value_counts()，以及核心数据矩阵 X 是否稀疏、极值、是否存在样本前缀 barcode 等）。

注意：
1. 本计划只供探索数据特征结构，不进行正式分析建模与绘图（因为数据还没看清结构，分析步骤必须在探查完成后再生成）。
2. 这是要交给 code_dev （即代码执行）的检查任务，不要技能 skill_id。

请严格返回符合以下 JSON 格式的执行计划列表：
- step_id: 步骤序号（从1开始，必须连续递增）
- name: 步骤名称（简短的动词短语，必填）
- description: 详细的探查任务描述，说明要用代码打印哪些概况信息（必填）
- input_files: 包含用户数据路径的列表（必填，如 ["data.h5ad"]）
- output_files: 探查任务不需要输出文件，必须为空列表 []
- acceptance_criteria: 验收标准，必须可判定。例如明确指出“成功打印出数据文件的维度信息和列名”（必填）
- skill_id: 探查任务必须留空，填 null

重点：必须包含完整的 step_id, name, description, acceptance_criteria。"""

    user_query = state.get("user_query", "")
    user_prompt = f"数据路径为: {data_path}\n原始用户查询: {user_query}\n请仅生成初步数据探查计划。"
    
    # 获取LLM
    llm_plan = get_llm(temperature=0.3)
    
    try:
        chain = llm_plan.with_structured_output(PlanResponse)
        messages = [
            SystemMessage(content=system_prompt),
            HumanMessage(content=user_prompt)
        ]
        
        print("\n" + "="*60)
        print("🔍 [Supervisor Debug] 发送给 LLM 的生成【初步数据探查】计划提示词：")
        print("[System Prompt]:")
        print(system_prompt)
        print("-" * 60)
        print("[User Prompt]:")
        print(user_prompt)
        print("="*60 + "\n")
        
        response = chain.invoke(messages)
        plan = response.plan
        
        # 为每个步骤添加 ID，并设置预定义的输入
        for i, step in enumerate(plan):
            step.step_id = i + 1
            if not step.input_files:
                step.input_files = [data_path]
            # 探索计划主要目标是执行探查代码，无需特定skill
            step.skill_id = None 
                
        state["plan"] = plan
        state["current_step_index"] = 0
        state["data_exploration_done"] = False # 探查还没执行完
        print(f"  --> 数据探查计划生成成功，共 {len(plan)} 个探查步骤。")
        
        print("\n" + "="*30 + " [初步数据探查计划结果] " + "="*30)
        for step in plan:
            print(f"  [步骤 {step.step_id}] {step.name}")
            print(f"  - 描述: {step.description}")
            print(f"  - 验收: {step.acceptance_criteria}")
            print(f"  - 输入: {step.input_files}")
            print(f"  - 输出: {step.output_files}")
            print(f"  - Skill: {step.skill_id}")
            print(f"  " + "-"*40)
        print("="*80 + "\n")
        
    except Exception as e:
        print(f"  --> [Error] 探查计划生成失败: {e}")
        if retry_count < max_retries:
            print(f"      重试 {retry_count + 1}/{max_retries}...")
            time.sleep(1)
            return generate_exploration_plan(state, data_path, retry_count + 1, max_retries)
        else:
            state["plan"] = []
            print("  --> 重试次数耗尽，无法生成探查计划。")
            
    return state
