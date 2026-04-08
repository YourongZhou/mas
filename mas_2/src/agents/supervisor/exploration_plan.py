from .state import SupervisorAgentState
from src.core.llm import get_llm
from src.core.state import PlanStep
from langchain_core.messages import SystemMessage, HumanMessage
import time
from pydantic import BaseModel, Field, model_validator
from typing import List
from .prompt import get_exploration_system_prompt, get_exploration_user_prompt

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
    
    system_prompt = get_exploration_system_prompt()
    user_query = state.get("user_query", "")
    user_prompt = get_exploration_user_prompt(data_path, user_query)
    
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
