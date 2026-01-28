"""
Code Developer Agent 状态定义
定义代码生成和执行相关的状态字段
"""
from typing import TypedDict, Optional
from src.core.state import GlobalState


class CodeAgentState(GlobalState):
    """
    Code Developer Agent 的状态
    继承自 GlobalState，添加代码生成相关的私有字段
    """
    
    # === 代码生成相关字段 ===
    # 当前任务描述
    task: str
    # 来自 Critic 的反馈（如果是驳回重做）
    feedback: Optional[str]
    # 数据路径
    data_path: str
    # 结果路径
    result_path: str
    # 生成的代码
    scanpy_code: str
    # 依赖文件内容
    requirements_txt: str
    # 执行结果文本
    analysis_result: str
    # 是否执行成功
    success: bool
    # 内部迭代计数（用于记录重试次数）
    internal_iteration_count: int

