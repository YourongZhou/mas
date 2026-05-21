"""
Code Developer Agent 状态定义
定义代码生成和执行相关的状态字段
"""
from typing import TypedDict, Optional, Any
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
    step_tool_context: Optional[str]
    # 复用的 Docker 容器 ID
    docker_container_id: Optional[str]
    # 当前 attempt 的细分耗时
    current_attempt_timing: Optional[dict[str, Any]]
    # 历史 attempt 耗时列表
    code_dev_attempt_timings: list[dict[str, Any]]
    sandbox_backend: Optional[str]
    sandbox_endpoint: Optional[str]
    sandbox_allowed_roots: Optional[list[str]]
    code_dev_step_key: Optional[str]
    tool_context_step_key: Optional[str]
    critic_reject_count: int
    step_blocked: bool
    step_block_reason: Optional[str]
    last_failure_fingerprint: Optional[str]
    failure_streak: int
