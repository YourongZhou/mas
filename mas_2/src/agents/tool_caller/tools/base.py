from typing import Callable, TypedDict

class ToolDefinition(TypedDict):
    """
    定义一个工具的标准结构
    """
    name: str           # 工具的唯一标识符，LLM 输出时用，例如 "mygene"
    description: str    # 给 LLM 看的说明书，告诉它什么时候用这个工具
    func: Callable      # 实际执行的 Python 函数