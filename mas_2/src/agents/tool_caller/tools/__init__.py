import pkgutil
import os
import importlib
from langchain_core.tools import BaseTool

AVAILABLE_TOOLS = []
TOOL_MAP = {}

package_dir = os.path.dirname(__file__)

# 1. 遍历当前目录下的所有文件
for _, module_name, _ in pkgutil.iter_modules([package_dir]):
    # 跳过 base.py 和私有文件
    if module_name == "base" or module_name.startswith("_"):
        continue

    try:
        # 2. 动态导入模块
        module = importlib.import_module(f".{module_name}", package=__name__)
        
        # 3. 扫描模块内的所有属性
        for attr_name in dir(module):
            attr = getattr(module, attr_name)
            
            # === 核心逻辑 ===
            # 检查这个属性是不是 LangChain 的 BaseTool 实例
            # 被 @tool 装饰的函数，类型就是 BaseTool (或 StructuredTool)
            if isinstance(attr, BaseTool):
                # 防止重复添加（如果同一个对象被引用了多次）
                if attr not in AVAILABLE_TOOLS:
                    AVAILABLE_TOOLS.append(attr)
                    
    except Exception as e:
        print(f"[Error] Loading module {module_name}: {e}")

# 4. 生成映射字典 (LangChain tool 对象自带 .name 属性)
TOOL_MAP = {t.name: t for t in AVAILABLE_TOOLS}

print(f"[System] Auto-registered tools: {list(TOOL_MAP.keys())}")