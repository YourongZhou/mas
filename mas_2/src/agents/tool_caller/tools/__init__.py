from .mygene import MYGENE_DEF
from .enrichment import ENRICHMENT_DEF
from .celltype import CELLTYPE_DEF

# 将所有可用工具放入一个列表
AVAILABLE_TOOLS = [
    MYGENE_DEF,
    ENRICHMENT_DEF,
    CELLTYPE_DEF
]

# 生成一个字典方便按名字查找工具
# 结果示例: {'mygene': <func>, 'add': <func>}
TOOL_MAP = {t["name"]: t for t in AVAILABLE_TOOLS}