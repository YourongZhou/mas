from langchain_openai import ChatOpenAI
from config import config

# 初始化一个全局可用的 LLM 实例
# 这样不仅 Agent 能用，工具内部如果需要 LLM 也能直接 import
llm = ChatOpenAI(
    api_key=config.API_KEY,
    base_url=config.BASE_URL,
    model=config.MODEL_NAME,
    temperature=config.TEMPERATURE
)