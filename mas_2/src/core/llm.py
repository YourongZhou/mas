"""
LLM 工厂模块
提供统一的 LLM 实例创建接口
"""
from langchain_openai import ChatOpenAI
from typing import Optional
from .config import config


def get_llm(
    model_name: Optional[str] = None,
    temperature: Optional[float] = None,
    base_url: Optional[str] = None,
    api_key: Optional[str] = None
) -> ChatOpenAI:
    """
    创建并返回配置好的 ChatOpenAI 实例
    
    Args:
        model_name: 模型名称，默认使用 config.MODEL_NAME
        temperature: 温度参数，默认使用 config.DEFAULT_TEMPERATURE
        base_url: API 基础 URL，默认使用 config.BASE_URL
        api_key: API 密钥，默认使用 config.API_KEY
    
    Returns:
        配置好的 ChatOpenAI 实例
    """
    return ChatOpenAI(
        api_key=api_key or config.API_KEY,
        base_url=base_url or config.BASE_URL,
        model=model_name or config.MODEL_NAME,
        temperature=temperature if temperature is not None else config.DEFAULT_TEMPERATURE
    )

