"""
LLM 工厂模块
提供统一的 LLM 实例创建接口
"""
import logging
from typing import Any, Optional

import httpx
from langchain_openai import ChatOpenAI

from .config import config


logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("openai").setLevel(logging.WARNING)


def _llm_extra_kwargs() -> dict[str, Any]:
    """可选 httpx 客户端：在 MAS_HTTPS_TRUST_ENV=false 时忽略环境代理。"""
    t = config.LLM_REQUEST_TIMEOUT
    timeout = httpx.Timeout(t, connect=min(60.0, t))
    out: dict[str, Any] = {"request_timeout": t}
    if not config.HTTPS_TRUST_ENV:
        out["http_client"] = httpx.Client(trust_env=False, timeout=timeout)
        out["http_async_client"] = httpx.AsyncClient(trust_env=False, timeout=timeout)
    return out


def get_llm(
    model_name: Optional[str] = None,
    temperature: Optional[float] = None,
    base_url: Optional[str] = None,
    api_key: Optional[str] = None,
    streaming: bool = False,
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
        temperature=temperature if temperature is not None else config.DEFAULT_TEMPERATURE,
        streaming=streaming,
        **_llm_extra_kwargs(),
    )
