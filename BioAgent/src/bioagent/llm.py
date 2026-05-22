from __future__ import annotations

import httpx
from langchain_openai import ChatOpenAI

from .config import AgentConfig


def build_llm(config: AgentConfig) -> ChatOpenAI:
    extra: dict = {
        "request_timeout": config.request_timeout,
    }
    if config.mimo_thinking_type in {"enabled", "disabled"}:
        extra["extra_body"] = {"thinking": {"type": config.mimo_thinking_type}}
    if (config.base_url or "").startswith("http"):
        timeout = httpx.Timeout(config.request_timeout, connect=min(60.0, config.request_timeout))
        extra["http_client"] = httpx.Client(timeout=timeout, trust_env=False)
        extra["http_async_client"] = httpx.AsyncClient(timeout=timeout, trust_env=False)

    return ChatOpenAI(
        api_key=config.api_key,
        base_url=config.base_url,
        model=config.model_name,
        temperature=config.temperature,
        streaming=False,
        **extra,
    )


def runtime_summary(config: AgentConfig) -> str:
    thinking = config.mimo_thinking_type or "(not set)"
    return (
        "LLM RUNTIME SUMMARY\n"
        f"model={config.model_name}\n"
        f"base_url={config.base_url}\n"
        f"api_key={config.mask_api_key()}\n"
        f"temperature={config.temperature}\n"
        f"request_timeout={config.request_timeout}\n"
        f"mimo_thinking_type={thinking}"
    )

