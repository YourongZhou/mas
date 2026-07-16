from __future__ import annotations

import httpx
from langchain_anthropic import ChatAnthropic
from langchain_core.language_models.chat_models import BaseChatModel
from langchain_openai import ChatOpenAI

from .config import AgentConfig


def message_text(message: object) -> str:
    content = getattr(message, "content", "")
    if isinstance(content, str):
        return content.strip()
    if not isinstance(content, list):
        return str(content or "").strip()
    texts: list[str] = []
    for block in content:
        if isinstance(block, dict):
            block_type = block.get("type")
            text = block.get("text")
        else:
            block_type = getattr(block, "type", "")
            text = getattr(block, "text", "")
        if block_type == "text" and text:
            texts.append(str(text).strip())
    return "\n".join(text for text in texts if text)


def build_llm(config: AgentConfig) -> BaseChatModel:
    if config.provider == "anthropic":
        return ChatAnthropic(
            api_key=config.api_key,
            base_url=_anthropic_base_url(config.base_url),
            model=config.model_name,
            temperature=config.temperature,
            timeout=config.request_timeout,
            streaming=False,
        )
    if config.provider != "openai_compatible":
        raise ValueError(f"Unsupported model provider: {config.provider}")

    extra: dict = {
        "request_timeout": config.request_timeout,
    }
    extra_body: dict = {}
    if config.mimo_thinking_type in {"enabled", "disabled"}:
        extra_body["thinking"] = {"type": config.mimo_thinking_type}
    if config.chat_template_enable_thinking is not None:
        extra_body["chat_template_kwargs"] = {
            "enable_thinking": config.chat_template_enable_thinking,
        }
    if extra_body:
        extra["extra_body"] = extra_body
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


def _anthropic_base_url(base_url: str) -> str:
    normalized = (base_url or "").rstrip("/")
    if normalized.endswith("/v1"):
        return normalized[:-3]
    return normalized


def runtime_summary(config: AgentConfig) -> str:
    thinking = config.mimo_thinking_type or "(not set)"
    chat_template_thinking = (
        "(not set)"
        if config.chat_template_enable_thinking is None
        else str(config.chat_template_enable_thinking)
    )
    return (
        "LLM RUNTIME SUMMARY\n"
        f"provider={config.provider}\n"
        f"model={config.model_name}\n"
        f"base_url={config.base_url}\n"
        f"api_key={config.mask_api_key()}\n"
        f"temperature={config.temperature}\n"
        f"request_timeout={config.request_timeout}\n"
        f"mimo_thinking_type={thinking}\n"
        f"chat_template_enable_thinking={chat_template_thinking}"
    )
