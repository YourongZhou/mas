"""
LLM 工厂模块
提供统一的 LLM 实例创建接口
"""
import inspect
import logging
from typing import Any, Optional

import httpx
from langchain_openai import ChatOpenAI

from .config import config


logging.getLogger("httpx").setLevel(logging.WARNING)
logging.getLogger("openai").setLevel(logging.WARNING)


_LLM_RUNTIME_REGISTRY: list[dict[str, Any]] = []


def _mask_secret(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return "(empty)"
    if len(text) <= 10:
        return text[:2] + "***"
    return f"{text[:6]}...{text[-4:]}"


def _safe_secret_value(value: Any) -> str:
    if value is None:
        return ""
    getter = getattr(value, "get_secret_value", None)
    if callable(getter):
        try:
            return str(getter() or "")
        except Exception:
            return ""
    return str(value or "")


def summarize_llm_instance(llm: ChatOpenAI) -> dict[str, Any]:
    base_url = (
        getattr(llm, "openai_api_base", None)
        or getattr(llm, "base_url", None)
        or ""
    )
    api_key = _safe_secret_value(getattr(llm, "openai_api_key", None))
    extra_body = getattr(llm, "extra_body", None)
    model_kwargs = getattr(llm, "model_kwargs", None)
    request_timeout = getattr(llm, "request_timeout", None)
    http_client = getattr(llm, "http_client", None)
    return {
        "model": getattr(llm, "model_name", None) or getattr(llm, "model", None) or "",
        "temperature": getattr(llm, "temperature", None),
        "streaming": bool(getattr(llm, "streaming", False)),
        "base_url": str(base_url or ""),
        "api_key_masked": _mask_secret(api_key),
        "extra_body": extra_body if isinstance(extra_body, dict) else {},
        "model_kwargs": model_kwargs if isinstance(model_kwargs, dict) else {},
        "request_timeout": request_timeout,
        "trust_env_disabled": http_client is not None,
    }


def _registry_signature(payload: dict[str, Any]) -> tuple[Any, ...]:
    return (
        payload.get("caller"),
        payload.get("model"),
        payload.get("temperature"),
        payload.get("streaming"),
        payload.get("base_url"),
        payload.get("api_key_masked"),
        repr(payload.get("extra_body", {})),
        repr(payload.get("model_kwargs", {})),
        payload.get("request_timeout"),
        payload.get("trust_env_disabled"),
    )


def _register_llm_instance(llm: ChatOpenAI) -> None:
    caller = "unknown"
    stack = inspect.stack()
    try:
        for frame_info in stack[2:]:
            filename = str(frame_info.filename or "").replace("\\", "/")
            if "/src/core/llm.py" in filename:
                continue
            caller = f"{filename}:{frame_info.lineno}"
            break
    finally:
        del stack

    summary = summarize_llm_instance(llm)
    payload = {
        "caller": caller,
        **summary,
    }
    sig = _registry_signature(payload)
    for item in _LLM_RUNTIME_REGISTRY:
        if _registry_signature(item) == sig:
            return
    _LLM_RUNTIME_REGISTRY.append(payload)


def format_llm_runtime_banner() -> str:
    lines = [
        "",
        "=" * 88,
        "========================= [LLM RUNTIME SUMMARY] =========================",
    ]
    if not _LLM_RUNTIME_REGISTRY:
        lines.append("  (no ChatOpenAI instances registered yet)")
    else:
        for idx, item in enumerate(_LLM_RUNTIME_REGISTRY, start=1):
            lines.append(
                f"[{idx}] model={item.get('model') or '(unknown)'} | "
                f"temp={item.get('temperature')} | "
                f"streaming={item.get('streaming')} | "
                f"base_url={item.get('base_url') or '(default)'}"
            )
            lines.append(
                f"    caller={item.get('caller')} | "
                f"api_key={item.get('api_key_masked')} | "
                f"timeout={item.get('request_timeout')} | "
                f"trust_env_disabled={item.get('trust_env_disabled')}"
            )
            extra_body = item.get("extra_body") or {}
            model_kwargs = item.get("model_kwargs") or {}
            if extra_body:
                lines.append(f"    extra_body={extra_body}")
            if model_kwargs:
                lines.append(f"    model_kwargs={model_kwargs}")
    lines.append("=" * 88)
    return "\n".join(lines)


def _llm_extra_kwargs() -> dict[str, Any]:
    """可选 httpx 客户端：在 MAS_HTTPS_TRUST_ENV=false 时忽略环境代理。"""
    t = config.LLM_REQUEST_TIMEOUT
    timeout = httpx.Timeout(t, connect=min(60.0, t))
    out: dict[str, Any] = {"request_timeout": t}
    if config.MIMO_THINKING_TYPE in {"enabled", "disabled"}:
        out["extra_body"] = {"thinking": {"type": config.MIMO_THINKING_TYPE}}
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
    llm = ChatOpenAI(
        api_key=api_key or config.API_KEY,
        base_url=base_url or config.BASE_URL,
        model=model_name or config.MODEL_NAME,
        temperature=temperature if temperature is not None else config.DEFAULT_TEMPERATURE,
        streaming=streaming,
        **_llm_extra_kwargs(),
    )
    _register_llm_instance(llm)
    return llm
