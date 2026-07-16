from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Any, Callable

from langchain_core.tools import BaseTool

from .config import AgentConfig


LangMemImport = Callable[[], tuple[Callable[..., BaseTool], Callable[..., BaseTool]]]


@dataclass
class MemoryHarness:
    enabled: bool
    namespace: tuple[str, ...]
    tools: list[BaseTool]
    disabled_reason: str = ""

    def runtime_summary(self) -> str:
        if not self.enabled:
            return f"memory=disabled reason={self.disabled_reason or 'not configured'}"
        return f"memory=enabled namespace={'/'.join(self.namespace)} tools={','.join(tool.name for tool in self.tools)}"

    def recall_text(self, query: str, *, limit: int = 5, max_chars: int = 2000) -> str:
        if not self.enabled:
            return ""
        search_tool = next((tool for tool in self.tools if tool.name == "search_memory"), None)
        if search_tool is None:
            return ""
        try:
            raw = search_tool.invoke({"query": query, "limit": limit})
        except Exception:
            return ""
        memories = _extract_memory_contents(raw)
        if not memories:
            return ""
        lines = ["Relevant BioAgent memories:"]
        for memory in memories[:limit]:
            one_line = " ".join(str(memory).split())
            if one_line:
                lines.append(f"- {one_line[:500]}")
        return "\n".join(lines)[:max_chars]


_IN_PROCESS_STORES: dict[tuple[str, ...], Any] = {}


def build_memory_harness(
    config: AgentConfig,
    *,
    import_langmem_tools: LangMemImport | None = None,
) -> MemoryHarness:
    namespace = (config.memory_namespace, config.memory_user_id)
    if not config.memory_enabled:
        return MemoryHarness(False, namespace, [], "BIOAGENT_MEMORY_ENABLED=false")

    importer = import_langmem_tools or _import_langmem_tools
    try:
        create_manage_memory_tool, create_search_memory_tool = importer()
        from langgraph.store.memory import InMemoryStore
    except ModuleNotFoundError as exc:
        return MemoryHarness(False, namespace, [], f"{exc.name or 'langmem'} is not installed")
    except Exception as exc:
        return MemoryHarness(False, namespace, [], f"memory backend unavailable: {exc}")

    store = _IN_PROCESS_STORES.get(namespace)
    if store is None:
        store = InMemoryStore()
        _IN_PROCESS_STORES[namespace] = store

    tools = [
        create_manage_memory_tool(
            namespace=namespace,
            store=store,
            instructions=(
                "Record durable BioAgent context only when it will help future runs: "
                "user preferences, project conventions, resolved environment facts, "
                "or lessons learned from repeated workflow failures. Do not store raw "
                "private datasets, secrets, full logs, or transient stack traces."
            ),
        ),
        create_search_memory_tool(namespace=namespace, store=store),
    ]
    return MemoryHarness(True, namespace, tools, "")


def _import_langmem_tools() -> tuple[Callable[..., BaseTool], Callable[..., BaseTool]]:
    from langmem import create_manage_memory_tool, create_search_memory_tool

    return create_manage_memory_tool, create_search_memory_tool


def _extract_memory_contents(raw: Any) -> list[str]:
    if raw is None:
        return []
    if isinstance(raw, str):
        stripped = raw.strip()
        if not stripped:
            return []
        try:
            return _extract_memory_contents(json.loads(stripped))
        except Exception:
            return [stripped]
    if isinstance(raw, tuple):
        items = raw[0] if raw else []
        return _extract_memory_contents(items)
    if isinstance(raw, list):
        contents: list[str] = []
        for item in raw:
            contents.extend(_extract_memory_contents(item))
        return contents
    if isinstance(raw, dict):
        value = raw.get("value")
        if isinstance(value, dict) and value.get("content"):
            return [str(value.get("content"))]
        if raw.get("content"):
            return [str(raw.get("content"))]
        return []
    content = getattr(raw, "value", None)
    if isinstance(content, dict) and content.get("content"):
        return [str(content.get("content"))]
    return []
