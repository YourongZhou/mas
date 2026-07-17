from __future__ import annotations

import uuid
from datetime import datetime, timezone
from typing import Any


EVENT_VERSION = 1


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def event_id(prefix: str = "event") -> str:
    return f"{prefix}_{uuid.uuid4().hex[:16]}"


def envelope_event(
    event: dict[str, Any],
    *,
    session_id: str,
    current_run_id: str = "",
    prior_events: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    value = dict(event)
    run_id = str(value.get("runId") or value.get("run_id") or current_run_id or "")
    turn = value.get("turn")
    tool_call_id = str(value.get("toolCallId") or value.get("call_id") or "")
    parent_event_id = str(value.get("parentEventId") or value.get("parent_event_id") or "")
    if not parent_event_id and value.get("type") == "tool_result" and tool_call_id:
        for prior in reversed(prior_events or []):
            if prior.get("type") == "tool_call" and str(prior.get("toolCallId") or prior.get("call_id") or "") == tool_call_id:
                parent_event_id = str(prior.get("eventId") or "")
                break
    value.update(
        {
            "eventVersion": int(value.get("eventVersion") or EVENT_VERSION),
            "eventId": str(value.get("eventId") or value.get("event_id") or event_id()),
            "sessionId": str(value.get("sessionId") or session_id),
            "runId": run_id,
            "turnId": str(value.get("turnId") or (f"{run_id}:turn:{turn}" if run_id and turn is not None else "")),
            "parentEventId": parent_event_id,
            "toolCallId": tool_call_id,
            "createdAt": str(value.get("createdAt") or value.get("created_at") or utc_now()),
        }
    )
    return value
