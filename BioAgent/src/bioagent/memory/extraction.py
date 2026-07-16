from __future__ import annotations

import re
from hashlib import sha256
from datetime import datetime, timezone

from .schemas import MemoryEpisode, TaskState


def extract_verified_episode(
    state: TaskState,
    *,
    source_run_id: str,
    execution_outcome: str | None,
    execution_tool: str = "execution",
) -> MemoryEpisode | None:
    if execution_outcome not in {"success", "failure"}:
        return None

    root_cause = ""
    if state.resolved_errors:
        root_cause = state.resolved_errors[-1]
    elif state.active_errors:
        root_cause = state.active_errors[-1]
    root_cause = _sanitize(root_cause, max_chars=400)

    runtime = _sanitize(state.runtime_environment, max_chars=120)
    skill_id = _sanitize(state.selected_skill, max_chars=160)
    if execution_outcome == "success":
        verified_fix = _sanitize(
            f"A subsequent {execution_tool} execution completed successfully"
            + (f" in {runtime}." if runtime else "."),
            max_chars=300,
        )
        lesson_subject = skill_id or "the verified workflow"
        reusable_lesson = _sanitize(
            f"Use {lesson_subject}"
            + (f" with {runtime}" if runtime else "")
            + " for tasks matching this task and data signature; execution succeeded.",
            max_chars=400,
        )
    else:
        verified_fix = ""
        reusable_lesson = _sanitize(
            "Do not reuse this execution route without first addressing the verified root cause."
            if root_cause
            else "The execution route failed without a verified reusable fix.",
            max_chars=400,
        )

    return MemoryEpisode(
        task_signature=_task_signature(state.task),
        data_signature=_data_signature(state),
        skill_id=skill_id,
        runtime_environment=runtime,
        outcome=execution_outcome,
        verified_root_cause=root_cause,
        verified_fix=verified_fix,
        reusable_lesson=reusable_lesson,
        source_run_id=_sanitize(source_run_id, max_chars=160),
        timestamp=datetime.now(timezone.utc).isoformat(timespec="seconds"),
    )


def task_and_data_signatures(state: TaskState) -> tuple[str, str]:
    return _task_signature(state.task), _data_signature(state)


def _task_signature(task: str) -> str:
    normalized = _sanitize(" ".join(task.split()), max_chars=2000).lower()
    stopwords = {
        "a",
        "an",
        "and",
        "analyze",
        "for",
        "from",
        "in",
        "of",
        "on",
        "path",
        "please",
        "run",
        "the",
        "to",
        "use",
        "with",
    }
    tokens = re.findall(r"[a-z0-9_.-]+|[\u4e00-\u9fff]", normalized)
    terms = list(dict.fromkeys(token for token in tokens if token not in stopwords))[:20]
    digest = sha256(normalized.encode("utf-8")).hexdigest()[:12]
    return f"terms={','.join(terms) or 'unspecified'}; hash={digest}"


def _data_signature(state: TaskState) -> str:
    candidates = [*state.confirmed_inputs, state.task]
    for value in candidates:
        match = re.search(r"\.([A-Za-z0-9]{1,12})(?:\b|$)", value)
        if match:
            return f"format=.{match.group(1).lower()}"
    return ""


def _sanitize(value: str, *, max_chars: int) -> str:
    text = " ".join(str(value or "").split())
    text = re.sub(r"(?:/[A-Za-z0-9._~@%+,:=-]+)+", "<path>", text)
    text = re.sub(r"[A-Za-z]:\\(?:[^\s]+\\)*[^\s]*", "<path>", text)
    return text[:max_chars]
