from __future__ import annotations

import json
from pathlib import Path


class ExecutionAttemptBudget:
    def __init__(self, run_dir: Path, limit: int) -> None:
        self.limit = limit
        self.path = run_dir / "state" / "execution_attempts.json"

    @property
    def used(self) -> int:
        try:
            value = json.loads(self.path.read_text(encoding="utf-8"))
        except (FileNotFoundError, json.JSONDecodeError, OSError):
            return 0
        return max(0, int(value.get("used") or 0))

    def consume(self, *, tool: str) -> dict:
        used = self.used
        if used >= self.limit:
            return {
                "ok": False,
                "error": f"Reached global max_execution_attempts={self.limit}.",
                "attempts_used": used,
            }
        updated = used + 1
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.path.with_suffix(".tmp")
        temporary.write_text(
            json.dumps({"used": updated, "limit": self.limit, "last_tool": tool}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        temporary.replace(self.path)
        return {"ok": True, "attempts_used": updated, "attempts_remaining": self.limit - updated}
