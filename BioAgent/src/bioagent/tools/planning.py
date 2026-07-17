from __future__ import annotations

import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PLAN_STATUSES = {"draft", "awaiting_approval", "active", "completed", "cancelled"}
STEP_STATUSES = {"pending", "in_progress", "completed", "blocked", "cancelled"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _id(prefix: str) -> str:
    return f"{prefix}_{uuid.uuid4().hex[:12]}"


class PlanStore:
    def __init__(self, run_dir: Path) -> None:
        self.path = run_dir / "state" / "plan.json"

    def propose(
        self,
        *,
        goal: str,
        steps: list[dict[str, Any]],
        assumptions: list[str] | None = None,
        requires_approval: bool = False,
    ) -> dict[str, Any]:
        normalized_steps = []
        for item in steps:
            if hasattr(item, "model_dump"):
                item = item.model_dump()
            title = str(item.get("title") or "").strip()
            if not title:
                raise ValueError("Every plan step requires a title.")
            normalized_steps.append(
                {
                    "id": _id("step"),
                    "title": title,
                    "status": "pending",
                    "dependencies": [str(value) for value in item.get("dependencies") or [] if str(value).strip()],
                    "success_criteria": str(item.get("success_criteria") or "").strip(),
                    "note": "",
                }
            )
        if not normalized_steps:
            raise ValueError("A plan requires at least one step.")
        now = _now()
        question_id = _id("question") if requires_approval else ""
        plan = {
            "plan_id": _id("plan"),
            "goal": goal.strip(),
            "status": "awaiting_approval" if requires_approval else "active",
            "assumptions": [str(value) for value in assumptions or [] if str(value).strip()],
            "steps": normalized_steps,
            "approval_question_id": question_id,
            "created_at": now,
            "updated_at": now,
        }
        self._save(plan)
        if not requires_approval:
            return {"ok": True, "status": "active", "plan": plan}
        return {
            "ok": True,
            "status": "needs_user_input",
            "interaction_type": "plan_approval",
            "question_id": question_id,
            "question": "Please review the proposed plan. Approve it or describe the changes you need.",
            "reason": "Execution is waiting for plan approval.",
            "options": ["Approve plan", "Request changes"],
            "allow_free_text": True,
            "input_type": "single_choice",
            "required_fields": ["plan_approval"],
            "resume_hint": "If approved, call update_plan with plan_status='active' before continuing.",
            "plan": plan,
        }

    def update(
        self,
        *,
        step_id: str = "",
        step_status: str = "",
        plan_status: str = "",
        note: str = "",
    ) -> dict[str, Any]:
        plan = self.load()
        if not plan:
            return {"ok": False, "error": "No persisted plan exists for this run."}
        if plan_status:
            if plan_status not in PLAN_STATUSES:
                return {"ok": False, "error": f"Invalid plan status: {plan_status}"}
            plan["status"] = plan_status
        if step_id:
            step = next((item for item in plan.get("steps", []) if item.get("id") == step_id), None)
            if step is None:
                return {"ok": False, "error": f"Unknown plan step id: {step_id}"}
            if step_status:
                if step_status not in STEP_STATUSES:
                    return {"ok": False, "error": f"Invalid plan step status: {step_status}"}
                step["status"] = step_status
            if note:
                step["note"] = note.strip()
        steps = plan.get("steps") or []
        if steps and all(item.get("status") in {"completed", "cancelled"} for item in steps):
            plan["status"] = "completed"
        plan["updated_at"] = _now()
        self._save(plan)
        return {"ok": True, "status": plan.get("status"), "plan": plan}

    def load(self) -> dict[str, Any]:
        try:
            value = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return {}
        return value if isinstance(value, dict) else {}

    def context_text(self) -> str:
        plan = self.load()
        if not plan:
            return ""
        lines = [
            "PlanState (persisted semantic plan; keep it current with update_plan):",
            f"- plan_id: {plan.get('plan_id', '')}",
            f"- goal: {plan.get('goal', '')}",
            f"- status: {plan.get('status', '')}",
        ]
        for step in plan.get("steps") or []:
            if not isinstance(step, dict):
                continue
            detail = f" | success: {step.get('success_criteria')}" if step.get("success_criteria") else ""
            note = f" | note: {step.get('note')}" if step.get("note") else ""
            lines.append(f"- [{step.get('status', 'pending')}] {step.get('id')}: {step.get('title')}{detail}{note}")
        return "\n".join(lines)

    def _save(self, plan: dict[str, Any]) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        temp = self.path.with_suffix(".json.tmp")
        temp.write_text(json.dumps(plan, ensure_ascii=False, indent=2), encoding="utf-8")
        temp.replace(self.path)
