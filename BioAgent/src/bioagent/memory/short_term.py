from __future__ import annotations

from typing import Any

from .schemas import TaskState


EXECUTION_TOOLS = {"execute_python", "execute_r", "start_job", "poll_job", "run_code_workflow", "run_skill_workflow"}
SKILL_DETAIL_TOOLS = {"read_skill_script", "inspect_skill_script_symbols", "inspect_skill_function"}
JOB_MANAGEMENT_TOOLS = {"list_jobs", "get_job", "tail_job", "cancel_job"}


class ShortTermMemory:
    def __init__(self, state: TaskState) -> None:
        self.state = state

    @classmethod
    def start(cls, *, task: str, data_path: str | None = None, result_dir: str | None = None) -> "ShortTermMemory":
        confirmed_inputs = [f"Primary data: {data_path}"] if data_path else []
        expected_outputs = ["Write artifacts under /work/outputs for frontend display."] if result_dir else []
        return cls(
            TaskState(
                task=task.strip(),
                current_goal=task.strip(),
                confirmed_inputs=confirmed_inputs,
                expected_outputs=expected_outputs,
                next_action="Select a relevant skill or inspect the confirmed inputs.",
            )
        )

    @classmethod
    def restore(cls, value: dict[str, Any]) -> "ShortTermMemory":
        return cls(TaskState.from_dict(value))

    def update(self, **changes: Any) -> None:
        for name, value in changes.items():
            if not hasattr(self.state, name):
                raise AttributeError(f"Unknown TaskState field: {name}")
            setattr(self.state, name, value)

    def apply_user_input(self, user_input: str) -> None:
        value = _one_line(user_input)[:500]
        if not value:
            return
        confirmed_inputs = _unique([*self.state.confirmed_inputs, f"Latest user input: {value}"])
        self.update(
            current_goal=value,
            confirmed_inputs=confirmed_inputs,
            current_stage="resuming",
            blockers=[],
            next_action="Continue the task using the latest user input and retained verified state.",
        )

    def observe_tool_call(self, tool_name: str, args: dict[str, Any]) -> None:
        if tool_name == "list_workflow_skills":
            self.update(current_stage="selecting_skill", next_action="Inspect the best matching workflow skill.")
        elif tool_name == "inspect_workflow_skill":
            self.update(
                selected_skill=str(args.get("skill_id") or self.state.selected_skill),
                current_stage="inspecting_skill",
                next_action="Inspect the exact scripts and function signatures needed for the task.",
            )
        elif tool_name in SKILL_DETAIL_TOOLS:
            self.update(current_stage="inspecting_skill_details", next_action="Use verified signatures to prepare the analysis code.")
        elif tool_name in {"list_files", "read_file", "glob_search", "grep_text"}:
            self.update(current_stage="inspecting_files", next_action="Use the inspected input context to choose the next action.")
        elif tool_name in {"write_run_file", "edit_run_file", "validate_script"}:
            self.update(current_stage="preparing_script", next_action="Validate the persistent script, then start the selected runtime job.")
        elif tool_name == "inspect_image_catalog":
            self.update(current_stage="selecting_runtime", next_action="Choose the matching runtime environment.")
        elif tool_name in EXECUTION_TOOLS:
            runtime = str(args.get("env_profile") or args.get("runtime") or self.state.runtime_environment)
            skill_id = str(args.get("skill_id") or self.state.selected_skill)
            self.update(
                selected_skill=skill_id,
                runtime_environment=runtime,
                current_stage="executing",
                next_action="Wait for and inspect the execution result.",
            )
        elif tool_name == "request_user_input":
            self.update(current_stage="waiting_for_user", next_action="Wait for the requested user input.")
        elif tool_name == "inspect_artifact":
            self.update(current_stage="verifying_artifacts", next_action="Inspect all important outputs, then finish with their evidence ids.")
        elif tool_name == "finish_task":
            self.update(current_stage="finalizing", next_action="Return the evidence-grounded final answer.")
        elif tool_name in JOB_MANAGEMENT_TOOLS:
            self.update(current_stage="checking_job", next_action="Use the persisted job status and exact logs to choose the next action.")

    def observe_tool_result(self, tool_name: str, args: dict[str, Any], result: Any) -> None:
        if tool_name == "request_user_input" and isinstance(result, dict):
            question = _one_line(result.get("question") or result.get("reason") or "User input required")
            self.update(blockers=[question], current_stage="waiting_for_user", next_action="Wait for user input, then continue this task.")
            return
        if tool_name not in EXECUTION_TOOLS:
            if tool_name in JOB_MANAGEMENT_TOOLS and isinstance(result, dict):
                active_job_id = str(result.get("active_job_id") or result.get("job_id") or self.state.active_job_id)
                job_status = str(result.get("status") or self.state.job_status)
                self.update(active_job_id=active_job_id, job_status=job_status)
            if tool_name == "inspect_artifact" and _result_succeeded(result):
                artifacts = _unique([*self.state.generated_artifacts, *_extract_artifacts({"artifacts": [result.get("path")]})])
                self.update(generated_artifacts=artifacts, current_stage="artifact_verified")
            elif tool_name == "finish_task" and _result_succeeded(result):
                self.update(
                    active_errors=[],
                    current_stage="completed",
                    next_action="No further action required.",
                    blockers=[],
                )
            return

        if isinstance(result, dict):
            active_job_id = str(result.get("active_job_id") or result.get("job_id") or self.state.active_job_id)
            job_status = str(result.get("status") or self.state.job_status)
            self.update(active_job_id=active_job_id, job_status=job_status)

        if isinstance(result, dict) and result.get("status") == "running":
            self.update(current_stage="job_running", blockers=[], next_action="Poll the active job with an appropriate wait interval.")
            return

        if _result_succeeded(result):
            resolved = _unique([*self.state.resolved_errors, *self.state.active_errors])
            artifacts = _unique([*self.state.generated_artifacts, *_extract_artifacts(result)])
            self.update(
                active_errors=[],
                resolved_errors=resolved,
                generated_artifacts=artifacts,
                blockers=[],
                current_stage="execution_verified",
                next_action="Review verified outputs and report results.",
            )
            return

        error = _execution_error(result)
        self.update(
            active_errors=[error] if error else ["Execution failed without a structured error."],
            blockers=[error] if error else ["Execution failed without a structured error."],
            current_stage="repairing_error",
            next_action="Use the exact error and verified signatures to make the smallest repair.",
        )

    def mark_completed(self) -> None:
        self.update(current_stage="completed", next_action="No further action required.", blockers=[])

    def mark_failed(self) -> None:
        if self.state.active_errors:
            self.update(
                current_stage="repairing_error",
                next_action="Use the exact error and verified signatures to make the smallest repair.",
            )
            return
        self.update(current_stage="failed", next_action="Review the failed run and choose the next corrective action.")

    def mark_paused(self) -> None:
        self.update(current_stage="paused", next_action="Wait for user instructions before continuing.")

    def compact_summary(self, *, max_chars: int = 2000) -> str:
        state = self.state
        lines = [
            "TaskState (current task-local working state; this replaces older state):",
            f"- task: {_one_line(state.task)}",
            f"- current_goal: {_one_line(state.current_goal)}",
            f"- confirmed_inputs: {_format_list(state.confirmed_inputs)}",
            f"- expected_outputs: {_format_list(state.expected_outputs)}",
            f"- selected_skill: {state.selected_skill or '(not selected)'}",
            f"- runtime_environment: {state.runtime_environment or '(not selected)'}",
            f"- active_job_id: {state.active_job_id or '(none; call start_job before poll_job)'}",
            f"- job_status: {state.job_status or '(none)'}",
            f"- current_stage: {state.current_stage}",
            f"- active_errors: {_format_list(state.active_errors)}",
            f"- resolved_errors: {_format_list(state.resolved_errors)}",
            f"- generated_artifacts: {_format_list(state.generated_artifacts)}",
            f"- blockers: {_format_list(state.blockers)}",
            f"- next_action: {_one_line(state.next_action)}",
        ]
        return "\n".join(lines)[:max_chars]


def _result_succeeded(result: Any) -> bool:
    if not isinstance(result, dict):
        return False
    if "ok" in result:
        return bool(result.get("ok"))
    return not bool(result.get("error") or result.get("error_reason"))


def _execution_error(result: Any) -> str:
    if not isinstance(result, dict):
        return _one_line(result)
    value = result.get("error_reason") or result.get("error")
    if not value:
        stderr = str(result.get("stderr") or "").strip()
        value = stderr.splitlines()[-1] if stderr else ""
    return _one_line(value)[:500]


def _extract_artifacts(result: Any) -> list[str]:
    if not isinstance(result, dict):
        return []
    artifacts = result.get("artifacts") or []
    if not isinstance(artifacts, list):
        return []
    values: list[str] = []
    for artifact in artifacts:
        if isinstance(artifact, dict):
            value = artifact.get("path") or artifact.get("name")
        else:
            value = artifact
        if value:
            values.append(str(value))
    return values


def _format_list(values: list[str]) -> str:
    if not values:
        return "(none)"
    return "; ".join(_one_line(value)[:240] for value in values[-5:])


def _one_line(value: Any) -> str:
    return " ".join(str(value or "").split())


def _unique(values: list[str]) -> list[str]:
    return list(dict.fromkeys(value for value in values if value))
