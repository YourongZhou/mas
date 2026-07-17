from __future__ import annotations

import json
import re
import time
from pathlib import Path
from pathlib import PurePosixPath
from typing import Any

from langchain_core.messages import HumanMessage, SystemMessage

from bioagent.config import AgentConfig
from bioagent.docker_runner import DockerRunner
from bioagent.llm import build_llm
from bioagent.logging_utils import RunLogger
from bioagent.skills import image_for_profile, inspect_workflow_skill
from bioagent.skills.registry import parse_frontmatter

from .common import resolve_data_path, to_repo_mount_path


GENERATOR_SYSTEM_PROMPT = """You generate executable bioinformatics workflow code from local workflow skills.

Return only a JSON object with these fields:
- runtime: "python" or "r"
- env_profile: Docker profile to use
- code: complete executable script
- requirements: optional pip requirements for Python only, usually empty because profiles are prebuilt
- rationale: short explanation of how the code follows the skill

Rules:
- The repository is mounted read-only at /repo.
- The run workspace is mounted at /work.
- Write every output under /work/outputs.
- Use the workflow skill as the primary source of truth.
- Prefer importing the skill's scripts from /repo/.../scripts when the skill says to use them.
- For Python, add the provided skill_script_dir_in_container to sys.path before importing skill scripts.
- The provided context is intentionally compact. Use selected script previews and the full available script paths; do not assume omitted previews mean omitted files are unavailable.
- Keep code self-contained enough to run as a single script.
- End the script by printing a machine-readable or human-readable line starting with ===RESULT===.
- Do not ask clarification questions; choose conservative defaults and record assumptions in outputs.
- Prioritize a robust core workflow over an over-broad workflow. For broad "standard analysis" requests, first complete load/QC/normalization/reduction/clustering/visualization/export.
- Do not run optional DE, pseudobulk, trajectory, communication, annotation, or large PDF-report steps unless the user explicitly asks and required metadata exists.
- Optional plots/reports must be wrapped in try/except so a non-critical visualization error does not fail the entire workflow after core outputs are saved.
- Save core outputs before optional sections: summary JSON/TXT, processed object when applicable, key tables, and at least one diagnostic plot.
"""


CODE_GENERATOR_SYSTEM_PROMPT = """You generate executable, sandbox-safe bioinformatics or data-analysis code.

Return only a JSON object with these fields:
- runtime: "python" or "r"
- env_profile: Docker profile to use
- code: complete executable script
- requirements: optional pip requirements for Python only, usually empty because profiles are prebuilt
- rationale: short explanation

Rules:
- Use this route when no local workflow skill is necessary or no matching skill exists.
- The repository is mounted read-only at /repo.
- The run workspace is mounted at /work.
- Write every output under /work/outputs.
- Keep the script focused on the requested task.
- Prefer standard libraries and packages available in the selected Docker profile.
- End the script by printing a line starting with ===RESULT===.
- Do not ask clarification questions; choose conservative defaults and record assumptions in outputs.
- Save useful artifacts such as summary JSON/TXT, CSV tables, or plots when appropriate.
- Optional plotting/report steps must be fail-soft and must not break the core result.
"""


def run_skill_workflow_impl(
    config: AgentConfig,
    logger: RunLogger,
    run_dir: Path,
    skill_id: str,
    task: str,
    data_path: str = "",
    runtime: str = "",
    env_profile: str = "",
    max_attempts: int = 5,
    timeout_s: int = 1800,
) -> dict[str, Any]:
    skill = inspect_workflow_skill(config.workflows_root, skill_id, max_chars=30000)
    if skill.get("error"):
        logger.error_reason(f"读取 Skill 失败：{skill['error']}")
        return {"ok": False, "error": skill["error"], "skill_id": skill_id}

    preflight_request = _preflight_user_input_request(skill=skill, skill_id=skill_id, task=task, data_path=data_path)
    if preflight_request:
        logger.progress("Skill 工作流暂停等待用户输入", preflight_request["question"])
        return preflight_request

    metadata = skill.get("metadata") if isinstance(skill.get("metadata"), dict) else {}
    runtime = (runtime or str(metadata.get("runtime") or "")).strip().lower()
    env_profile = (env_profile or str(metadata.get("env_profile") or "")).strip()
    runtime, env_profile = _resolve_runtime_and_profile(config, runtime, env_profile)

    host_data_path = None
    container_data_path = ""
    if data_path.strip():
        try:
            host_data_path = resolve_data_path(config, run_dir, data_path)
            container_data_path = to_repo_mount_path(config, host_data_path)
        except Exception as exc:
            logger.error_reason(f"解析输入数据路径失败：{exc}")
            return {"ok": False, "error": f"Could not resolve data_path: {exc}", "data_path": data_path}

    logger.node("Skill 工作流启动", f"skill_id={skill_id} runtime={runtime} env_profile={env_profile}")
    logger.progress("任务说明", task)
    if host_data_path:
        logger.progress("输入数据", f"host={host_data_path} container={container_data_path}")
    logger.progress("输出目录", f"container=/work/outputs host={run_dir / 'outputs'}")
    required_outputs = _required_output_paths_from_task(task)
    if required_outputs:
        logger.progress("要求输出文件", ", ".join(required_outputs))

    docker = DockerRunner(config=config, logger=logger, run_dir=run_dir)
    attempts: list[dict[str, Any]] = []
    previous_code = ""
    previous_result: dict[str, Any] | None = None

    for attempt in range(1, max(1, max_attempts) + 1):
        logger.node(
            f"Skill 工作流尝试 {attempt}/{max(1, max_attempts)}",
            "首次根据 Skill 生成代码。" if previous_result is None else "根据上一轮代码、报错和输出进行最小修复。",
        )
        if previous_result is None:
            logger.progress("生成代码", "首次生成，将携带 Skill 文档、相关脚本摘要、输入数据路径和 Docker 环境。")
        else:
            reason = _extract_error_reason(previous_result.get("stderr"), previous_result.get("stdout"))
            logger.progress("修复代码", f"previous_exit_code={previous_result.get('exit_code')}")
            logger.error_reason(reason)
        prompt = _build_generation_prompt(
            config=config,
            run_dir=run_dir,
            skill=skill,
            skill_id=skill_id,
            task=task,
            runtime=runtime,
            env_profile=env_profile,
            host_data_path=str(host_data_path) if host_data_path else "",
            container_data_path=container_data_path,
            previous_code=previous_code,
            previous_result=previous_result,
            attempt=attempt,
        )
        started = time.perf_counter()
        logger.progress("调用代码生成模型", f"attempt={attempt}")
        generated = _invoke_code_generator(config, prompt)
        generation_elapsed = time.perf_counter() - started
        logger.progress("代码生成完成", f"attempt={attempt} elapsed={generation_elapsed:.2f}s")
        logger.log(f"SKILL_CODE_GENERATION attempt={attempt} elapsed={generation_elapsed:.2f}s")
        logger.preview(
            f"SKILL_GENERATED_SUMMARY attempt={attempt}",
            json.dumps(_generated_log_summary(generated), ensure_ascii=False, indent=2),
            max_chars=5000,
        )

        code = str(generated.get("code") or "").strip()
        if not code:
            return {"ok": False, "error": "Code generator returned empty code.", "attempts": attempts}

        if _looks_like_malformed_generator_json(code):
            error_reason = "Code generator returned malformed JSON-like content instead of executable code."
            result = {
                "ok": False,
                "runtime": runtime,
                "env_profile": env_profile,
                "script_path": "",
                "exit_code": 2,
                "stdout": "",
                "stderr": error_reason,
                "error_reason": error_reason,
            }
            attempt_record = {
                "attempt": attempt,
                "ok": False,
                "runtime": runtime,
                "env_profile": env_profile,
                "script_path": "",
                "exit_code": 2,
                "error_reason": error_reason,
                "stdout_tail": "",
                "stderr_tail": error_reason,
                "rationale": generated.get("rationale", ""),
            }
            attempts.append(attempt_record)
            logger.preview(f"SKILL_ATTEMPT_RESULT attempt={attempt}", json.dumps(attempt_record, ensure_ascii=False, indent=2), max_chars=8000)
            previous_code = code[-4000:]
            previous_result = result
            logger.progress("本次生成无效", f"attempt={attempt} will_retry={attempt < max(1, max_attempts)}")
            logger.error_reason(error_reason)
            continue

        generated_runtime = str(generated.get("runtime") or runtime).strip().lower()
        generated_profile = str(generated.get("env_profile") or env_profile).strip()
        runtime, env_profile = _resolve_runtime_and_profile(config, generated_runtime, generated_profile)

        if runtime == "r":
            logger.progress("执行生成代码", f"runtime=r env_profile={env_profile}")
            result = docker.execute_r(code, env_profile=env_profile, timeout_s=timeout_s).to_dict()
        else:
            logger.progress("执行生成代码", f"runtime=python env_profile={env_profile}")
            result = docker.execute_python(
                code,
                env_profile=env_profile,
                requirements=str(generated.get("requirements") or ""),
                timeout_s=timeout_s,
            ).to_dict()

        missing_outputs = _missing_required_outputs(run_dir, required_outputs) if result.get("ok") else []
        if missing_outputs:
            message = "Missing required output files after successful exit: " + ", ".join(missing_outputs)
            result = dict(result)
            result["ok"] = False
            result["error_reason"] = message
            result["stderr"] = (str(result.get("stderr") or "").rstrip() + "\n" + message).strip()
            logger.error_reason(message)

        error_reason = str(result.get("error_reason") or _extract_error_reason(result.get("stderr"), result.get("stdout")))
        attempt_record = {
            "attempt": attempt,
            "ok": bool(result.get("ok")),
            "runtime": runtime,
            "env_profile": env_profile,
            "script_path": result.get("script_path"),
            "exit_code": result.get("exit_code"),
            "error_reason": "" if result.get("ok") else error_reason,
            "stdout_tail": str(result.get("stdout") or "")[-2000:],
            "stderr_tail": str(result.get("stderr") or "")[-4000:],
            "rationale": generated.get("rationale", ""),
        }
        attempts.append(attempt_record)
        logger.preview(f"SKILL_ATTEMPT_RESULT attempt={attempt}", json.dumps(attempt_record, ensure_ascii=False, indent=2), max_chars=8000)

        if result.get("ok"):
            logger.progress("Skill 工作流完成", f"attempts={attempt} output_dir={run_dir / 'outputs'}")
            output_dir = str(run_dir / "outputs")
            return {
                "ok": True,
                "skill_id": skill_id,
                "runtime": runtime,
                "env_profile": env_profile,
                "host_data_path": str(host_data_path) if host_data_path else "",
                "container_data_path": container_data_path,
                "output_dir": output_dir,
                "attempts": attempts,
                "final_code": code,
                "execution": result,
            }

        previous_code = code
        previous_result = result
        logger.progress("本次执行失败", f"attempt={attempt} exit_code={result.get('exit_code')} will_retry={attempt < max(1, max_attempts)}")
        logger.error_reason(error_reason)

    logger.progress("Skill 工作流失败终止", f"attempts={len(attempts)} 已达到最大尝试次数。")
    return {
        "ok": False,
        "skill_id": skill_id,
        "runtime": runtime,
        "env_profile": env_profile,
        "host_data_path": str(host_data_path) if host_data_path else "",
        "container_data_path": container_data_path,
        "attempts": attempts,
        "final_code": previous_code,
        "last_execution": previous_result,
        "error": f"Skill workflow failed after {len(attempts)} attempt(s).",
    }


def run_code_workflow_impl(
    config: AgentConfig,
    logger: RunLogger,
    run_dir: Path,
    task: str,
    data_path: str = "",
    runtime: str = "python",
    env_profile: str = "",
    max_attempts: int = 5,
    timeout_s: int = 900,
) -> dict[str, Any]:
    runtime, env_profile = _resolve_runtime_and_profile(config, runtime, env_profile)

    host_data_path = None
    container_data_path = ""
    if data_path.strip():
        try:
            host_data_path = resolve_data_path(config, run_dir, data_path)
            container_data_path = to_repo_mount_path(config, host_data_path)
        except Exception as exc:
            logger.error_reason(f"解析输入数据路径失败：{exc}")
            return {"ok": False, "error": f"Could not resolve data_path: {exc}", "data_path": data_path}

    logger.node("通用代码工作流启动", f"runtime={runtime} env_profile={env_profile}")
    logger.progress("任务说明", task)
    if host_data_path:
        logger.progress("输入数据", f"host={host_data_path} container={container_data_path}")
    logger.progress("输出目录", f"container=/work/outputs host={run_dir / 'outputs'}")

    docker = DockerRunner(config=config, logger=logger, run_dir=run_dir)
    attempts: list[dict[str, Any]] = []
    previous_code = ""
    previous_result: dict[str, Any] | None = None

    for attempt in range(1, max(1, max_attempts) + 1):
        logger.node(
            f"通用代码工作流尝试 {attempt}/{max(1, max_attempts)}",
            "首次按任务直接生成代码。" if previous_result is None else "根据上一轮代码、报错和输出进行最小修复。",
        )
        if previous_result is None:
            logger.progress("生成代码", "首次生成，将携带任务、输入数据路径、运行目录和 Docker 环境。")
        else:
            reason = _extract_error_reason(previous_result.get("stderr"), previous_result.get("stdout"))
            logger.progress("修复代码", f"previous_exit_code={previous_result.get('exit_code')}")
            logger.error_reason(reason)
        prompt = _build_code_generation_prompt(
            config=config,
            run_dir=run_dir,
            task=task,
            runtime=runtime,
            env_profile=env_profile,
            host_data_path=str(host_data_path) if host_data_path else "",
            container_data_path=container_data_path,
            previous_code=previous_code,
            previous_result=previous_result,
            attempt=attempt,
        )
        started = time.perf_counter()
        logger.progress("调用代码生成模型", f"attempt={attempt}")
        generated = _invoke_code_generator(config, prompt, system_prompt=CODE_GENERATOR_SYSTEM_PROMPT)
        generation_elapsed = time.perf_counter() - started
        logger.progress("代码生成完成", f"attempt={attempt} elapsed={generation_elapsed:.2f}s")
        logger.log(f"CODE_GENERATION attempt={attempt} elapsed={generation_elapsed:.2f}s")
        logger.preview(
            f"CODE_GENERATED_SUMMARY attempt={attempt}",
            json.dumps(_generated_log_summary(generated), ensure_ascii=False, indent=2),
            max_chars=5000,
        )

        code = str(generated.get("code") or "").strip()
        if not code:
            return {"ok": False, "error": "Code generator returned empty code.", "attempts": attempts}

        runtime, env_profile = _resolve_runtime_and_profile(
            config,
            str(generated.get("runtime") or runtime).strip().lower(),
            str(generated.get("env_profile") or env_profile).strip(),
        )
        if runtime == "r":
            logger.progress("执行生成代码", f"runtime=r env_profile={env_profile}")
            result = docker.execute_r(code, env_profile=env_profile, timeout_s=timeout_s).to_dict()
        else:
            logger.progress("执行生成代码", f"runtime=python env_profile={env_profile}")
            result = docker.execute_python(
                code,
                env_profile=env_profile,
                requirements=str(generated.get("requirements") or ""),
                timeout_s=timeout_s,
            ).to_dict()

        error_reason = str(result.get("error_reason") or _extract_error_reason(result.get("stderr"), result.get("stdout")))
        attempt_record = {
            "attempt": attempt,
            "ok": bool(result.get("ok")),
            "runtime": runtime,
            "env_profile": env_profile,
            "script_path": result.get("script_path"),
            "exit_code": result.get("exit_code"),
            "error_reason": "" if result.get("ok") else error_reason,
            "stdout_tail": str(result.get("stdout") or "")[-2000:],
            "stderr_tail": str(result.get("stderr") or "")[-4000:],
            "rationale": generated.get("rationale", ""),
        }
        attempts.append(attempt_record)
        logger.preview(f"CODE_ATTEMPT_RESULT attempt={attempt}", json.dumps(attempt_record, ensure_ascii=False, indent=2), max_chars=8000)

        if result.get("ok"):
            logger.progress("通用代码工作流完成", f"attempts={attempt} output_dir={run_dir / 'outputs'}")
            return {
                "ok": True,
                "workflow_type": "code",
                "runtime": runtime,
                "env_profile": env_profile,
                "host_data_path": str(host_data_path) if host_data_path else "",
                "container_data_path": container_data_path,
                "output_dir": str(run_dir / "outputs"),
                "attempts": attempts,
                "final_code": code,
                "execution": result,
            }

        previous_code = code
        previous_result = result
        logger.progress("本次执行失败", f"attempt={attempt} exit_code={result.get('exit_code')} will_retry={attempt < max(1, max_attempts)}")
        logger.error_reason(error_reason)

    logger.progress("通用代码工作流失败终止", f"attempts={len(attempts)} 已达到最大尝试次数。")
    return {
        "ok": False,
        "workflow_type": "code",
        "runtime": runtime,
        "env_profile": env_profile,
        "host_data_path": str(host_data_path) if host_data_path else "",
        "container_data_path": container_data_path,
        "attempts": attempts,
        "final_code": previous_code,
        "last_execution": previous_result,
        "error": f"Code workflow failed after {len(attempts)} attempt(s).",
    }


def _resolve_runtime_and_profile(config: AgentConfig, runtime: str, env_profile: str) -> tuple[str, str]:
    runtime = (runtime or "").strip().lower()
    env_profile = (env_profile or "").strip()
    if not runtime:
        runtime = "python"
    if runtime == "mixed":
        runtime = "python"
    if not env_profile:
        env_profile = "py-general-v1" if runtime == "python" else "r-bioc-v1"
    entry = image_for_profile(config.docker_root, env_profile)
    if entry:
        entry_runtime = str(entry.get("runtime") or "").strip().lower()
        if entry_runtime in {"python", "r"}:
            runtime = entry_runtime
    return runtime, env_profile


def _generated_log_summary(generated: dict[str, Any]) -> dict[str, Any]:
    code = str(generated.get("code") or "")
    lines = code.splitlines()
    return {
        "runtime": str(generated.get("runtime") or ""),
        "env_profile": str(generated.get("env_profile") or ""),
        "requirements": str(generated.get("requirements") or "").strip(),
        "code_lines": len(lines),
        "rationale": str(generated.get("rationale") or ""),
        "code_preview_head": "\n".join(lines[:40]),
        "code_preview_tail": "\n".join(lines[-20:]) if len(lines) > 60 else "",
        "trace_note": "完整代码已写入本次 attempt 的 script_path，日志只保留首尾摘要以便阅读。",
    }


def _extract_error_reason(stderr: Any, stdout: Any = "") -> str:
    text = str(stderr or "").strip() or str(stdout or "").strip()
    if not text:
        return "进程返回失败，但 stdout/stderr 为空。"
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    context = ""
    for line in reversed(lines):
        if not context and line.startswith("Error raised while"):
            context = line
        if _looks_like_exception_line(line):
            return f"{line}；{context}" if context and context != line else line
    for line in reversed(lines):
        lowered = line.lower()
        if "error" in lowered or "exception" in lowered or "timed out" in lowered:
            return line
    return lines[-1][-500:] if lines else "未能自动提取明确错误原因。"


def _looks_like_exception_line(line: str) -> bool:
    if line.startswith(("Traceback ", "During handling ")):
        return False
    return bool(
        line.startswith(("Error: ", "Exception: "))
        or "Error:" in line[:120]
        or line.startswith(
            (
                "AssertionError:",
                "ImportError:",
                "KeyError:",
                "ModuleNotFoundError:",
                "RuntimeError:",
                "TypeError:",
                "ValueError:",
            )
        )
    )


def _looks_like_malformed_generator_json(code: str) -> bool:
    stripped = code.lstrip()
    if not stripped.startswith("{"):
        return False
    return '"code"' in stripped[:2000] or "'code'" in stripped[:2000]


def _build_generation_prompt(
    *,
    config: AgentConfig,
    run_dir: Path,
    skill: dict[str, Any],
    skill_id: str,
    task: str,
    runtime: str,
    env_profile: str,
    host_data_path: str,
    container_data_path: str,
    previous_code: str,
    previous_result: dict[str, Any] | None,
    attempt: int,
) -> str:
    skill_path = Path(str(skill.get("path") or (config.workflows_root / skill_id))).resolve()
    skill_body = str(skill.get("body_preview") or "")
    context_query = _context_query(task, previous_result)
    scripts = _script_inventory(skill_path, query=context_query, skill_body=skill_body)
    references = _reference_inventory(skill_path, query=context_query, skill_body=skill_body)
    script_dir = "/repo/" + skill_path.relative_to(config.repo_root).as_posix() + "/scripts"
    payload = {
        "attempt": attempt,
        "task": task,
        "skill_id": skill_id,
        "skill_path_host": str(skill_path),
        "skill_script_dir_in_container": script_dir,
        "runtime": runtime,
        "env_profile": env_profile,
        "host_data_path": host_data_path,
        "container_data_path": container_data_path,
        "repo_mount": "/repo",
        "work_mount": "/work",
        "output_dir": "/work/outputs",
        "required_output_paths": _required_output_paths_from_task(task),
        "run_dir_host": str(run_dir),
        "skill_metadata": skill.get("metadata", {}),
        "context_policy": {
            "strategy": "compact_relevance_ranked_context",
            "script_previews": "Only selected scripts include signatures/previews. All script paths remain listed.",
            "reference_previews": "Only selected references include previews. Load behavior should be inferred conservatively from Skill instructions.",
        },
        "skill_body": skill_body[:10000],
        "available_script_paths": _relative_files(skill_path / "scripts", skill_path),
        "available_reference_paths": _relative_files(skill_path / "references", skill_path),
        "selected_script_context": scripts,
        "selected_reference_context": references,
    }
    if previous_result is not None:
        payload["repair_context"] = {
            "previous_code": previous_code[-8000:],
            "error_reason": _extract_error_reason(previous_result.get("stderr"), previous_result.get("stdout")),
            "exit_code": previous_result.get("exit_code"),
            "stdout_tail": str(previous_result.get("stdout") or "")[-2000:],
            "stderr_tail": str(previous_result.get("stderr") or "")[-3000:],
        }
    return json.dumps(payload, ensure_ascii=False, indent=2, default=str)


def _preflight_user_input_request(
    *,
    skill: dict[str, Any],
    skill_id: str,
    task: str,
    data_path: str,
) -> dict[str, Any] | None:
    if data_path.strip():
        return None
    lowered_task = task.lower()
    if any(token in lowered_task for token in ("example", "demo", "示例", "演示")):
        return None

    body = str(skill.get("body_preview") or "")
    body_lower = body.lower()
    requires_first_question = "always ask question 1 first" in body_lower or "ask this first" in body_lower
    asks_for_input_files = "input files" in body_lower or "输入文件" in body_lower
    if not (requires_first_question and asks_for_input_files):
        return None

    question = _extract_first_clarification_question(body) or (
        "Which input data file should this workflow analyze? "
        "Please provide a repo-relative or absolute path, or say that you want to use the example/demo data."
    )
    return {
        "status": "needs_user_input",
        "skill_id": skill_id,
        "question": question,
        "reason": "The selected workflow skill explicitly requires the input data choice before execution.",
        "required_fields": ["data_path"],
        "resume_hint": "Use the user's answer as data_path if it is a path; if they request demo/example data, rerun with the skill's documented example defaults.",
    }


def _extract_first_clarification_question(body: str) -> str:
    lines = [line.strip(" -*\t") for line in body.splitlines()]
    in_questions = False
    for line in lines:
        lowered = line.lower()
        if "clarification questions" in lowered:
            in_questions = True
            continue
        if not in_questions:
            continue
        if line.endswith("?") or line.endswith("？"):
            return line
        if "do you have specific" in lowered and "data" in lowered:
            return line.rstrip(".") + "?"
    return ""


def _build_code_generation_prompt(
    *,
    config: AgentConfig,
    run_dir: Path,
    task: str,
    runtime: str,
    env_profile: str,
    host_data_path: str,
    container_data_path: str,
    previous_code: str,
    previous_result: dict[str, Any] | None,
    attempt: int,
) -> str:
    payload: dict[str, Any] = {
        "attempt": attempt,
        "task": task,
        "runtime": runtime,
        "env_profile": env_profile,
        "host_data_path": host_data_path,
        "container_data_path": container_data_path,
        "repo_root_host": str(config.repo_root),
        "repo_mount": "/repo",
        "work_mount": "/work",
        "output_dir": "/work/outputs",
        "run_dir_host": str(run_dir),
        "available_profiles_hint": [
            "py-general-v1: general Python, pandas/numpy/scipy/sklearn/matplotlib/reportlab",
            "py-scverse-v1: Python single-cell/scverse/scanpy",
            "r-bioc-v1: R Bioconductor/DESeq2/survival",
            "r-singlecell-v1: R Seurat/single-cell",
            "mixed-genetics-v1: mixed Python/R genetics tooling",
        ],
    }
    if previous_result is not None:
        payload["repair_context"] = {
            "previous_code": previous_code[-20000:],
            "exit_code": previous_result.get("exit_code"),
            "stdout_tail": str(previous_result.get("stdout") or "")[-5000:],
            "stderr_tail": str(previous_result.get("stderr") or "")[-8000:],
        }
    return json.dumps(payload, ensure_ascii=False, indent=2, default=str)


def _script_inventory(skill_path: Path, *, query: str, skill_body: str, max_items: int = 10) -> list[dict[str, str]]:
    scripts_dir = skill_path / "scripts"
    if not scripts_dir.is_dir():
        return []
    scored: list[tuple[int, Path, str]] = []
    for path in sorted(scripts_dir.rglob("*")):
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        score = _context_score(path, text, query)
        rel = str(path.relative_to(skill_path)).replace("\\", "/")
        if rel in skill_body:
            score += 3
        scored.append((score, path, text))
    selected = sorted(scored, key=lambda item: (-item[0], str(item[1])))[:max_items]
    return [
        {
            "path": str(path.relative_to(skill_path)).replace("\\", "/"),
            "preview": _compact_source_preview(text),
            "selection_score": str(score),
        }
        for score, path, text in selected
    ]


def _reference_inventory(skill_path: Path, *, query: str, skill_body: str, max_items: int = 2) -> list[dict[str, str]]:
    references_dir = skill_path / "references"
    if not references_dir.is_dir():
        return []
    scored: list[tuple[int, Path, str]] = []
    for path in sorted(references_dir.rglob("*")):
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        _meta, body = parse_frontmatter(text)
        score = _context_score(path, body, query)
        rel = str(path.relative_to(skill_path)).replace("\\", "/")
        if rel in skill_body:
            score += 2
        scored.append((score, path, body))
    selected = sorted(scored, key=lambda item: (-item[0], str(item[1])))[:max_items]
    return [
        {
            "path": str(path.relative_to(skill_path)).replace("\\", "/"),
            "preview": body[:1200],
            "selection_score": str(score),
        }
        for score, path, body in selected
    ]


def _relative_files(root: Path, base: Path) -> list[str]:
    if not root.is_dir():
        return []
    return [str(path.relative_to(base)).replace("\\", "/") for path in sorted(root.rglob("*")) if path.is_file()]


def _context_query(task: str, previous_result: dict[str, Any] | None) -> str:
    parts = [task]
    parts.extend(_context_query_expansions(task))
    if previous_result is not None:
        parts.append(str(previous_result.get("stderr") or "")[-2000:])
        parts.append(str(previous_result.get("stdout") or "")[-1200:])
    return "\n".join(parts)


def _context_query_expansions(task: str) -> list[str]:
    lowered = task.lower()
    expansions: list[str] = []
    if any(token in lowered for token in ("qc", "quality", "质量")):
        expansions.append("qc quality metrics filter cells mitochondrial ribosomal violin")
    if any(token in lowered for token in ("hvg", "variable", "高变")):
        expansions.append("hvg highly variable genes find_variable_genes")
    if any(token in lowered for token in ("normalize", "normalization", "标准化")):
        expansions.append("normalize normalization log1p normalize_data")
    if any(token in lowered for token in ("pca", "principal")):
        expansions.append("pca scale scale_and_pca principal components")
    if any(token in lowered for token in ("neighbor", "neighbors", "邻居", "leiden", "cluster", "聚类")):
        expansions.append("neighbors leiden cluster cluster_cells")
    if any(token in lowered for token in ("umap", "降维")):
        expansions.append("umap dimreduction run_umap plot_dimreduction")
    if any(token in lowered for token in ("plot", "figure", "visual", "图", "可视化")):
        expansions.append("plot figure visualization violin scatter umap plot_qc plot_dimreduction")
    if any(token in lowered for token in ("summary", "json", "processed", "h5ad", "输出", "export")):
        expansions.append("summary json processed h5ad output export_results")
    return expansions


def _context_score(path: Path, text: str, query: str) -> int:
    haystack = (str(path.name) + "\n" + str(path) + "\n" + text[:4000]).lower()
    tokens = _query_tokens(query)
    score = sum(1 for token in tokens if token in haystack)
    stem = path.stem.lower().replace("_", " ")
    score += sum(2 for token in tokens if token in stem)
    return score


def _query_tokens(text: str) -> list[str]:
    raw = re.findall(r"[A-Za-z][A-Za-z0-9_+-]{2,}|[\u4e00-\u9fff]{2,}", text.lower())
    stop = {
        "the",
        "and",
        "for",
        "with",
        "from",
        "under",
        "output",
        "outputs",
        "analysis",
        "workflow",
        "数据",
        "分析",
        "输出",
    }
    seen: set[str] = set()
    tokens: list[str] = []
    for token in raw:
        if token in stop or token in seen:
            continue
        seen.add(token)
        tokens.append(token)
    return tokens[:80]


def _compact_source_preview(text: str, max_chars: int = 2500) -> str:
    lines = text.splitlines()
    selected: list[str] = []
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        if stripped.startswith(("import ", "from ")):
            selected.append(line)
            i += 1
            continue
        if stripped.startswith(("def ", "class ")):
            block: list[str] = [line]
            paren_balance = line.count("(") - line.count(")")
            i += 1
            while i < len(lines) and paren_balance > 0:
                block.append(lines[i])
                paren_balance += lines[i].count("(") - lines[i].count(")")
                i += 1
            while i < len(lines) and len(block) < 18:
                next_stripped = lines[i].strip()
                if not next_stripped:
                    block.append(lines[i])
                    i += 1
                    continue
                if next_stripped.startswith(('"""', "'''")):
                    block.append(lines[i])
                    quote = next_stripped[:3]
                    i += 1
                    while i < len(lines) and quote not in lines[i] and len(block) < 18:
                        block.append(lines[i])
                        i += 1
                    if i < len(lines) and len(block) < 18:
                        block.append(lines[i])
                        i += 1
                break
            selected.extend(block)
            continue
        i += 1
    preview = "\n".join(selected[:120])
    if not preview:
        preview = "\n".join(lines[:80])
    return preview[:max_chars]


def _required_output_paths_from_task(task: str) -> list[str]:
    seen: set[str] = set()
    paths: list[str] = []
    for match in re.findall(r"/work/outputs/[^\s`'\"<>\])},;]+", task):
        path = match.rstrip(".:")
        if path in {"/work/outputs", "/work/outputs/"}:
            continue
        rel = path.removeprefix("/work/outputs/").strip("/")
        if not rel or _unsafe_relative_output_path(rel):
            continue
        if path not in seen:
            seen.add(path)
            paths.append(path)
    if "/work/outputs" in task:
        output_section = task.split("/work/outputs", 1)[1]
        for filename in re.findall(
            r"(?<![/\w.-])([A-Za-z0-9][A-Za-z0-9_.-]*\.(?:h5ad|h5|json|png|svg|pdf|csv|tsv|txt|xlsx|rds|mtx))(?=$|[\s,;:.)\]}。；，、])",
            output_section,
            flags=re.I,
        ):
            path = f"/work/outputs/{filename.rstrip('.:')}"
            rel = path.removeprefix("/work/outputs/").strip("/")
            if not rel or _unsafe_relative_output_path(rel):
                continue
            if path not in seen:
                seen.add(path)
                paths.append(path)
    return paths


def _missing_required_outputs(run_dir: Path, required_outputs: list[str]) -> list[str]:
    missing: list[str] = []
    for output_path in required_outputs:
        rel = output_path.removeprefix("/work/outputs/").strip("/")
        if not rel or _unsafe_relative_output_path(rel):
            continue
        host_path = run_dir / "outputs" / rel
        if not host_path.exists() or (host_path.is_file() and host_path.stat().st_size == 0):
            missing.append(output_path)
    return missing


def _unsafe_relative_output_path(path: str) -> bool:
    return PurePosixPath(path).is_absolute() or ".." in PurePosixPath(path).parts


def _invoke_code_generator(config: AgentConfig, prompt: str, *, system_prompt: str = GENERATOR_SYSTEM_PROMPT) -> dict[str, Any]:
    llm = build_llm(config)
    response = llm.invoke(
        [
            SystemMessage(content=system_prompt),
            HumanMessage(content=prompt),
        ]
    )
    content = str(response.content or "")
    parsed = _parse_json_object(content)
    if parsed:
        return parsed
    code = _extract_fenced_code(content)
    return {
        "runtime": "python",
        "env_profile": "",
        "code": code or content,
        "requirements": "",
        "rationale": "Generator did not return JSON; used raw/fenced content as code.",
    }


def _parse_json_object(text: str) -> dict[str, Any]:
    text = text.strip()
    if text.startswith("```"):
        text = re.sub(r"^```(?:json)?\s*", "", text)
        text = re.sub(r"\s*```$", "", text)
    try:
        parsed = json.loads(text)
        return parsed if isinstance(parsed, dict) else {}
    except Exception:
        pass
    match = re.search(r"\{.*\}", text, flags=re.S)
    if not match:
        return _parse_json_like_fields(text)
    try:
        parsed = json.loads(match.group(0))
    except Exception:
        return _parse_json_like_fields(text)
    return parsed if isinstance(parsed, dict) else {}


def _parse_json_like_fields(text: str) -> dict[str, Any]:
    code = _extract_json_string_field(text, "code")
    if not code:
        return {}
    return {
        "runtime": _extract_json_string_field(text, "runtime") or "python",
        "env_profile": _extract_json_string_field(text, "env_profile"),
        "code": code,
        "requirements": _extract_json_string_field(text, "requirements"),
        "rationale": _extract_json_string_field(text, "rationale") or "Recovered code from a malformed JSON-like generator response.",
    }


def _extract_json_string_field(text: str, field: str) -> str:
    match = re.search(rf'"{re.escape(field)}"\s*:\s*"((?:\\.|[^"\\])*)"', text, flags=re.S)
    if not match:
        return ""
    try:
        value = json.loads('"' + match.group(1) + '"')
    except Exception:
        return ""
    return value if isinstance(value, str) else ""


def _extract_fenced_code(text: str) -> str:
    match = re.search(r"```(?:python|r|R)?\s*(.*?)```", text, flags=re.S)
    return match.group(1).strip() if match else ""
