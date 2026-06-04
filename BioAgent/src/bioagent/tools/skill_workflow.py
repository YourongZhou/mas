from __future__ import annotations

import json
import re
import time
from pathlib import Path
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
        "run_dir_host": str(run_dir),
        "skill_metadata": skill.get("metadata", {}),
        "context_policy": {
            "strategy": "compact_relevance_ranked_context",
            "script_previews": "Only selected scripts include signatures/previews. All script paths remain listed.",
            "reference_previews": "Only selected references include previews. Load behavior should be inferred conservatively from Skill instructions.",
        },
        "skill_body": skill_body[:18000],
        "available_script_paths": _relative_files(skill_path / "scripts", skill_path),
        "available_reference_paths": _relative_files(skill_path / "references", skill_path),
        "selected_script_context": scripts,
        "selected_reference_context": references,
    }
    if previous_result is not None:
        payload["repair_context"] = {
            "previous_code": previous_code[-20000:],
            "exit_code": previous_result.get("exit_code"),
            "stdout_tail": str(previous_result.get("stdout") or "")[-5000:],
            "stderr_tail": str(previous_result.get("stderr") or "")[-8000:],
        }
    return json.dumps(payload, ensure_ascii=False, indent=2, default=str)


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


def _script_inventory(skill_path: Path, *, query: str, skill_body: str, max_items: int = 12) -> list[dict[str, str]]:
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


def _reference_inventory(skill_path: Path, *, query: str, skill_body: str, max_items: int = 4) -> list[dict[str, str]]:
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
            "preview": body[:1800],
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
    if previous_result is not None:
        parts.append(str(previous_result.get("stderr") or "")[-2000:])
        parts.append(str(previous_result.get("stdout") or "")[-1200:])
    return "\n".join(parts)


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


def _compact_source_preview(text: str) -> str:
    lines = text.splitlines()
    selected: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith(("def ", "class ", "import ", "from ")):
            selected.append(line)
    preview = "\n".join(selected[:120])
    if not preview:
        preview = "\n".join(lines[:80])
    return preview[:5000]


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
        return {}
    try:
        parsed = json.loads(match.group(0))
    except Exception:
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _extract_fenced_code(text: str) -> str:
    match = re.search(r"```(?:python|r|R)?\s*(.*?)```", text, flags=re.S)
    return match.group(1).strip() if match else ""
