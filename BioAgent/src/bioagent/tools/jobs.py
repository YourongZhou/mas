from __future__ import annotations

import json
import hashlib
import mimetypes
import re
import shlex
import shutil
import subprocess
import threading
import time
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from bioagent.config import AgentConfig
from bioagent.logging_utils import RunLogger
from bioagent.skills import image_for_profile

from .attempts import ExecutionAttemptBudget
from .workspace import resolve_run_path, run_relative_path, validate_script_impl


JOB_ID_PATTERN = re.compile(r"^job_\d{8}_\d{6}_[0-9a-f]{6}$")
TERMINAL_JOB_STATUSES = {"completed", "failed", "cancelled", "timed_out", "failed_to_start"}
JOB_METADATA_LOCK = threading.RLock()


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class DockerJobManager:
    def __init__(
        self,
        config: AgentConfig,
        logger: RunLogger | None,
        run_dir: Path,
        *,
        prior_run_dirs: list[Path] | None = None,
    ) -> None:
        self.config = config
        self.logger = logger
        self.run_dir = run_dir
        self.run_dirs = _session_run_dirs(config.runs_dir, run_dir, prior_run_dirs or [])
        self.jobs_dir = run_dir / "jobs"
        self.job_state_path = run_dir / "state" / "job_state.json"
        self.budget = ExecutionAttemptBudget(run_dir, config.max_execution_attempts)

    def start_job(
        self,
        *,
        runtime: str,
        script_path: str,
        env_profile: str,
        requirements: str = "",
        timeout_s: int = 1800,
    ) -> dict:
        active = self._active_running_job()
        if active is not None:
            return {
                "ok": False,
                "status": "job_already_running",
                "error": "A job is already running in this run. Poll or cancel it before starting another job.",
                "active_job_id": active["job_id"],
                "job": active,
            }
        normalized_runtime = runtime.strip().lower()
        if normalized_runtime not in {"python", "r"}:
            return {"ok": False, "error": f"Unsupported runtime: {runtime}. Use python or r."}
        script = resolve_run_path(self.run_dir, script_path)
        validation = validate_script_impl(self.run_dir, path=str(script), runtime=normalized_runtime)
        if not validation.get("ok"):
            return validation
        image_entry = image_for_profile(self.config.docker_root, env_profile)
        if not image_entry:
            return {"ok": False, "error": f"Unknown Docker env_profile: {env_profile}."}
        image = str(image_entry.get("image_tag") or "").strip()
        profile_runtime = str(image_entry.get("runtime") or "").strip()
        if profile_runtime and profile_runtime not in {normalized_runtime, "mixed"}:
            return {
                "ok": False,
                "error": f"Profile {env_profile} is runtime={profile_runtime}, not {normalized_runtime}.",
            }
        if shutil.which("docker") is None:
            return {"ok": False, "error": "docker CLI not found."}
        image_check = _run(["docker", "image", "inspect", image], timeout=30)
        if image_check.returncode != 0:
            return {"ok": False, "error": f"Docker image not found locally: {image}.", "stderr": image_check.stderr}

        attempt = self.budget.consume(tool="start_job")
        if not attempt.get("ok"):
            return attempt
        job_id = f"job_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:6]}"
        (self.run_dir / "state" / "final_verification.json").unlink(missing_ok=True)
        job_dir = self.jobs_dir / job_id
        job_dir.mkdir(parents=True, exist_ok=False)
        cidfile = job_dir / "container.cid"
        relative_script = run_relative_path(self.run_dir, script)
        script_sha256 = _sha256(script)
        docker_script = f"/work/{relative_script}"
        inner = self._inner_command(
            runtime=normalized_runtime,
            docker_script=docker_script,
            requirements=requirements,
            job_dir=job_dir,
            timeout_s=max(1, timeout_s),
        )
        command = [
            "docker",
            "run",
            "-d",
            "--cidfile",
            str(cidfile),
            "--label",
            f"bioagent.run_id={self.run_dir.name}",
            "--label",
            f"bioagent.job_id={job_id}",
            "-v",
            f"{self.config.repo_root}:/repo:ro",
            "-v",
            f"{self.run_dir}:/work",
            "-w",
            "/work",
            image,
            "sh",
            "-lc",
            inner,
        ]
        self._progress("Docker job starting", f"job_id={job_id} profile={env_profile} script={script}")
        started = _run(command, timeout=60)
        if started.returncode != 0:
            self._write_job(
                job_dir,
                {
                    "job_id": job_id,
                    "status": "failed_to_start",
                    "runtime": normalized_runtime,
                    "env_profile": env_profile,
                    "image": image,
                    "script_path": str(script),
                    "exit_code": started.returncode,
                    "stderr": started.stderr,
                    "created_at": _now(),
                    "callback_state": "pending",
                },
            )
            self._set_active_job(job_id)
            return {
                "ok": False,
                "status": "failed_to_start",
                "job_id": job_id,
                "exit_code": started.returncode,
                "error": started.stderr.strip() or "docker run failed",
                **attempt,
            }
        container_id = started.stdout.strip()
        cidfile.write_text(container_id + "\n", encoding="utf-8")
        metadata = {
            "job_id": job_id,
            "run_id": self.run_dir.name,
            "run_dir": str(self.run_dir),
            "status": "running",
            "runtime": normalized_runtime,
            "env_profile": env_profile,
            "image": image,
            "script_path": str(script),
            "script_sha256": script_sha256,
            "workspace_script_path": docker_script,
            "container_id": container_id,
            "timeout_s": max(1, timeout_s),
            "created_at": _now(),
            "updated_at": _now(),
            "attempts_used": attempt["attempts_used"],
            "callback_state": "pending",
        }
        self._write_job(job_dir, metadata)
        metadata["execution_manifest_path"] = str(_write_execution_manifest(self.run_dir, metadata))
        self._write_job(job_dir, metadata)
        self._set_active_job(job_id)
        return {"ok": True, **metadata, "attempts_remaining": attempt["attempts_remaining"]}

    def list_jobs(self, *, status: str = "") -> dict:
        jobs: list[dict[str, Any]] = []
        for run_dir, path in self._job_paths():
            try:
                metadata = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if status and metadata.get("status") != status:
                continue
            jobs.append(_compact_job(metadata, run_dir=run_dir))
        jobs.sort(key=lambda item: str(item.get("created_at") or ""), reverse=True)
        return {
            "ok": True,
            "active_job_id": self._running_job_id() or self.active_job_id(),
            "job_count": len(jobs),
            "jobs": jobs,
        }

    def get_job(self, job_id: str = "", *, refresh: bool = False) -> dict:
        loaded = self._load_job_or_error(job_id)
        if isinstance(loaded, dict):
            return loaded
        job_dir, metadata = loaded
        if refresh and metadata.get("status") == "running":
            return self._poll_once(job_dir, metadata)
        return self._public_result(metadata)

    def poll_job(self, job_id: str = "", *, wait_s: int = 0) -> dict:
        loaded = self._load_job_or_error(job_id)
        if isinstance(loaded, dict):
            return loaded
        job_dir, metadata = loaded
        deadline = time.monotonic() + max(0, min(wait_s, 300))
        while True:
            result = self._poll_once(job_dir, metadata)
            if result.get("status") != "running" or time.monotonic() >= deadline:
                return result
            time.sleep(min(2.0, max(0.1, deadline - time.monotonic())))

    def tail_job(self, job_id: str = "", *, lines: int = 200) -> dict:
        loaded = self._load_job_or_error(job_id)
        if isinstance(loaded, dict):
            return loaded
        job_dir, metadata = loaded
        resolved_job_id = str(metadata["job_id"])
        if metadata.get("status") in TERMINAL_JOB_STATUSES:
            return {
                "ok": True,
                "job_ok": metadata.get("status") == "completed",
                "job_id": resolved_job_id,
                "status": metadata.get("status"),
                "stdout_tail": str(metadata.get("stdout_tail") or "")[-12000:],
                "stderr_tail": str(metadata.get("stderr_tail") or metadata.get("stderr") or "")[-4000:],
                "error_reason": metadata.get("error") or "",
            }
        container_id = str(metadata.get("container_id") or "")
        logs = _run(["docker", "logs", "--tail", str(max(1, min(lines, 2000))), container_id], timeout=30)
        metadata.update(
            stdout_tail=logs.stdout[-12000:],
            stderr_tail=logs.stderr[-4000:],
            updated_at=_now(),
        )
        self._write_job(job_dir, metadata)
        return {
            "ok": logs.returncode == 0,
            "job_id": resolved_job_id,
            "status": metadata.get("status", "unknown"),
            "stdout_tail": metadata["stdout_tail"],
            "stderr_tail": metadata["stderr_tail"],
        }

    def cancel_job(self, job_id: str = "") -> dict:
        loaded = self._load_job_or_error(job_id)
        if isinstance(loaded, dict):
            return loaded
        job_dir, metadata = loaded
        resolved_job_id = str(metadata["job_id"])
        if metadata.get("status") in TERMINAL_JOB_STATUSES:
            return {
                "ok": False,
                "job_id": resolved_job_id,
                "status": metadata.get("status"),
                "error": "Job is already in a terminal state.",
            }
        container_id = str(metadata.get("container_id") or "")
        removed = _run(["docker", "rm", "-f", container_id], timeout=30)
        inspected = _run(["docker", "inspect", container_id], timeout=30) if removed.returncode == 0 else removed
        verified_exit = removed.returncode == 0 and inspected.returncode != 0
        metadata.update(
            status="cancelled" if verified_exit else metadata.get("status", "running"),
            updated_at=_now(),
            cancelled_at=_now() if verified_exit else "",
            cancel_verified=verified_exit,
        )
        self._write_job(job_dir, metadata)
        return {
            "ok": verified_exit,
            "job_id": resolved_job_id,
            "status": metadata["status"],
            "verified_exit": verified_exit,
            "error": "" if verified_exit else (removed.stderr.strip() or "Container still exists after cancellation request."),
        }

    def mark_callback_delivered(self, job_id: str, callback_id: str) -> None:
        with JOB_METADATA_LOCK:
            loaded = self._load_job_or_error(job_id)
            if isinstance(loaded, dict):
                return
            job_dir, metadata = loaded
            metadata.update(
                callback_state="delivered",
                callback_id=callback_id,
                callback_delivered_at=_now(),
                updated_at=_now(),
            )
            self._write_job(job_dir, metadata)
            run_dir = job_dir.parent.parent
            metadata["execution_manifest_path"] = str(_write_execution_manifest(run_dir, metadata))
            self._write_job(job_dir, metadata)

    def _poll_once(self, job_dir: Path, metadata: dict[str, Any]) -> dict:
        with JOB_METADATA_LOCK:
            return self._poll_once_locked(job_dir, metadata)

    def _poll_once_locked(self, job_dir: Path, metadata: dict[str, Any]) -> dict:
        try:
            persisted = json.loads((job_dir / "job.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            persisted = None
        if isinstance(persisted, dict):
            metadata = persisted
        if metadata.get("status") in TERMINAL_JOB_STATUSES:
            return self._public_result(metadata)
        container_id = str(metadata.get("container_id") or "")
        inspected = _run(["docker", "inspect", "--format", "{{json .State}}", container_id], timeout=30)
        if inspected.returncode != 0:
            metadata.update(status="cancelled", updated_at=_now(), error=inspected.stderr.strip() or "Container no longer exists.")
            self._write_job(job_dir, metadata)
            return self._public_result(metadata)
        try:
            state = json.loads(inspected.stdout)
        except json.JSONDecodeError:
            return {"ok": False, "job_id": metadata["job_id"], "status": "unknown", "error": "Invalid docker inspect output."}
        if bool(state.get("Running")) or state.get("Status") in {"created", "running", "restarting"}:
            logs = _run(["docker", "logs", "--tail", "100", container_id], timeout=30)
            metadata.update(
                status="running",
                stdout_tail=logs.stdout[-12000:],
                stderr_tail=logs.stderr[-4000:],
                updated_at=_now(),
            )
            self._write_job(job_dir, metadata)
            return self._public_result(metadata)
        exit_code = int(state.get("ExitCode") or 0)
        logs = _run(["docker", "logs", "--tail", "400", container_id], timeout=30)
        removed = _run(["docker", "rm", "-f", container_id], timeout=30)
        status = "timed_out" if exit_code == 124 else ("completed" if exit_code == 0 else "failed")
        job_run_dir = job_dir.parent.parent
        metadata.update(
            status=status,
            exit_code=exit_code,
            stdout_tail=logs.stdout[-12000:],
            stderr_tail=logs.stderr[-4000:],
            error=str(state.get("Error") or ""),
            updated_at=_now(),
            finished_at=_now(),
            artifacts=_output_artifacts(job_run_dir, metadata),
            container_removed=removed.returncode == 0,
        )
        self._write_job(job_dir, metadata)
        metadata["execution_manifest_path"] = str(_write_execution_manifest(job_run_dir, metadata))
        self._write_job(job_dir, metadata)
        self._progress("Docker job finished", f"job_id={metadata['job_id']} status={status} exit_code={exit_code}")
        return self._public_result(metadata)

    def _inner_command(
        self,
        *,
        runtime: str,
        docker_script: str,
        requirements: str,
        job_dir: Path,
        timeout_s: int,
    ) -> str:
        quoted_script = shlex.quote(docker_script)
        timeout = f"timeout --signal=TERM --kill-after=30s {timeout_s}s"
        if runtime == "r":
            return f"set -e\nRscript -e {shlex.quote(f'parse(file={docker_script!r})')}\n{timeout} Rscript {quoted_script}\n"
        lines = ["set -e"]
        if requirements.strip():
            requirements_path = job_dir / "requirements.txt"
            requirements_path.write_text(requirements.strip() + "\n", encoding="utf-8")
            docker_requirements = f"/work/{run_relative_path(self.run_dir, requirements_path)}"
            lines.append(f"python -m pip install -r {shlex.quote(docker_requirements)}")
        lines.append(f"{timeout} python {quoted_script}")
        return "\n".join(lines) + "\n"

    def active_job_id(self) -> str:
        for run_dir in self.run_dirs:
            try:
                state = json.loads((run_dir / "state" / "job_state.json").read_text(encoding="utf-8"))
            except (FileNotFoundError, OSError, json.JSONDecodeError):
                continue
            job_id = str(state.get("active_job_id") or "")
            if JOB_ID_PATTERN.fullmatch(job_id) and (run_dir / "jobs" / job_id / "job.json").is_file():
                return job_id
        candidates = sorted(
            (path for _, path in self._job_paths()),
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
        if not candidates:
            return ""
        fallback = candidates[0].parent.name
        if JOB_ID_PATTERN.fullmatch(fallback):
            return fallback
        return ""

    def _set_active_job(self, job_id: str) -> None:
        self.job_state_path.parent.mkdir(parents=True, exist_ok=True)
        temporary = self.job_state_path.with_suffix(".tmp")
        temporary.write_text(
            json.dumps({"active_job_id": job_id, "updated_at": _now()}, ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        temporary.replace(self.job_state_path)

    def _load_job_or_error(self, job_id: str) -> tuple[Path, dict[str, Any]] | dict[str, Any]:
        requested = str(job_id or "").strip()
        active = self._running_job_id() or self.active_job_id()
        resolved = requested or active
        available = [item["job_id"] for item in self.list_jobs()["jobs"]]
        if not resolved:
            return {
                "ok": False,
                "status": "no_active_job",
                "error": "No job has been started in this session. Call start_job first.",
                "active_job_id": "",
                "available_job_ids": available,
            }
        if not JOB_ID_PATTERN.fullmatch(resolved):
            return self._job_id_error(requested, active, available, "Job id was not issued by the harness.")
        job_dir = self._find_job_dir(resolved)
        if job_dir is None:
            return self._job_id_error(requested, active, available, "Job id is not registered in this session.")
        path = job_dir / "job.json"
        try:
            metadata = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return self._job_id_error(requested, active, available, "Persisted job metadata is unreadable.")
        if metadata.get("job_id") != resolved:
            return self._job_id_error(requested, active, available, "Persisted job id does not match its registry entry.")
        metadata.setdefault("run_id", job_dir.parent.parent.name)
        metadata.setdefault("run_dir", str(job_dir.parent.parent))
        return job_dir, metadata

    @staticmethod
    def _job_id_error(requested: str, active: str, available: list[str], reason: str) -> dict[str, Any]:
        return {
            "ok": False,
            "status": "invalid_job_id",
            "error": f"Invalid job_id={requested!r}. {reason} Use active_job_id or call list_jobs.",
            "active_job_id": active,
            "available_job_ids": available,
        }

    def _active_running_job(self) -> dict[str, Any] | None:
        for run_dir, path in self._job_paths():
            try:
                metadata = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if metadata.get("status") != "running":
                continue
            metadata.setdefault("run_id", run_dir.name)
            metadata.setdefault("run_dir", str(run_dir))
            refreshed = self._poll_once(path.parent, metadata)
            if refreshed.get("status") == "running":
                return refreshed
        return None

    def _running_job_id(self) -> str:
        for _, path in self._job_paths():
            try:
                metadata = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if metadata.get("status") == "running":
                return str(metadata.get("job_id") or "")
        return ""

    def _job_paths(self) -> list[tuple[Path, Path]]:
        paths = [
            (run_dir, path)
            for run_dir in self.run_dirs
            for path in (run_dir / "jobs").glob("*/job.json")
            if path.is_file()
        ]
        return sorted(paths, key=lambda item: item[1].stat().st_mtime, reverse=True)

    def _find_job_dir(self, job_id: str) -> Path | None:
        for run_dir in self.run_dirs:
            candidate = run_dir / "jobs" / job_id
            if (candidate / "job.json").is_file():
                return candidate
        return None

    @staticmethod
    def _write_job(job_dir: Path, metadata: dict[str, Any]) -> None:
        with JOB_METADATA_LOCK:
            path = job_dir / "job.json"
            temporary = path.with_name(f"{path.name}.{uuid.uuid4().hex[:8]}.tmp")
            temporary.write_text(json.dumps(metadata, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
            temporary.replace(path)

    @staticmethod
    def _public_result(metadata: dict[str, Any]) -> dict:
        status = str(metadata.get("status") or "unknown")
        return {
            "ok": status in {"running", "completed"},
            **metadata,
            "error_reason": metadata.get("error") or (f"Job ended with status={status}" if status not in {"running", "completed"} else ""),
        }

    def _progress(self, title: str, detail: str) -> None:
        if self.logger is not None:
            self.logger.progress(title, detail)


def cancel_run_jobs(run_dir: Path) -> None:
    jobs_dir = run_dir / "jobs"
    if not jobs_dir.is_dir():
        return
    for cidfile in jobs_dir.glob("*/container.cid"):
        container_id = cidfile.read_text(encoding="utf-8").strip()
        if container_id:
            _run(["docker", "rm", "-f", container_id], timeout=30)
        metadata_path = cidfile.parent / "job.json"
        if metadata_path.is_file():
            try:
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                continue
            if metadata.get("status") == "running":
                metadata.update(status="cancelled", cancelled_at=_now(), updated_at=_now(), error="Cancelled by user pause.")
                DockerJobManager._write_job(cidfile.parent, metadata)


def _run(command: list[str], *, timeout: int) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
        check=False,
    )


def _output_artifacts(run_dir: Path, job: dict[str, Any]) -> list[dict[str, Any]]:
    output_dir = run_dir / "outputs"
    if not output_dir.is_dir():
        return []
    result = []
    for path in sorted(output_dir.rglob("*")):
        if not path.is_file():
            continue
        digest = _sha256(path)
        artifact_id = "artifact_" + hashlib.sha256(
            f"{job.get('job_id')}:{path.resolve()}:{digest}".encode("utf-8")
        ).hexdigest()[:16]
        result.append(
            {
                "artifact_id": artifact_id,
                "path": str(path.resolve()),
                "size_bytes": path.stat().st_size,
                "mime_type": mimetypes.guess_type(path.name)[0] or "application/octet-stream",
                "sha256": digest,
                "created_by_job_id": str(job.get("job_id") or ""),
                "created_by_tool_call_id": str(job.get("created_by_tool_call_id") or ""),
                "created_by_event_id": str(job.get("created_by_event_id") or ""),
                "script_path": str(job.get("script_path") or ""),
                "script_sha256": str(job.get("script_sha256") or ""),
                "runtime": str(job.get("runtime") or ""),
                "env_profile": str(job.get("env_profile") or ""),
                "image": str(job.get("image") or ""),
                "input_artifact_ids": [],
                "created_at": _now(),
            }
        )
    return result


def attach_job_provenance(
    run_dirs: list[Path],
    *,
    job_id: str,
    tool_call_id: str,
    event_id: str,
) -> None:
    with JOB_METADATA_LOCK:
        for run_dir in run_dirs:
            metadata_path = run_dir / "jobs" / job_id / "job.json"
            if not metadata_path.is_file():
                continue
            try:
                metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                return
            metadata.update(created_by_tool_call_id=tool_call_id, created_by_event_id=event_id, updated_at=_now())
            DockerJobManager._write_job(metadata_path.parent, metadata)
            metadata["execution_manifest_path"] = str(_write_execution_manifest(run_dir, metadata))
            DockerJobManager._write_job(metadata_path.parent, metadata)
            return


def _write_execution_manifest(run_dir: Path, job: dict[str, Any]) -> Path:
    with JOB_METADATA_LOCK:
        path = run_dir / "state" / "execution_manifest.json"
        try:
            manifest = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            manifest = {"version": 1, "jobs": {}, "artifacts": {}}
        if not isinstance(manifest, dict):
            manifest = {"version": 1, "jobs": {}, "artifacts": {}}
        jobs = manifest.setdefault("jobs", {})
        artifacts = manifest.setdefault("artifacts", {})
        jobs[str(job.get("job_id") or "unknown")] = {
            key: job.get(key)
            for key in (
                "job_id",
                "run_id",
                "status",
                "runtime",
                "env_profile",
                "image",
                "script_path",
                "script_sha256",
                "created_by_tool_call_id",
                "created_by_event_id",
                "created_at",
                "finished_at",
                "exit_code",
            )
            if job.get(key) is not None
        }
        for artifact in job.get("artifacts") or []:
            if isinstance(artifact, dict) and artifact.get("artifact_id"):
                artifacts[str(artifact["artifact_id"])] = artifact
        manifest["updated_at"] = _now()
        path.parent.mkdir(parents=True, exist_ok=True)
        temp = path.with_name(f"{path.name}.{uuid.uuid4().hex[:8]}.tmp")
        temp.write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
        temp.replace(path)
        return path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _compact_job(metadata: dict[str, Any], *, run_dir: Path) -> dict[str, Any]:
    compact = {
        key: metadata.get(key)
        for key in (
            "job_id",
            "status",
            "runtime",
            "env_profile",
            "script_path",
            "workspace_script_path",
            "created_at",
            "updated_at",
            "finished_at",
            "exit_code",
            "attempts_used",
            "callback_state",
        )
        if metadata.get(key) is not None
    }
    compact["run_id"] = str(metadata.get("run_id") or run_dir.name)
    compact["run_dir"] = str(metadata.get("run_dir") or run_dir)
    return compact


def _session_run_dirs(runs_root: Path, current_run: Path, prior_run_dirs: list[Path]) -> list[Path]:
    root = runs_root.resolve()
    result: list[Path] = []
    for candidate in [current_run, *prior_run_dirs]:
        resolved = Path(candidate).resolve()
        if resolved != root and root not in resolved.parents:
            continue
        if resolved not in result:
            result.append(resolved)
    return result
