from __future__ import annotations

import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

from .config import AgentConfig
from .logging_utils import RunLogger
from .skills import image_for_profile


@dataclass
class ExecutionResult:
    ok: bool
    runtime: str
    env_profile: str
    image: str
    run_dir: str
    script_path: str
    exit_code: int
    stdout: str
    stderr: str
    command: list[str]

    def to_dict(self) -> dict:
        return {
            "ok": self.ok,
            "runtime": self.runtime,
            "env_profile": self.env_profile,
            "image": self.image,
            "run_dir": self.run_dir,
            "script_path": self.script_path,
            "exit_code": self.exit_code,
            "stdout": self.stdout[-8000:],
            "stderr": self.stderr[-8000:],
            "error_reason": "" if self.ok else _summarize_process_error(self.stderr, self.stdout),
            "command": self.command,
        }


class DockerRunner:
    def __init__(self, config: AgentConfig, logger: RunLogger, run_dir: Path) -> None:
        self.config = config
        self.logger = logger
        self.run_dir = run_dir
        self.scripts_dir = run_dir / "scripts"
        self.outputs_dir = run_dir / "outputs"
        self.scripts_dir.mkdir(parents=True, exist_ok=True)
        self.outputs_dir.mkdir(parents=True, exist_ok=True)

    def _image(self, env_profile: str, runtime: str) -> str:
        entry = image_for_profile(self.config.docker_root, env_profile)
        if not entry:
            raise ValueError(f"Unknown Docker env_profile: {env_profile}. Check mas_2/docker/image_catalog.json.")
        image = str(entry.get("image_tag") or "").strip()
        entry_runtime = str(entry.get("runtime") or "").strip()
        if entry_runtime and runtime != "mixed" and entry_runtime not in {runtime, "mixed"}:
            raise ValueError(
                f"Profile {env_profile} is runtime={entry_runtime}, not requested runtime={runtime}."
            )
        if not image:
            raise ValueError(f"Docker env_profile {env_profile} has no image_tag.")
        return image

    def _docker_available(self) -> bool:
        return shutil.which("docker") is not None

    def _image_exists(self, image: str) -> bool:
        result = subprocess.run(
            ["docker", "image", "inspect", image],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=30,
            check=False,
        )
        return result.returncode == 0

    def _run(
        self,
        command: list[str],
        timeout_s: int,
        *,
        cidfile: Path | None = None,
    ) -> subprocess.CompletedProcess[str]:
        self.logger.progress("Docker 命令启动", f"timeout={timeout_s}s command={' '.join(command)}")
        started = time.perf_counter()
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                encoding="utf-8",
                errors="replace",
                timeout=timeout_s,
                check=False,
            )
        except subprocess.TimeoutExpired:
            self._remove_timed_out_container(cidfile)
            raise
        finally:
            if cidfile is not None:
                cidfile.unlink(missing_ok=True)
        elapsed = time.perf_counter() - started
        self.logger.progress(
            "Docker 命令结束",
            f"exit_code={result.returncode} elapsed={elapsed:.2f}s stdout_chars={len(result.stdout or '')} stderr_chars={len(result.stderr or '')}",
        )
        if result.returncode != 0:
            self.logger.error_reason(_summarize_process_error(result.stderr or "", result.stdout or ""))
        return result

    def _remove_timed_out_container(self, cidfile: Path | None) -> None:
        if cidfile is None or not cidfile.is_file():
            return
        container_id = cidfile.read_text(encoding="utf-8").strip()
        if not container_id:
            return
        self.logger.progress("清理超时 Docker 容器", f"container_id={container_id}")
        subprocess.run(
            ["docker", "rm", "-f", container_id],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=30,
            check=False,
        )

    def execute_python(
        self,
        code: str,
        *,
        env_profile: str = "py-general-v1",
        requirements: str = "",
        timeout_s: int = 1800,
    ) -> ExecutionResult:
        script = self.scripts_dir / "code.py"
        script.write_text(code, encoding="utf-8")
        req_file = self.scripts_dir / "requirements.txt"
        req_file.write_text(requirements.strip() + "\n" if requirements.strip() else "", encoding="utf-8")
        self.logger.progress(
            "准备 Python 执行脚本",
            f"script={script} env_profile={env_profile} requirements={'yes' if requirements.strip() else 'no'}",
        )

        if not self._docker_available():
            self.logger.error_reason("docker CLI not found，无法启动沙箱执行。")
            return ExecutionResult(False, "python", env_profile, "", str(self.run_dir), str(script), 127, "", "docker CLI not found", [])
        try:
            image = self._image(env_profile, "python")
        except ValueError as exc:
            self.logger.error_reason(str(exc))
            return ExecutionResult(False, "python", env_profile, "", str(self.run_dir), str(script), 2, "", str(exc), [])
        self.logger.progress("选择 Docker 镜像", f"runtime=python env_profile={env_profile} image={image}")
        if not self._image_exists(image):
            self.logger.error_reason(f"Docker image not found locally: {image}. 请先构建 mas_2/docker 中对应镜像。")
            return ExecutionResult(
                False,
                "python",
                env_profile,
                image,
                str(self.run_dir),
                str(script),
                125,
                "",
                f"Docker image not found locally: {image}. Build it from mas_2/docker before running this profile.",
                [],
            )

        inner = "set -e\n"
        if requirements.strip():
            inner += "python -m pip install -r /work/scripts/requirements.txt\n"
        inner += "python /work/scripts/code.py\n"

        cidfile = self.run_dir / ".docker-python.cid"
        cidfile.unlink(missing_ok=True)
        command = [
            "docker",
            "run",
            "--rm",
            "--cidfile",
            str(cidfile),
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
        try:
            result = self._run(command, timeout_s, cidfile=cidfile)
            return ExecutionResult(
                result.returncode == 0,
                "python",
                env_profile,
                image,
                str(self.run_dir),
                str(script),
                result.returncode,
                result.stdout or "",
                result.stderr or "",
                command,
            )
        except subprocess.TimeoutExpired as exc:
            self.logger.error_reason(f"Docker Python 执行超过 timeout={timeout_s}s。")
            stdout = _decode_process_output(exc.stdout)
            stderr = _decode_process_output(exc.stderr)
            return ExecutionResult(False, "python", env_profile, image, str(self.run_dir), str(script), -1, stdout, stderr + "\nTimed out", command)

    def execute_r(
        self,
        code: str,
        *,
        env_profile: str = "r-bioc-v1",
        timeout_s: int = 1800,
    ) -> ExecutionResult:
        script = self.scripts_dir / "code.R"
        script.write_text(code, encoding="utf-8")
        self.logger.progress("准备 R 执行脚本", f"script={script} env_profile={env_profile}")

        if not self._docker_available():
            self.logger.error_reason("docker CLI not found，无法启动沙箱执行。")
            return ExecutionResult(False, "r", env_profile, "", str(self.run_dir), str(script), 127, "", "docker CLI not found", [])
        try:
            image = self._image(env_profile, "r")
        except ValueError as exc:
            self.logger.error_reason(str(exc))
            return ExecutionResult(False, "r", env_profile, "", str(self.run_dir), str(script), 2, "", str(exc), [])
        self.logger.progress("选择 Docker 镜像", f"runtime=r env_profile={env_profile} image={image}")
        if not self._image_exists(image):
            self.logger.error_reason(f"Docker image not found locally: {image}. 请先构建 mas_2/docker 中对应镜像。")
            return ExecutionResult(
                False,
                "r",
                env_profile,
                image,
                str(self.run_dir),
                str(script),
                125,
                "",
                f"Docker image not found locally: {image}. Build it from mas_2/docker before running this profile.",
                [],
            )

        cidfile = self.run_dir / ".docker-r.cid"
        cidfile.unlink(missing_ok=True)
        command = [
            "docker",
            "run",
            "--rm",
            "--cidfile",
            str(cidfile),
            "-v",
            f"{self.config.repo_root}:/repo:ro",
            "-v",
            f"{self.run_dir}:/work",
            "-w",
            "/work",
            image,
            "Rscript",
            "/work/scripts/code.R",
        ]
        try:
            result = self._run(command, timeout_s, cidfile=cidfile)
            return ExecutionResult(
                result.returncode == 0,
                "r",
                env_profile,
                image,
                str(self.run_dir),
                str(script),
                result.returncode,
                result.stdout or "",
                result.stderr or "",
                command,
            )
        except subprocess.TimeoutExpired as exc:
            self.logger.error_reason(f"Docker R 执行超过 timeout={timeout_s}s。")
            stdout = _decode_process_output(exc.stdout)
            stderr = _decode_process_output(exc.stderr)
            return ExecutionResult(False, "r", env_profile, image, str(self.run_dir), str(script), -1, stdout, stderr + "\nTimed out", command)


def _decode_process_output(value: str | bytes | None) -> str:
    if value is None:
        return ""
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return value


def _summarize_process_error(stderr: str, stdout: str) -> str:
    text = (stderr or "").strip() or (stdout or "").strip()
    if not text:
        return "进程返回非零退出码，但 stdout/stderr 为空。"
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
