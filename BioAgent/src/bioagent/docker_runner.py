from __future__ import annotations

import shutil
import subprocess
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

    def _run(self, command: list[str], timeout_s: int) -> subprocess.CompletedProcess[str]:
        self.logger.log("DOCKER COMMAND " + " ".join(command))
        return subprocess.run(
            command,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=timeout_s,
            check=False,
        )

    def execute_python(
        self,
        code: str,
        *,
        env_profile: str = "py-general-v1",
        requirements: str = "",
        timeout_s: int = 900,
    ) -> ExecutionResult:
        script = self.scripts_dir / "code.py"
        script.write_text(code, encoding="utf-8")
        req_file = self.scripts_dir / "requirements.txt"
        req_file.write_text(requirements.strip() + "\n" if requirements.strip() else "", encoding="utf-8")

        if not self._docker_available():
            return ExecutionResult(False, "python", env_profile, "", str(self.run_dir), str(script), 127, "", "docker CLI not found", [])
        try:
            image = self._image(env_profile, "python")
        except ValueError as exc:
            return ExecutionResult(False, "python", env_profile, "", str(self.run_dir), str(script), 2, "", str(exc), [])
        if not self._image_exists(image):
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

        command = [
            "docker",
            "run",
            "--rm",
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
            result = self._run(command, timeout_s)
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
            return ExecutionResult(False, "python", env_profile, image, str(self.run_dir), str(script), -1, exc.stdout or "", (exc.stderr or "") + "\nTimed out", command)

    def execute_r(
        self,
        code: str,
        *,
        env_profile: str = "r-bioc-v1",
        timeout_s: int = 900,
    ) -> ExecutionResult:
        script = self.scripts_dir / "code.R"
        script.write_text(code, encoding="utf-8")

        if not self._docker_available():
            return ExecutionResult(False, "r", env_profile, "", str(self.run_dir), str(script), 127, "", "docker CLI not found", [])
        try:
            image = self._image(env_profile, "r")
        except ValueError as exc:
            return ExecutionResult(False, "r", env_profile, "", str(self.run_dir), str(script), 2, "", str(exc), [])
        if not self._image_exists(image):
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

        command = [
            "docker",
            "run",
            "--rm",
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
            result = self._run(command, timeout_s)
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
            return ExecutionResult(False, "r", env_profile, image, str(self.run_dir), str(script), -1, exc.stdout or "", (exc.stderr or "") + "\nTimed out", command)

