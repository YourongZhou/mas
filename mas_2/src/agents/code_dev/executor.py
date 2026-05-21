"""
代码执行器模块
使用 Docker 容器安全执行生成的代码
"""
import os
import shlex
import tempfile
import docker
import shutil
import logging
import json
import re
import time
from datetime import datetime, timezone
from pathlib import Path

from src.utils.execution_environment import get_image_catalog_entry, image_catalog_root

logging.basicConfig(level=logging.INFO)


def _env_truthy(name: str) -> bool:
    return os.environ.get(name, "").strip().lower() in ("1", "true", "yes", "on")


_ENV_CACHE_ROOT = "/opt/mas_env_cache"
_ENV_CACHE_DEPS_DIR = f"{_ENV_CACHE_ROOT}/pydeps"
_ENV_CACHE_MANIFEST = f"{_ENV_CACHE_ROOT}/manifest.json"
_ENV_CACHE_LOCK = f"{_ENV_CACHE_ROOT}/install.lock"


class CodeExecutor:
    """使用 Docker 容器执行代码的执行器"""
    
    def __init__(self, docker_path: str = None, data_dir: str = None, data_dirs: list = None,
                 output_dir: str = None, input_files: list = None,
                 workflow_host_path: str = None, container_id: str = None,
                 docker_image: str | None = None, runtime: str = "python",
                 env_profile: str | None = None, env_signature: str | None = None,
                 env_mode: str = "isolated", asset_cache_host_path: str | None = None,
                 env_cache_key: str | None = None, env_cache_volume: str | None = None,
                 required_input_container_paths: list[str] | None = None):
        """
        初始化代码执行器

        workflow_host_path: 若提供，将宿主目录只读挂载到容器内 /app/workflow（供 SKILL 内 scripts 导入）。
        """
        self.logger = logging.getLogger(__name__)
        self.docker_available = self._check_docker_availability()
        self.code_path = f"{docker_path}/code.py" if docker_path else None
        self.requirements_path = f"{docker_path}/requirements.txt" if docker_path else None
        self.output_dir = output_dir if output_dir else "/tmp/output"  # 默认输出目录
        self.docker_image = docker_image or "python:3.11"
        self.container_id = container_id
        self.runtime = (runtime or "python").strip() or "python"
        self.env_profile = (env_profile or "").strip() or None
        self.env_signature = (env_signature or "").strip() or None
        self.env_mode = (env_mode or "isolated").strip() or "isolated"
        self.asset_cache_host_path = asset_cache_host_path
        self.env_cache_key = (env_cache_key or "").strip() or None
        self.env_cache_volume = (env_cache_volume or "").strip() or None
        self.required_input_container_paths = [
            str(path).strip()
            for path in (required_input_container_paths or [])
            if str(path).strip()
        ]
        
        # --- 关键修改 1：使用列表存储挂载信息，避免字典 Key 冲突 ---
        self.volume_mounts = []
        
        # 处理数据目录：支持单个或多个目录
        # 确保 self.data_dirs 是有序列表，且 data_dir 始终在第一个
        self.data_dirs = []
        if data_dirs:
            self.data_dirs = [d for d in data_dirs if d and os.path.exists(d)]
        elif data_dir:
            self.data_dirs = [data_dir] if os.path.exists(data_dir) else []
        
        if input_files:
            self._determine_data_dirs_from_input_files(input_files)
        
        self.data_mount_path = self.data_dirs[0] if self.data_dirs else None
        self.data_file_name = None
        
        # (不再挂载 code.py 和 requirements.txt，改为运行时推送)

        # 挂载数据目录（支持多个）
        for idx, data_dir_path in enumerate(self.data_dirs):
            if not os.path.exists(data_dir_path):
                continue
                
            data_mount_path = os.path.abspath(data_dir_path)
            if os.path.isfile(data_mount_path):
                data_mount_path = os.path.dirname(data_mount_path)
                if idx == 0:
                    self.data_file_name = os.path.basename(data_dir_path)
            
            if idx == 0:
                bind_path = '/app/data'
            else:
                bind_path = f'/app/data{idx}'
            
            # --- 关键修改 2：追加到列表，允许同一个 host path 挂载到多个位置 ---
            # 使用 rw 模式
            mount_str = f"{data_mount_path}:{bind_path}:rw"
            self.volume_mounts.append(mount_str)
            if not self.container_id:
                self.logger.info(f"数据目录挂载: {data_mount_path} -> {bind_path}")

        # 挂载输出目录
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
            
        # --- 关键修改 3：追加输出目录挂载 ---
        mount_str = f"{os.path.abspath(self.output_dir)}:/app/output:rw"
        self.volume_mounts.append(mount_str)

        self.workflow_host_path = workflow_host_path
        if workflow_host_path and os.path.isdir(workflow_host_path):
            wf = os.path.abspath(workflow_host_path)
            self.volume_mounts.append(f"{wf}:/app/workflows:ro")
            if not self.container_id:
                self.logger.info(f"Workflow 目录挂载: {wf} -> /app/workflows")

        if asset_cache_host_path and os.path.isdir(asset_cache_host_path):
            cache_path = os.path.abspath(asset_cache_host_path)
            self.volume_mounts.append(f"{cache_path}:/app/assets:rw")
            if not self.container_id:
                self.logger.info(f"Asset cache 挂载: {cache_path} -> /app/assets")

        if self.env_cache_volume:
            self.volume_mounts.append(f"{self.env_cache_volume}:{_ENV_CACHE_ROOT}:rw")
            if not self.container_id:
                self.logger.info(f"Env cache 挂载: {self.env_cache_volume} -> {_ENV_CACHE_ROOT}")
    
    def _determine_data_dirs_from_input_files(self, input_files: list):
        """
        根据输入文件列表智能确定需要挂载的数据目录
        """
        if not input_files:
            return
        
        existing_set = {os.path.abspath(d) for d in self.data_dirs}
        dirs_to_add = []
        
        for input_file in input_files:
            if not input_file:
                continue
            
            found_dir = None
            if os.path.isabs(input_file) and os.path.exists(input_file):
                file_dir = os.path.dirname(input_file) if os.path.isfile(input_file) else input_file
                found_dir = os.path.abspath(file_dir)
            elif not os.path.isabs(input_file):
                for existing_dir in self.data_dirs:
                    candidate = os.path.join(existing_dir, input_file)
                    if os.path.exists(candidate):
                        found_dir = None
                        break
                    if os.path.isfile(existing_dir):
                        existing_dir = os.path.dirname(existing_dir)
                    candidate = os.path.join(existing_dir, os.path.basename(input_file))
                    if os.path.exists(candidate):
                        found_dir = None
                        break
            
            if found_dir and found_dir not in existing_set:
                dirs_to_add.append(found_dir)
                existing_set.add(found_dir)
        
        if dirs_to_add:
            self.data_dirs.extend(dirs_to_add)
            self.logger.info(f"根据输入文件追加挂载目录: {dirs_to_add}")

    def _check_docker_availability(self) -> bool:
        try:
            self.client = docker.from_env()
            return True
        except ImportError:
            self.logger.error("FAILED.Docker模块不可用，请安装: pip install docker")
            return False
        except Exception as e:
            self.logger.error(f"FAILED.Docker客户端初始化失败: {e}")
            return False

    def _shared_container_labels(self) -> dict[str, str]:
        labels = {
            "mas.managed": "1",
            "mas.runtime": self.runtime,
            "mas.env_mode": self.env_mode,
            "mas.last_used_at": datetime.now(timezone.utc).isoformat(),
        }
        if self.env_profile:
            labels["mas.env_profile"] = self.env_profile
        if self.env_signature:
            labels["mas.env_signature"] = self.env_signature
        return labels

    def _ensure_image_available(self) -> None:
        try:
            self.client.images.get(self.docker_image)
            return
        except docker.errors.ImageNotFound:
            pass

        entry = get_image_catalog_entry(self.env_profile)
        if not _env_truthy("MAS_DOCKER_DEV_BOOTSTRAP"):
            raise RuntimeError(
                f"所需 Docker 镜像不存在: {self.docker_image}。"
                "如为本地开发，请设置 MAS_DOCKER_DEV_BOOTSTRAP=1 后重试。"
            )
        if entry is None:
            raise RuntimeError(f"镜像缺失且无法 bootstrap：未找到 env_profile={self.env_profile} 的 image catalog")

        docker_root = image_catalog_root()
        dockerfile_path = (docker_root / entry.dockerfile).resolve()
        context_path = (docker_root / entry.build_context).resolve()
        if not dockerfile_path.is_file():
            raise RuntimeError(f"Dockerfile 不存在，无法 bootstrap: {dockerfile_path}")
        dockerfile_arg = os.path.relpath(dockerfile_path, start=context_path)
        self.logger.info(
            "ENV IMAGE BOOTSTRAP profile=%s image=%s dockerfile=%s",
            self.env_profile or "unknown",
            self.docker_image,
            dockerfile_path,
        )
        self.client.images.build(
            path=str(context_path),
            dockerfile=dockerfile_arg,
            tag=self.docker_image,
            rm=True,
        )

    def _ensure_volume_exists(self, volume_name: str) -> None:
        try:
            self.client.volumes.get(volume_name)
        except docker.errors.NotFound:
            self.client.volumes.create(name=volume_name)

    def _create_container_env_vars(self, environment_vars: dict | None = None) -> dict[str, str]:
        env_vars = {
            'PYTHONUNBUFFERED': '1',
            'DEBIAN_FRONTEND': 'noninteractive',
            'MPLCONFIGDIR': '/tmp/mas_runtime/.mas_mpl',
            'NUMBA_CACHE_DIR': '/tmp/mas_runtime/.mas_numba',
            'HOME': '/tmp',
            'PYTHONPYCACHEPREFIX': '/tmp/mas_runtime/.mas_pycache',
            'PIP_DISABLE_PIP_VERSION_CHECK': '1',
        }
        if environment_vars:
            env_vars.update(environment_vars)
        return env_vars

    def _runtime_env_prelude(self, deps_dir: str) -> str:
        return f"""set -e
export PYTHONPATH="/app/workflows${{PYTHONPATH:+:$PYTHONPATH}}"
RUNTIME_ROOT=/tmp/mas_runtime
DEPS="{deps_dir}"
export TMPDIR="$RUNTIME_ROOT/.mas_tmp"
mkdir -p "$TMPDIR" "$DEPS" "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$PYTHONPYCACHEPREFIX"
"""

    def ensure_env_cache_ready(
        self,
        requirements_str: str | None = None,
        environment_vars: dict | None = None,
        mem_limit: str | None = '4g',
        timeout: int | None = 900,
    ) -> dict:
        if not self.docker_available:
            return {"success": False, "error": "Docker不可用", "output": ""}
        if self.runtime != "python":
            return {
                "success": False,
                "error": f"当前执行器仅支持 Python runtime，收到 runtime={self.runtime}",
                "output": "",
            }
        if not self.env_cache_volume:
            return {
                "success": True,
                "output": "env cache disabled",
                "timing": {
                    "cache_prepare_elapsed_seconds": 0.0,
                    "cache_hit": True,
                    "cache_install_performed": False,
                },
            }

        req_lines = []
        for raw_line in str(requirements_str or "").splitlines():
            line = raw_line.strip()
            if line and not line.startswith("#"):
                req_lines.append(line)

        self._ensure_image_available()
        self._ensure_volume_exists(self.env_cache_volume)

        total_start = time.perf_counter()
        bootstrap = None
        try:
            env_vars = self._create_container_env_vars(environment_vars)
            run_kwargs = {
                'image': self.docker_image,
                'command': ['sleep', 'infinity'],
                'volumes': [f"{self.env_cache_volume}:{_ENV_CACHE_ROOT}:rw"],
                'environment': env_vars,
                'mem_limit': mem_limit,
                'network_mode': 'bridge',
                'detach': True,
                'auto_remove': False,
                'labels': {
                    'mas.managed': '1',
                    'mas.runtime': self.runtime,
                    'mas.kind': 'env-bootstrap',
                    'mas.env_profile': self.env_profile or '',
                    'mas.env_signature': self.env_signature or '',
                    'mas.env_cache_key': self.env_cache_key or '',
                },
            }

            bootstrap = self.client.containers.run(**run_kwargs)
            payload = {
                "env_signature": self.env_signature or "",
                "env_profile": self.env_profile or "",
                "env_cache_key": self.env_cache_key or "",
                "requirements": req_lines,
                "manifest_path": _ENV_CACHE_MANIFEST,
                "lock_path": _ENV_CACHE_LOCK,
                "deps_dir": _ENV_CACHE_DEPS_DIR,
            }
            bootstrap_script = self._runtime_env_prelude(_ENV_CACHE_DEPS_DIR) + (
                "python - <<'PY'\n"
                "import fcntl, json, os, subprocess, sys\n"
                f"payload = json.loads({json.dumps(json.dumps(payload, ensure_ascii=True))})\n"
                "cache_root = os.path.dirname(payload['manifest_path'])\n"
                "os.makedirs(cache_root, exist_ok=True)\n"
                "os.makedirs(payload['deps_dir'], exist_ok=True)\n"
                "os.chmod(cache_root, 0o777)\n"
                "os.chmod(payload['deps_dir'], 0o777)\n"
                "open(payload['lock_path'], 'a', encoding='utf-8').close()\n"
                "os.chmod(payload['lock_path'], 0o666)\n"
                "with open(payload['lock_path'], 'r+', encoding='utf-8') as fh:\n"
                "    fcntl.flock(fh.fileno(), fcntl.LOCK_EX)\n"
                "    manifest = {}\n"
                "    if os.path.exists(payload['manifest_path']):\n"
                "        try:\n"
                "            with open(payload['manifest_path'], 'r', encoding='utf-8') as mf:\n"
                "                manifest = json.load(mf)\n"
                "        except Exception:\n"
                "            manifest = {}\n"
                "    expected = {\n"
                "        'env_signature': payload['env_signature'],\n"
                "        'env_profile': payload['env_profile'],\n"
                "        'env_cache_key': payload['env_cache_key'],\n"
                "        'installed_requirements': payload['requirements'],\n"
                "    }\n"
                "    hit = (\n"
                "        manifest.get('env_signature') == expected['env_signature']\n"
                "        and manifest.get('env_cache_key') == expected['env_cache_key']\n"
                "        and manifest.get('installed_requirements') == expected['installed_requirements']\n"
                "    )\n"
                "    if not hit and payload['requirements']:\n"
                "        cmd = [\n"
                "            sys.executable, '-m', 'pip', 'install',\n"
                "            '-i', 'https://pypi.tuna.tsinghua.edu.cn/simple',\n"
                "            '--target', payload['deps_dir'],\n"
                "            *payload['requirements'],\n"
                "        ]\n"
                "        subprocess.check_call(cmd)\n"
                "    manifest = {\n"
                "        **expected,\n"
                "        'last_used_at': __import__('datetime').datetime.utcnow().isoformat() + 'Z',\n"
                "    }\n"
                "    with open(payload['manifest_path'], 'w', encoding='utf-8') as mf:\n"
                "        json.dump(manifest, mf, ensure_ascii=False, indent=2)\n"
                "    os.chmod(payload['manifest_path'], 0o666)\n"
                "    print('ENV CACHE HIT' if hit else 'ENV CACHE READY')\n"
                "PY"
            )
            command = "bash -lc " + shlex.quote(bootstrap_script)
            if timeout is not None:
                command = f"timeout {int(timeout)} " + command
            exec_result = bootstrap.exec_run(cmd=command, user='', demux=False)
            output = exec_result.output.decode('utf-8', errors='replace').strip() if exec_result.output else ""
            elapsed = round(time.perf_counter() - total_start, 3)
            return {
                "success": exec_result.exit_code == 0,
                "error": "" if exec_result.exit_code == 0 else self._extract_error_summary(output) or "环境缓存准备失败",
                "output": output,
                "exit_code": exec_result.exit_code,
                "timing": {
                    "cache_prepare_elapsed_seconds": elapsed,
                    "cache_hit": "ENV CACHE HIT" in output,
                    "cache_install_performed": "ENV CACHE READY" in output,
                },
            }
        except Exception as exc:
            return {
                "success": False,
                "error": str(exc),
                "output": "",
                "timing": {
                    "cache_prepare_elapsed_seconds": round(time.perf_counter() - total_start, 3),
                    "cache_hit": False,
                    "cache_install_performed": False,
                },
            }
        finally:
            if bootstrap is not None:
                try:
                    bootstrap.remove(force=True)
                except Exception:
                    pass

    def _run_required_inputs_preflight(self, container, user: str, timeout: int | None) -> dict | None:
        if not self.required_input_container_paths:
            return None
        checks = "\n".join(
            [
                f'if [ ! -r {shlex.quote(path)} ]; then echo "MISSING_INPUT::{path}"; exit 31; fi'
                for path in self.required_input_container_paths
            ]
        )
        script = self._runtime_env_prelude(
            _ENV_CACHE_DEPS_DIR if self.env_cache_volume else "/tmp/mas_runtime/.mas_pydeps"
        ) + (
            'echo "===PRE_EXEC_PREFLIGHT==="\n'
            'if [ -d /app/data ]; then ls -1 /app/data || true; fi\n'
            f"{checks}\n"
            'echo "PREFLIGHT_OK"\n'
        )
        command = "bash -lc " + shlex.quote(script)
        if timeout is not None:
            command = f"timeout {int(timeout)} " + command
        result = container.exec_run(cmd=command, user=user, demux=False)
        output = result.output.decode('utf-8', errors='replace').strip() if result.output else ""
        if result.exit_code == 0:
            return None
        missing_lines = [
            line.split("MISSING_INPUT::", 1)[1].strip()
            for line in output.splitlines()
            if "MISSING_INPUT::" in line
        ]
        err = "挂载预检查失败"
        if missing_lines:
            err = "挂载预检查失败，容器内缺少输入文件: " + ", ".join(missing_lines)
        elif output:
            err = self._extract_error_summary(output) or err
        return {
            "success": False,
            "error": err,
            "output": output,
            "exit_code": result.exit_code,
        }

    def _prepare_temp_directory(self, temp_dir: str) -> None:
        shutil.copy2(self.code_path, os.path.join(temp_dir, 'code.py'))
        if os.path.exists(self.requirements_path):
            shutil.copy2(self.requirements_path, os.path.join(temp_dir, 'requirements.txt'))

    @staticmethod
    def _extract_error_summary(logs: str, max_chars: int = 1200) -> str:
        """从容器日志提取末尾最关键错误，避免安装过程噪声。"""
        text = (logs or "").strip()
        if not text:
            return ""

        tb_matches = list(
            re.finditer(
                r"Traceback \(most recent call last\):[\s\S]*?(?=\n\n|\Z)",
                text,
                re.DOTALL,
            )
        )
        if tb_matches:
            return tb_matches[-1].group(0).strip()[-max_chars:]

        lines = text.splitlines()
        err_keys = (
            "error:",
            "exception",
            "failed",
            "module not found",
            "modulenotfounderror",
            "importerror",
            "typeerror",
            "valueerror",
            "attributeerror",
            "nameerror",
            "runtimeerror",
            "syntaxerror",
            "indentationerror",
            "assertionerror",
            "unable to",
        )
        for i in range(len(lines) - 1, -1, -1):
            line_low = lines[i].lower()
            if any(k in line_low for k in err_keys):
                start = max(0, i - 6)
                return "\n".join(lines[start:i + 1]).strip()[-max_chars:]

        return "\n".join(lines[-20:]).strip()[-max_chars:]

    def execute(self,
                code_str: str = None,
                requirements_str: str = None,
                environment_vars: dict | None = None,
                mem_limit: str | None = '4g',
                timeout: int | None = 600) -> dict:
        """
        使用复用或新建的容器执行代码
        """
        if not self.docker_available:
            return {
                'success': False,
                'error': 'Docker不可用',
                'output': '',
                'files': []
            }
        if self.runtime != "python":
            return {
                'success': False,
                'error': f'当前执行器仅支持 Python runtime，收到 runtime={self.runtime}',
                'output': '',
                'files': [],
                'container_id': self.container_id,
            }

        container = None
        total_start = time.perf_counter()
        timing: dict[str, float | int | bool] = {
            'container_setup_elapsed_seconds': 0.0,
            'push_to_container_elapsed_seconds': 0.0,
            'pip_install_elapsed_seconds': 0.0,
            'python_exec_elapsed_seconds': 0.0,
            'collect_outputs_elapsed_seconds': 0.0,
            'total_elapsed_seconds': 0.0,
            'pip_install_skipped': True,
            'install_exit_code': None,
            'python_exit_code': None,
        }

        def _remaining_timeout_seconds() -> int | None:
            if timeout is None:
                return None
            elapsed = time.perf_counter() - total_start
            remaining = timeout - elapsed
            if remaining <= 0:
                return 0
            return max(1, int(remaining))

        def _exec_phase(script: str, phase_name: str) -> dict:
            phase_start = time.perf_counter()
            phase_timeout = _remaining_timeout_seconds()
            if timeout is not None and phase_timeout == 0:
                return {
                    'exit_code': 124,
                    'output': f'[{phase_name}] timeout budget exhausted before phase started',
                    'elapsed_seconds': time.perf_counter() - phase_start,
                }

            command = "bash -lc " + shlex.quote(script)
            if phase_timeout is not None:
                command = f"timeout {phase_timeout} " + command

            exec_result = container.exec_run(
                cmd=command,
                user=user,
                demux=False
            )
            phase_output = exec_result.output.decode('utf-8', errors='replace').strip() if exec_result.output else ""
            return {
                'exit_code': exec_result.exit_code,
                'output': phase_output,
                'elapsed_seconds': time.perf_counter() - phase_start,
            }

        def _collect_output_files() -> list[dict]:
            collect_start = time.perf_counter()
            output_files: list[dict] = []
            if self.output_dir:
                output_path = Path(self.output_dir)
                for file_path in output_path.rglob('*'):
                    if file_path.is_file():
                        if any(part.startswith('.mas_') for part in file_path.parts):
                            continue
                        output_files.append({
                            'path': str(file_path),
                            'name': file_path.name,
                            'size': file_path.stat().st_size,
                            'size_mb': file_path.stat().st_size / (1024 * 1024)
                        })
            timing['collect_outputs_elapsed_seconds'] = round(time.perf_counter() - collect_start, 3)
            return output_files

        try:
            # 1. 尝试获取已有容器
            container_setup_start = time.perf_counter()
            if self.container_id:
                try:
                    container = self.client.containers.get(self.container_id)
                    labels = (container.attrs.get("Config", {}) or {}).get("Labels", {}) or {}
                    if self.env_signature:
                        current_sig = (labels.get("mas.env_signature") or "").strip()
                        current_profile = (labels.get("mas.env_profile") or "legacy").strip() or "legacy"
                        if current_sig != self.env_signature:
                            self.logger.info(
                                "ENV SWITCH from=%s to=%s",
                                current_profile,
                                self.env_profile or "legacy",
                            )
                            container = None
                        else:
                            self.logger.info(
                                "ENV REUSE profile=%s container=%s",
                                self.env_profile or current_profile,
                                self.container_id[:12],
                            )
                    if container is not None and container.status != 'running':
                        container.start()
                except docker.errors.NotFound:
                    container = None

            # 2. 如果无复用容器则新建
            if container is None:
                self._ensure_image_available()
                if self.env_cache_volume:
                    self._ensure_volume_exists(self.env_cache_volume)
                env_vars = self._create_container_env_vars(environment_vars)

                run_kwargs = {
                    'image': self.docker_image,
                    'command': ['sleep', 'infinity'],
                    'volumes': self.volume_mounts,
                    'environment': env_vars,
                    'mem_limit': mem_limit,
                    'network_mode': 'bridge',
                    'detach': True,
                    'auto_remove': False,
                    'labels': self._shared_container_labels(),
                }
                if hasattr(os, 'getuid') and hasattr(os, 'getgid'):
                    try:
                        run_kwargs['user'] = f"{os.getuid()}:{os.getgid()}"
                    except Exception as e:
                        self.logger.warning(f"无法获取当前用户 UID/GID: {e}")

                self.logger.info(
                    "ENV CREATE profile=%s image=%s mode=%s",
                    self.env_profile or "legacy",
                    self.docker_image,
                    self.env_mode,
                )
                container = self.client.containers.run(**run_kwargs)
                self.container_id = container.id
            timing['container_setup_elapsed_seconds'] = round(time.perf_counter() - container_setup_start, 3)

            # 3. 将代码和 requirements 放入容器内
            import tarfile
            import io
            pw_tarstream = io.BytesIO()
            with tarfile.open(fileobj=pw_tarstream, mode='w') as tar:
                # 写入 code.py
                code_content = code_str or ""
                if not code_content and self.code_path and os.path.exists(self.code_path):
                    with open(self.code_path, 'r', encoding='utf-8') as f:
                        code_content = f.read()
                if code_content:
                    c_bytes = code_content.encode('utf-8')
                    tinfo = tarfile.TarInfo(name='code.py')
                    tinfo.size = len(c_bytes)
                    tar.addfile(tinfo, io.BytesIO(c_bytes))

                # 写入 requirements.txt
                req_content = requirements_str or ""
                if not req_content and self.requirements_path and os.path.exists(self.requirements_path):
                    with open(self.requirements_path, 'r', encoding='utf-8') as f:
                        req_content = f.read()
                if req_content:
                    r_bytes = req_content.encode('utf-8')
                    tinfo = tarfile.TarInfo(name='requirements.txt')
                    tinfo.size = len(r_bytes)
                    tar.addfile(tinfo, io.BytesIO(r_bytes))

            pw_tarstream.seek(0)
            push_start = time.perf_counter()
            self.logger.info(f"正在将最新代码与依赖推送到容器 {self.container_id[:12]} 内的 /app 目录...")
            container.put_archive('/app', pw_tarstream.read())
            timing['push_to_container_elapsed_seconds'] = round(time.perf_counter() - push_start, 3)

            user = f"{os.getuid()}:{os.getgid()}" if hasattr(os, 'getuid') else ''
            deps_dir = _ENV_CACHE_DEPS_DIR if self.env_cache_volume else "/tmp/mas_runtime/.mas_pydeps"

            combined_logs: list[str] = []

            preflight = self._run_required_inputs_preflight(
                container,
                user=user,
                timeout=_remaining_timeout_seconds(),
            )
            if preflight is not None:
                combined_logs.append(preflight.get("output", ""))
                output_files = _collect_output_files()
                timing['total_elapsed_seconds'] = round(time.perf_counter() - total_start, 3)
                return {
                    'success': False,
                    'error': preflight.get('error', '挂载预检查失败'),
                    'output': "\n".join(part for part in combined_logs if part).strip(),
                    'files': output_files,
                    'container_id': self.container_id,
                    'exit_code': preflight.get('exit_code'),
                    'timing': timing,
                }

            # 4. 先安装依赖（如果 requirements.txt 非空），再执行代码，分别计时
            req_content = requirements_str or ""
            if not req_content and self.requirements_path and os.path.exists(self.requirements_path):
                with open(self.requirements_path, 'r', encoding='utf-8') as f:
                    req_content = f.read()

            if req_content.strip():
                timing['pip_install_skipped'] = False
                install_script = self._runtime_env_prelude(deps_dir) + (
                    'python -m pip install -i https://pypi.tuna.tsinghua.edu.cn/simple '
                    '-r /app/requirements.txt --target "$DEPS"'
                )
                self.logger.info(f"正在容器 {self.container_id[:12]} 内安装 Python 依赖...")
                install_result = _exec_phase(install_script, "pip_install")
                timing['pip_install_elapsed_seconds'] = round(install_result['elapsed_seconds'], 3)
                timing['install_exit_code'] = install_result['exit_code']
                if install_result['output']:
                    combined_logs.append(install_result['output'])

                if install_result['exit_code'] == 124:
                    output_files = _collect_output_files()
                    timing['total_elapsed_seconds'] = round(time.perf_counter() - total_start, 3)
                    return {
                        'success': False,
                        'error': f'依赖安装超时 (总超时限制 {timeout}s)',
                        'output': "\n".join(combined_logs),
                        'files': output_files,
                        'container_id': self.container_id,
                        'exit_code': install_result['exit_code'],
                        'timing': timing,
                    }

                if install_result['exit_code'] != 0:
                    output_files = _collect_output_files()
                    install_err = self._extract_error_summary(install_result['output'])
                    if not install_err:
                        install_err = f"依赖安装失败，退出码: {install_result['exit_code']}"
                    timing['total_elapsed_seconds'] = round(time.perf_counter() - total_start, 3)
                    return {
                        'success': False,
                        'error': install_err,
                        'output': "\n".join(combined_logs),
                        'files': output_files,
                        'container_id': self.container_id,
                        'exit_code': install_result['exit_code'],
                        'timing': timing,
                    }

            exec_script = self._runtime_env_prelude(deps_dir) + 'export PYTHONPATH="$DEPS:$PYTHONPATH"\npython /app/code.py'
            self.logger.info(f"正在容器 {self.container_id[:12]} 内执行 Python 代码 (超时设置: {timeout}s)...")
            python_result = _exec_phase(exec_script, "python_exec")
            timing['python_exec_elapsed_seconds'] = round(python_result['elapsed_seconds'], 3)
            timing['python_exit_code'] = python_result['exit_code']
            if python_result['output']:
                combined_logs.append(python_result['output'])

            status_code = python_result['exit_code']
            logs = "\n".join(part for part in combined_logs if part).strip()

            if status_code == 124: # timeout 命令的超时退出码
                output_files = _collect_output_files()
                timing['total_elapsed_seconds'] = round(time.perf_counter() - total_start, 3)
                return {
                    'success': False,
                    'error': f'执行超时 ({timeout}s)',
                    'output': logs,
                    'files': output_files,
                    'container_id': self.container_id,
                    'timing': timing,
                }

            output_files = _collect_output_files()

            if status_code != 0:
                # 提取去除了 pip 等前置噪音的日志进行错误分析
                exec_logs = logs.split("===MAS_EXEC_START===")[-1] if "===MAS_EXEC_START===" in logs else logs
                err = self._extract_error_summary(exec_logs)
                if not err:
                    err = f"容器退出码非0: {status_code}"
                timing['total_elapsed_seconds'] = round(time.perf_counter() - total_start, 3)
                return {
                    'success': False,
                    'error': err,
                    'output': logs,
                    'files': output_files,
                    'container_id': self.container_id,
                    'exit_code': status_code,
                    'timing': timing,
                }

            timing['total_elapsed_seconds'] = round(time.perf_counter() - total_start, 3)
            return {
                'success': True,
                'output': logs,
                'files': output_files,
                'container_id': self.container_id,
                'exit_code': status_code,
                'timing': timing,
            }
        except Exception as e:
            self.logger.error(f"执行失败: {e}")
            timing['total_elapsed_seconds'] = round(time.perf_counter() - total_start, 3)
            return {
                'success': False,
                'error': str(e),
                'output': '',
                'files': [],
                'container_id': self.container_id,
                'timing': timing,
            }
