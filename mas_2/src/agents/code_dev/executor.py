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
from pathlib import Path

logging.basicConfig(level=logging.INFO)


class CodeExecutor:
    """使用 Docker 容器执行代码的执行器"""
    
    def __init__(self, docker_path: str, data_dir: str = None, data_dirs: list = None,
                 output_dir: str = None, input_files: list = None,
                 workflow_host_path: str = None):
        """
        初始化代码执行器

        workflow_host_path: 若提供，将宿主目录只读挂载到容器内 /app/workflow（供 SKILL 内 scripts 导入）。
        """
        self.logger = logging.getLogger(__name__)
        self.docker_available = self._check_docker_availability()
        self.code_path = f"{docker_path}/code.py"
        self.requirements_path = f"{docker_path}/requirements.txt"
        self.output_dir = output_dir if output_dir else "/tmp/output"  # 默认输出目录
        self.docker_image = "python:3.13-slim"
        
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
        
        # 挂载代码文件 (将宿主机的 code.py 挂载到容器的 /app/code.py)
        if os.path.exists(self.code_path):
            # 格式: host_path:container_path:mode
            mount_str = f"{os.path.abspath(self.code_path)}:/app/code.py:ro"
            self.volume_mounts.append(mount_str)

        # 挂载 requirements 文件（如存在）
        if os.path.exists(self.requirements_path):
            mount_str = f"{os.path.abspath(self.requirements_path)}:/app/requirements.txt:ro"
            self.volume_mounts.append(mount_str)

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
            self.volume_mounts.append(f"{wf}:/app/workflow:ro")
            self.logger.info(f"Workflow 目录挂载: {wf} -> /app/workflow")
    
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

    def _prepare_temp_directory(self, temp_dir: str) -> None:
        shutil.copy2(self.code_path, os.path.join(temp_dir, 'code.py'))
        if os.path.exists(self.requirements_path):
            shutil.copy2(self.requirements_path, os.path.join(temp_dir, 'requirements.txt'))

    def execute(self,
                environment_vars: dict | None = None,
                mem_limit: str | None = '4g',
                timeout: int | None = 300) -> dict:
        """
        使用卷挂载执行代码
        """
        if not self.docker_available:
            return {
                'success': False,
                'error': 'Docker不可用',
                'output': '',
                'files': []
            }

        with tempfile.TemporaryDirectory() as temp_dir:
            container = None
            try:
                self._prepare_temp_directory(temp_dir)

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

                # set -e：pip 失败则不再执行 python，避免「装包失败仍跑代码」导致 ModuleNotFoundError 误导。
                # 运行时依赖与缓存写入容器 /tmp，避免污染宿主输出目录并减少文件监听噪声。
                _container_script = """set -e
export PYTHONPATH="/app/workflow${PYTHONPATH:+:$PYTHONPATH}"
RUNTIME_ROOT=/tmp/mas_runtime
DEPS="$RUNTIME_ROOT/.mas_pydeps"
export TMPDIR="$RUNTIME_ROOT/.mas_tmp"
mkdir -p "$TMPDIR" "$DEPS" "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR" "$PYTHONPYCACHEPREFIX"
if [ -f /app/requirements.txt ] && [ -s /app/requirements.txt ]; then
  python -m pip install --no-cache-dir --upgrade -r /app/requirements.txt --target "$DEPS"
fi
export PYTHONPATH="$DEPS:$PYTHONPATH"
python /app/code.py"""
                command = "bash -lc " + shlex.quote(_container_script)

                self.logger.info(f"运行容器，镜像: {self.docker_image}")
                
                run_kwargs = {
                    'image': self.docker_image,
                    'command': command,
                    'volumes': self.volume_mounts,  # <--- 直接传列表
                    'environment': env_vars,
                    'mem_limit': mem_limit,
                    'network_mode': 'bridge',
                    'detach': True,
                    'auto_remove': False,
                }
                if hasattr(os, 'getuid') and hasattr(os, 'getgid'):
                    try:
                        run_kwargs['user'] = f"{os.getuid()}:{os.getgid()}"
                    except Exception as e:
                        self.logger.warning(f"无法获取当前用户 UID/GID，将使用 Docker 默认用户: {e}")

                container = self.client.containers.run(**run_kwargs)

                try:
                    container.wait(timeout=timeout)
                except Exception as e:
                    self.logger.warning(f"容器等待超时或出错: {e}")
                    try:
                        container.stop(timeout=10)
                    except Exception as stop_err:
                        self.logger.warning(f"容器停止失败: {stop_err}")
                    return {
                        'success': False,
                        'error': f'执行超时: {e}',
                        'output': '',
                        'files': []
                    }

                logs = container.logs().decode('utf-8').strip()

                output_files = []
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

                return {
                    'success': True,
                    'output': logs,
                    'files': output_files,
                    'container_id': container.id[:12]
                }
            except Exception as e:
                self.logger.error(f"执行失败: {e}")
                return {
                    'success': False,
                    'error': str(e),
                    'output': '',
                    'files': []
                }
            finally:
                if container is not None:
                    try:
                        container.remove(force=True)
                    except Exception as remove_container_err:
                        self.logger.warning(f"容器清理失败: {remove_container_err}")
