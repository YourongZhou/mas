"""
代码执行器模块
使用 Docker 容器安全执行生成的代码
"""
import os
import tempfile
import docker
import shutil
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)


class CodeExecutor:
    """使用 Docker 容器执行代码的执行器"""
    
    def __init__(self, docker_path: str, data_dir: str = None, data_dirs: list = None, 
                 output_dir: str = None, input_files: list = None):
        """
        初始化代码执行器
        """
        self.docker_available = self._check_docker_availability()
        self.code_path = f"{docker_path}/code.py"
        self.requirements_path = f"{docker_path}/requirements.txt"
        self.output_dir = output_dir if output_dir else "/tmp/output"  # 默认输出目录
        self.docker_image = "python:3.13-slim"
        self.logger = logging.getLogger(__name__)
        
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

    def _create_dockerfile(self, temp_dir: str, output_mount_path: str | None = None):
        dockerfile_lines = [
            f"FROM {self.docker_image}",
            "WORKDIR /app",
            "RUN mkdir -p /app/data /app/output"
        ]

        if os.path.exists(self.requirements_path):
            shutil.copy2(self.requirements_path, os.path.join(temp_dir, "requirements.txt"))
            dockerfile_lines.extend([
                "COPY requirements.txt .",
                "RUN pip install --no-cache-dir -r requirements.txt",
            ])

        dockerfile_lines.append('CMD ["python", "code.py"]')

        dockerfile_content = "\n".join(dockerfile_lines)
        dockerfile_path = os.path.join(temp_dir, 'Dockerfile')
        with open(dockerfile_path, 'w') as f:
            f.write(dockerfile_content)

        return dockerfile_path

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
            try:
                self._prepare_temp_directory(temp_dir)
                import uuid
                image_tag = f"code-executor-{uuid.uuid4().hex[:8]}"

                self.logger.info(f"构建Docker镜像: {image_tag}")
                self._create_dockerfile(temp_dir)

                self.client.images.build(
                    path=temp_dir,
                    tag=image_tag,
                    rm=True,
                    forcerm=True
                )

                # --- 关键修改 4：不需要解析 volumes 字典，直接使用列表 ---
                # self.volume_mounts 已经是符合 Docker SDK 格式的列表了

                env_vars = {
                    'PYTHONUNBUFFERED': '1',
                    'DEBIAN_FRONTEND': 'noninteractive',
                    'MPLCONFIGDIR': '/tmp/matplotlib',
                    'NUMBA_CACHE_DIR': '/tmp/numba',
                    'HOME': '/tmp',
                    'PYTHONPYCACHEPREFIX': '/tmp/__pycache__',
                }
                if environment_vars:
                    env_vars.update(environment_vars)

                self.logger.info(f"运行容器，镜像: {image_tag}")
                
                user = f"{os.getuid()}:{os.getgid()}"
                
                container = self.client.containers.run(
                    image=image_tag,
                    volumes=self.volume_mounts,  # <--- 直接传列表
                    environment=env_vars,
                    user=user,
                    mem_limit=mem_limit,
                    network_mode='bridge',
                    detach=True,
                    auto_remove=False
                )

                try:
                    container.wait(timeout=timeout)
                except Exception as e:
                    self.logger.warning(f"容器等待超时或出错: {e}")
                    container.stop(timeout=10)
                    container.remove(force=True)
                    self.client.images.remove(image=image_tag, force=True)
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
                            output_files.append({
                                'path': str(file_path),
                                'name': file_path.name,
                                'size': file_path.stat().st_size,
                                'size_mb': file_path.stat().st_size / (1024 * 1024)
                            })

                container.remove(force=True)
                self.client.images.remove(image=image_tag, force=True)

                return {
                    'success': True,
                    'output': logs,
                    'files': output_files,
                    'image_tag': image_tag,
                    'container_id': container.id[:12]
                }

            except docker.errors.BuildError as e:
                self.logger.error(f"镜像构建失败: {e}")
                return {
                    'success': False,
                    'error': f'镜像构建失败: {e}',
                    'output': '',
                    'files': []
                }
            except Exception as e:
                self.logger.error(f"执行失败: {e}")
                return {
                    'success': False,
                    'error': str(e),
                    'output': '',
                    'files': []
                }