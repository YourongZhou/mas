import os
import tempfile
import docker
import shutil
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO)


class CodeExecutor:
    def __init__(self, docker_path: str, data_dir: str = None, output_dir: str = None):
        self.docker_available = self._check_docker_availability()
        self.code_path = f"{docker_path}/code.py"
        self.requirements_path = f"{docker_path}/requirements.txt"
        self.data_dir = data_dir
        self.output_dir = output_dir if output_dir else "/tmp/output"  # Default output directory
        self.docker_image = "python:3.13-slim"
        self.logger = logging.getLogger(__name__)

        # Set up volume mounts
        self.volume_mounts = {}
        if self.data_dir:
            self.volume_mounts[self.data_dir] = {
                'bind': '/app/data',
                'mode': 'ro'  # Read-only for data input
            }
        # Always mount output directory as writable
        self.volume_mounts[self.output_dir] = {
            'bind': '/app/output',
            'mode': 'rw'  # Read-write for output
        }



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

        # build Dockerfile
        dockerfile_lines = [
            f"FROM {self.docker_image}",
            "WORKDIR /app",
        ]

        # Add requirements installation
        if os.path.exists(self.requirements_path):
            dockerfile_lines.extend([
                f"COPY requirements.txt .",
                "RUN pip install --no-cache-dir -r requirements.txt",
            ])

        # Copy code file
        dockerfile_lines.append(f"COPY code.py .")

        # Copy data directory
        if hasattr(self, 'data_dir') and self.data_dir and os.path.exists(self.data_dir):
        # Copy data directory contents to container
            dest_data_path = os.path.join(temp_dir, 'data')
            if os.path.isdir(self.data_dir):
                shutil.copytree(self.data_dir, dest_data_path, dirs_exist_ok=True)
            else:
                # If it's a single file, copy it to a data directory
                os.makedirs(dest_data_path, exist_ok=True)
                shutil.copy2(self.data_dir, dest_data_path)

            dockerfile_lines.append("COPY data /app/data")

        # Ensure both data and output directories exist in the container
        dockerfile_lines.append("RUN mkdir -p /app/data")
        dockerfile_lines.append("RUN mkdir -p /app/output")

        # Set execution command
        dockerfile_lines.append('CMD ["python", "code.py"]')

        dockerfile_content = "\n".join(dockerfile_lines)
        dockerfile_path = os.path.join(temp_dir, 'Dockerfile')
        with open(dockerfile_path, 'w') as f:
            f.write(dockerfile_content)

        return dockerfile_path



    def _prepare_temp_directory(self, temp_dir: str) -> None:
        """准备临时目录内容"""
        # 复制代码文件
        shutil.copy2(self.code_path, os.path.join(temp_dir, 'code.py'))

        # 复制requirements.txt（如果存在）
        if os.path.exists(self.requirements_path):
            shutil.copy2(self.requirements_path, os.path.join(temp_dir, 'requirements.txt'))



    def execute(self,
        environment_vars: dict | None = None,
        mem_limit: str | None = '4g',
        timeout: int | None = 300) -> dict:
        """
        使用卷挂载执行代码

        Args:
            environment_vars: 环境变量
            mem_limit: 内存限制
            timeout: 执行超时时间（秒）

        Returns:
            包含执行结果的字典
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
                # 准备文件
                self._prepare_temp_directory(temp_dir)

                # 生成唯一的镜像标签
                import uuid
                image_tag = f"code-executor-{uuid.uuid4().hex[:8]}"

                # 构建Docker镜像
                self.logger.info(f"构建Docker镜像: {image_tag}")

                dockerfile_path = self._create_dockerfile(temp_dir)

                build_logs = self.client.images.build(
                    path=temp_dir,
                    tag=image_tag,
                    rm=True,
                    forcerm=True
                )

                # 解析卷挂载
                volumes = {}
                if self.volume_mounts:
                    for host_path, container_path in self.volume_mounts.items():
                        if isinstance(container_path, dict):
                            # 已经是Docker SDK格式
                            volumes[host_path] = container_path
                        else:
                            # 转换为Docker SDK格式
                            volumes[host_path] = {
                                'bind': container_path,
                                'mode': 'rw'
                            }

                # 准备环境变量
                env_vars = {
                    'PYTHONUNBUFFERED': '1',  # 实时输出
                    'DEBIAN_FRONTEND': 'noninteractive',  # 非交互模式
                    'MPLCONFIGDIR': '/tmp/matplotlib',  # Matplotlib缓存
                    'NUMBA_CACHE_DIR': '/tmp/numba',    # Numba缓存
                    'HOME': '/tmp',                     # HOME目录
                    'PYTHONPYCACHEPREFIX': '/tmp/__pycache__',  # Python缓存
                }
                if environment_vars:
                    env_vars.update(environment_vars)

                # 运行容器
                self.logger.info(f"运行容器，镜像: {image_tag}")

                user = f"{os.getuid()}:{os.getgid()}"
                print(f"Using user: {user}")

                container = self.client.containers.run(
                    image=image_tag,
                    volumes=volumes,
                    environment=env_vars,
                    user=user,
                    mem_limit=mem_limit,
                    network_mode='bridge',
                    detach=True,
                    auto_remove=False  # 手动控制删除
                )

                # 等待容器完成
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

                # 获取容器输出
                logs = container.logs().decode('utf-8').strip()

                # 收集输出文件信息
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

                # 清理
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