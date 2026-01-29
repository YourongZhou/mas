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
    
    def __init__(self, docker_path: str, data_dir: str = None, output_dir: str = None):
        """
        初始化代码执行器
        
        Args:
            docker_path: 包含 code.py 和 requirements.txt 的目录路径
            data_dir: 数据目录路径（将被挂载到容器的 /app/data）
            output_dir: 输出目录路径（将被挂载到容器的 /app/output）
        """
        self.docker_available = self._check_docker_availability()
        self.code_path = f"{docker_path}/code.py"
        self.requirements_path = f"{docker_path}/requirements.txt"
        self.data_dir = data_dir
        self.output_dir = output_dir if output_dir else "/tmp/output"  # 默认输出目录
        self.docker_image = "python:3.13-slim"
        self.logger = logging.getLogger(__name__)
        
        # 处理数据路径：如果是文件，使用文件所在目录；如果是目录，直接使用
        self.data_mount_path = None
        self.data_file_name = None
        if self.data_dir and os.path.exists(self.data_dir):
            if os.path.isfile(self.data_dir):
                # 如果是文件，挂载文件所在的目录，并记录文件名
                self.data_mount_path = os.path.dirname(os.path.abspath(self.data_dir))
                self.data_file_name = os.path.basename(self.data_dir)
                self.logger.info(f"数据路径是文件，将挂载目录: {self.data_mount_path}, 文件名: {self.data_file_name}")
            elif os.path.isdir(self.data_dir):
                # 如果是目录，直接挂载目录
                self.data_mount_path = os.path.abspath(self.data_dir)
                self.logger.info(f"数据路径是目录，将挂载: {self.data_mount_path}")

        # 设置卷挂载
        self.volume_mounts = {}
        # 挂载代码文件 (将宿主机的 code.py 挂载到容器的 /app/code.py)
        if os.path.exists(self.code_path):
            self.volume_mounts[self.code_path] = {'bind': '/app/code.py', 'mode': 'ro'}

        # 挂载数据目录
        if self.data_dir and os.path.exists(self.data_dir):
            data_mount_path = os.path.abspath(self.data_dir)
            if os.path.isfile(data_mount_path):
                data_mount_path = os.path.dirname(data_mount_path)
            self.volume_mounts[data_mount_path] = {'bind': '/app/data', 'mode': 'ro'}

        # 挂载输出目录
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)
        self.volume_mounts[os.path.abspath(self.output_dir)] = {'bind': '/app/output', 'mode': 'rw'}

    def _check_docker_availability(self) -> bool:
        """检查 Docker 是否可用"""
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
        """创建 Dockerfile"""
        # 构建 Dockerfile
        dockerfile_lines = [
            f"FROM {self.docker_image}",
            "WORKDIR /app",
            "RUN mkdir -p /app/data /app/output"
        ]

        # 只有当 requirements 变化时，build 才会耗时
        if os.path.exists(self.requirements_path):
            # 临时复制一份用于 build
            shutil.copy2(self.requirements_path, os.path.join(temp_dir, "requirements.txt"))
            dockerfile_lines.extend([
                "COPY requirements.txt .",
                "RUN pip install --no-cache-dir -r requirements.txt",
            ])


        # To zhenwei：我改了一下，数据通过 volume mount 挂载，不需要在 Dockerfile 中 COPY
        # 这样可以避免数据重复，并且支持大文件

        # 设置执行命令
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
                if self.data_mount_path:
                    self.logger.info(f"数据目录挂载: {self.data_mount_path} -> /app/data")
                    if self.data_file_name:
                        self.logger.info(f"数据文件名: {self.data_file_name} (在容器内路径: /app/data/{self.data_file_name})")

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

