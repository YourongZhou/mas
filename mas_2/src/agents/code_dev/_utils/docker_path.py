"""
Docker 路径转换工具
将本地文件路径转换为 Docker 容器内的路径
"""
import os


def convert_to_docker_path(local_path: str, mode: str = "data") -> str:
    """
    将本地文件路径转换为Docker容器内的路径

    Args:
        local_path: 本地文件路径，例如 '/home/user/Data/scHetero/adata_raw_new.h5ad'
        mode: 模式，'data' 表示数据路径，'output' 表示输出路径

    Returns:
        Docker容器内的路径
    """
    if mode == "data":
        base_docker_path = "/app/data"
    elif mode == "output":
        base_docker_path = "/app/output"
    else:
        raise ValueError("mode参数必须为'data'或'output'")
    
    # 检查路径是否为目录
    if os.path.isdir(local_path):
        # 如果是目录，整个目录挂载到 base_docker_path
        return base_docker_path
    else:
        # 如果是文件，将文件名追加到 base_docker_path
        filename = os.path.basename(local_path)
        docker_path = os.path.join(base_docker_path, filename)
        return docker_path

