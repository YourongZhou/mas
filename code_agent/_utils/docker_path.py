import os

def convert_to_docker_path(local_path, mode = "data"):
    """
    将本地文件路径转换为Docker容器内的路径

    Args:
        local_path (str): 本地文件路径，例如 '/home/zwhuang/Data/scHetero/adata_raw_new.h5ad'
        base_docker_path (str): Docker容器内的基础路径，默认为 '/app/data'

    Returns:
        str: Docker容器内的路径
    """
    if mode == "data":
        base_docker_path = "/app/data"
    elif mode == "output":
        base_docker_path = "/app/output"
    else:
        raise ValueError("mode参数必须为'data'或'output'")
    # Check if the path is a directory
    if os.path.isdir(local_path):
        # If it's a directory, the whole directory is mounted to base_docker_path
        return base_docker_path
    else:
        # If it's a file, append the filename to base_docker_path
        filename = os.path.basename(local_path)
        docker_path = os.path.join(base_docker_path, filename)
        return docker_path