import os
from langchain_core.tools import tool

@tool
def read_local_file(file_path: str, offset: int = 0, limit: int = 4000) -> str:
    """
    读取指定路径的文本文件内容（仅限只读操作，不可用于读取二进制文件）。
    用于查看当前系统里的本地脚本（如 workflows/*/scripts/*.py）或者长篇 markdown 文档。
    
    Args:
        file_path (str): 需要读取的文件的绝对路径（物理机上的绝对路径，如果你只知道相对路径需要先推测出绝对路径）。
        offset (int): 从哪个字符开始读取（支持大文件分块读取）。默认为 0。
        limit (int): 最大读取的字符数。默认为 4000，如果文件过长建议通过 offset 多次获取。

    Returns:
        str: 读取到的文件文本内容，或者在读取失败时的错误信息。如果内容被截断，会带有提示。
    """
    if not os.path.exists(file_path):
        return f"Error: 文件 {file_path} 不存在。"
    
    if not os.path.isfile(file_path):
        return f"Error: 路径 {file_path} 不是一个常规文件。"
    
    try:
        # 支持处理不同编码，首先尝试 utf-8
        try:
            with open(file_path, "r", encoding="utf-8") as f:
                f.seek(offset)
                content = f.read(limit)
                
                current_pos = f.tell()
                f.seek(0, os.SEEK_END)
                total_size = f.tell()
        except UnicodeDecodeError:
            with open(file_path, "r", encoding="gbk") as f:
                f.seek(offset)
                content = f.read(limit)
                
                current_pos = f.tell()
                f.seek(0, os.SEEK_END)
                total_size = f.tell()

        if current_pos < total_size:
            content += f"\n\n[提示: 文件未读取完。已读取到了 {current_pos}/{total_size} 字节。你可以增加 offset 继续读取之后的内容。]"
            
        return content
    except Exception as e:
        return f"Error: 读取文件时发生异常: {str(e)}"

@tool
def list_directory(dir_path: str) -> str:
    """
    列出指定目录下的文件和文件夹结构，用于了解可用脚本和资源。
    用于浏览本地技能路径(workflow)内部存在哪些文件（如 scripts/ 等子目录）。
    
    Args:
        dir_path (str): 需要检查的目录的绝对路径。
        
    Returns:
        str: 该目录下文件和文件夹的列表信息。
    """
    if not os.path.exists(dir_path):
        return f"Error: 目录 {dir_path} 不存在。"
    
    if not os.path.isdir(dir_path):
        return f"Error: 路径 {dir_path} 不是一个目录。"
    
    try:
        items = os.listdir(dir_path)
        if not items:
            return "该目录为空。"
        
        res = f"目录 {dir_path} 内容:\n"
        for item in sorted(items):
            full_path = os.path.join(dir_path, item)
            item_type = "DIR " if os.path.isdir(full_path) else "FILE"
            if os.path.isfile(full_path):
                size = os.path.getsize(full_path)
                res += f"[{item_type}] {item} (大小: {size} 字节)\n"
            else:
                res += f"[{item_type}] {item}\n"
        return res
    except Exception as e:
        return f"Error: 读取目录时发生异常: {str(e)}"
