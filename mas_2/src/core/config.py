"""
核心配置模块
负责加载环境变量和定义全局常量
"""
import os
from pathlib import Path

from dotenv import load_dotenv

# 始终从 mas_2 项目根目录加载 .env（不依赖进程 cwd，避免 notebook 在子目录、
# 或从其它路径启动 chainlit 时漏载导致与终端/Jupyter 行为不一致）
_MAS_ROOT = Path(__file__).resolve().parents[2]
load_dotenv(_MAS_ROOT / ".env", override=True)


def _env_bool(name: str, default: bool = True) -> bool:
    v = os.getenv(name)
    if v is None or str(v).strip() == "":
        return default
    return str(v).strip().lower() not in ("0", "false", "no", "off")

def _first_nonempty_env(*names: str) -> str | None:
    """按顺序取第一个非空环境变量；值会 strip，避免密钥首尾空格导致鉴权失败。"""
    for name in names:
        v = os.getenv(name)
        if v is not None and str(v).strip():
            return str(v).strip()
    return None


class Config:
    """全局配置类"""
    
    API_KEY = _first_nonempty_env("OPENAI_API_KEY", "DASHSCOPE_API_KEY")
    BASE_URL = os.getenv("BASE_URL", "https://dashscope.aliyuncs.com/compatible-mode/v1")
    
    # 模型配置
    MODEL_NAME = os.getenv("MODEL_NAME", "qwen-turbo")
    DEFAULT_TEMPERATURE = float(os.getenv("TEMPERATURE", "0.5"))
    MIMO_THINKING_TYPE = os.getenv("MIMO_THINKING_TYPE", "").strip().lower()

    # LLM HTTP：httpx 默认会读取 HTTP_PROXY/HTTPS_PROXY；若代理未开会导致 Connection error，
    # 而 curl 未走代理时仍可能成功。设为 false 可强制直连（与多数「关掉系统代理」效果类似）。
    HTTPS_TRUST_ENV = _env_bool("MAS_HTTPS_TRUST_ENV", True)
    LLM_REQUEST_TIMEOUT = float(os.getenv("LLM_REQUEST_TIMEOUT", "120"))
    
    # LangChain 配置
    LANGCHAIN_TRACING_V2 = os.getenv("LANGCHAIN_TRACING_V2", "false").lower()
    LANGCHAIN_API_KEY = os.getenv("LANGCHAIN_API_KEY")
    LANGCHAIN_PROJECT = os.getenv("LANGCHAIN_PROJECT", "mas_2")


# 全局配置实例
config = Config()

