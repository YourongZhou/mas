"""
核心配置模块
负责加载环境变量和定义全局常量
"""
import os
from dotenv import load_dotenv

# 加载 .env 文件中的变量
load_dotenv()


class Config:
    """全局配置类"""
    
    # API 配置
    API_KEY = "sk-3e43aba7e80343bc96fb5e7d549837ac"#os.getenv("OPENAI_API_KEY")
    BASE_URL = os.getenv("BASE_URL", "https://dashscope.aliyuncs.com/compatible-mode/v1")
    
    # 模型配置
    MODEL_NAME = "qwen-plus"#os.getenv("MODEL_NAME", "qwen-plus")
    DEFAULT_TEMPERATURE = float(os.getenv("TEMPERATURE", "0.5"))
    
    # LangChain 配置
    LANGCHAIN_TRACING_V2 = "true"#os.getenv("LANGCHAIN_TRACING_V2", "false").lower() == "true"
    LANGCHAIN_API_KEY = "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"#os.getenv("LANGCHAIN_API_KEY")
    LANGCHAIN_PROJECT = "mas_2"#os.getenv("LANGCHAIN_PROJECT", "mas_2")


# 全局配置实例
config = Config()

