# 测试说明

## 目录结构

```
tests/
├── conftest.py              # Pytest 配置和 fixtures
├── integration/             # 集成测试
│   └── test_code_agent.py  # Code Dev Agent 集成测试
└── unit/                    # 单元测试
    └── test_supervisor.py  # Supervisor Agent 单元测试
```

## 运行测试

### 运行所有测试
```bash
pytest
```

### 运行单元测试（使用 Mock，不需要真实 API）
```bash
pytest tests/unit/
```

### 运行集成测试（需要真实的 LLM API，需要配置 .env）
```bash
pytest tests/integration/
```

### 运行特定测试文件
```bash
pytest tests/unit/test_supervisor.py
pytest tests/integration/test_code_agent.py
```

### 运行特定测试函数
```bash
pytest tests/unit/test_supervisor.py::test_supervisor_decision_with_mock
```

### 显示详细输出
```bash
pytest -v
```

### 显示打印输出
```bash
pytest -s
```

## 前置条件

### 单元测试
- 无需额外配置，使用 Mock 模拟 LLM 调用

### 集成测试
- 需要配置 `.env` 文件，包含：
  ```
  OPENAI_API_KEY=your_api_key
  BASE_URL=https://dashscope.aliyuncs.com/compatible-mode/v1
  MODEL_NAME=qwen-plus
  ```
- 需要网络连接以调用 LLM API

## Fixtures

### `basic_global_state`
基础的 GlobalState 字典，包含所有必需字段的默认值。

### `code_dev_state`
专门用于 Code Dev Agent 测试的状态，包含 Code Dev 特有的字段。

### `supervisor_state`
专门用于 Supervisor Agent 测试的状态。

## 编写新测试

### 单元测试示例
```python
from unittest.mock import patch, MagicMock

def test_my_agent(supervisor_state):
    with patch('src.agents.my_agent.graph.llm') as mock_llm:
        # 设置 Mock
        mock_llm.invoke.return_value = MagicMock(content="mocked response")
        
        # 执行测试
        result = my_agent_graph.invoke(supervisor_state)
        
        # 断言
        assert result["some_field"] == "expected_value"
```

### 集成测试示例
```python
def test_my_agent_integration(basic_global_state):
    state = basic_global_state.copy()
    state["user_query"] = "测试查询"
    
    result = my_agent_graph.invoke(state)
    
    assert result.get("pending_contribution") is not None
```

