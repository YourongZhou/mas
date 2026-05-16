import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from src.main import finalize_step


def test_finalize_step_always_removes_task_container(basic_global_state):
    state = basic_global_state.copy()
    state["docker_container_id"] = "cid-1234567890ab"
    state["docker_env_mode"] = "shared"

    fake_container = MagicMock()
    fake_client = SimpleNamespace(
        containers=SimpleNamespace(get=lambda container_id: fake_container)
    )
    fake_docker_module = SimpleNamespace(from_env=lambda: fake_client)
    fake_llm = MagicMock()
    fake_llm.invoke.return_value = SimpleNamespace(content="final answer")

    with patch("src.main.get_llm", return_value=fake_llm), patch.dict(
        sys.modules, {"docker": fake_docker_module}
    ):
        result = finalize_step(state)

    fake_container.remove.assert_called_once_with(force=True)
    assert result["docker_container_id"] is None
    assert result["final_answer"] == "final answer"
