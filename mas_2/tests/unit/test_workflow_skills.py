"""workflow_skills：挂载策略与注册表。"""
import pytest

from src.utils.execution_environment import (
    filter_profiled_requirements,
    get_image_catalog_entry,
    resolve_environment_for_state,
    resolve_execution_environment,
)
from src.utils.workflow_skills import (
    SCANPY_CORE_SKILL_ID,
    get_skill,
    list_skills,
    should_mount_workflow_in_docker,
    use_scanpy_code_style,
)


def test_use_scanpy_style_only_for_core_skill():
    assert use_scanpy_code_style(None) is False
    assert use_scanpy_code_style("") is False
    assert use_scanpy_code_style("cell-cell-communication") is False
    assert use_scanpy_code_style(SCANPY_CORE_SKILL_ID) is True


def test_should_mount_false_when_empty_or_unknown():
    assert should_mount_workflow_in_docker(None) is False
    assert should_mount_workflow_in_docker("") is False
    assert should_mount_workflow_in_docker("definitely-not-a-workflow-skill-xyz") is False


def test_should_mount_true_for_registered_skill():
    ids = [r.skill_id for r in list_skills()]
    if not ids:
        pytest.skip("无 workflows/*/SKILL.md")
    assert should_mount_workflow_in_docker(ids[0]) is True


def test_scanpy_skill_has_structured_environment_metadata():
    rec = get_skill(SCANPY_CORE_SKILL_ID)
    assert rec is not None
    assert rec.runtime == "python"
    assert rec.env_profile == "py-scverse-v1"
    assert rec.env_image == "mas/py-scverse:v1"
    assert rec.env_mode == "shared"


def test_resolve_execution_environment_for_scanpy_skill():
    env = resolve_execution_environment(SCANPY_CORE_SKILL_ID)
    assert env is not None
    assert env.runtime == "python"
    assert env.env_profile == "py-scverse-v1"
    assert env.env_image == "mas/py-scverse:v1"
    assert env.env_signature


def test_image_catalog_contains_scanpy_profile():
    entry = get_image_catalog_entry("py-scverse-v1")
    assert entry is not None
    assert entry.runtime == "python"
    assert entry.image_tag == "mas/py-scverse:v1"
    assert "openpyxl" in entry.allowed_runtime_extras


def test_resolve_environment_for_state_falls_back_to_selected_skill():
    env = resolve_environment_for_state(
        {
            "current_step_skill_id": None,
            "selected_skill_id": SCANPY_CORE_SKILL_ID,
        }
    )
    assert env is not None
    assert env.env_profile == "py-scverse-v1"


def test_profiled_requirements_filter_blocks_core_and_non_allowlisted_packages():
    env = resolve_execution_environment(SCANPY_CORE_SKILL_ID)
    assert env is not None

    filtered = filter_profiled_requirements(
        env,
        "\n".join(
            [
                "scanpy==1.10.1",
                "openpyxl==3.1.5",
                "requests==2.32.0",
                "pyarrow",
            ]
        ),
    )

    assert "scanpy==1.10.1" in filtered.blocked_core_lines
    assert "requests==2.32.0" in filtered.blocked_lines
    assert "openpyxl==3.1.5" in filtered.allowed_lines
    assert "pyarrow" in filtered.allowed_lines
    assert filtered.env_cache_key
    assert filtered.env_cache_volume.startswith("mas-env-cache-")


def test_profiled_requirements_cache_key_changes_with_allowed_extras():
    env = resolve_execution_environment(SCANPY_CORE_SKILL_ID)
    assert env is not None

    filtered_a = filter_profiled_requirements(env, "openpyxl==3.1.5")
    filtered_b = filter_profiled_requirements(env, "pyarrow")

    assert filtered_a.env_cache_key != filtered_b.env_cache_key
    assert filtered_a.env_cache_volume != filtered_b.env_cache_volume
