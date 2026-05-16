from src.sandbox.tools import build_code_dev_sandbox_tool_bundle


def _build_state() -> dict:
    return {
        "current_step_skill_id": "scrnaseq-scanpy-core-analysis",
        "data_path": "./data",
        "result_path": "./result/sandbox_unit",
        "current_step_file_paths": {
            "input_files": ["README.md"],
            "output_files": ["./result/sandbox_unit/out.txt"],
        },
    }


def test_build_code_dev_sandbox_tool_bundle_exposes_expected_tools():
    bundle = build_code_dev_sandbox_tool_bundle(_build_state())

    assert bundle.backend == "local"
    assert [tool.name for tool in bundle.tools] == [
        "read_file",
        "list_files",
        "glob_search",
        "grep_text",
        "file_exists",
        "exec_command",
    ]
    assert any(root.endswith("/mas_2") for root in bundle.allowed_roots)


def test_sandbox_tool_bundle_can_read_and_list_files():
    bundle = build_code_dev_sandbox_tool_bundle(_build_state())

    listing = bundle.tool_map["list_files"].invoke({"path": ".", "recursive": False, "max_entries": 5})
    readme = bundle.tool_map["read_file"].invoke({"path": "README.md", "line_offset": 1, "max_lines": 3})
    exists = bundle.tool_map["file_exists"].invoke({"path": "README.md"})

    assert "Directory:" in listing
    assert "README.md" in readme
    assert "exists=True" in exists


def test_sandbox_tool_bundle_can_search_and_exec():
    bundle = build_code_dev_sandbox_tool_bundle(_build_state())

    glob_result = bundle.tool_map["glob_search"].invoke(
        {"pattern": "src/sandbox/*.py", "path": ".", "max_entries": 10}
    )
    grep_result = bundle.tool_map["grep_text"].invoke(
        {"pattern": "SandboxClient", "path": "src/sandbox", "head_limit": 20}
    )
    exec_result = bundle.tool_map["exec_command"].invoke(
        {"argv": ["pwd"], "cwd": ".", "timeout_s": 5}
    )

    assert "src/sandbox/client.py" in glob_result
    assert "SandboxClient" in grep_result
    assert "/mas_2" in exec_result
