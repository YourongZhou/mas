from types import SimpleNamespace

from src.agents.code_dev.executor import CodeExecutor


class _FakeContainer:
    def __init__(self, container_id: str = "cid-123", preflight_ok: bool = True):
        self.id = container_id
        self.status = "running"
        self.attrs = {"Config": {"Labels": {}}}
        self.preflight_ok = preflight_ok
        self.exec_commands: list[str] = []
        self.put_archive_calls = 0

    def start(self):
        self.status = "running"

    def put_archive(self, path, data):
        self.put_archive_calls += 1

    def exec_run(self, cmd, user="", demux=False):
        self.exec_commands.append(cmd)
        if "===PRE_EXEC_PREFLIGHT===" in cmd:
            if self.preflight_ok:
                return SimpleNamespace(exit_code=0, output=b"PREFLIGHT_OK\n")
            return SimpleNamespace(
                exit_code=31,
                output=b"MISSING_INPUT::/app/data/example.txt\n",
            )
        return SimpleNamespace(
            exit_code=0,
            output=b"===MAS_EXEC_START===\n===RESULT===ok===\n",
        )

    def remove(self, force=False):
        return None


class _FakeContainersAPI:
    def __init__(self, existing=None):
        self.existing = existing or {}
        self.run_calls = 0
        self.last_run_kwargs = None

    def get(self, container_id):
        if container_id not in self.existing:
            raise Exception(f"missing container: {container_id}")
        return self.existing[container_id]

    def run(self, **kwargs):
        self.run_calls += 1
        self.last_run_kwargs = kwargs
        container = _FakeContainer(container_id=f"new-{self.run_calls}")
        self.existing[container.id] = container
        return container


class _FakeVolumesAPI:
    def __init__(self):
        self.created = []

    def get(self, name):
        raise Exception(name)

    def create(self, name):
        self.created.append(name)
        return SimpleNamespace(name=name)


class _FakeDockerClient:
    def __init__(self, existing=None):
        self.containers = _FakeContainersAPI(existing=existing)
        self.volumes = _FakeVolumesAPI()
        self.images = SimpleNamespace(get=lambda image: True)


def _patch_docker(monkeypatch, client):
    def _fake_check(self):
        self.client = client
        return True

    monkeypatch.setattr(CodeExecutor, "_check_docker_availability", _fake_check)


def test_executor_reuses_same_experiment_container(monkeypatch, tmp_path):
    existing = {"cid-123": _FakeContainer(container_id="cid-123")}
    client = _FakeDockerClient(existing=existing)
    _patch_docker(monkeypatch, client)

    executor = CodeExecutor(
        docker_path=None,
        container_id="cid-123",
        output_dir=str(tmp_path / "output"),
    )
    result = executor.execute(code_str="print('x')", requirements_str="")

    assert result["success"] is True
    assert result["container_id"] == "cid-123"
    assert client.containers.run_calls == 0


def test_executor_does_not_reuse_cross_experiment_container_by_signature(monkeypatch, tmp_path):
    client = _FakeDockerClient()
    _patch_docker(monkeypatch, client)

    executor = CodeExecutor(
        docker_path=None,
        container_id=None,
        output_dir=str(tmp_path / "output"),
        docker_image="mas/py-scverse:v1",
        runtime="python",
        env_profile="py-scverse-v1",
        env_signature="sig-123",
        env_mode="shared",
    )
    result = executor.execute(code_str="print('x')", requirements_str="")

    assert result["success"] is True
    assert client.containers.run_calls == 1
    assert result["container_id"].startswith("new-")


def test_executor_preflight_fails_when_required_input_missing(monkeypatch, tmp_path):
    existing = {"cid-preflight": _FakeContainer(container_id="cid-preflight", preflight_ok=False)}
    client = _FakeDockerClient(existing=existing)
    _patch_docker(monkeypatch, client)

    executor = CodeExecutor(
        docker_path=None,
        container_id="cid-preflight",
        output_dir=str(tmp_path / "output"),
        required_input_container_paths=["/app/data/example.txt"],
    )
    result = executor.execute(code_str="print('x')", requirements_str="")

    assert result["success"] is False
    assert "缺少输入文件" in result["error"]
