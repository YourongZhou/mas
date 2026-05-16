"""Structured execution environment metadata and requirement filtering."""
from __future__ import annotations

import hashlib
import json
import logging
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Optional

from src.utils.project_paths import get_mas2_project_root
from src.utils.workflow_skills import WorkflowSkillRecord, get_skill

_logger = logging.getLogger(__name__)

_CATALOG_PATH = get_mas2_project_root() / "docker" / "image_catalog.json"
_ENV_MANIFESTS_ROOT = get_mas2_project_root() / "docker" / "env_manifests"
_ASSET_CACHE_ROOT = get_mas2_project_root() / ".cache" / "workflow_assets"
_DEFAULT_ALLOWED_RUNTIME_EXTRAS = (
    "openpyxl",
    "xlrd",
    "pyarrow",
    "fastparquet",
    "tables",
    "xlsxwriter",
)
_ENV_CACHE_VOLUME_PREFIX = "mas-env-cache-"


def _normalize_items(values: Iterable[str]) -> tuple[str, ...]:
    out = []
    for item in values:
        text = str(item).strip()
        if text:
            out.append(text)
    return tuple(sorted(dict.fromkeys(out)))


def canonicalize_requirement_name(name: str) -> str:
    text = str(name or "").strip()
    if not text:
        return ""
    return re.sub(r"[-_.]+", "-", text).lower()


def requirement_name_from_line(line: str) -> str:
    text = str(line or "").strip()
    if not text or text.startswith("#") or text.startswith("-"):
        return ""
    match = re.match(r"^\s*([A-Za-z0-9_.-]+)", text)
    if not match:
        return ""
    return canonicalize_requirement_name(match.group(1))


def requirement_lines_from_text(text: str) -> list[str]:
    out: list[str] = []
    for raw_line in str(text or "").splitlines():
        line = raw_line.strip()
        if line and not line.startswith("#"):
            out.append(line)
    return out


def _compute_env_signature(
    *,
    runtime: str,
    env_profile: str,
    env_image: str,
    env_extra_requirements: Iterable[str],
    required_assets: Iterable[str],
) -> str:
    payload = {
        "runtime": runtime,
        "env_profile": env_profile,
        "env_image": env_image,
        "env_extra_requirements": list(_normalize_items(env_extra_requirements)),
        "required_assets": list(_normalize_items(required_assets)),
    }
    blob = json.dumps(payload, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:16]


def _compute_env_cache_key(env_signature: str, runtime_extras: Iterable[str]) -> str:
    payload = {
        "env_signature": env_signature,
        "runtime_extras": list(_normalize_items(runtime_extras)),
    }
    blob = json.dumps(payload, sort_keys=True, ensure_ascii=True).encode("utf-8")
    return hashlib.sha256(blob).hexdigest()[:20]


def env_cache_volume_name(cache_key: str) -> str:
    return f"{_ENV_CACHE_VOLUME_PREFIX}{cache_key}"


def _default_package_versions_relpath(env_profile: str, runtime: str) -> str:
    suffix = "packages.txt"
    if runtime == "python":
        suffix = "pip-freeze.txt"
    elif runtime == "r":
        suffix = "installed-packages.txt"
    return f"env_manifests/{env_profile}.{suffix}"


def _resolve_package_versions_path(relpath: Optional[str]) -> Optional[Path]:
    text = str(relpath or "").strip()
    if not text:
        return None
    return (image_catalog_root() / text).resolve()


@dataclass(frozen=True)
class ImageCatalogEntry:
    env_profile: str
    image_tag: str
    runtime: str
    base_image: str
    dockerfile: str
    build_context: str
    smoke_test: str
    package_versions_path: Optional[str]
    core_dependencies: tuple[str, ...]
    allowed_runtime_extras: tuple[str, ...]


@dataclass(frozen=True)
class FilteredRequirements:
    allowed_lines: tuple[str, ...]
    blocked_lines: tuple[str, ...]
    blocked_core_lines: tuple[str, ...]
    metadata_lines: tuple[str, ...]
    runtime_extra_lines: tuple[str, ...]
    env_cache_key: str
    env_cache_volume: str

    @property
    def allowed_text(self) -> str:
        return "\n".join(self.allowed_lines)


@dataclass(frozen=True)
class ResolvedExecutionEnvironment:
    skill_id: Optional[str]
    runtime: str
    env_profile: str
    env_image: str
    env_mode: str
    env_extra_requirements: tuple[str, ...]
    required_assets: tuple[str, ...]
    env_signature: str
    dockerfile: Optional[str]
    build_context: Optional[str]
    smoke_test: Optional[str]
    package_versions_path: Optional[str]
    core_dependencies: tuple[str, ...]
    allowed_runtime_extras: tuple[str, ...]
    env_cache_key: str
    env_cache_volume: str

    @property
    def is_python(self) -> bool:
        return self.runtime == "python"

    @property
    def is_shared(self) -> bool:
        return self.env_mode == "shared"


@lru_cache(maxsize=1)
def load_image_catalog() -> dict[str, ImageCatalogEntry]:
    if not _CATALOG_PATH.is_file():
        _logger.warning("image catalog missing: %s", _CATALOG_PATH)
        return {}
    data = json.loads(_CATALOG_PATH.read_text(encoding="utf-8"))
    out: dict[str, ImageCatalogEntry] = {}
    for item in data:
        if not isinstance(item, dict):
            continue
        profile = str(item.get("env_profile") or "").strip()
        if not profile:
            continue
        out[profile] = ImageCatalogEntry(
            env_profile=profile,
            image_tag=str(item.get("image_tag") or "").strip(),
            runtime=str(item.get("runtime") or "").strip(),
            base_image=str(item.get("base_image") or "").strip(),
            dockerfile=str(item.get("dockerfile") or "").strip(),
            build_context=str(item.get("build_context") or ".").strip() or ".",
            smoke_test=str(item.get("smoke_test") or "").strip(),
            package_versions_path=str(item.get("package_versions_path") or "").strip() or None,
            core_dependencies=_normalize_items(item.get("core_dependencies") or []),
            allowed_runtime_extras=_normalize_items(
                item.get("allowed_runtime_extras") or _DEFAULT_ALLOWED_RUNTIME_EXTRAS
            ),
        )
    return out


def get_image_catalog_entry(env_profile: Optional[str]) -> Optional[ImageCatalogEntry]:
    if not env_profile:
        return None
    return load_image_catalog().get(env_profile)


def resolve_execution_environment(skill_id: Optional[str]) -> Optional[ResolvedExecutionEnvironment]:
    rec = get_skill(skill_id or "")
    if rec is None:
        return None
    return resolve_execution_environment_for_record(rec)


def resolve_execution_environment_for_record(
    rec: WorkflowSkillRecord,
) -> Optional[ResolvedExecutionEnvironment]:
    runtime = (rec.runtime or "").strip()
    env_profile = (rec.env_profile or "").strip()
    env_image = (rec.env_image or "").strip()
    if not runtime or not env_profile:
        return None

    catalog_entry = get_image_catalog_entry(env_profile)
    if catalog_entry:
        if catalog_entry.runtime and runtime and catalog_entry.runtime != runtime:
            raise ValueError(
                f"workflow `{rec.skill_id}` runtime={runtime} 与 image catalog runtime={catalog_entry.runtime} 不一致"
            )
        if not env_image:
            env_image = catalog_entry.image_tag

    if not env_image:
        raise ValueError(f"workflow `{rec.skill_id}` 缺少 env_image，且 image catalog 无默认值")

    extra_requirements = _normalize_items(rec.env_extra_requirements)
    required_assets = _normalize_items(rec.required_assets)
    signature = _compute_env_signature(
        runtime=runtime,
        env_profile=env_profile,
        env_image=env_image,
        env_extra_requirements=extra_requirements,
        required_assets=required_assets,
    )
    core_dependencies = catalog_entry.core_dependencies if catalog_entry else ()
    allowed_runtime_extras = catalog_entry.allowed_runtime_extras if catalog_entry else _normalize_items(_DEFAULT_ALLOWED_RUNTIME_EXTRAS)
    package_versions_path = (
        catalog_entry.package_versions_path
        if catalog_entry and catalog_entry.package_versions_path
        else _default_package_versions_relpath(env_profile, runtime)
    )
    cache_key = _compute_env_cache_key(signature, extra_requirements)
    return ResolvedExecutionEnvironment(
        skill_id=rec.skill_id,
        runtime=runtime,
        env_profile=env_profile,
        env_image=env_image,
        env_mode=(rec.env_mode or "isolated").strip() or "isolated",
        env_extra_requirements=extra_requirements,
        required_assets=required_assets,
        env_signature=signature,
        dockerfile=catalog_entry.dockerfile if catalog_entry else None,
        build_context=catalog_entry.build_context if catalog_entry else None,
        smoke_test=catalog_entry.smoke_test if catalog_entry else None,
        package_versions_path=package_versions_path,
        core_dependencies=core_dependencies,
        allowed_runtime_extras=allowed_runtime_extras,
        env_cache_key=cache_key,
        env_cache_volume=env_cache_volume_name(cache_key),
    )


def resolve_environment_for_state(state: dict[str, Any], skill_id: Optional[str] = None) -> Optional[ResolvedExecutionEnvironment]:
    preferred = (skill_id or state.get("current_step_skill_id") or state.get("selected_skill_id") or "").strip()
    if preferred:
        return resolve_execution_environment(preferred)
    return None


def filter_profiled_requirements(
    env: ResolvedExecutionEnvironment,
    runtime_requirements_text: str,
) -> FilteredRequirements:
    metadata_lines = list(env.env_extra_requirements)
    metadata_names = {
        requirement_name_from_line(line)
        for line in metadata_lines
        if requirement_name_from_line(line)
    }
    core_names = {
        requirement_name_from_line(line)
        for line in env.core_dependencies
        if requirement_name_from_line(line)
    }
    allowlist_names = {
        canonicalize_requirement_name(name)
        for name in env.allowed_runtime_extras
        if canonicalize_requirement_name(name)
    }

    allowed_lines = list(metadata_lines)
    blocked_lines: list[str] = []
    blocked_core_lines: list[str] = []
    seen_allowed_names = set(metadata_names)

    for line in requirement_lines_from_text(runtime_requirements_text):
        req_name = requirement_name_from_line(line)
        if not req_name:
            blocked_lines.append(line)
            continue
        if req_name in core_names:
            blocked_lines.append(line)
            blocked_core_lines.append(line)
            continue
        if req_name in metadata_names:
            continue
        if req_name not in allowlist_names:
            blocked_lines.append(line)
            continue
        if req_name in seen_allowed_names:
            continue
        allowed_lines.append(line)
        seen_allowed_names.add(req_name)

    cache_key = _compute_env_cache_key(env.env_signature, allowed_lines)
    metadata_set = set(metadata_lines)
    runtime_extra_lines = [line for line in allowed_lines if line not in metadata_set]
    return FilteredRequirements(
        allowed_lines=tuple(allowed_lines),
        blocked_lines=tuple(blocked_lines),
        blocked_core_lines=tuple(blocked_core_lines),
        metadata_lines=tuple(metadata_lines),
        runtime_extra_lines=tuple(runtime_extra_lines),
        env_cache_key=cache_key,
        env_cache_volume=env_cache_volume_name(cache_key),
    )


def sync_environment_to_state(
    state: dict[str, Any],
    env: ResolvedExecutionEnvironment,
    *,
    container_id: Optional[str] = None,
    env_cache_key: Optional[str] = None,
    env_cache_volume: Optional[str] = None,
) -> None:
    state["docker_env_profile"] = env.env_profile
    state["docker_env_signature"] = env.env_signature
    state["docker_env_image"] = env.env_image
    state["docker_env_runtime"] = env.runtime
    state["docker_env_mode"] = env.env_mode
    state["docker_env_package_versions_path"] = env.package_versions_path
    state["docker_required_assets"] = list(env.required_assets)
    state["docker_env_cache_key"] = env_cache_key or env.env_cache_key
    state["docker_env_cache_volume"] = env_cache_volume or env.env_cache_volume
    if container_id:
        state["docker_container_id"] = container_id


def clear_environment_from_state(state: dict[str, Any]) -> None:
    for key in (
        "docker_env_profile",
        "docker_env_signature",
        "docker_env_image",
        "docker_env_runtime",
        "docker_env_mode",
        "docker_env_package_versions_path",
        "docker_required_assets",
        "docker_env_cache_key",
        "docker_env_cache_volume",
    ):
        state[key] = None


def describe_environment(env: Optional[ResolvedExecutionEnvironment]) -> str:
    if env is None:
        return "legacy"
    return (
        f"profile={env.env_profile} runtime={env.runtime} image={env.env_image} "
        f"mode={env.env_mode} signature={env.env_signature} cache={env.env_cache_key}"
    )


def ensure_asset_cache(env: Optional[ResolvedExecutionEnvironment]) -> Optional[Path]:
    if env is None or not env.required_assets:
        return None
    _ASSET_CACHE_ROOT.mkdir(parents=True, exist_ok=True)
    profile_root = _ASSET_CACHE_ROOT / env.env_profile
    profile_root.mkdir(parents=True, exist_ok=True)
    for asset in env.required_assets:
        (profile_root / asset).mkdir(parents=True, exist_ok=True)
    return profile_root


def asset_cache_status(env: Optional[ResolvedExecutionEnvironment]) -> list[tuple[str, bool, Path]]:
    if env is None or not env.required_assets:
        return []
    profile_root = _ASSET_CACHE_ROOT / env.env_profile
    out = []
    for asset in env.required_assets:
        asset_path = profile_root / asset
        out.append((asset, asset_path.exists(), asset_path))
    return out


def image_catalog_root() -> Path:
    return _CATALOG_PATH.parent


def env_manifests_root() -> Path:
    return _ENV_MANIFESTS_ROOT


def load_environment_package_versions_text(env: Optional[ResolvedExecutionEnvironment]) -> str:
    if env is None:
        return ""
    path = _resolve_package_versions_path(env.package_versions_path)
    if path is None or not path.is_file():
        return ""
    try:
        return path.read_text(encoding="utf-8").strip()
    except OSError:
        return ""


def format_environment_package_versions_for_prompt(env: Optional[ResolvedExecutionEnvironment]) -> str:
    if env is None:
        return ""
    package_versions = load_environment_package_versions_text(env)
    if not package_versions:
        return ""
    source_path = _resolve_package_versions_path(env.package_versions_path)
    source_note = str(source_path) if source_path else str(env.package_versions_path or "")
    return (
        "\n【当前 Docker 环境完整包版本清单】\n"
        f"env_profile: {env.env_profile}\n"
        f"env_image: {env.env_image}\n"
        f"source: {source_note}\n"
        "以下内容是当前容器环境中已安装包的完整版本列表。生成代码时必须以此为准，不要猜测版本，"
        "不要调用与该环境不兼容或已不存在的 API。\n"
        "```text\n"
        f"{package_versions}\n"
        "```\n"
    )
