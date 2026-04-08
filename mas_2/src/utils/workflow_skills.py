"""
工作流技能（workflows/*/SKILL.md）注册表。

用于 Supervisor 选型、Code Dev 提示词与 Docker 挂载路径解析。
"""
from __future__ import annotations

import logging
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict, List, Optional

_logger = logging.getLogger(__name__)

# 第一期试点：历史端到端验收用例；目录展示仍标【试点】。凡已注册技能均会挂载 workflows 目录（只读 /app/workflow）。
PILOT_SKILL_IDS = frozenset(
    {
        "gwas-to-function-twas",
        "scrna-trajectory-inference",
        "scrnaseq-scanpy-core-analysis",
    }
)

SCANPY_CORE_SKILL_ID = "scrnaseq-scanpy-core-analysis"

_UTILS_DIR = Path(__file__).resolve().parent
_MAS2_ROOT = _UTILS_DIR.parent.parent
_WORKFLOWS_ROOT = _MAS2_ROOT / "workflows"


@dataclass(frozen=True)
class WorkflowSkillRecord:
    skill_id: str
    name: str
    category: str
    short_description: str
    root_path: Path
    raw_meta: Dict[str, Any]


def _split_frontmatter(md_text: str) -> tuple[Dict[str, Any], str]:
    """解析 SKILL.md 的 YAML frontmatter。"""
    try:
        import yaml
    except ImportError:
        yaml = None  # type: ignore

    if not md_text.startswith("---"):
        return {}, md_text

    parts = md_text.split("---", 2)
    if len(parts) < 3:
        return {}, md_text

    raw = parts[1].strip()
    body = parts[2].lstrip("\n")
    if yaml is None:
        meta: Dict[str, Any] = {}
        for line in raw.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if ":" in line:
                k, v = line.split(":", 1)
                meta[k.strip()] = v.strip().strip('"').strip("'")
        return meta, body

    loaded = yaml.safe_load(raw)
    if not isinstance(loaded, dict):
        return {}, body
    return loaded, body


def _record_from_skill_path(skill_md: Path) -> Optional[WorkflowSkillRecord]:
    try:
        text = skill_md.read_text(encoding="utf-8")
    except OSError:
        return None
    meta, _ = _split_frontmatter(text)
    skill_id = str(meta.get("id") or skill_md.parent.name).strip()
    if not skill_id:
        return None
    name = str(meta.get("name") or skill_id)
    category = str(meta.get("category") or "")
    short_description = str(meta.get("short-description") or meta.get("short_description") or "")
    root = skill_md.parent
    return WorkflowSkillRecord(
        skill_id=skill_id,
        name=name,
        category=category,
        short_description=short_description,
        root_path=root,
        raw_meta=meta,
    )


@lru_cache(maxsize=1)
def list_skills() -> List[WorkflowSkillRecord]:
    """扫描 workflows/*/SKILL.md，返回全部技能元数据。"""
    out: List[WorkflowSkillRecord] = []
    if not _WORKFLOWS_ROOT.is_dir():
        return out
    for skill_md in sorted(_WORKFLOWS_ROOT.glob("*/SKILL.md")):
        rec = _record_from_skill_path(skill_md)
        if rec:
            out.append(rec)
    return out


def get_skill(skill_id: str) -> Optional[WorkflowSkillRecord]:
    if not skill_id:
        return None
    for rec in list_skills():
        if rec.skill_id == skill_id:
            return rec
    return None


def resolve_workflow_root(skill_id: str) -> Optional[str]:
    """返回宿主机上该技能根目录绝对路径；不存在则 None。"""
    rec = get_skill(skill_id)
    if rec and rec.root_path.is_dir():
        return str(rec.root_path.resolve())
    candidate = _WORKFLOWS_ROOT / skill_id
    if candidate.is_dir():
        return str(candidate.resolve())
    return None

def get_workflows_root() -> str:
    """返回宿主机 workflows 根目录绝对路径。"""
    return str(_WORKFLOWS_ROOT.resolve())


def format_skills_catalog_for_prompt(max_items: int = 64) -> str:
    """供 Supervisor 注入的简短目录（避免全文）。"""
    lines: List[str] = []
    for rec in list_skills()[:max_items]:
        pilot = "【试点】" if rec.skill_id in PILOT_SKILL_IDS else ""
        lines.append(
            f"- id=`{rec.skill_id}` {pilot}| {rec.name} ({rec.category}) — {rec.short_description}"
        )
    if not lines:
        return "（未找到 workflows/*/SKILL.md，请检查 mas_2/workflows 目录）"
    return "\n".join(lines)


def format_skill_injection_for_code_dev(skill_id: Optional[str], max_chars: int = 4000) -> str:
    """为 Code Dev 注入：提供元数据和物理路径，指示 Agent 自主调用工具探索目录和读取文件。"""
    if not skill_id:
        return ""
    rec = get_skill(skill_id)
    if not rec:
        return f"\n【workflow skill】skill_id={skill_id}（未在注册表中找到，请仅依赖任务描述与 RAG）\n"
    
    root_path_str = str(rec.root_path.resolve()).replace('\\', '/')
    
    meta_lines = [
        f"skill_id: {rec.skill_id}",
        f"name: {rec.name}",
        f"category: {rec.category}",
        f"short-description: {rec.short_description}",
        f"host_path: {root_path_str} (宿主机绝对路径，用于工具调用读取本地文件)"
    ]
    
    _logger.debug("skill_injection for code_dev tool exploration skill_id=%s meta=%s", skill_id, meta_lines)
    
    return (
        "\n【当前步骤绑定的 Workflow Skill】\n"
        + "\n".join(meta_lines)
        + "\n"
    )


def format_skill_for_critic(skill_id: Optional[str]) -> str:
    if not skill_id:
        return ""
    rec = get_skill(skill_id)
    if not rec:
        return f"\n【技能上下文】用户计划指定 skill_id={skill_id}（未在注册表中找到）。\n"
    return (
        f"\n【技能上下文】当前步骤 skill_id=`{rec.skill_id}`（{rec.name}）。"
        "若技能文档要求使用其目录内参考脚本，容器内通常只读挂载为 /app/workflow，属正常。\n"
        f"验收时若 acceptance 与技能文档一致，且执行成功，应通过；不要按不相关的 Scanpy 默认流程驳回。\n"
    )


def should_mount_workflow_in_docker(skill_id: Optional[str]) -> bool:
    """凡已注册且目录存在的技能均挂载 workflows 根目录到容器 /app/workflow（只读）。"""
    if not skill_id:
        return False
    return resolve_workflow_root(skill_id) is not None
